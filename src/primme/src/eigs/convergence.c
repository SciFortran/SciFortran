/*******************************************************************************
 * Copyright (c) 2018, College of William & Mary
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the College of William & Mary nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * PRIMME: https://github.com/primme/primme
 * Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
 *******************************************************************************
 * File: convergence.c
 *
 * Purpose - Checks for convergence of the block vectors.
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/convergence.c"
#endif

#include "numerical.h"
#include "template_normal.h"
#include "common_eigs.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "convergence.h"
#include "ortho.h"
#include "auxiliary_eigs.h"
#include "auxiliary_eigs_normal.h"
#endif

#ifdef SUPPORTED_TYPE

/*******************************************************************************
 * Subroutine check_convergence - This procedure checks the block vectors for  
 *    convergence.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * X              The Ritz vectors in the block
 * nLocal, ldX    The number of rows and the leading dimension of X
 * givenX         Whether X is provided
 * ldR            The leading dimension of R
 * givenR        Whether R is provided
 * evecs          The locked vectors
 * numLocked      The number of columns in evecs besides numOrthoConst
 * ldevecs        The leading dimension of evecs
 * left, right    Range of vectors to be checked for convergence
 * blockNorms     Residual norms of the Ritz vectors starting from left
 * hVals          The Ritz values
 * practConvCheck Disable (-1) or enforce (1) the practically convergence checking
 * VtBV           evecs'*B*evecs
 * ctx            Structure containing various solver parameters
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * R              The residual vectors of the Ritz vectors in the block
 *                (this routine may remove the evecs directions in some vectors)
 * flags          Array indicating which eigenvectors have converged     
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * reset          flag to reset V and W in the next restart
 ******************************************************************************/

TEMPLATE_PLEASE
int check_convergence_Sprimme(SCALAR *X, PRIMME_INT ldX, int givenX, SCALAR *R,
      PRIMME_INT ldR, int givenR, SCALAR *evecs, int numLocked,
      PRIMME_INT ldevecs, SCALAR *Bevecs, PRIMME_INT ldBevecs, HSCALAR *VtBV,
      int ldVtBV, int left, int right, int *flags, HREAL *blockNorms,
      HEVAL *hVals, int *reset, int practConvCheck, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i;                  /* Loop variable                                      */
   int numToProject;       /* Number of vectors with potential accuracy problem  */
   int *toProject = NULL; /* Indices from left with potential accuracy problem  */
   double tol;             /* Residual tolerance                                 */
   double attainableTol=0; /* Used in locking to check near convergence problem  */
   int isConv;             /* return of convTestFun                              */

   CHKERR(Num_malloc_iprimme(right-left, &toProject, ctx));

   /* -------------------------------------------- */
   /* Tolerance based on our dynamic norm estimate */
   /* -------------------------------------------- */

   double eps_matrix;
   CHKERR(machineEpsMatrix_Sprimme(&eps_matrix, ctx));
   tol = max(
         eps_matrix * problemNorm_Sprimme(1, primme), primme->stats.maxConvTol);

   /* ---------------------------------------------------------------------- */
   /* If locking, set tol beyond which we need to check for accuracy problem */
   /* ---------------------------------------------------------------------- */
   if (primme->locking) {
      attainableTol = sqrt((double)(primme->numOrthoConst+numLocked))*tol;
   }

   /* ----------------------------------------------------------------- */
   /* Determine which Ritz vectors have converged < tol and flag them.  */
   /* ----------------------------------------------------------------- */

   numToProject = 0;
   for (i=left; i < right; i++) {
       
      /* Refine doesn't order the pairs considering closest_leq/gep. */
      /* Then ignore values so that value +-residual is completely   */
      /* outside of the desired region.                              */

#ifdef USE_HERMITIAN
      double targetShift = primme->numTargetShifts > 0
                                 ? primme->targetShifts[min(primme->initSize,
                                         primme->numTargetShifts - 1)]
                                 : 0.0;

      if ((primme->target == primme_closest_leq
               && hVals[i]-blockNorms[i-left] > targetShift) ||
            (primme->target == primme_closest_geq
             && hVals[i]+blockNorms[i-left] < targetShift)) {
         flags[i] = UNCONVERGED;
         continue;
      }
#endif

      if (blockNorms[i-left] <= primme->stats.maxConvTol) {
         flags[i] = CONVERGED;
         continue;
      }

      CHKERR(convTestFun_Sprimme(hVals[i], X ? &X[ldX * (i - left)] : NULL,
            givenX, blockNorms[i - left], &isConv, ctx));

      if (isConv) {
         flags[i] = CONVERGED;
      }

      /* ----------------------------------------------------------------- */
      /* If residual norm is around the bound of the error in the          */
      /* residual norm, then stop converging this value and force reset    */
      /* of V and W in the next restart.                                   */
      /* ----------------------------------------------------------------- */

      else if (blockNorms[i-left] <= primme->stats.estimateResidualError && reset) {
         flags[i] = SKIP_UNTIL_RESTART;
         *reset = 1;
      }

      /* ----------------------------------------------------------------- */
      /* If locking there may be an accuracy problem close to convergence. */
      /* Check if there is danger if R is provided. If the Ritz vector was */
      /* flagged practically converged before and R is not provided then   */
      /* consider converged still.                                         */
      /* ----------------------------------------------------------------- */

      else if (primme->locking && numLocked > 0 && practConvCheck >= 0) {
         if (givenR && blockNorms[i - left] < attainableTol) {
            toProject[numToProject++] = i-left;
         }
         else if (flags[i] != PRACTICALLY_CONVERGED) {
            flags[i] = UNCONVERGED;
         }
      }

      else {
         flags[i] = UNCONVERGED;
      }
   }

   /* --------------------------------------------------------------- */
   /* Project the TO_BE_PROJECTED residuals and check for practical   */
   /* convergence among them.                                         */
   /* --------------------------------------------------------------- */

   if (numToProject > 0) {
      CHKERR(check_practical_convergence(R, ldR, evecs,
            primme->numOrthoConst + numLocked, ldevecs, Bevecs, ldBevecs, left,
            toProject, numToProject, flags, blockNorms, tol, VtBV, ldVtBV,
            ctx));
   }

   CHKERR(Num_free_iprimme(toProject, ctx));

   return 0;

}

/*******************************************************************************
 * Subroutine check_practical_convergence(): for pairs whose residual norm is
 *    less than tol*sqrt(numConverged) but greater than tol, they will be
 *    flagged as converged if || (I-BQQ')R(i) || < tol, where Q are the locked
 *    vectors.
 *
 *    NOTE: the routine removes evecs directions from the residual vectors,
 *          but blockNorms isn't changed.
 *
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * nLocal          The local length of the vectors in the basis
 * evecs           The locked eigenvectors
 * evecsSize       The number of locked eigenvectors
 * ldevecs         The leading dimension of evecs
 * left            Base index indicating which flags are to be recomputed
 * iev             Indices of flags to recompute
 * numToProject    Size of iev
 * blockNorms      The norms of the residual vectors starting by index 'left'
 * tol             The required convergence tolerance
 * VtBV            evecs'*B*evecs
 * ctx             Structure containing various solver parameters
 *
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * R               The residual vectors
 * ldR             The leading dimension of R
 * flags           Indicates which Ritz pairs have converged
 ******************************************************************************/

STATIC int check_practical_convergence(SCALAR *R, PRIMME_INT ldR, SCALAR *evecs,
      int evecsSize, PRIMME_INT ldevecs, SCALAR *Bevecs, PRIMME_INT ldBevecs,
      int left, int *iev, int numToProject, int *flags, HREAL *blockNorms,
      double tol, HSCALAR *VtBV, int ldVtBV, primme_context ctx) {

   int i;
   HREAL *norms;

   /* R = (I-Bevecs*evecs')*R */
   /* newResiduals(i) = || (I-Bevecs*evecs')*R(i) || */

   CHKERR(Num_malloc_RHprimme(numToProject, &norms, ctx));

   CHKERR(ortho_single_iteration_Sprimme(evecs, evecsSize, ldevecs,
         Bevecs ? Bevecs : evecs, Bevecs ? ldBevecs : ldevecs, VtBV, ldVtBV, R,
         iev, numToProject, ldR, norms, ctx));

   for (i=0; i < numToProject; i++) {

      /* ------------------------------------------------------------------ */
      /* NOTE: previous versions than 2.0 used the next criterion instead:  */
      /*                                                                    */
      /* if (normDiff >= tol && normPr < tol*tol/blockNorm/2)               */
      /*                                                                    */
      /* Due to rounding-off errors in computing normDiff, the test         */
      /* normDiff >= tol may fail even it is satisfied in exact arithmetic. */
      /* ------------------------------------------------------------------ */

      blockNorms[iev[i]] = norms[i];
     
      if (norms[i] <= tol) {
         PRINTF(5, " PRACTICALLY_CONVERGED %d norm(I-BQQt)r %e",
               left + iev[i], (double)blockNorms[i]);
         flags[left+iev[i]] = PRACTICALLY_CONVERGED;
      }
      else {
         flags[left+iev[i]] = UNCONVERGED;
      }
   }

   CHKERR(Num_free_RHprimme(norms, ctx));

   return 0;
}

#endif /* SUPPORTED_TYPE */
