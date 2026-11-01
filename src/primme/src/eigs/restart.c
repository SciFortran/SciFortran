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
 * File: restart.c
 *
 * Purpose - Restart V and related matrices (eg. W, H, Q, R, QtV...).
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/restart.c"
#endif

#include "numerical.h"
#include "template_normal.h"
#include "common_eigs.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "auxiliary_eigs.h"
#include "auxiliary_eigs_normal.h"
#include "restart.h"
#include "ortho.h"
#include "solve_projection.h"
#include "factorize.h"
#include "update_projection.h"
#include "update_W.h"
#include "convergence.h"
#endif

#ifdef SUPPORTED_TYPE

/*******************************************************************************
 * Subroutine: restart - This routine replaces V with V*c, some subset
 *             of the Ritz vectors of the current and the previous iteration.
 *             Related bases and matrices are updated accordingly.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 *
 * nLocal           Number of rows of V, W, Q, evecs and evecsHat assigned to the node
 *
 * ldV              The leading dimension of V, W, Q, evecs and evecsHat
 *
 * basisSize        The number of columns in V, W and Q
 *
 * resNorms         The residual norms of the converged eigenpairs
 *
 * prevhVecs        Coefficient vectors retained from the previous iteration
 *
 * ldprevhVecs      The leading dimension of prevhVecs
 *
 * numGuesses       Number of remaining initial guesses
 *
 * rwork            Real work array
 *
 * rworkSize        Must be of size 
 *
 * iwork            Integer work array
 *                  
 * primme           Structure containing various solver parameters
 *
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * restartSizeOutput Output the number of columns of the restarted V.
 *
 * V                The orthonormal basis. After restart, contains Ritz vectors
 *                  plus the orthogonal components from prevhVecs 
 *
 * W                A*V
 *
 * BV               B*V
 *
 * hU               The left singular vectors of R or the eigenvectors of QtV/R
 *
 * ldhU             The leading dimension of the input hU
 *
 * newldhU          The leading dimension of the output hU
 *
 * hVecs            The coefficient vectors
 *
 * ldhVecs          The leading dimension of the input hVecs
 *
 * newldhVecs       The leading dimension of the output hVecs
 *
 * hVals            The Rayleigh quotient of candidates
 *
 * hSVals           The singular values of R
 *
 * flags            Array indicating the convergence of the Ritz vectors
 *
 * iev              Array of size blockSize indicating which Ritz vectors are
 *                  targeted in the block
 *
 * ievSize          The length of iev
 *
 * blockNorms       The residual norms of the eigenpairs in the block
 *
 * evecs            The converged Ritz vectors. Without locking, all converged
 *                  eigenvectors are copied from V to evecs if skew projections
 *                  are required
 *
 * Bevecs           B*evecs
 *
 * evals            The converged Ritz values
 *
 * resNorms         The residual norms of the converged eigenpairs
 *
 * evecsHat         K^{-1}evecs
 *
 * prevRitzVals     Projected values retained from the previous iteration
 *
 * numPrevRitzVals  Length of the vector prevRitzVals
 *
 * ldevecsHat       The leading dimension of evecsHat
 *
 * M, ldM           evecs'*evecsHat and the leading dimension of M
 *
 * Mfact            The factorization of M
 *
 * ldMfact          The leading dimension of Mfact
 *
 * ipivot           The pivot array of the Mfact factorization
 *
 * H                The projection V'*A*V
 *
 * ldH              The leading dimension of H
 *
 * VtBV             The projection V'*B*V
 *
 * ldVtBV           The leading dimension of VtBV
 *
 * fVtBV            The Cholesky factor of VtBV
 *
 * ldfVtBV          The leading dimension of fVtBV
 * 
 * Q, R             The factors of the QR decomposition of (A - targetShift*B)*V
 *
 * ldQ, ldR         The leading dimension of Q and R
 *
 * fQtQ             The Cholesky factor of QtQ
 *
 * ldfQtQ           The leading dimension of fQtQ
 * 
 * numConverged     The number of converged eigenpairs
 *
 * numLocked        The number of locked eigenpairs
 *
 * lockedFlags      The flags of the locked pairs
 *
 * numConvergedStored The # of converged vectors copied to evecs
 *
 * nprevhVecs       The number of columns of prevhVecs
 *
 * targetShiftIndex The target shift used in (A - targetShift*B) = Q*R
 *
 * numArbitraryVecs The number of columns of hVecsRot
 *
 * hVecsRot         hVecs = hV*hVecsRot, where hV are the original coefficient
 *                  vectors returned by solve_H
 *
 * ldhVecsRot       The leading dimension of hVecsRot
 *
 * restartsSinceReset Number of restarts since last reset of V and W
 *
 *
 * Return value
 * ------------
 * error code
 *       
 ******************************************************************************/
 
TEMPLATE_PLEASE
int restart_Sprimme(SCALAR *V, SCALAR *W, SCALAR *BV, PRIMME_INT nLocal,
      int basisSize, PRIMME_INT ldV, HEVAL *hVals, HREAL *hSVals, int *flags,
      int *iev, int *ievSize, HREAL *blockNorms, SCALAR *evecs,
      PRIMME_INT ldevecs, SCALAR *Bevecs, PRIMME_INT ldBevecs, int *evecsPerm,
      HEVAL *evals, HREAL *resNorms, SCALAR *evecsHat, PRIMME_INT ldevecsHat,
      HSCALAR *M, int ldM, HSCALAR *Mfact, int ldMfact, int *ipivot,
      int *numConverged, int *numLocked, int *lockedFlags,
      int *numConvergedStored, HSCALAR *prevhVecs, int nprevhVecs,
      int ldprevhVecs, int numGuesses, HEVAL *prevRitzVals,
      int *numPrevRitzVals, HSCALAR *H, int ldH, HSCALAR *VtBV, int ldVtBV,
      HSCALAR *fVtBV, int ldfVtBV, SCALAR *Q, PRIMME_INT ldQ, HSCALAR *R,
      int ldR, HSCALAR *QtV, int ldQtV, HSCALAR *QtQ, int ldQtQ, HSCALAR *fQtQ,
      int ldfQtQ, HSCALAR *hU, int ldhU, int newldhU, HSCALAR *hVecs,
      int ldhVecs, int newldhVecs, int *restartSizeOutput,
      int *targetShiftIndex, int *numArbitraryVecs, HSCALAR *hVecsRot,
      int ldhVecsRot, int *restartsSinceReset, double startTime,
      primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i;                   /* Loop indices */
   int restartSize;         /* Basis size after restarting                   */
   int indexOfPreviousVecsBeforeRestart=0;/* descriptive enough name, isn't? */

   /* ----------------------------------------------------------- */
   /* Remove the SKIP_UNTIL_RESTART flags.                        */
   /* ----------------------------------------------------------- */

   for (i=0, *numConverged=*numLocked; i<basisSize; i++) {
      if (flags[i] == SKIP_UNTIL_RESTART) {
         flags[i] = UNCONVERGED;
      }
      else if (flags[i] != UNCONVERGED &&
            /* Don't check more than numEvals */
               *numConverged < primme->numEvals &&
            /* Check only the first pairs, except if finding closest_leq/geq  */
            /* because with refined extraction the pairs may not be ordered   */
            /* by this criterion.                                             */
               (i < primme->numEvals-*numLocked
                || primme->target == primme_closest_geq
                || primme->target == primme_closest_leq)) {
         (*numConverged)++;
      }
   }

   /* ----------------------------------------------------------- */
   /* Special case: If (basisSize+numLocked) is the entire space, */
   /* then everything should be converged. Do not test, just flag */
   /* everything as converged. But only do that when multiple     */
   /* shifts are not involved.                                    */
   /* ----------------------------------------------------------- */

   int numPrevRetained = primme->restartingParams.maxPrevRetain;
   if (!primme->locking &&
         basisSize + *numLocked + primme->numOrthoConst >= primme->n) {
      restartSize = basisSize;
      numPrevRetained = 0;
   }
   /* --------------------------------------------------------------------- */
   /* If basis isn't full, restart with the current basis size.             */
   /* --------------------------------------------------------------------- */
   else if (basisSize <= primme->maxBasisSize - primme->maxBlockSize) {
      restartSize = basisSize;
      numPrevRetained = 0;
   }
   else {
      restartSize = min(basisSize, primme->minRestartSize);
   }

   /* ----------------------------------------------------------------------- */
   /* Insert as many initial guesses as eigenpairs have converged.            */
   /* Leave sufficient restarting room in the restarted basis so that to      */
   /* insert (in main_iter) as many initial guesses as the number of          */
   /* eigenpairs that converged.                                              */
   /* ----------------------------------------------------------------------- */

   restartSize -= min(min(numGuesses, *numConverged-*numLocked), restartSize);

   /* ----------------------------------------------------------------------- */
   /* Limit restartSize so that it plus 'to be locked' plus previous Ritz     */
   /* vectors do not exceed basisSize.                                        */
   /* ----------------------------------------------------------------------- */

   if (primme->locking)
      restartSize = min(restartSize, basisSize-(*numConverged-*numLocked));

   /* ----------------------------------------------------------------------- */
   /* Limit the number of previous retained vectors such that the final basis */
   /* size isn't larger than the current basis size.                          */
   /* ----------------------------------------------------------------------- */

   numPrevRetained = max(0, min(min(
         numPrevRetained,
         primme->maxBasisSize - restartSize - 1),
         primme->n - restartSize - *numConverged - primme->numOrthoConst));

   /* ----------------------------------------------------------------------- */
   /* Restarting with a small number of coefficient vectors from the previous */
   /* iteration can accelerate convergence.  The previous                     */
   /* coefficient vectors must be combined with the current coefficient       */
   /* vectors by first orthogonalizing the previous ones versus the current   */
   /* restartSize ones.  The orthogonalized previous vectors are then         */
   /* inserted into the hVecs array at hVecs(:,indexOfPreviousVecs).          */
   /* When using refined, avoid to replace linear dependent previous coeffi-  */
   /* cient vectors by random vectors, it may increase the norm of R after    */
   /* restarting.                                                             */
   /* ----------------------------------------------------------------------- */

   int indexOfPreviousVecs; /* Column index in hVecs with previous vecs      */
   if (primme->locking)
      indexOfPreviousVecs = restartSize+*numConverged-*numLocked;
   else
      indexOfPreviousVecs = restartSize;
   indexOfPreviousVecsBeforeRestart = indexOfPreviousVecs;

   int nLocked = primme->numOrthoConst + *numLocked;
   {
      CHKERR(ortho_coefficient_vectors_Sprimme(hVecs, basisSize, ldhVecs,
            indexOfPreviousVecs,
            VtBV ? &VtBV[nLocked * ldVtBV + nLocked] : NULL, ldVtBV, prevhVecs,
            nprevhVecs, ldprevhVecs, flags, &numPrevRetained, ctx));
   }

   /* Restart V and W, and compute X and residual vectors for next candidates
    */

   int *restartPerm; /* Permutation of hVecs used to restart V        */
   int *hVecsPerm;   /* Permutation of hVecs to sort as primme.target */
   CHKERR(Num_malloc_iprimme(basisSize, &restartPerm, ctx));
   CHKERR(Num_malloc_iprimme(basisSize, &hVecsPerm, ctx));

   if (!primme->locking) {
      CHKERR(restart_soft_locking_Sprimme(&restartSize, V, W, BV, nLocal,
            basisSize, ldV, hVecs, ldhVecs, restartPerm, hVals, flags, iev,
            ievSize, blockNorms, evecs, ldevecs, Bevecs, ldBevecs, evals,
            resNorms, numConverged, numPrevRetained, &indexOfPreviousVecs,
            hVecsPerm, H, ldH, VtBV, ldVtBV, ctx));
   } else {
      CHKERR(restart_locking_Sprimme(&restartSize, V, W, BV, nLocal, basisSize,
            ldV, hVecs, ldhVecs, restartPerm, hVals, flags, iev, ievSize,
            blockNorms, evecs, ldevecs, Bevecs, ldBevecs, evals, numConverged,
            numLocked, resNorms, lockedFlags, evecsPerm, numPrevRetained,
            &indexOfPreviousVecs, hVecsPerm, H, ldH, VtBV, ldVtBV, startTime,
            ctx));
   }

   /* Update fVtBV */

   if (fVtBV) {
      int newnLocked = primme->numOrthoConst + *numLocked;
      CHKERR(update_cholesky_Sprimme(VtBV, ldVtBV, fVtBV, ldfVtBV, nLocked,
            newnLocked + restartSize, ctx));
   }

   /* Rearrange prevRitzVals according to restartPerm */

   if (primme->target != primme_smallest && primme->target != primme_largest) {
      if (*numPrevRitzVals > 0) {
         for (i = *numPrevRitzVals; i < basisSize; i++)
            prevRitzVals[i] = prevRitzVals[*numPrevRitzVals - 1];

         CHKERR(KIND(permute_vecs_RHprimme, permute_vecs_SHprimme)(
               prevRitzVals, 1, basisSize, 1, restartPerm, ctx));
      }
      for (i = 0; i < restartSize; i++) {
         if (restartPerm[i] >= *numPrevRitzVals) prevRitzVals[i] = hVals[i];
      }
      CHKERR(KIND(permute_vecs_RHprimme, permute_vecs_SHprimme)(
            prevRitzVals, 1, restartSize, 1, hVecsPerm, ctx));
      *numPrevRitzVals = restartSize;
   }

   if (newldhVecs == 0) newldhVecs = restartSize;
   if (newldhU == 0) newldhU = restartSize;
   CHKERR(restart_projection_Sprimme(V, ldV, W, ldV, BV, ldV, H, ldH, VtBV,
         ldVtBV, Q, ldQ, nLocal, R, ldR, QtV, ldQtV, QtQ, ldQtQ, fQtQ, ldfQtQ,
         hU, ldhU, newldhU, indexOfPreviousVecsBeforeRestart, hVecs, ldhVecs,
         newldhVecs, hVals, hSVals, restartPerm, hVecsPerm, restartSize,
         basisSize, numPrevRetained, indexOfPreviousVecs, evecs,
         numConvergedStored, ldevecs, Bevecs, ldBevecs, evecsHat, ldevecsHat, M,
         ldM, Mfact, ldMfact, ipivot, targetShiftIndex, *numConverged,
         numArbitraryVecs, hVecsRot, ldhVecsRot, ctx));
   CHKERR(Num_free_iprimme(restartPerm, ctx));

   /* If all request eigenpairs converged, move the converged vectors at the  */
   /* beginning of V                                                          */

   if (*numConverged >= primme->numEvals && !primme->locking) {
      CHKERR(permute_vecs_Sprimme(V, nLocal, restartSize, ldV, hVecsPerm, ctx));
      CHKERR(permute_vecs_Sprimme(W, nLocal, restartSize, ldV, hVecsPerm, ctx));
      if (BV) {
         CHKERR(permute_vecs_Sprimme(
               BV, nLocal, restartSize, ldV, hVecsPerm, ctx));
      }
   }
   CHKERR(Num_free_iprimme(hVecsPerm, ctx));

   *restartSizeOutput = restartSize; 

   /* If VtBV is computed, estimate the orthogonality error as                */
   /* ||I - V'*B*V||_F * problemNorm. Use this value to estimate the residual error */
   /* and to bound the minimum residual error norm the solver can converge    */

   HREAL fn = 0.0; /* = max(|I - VtBV|_F, |I - QtQ|_F) * problemNorm */

   if (ctx.procID == 0) {
      if (VtBV) {
         HREAL n = 0;
         int i,j, nVtBV = primme->numOrthoConst + *numLocked + restartSize;
         for (i = 0; i < nVtBV; i++) {
            for (j = 0; j < i; j++)
               n += 2 * REAL_PART(CONJ(VtBV[i * ldVtBV + j]) *
                                  VtBV[i * ldVtBV + j]) /
                    ABS(VtBV[i * ldVtBV + i]) / ABS(VtBV[j * ldVtBV + j]);
         }
         fn = sqrt(n);
      }

      if (QtQ) {
         HREAL n = 0;
         int i,j;
         for (i = 0; i < restartSize; i++) {
            for (j = 0; j < i; j++)
               n += 2 * REAL_PART(CONJ(QtQ[i * ldQtQ + j]) *
                                  QtQ[i * ldQtQ + j]) /
                    ABS(QtQ[i * ldQtQ + i]) / ABS(QtQ[j * ldQtQ + j]);
         }
         fn = max(fn, sqrt(n));
      }
   }
   CHKERR(broadcast_RHprimme(&fn, 1, ctx));
      
   if (fn > 0.0) {
      if (*restartsSinceReset <= 1) {
         PRINTF(5, "Orthogonalization level: %g", (double)fn);
         primme->stats.maxConvTol = max(primme->stats.maxConvTol,
               fn * primme->stats.estimateLargestSVal);
      }
      primme->stats.estimateResidualError = sqrt((double)*restartsSinceReset) *
                                            fn * problemNorm_Sprimme(1, primme);
   }
   else {
      double eps_orth;
      CHKERR(machineEpsOrth_Sprimme(&eps_orth, ctx));
      primme->stats.estimateResidualError =
            2 * sqrt((double)*restartsSinceReset) * eps_orth *
            problemNorm_Sprimme(1, primme);
   }

   /* Check VtBV */

   if (0) {
      HSCALAR *VtBV0 = NULL;
      int n = primme->numOrthoConst + *numLocked + restartSize;
      CHKERR(Num_malloc_SHprimme(n * n, &VtBV0, ctx));
      CHKERR(Num_gemm_ddh_Sprimme("C", "N", primme->numOrthoConst + *numLocked,
            primme->numOrthoConst + *numLocked, primme->nLocal, 1.0, evecs,
            ldevecs, evecs, ldevecs, 0.0, VtBV0, n, ctx));
      CHKERR(Num_gemm_ddh_Sprimme("C", "N", primme->numOrthoConst + *numLocked,
            restartSize, primme->nLocal, 1.0, evecs, ldevecs, V, ldV, 0.0,
            &VtBV0[n * (primme->numOrthoConst + *numLocked)], n, ctx));
      CHKERR(Num_gemm_ddh_Sprimme("C", "N", restartSize, restartSize,
            primme->nLocal, 1.0, V, ldV, V, ldV, 0.0,
            &VtBV0[n * (primme->numOrthoConst + *numLocked) +
                   primme->numOrthoConst + *numLocked],
            n, ctx));
      CHKERR(globalSum_SHprimme(VtBV0, n*n, ctx));
      int i,j;
      for (i = 0; i < n; i++) {
         for (j = 0; j <= i; j++) {
            assert(ABS(VtBV[ldVtBV * i + j] - VtBV0[n * i + j]) <
                   MACHINE_EPSILON * 10);
         }
      }
      CHKERR(Num_free_SHprimme(VtBV0, ctx));
   }

   /* Check H */

   if (0) {
      HSCALAR *H0 = NULL;
      CHKERR(Num_malloc_SHprimme(restartSize * restartSize, &H0, ctx));
      CHKERR(Num_gemm_ddh_Sprimme("C", "N", restartSize, restartSize,
            primme->nLocal, 1.0, W, ldV, V, ldV, 0.0, H0, restartSize, ctx));
      CHKERR(globalSum_SHprimme(H0, restartSize * restartSize, ctx));
      int i, j;
      for (i = 0; i < restartSize; i++) {
         for (j = 0; j <= i; j++) {
            assert(ABS(H[ldH * i + j] - H0[restartSize * i + j]) <
                   problemNorm_Sprimme(1, primme) * MACHINE_EPSILON * 10);
         }
      }
      CHKERR(Num_free_SHprimme(H0, ctx));
   }
   return 0;
}

/*******************************************************************************
 * Subroutine: restart_soft_locking - This routine replaces V with V*c, some
 *             subset of hVecs. Also it may include components from vectors 
 *             from the previous iteration.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 *
 * nLocal           Number of rows of V assigned to the node
 *
 * ldV              The leading dimension of V and W
 *
 * hR               The factors of the QR decomposition of (A - targetShift*B)*V
 *
 * ldhR             The leading dimension of Q and R
 *
 * hU               The eigenvectors of QtV/R
 *
 * ldhU             The leading dimension of hU
 *
 * basisSize        Size of the basis V
 *
 * resNorms         The residual norms of the converged eigenpairs
 *
 * previousHVecs    Coefficient vectors retained from the previous iteration
 *
 * rwork            Real work array
 *
 * rworkSize        Must be of size 
 *
 * iwork            Integer work array
 *                  
 * primme           Structure containing various solver parameters
 *
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * restartSize      The number of vectors to restart with.
 *
 * V                The orthonormal basis. After restart, contains Ritz vectors
 *                  plus the orthogonal components from numPrevRetained Ritz 
 *                  vectors from the penultimate step.
 *
 * W                A*V
 *
 * BV               B*V
 *
 * X                Reference to the Ritz vectors of the eigenpairs in the block
 *
 * R                Reference to the residual vectors of the eigenpairs in the block
 *
 * hVecs, ldhVecs   The eigenvectors of H and the leading dimension of H
 *
 * restartPerm      The permutation applied to the columns of hVecs before restarting
 *
 * hVals            The eigenvalues of H
 *
 * flags            Array indicating the convergence of the Ritz vectors
 *
 * iev              Array of size blockSize indicating which Ritz vectors are
 *                  targeted in the block
 *
 * ievSize          The length of iev
 *
 * blockNorms       The residual norms of the eigenpairs in the block
 *
 * evecs            The converged Ritz vectors. Without locking, all converged
 *                  eigenvectors are copied from V to evecs if skew projections
 *                  are required
 *
 * evals            The converged Ritz values
 *
 * resNorms         The residual norms of the converged eigenpairs
 *
 * evecsHat         K^{-1}evecs
 *
 * ldevecsHat       The leading dimension of evecsHat
 *
 * M, ldM           evecs'*evecsHat and the leading dimension of M
 *
 * numConverged     The number of converged eigenpairs
 *
 * numConvergedStored The # of converged vectors copied to evecs
 *
 * numPrevRetained  As input the number of columns of previousHVecs. As output the
 *                  number of columns added to V
 *
 * indexOfPreviousVecs The first column in the output V that has a vector from previousHVecs
 *
 * hVecsPerm        The permutation that orders the output hVals and hVecs as primme.target
 *
 * 
 *
 * Return value
 * ------------
 * int   > 0 the restart size   
 *        -2 restart_H failed
 *        -4 factorization of M failed
 *        -5 flags do not correspond to converged pairs in pseudolocking
 *       
 ******************************************************************************/

STATIC int restart_soft_locking_Sprimme(int *restartSize, SCALAR *V, SCALAR *W,
      SCALAR *BV, PRIMME_INT nLocal, int basisSize, PRIMME_INT ldV,
      HSCALAR *hVecs, int ldhVecs, int *restartPerm, HEVAL *hVals, int *flags,
      int *iev, int *ievSize, HREAL *blockNorms, SCALAR *evecs,
      PRIMME_INT ldevecs, SCALAR *Bevecs, PRIMME_INT ldBevecs, HEVAL *evals,
      HREAL *resNorms, int *numConverged, int numPrevRetained,
      int *indexOfPreviousVecs, int *hVecsPerm, HSCALAR *H, int ldH,
      HSCALAR *VtBV, int ldVtBV, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, j, k;               /* loop indices */


   /* -------------------------------------------------------------------------- */ 
   /* Check if any of the previous flagged converged eigenvalues seems           */
   /* to have become unconverged by checking hVals[i]-evals[i] < tol.            */
   /* If it fails, flag it UNCONVERGED and let it be targeted again. This avoids */  
   /* early converged but unwanted evs preventing wanted from being targeted.    */
   /* Update maxConvTol for the remaining converged pairs.                       */
   /* Update the number of converged values also.                                */
   /* -------------------------------------------------------------------------- */

   *numConverged = 0;
   for (i=0; i<primme->numEvals; i++) {
      if (flags[i] != UNCONVERGED &&
            EVAL_ABS(hVals[i] - evals[i]) > resNorms[i]) {
         flags[i] = UNCONVERGED;
      } else if (flags[i] != UNCONVERGED && i < primme->numEvals) {
         if (flags[i] == CONVERGED) {
            if (*numConverged == 0) primme->stats.maxConvTol = 0.0;
            primme->stats.maxConvTol =
               max(primme->stats.maxConvTol, resNorms[i]);
         }
         (*numConverged)++;
      }
   }

   /* -------------------------------------------------------------- */
   /* Restart V and W by replacing it with the current Ritz vectors. */
   /* Result of restartPerm (modified part indicated with ---)       */
   /*                                                                */
   /*      converged | non-conv | prevRitzVecs | X & R               */
   /* V: [-----------|----------|--------------|- X ---|    )        */
   /* W: [-----------|----------|--------------|- R ---|    )        */
   /* BV:[-----------|----------|--------------|- BX --|    )        */
   /*                           ^ indexOfPreviousVecs                */
   /*    [-----------) numConverged            ^ restartSize         */
   /*                           [--------------) numPrevRetained     */
   /*                                  ievSize [-------)             */
   /*                                                                */
   /* X & R have the eigenvectors and residual vectors of the        */
   /* first ievSize candidates pairs to be targeted after restart.   */
   /* Their computation is performed more efficiently here at the    */
   /* same time as V, W and BV are updated.                          */
   /* -------------------------------------------------------------- */

   *indexOfPreviousVecs = *restartSize;

   *restartSize += numPrevRetained;

   *ievSize = max(0, min(min(min(min(min(*ievSize, primme->maxBlockSize),
                                       primme->numEvals - *numConverged + 1),
                                   primme->maxBasisSize - *restartSize),
                               basisSize - *numConverged),
                           primme->minRestartSize - *numConverged));

   /* Generate restartPerm */

   for (i=j=k=0; i<basisSize; i++) {
      if (k >= *numConverged || flags[i] == UNCONVERGED) {
         restartPerm[*numConverged + j++] = i;
      }
      else {
         restartPerm[k++] = i;
      }
   }
   assert(j+k == basisSize);
 
   /* Permute hVals and hVecs */

   CHKERR(KIND(permute_vecs_RHprimme, permute_vecs_SHprimme)(
         hVals, 1, basisSize, 1, restartPerm, ctx));
   CHKERR(permute_vecs_SHprimme(hVecs, basisSize, basisSize, ldhVecs,
                                restartPerm, ctx));

   /* -------------------------------------------------------------- */
   /* Restart V and W by replacing it with the current Ritz vectors. */
   /* Compute X, R, blockNorms for the next values in the block.     */
   /* -------------------------------------------------------------- */

   SCALAR *X = &V[*restartSize*ldV];
   SCALAR *R = &W[*restartSize*ldV];
   SCALAR *BX = BV ? &BV[*restartSize*ldV] : NULL;

   CHKERR(Num_aux_update_VWXR_Sprimme(V, W, BV, nLocal, basisSize, ldV,
            hVecs, *restartSize, ldhVecs, hVals,
            V, 0, *restartSize, ldV,
            X, *numConverged, *numConverged+*ievSize, ldV,
            evecs, primme->numOrthoConst, 0, *numConverged, ldevecs,
            W, 0, *restartSize, ldV,
            R, *numConverged, *numConverged+*ievSize, ldV, blockNorms,
            BV, 0, *restartSize, ldV,
            BX, *numConverged, *numConverged+*ievSize, ldV,
            Bevecs, primme->numOrthoConst, 0, *numConverged, ldBevecs,
            NULL, 0, 0,
            VtBV, *restartSize, ldVtBV,
            H, *restartSize, ldH,
            ctx));

   /* ----------------------------------------------------------------- */
   /* Generate the permutation hVecsPerm that undoes restartPerm        */
   /* ----------------------------------------------------------------- */

   for (i=0; i<basisSize; i++)
      hVecsPerm[restartPerm[i]] = i;

   /* Compute iev */

   for (i=0; i<*ievSize; i++)
      for (j=0; j<*restartSize; j++)
         if (hVecsPerm[j] == *numConverged+i)
            iev[i] = j;

   return 0;
}

/*******************************************************************************
 * Subroutine: restart_locking - This routine is only called when locking and
 *    replaces V (W) with V*c (W*c), where c is some subset of hVecs 
 *    and prevRitzVecs.
 *
 *    This subroutine locks converged Ritz pairs. The converged Ritz vectors
 *    are copied to evecs, and the converged Ritz values are copied to evals.
 *    The evals array is maintained in sorted order; however the evecs array
 *    is not.  Instead, a permutation array is maintained so that the locked
 *    Ritz vectors may be sorted upon return to the user.  
 *
 *    Vectors that have remained converged are locked to the evecs array and
 *    may be replaced afterward by initial guesses if there are any remaining.
 *
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * nLocal           Number of rows of V assigned to the node
 *
 * ldV              The leading dimension of V and W
 *
 * ldhR             The leading dimension of Q and R
 *
 * hU               The eigenvectors of QtV/R
 *
 * ldhU             The leading dimension of hU
 *
 * basisSize        Size of the basis V
 *
 * numGuesses       Number of remaining initial guesses
 *
 * previousHVecs    Coefficient vectors retained from the previous iteration
 *
 * ldpreviousHVecs  The leading dimension of previousHVecs
 *
 * rwork            Real work array
 *
 * rworkSize        Must be of size 
 *                  
 * iwork            Integer work array
 *                  
 * primme           Structure containing various solver parameters
 *
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * restartSize      The number of vectors to restart with
 *
 * V                The orthonormal basis. After restart, contains Ritz vectors
 *                  plus the orthogonal components from numPrevRetained Ritz 
 *                  vectors from the penultimate step.
 *
 * W                A*V
 *
 * BV               B*V
 *
 * X                Reference to the Ritz vectors of the eigenpairs in the block
 *
 * R                Reference to the residual vectors of the eigenpairs in the block
 *
 * hVecs, ldhVecs   The eigenvectors of H and the leading dimension of H
 *
 * restartPerm      The permutation applied to the columns of hVecs before restarting
 *
 * hVals            The eigenvalues of H
 *
 * flags            Array indicating the convergence of the Ritz vectors
 *
 * iev              Array of size blockSize indicating which Ritz vectors are
 *                  targeted in the block
 *
 * ievSize          The length of iev
 *
 * blockNorms       The residual norms of the eigenpairs in the block
 *
 * evecs            The converged Ritz vectors
 *
 * numConverged     The number of converged eigenpairs
 *
 * numLocked        The number of eigenpairs to be locked
 *
 * resNorms         The residual norms of the converged eigenpairs
 *
 * evecsperm        The permutation that orders the converged pairs as primme.target
 *
 * lockedFlags      The flags of the locked pairs
 *
 * numPrevRetained  As input the number of columns of previousHVecs. As output the
 *                  number of columns added to V
 *
 * indexOfPreviousVecs The first column in the output V that has a vector from previousHVecs
 *
 * hVecsPerm        The permutation that orders the output hVals and hVecs as primme.target
 *
 * numArbitraryVecs On input, the number of leading coefficient vectors in
 *                  hVecs that do not come from solving the projected problem.
 *                  On output, the number of such vectors that are in the
 *                  restarted basis
 *
 * Return value
 * ------------
 * int   > 0 the restart size   
 *        -2 restart_H failed
 *        -4 factorization of M failed
 *        -5 flags do not correspond to converged pairs in pseudolocking
 *       
 ******************************************************************************/

STATIC int restart_locking_Sprimme(int *restartSize, SCALAR *V, SCALAR *W,
      SCALAR *BV, PRIMME_INT nLocal, int basisSize, PRIMME_INT ldV,
      HSCALAR *hVecs, int ldhVecs, int *restartPerm, HEVAL *hVals, int *flags,
      int *iev, int *ievSize, HREAL *blockNorms, SCALAR *evecs,
      PRIMME_INT ldevecs, SCALAR *Bevecs, PRIMME_INT ldBevecs, HEVAL *evals,
      int *numConverged, int *numLocked, HREAL *resNorms, int *lockedFlags,
      int *evecsperm, int numPrevRetained, int *indexOfPreviousVecs,
      int *hVecsPerm, HSCALAR *H, int ldH, HSCALAR *VtBV, int ldVtBV,
      double startTime, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, j, k;             /* Loop variables                                 */
   int numPacked;           /* The number of coefficient vectors moved to the */
                            /* end of the hVecs array.                        */
   int maxBlockSize;        /* max block size for the next iteration          */
   int sizeBlockNorms;      /* Size of the block                              */
   int failed;              /* Number of vectors to be locked that didn't pass the conv. test */
   int left;                /* Index of the first column with a vector to be locked */
   HREAL *lockedResNorms;  /* residual norms for locked vectors              */
   int numLocked0 = *numLocked; /* aux variables                         */


   /* -------------------------------------------------------------- */
   /* Rearrange hVals and hVecs so that V*hVecs and W*hVecs will     */
   /* look like this (modified part indicated with ---):             */
   /*                                                                */
   /*      restarted  | prevRitzVecs | to be locked | X & R          */
   /* V: [------------|--------------|--------------|- X -|    )     */
   /* W: [------------|--------------|--------------|- R -|    )     */
   /* BV:[------------|--------------|--------------|- BX-|    )     */
   /*                                ^ left         ^ restartSize    */
   /*                 ^ indexOfPreviousVecs                          */
   /*          [--------------) numPrevRetained                      */
   /*                                [--------------) numPacked      */
   /*                         [----) sizeBlockNorms [-----)          */
   /*                                maxBlockSize   [-------)        */
   /*                                                                */
   /*                           to be locked                         */
   /* evecs: [----------------|--------------)                       */
   /* Bevecs:[----------------|--------------)                       */
   /*                         ^ numLocked                            */
   /*                                                                */
   /* X & R have the eigenvectors and residual vectors of the        */
   /* first sizeBlockNorms candidates pairs to be targeted after     */
   /* restart. Their computation is performed more efficiently here  */
   /* together with the V, W                                         */
   /* -------------------------------------------------------------- */

   maxBlockSize = max(0, min(min(
         *restartSize,
         primme->maxBlockSize),
         primme->numEvals-*numConverged+1));

   sizeBlockNorms = max(0, min(
         maxBlockSize,
         primme->maxBasisSize-*restartSize-numPrevRetained
               -*numConverged+*numLocked));

   *indexOfPreviousVecs = *restartSize;

   left = *restartSize + numPrevRetained;

   for (i=k=numPacked=0; i<basisSize; i++) {
      if (flags[i] != UNCONVERGED
            /* Otherwise don't check more than numEvals-numLocked             */
                  && numPacked < primme->numEvals-*numLocked
            /* Check only the first pairs, except if finding closest_leq/geq  */
            /* because with refined extraction the pairs may not be ordered   */
            /* by this criterion.                                             */
                  && (i < primme->numEvals-*numLocked
                     || primme->target == primme_closest_geq
                     || primme->target == primme_closest_leq)) {
         restartPerm[left + numPacked++] = i;
      }
      else if (k < left) {
         restartPerm[k++] = i;
      }
      else {
         restartPerm[min(*numConverged, primme->numEvals)-*numLocked + k++] = i;
      }
   }

   *restartSize = left + numPacked;

   KIND(permute_vecs_RHprimme, permute_vecs_SHprimme)
   (hVals, 1, basisSize, 1, restartPerm, ctx);
   permute_vecs_SHprimme(hVecs, basisSize, basisSize, ldhVecs, restartPerm,
                         ctx);

   /* -------------------------------------------------------------- */
   /* Restart V and W by replacing it V*hVecs and W*hVecs, copy      */
   /* vectors to be locked into evecs and compute X and R.           */
   /* -------------------------------------------------------------- */

   lockedResNorms = &resNorms[*numLocked];
   SCALAR *X = &V[*restartSize*ldV];
   SCALAR *R = &W[*restartSize*ldV];
   SCALAR *BX = BV ? &BV[*restartSize*ldV] : NULL;
   assert(*restartSize <= primme->maxBasisSize);
   CHKERR(Num_aux_update_VWXR_Sprimme(V, W, BV, nLocal, basisSize, ldV,
            hVecs, *restartSize, ldhVecs, hVals,
            V, 0, *restartSize, ldV,
            X, 0, sizeBlockNorms, ldV,
            evecs, *numLocked+primme->numOrthoConst, left,
                left+numPacked, ldevecs,
            W, 0, *restartSize, ldV,
            R, 0, sizeBlockNorms, ldV, blockNorms,
            BV, 0, *restartSize, ldV,
            BX, 0, sizeBlockNorms, ldV,
            Bevecs, *numLocked+primme->numOrthoConst, left,
                left+numPacked, ldBevecs,
            lockedResNorms, left, *restartSize,
            VtBV, *restartSize, ldVtBV,
            H, *restartSize, ldH,
            ctx));
 
   /* -------------------------------------------------------------- */
   /* Recompute flags for the vectors to be locked.                  */
   /* NOTE: the eigenvals flagged as practically converged will keep */
   /*       it if their residual norm is still less than             */
   /*       sqrt(numLocked)*tol                                      */
   /* -------------------------------------------------------------- */

   permute_vecs_iprimme(flags, basisSize, restartPerm, ctx);
   CHKERR(check_convergence_Sprimme(&V[ldV * left], ldV, 1 /* X given */, NULL,
         0, 0 /* R not given */, NULL, *numLocked, 0, NULL, 0, NULL, 0, left,
         left + numPacked, flags, lockedResNorms, hVals, NULL, 0, ctx));

   /* -------------------------------------------------------------- */
   /* Copy the values for the converged values into evals            */
   /* -------------------------------------------------------------- */

   for (i=left, j=0; i < left+numPacked; i++) {
      if (flags[i] != UNCONVERGED && *numLocked+j < primme->numEvals) {
         evals[*numLocked+j++] = hVals[i];
      }
      else {
         flags[i] = UNCONVERGED;
      }
   }

   /* -------------------------------------------------------------- */
   /* Merge back the pairs that failed to be locked into the         */
   /* restarted basis.                                               */
   /* -------------------------------------------------------------- */

   /* Generate the permutation that separates the pairs that failed  */
   /* to be locked to the ones that passed.                          */

   int *ifailed;            /* Indices of vectors to be locked vectors that failed */
   CHKERR(Num_malloc_iprimme(numPacked, &ifailed, ctx));
   for (i=left, failed=0; i < left+numPacked; i++)
      if (flags[i] == UNCONVERGED) ifailed[failed++] = i-left;
   for (i=left, j=0; i < left+numPacked; i++)
      if (flags[i] != UNCONVERGED) ifailed[failed+j++] = i-left;

   if (1 /* Put zero to disable the new feature */) {
      /* New feature: the pairs that failed to be locked are         */
      /* rearranged with the rest of the restarted vectors and they  */
      /* are considered to be included in the block.                 */

      /* Generate hVecsPerm merging candidate vectors and vectors to */
      /* be locked that failed to pass the convergence test with the */
      /* restarted pairs. The order is determined by the original    */
      /* position of the pairs.                                      */

      maxBlockSize = max(0, min(
            maxBlockSize,
            primme->maxBasisSize-*restartSize-numPrevRetained
                  -*numConverged+*numLocked));

      HREAL *blockNorms0;
      CHKERR(Num_malloc_RHprimme(sizeBlockNorms, &blockNorms0, ctx));
      for (i=0; i<sizeBlockNorms; i++) blockNorms0[i] = blockNorms[i];
      for (i=j=k=0; i<*indexOfPreviousVecs || j<failed; k++) {

         /* If there is no more failed vectors or the candidate      */
         /* vector was before the failed, take the candidate         */

         if (i < *indexOfPreviousVecs && (j >= failed 
                  || restartPerm[i] < restartPerm[left+ifailed[j]])) {

            if (k < maxBlockSize && i < sizeBlockNorms)
               blockNorms[k] = blockNorms0[i];
            hVecsPerm[k] = i++;

         }

         /* Otherwise take the failed */

         else {
            if (k < maxBlockSize)
               blockNorms[k] = resNorms[numLocked0+ifailed[j]];
            hVecsPerm[k] = left + j++;
         }
      }
      CHKERR(Num_free_RHprimme(blockNorms0, ctx));

      /* Generate the rest of the permutation of hVecsPerm           */

      for (i=0; i<numPrevRetained; i++) hVecsPerm[k++] = i+*indexOfPreviousVecs;
      assert(k == *indexOfPreviousVecs + numPrevRetained + failed);
      for (; k < basisSize; k++) hVecsPerm[k] = -1;

      /* -------------------------------------------------------------- */
      /* Update X and R with the failed vectors and place the failed    */
      /* vectors after the candidate vectors (modified part indicated   */
      /* with ---):                                                     */
      /*                                                                */
      /*      restarted  | prevRitzVecs | failed locked | X & R         */
      /* V: [            |              |---------------|- X ---|    )  */
      /* W: [            |              |---------------|- R ---|    )  */
      /* BV:[            |              |---------------|- BX --|    )  */
      /*                                ^ left                          */
      /*                         failed [---------------)               */
      /*                                 maxBlockSize   [---------)     */
      /*                                                                */
      /* -------------------------------------------------------------- */

      CHKERR(compute_residual_columns(nLocal, &hVals[left], &V[left * ldV],
            BV ? &BV[left * ldV] : NULL, failed, ifailed, ldV, &W[left * ldV],
            ldV, X, BX, sizeBlockNorms, ldV, 0, R, ldV,
            &V[(left + failed) * ldV], BV ? &BV[(left + failed) * ldV] : NULL,
            maxBlockSize, hVecsPerm, ldV, &W[(left + failed) * ldV], ldV, ctx));
   }
   else {
      /* The failed pairs are not rearranged with the rest of           */
      /* non-converged pairs and they may not be in the block in the    */
      /* next iteration. This was the behaviour of locking in previous  */
      /* versions than 2.0.                                             */

      for (i=0; i<basisSize; i++) hVecsPerm[i] = i;

      /* Pack together the failed vectors and place them after          */
      /* candidates.                                                    */

      Num_compact_vecs_Sprimme(&V[left*ldV], nLocal, failed, ldV,
            ifailed, &V[left*ldV], ldV, 0, ctx);
      Num_compact_vecs_Sprimme(&W[left*ldV], nLocal, failed, ldV,
            ifailed, &W[left*ldV], ldV, 0, ctx);

      /* Copy X and R after the failed vectors */

      CHKERR(Num_copy_matrix_Sprimme(X, nLocal, sizeBlockNorms, ldV,
            &V[(left + failed) * ldV], ldV, ctx));
      CHKERR(Num_copy_matrix_Sprimme(R, nLocal, sizeBlockNorms, ldV,
            &W[(left + failed) * ldV], ldV, ctx));
      CHKERR(Num_copy_matrix_Sprimme(BX, nLocal, sizeBlockNorms, ldV,
            &BV[(left + failed) * ldV], ldV, ctx));
   }

   /* Modify hVals, hVecs and restartPerm to add the pairs failed to */
   /* be locked just after the candidates                            */

   Num_compact_vecs_SHprimme(&hVecs[left * ldhVecs], basisSize, failed, ldhVecs,
         ifailed, &hVecs[left * ldhVecs], ldhVecs, 0, ctx);
   KIND(Num_compact_vecs_RHprimme, Num_compact_vecs_SHprimme)
   (&hVals[left], 1, failed, 1, ifailed, &hVals[left], 1, 0, ctx);
   CHKERR(permute_vecs_iprimme(&restartPerm[left], numPacked, ifailed, ctx));

   /* Update VtBV */

   if (VtBV) {
      int *iV;
      CHKERR(Num_malloc_iprimme(left + numPacked, &iV, ctx));
      for (i = 0; i < numPacked - failed; i++)
         iV[i] = ifailed[failed + i] + left;
      for (i = 0; i < left; i++) iV[i + numPacked - failed] = i;
      for (i = 0; i < failed; i++)
         iV[i + left + numPacked - failed] = ifailed[i] + left;
      int nLocked = primme->numOrthoConst + *numLocked;
      int nVtBV = nLocked + left + numPacked;
      HSCALAR *rwork;
      CHKERR(Num_malloc_SHprimme(nVtBV * (left + numPacked), &rwork, ctx));
      CHKERR(Num_copy_matrix_columns_SHprimme(&VtBV[nLocked * ldVtBV], nVtBV,
            iV, left + numPacked, ldVtBV, rwork, NULL, nVtBV, ctx));
      CHKERR(Num_copy_matrix_SHprimme(rwork, nLocked, left + numPacked, nVtBV,
            &VtBV[nLocked * ldVtBV], ldVtBV, ctx));
      CHKERR(Num_copy_matrix_rows_SHprimme(rwork + nLocked, iV,
            left + numPacked, left + numPacked, nVtBV,
            &VtBV[ldVtBV * nLocked + nLocked], NULL, ldVtBV, ctx));
      CHKERR(Num_free_SHprimme(rwork, ctx));
      CHKERR(Num_free_iprimme(iV, ctx));
   }

   /* Update H */

   if (primme->orth == primme_orth_explicit_I) {
      CHKERR(Num_copy_matrix_columns_SHprimme(&H[left * ldH], left + numPacked,
            ifailed, failed, ldH, &H[left * ldH], NULL, ldH, ctx));
      CHKERR(Num_copy_matrix_rows_SHprimme(H + left, ifailed, failed,
            left + failed, ldH, &H[left], NULL, ldH, ctx));
   }

   CHKERR(Num_free_iprimme(ifailed, ctx));

   X = &V[(left+failed)*ldV];
   R = &W[(left+failed)*ldV];

   /* Pack those Ritz vectors that are actually converged from */
   /* those that have been copied into evecs. Also insertSort  */
   /* the converged Ritz value within the evals array.         */ 

   for (i=left; i < left+numPacked; i++) {
       if (flags[i] != UNCONVERGED && *numLocked < primme->numEvals) {
         HREAL resNorm = resNorms[*numLocked] = lockedResNorms[i-left];
         HEVAL eval = evals[*numLocked];
         CHKERR(Num_copy_matrix_Sprimme(
               &evecs[(numLocked0 + i - left + primme->numOrthoConst) *
               ldevecs],
               nLocal, 1, ldevecs,
               &evecs[(*numLocked + primme->numOrthoConst) * ldevecs], ldevecs,
               ctx));
         if (Bevecs) {
            CHKERR(Num_copy_matrix_Sprimme(
                  &Bevecs[(numLocked0 + i - left + primme->numOrthoConst) *
                  ldBevecs],
                  nLocal, 1, ldBevecs,
                  &Bevecs[(*numLocked + primme->numOrthoConst) * ldBevecs],
                  ldBevecs, ctx));
         }

         (*numLocked)++;

         /* Report a pair was hard locked */
         /* NOTE: do this before sorting evals */
         if (primme->monitorFun) {
            lockedFlags[*numLocked-1] = flags[i];
            CHKERR(monitorFun_Sprimme(NULL, 0, NULL, NULL, 0, NULL, 0, evals,
                  *numLocked, lockedFlags, resNorms, -1, 0.0, NULL, 0.0,
                  primme_event_locked, startTime, ctx));
         }

         CHKERR(insertionSort_Sprimme(eval, evals, resNorm, resNorms, flags[i],
               lockedFlags, evecsperm, *numLocked - 1, 0, primme));

         /* Update maxConvTol if it wasn't practically converged */
         if (flags[i] == CONVERGED) {
            primme->stats.maxConvTol = max(primme->stats.maxConvTol, resNorm);
         }
      }
   }

   *restartSize = left + failed;
   *ievSize = min(maxBlockSize, sizeBlockNorms+failed);
   *numConverged = *numLocked;
   for (i=0; i<*ievSize; i++) iev[i] = i;

   /* ----------------------------------------------------------------- */
   /* There is no converged pair in the restarted basis                 */
   /* ----------------------------------------------------------------- */
   
   for (i=0; i<basisSize; i++) flags[i] = UNCONVERGED;

   return 0;
}


/******************************************************************************
 * Function Num_aux_update_VWXR - This subroutine performs the next operations:
 *
 *    X0 = V*h(nX0b+1:nX0e), X1 = V*h(nX1b+1:nX1e), evecs(evecsSize:) = V*h(nX2b+1:nX2e)
 *    Wo = A*X0(nWob+1-nX0b:nWoe-nX0b),
 *    R = Wo(nRb+1-nWob:nRe-nWob) - X0(nRb+1-nX0b:nRe-nX0b)*diag(hVals(nRb+1:nRe)),
 *    Rnorms = norms(R),
 *    rnorms = norms(Wo(nrb+1-nWob:nre-nWob) - X0(nrb+1-nX0b:nre-nX0b)*diag(hVals(nrb+1:nre))),
 *
 * NOTE: if Rnorms and rnorms are requested, nRb-nRe+nrb-nre < mV
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V, W        input basis
 * mV,nV,ldV   number of rows and columns and leading dimension of V and W
 * h           input rotation matrix
 * nh          Number of columns of h
 * ldh         The leading dimension of h
 * hVals       Array of values
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * X0          Output matrix V*h(nX0b:nX0e-1) (optional)
 * nX0b, nX0e  Range of columns of h
 * X1          Output matrix V*h(nX1b:nX1e-1) (optional)
 * nX1b, nX1e  Range of columns of h
 * evecs       Output matrix V*h(nX2b:nX2e-1) (optional)
 * evecsSize   First column where to copy V*h(nX2b:nX2e-1) (optional)
 * nX2b, nX2e  Range of columns of h
 * Wo          Output matrix W*h(nWob:nWoe-1) (optional)
 * nWob, nWoe  Range of columns of h
 * R           Output matrix (optional)
 * nRb, nRe    Range of columns of h and hVals
 * Rnorms      Output array with the norms of R (optional)
 * rnorms      Output array with the extra residual vector norms (optional)
 * nrb, nre    Columns of residual vector to compute the norm
 * 
 * NOTE: n*e, n*b are zero-base indices of ranges where the first value is
 *       included and the last isn't.
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_aux_update_VWXR_Sprimme(SCALAR *V, SCALAR *W, SCALAR *BV,
   PRIMME_INT mV, int nV, PRIMME_INT ldV,
   HSCALAR *h, int nh, int ldh, HEVAL *hVals,
   SCALAR *X0, int nX0b, int nX0e, PRIMME_INT ldX0,
   SCALAR *X1, int nX1b, int nX1e, PRIMME_INT ldX1,
   SCALAR *evecs, int evecsSize, int nX2b, int nX2e, PRIMME_INT ldevecs,
   SCALAR *Wo, int nWob, int nWoe, PRIMME_INT ldWo,
   SCALAR *R, int nRb, int nRe, PRIMME_INT ldR, HREAL *Rnorms,
   SCALAR *BX0, int nBX0b, int nBX0e, PRIMME_INT ldBX0,
   SCALAR *BX1, int nBX1b, int nBX1e, PRIMME_INT ldBX1,
   SCALAR *Bevecs, int BevecsSize, int nBX2b, int nBX2e, PRIMME_INT ldBevecs,
   HREAL *rnorms, int nrb, int nre,
   HSCALAR *VtBV, int nVtBV, int ldVtBV,
   HSCALAR *H, int nH, int ldH, 
   primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i;         /* Loop variables */

   CHKERR(Num_update_VWXR_Sprimme(
            V, W, BV, mV, nV, ldV, h, nh, ldh, hVals,
            X0, nX0b, nX0e, ldX0,
            X1, nX1b, nX1e, ldX1,
            evecs?&evecs[ldevecs*evecsSize]:NULL, nX2b, nX2e, ldevecs,
            Wo, nWob, nWoe, ldWo,
            R, nRb, nRe, ldR, Rnorms,
            BX0, nBX0b, nBX0e, ldBX0,
            BX1, nBX1b, nBX1e, ldBX1,
            Bevecs?&Bevecs[ldBevecs*BevecsSize]:NULL, nBX2b, nBX2e, ldBevecs,
            rnorms, nrb, nre,
            VtBV?&VtBV[evecsSize*ldVtBV+evecsSize]:NULL, nVtBV, ldVtBV,
            KIND(primme->orth == primme_orth_explicit_I,1)?H:NULL, nH, ldH,
            NULL, 0, 0,
            ctx));

   /* Don't trust residual norm smaller than the error in the residual norm */

   for (i = 0; i < nre - nrb; i++) {
      rnorms[i] = max(rnorms[i], primme->stats.estimateResidualError);
   }

   /* Update VtBV */

   if (VtBV) {
      HSCALAR *work;
      assert(nVtBV == nX0e - nX0b);
      CHKERR(Num_malloc_SHprimme((nX0e - nX0b) * evecsSize, &work, ctx));
      if (ctx.procID == 0) {
         Num_zero_matrix_SHprimme(work, evecsSize, nX0e - nX0b, evecsSize, ctx);
         CHKERR(Num_gemm_SHprimme("N", "N", evecsSize, nX0e - nX0b, nV, 1.0,
               &VtBV[evecsSize * ldVtBV], ldVtBV, h, ldh, 0.0, work, evecsSize,
               ctx));
      }
      CHKERR(broadcast_SHprimme(work, evecsSize * (nX0e - nX0b), ctx));
      CHKERR(Num_copy_matrix_SHprimme(work, evecsSize, nX0e - nX0b, evecsSize,
            &VtBV[evecsSize * ldVtBV], ldVtBV, ctx));
      CHKERR(Num_free_SHprimme(work, ctx));
   }

   return 0;

}


/*******************************************************************************
 * Subroutine: restart_projection - This routine updates Q, R, H and QtV to match the
 *    changes in V.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * nLocal           Number of rows of V assigned to the node
 *
 * basisSize        Size of the basis V
 *
 * numConverged     The number of converged eigenpairs without locking
 *
 * numLocked        The number of Ritz vectors that have been locked 
 *
 * numPrevRetained  The number of coefficient vectors from previousHVecs in hVecs
 *  
 * indexOfPreviousVecs The first column in hVecs that has a vector from previousHVecs
 *
 * indexOfPreviousVecsBeforeRestart The number of columns in hVecs that previousHVecs
 *                  where orthogonalized against
 *
 * rwork            Real work array
 *
 * rworkSize        Size of rwork
 *                  
 * iwork            Integer work array
 *
 * primme           Structure containing various solver parameters
 *
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * restartSize      The number of vectors to restart with. If negative,
 *                  use dynamic thick restart.
 *
 * V                The orthonormal basis. After restart, contains Ritz vectors
 *                  plus the orthogonal components from numPrevRetained Ritz 
 *                  vectors from the penultimate step.
 *
 * ldV              The leading dimension of V
 *
 * W                A*V
 *
 * ldW              The leading dimension of W
 *
 * H                The projection V'*A*V
 *
 * ldH              The leading dimension of H
 *
 * VtBV             The projection V'*B*V
 *
 * ldVtBV           The leading dimension of VtBV
 *
 * Q, R             The factors of the QR decomposition of (A - targetShift*B)*V
 *
 * ldQ, ldR         The leading dimension of Q and R
 *
 * fQtQ             The Cholesky factor of QtQ
 *
 * ldfQtQ           The leading dimension of fQtQ
 * 
 * QtV              = Q'*V
 *
 * ldQtV            The leading dimension of QtV
 *
 * hU               The left coefficient vectors or the eigenvectors of QtV/R
 *
 * ldhU             The leading dimension of the input hU
 *
 * newldhU          The leading dimension of the output hU
 *
 * hVecs            The coefficient vectors
 *
 * ldhVecs          The leading dimension of the input hVecs
 *
 * newldhVecs       The leading dimension of the output hVecs
 *
 * hVals            The eigenvalues of H
 *
 * hSVals           The singular values of R
 *
 * hVecsPerm        The permutation that orders hVals and hVecs as primme.target
 *
 * restartPerm      The permutation applied to the columns of hVecs before restarting
 *
 * flags            Array indicating the convergence of the Ritz vectors
 *
 * evecs            The converged Ritz vectors. Without locking, all converged
 *                  eigenvectors are copied from V to evecs if skew projections
 *                  are required.
 *
 * ldevecs          The leading dimension of evecs
 *
 * evecsHat         K^{-1}evecs
 *
 * ldevecsHat       The leading dimension of evecsHat
 *
 * M                evecs'*evecsHat
 *
 * ldM              The leading dimension of M
 *
 * Mfact            The factorization of M
 *
 * ldMfact          The leading dimension of Mfact
 *
 * ipivot           The pivot array of the Mfact factorization
 *
 * targetShiftIndex The target shift used in (A - targetShift*B) = Q*R
 *
 * numArbitraryVecs On input, the number of coefficients vectors that do
 *                  not correspond to solutions of the projected problem.
 *                  They appear in the first numArbitraryVecs positions of hVecs.
 *                  On output, the number of such vectors that are in the
 *                  restarted basis.
 *
 * Return value
 * ------------
 * int   > 0 the restart size   
 *        -2 restart_H failed
 *        -4 factorization of M failed
 *        -5 flags do not correspond to converged pairs in pseudolocking
 *       
 ******************************************************************************/

STATIC int restart_projection_Sprimme(SCALAR *V, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, SCALAR *BV, PRIMME_INT ldBV, HSCALAR *H, int ldH,
      HSCALAR *VtBV, int ldVtBV, SCALAR *Q, PRIMME_INT ldQ, PRIMME_INT nLocal,
      HSCALAR *R, int ldR, HSCALAR *QtV, int ldQtV, HSCALAR *QtQ, int ldQtQ,
      HSCALAR *fQtQ, int ldfQtQ, HSCALAR *hU, int ldhU, int newldhU,
      int indexOfPreviousVecsBeforeRestart, HSCALAR *hVecs, int ldhVecs,
      int newldhVecs, HEVAL *hVals, HREAL *hSVals, int *restartPerm,
      int *hVecsPerm, int restartSize, int basisSize, int numPrevRetained,
      int indexOfPreviousVecs, SCALAR *evecs, int *evecsSize,
      PRIMME_INT ldevecs, SCALAR *Bevecs, PRIMME_INT ldBevecs, SCALAR *evecsHat,
      PRIMME_INT ldevecsHat, HSCALAR *M, int ldM, HSCALAR *Mfact, int ldMfact,
      int *ipivot, int *targetShiftIndex, int numConverged,
      int *numArbitraryVecs, HSCALAR *hVecsRot, int ldhVecsRot,
      primme_context ctx) {

   primme_params *primme = ctx.primme;

   /* -------------------------------------------------------- */
   /* Restart projected problem matrices H and R               */
   /* -------------------------------------------------------- */

   switch (primme->projectionParams.projection) {
   case primme_proj_RR:
      CHKERR(restart_RR(H, ldH, VtBV, ldVtBV, hVecs, ldhVecs, newldhVecs, hVals,
            restartSize, basisSize, numConverged, numPrevRetained,
            indexOfPreviousVecs, hVecsPerm, targetShiftIndex, ctx));
      break;

   case primme_proj_harmonic:
      CHKERR(restart_harmonic(V, ldV, W, ldW, BV, ldBV, H, ldH, Q, ldQ, nLocal,
            R, ldR, QtV, ldQtV, QtQ, ldQtQ, fQtQ, ldfQtQ, VtBV, ldVtBV, hU,
            ldhU, newldhU, hVecs, ldhVecs, newldhVecs, hVals, hSVals,
            restartPerm, hVecsPerm, restartSize, basisSize, numPrevRetained,
            indexOfPreviousVecs, targetShiftIndex, numConverged,
            numArbitraryVecs, hVecsRot, ldhVecsRot, ctx));
      break;

   case primme_proj_refined:
      CHKERR(restart_refined(V, ldV, W, ldW, BV, ldBV, H, ldH, Q, ldQ, nLocal,
            R, ldR, QtQ, ldQtQ, fQtQ, ldfQtQ, VtBV, ldVtBV, hU, ldhU, newldhU,
            indexOfPreviousVecsBeforeRestart, hVecs, ldhVecs, newldhVecs, hVals,
            hSVals, restartPerm, hVecsPerm, restartSize, basisSize,
            numPrevRetained, indexOfPreviousVecs, targetShiftIndex,
            numConverged, numArbitraryVecs, hVecsRot, ldhVecsRot, ctx));
      break;

   default:
      assert(0);
   }

   if (evecsHat) {
   /* --------------------------------------------------------------------- */
   /* If the user requires (I-QQ') projectors in JDQMR without locking,     */
   /* the converged eigenvectors are copied temporarily to evecs. There     */
   /* they stay locked  for use in (I-QQ') and (I-K^{-1}Q () Q') projectors.*/
   /* The Ritz vectors remain in the basis, and they will overwrite evecs   */
   /* the end.                                                              */
   /* We recommend against this type of usage. It's better to use locking.  */
   /* --------------------------------------------------------------------- */

   /* Andreas NOTE: is done inefficiently for the moment. We should only */
   /* add the recently converged. But we need to differentiate them      */
   /* from flags...                                                      */
   /* Eloy NOTE: that will not be enough, the Ritz vectors change from   */
   /* one restart to the next even if they weren't targeted. Probably    */
   /* the best we can do is to always copy all converged vectors.        */

      if (!primme->locking) *evecsSize = 0;

      int numRecentlyConverged = numConverged - *evecsSize;

      /* Compute K^{-1}x for all newly locked eigenvectors */

      double *shiftsForPreconditioner;
      if (numConverged <= primme->numTargetShifts) {
         shiftsForPreconditioner = &primme->targetShifts[*evecsSize];
      } else if (primme->numTargetShifts > 0) {
         CHKERR(
               Num_malloc_dprimme(numConverged, &shiftsForPreconditioner, ctx));
         int i;
         for (i = 0; i < numRecentlyConverged; i++) {
            shiftsForPreconditioner[i] = primme->targetShifts[min(
                  i + *evecsSize, primme->numTargetShifts - 1)];
         }
      } else {
         shiftsForPreconditioner = NULL;
      }

      primme->ShiftsForPreconditioner = shiftsForPreconditioner;

      CHKERR(applyPreconditioner_Sprimme(
               &Bevecs[ldBevecs*(*evecsSize+primme->numOrthoConst)], nLocal,
               ldBevecs,
               &evecsHat[ldevecsHat*(*evecsSize+primme->numOrthoConst)],
               ldevecsHat, numRecentlyConverged, ctx));

      if (numConverged > primme->numTargetShifts) {
         CHKERR(Num_free_dprimme(shiftsForPreconditioner, ctx));
      }

      /* Update the projection evecs'*evecsHat now that evecs and evecsHat   */
      /* have been expanded by numRecentlyConverged columns.  Required       */
      /* workspace is numLocked*numEvals.  The most ever needed would be     */
      /* maxBasisSize*numEvals.                                              */

      CHKERR(update_XKinvBX_Sprimme(evecs, ldevecs, evecsHat, ldevecsHat, M,
            ldM, *evecsSize, numConverged - *evecsSize, Mfact, ldMfact, ipivot,
            ctx));

      *evecsSize = numConverged;
   }

   return 0;
}

/*******************************************************************************
 * Function restart_RR - This routine is used to recompute H = V'*A*V once V 
 *   has been restarted.  If no coefficient vectors from the previous iteration
 *   have been retained, then the restarted H will be diagonal and the 
 *   new eigenvectors (coefficient vectors) of H will be the standard basis
 *   vectors.  If previous coefficient vectors have been retained, then H will 
 *   contain the numPrevRetained x numPrevRetained submatrix 
 *   previousHVecs'*H*previousHvecs and the rest of the elements of H will be
 *   along the diagonal.
 *   
 *
 * INPUT PARAMETERS
 * ----------------
 * restartSize   Number of vectors the basis was restarted with
 * 
 * basisSize     Maximum size of the basis V
 *
 * numLocked     The number of Ritz vectors that have been locked 
 *
 * numConverged  The number of converged eigenpairs
 *
 * numPrevRetained The number of vectors retained from the previous iteration
 *
 * indexOfPreviousVecs  The index within hVecs where the previous vectors were
 *                      inserted.  Its also where the overlap matrix
 *                      previousHVecs'*H*previousHvecs will be inserted within
 *                      the restarted H.
 *
 * hVecsPerm     The permutation that orders hVals and hVecs as primme.target
 *
 * iwork         Integer work array
 *
 * rwork         Work array.  Necessary only when coefficient vectors from the
 *               previous iteration are retained.
 *
 * rworkSize     Size of rwork
 *
 *
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * H      Will contain H = V'*A*V given the restarted V.  H will be
 *        diagonal since V will contain only Ritz vectors, unless previous
 *        vectors are retained.  In that case, it will be diagonal except for
 *        a numPrevRetained x numPrevRetained submatrix inserted at
 *        H(numPrevRetained, numPrevRetained)
 *
 * ldH    The leading dimension of H
 *
 * VtBV   Will contain VtBV = V'*B*V given the restarted V. VtBV will be
 *        diagonal since V will contain only Ritz vectors, unless previous
 *        vectors are retained.  In that case, it will be diagonal except for
 *        a numPrevRetained x numPrevRetained submatrix inserted at
 *        VtBV(numPrevRetained, numPrevRetained)
 *
 * ldVtBV The leading dimension of VtBV
 *
 * hVecs  If the new H is diagonal, then it will contain the standard basis
 *        vectors.  If previous coefficient vectors are retained, then 
 *        restartSize - numPrevRetained of the vectors will be standard basis
 *        vectors.  The remaining numPrevRetained vectors will contain
 *        numPrevRetained non-zero elements corresponding to the 
 *        numPrevRetined x numPrevRetained submatrix.
 * 
 * ldhVecs  The leading dimension of the input hVecs
 *
 * newldhVecs  The leading dimension of the output hVecs
 *
 * hVals  The eigenvalues of the restarted H
 * 
 * targetShiftIndex The target shift used in (A - targetShift*B) = Q*R
 * 
 * Return value
 * ------------
 * Error code: 0 upon success
 *            -1 eigenvalues of submatrix could not be computed
 *
 ******************************************************************************/

STATIC int restart_RR(HSCALAR *H, int ldH, HSCALAR *VtBV, int ldVtBV,
      HSCALAR *hVecs, int ldhVecs, int newldhVecs, HEVAL *hVals, int restartSize,
      int basisSize, int numConverged, int numPrevRetained,
      int indexOfPreviousVecs, int *hVecsPerm, int *targetShiftIndex,
      primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, j;          /* Loop variables                                       */
   int orderedIndexOfPreviousVecs;  /* index of prev. vecs after applying hVecsPerm */
   double aNorm = primme?max(primme->aNorm, primme->stats.estimateLargestSVal):0.0;

#ifdef USE_HERMITIAN
   if (primme->orth == primme_orth_implicit_I) {

      // If coefficient vectors from the previous iteration were retained, then
      // insert the computed overlap matrix into the restarted H

      CHKERR(compute_submatrix_SHprimme(&hVecs[ldhVecs * indexOfPreviousVecs],
            numPrevRetained, ldhVecs, H, basisSize, ldH, 1 /* Hermitian */,
            &H[ldH * indexOfPreviousVecs + indexOfPreviousVecs], ldH, ctx));

      /* ----------------------------------------------------------------- */
      /* Y = V*hVecs([0:indexOfPreviousVecs-1 \                            */
      /*                 indexOfPreviousVecs+numPrevRetained:restartSize]) */
      /* are the new Ritz vectors so the new H and VtBV will be diag(hVals)*/
      /* and I respectively. Note it is only stored the upper triangular   */
      /* part of H and VtBV.                                               */
      /* ----------------------------------------------------------------- */

      CHKERR(Num_zero_matrix_SHprimme(
            H, restartSize, indexOfPreviousVecs, ldH, ctx));
      CHKERR(Num_zero_matrix_SHprimme(&H[ldH * indexOfPreviousVecs],
            indexOfPreviousVecs, numPrevRetained, ldH, ctx));
      CHKERR(Num_zero_matrix_SHprimme(
            &H[ldH * (indexOfPreviousVecs + numPrevRetained)], restartSize,
            restartSize - indexOfPreviousVecs - numPrevRetained, ldH, ctx));

      for (j = 0; j < indexOfPreviousVecs; j++) {
         H[ldH * j + j] = hVals[j];
      }
      for (j = indexOfPreviousVecs + numPrevRetained; j < restartSize; j++) {
         H[ldH * j + j] = hVals[j];
      }
   }
#else
   (void)ldhVecs;
   (void)basisSize;
#endif /* USE_HERMITIAN */

   /* ---------------------------------------------------------------------- */
   /* Solve the whole matrix when the targetShift has to change              */ 
   /* ---------------------------------------------------------------------- */

   int nLocked = primme->numOrthoConst + (primme->locking ? numConverged : 0);
   double eps_matrix;
   CHKERR(machineEpsMatrix_Sprimme(&eps_matrix, ctx));
   if (targetShiftIndex && primme->targetShifts &&
         (*targetShiftIndex < 0 ||
               fabs(primme->targetShifts[*targetShiftIndex] -
                     primme->targetShifts[min(primme->numTargetShifts - 1,
                           numConverged)])) > eps_matrix * aNorm) {

      *targetShiftIndex = min(primme->numTargetShifts-1, numConverged);

      CHKERR(solve_H_SHprimme(H, restartSize, ldH,
            VtBV ? &VtBV[nLocked * ldVtBV + nLocked] : NULL, ldVtBV, NULL, 0,
            NULL, 0, NULL, 0, NULL, 0, hVecs, newldhVecs, hVals, NULL,
            numConverged, ctx));

      return 0;
   }


   /* --------------------------------------------------------------------- */
   /* Find the index of the first column with a retained vector. All        */
   /* retained vectors should be together.                                  */
   /* --------------------------------------------------------------------- */

   for (i=0, orderedIndexOfPreviousVecs=restartSize; i<restartSize; i++) {
      if (hVecsPerm[i] == indexOfPreviousVecs) {
         orderedIndexOfPreviousVecs = i;
         break;
      }
   }
   assert(orderedIndexOfPreviousVecs != restartSize
         || indexOfPreviousVecs >= restartSize);

   for (i=orderedIndexOfPreviousVecs+1;
         i<orderedIndexOfPreviousVecs+numPrevRetained; i++)
      assert(hVecsPerm[i-1]+1 == hVecsPerm[i]);

   /* --------------------------------------------------------------------- */
   /* Given the above H, we know the eigenvectors of H will be the          */
   /* canonical basis except for the retained vectors.                      */
   /* --------------------------------------------------------------------- */

   for (j=0; j < restartSize; j++) {
      for (i=0; i < restartSize; i++) {
          hVecs[newldhVecs*j+i] = 0.0;
      }
      hVecs[newldhVecs*j+hVecsPerm[j]] = 1.0;
   }

   /* Apply permutation hVecsPerm to hVals */
   CHKERR(KIND(permute_vecs_RHprimme, permute_vecs_SHprimme)(
         hVals, 1, restartSize, 1, hVecsPerm, ctx));

   /* ---------------------------------------------------------------------- */
   /* Solve the overlap matrix corresponding for the retained vectors to     */ 
   /* compute the coefficient vectors.                                       */
   /* ---------------------------------------------------------------------- */

   CHKERR(solve_H_SHprimme(&H[ldH * indexOfPreviousVecs + indexOfPreviousVecs],
         numPrevRetained, ldH,
         VtBV ? &VtBV[ldVtBV * (nLocked + indexOfPreviousVecs) + nLocked +
                      indexOfPreviousVecs]
              : NULL,
         ldVtBV, NULL, 0, NULL, 0, NULL, 0, NULL, 0,
         &hVecs[newldhVecs * orderedIndexOfPreviousVecs + indexOfPreviousVecs],
         newldhVecs, &hVals[orderedIndexOfPreviousVecs], NULL, numConverged, ctx));
   return 0;
}

/*******************************************************************************
 * Function restart_refined - This routine is used to recompute the QR decomposition
 *    of W (=A*V), after V being replaced by V*hVecs. The coefficient vectors
 *    except hVecs(indexOfPrevRetained:indexOfPrevRetained+numPrevRetained) are
 *    right singular vectors of R. The output R for these columns will be
 *    a diagonal matrix with the singular values as diagonal elements. The
 *    output Q for these columns will be the update with the left singular
 *    vectors. The other columns of R and Q are recomputed properly without
 *    requiring explicitly recompute the QR factorization.
 *
 *    Also H = V'*A*V and QtV = Q'*V are recomputed properly.
 *   
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * nLocal           Number of rows of V assigned to the node
 *
 * basisSize        Size of the basis V
 *
 * numConverged     The number of converged eigenpairs without locking
 *
 * numPrevRetained  The number of coefficient vectors from previousHVecs in hVecs
 *  
 * indexOfPreviousVecs The first column in hVecs that has a vector from previousHVecs
 *
 * indexOfPreviousVecsBeforeRestart The number of columns in hVecs that previousHVecs
 *                  where orthogonalized against
 *
 * rwork            Real work array
 *
 * rworkSize        Size of rwork
 *                  
 * iwork            Integer work array
 *
 * primme           Structure containing various solver parameters
 *
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * restartSize      The number of vectors to restart with. If negative,
 *                  use dynamic thick restart.
 *
 * V                The orthonormal basis. After restart, contains Ritz vectors
 *                  plus the orthogonal components from numPrevRetained Ritz 
 *                  vectors from the penultimate step.
 *
 * ldV              The leading dimension of V
 *
 * W                A*V
 *
 * ldW              The leading dimension of W
 *
 * H                The projection V'*A*V
 *
 * ldH              The leading dimension of H
 *
 * Q, R             The factors of the QR decomposition of (A - targetShift*B)*V
 *
 * ldQ, ldR         The leading dimension of Q and R
 *
 * hU               The left singular vectors of R
 *
 * fQtQ             The Cholesky factor of QtQ
 *
 * ldfQtQ           The leading dimension of fQtQ
 *
 * ldhU             The leading dimension of the input hU
 *
 * newldhU          The leading dimension of the output hU
 *
 * hVecs            The eigenvectors of H
 *
 * ldhVecs          The leading dimension of the input hVecs
 *
 * newldhVecs       The leading dimension of the output hVecs
 *
 * hVals            The eigenvalues of H
 *
 * hSVals           The singular values of R
 *
 * hVecsPerm        The permutation that orders hVals and hVecs as primme.target
 *
 * restartPerm      The permutation applied to the columns of hVecs before restarting
 *
 * flags            Array indicating the convergence of the Ritz vectors
 *
 * targetShiftIndex The target shift used in (A - targetShift*B) = Q*R
 *
 * numArbitraryVecs On input, the number of coefficients vectors that do
 *                  not correspond to solutions of the projected problem.
 *                  They appear in the first numArbitraryVecs positions of hVecs.
 *                  On output, the number of such vectors that are in the
 *                  restarted basis.
 * 
 * Return value
 * ------------
 * Error code: 0 upon success
 *            -1 eigenvalues of submatrix could not be computed
 *
 ******************************************************************************/

STATIC int restart_refined(SCALAR *V, PRIMME_INT ldV, SCALAR *W, PRIMME_INT ldW,
      SCALAR *BV, PRIMME_INT ldBV, HSCALAR *H, int ldH, SCALAR *Q,
      PRIMME_INT ldQ, PRIMME_INT nLocal, HSCALAR *R, int ldR, HSCALAR *QtQ,
      int ldQtQ, HSCALAR *fQtQ, int ldfQtQ, HSCALAR *VtBV, int ldVtBV,
      HSCALAR *hU, int ldhU, int newldhU, int indexOfPreviousVecsBeforeRestart,
      HSCALAR *hVecs, int ldhVecs, int newldhVecs, HEVAL *hVals, HREAL *hSVals,
      int *restartPerm, int *hVecsPerm, int restartSize, int basisSize,
      int numPrevRetained, int indexOfPreviousVecs, int *targetShiftIndex,
      int numConverged, int *numArbitraryVecs, HSCALAR *hVecsRot,
      int ldhVecsRot, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, j;          /* Loop variables                                       */
   int mhVecsRot0;    /* Number of rows of hVecsRot0                          */
   int newNumArbitraryVecs; /* arbitrary vectors after restarting             */
   int nRegular;      /* restarted vectors that weren't retained              */
   HSCALAR *hVecsRot0; /* Part of hVecsRot                                   */
   HSCALAR *RPrevhVecs;
   double aNorm = primme?max(primme->aNorm, primme->stats.estimateLargestSVal):0.0;

    /* ------------------------------- */
   /* Replace H by hVecs' * H * hVecs */
   /* ------------------------------- */

#ifdef USE_HERMITIAN
   if (H && primme->orth == primme_orth_implicit_I) {
      CHKERR(compute_submatrix_SHprimme(hVecs, restartSize, ldhVecs, H,
            basisSize, ldH, 1 /* Hermitian */, H, ldH, ctx));
   }
#endif /* USE_HERMITIAN */

   /* -------------------------------------- */
   /* Quick exit if the target has changed   */
   /* -------------------------------------- */

   /* NOTE: Force to pass the next condition if you want to rebuild the QR    */
   /* factorization at every restart.                                         */

   /* NOTE: keep the same condition here as in main_iter */

   int nLocked = primme->numOrthoConst + (primme->locking ? numConverged : 0);
   double eps_matrix;
   CHKERR(machineEpsMatrix_Sprimme(&eps_matrix, ctx));
   if (*targetShiftIndex < 0 ||
         fabs(primme->targetShifts[*targetShiftIndex] -
               primme->targetShifts[min(primme->numTargetShifts - 1,
                     numConverged)]) > eps_matrix * aNorm) {

      *targetShiftIndex = min(primme->numTargetShifts-1, numConverged);

      int nQ = 0;
      CHKERR(update_Q_Sprimme(BV ? BV : V, nLocal, BV ? ldBV : ldV, W, ldW, Q,
            ldQ, R, ldR, QtQ, ldQtQ, fQtQ, ldfQtQ,
            primme->targetShifts[*targetShiftIndex], 0, restartSize, &nQ, ctx));
      CHKERRM(restartSize != nQ, -1, "Not supported deficient QR");

      CHKERR(solve_H_SHprimme(H, restartSize, ldH,
            VtBV ? &VtBV[ldVtBV * nLocked + nLocked] : NULL, ldVtBV, R, ldR,
            NULL, 0, QtQ, ldQtQ, hU, newldhU, hVecs, newldhVecs, hVals,
            hSVals, numConverged, ctx));

      *numArbitraryVecs = 0;

      return 0;
   }

   /* restartPerm0 = hVecsPerm(restartPerm)                             */

   int *restartPerm0;
   CHKERR(Num_malloc_iprimme(restartSize, &restartPerm0, ctx));
   for (i=0; i<restartSize; i++) {
      restartPerm0[i] = restartPerm[hVecsPerm[i]];
      assert(restartPerm0[i] >= 0 && restartPerm0[i] < basisSize);
   }

   /* ----------------------------------------------------------------- */
   /* Update numArbitraryVecs as the number of arbitrary vectors in     */
   /* the restarted basis.                                              */
   /* ----------------------------------------------------------------- */

   for (i=newNumArbitraryVecs=0; i < restartSize-numPrevRetained; i++) {
      if (restartPerm0[i] < *numArbitraryVecs) {
         newNumArbitraryVecs++;
      }
   }

   /* Check all arbitrary vectors will be together and at the beginning after restart */

   for (i=0; i<newNumArbitraryVecs; i++)
      assert(restartPerm0[i] < *numArbitraryVecs);

   /* -------------------------------------------------------------------- */
   /* Note that hVecs=[hV*hVecsRot prevhVecs], where hV are the right      */
   /* singular vectors of R and prevhVecs are the retained coefficient     */
   /* vectors after orthogonalizing them against the first                 */
   /* indexOfPreviousVecsBeforeRestart coefficient vectors in hVecs.       */
   /*                                                                      */
   /* In restart V is replaced by by V*Y and W by W*Y.                     */
   /* To restart the QR decomposition of (A-\tau I)*V, Q is replaced by    */
   /* Q*Z and R by T, where Z and T are the QR decomposition of R*Y.       */
   /* We avoid computing R*Y because the errors go into the direction of   */
   /* of the largest singular values. Note that                            */
   /*                                                                      */
   /* R*Y = R*[ hV*hVecsRot prevhVecs ]                                    */
   /*     = [ hU*diag(hSVals)*hVecsRot R*prevhVecs ].                      */
   /*                                                                      */
   /* Note that the code works with R(:,restartPerm0) instead of R,        */
   /* because the new arbitrary columns are all together at the leading    */
   /* columns of R(:,restartPerm0).                                        */
   /* -------------------------------------------------------------------- */


   /* RPrevhVecs = R * hVecs(indexOfPreviousVecs:indexOfPreviousVecs+numPrevRetained) */

   CHKERR(Num_malloc_SHprimme((size_t)numPrevRetained * basisSize, &RPrevhVecs,
                              ctx));
   CHKERR(Num_gemm_SHprimme("N", "N", basisSize, numPrevRetained, basisSize,
         1.0, R, ldR, &hVecs[ldhVecs * indexOfPreviousVecs], ldhVecs, 0.0,
         RPrevhVecs, basisSize, ctx));

   /* Do hVecsRot0 = diag(hSvals)*hVecsRot(hVecsPerm(restartPerm)) limited */
   /* to the columns 0:newNumArbitraryVecs-1 (in 3 steps)                  */
   
   nRegular = restartSize - numPrevRetained;
   for (i=0, mhVecsRot0=*numArbitraryVecs; i < nRegular; i++)
      mhVecsRot0 = max(mhVecsRot0, restartPerm0[i]+1);
   CHKERR(Num_malloc_SHprimme((size_t)mhVecsRot0 * nRegular, &hVecsRot0, ctx));

   /* 1) hVecsRot0 = hVecsRot(hVecsPerm(restartPerm(0:newNumArbitraryVecs-1))) */
   CHKERR(Num_zero_matrix_SHprimme(
         hVecsRot0, mhVecsRot0, nRegular, mhVecsRot0, ctx));
   CHKERR(Num_copy_matrix_columns_SHprimme(hVecsRot, *numArbitraryVecs,
         restartPerm0, newNumArbitraryVecs, ldhVecsRot, hVecsRot0, NULL,
         mhVecsRot0, ctx));

   /* 2) hVecsRot0 = diag(hSVals)*hVecsRot0 */
   for (i=0; i<newNumArbitraryVecs; i++) {
      for (j=0; j<*numArbitraryVecs; j++) {
         hVecsRot0[mhVecsRot0*i+j] *= hSVals[j];
      }
   }

   /* 3) hVecsRot0(:,c) = diag(hSVals)*I(:,restartPerm0(c))  */
   /*  for c = newNumArbVecs:nRegular-1                   */
   for (i=newNumArbitraryVecs; i<nRegular; i++) {
      hVecsRot0[mhVecsRot0*i+restartPerm0[i]] = hSVals[restartPerm0[i]];
   }

   /* R = zeros() */
   Num_zero_matrix_SHprimme(R, primme->maxBasisSize, primme->maxBasisSize,
         ldR, ctx);
   
   /* [hVecsRot0, R] = ortho(hVecsRot0) */
   CHKERR(Bortho_local_SHprimme(hVecsRot0, mhVecsRot0, R, ldR, 0, nRegular - 1,
         NULL, 0, 0, mhVecsRot0, QtQ, ldQtQ, primme->iseed, ctx));

   /* hU = hU * hVecsRot0 */
   HSCALAR *rwork;
   CHKERR(Num_malloc_SHprimme((size_t)basisSize * nRegular, &rwork, ctx));
   CHKERR(Num_gemm_SHprimme("N", "N", basisSize, nRegular, mhVecsRot0, 1.0, hU,
         ldhU, hVecsRot0, mhVecsRot0, 0.0, rwork, basisSize, ctx));
   CHKERR(Num_copy_matrix_SHprimme(
         rwork, basisSize, nRegular, basisSize, hU, basisSize, ctx));
   CHKERR(Num_free_SHprimme(rwork, ctx));
   CHKERR(Num_free_SHprimme(hVecsRot0, ctx));

   /* hU = C'\hU, where C is the Cholesky factor of QtQ */

   if (QtQ) {
      CHKERR(Num_trsm_SHprimme("R", "U", "N", "N", basisSize, nRegular, 1.0,
            fQtQ, ldfQtQ, hU, ldhU, ctx));
   }

   /* hU = [hU RPrevhVecs] */
   CHKERR(Num_copy_matrix_SHprimme(RPrevhVecs, basisSize, numPrevRetained,
         basisSize, &hU[basisSize * nRegular], basisSize, ctx));
   CHKERR(Num_free_SHprimme(RPrevhVecs, ctx));

   /* [hU, R] = ortho(hU, nRegular:nRegular+numPrevRetained-1) */
   CHKERR(Bortho_local_SHprimme(hU, basisSize, R, ldR, nRegular,
         nRegular + numPrevRetained - 1, NULL, 0, 0, basisSize, QtQ, ldQtQ,
         primme->iseed, ctx));

   /* Zero upper triangular part of R(:,newNumArbVecs:nRegular-1) for the     */
   /* columns that restartPerm0[i] >= numArbVecs                              */
   for (i=newNumArbitraryVecs; i<nRegular; i++) {
      if (restartPerm0[i] >= *numArbitraryVecs) {
         for (j=0; j<=i; j++) {
            R[ldR*i+j] = 0.0;
         }
         R[ldR*i+i] = hSVals[restartPerm0[i]];
      }
   }
   CHKERR(Num_free_iprimme(restartPerm0, ctx));

   /* When the retained coefficient vectors were orthogonalized against all   */
   /* arbitrary vectors then R*prevhVecs is orthogonal to                     */
   /* hU(0:indexOfPreviousVecsBeforeRestart) and that rows of R corresponding */
   /* to the retained vectors should be zero.                                 */

    if (*numArbitraryVecs <= indexOfPreviousVecsBeforeRestart) {
      /* Zero R(0:nRegular-1,nRegular:restartSize) */
      Num_zero_matrix_SHprimme(&R[ldR*nRegular], nRegular, numPrevRetained,
            ldR, ctx);
   }

   /* ----------------------------------- */
   /* Restart Q by replacing it with Q*hU */
   /* ----------------------------------- */

   CHKERR(Num_update_VWXR_Sprimme(Q, NULL, NULL, nLocal, basisSize, ldQ, hU,
            restartSize,
            basisSize, NULL,
            Q, 0, restartSize, ldQ,
            NULL, 0, 0, 0,
            NULL, 0, 0, 0,
            NULL, 0, 0, 0,
            NULL, 0, 0, 0, NULL,
            NULL, 0, 0, 0,
            NULL, 0, 0, 0,
            NULL, 0, 0, 0,
            NULL, 0, 0,
            QtQ, restartSize,ldQtQ,
            NULL, 0, 0,
            NULL, 0, 0,
            ctx));

   /* Update fQtQ */

   if (QtQ) {
      CHKERR(update_cholesky_Sprimme(
            QtQ, ldQtQ, fQtQ, ldfQtQ, 0, restartSize, ctx));
   }

   /* ---------------------------------------------------------------------- */
   /* R may lost the triangular structure after the previous permutation, so */
   /* the easiest way to recompute the projected matrices hVecs, hU and      */
   /* and hSVals is to call solve_H. For LOBPCG style, it means to           */
   /* solve twice the projected problem, and it may be a performance issue.  */
   /*                                                                        */
   /* TODO: solve only the columns of R corresponding to the new             */
   /*       arbitrary vectors and the retained vectors.                      */
   /* ---------------------------------------------------------------------- */

   /* Compute the singular value decomposition of R.                         */
   /* Notice that hVals aren't changed. They are updated as a permutation    */
   /* over the original hVals.                                               */
   /* NOTE: the eigenvalues of the retained pairs aren't correctly computed. */

   HEVAL *dummyhVals;
   CHKERR(KIND(Num_malloc_RHprimme, Num_malloc_SHprimme)(
         restartSize, &dummyhVals, ctx));
   CHKERR(solve_H_SHprimme(H, restartSize, ldH,
         VtBV ? &VtBV[nLocked * ldVtBV + nLocked] : NULL, ldVtBV, R, ldR, NULL,
         0, QtQ, ldQtQ, hU, newldhU, hVecs, newldhVecs, dummyhVals, hSVals,
         numConverged, ctx));
   CHKERR(KIND(Num_free_RHprimme,Num_free_SHprimme)(dummyhVals, ctx));

   CHKERR(KIND(permute_vecs_RHprimme, permute_vecs_SHprimme)(
         hVals, 1, restartSize, 1, hVecsPerm, ctx));

   /* ---------------------------------------------------------------------- */
   /* Permute back the columns of R                                          */
   /* ---------------------------------------------------------------------- */

   /* hVecsPerm(invhVecsPerm) = 1:n */
   int *invhVecsPerm;
   CHKERR(Num_malloc_iprimme(restartSize, &invhVecsPerm, ctx));
   for (i=0; i<restartSize; i++) {
      invhVecsPerm[hVecsPerm[i]] = i;
   }

   /* R = R(:, invhVecsPerm) */
   CHKERR(permute_vecs_SHprimme(
         R, restartSize, restartSize, ldR, invhVecsPerm, ctx));
   CHKERR(Num_free_iprimme(invhVecsPerm, ctx));

   /* ----------------------------------------------------------------- */
   /* After all the changes in hVecs and R new arbitrary vectors may    */
   /* have been introduced. When the retained coefficient vectors are   */
   /* orthogonalized against all arbitrary vectors then the new number  */
   /* of arbitrary vectors is at most the largest index out of order    */
   /* in hVecsPerm, plus one. Otherwise the easiest way is to consider  */
   /* all restarted vectors as arbitrary vectors.                       */
   /* ----------------------------------------------------------------- */

   if (*numArbitraryVecs <= indexOfPreviousVecsBeforeRestart) {
      for (i=*numArbitraryVecs=newNumArbitraryVecs; i<restartSize; i++)
         if (hVecsPerm[i] != i) *numArbitraryVecs=i+1;
   }
   else {
      *numArbitraryVecs = restartSize;
   }

   /* ----------------------------------------------------------------------- */
   /* We want to compute hVecsRot so that it corresponds to the rotation of   */
   /* the right singular vectors of R that produces I(:,hVecsPerm), because   */
   /* the desired coefficient vectors were already decided in the last        */
   /* iteration. The current hVecs is the right singular vectors of           */
   /* R*I(:,hVecsPerm). Notice R*I(:,hVecsPerm) = hU * S * hV'*I(:,hVecsPerm).*/
   /* Setting hVecsRot as current hVecs'=hV'*I(:,hVecsPerm) will satisfy that */
   /* hVecs = hV * hVecsRot, with hVecs=I(:,hVecsPerm).                       */
   /* ----------------------------------------------------------------------- */

   /* hVecsRot <- hVecs' for arbitrary vectors */

   Num_zero_matrix_SHprimme(hVecsRot, primme->maxBasisSize, primme->maxBasisSize,
         ldhVecsRot, ctx);
   for (j=0; j < *numArbitraryVecs; j++) {
      for (i=0; i<restartSize; i++) {
         hVecsRot[ldhVecsRot*j+i] = CONJ(hVecs[newldhVecs*i+j]);
      }
   }
 
   /* hVecs <- I for arbitrary vectors */

   Num_zero_matrix_SHprimme(hVecs, restartSize, *numArbitraryVecs, newldhVecs, ctx);
   for (j=0; j < *numArbitraryVecs; j++) {
      hVecs[newldhVecs*j+hVecsPerm[j]] = 1.0;
   }

   return 0;
}

/*******************************************************************************
 * Function restart_harmonic - This routine is used to recompute the QR decomposition
 *    of W (=A*V), after V being replaced by V*hVecs. 
 *    Also H = V'*A*V and QtV = Q'*V are recomputed properly.
 *   
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * nLocal           Number of rows of V assigned to the node
 *
 * basisSize        Size of the basis V
 *
 * numConverged     The number of converged eigenpairs without locking
 *
 * numPrevRetained  The number of coefficient vectors from previousHVecs in hVecs
 *  
 * indexOfPreviousVecs The first column in hVecs that has a vector from previousHVecs
 *
 * rwork            Real work array
 *
 * rworkSize        Size of rwork
 *                  
 * iwork            Integer work array
 *
 * primme           Structure containing various solver parameters
 *
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * restartSize      The number of vectors to restart with. If negative,
 *                  use dynamic thick restart.
 *
 * V                The orthonormal basis. After restart, contains Ritz vectors
 *                  plus the orthogonal components from numPrevRetained Ritz 
 *                  vectors from the penultimate step.
 *
 * ldV              The leading dimension of V
 *
 * W                A*V
 *
 * ldW              The leading dimension of W
 *
 * H                The projection V'*A*V
 *
 * ldH              The leading dimension of H
 *
 * Q, R             The factors of the QR decomposition of (A - targetShift*B)*V
 *
 * ldQ, ldR         The leading dimension of Q and R
 *
 * fQtQ             The Cholesky factor of QtQ
 *
 * ldfQtQ           The leading dimension of fQtQ
 *
 * QtV              = Q'*V
 *
 * ldQtV            The leading dimension of QtV
 *
 * hU               The eigenvectors of QtV/R
 *
 * ldhU             The leading dimension of the input hU
 *
 * newldhU          The leading dimension of the output hU
 *
 * hVecs            The eigenvectors of H
 *
 * ldhVecs          The leading dimension of the input hVecs
 *
 * newldhVecs       The leading dimension of the output hVecs
 *
 * hVals            The eigenvalues of H
 *
 * hSVals           The singular values of R
 *
 * hVecsPerm        The permutation that orders hVals and hVecs as primme.target
 *
 * restartPerm      The permutation applied to the columns of hVecs before restarting
 *
 * flags            Array indicating the convergence of the Ritz vectors
 *
 * targetShiftIndex The target shift used in (A - targetShift*B) = Q*R
 *
 * numArbitraryVecs On input, the number of coefficients vectors that do
 *                  not correspond to solutions of the projected problem.
 *                  They appear in the first numArbitraryVecs positions of hVecs.
 *                  On output, the number of such vectors that are in the
 *                  restarted basis.
 * 
 * Return value
 * ------------
 * Error code: 0 upon success
 *            -1 eigenvalues of submatrix could not be computed
 *
 ******************************************************************************/

STATIC int restart_harmonic(SCALAR *V, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, SCALAR *BV, PRIMME_INT ldBV, HSCALAR *H, int ldH,
      SCALAR *Q, PRIMME_INT ldQ, PRIMME_INT nLocal, HSCALAR *R, int ldR,
      HSCALAR *QtV, int ldQtV, HSCALAR *QtQ, int ldQtQ, HSCALAR *fQtQ,
      int ldfQtQ, HSCALAR *VtBV, int ldVtBV, HSCALAR *hU, int ldhU, int newldhU,
      HSCALAR *hVecs, int ldhVecs, int newldhVecs, HEVAL *hVals, HREAL *hSVals,
      int *restartPerm, int *hVecsPerm, int restartSize, int basisSize,
      int numPrevRetained, int indexOfPreviousVecs, int *targetShiftIndex,
      int numConverged, int *numArbitraryVecs, HSCALAR *hVecsRot,
      int ldhVecsRot, primme_context ctx) {

   primme_params *primme = ctx.primme;

   (void)ldhU; /* unused parameter */
   (void)restartPerm; /* unused parameter */
   (void)hVecsPerm; /* unused parameter */
   (void)numPrevRetained; /* unused parameter */
   (void)indexOfPreviousVecs; /* unused parameter */
   (void)hVecsRot; /* unused parameter */
   (void)ldhVecsRot; /* unused parameter */

   /* ------------------------------- */
   /* Replace H by hVecs' * H * hVecs */
   /* ------------------------------- */

#ifdef USE_HERMITIAN
   if (primme->orth == primme_orth_implicit_I) {
      CHKERR(compute_submatrix_SHprimme(hVecs, restartSize, ldhVecs, H,
            basisSize, ldH, 1 /* Hermitian */, H, ldH, ctx));
   }
#else
   (void)ldhVecs;
   (void)basisSize;
#endif /* USE_HERMITIAN */

   /* ------------------------------- */
   /* Update targetShiftIndex         */
   /* ------------------------------- */

   *targetShiftIndex = min(primme->numTargetShifts-1, numConverged);

   /* ------------------------------- */
   /* Compute QR                      */
   /* ------------------------------- */

   int nQ = 0;
   CHKERR(update_Q_Sprimme(BV ? BV : V, nLocal, BV ? ldBV : ldV, W, ldW, Q, ldQ,
         R, ldR, QtQ, ldQtQ, fQtQ, ldfQtQ,
         primme->targetShifts[*targetShiftIndex], 0, restartSize, &nQ, ctx));
   CHKERRM(restartSize != nQ, -1, "Not supported deficient QR");

   /* ------------------------------- */
   /* Update QtV                      */
   /* ------------------------------- */

   CHKERR(update_projection_Sprimme(Q, ldQ, V, ldV, QtV, ldQtV, nLocal, 0,
                                    restartSize, 0 /*unsymmetric*/, ctx));

   /* ------------------------------- */
   /* Solve the projected problem     */
   /* ------------------------------- */

   int nLocked = primme->numOrthoConst + (primme->locking ? numConverged : 0);
   CHKERR(solve_H_SHprimme(H, restartSize, ldH,
         VtBV ? &VtBV[ldVtBV * nLocked + nLocked] : NULL, ldVtBV, R, ldR, QtV,
         ldQtV, QtQ, ldQtQ, hU, newldhU, hVecs, newldhVecs, hVals, hSVals,
         numConverged, ctx));

   *numArbitraryVecs = 0;

   return 0;
}


/*******************************************************************************
 * Subroutine ortho_coefficient_vectors - Orthogonalize properly the columns of
 *    hVecs from indexOfPreviousVecs plus numPrevRetained.
 *    
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * hVecs            The eigenvectors of H
 * ldhVecs          The leading dimension of the input hVecs
 * basisSize        Rows in hVecs and size of R and hU
 * R                The factors of the QR decomposition of (A - targetShift*B)*V
 * ldR              The leading dimension of R
 * outR             Rotations of orthogonalizing hVecs
 * iwork            Integer work array
 * rwork            Work array
 * rworkSize        Length of the work array
 *
 ******************************************************************************/

STATIC int ortho_coefficient_vectors_Sprimme(HSCALAR *hVecs, int basisSize,
      int ldhVecs, int indexOfPreviousVecs, HSCALAR *VtBV, int ldVtBV,
      HSCALAR *prevhVecs, int nprevhVecs, int ldprevhVecs, int *flags,
      int *numPrevRetained, primme_context ctx) {

   HSCALAR *rwork; /* auxiliary space to broadcast the retained vectors */
   CHKERR(Num_malloc_SHprimme(basisSize*(*numPrevRetained)+1, &rwork, ctx));

   if (ctx.procID == 0) {

      /* First, retain coefficient vectors corresponding to current block */
      /* vectors.  If all of those have been retained, then retain the    */ 
      /* the next unconverged coefficient vectors beyond iev[blockSize-1].*/

      int retained = 0; 
      int i;
      for (i = 0; i < nprevhVecs && retained < *numPrevRetained &&
                  indexOfPreviousVecs + retained < basisSize;
            i++) {
         /* Skip converged pairs in soft locking */

         if (ctx.primme->locking == 0 && flags[i] != UNCONVERGED) continue;

         /* Avoid that ortho replaces linear dependent vectors by random
          * vectors. The random vectors make the restarting W with large
          * singular value. This change has shown benefit finding the smallest
          * singular values.  */

         HSCALAR R = 0.0;
         CHKERR(Bortho_local_SHprimme(&prevhVecs[ldprevhVecs * i], ldprevhVecs,
               &R, 1, 0, 0, hVecs, ldhVecs, indexOfPreviousVecs + retained,
               basisSize, VtBV, ldVtBV, ctx.primme->iseed, ctx));

         if (ABS(R) < MACHINE_EPSILON * sqrt(retained + 1.0)) continue;

         CHKERR(Num_copy_matrix_SHprimme(&prevhVecs[ldprevhVecs * i], basisSize,
               1, ldprevhVecs,
               &hVecs[ldhVecs * (indexOfPreviousVecs + retained)], ldhVecs,
               ctx));
         retained++;
      }

      PRINTF(5, "retain_previous: numPrevRetained: %d", retained);

      rwork[0] = (HSCALAR)retained;
      CHKERR(Num_copy_matrix_SHprimme(&hVecs[ldhVecs * indexOfPreviousVecs],
            basisSize, retained, ldhVecs, rwork + 1, basisSize, ctx));
      CHKERR(Num_zero_matrix_SHprimme(&rwork[1 + basisSize * retained],
            basisSize, *numPrevRetained - retained, basisSize, ctx));
   }

   /* Broadcast hVecs(indexOfPreviousVecs:indexOfPreviousVecs+numPrevRetained) */

   CHKERR(broadcast_SHprimme(rwork, basisSize * (*numPrevRetained) + 1, ctx));
   *numPrevRetained = (int)REAL_PART(rwork[0]);
   CHKERR(Num_copy_matrix_SHprimme(rwork + 1, basisSize, *numPrevRetained,
         basisSize, &hVecs[ldhVecs * indexOfPreviousVecs], ldhVecs, ctx));

   CHKERR(Num_free_SHprimme(rwork, ctx));

  return 0;
}

/******************************************************************************
 * Function compute_residual_columns - This subroutine performs the next
 *    operations in a cache-friendly way:
 *
 *    X = X(p); AX = AX(p); BX = BX(p)
 *    j = k = 0; XD = RD = []
 *    for i=0:nd-1
 *       if pd(i) == j+io0
 *          XD = [XD XO(j)]; RD = [RD RO(j)]; BXD = [BXD BXO(j)]; j++
 *       else
 *          XD = [XD X(k]; RD = [RD AX(k) - evals(p(k))*(BX?BX:X)(k)]; k++
 *       end if
 *    end for
 *
 *           n        nd             no
 * X:  [-----|-------------]           (input/output)
 * XO:                      [--------| (input)
 * XD:       [--------|                (output)
 *
 * NOTE: X and XD *can* overlap, but X(0:n-1) and XD *cannot* overlap (same for R and RD)
 *       XO and XD *can* overlap (same for RO and RD)
 *       p should be a list of increasing indices
 *       pd should be a merge of two increasing lists
 *
 * PARAMETERS
 * ---------------------------
 * m           The number of rows of matrices x, Ax, xo, ro, xd and rd
 * evals       The values to compute the residual
 * x           The matrix that does x = x(p)
 * Bx          The matrix that does Bx = Bx(p)
 * n           The number of columns of the output x
 * p           The columns to copy back to x and Ax
 * ldx         The leading dimension of x
 * Ax          The matrix that does Ax = Ax(p)
 * ldAx        The leading dimension of Ax
 * xo          Alternative source of columns for xd
 * Bxo         Alternative source of columns for xd
 * no          The number of columns in xo
 * ldxo        The leading dimension of xo
 * io0         The index of the first column in xo
 * ro          Alternative source of columns for rd
 * ldro        The leading dimension of ro
 * xd          The matrix that will have columns from x and xo
 * Bxd         The matrix that will have columns from x and xo
 * nd          The maximum size of xd
 * pd          The indices of the columns to generate xd and rd
 * ldxd        The leading dimension of xd
 * rd          The matrix that will have columns from r and ro
 * ldrd        The leading dimension of rd
 * rwork       Workspace
 * lrwork      The size of rwork
 *
 ******************************************************************************/

STATIC int compute_residual_columns(PRIMME_INT m, HEVAL *evals, SCALAR *x,
      SCALAR *Bx, int n, int *p, PRIMME_INT ldx, SCALAR *Ax, PRIMME_INT ldAx,
      SCALAR *xo, SCALAR *Bxo, int no, PRIMME_INT ldxo, int io0, SCALAR *ro,
      PRIMME_INT ldro, SCALAR *xd, SCALAR *Bxd, int nd, int *pd,
      PRIMME_INT ldxd, SCALAR *rd, PRIMME_INT ldrd, primme_context ctx) {

   int i, id, k, io, M=min(m,PRIMME_BLOCK_SIZE);

   /* Quick exit */

   if (n == 0) {
      CHKERR(Num_copy_matrix_Sprimme(xo, m, min(no,nd), ldxo, xd, ldxd, ctx));
      if (Bxd) {
         CHKERR(Num_copy_matrix_Sprimme(
               Bxo, m, min(no, nd), ldxo, Bxd, ldxd, ctx));
      }
      CHKERR(Num_copy_matrix_Sprimme(ro, m, min(no, nd), ldro, rd, ldrd, ctx));
      return 0;
   }

   SCALAR *X0, *R0, *BX0 = NULL;
   CHKERR(Num_malloc_Sprimme(nd*M, &X0, ctx));
   CHKERR(Num_malloc_Sprimme(nd*M, &R0, ctx));
   if (Bxo) CHKERR(Num_malloc_Sprimme(nd*M, &BX0, ctx));

   for (k=0; k<m; k+=M, M=min(M,m-k)) {
      for (i=id=io=0; i < n || id < nd; id++) {
         if (id < nd && io < no && pd[id] == io+io0) {
            CHKERR(Num_copy_matrix_Sprimme(
                  &xo[io * ldxo + k], M, 1, ldxo, &X0[id * M], M, ctx));
            if (Bxo)
               CHKERR(Num_copy_matrix_Sprimme(
                     &Bxo[io * ldxo + k], M, 1, ldxo, &BX0[id * M], M, ctx));
            CHKERR(Num_copy_matrix_Sprimme(
                  &ro[io * ldro + k], M, 1, ldro, &R0[id * M], M, ctx));
            io++;
         }
         else {
            assert(id >= nd || i < n);
            CHKERR(Num_copy_matrix_Sprimme(
                  &x[p[i] * ldx + k], M, 1, ldx, &x[i * ldx + k], ldx, ctx));
            if (Bx)
               CHKERR(Num_copy_matrix_Sprimme(&Bx[p[i] * ldx + k], M, 1, ldx,
                     &Bx[i * ldx + k], ldx, ctx));
            CHKERR(Num_copy_matrix_Sprimme(&Ax[p[i] * ldAx + k], M, 1, ldAx,
                  &Ax[i * ldAx + k], ldAx, ctx));
            if (id < nd) {
               CHKERR(Num_copy_matrix_Sprimme(
                     &x[p[i] * ldx + k], M, 1, ldx, &X0[id * M], M, ctx));
               if (Bx)
                  CHKERR(Num_copy_matrix_Sprimme(
                        &Bx[p[i] * ldx + k], M, 1, ldx, &BX0[id * M], M, ctx));
               CHKERR(Num_compute_residuals_Sprimme(M, 1, &evals[p[i]],
                     Bx ? &Bx[p[i] * ldx + k] : &x[p[i] * ldx + k], ldx,
                     &Ax[p[i] * ldAx + k], ldAx, &R0[id * M], M, ctx));
            }
            i++;
         }
      }
      assert(id >= nd);
      CHKERR(Num_copy_matrix_Sprimme(X0, M, nd, M, &xd[k], ldxd, ctx));
      if (Bxd) {
         CHKERR(Num_copy_matrix_Sprimme(BX0, M, nd, M, &Bxd[k], ldxd, ctx));
      }
      CHKERR(Num_copy_matrix_Sprimme(R0, M, nd, M, &rd[k], ldrd, ctx));
   }

   CHKERR(Num_free_Sprimme(X0, ctx));
   CHKERR(Num_free_Sprimme(R0, ctx));
   if (Bxo) CHKERR(Num_free_Sprimme(BX0, ctx));

   return 0;
}

#endif /* SUPPORTED_TYPE */
