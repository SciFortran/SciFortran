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
 * File: main_iter.c
 *
 * Purpose - This is the main Davidson-type iteration 
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/main_iter.c"
#endif


#include "numerical.h"
#include "template_normal.h"
#include "common_eigs.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "main_iter.h"
#include "convergence.h"
#include "correction.h"
#include "init.h"
#include "ortho.h"
#include "restart.h"
#include "solve_projection.h"
#include "update_projection.h"
#include "update_W.h"
#include "auxiliary_eigs.h"
#include "auxiliary_eigs_normal.h"
#endif

#ifdef SUPPORTED_TYPE

/*----------------------------------------------------------------------------*
 * The following are needed for the Dynamic Method Switching
 *----------------------------------------------------------------------------*/

#ifndef MAIN_ITER_PRIVATE_H
#define MAIN_ITER_PRIVATE_H
typedef struct {
   /* Time measurements for various components of the solver */
   double MV_PR;          /* OPp operator MV+PR time                          */
   double MV;             /* MV time only                                     */
   double PR;             /* PRecond only                                     */
   double qmr_only;       /* c_q only the QMR iteration                       */
   double qmr_plus_MV_PR; /* a   QMR plus operators                           */
   double gdk_plus_MV_PR; /* b   GD plus operators                            */
   double gdk_plus_MV;    /* b_m GD plus MV (i.e., GD without correction)     */
                          /* The following two are not currently updated:     */
   double project_locked; /* c_p projection time per locked eigenvector in QMR*/
   double reortho_locked; /* c_g (usu 2*c_p) ortho per locked vector in outer */

   /* Average convergence estimates. Updated at restart/switch/convergence    */
   double gdk_conv_rate;  /* convergence rate of all (|r|/|r0|) seen for GD+k */
   double jdq_conv_rate;  /* convergence rate of all (|r|/|r0|) seen for JDQMR*/
   double JDQMR_slowdown; /* log(gdRate)/log(jdRate) restricted in [1.1, 2.5] */
   double ratio_MV_outer; /* TotalMV/outerIts (average of last two updates)   */

   /* To maintain a convergence rate of all residual reductions seen we need: */
   int    nextReset;      /* When to reset the following averaging sums       */
   double gdk_sum_logResReductions;/* Sumof all log(residual reductions) in GD*/
   double jdq_sum_logResReductions;/* Sumof all log(residual reductions) in JD*/
   double gdk_sum_MV;     /* Total num of MV resulting in these GD reductions */
   double jdq_sum_MV;     /* Total num of MV resulting in these JD reductions */
   int nevals_by_gdk;     /* Number of evals found by GD+k since last reset   */
   int nevals_by_jdq;     /* Number of evals found by JDQMR since last reset  */

   /* Variables to remember MV/its/time/resNorm/etc since last update/switch */
   int numIt_0;           /*Remembers starting outer its/MVs since switched   */
   int numMV_0;           /*   to current method, or since an epair converged */
   double timer_0;        /*Remembers starting time since switched to a method*/
                          /*   or since an epair converged with that method   */
   double time_in_inner;  /*Accumulative time spent in inner iterations       */
                          /*   since last switch or since an epair converged  */
   double resid_0;        /*First residual norm of the convergence of a method*/
                          /*   since last switch or since an epair converged */

   /* Weighted ratio of expected times, for final method recommendation */
   double accum_jdq_gdk;  /*Expected ratio of accumulative times of JDQMR/GD+k*/
   double accum_jdq;      /* Accumulates jdq_times += ratio*(gdk+MV+PR)       */
   double accum_gdk;      /* Accumulates gdk_times += gdk+MV+PR               */

} primme_CostModel;
#endif /* MAIN_ITER_PRIVATE_H */


#if 0
STATIC void displayModel(primme_CostModel *model);
#endif

/******************************************************************************
 * Subroutine main_iter - This routine implements a more general, parallel, 
 *    block (Jacobi)-Davidson outer iteration with a variety of options.
 *
 * A coarse outline of the algorithm performed is as follows:
 *
 *  1. Initialize basis V
 *  2. while (converged Ritz vectors have become unconverged) do
 *  3.    Update W=A*V, and H=V'*A*V, and solve the H eigenproblem
 *  4.    while (not all Ritz pairs have converged) do
 *  5.       while (maxBasisSize has not been reached) 
 *  6.          Adjust the block size if necessary
 *  7.          Compute the residual norms of each block vector and 
 *                check for convergence.  If a Ritz vector has converged,
 *                target an unconverged Ritz vector to replace it.
 *  8.          Solve the correction equation for each residual vector
 *  9.          Insert corrections in V, orthogonalize, update W and H 
 * 10.          Solve the new H eigenproblem to obtain the next Ritz pairs    
 * 11.       endwhile
 * 12.       Restart V with appropriate vectors 
 * 13.       If (locking) Lock out Ritz pairs that converged since last restart
 * 14.    endwhile
 * 15.    if Ritz vectors have become unconverged, reset convergence flags
 * 16. endwhile
 *    
 *
 * OUTPUT arrays and parameters
 * ----------------------------
 * evals    The converged Ritz values.  It is maintained in sorted order.
 *          If locking is not engaged, values are copied to this array only 
 *          upon return.  If locking is engaged, then values are copied to this 
 *          array as they converge.
 *
 * resNorms The residual norms corresponding to the converged Ritz vectors
 *
 * INPUT/OUTPUT arrays and parameters
 * ----------------------------------
 * evecs    Stores initial guesses. Upon return, it contains the converged Ritz
 *          vectors.  If locking is not engaged, then converged Ritz vectors 
 *          are copied to this array just before return.  
 *
 * primme.initSize: On output, it stores the number of converged eigenvectors. 
 *           If smaller than numEvals and locking is used, there are
 *              only primme.initSize vectors in evecs.
 *           Without locking all numEvals approximations are in evecs
 *              but only the initSize ones are converged.
 *           During the execution, access to primme.initSize gives 
 *              the number of converged pairs up to that point. The pairs
 *              are available in evals and evecs, but only when locking is used
 * ret      successful state
 * numRet   number of returned pairs in evals, evecs, and resNorms
 *
 * Return Value
 * ------------
 * error code      
 ******************************************************************************/

TEMPLATE_PLEASE
int main_iter_Sprimme(HEVAL *evals, SCALAR *evecs, PRIMME_INT ldevecs,
      HREAL *resNorms, double startTime, int *ret, int *numRet,
      primme_context ctx) {

   /* Default error is something is wrong in this function */
   *ret = PRIMME_MAIN_ITER_FAILURE;

   primme_params *primme = ctx.primme;
                            /* primme parameters */
   int i;                   /* Loop variable                                 */
   int blockSize;           /* Current block size                            */
   int availableBlockSize;  /* There is enough space in basis for this block */
   int basisSize;           /* Current size of the basis V                   */
   int numLocked;           /* Number of locked Ritz vectors                 */
   int numGuesses;          /* Current number of initial guesses in evecs    */
   int nextGuess;           /* Index of the next initial guess in evecs      */ 
   int numConverged;        /* Number of converged Ritz pairs                */
   int targetShiftIndex;    /* Target shift used in the QR factorization     */
   int recentlyConverged;   /* Number of target Ritz pairs that have         */
                            /*    converged during the current iteration.    */
   int maxRecentlyConverged;/* Maximum value allowed for recentlyConverged   */
   int numConvergedStored;  /* Numb of Ritzvecs temporarily stored in evecs  */
                            /*    to allow for skew projectors w/o locking   */
   int LockingProblem;      /* Flag==1 if practically converged pairs locked */
   int restartLimitReached; /* True when maximum restarts performed          */
   int nprevhVecs;          /* Number of vectors stored in prevhVecs         */
   int numArbitraryVecs;    /* Columns in hVecs computed with RR instead of  */
                            /* the current extraction method.                */
   int maxEvecsSize;        /* Maximum capacity of evecs array               */
   int numPrevRitzVals = 0; /* Size of the prevRitzVals updated in correction*/
   int touch=0;             /* param used in inner solver stopping criteria  */

   int *perm=NULL;          /* Permutation on the locked pairs */
   int *flags;              /* Indicates which Ritz values have converged    */
   int *map;                /* Indices of the closest vector from prevhVecs  */
   int *lockedFlags=NULL;   /* Flags for the locked pairs                    */
   int *ipivot;             /* The pivot for the Mfact factorization of M      */
   int *iev;                /* Evalue index each block vector corresponds to */

   SCALAR *V;               /* Basis vectors                                 */
   PRIMME_INT ldV;          /* The leading dimension of V                    */
   SCALAR *W;               /* Work space storing A*V                        */
   PRIMME_INT ldW;          /* The leading dimension of W                    */
   SCALAR *BV=NULL;         /* Work space storing B*V                        */
   PRIMME_INT ldBV;         /* The leading dimension of BV                   */
   HSCALAR *H;              /* Upper triangular portion of V'*A*V            */
   HSCALAR *VtBV = NULL;    /* Upper triangular portion of V'*B*V            */
   HSCALAR *fVtBV = NULL;   /* Cholesky factor of VtBV                       */
   HSCALAR *QtQ = NULL;     /* Upper triangular portion of Q'*Q              */
   HSCALAR *fQtQ = NULL;    /* Cholesky factor of QtQ                        */
   HSCALAR *M = NULL;       /* The projection Q'*K*B*Q, where Q = [evecs, x]   */
                            /* x is the current Ritz vector and K is a       */
                            /* hermitian preconditioner.                     */
   HSCALAR *Mfact = NULL;     /* The factorization of M=Q'KBQ                   */
   SCALAR *Bevecs = NULL;   /* B*evecs                                       */
   PRIMME_INT ldBevecs;     /* Leading dimension of Bevecs                   */
   SCALAR *evecsHat = NULL; /* K^{-1}*B*evecs                                   */
   PRIMME_INT ldevecsHat=0; /* The leading dimension of evecsHat             */
   HSCALAR *hVecs;          /* Eigenvectors of H                             */
   HSCALAR *hU = NULL;      /* Left singular vectors of R                    */
   HSCALAR *prevhVecs=NULL; /* hVecs from previous iteration                 */

   int numQR;               /* Maximum number of QR factorizations           */
   SCALAR *Q = NULL;        /* QR decompositions for harmonic or refined     */
   PRIMME_INT ldQ;          /* The leading dimension of Q                    */
   HSCALAR *R = NULL;       /* projection: (A-target[i])*V = QR              */
   HSCALAR *QtV = NULL;     /* Q'*V                                          */
   HSCALAR *hVecsRot = NULL; /* transformation of hVecs in arbitrary vectors  */

   HEVAL *hVals;           /* Eigenvalues of H                              */
   HREAL *hSVals = NULL;    /* Singular values of R                          */
   HEVAL *prevRitzVals;     /* Eigenvalues of H at previous outer iteration  */
                            /* by robust shifting algorithm in correction.c  */
   HREAL *basisNorms;       /* Residual norms of basis at pairs              */
   HREAL *blockNorms;       /* Residual norms corresponding to current block */
                            /* vectors.                                      */
   double smallestResNorm;  /* the smallest residual norm in the block       */
   int reset=0;             /* Flag to reset V and W                         */
   int restartsSinceReset=0;/* Restart since last reset of V and W           */
   int wholeSpace=0;        /* search subspace reach max size                */

   /* Runtime measurement variables for dynamic method switching             */
   primme_CostModel CostModel; /* Structure holding the runtime estimates of */
                            /* the parameters of the model.Only visible here */
   double tstart=0.0;       /* Timing variable for accumulative time spent   */
   int maxNumRandoms = 10;  /* We do not allow more than 10 randomizations */

   /* -------------------------------------------------------------- */
   /* Allocate objects                                               */
   /* -------------------------------------------------------------- */

   maxEvecsSize = primme->numOrthoConst + primme->numEvals;
   if (primme->projectionParams.projection != primme_proj_RR) {
      numQR = 1;
   }
   else {
      numQR = 0;
   }

   /* Use leading dimension ldOPs for the large dimension mats: V, W and Q */

   ldV = ldW = ldBV = ldQ = primme->ldOPs;
   ldBevecs = primme->massMatrixMatvec ? primme->ldOPs : ldevecs;
   if (primme->massMatrixMatvec) {
      CHKERR(Num_malloc_Sprimme(ldBV*primme->maxBasisSize, &BV, ctx));
   }
   CHKERR(Num_malloc_Sprimme(primme->ldOPs*primme->maxBasisSize, &V, ctx));
   CHKERR(Num_malloc_Sprimme(primme->ldOPs*primme->maxBasisSize, &W, ctx));
   if (numQR > 0) {
      CHKERR(Num_malloc_Sprimme(primme->ldOPs*primme->maxBasisSize*numQR, &Q, ctx));
      CHKERR(Num_malloc_SHprimme(
          primme->maxBasisSize * primme->maxBasisSize * numQR, &R, ctx));
      CHKERR(Num_malloc_SHprimme(
          primme->maxBasisSize * primme->maxBasisSize * numQR, &hU, ctx));
   }
   if (primme->projectionParams.projection == primme_proj_harmonic) {
      CHKERR(Num_malloc_SHprimme(
          primme->maxBasisSize * primme->maxBasisSize * numQR, &QtV, ctx));
   }
   CHKERR(Num_malloc_SHprimme(primme->maxBasisSize * primme->maxBasisSize, &H,
                              ctx));
   CHKERR(Num_malloc_SHprimme(primme->maxBasisSize * primme->maxBasisSize,
                              &hVecs, ctx));
   int maxRank; /* maximum size of the main space being orthonormalize */
   if (primme->locking) {
      maxRank = primme->numOrthoConst + primme->numEvals + primme->maxBasisSize;
   } else {
      maxRank = primme->numOrthoConst + primme->maxBasisSize;
   }
   if (primme->orth == primme_orth_explicit_I) {
      CHKERR(Num_malloc_SHprimme(maxRank * maxRank, &VtBV, ctx));
      CHKERR(Num_malloc_SHprimme(maxRank * maxRank, &fVtBV, ctx));
   }
   int ldVtBV = maxRank;
   int ldfVtBV = maxRank;
   if (primme->orth == primme_orth_explicit_I && numQR) {
      CHKERR(Num_malloc_SHprimme(
            primme->maxBasisSize * primme->maxBasisSize * numQR, &QtQ, ctx));
      CHKERR(Num_malloc_SHprimme(
            primme->maxBasisSize * primme->maxBasisSize * numQR, &fQtQ, ctx));
   }
   int ldQtQ = primme->maxBasisSize;
   int ldfQtQ = primme->maxBasisSize;
   if (primme->projectionParams.projection == primme_proj_refined
       || primme->projectionParams.projection == primme_proj_harmonic) {
     CHKERR(Num_malloc_SHprimme(
         primme->maxBasisSize * primme->maxBasisSize * numQR, &hVecsRot, ctx));
   }

   if (primme->correctionParams.precondition && 
         primme->correctionParams.maxInnerIterations != 0 &&
         primme->correctionParams.projectors.RightQ &&
         primme->correctionParams.projectors.SkewQ           ) {
      ldevecsHat = primme->nLocal;
      CHKERR(Num_malloc_Sprimme(ldevecsHat * maxEvecsSize, &evecsHat, ctx));
      CHKERR(Num_malloc_SHprimme(maxEvecsSize * maxEvecsSize, &M, ctx));
      CHKERR(Num_malloc_SHprimme(maxEvecsSize * maxEvecsSize, &Mfact, ctx));
   }
   if (primme->massMatrixMatvec && primme->locking) {
      CHKERR(Num_malloc_Sprimme(ldBevecs * maxEvecsSize, &Bevecs, ctx));
   }
   CHKERR(Num_malloc_SHprimme(
         primme->maxBasisSize * primme->maxBasisSize, &prevhVecs, ctx));

   CHKERR(KIND(Num_malloc_RHprimme, Num_malloc_SHprimme)(primme->maxBasisSize,
         &hVals, ctx));
   if (numQR > 0) {
      CHKERR(Num_malloc_RHprimme(primme->maxBasisSize, &hSVals, ctx));
   }
   CHKERR(KIND(Num_malloc_RHprimme, Num_malloc_SHprimme)(
         primme->maxBasisSize + primme->numEvals, &prevRitzVals, ctx));
   CHKERR(Num_malloc_RHprimme(primme->maxBlockSize, &blockNorms, ctx));
   CHKERR(Num_malloc_RHprimme(primme->maxBasisSize, &basisNorms, ctx));
   CHKERR(
         Num_zero_matrix_RHprimme(basisNorms, 1, primme->maxBasisSize, 1, ctx));

   /* Integer objects */

   CHKERR(Num_malloc_iprimme(primme->numEvals, &perm, ctx));
   if (primme->locking) {
      CHKERR(Num_malloc_iprimme(primme->numEvals, &lockedFlags, ctx));
   }
   CHKERR(Num_malloc_iprimme(primme->maxBasisSize, &flags, ctx));
   CHKERR(Num_malloc_iprimme(primme->maxBasisSize, &map, ctx));
   CHKERR(Num_malloc_iprimme(primme->maxBlockSize, &iev, ctx));
   CHKERR(Num_malloc_iprimme(maxEvecsSize, &ipivot, ctx));

   /* -------------------------------------------------------------- */
   /* Initialize counters and flags                                  */
   /* -------------------------------------------------------------- */

   primme->stats.numOuterIterations            = 0; 
   primme->stats.numRestarts                   = 0;
   primme->stats.numMatvecs                    = 0;
   primme->stats.numPreconds                   = 0;
   primme->stats.numGlobalSum                  = 0;
   primme->stats.numBroadcast                  = 0;
   primme->stats.volumeGlobalSum               = 0;
   primme->stats.volumeBroadcast               = 0;
   primme->stats.flopsDense                    = 0;
   primme->stats.numOrthoInnerProds            = 0.0;
   primme->stats.elapsedTime                   = 0.0;
   primme->stats.timeMatvec                    = 0.0;
   primme->stats.timePrecond                   = 0.0;
   primme->stats.timeOrtho                     = 0.0;
   primme->stats.timeGlobalSum                 = 0.0;
   primme->stats.timeBroadcast                 = 0.0;
   primme->stats.timeDense                     = 0.0;
   primme->stats.estimateMinEVal               = HUGE_VAL;
   primme->stats.estimateMaxEVal               = -HUGE_VAL;
   primme->stats.estimateLargestSVal           = -HUGE_VAL;
   primme->stats.estimateBNorm                 = primme->massMatrixMatvec ? -HUGE_VAL : 1.0;
   primme->stats.estimateInvBNorm              = primme->massMatrixMatvec ? -HUGE_VAL : 1.0;
   primme->stats.maxConvTol                    = 0.0;
   primme->stats.estimateResidualError         = 0.0;
   primme->stats.lockingIssue                  = 0;

   numLocked = 0;
   LockingProblem = 0;

   blockSize = 0; 
   *numRet = 0;

   for (i=0; i<primme->numEvals; i++) perm[i] = i;

   /* -------------------------------------- */
   /* Quick return for matrix of dimension 1 */
   /* -------------------------------------- */

   if (primme->numEvals == 0) {
      primme->initSize = 0;
      *ret = 0;
      goto clean;
   }

   /* -------------------- */
   /* Initialize the basis */
   /* -------------------- */

   CHKERR(init_basis_Sprimme(V, primme->nLocal, ldV, W, ldW, BV, ldBV, evecs,
         ldevecs, Bevecs, ldBevecs, evecsHat, primme->nLocal, M, maxEvecsSize,
         Mfact, 0, ipivot, VtBV, ldVtBV, fVtBV, ldfVtBV, maxRank, &basisSize,
         &nextGuess, &numGuesses, ctx));

   /* Now initSize will store the number of converged pairs */
   primme->initSize = 0;

   /* ----------------------------------------------------------- */
   /* Dynamic method switch means we need to decide whether to    */
   /* allow inner iterations based on runtime timing measurements */
   /* ----------------------------------------------------------- */
   if (primme->dynamicMethodSwitch > 0) {
      initializeModel(&CostModel, primme);
      CostModel.MV = primme->stats.timeMatvec/primme->stats.numMatvecs;
      if (primme->numEvals < 5 ||
            primme->maxBasisSize + (primme->locking ? primme->numEvals : 0) >=
                  primme->n)
         primme->dynamicMethodSwitch = 1;   /* Start tentatively GD+k */
      else
         primme->dynamicMethodSwitch = 3;   /* Start GD+k for 1st pair */
      primme->correctionParams.maxInnerIterations = 0; 
   }

   /* ---------------------------------------------------------------------- */
   /* Outer most loop                                                        */
   /* Without locking, restarting can cause converged Ritz values to become  */
   /* unconverged. Keep performing JD iterations until they remain converged */
   /* ---------------------------------------------------------------------- */
   while (
         primme->stats.numMatvecs < primme->maxMatvecs &&
         (primme->maxOuterIterations == 0 ||
               primme->stats.numOuterIterations < primme->maxOuterIterations)) {

      if (reset > 0) PRINTF(5, "Resetting V, W and QR");

      /* Reset convergence flags. This may only reoccur without locking */

      primme->initSize = numConverged = numConvergedStored = numLocked;
      reset = 0;
      for (i=0; i<primme->maxBasisSize; i++)
         flags[i] = UNCONVERGED;

      /* Compute the initial H and solve for its eigenpairs */

      targetShiftIndex = 0;
      if (numQR) {
         int nQ = 0;
         CHKERR(update_Q_Sprimme(BV ? BV : V, primme->nLocal, ldBV, W, ldW, Q,
               ldQ, R, primme->maxBasisSize, QtQ, ldQtQ, fQtQ, ldfQtQ,
               primme->targetShifts[targetShiftIndex], 0, basisSize, &nQ, ctx));
         CHKERRM(numQR && basisSize != nQ, -1, "Not supported deficient QR");
      }

      if (H)
         CHKERR(update_projection_Sprimme(V, ldV, W, ldW, H,
               primme->maxBasisSize, primme->nLocal, 0, basisSize,
               KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx));

      if (QtV) {
        CHKERR(update_projection_Sprimme(
            Q, ldQ, V, ldV, QtV, primme->maxBasisSize, primme->nLocal, 0,
            basisSize, 0 /*unsymmetric*/, ctx));
      }

      CHKERR(solve_H_SHprimme(H, basisSize, primme->maxBasisSize,
            VtBV
                  ? &VtBV[(primme->numOrthoConst + numLocked) * ldVtBV +
                                    primme->numOrthoConst + numLocked]
                  : NULL,
            ldVtBV, R, primme->maxBasisSize, QtV, primme->maxBasisSize, QtQ,
            ldQtQ, hU, basisSize, hVecs, basisSize, hVals, hSVals,
            numConverged, ctx));

      numArbitraryVecs = 0;
      maxRecentlyConverged = availableBlockSize = blockSize = 0;
      smallestResNorm = HUGE_VAL;
      primme->stats.estimateResidualError = 0.0;
      if (!primme->locking) primme->stats.maxConvTol = 0.0;
      blockSize = 0;
      restartsSinceReset = 0;

      /* -------------------------------------------------------------- */
      /* Begin the iterative process.  Keep restarting until all of the */
      /* required eigenpairs have been found (no verification)          */
      /* -------------------------------------------------------------- */
      while (numConverged < primme->numEvals &&
             primme->stats.numMatvecs < primme->maxMatvecs &&
             (primme->maxOuterIterations == 0 ||
                   primme->stats.numOuterIterations <
                         primme->maxOuterIterations) &&
             !wholeSpace) {

         nprevhVecs = 0;
         int candidates_prepared = 0;  /* Has prepare_candidates been called? */

         /* ----------------------------------------------------------------- */
         /* Main block Davidson loop.                                         */
         /* Keep adding vectors to the basis V until the basis has reached    */
         /* maximum size or the basis plus the locked vectors span the entire */
         /* space. Once this happens, restart with a smaller basis.           */
         /* ----------------------------------------------------------------- */
         while (basisSize < primme->maxBasisSize &&
                primme->stats.numMatvecs < primme->maxMatvecs &&
                ( primme->maxOuterIterations == 0 ||
                  primme->stats.numOuterIterations < primme->maxOuterIterations) ) {

            primme->stats.numOuterIterations++;

            /* When QR are computed and there are more than one target shift, */
            /* limit blockSize and the converged values to one.               */

            if (primme->numTargetShifts > numConverged+1 && numQR) {
               availableBlockSize = 1;
               maxRecentlyConverged = numConverged-numLocked+1;
            }
            else {
               availableBlockSize = primme->maxBlockSize;
               maxRecentlyConverged = max(0, primme->numEvals-numConverged);
            }

            /* Limit blockSize to vacant vectors in the basis */

            availableBlockSize = min(availableBlockSize, primme->maxBasisSize-basisSize);

            /* Limit blockSize to remaining values to converge plus one */

            availableBlockSize = min(availableBlockSize, maxRecentlyConverged+1);

            /* Set the block with the first unconverged pairs */
            if (availableBlockSize > 0) {
               /* If locking with GD and no preconditioning is running, don't */
               /* check practically convergence now, unless whole space.      */

               int practConvCheck = 0;
               if (primme->n <= basisSize + numLocked + primme->numOrthoConst) {
                  practConvCheck = 1;
               } else if (primme->locking && !primme->massMatrixMatvec &&
                          !primme->correctionParams.precondition &&
                          primme->correctionParams.maxInnerIterations == 0) {
                  practConvCheck = -1;
               }
               CHKERR(prepare_candidates(V, ldV, W, ldW, BV, ldBV,
                     primme->nLocal, H, primme->maxBasisSize, basisSize,
                     &V[basisSize * ldV], &W[basisSize * ldW],
                     BV ? &BV[basisSize * ldBV] : NULL,
                     1 /* compute approx vectors and residuals */, hVecs,
                     basisSize, hVals, hSVals, flags, maxRecentlyConverged,
                     blockNorms, blockSize, availableBlockSize, evecs,
                     numLocked, ldevecs, Bevecs, ldBevecs, evals, resNorms,
                     targetShiftIndex, iev, &blockSize, &recentlyConverged,
                     &numArbitraryVecs, &smallestResNorm, hVecsRot,
                     primme->maxBasisSize, numConverged, basisNorms, &reset,
                     VtBV, ldVtBV, prevhVecs, nprevhVecs, primme->maxBasisSize,
                     practConvCheck, map, startTime, ctx));
               assert(recentlyConverged >= 0);
               candidates_prepared = 1;
            }
            else {
               blockSize = recentlyConverged = 0;
            }

            /* If the total number of converged pairs, including the     */
            /* recentlyConverged ones, are greater than or equal to the  */
            /* target number of eigenvalues, attempt to restart, verify  */
            /* their convergence, lock them if necessary, and return.    */
            /* For locking interior, restart and lock now any converged. */
            /* If Q, restart after an eigenpair converged to recompute   */
            /* QR with a different shift.                                */
            /* Also if it has been converged as many pairs as initial    */
            /* guesses has introduced in V, then restart and introduce   */
            /* new guesses.                                              */

            numConverged += recentlyConverged;

            /* Report iteration */

            CHKERR(monitorFun_Sprimme(hVals, basisSize, flags, iev, blockSize,
                  basisNorms, numConverged, evals, numLocked, lockedFlags,
                  resNorms, -1, -1.0, NULL, 0.0, primme_event_outer_iteration,
                  startTime, ctx));

            /* Reset touch every time an eigenpair converges */

            if (recentlyConverged > 0) touch = 0;

            if (primme->dynamicMethodSwitch > 0) {
               if (CostModel.resid_0 == -1.0L)       /* remember the very */
                  CostModel.resid_0 = blockNorms[0]; /* first residual */

               /* If some pairs converged OR we evaluate JDQMR at every step, */
               /* update convergence statistics and consider switching        */
               if (recentlyConverged > 0 || primme->dynamicMethodSwitch == 2)
               {
                  CostModel.MV =
                     primme->stats.timeMatvec/primme->stats.numMatvecs;
                  int ret0 = update_statistics(&CostModel, primme, tstart,
                        recentlyConverged, 0, numConverged, blockNorms[0]);

                  if (ret0) switch (primme->dynamicMethodSwitch) {
                     /* for few evals (dyn=1) evaluate GD+k only at restart*/
                     case 3:
                        CHKERR(switch_from_GDpk(&CostModel,ctx));
                        break;
                     case 2: case 4:
                        CHKERR(switch_from_JDQMR(&CostModel,ctx));
                  }
               }
            }

            if (numConverged >= primme->numEvals ||
                  (primme->locking && numConverged > numLocked &&
                        primme->target != primme_smallest &&
                        primme->target != primme_largest &&
                        (primme->projectionParams.projection ==
                                    primme_proj_RR ||
                              primme->target == primme_closest_geq ||
                              primme->target == primme_closest_leq)) ||
                  targetShiftIndex < 0 ||
                  (blockSize == 0 && recentlyConverged > 0) ||
                  /* NOTE: use the same condition as in restart_refined */
                  (numQR && fabs(primme->targetShifts[targetShiftIndex] -
                                  primme->targetShifts[min(
                                        primme->numTargetShifts - 1,
                                        numConverged)]) >=
                                  max(primme->aNorm,
                                        primme->stats.estimateLargestSVal)) ||
                  (numConverged >= nextGuess - primme->numOrthoConst &&
                        numGuesses > 0)) {

               break;

            }

            if (blockSize > 0) {
               /* Solve the correction equations with the new blockSize Ritz */
               /* vectors and residuals.                                     */

               /* If dynamic method switching, time the inner method     */
               if (primme->dynamicMethodSwitch > 0) {
                  tstart = primme_wTimer(); /* accumulate correction time */

               }

               CHKERR(solve_correction_Sprimme(V, ldV, W, ldW, BV ? BV : V,
                     ldBV, evecs, ldevecs, Bevecs ? Bevecs : evecs,
                     Bevecs ? ldBevecs : ldevecs, evecsHat, ldevecsHat, Mfact,
                     ipivot, evals, numLocked, numConvergedStored, hVals,
                     prevRitzVals, &numPrevRitzVals, flags, basisSize,
                     blockNorms, iev, blockSize, &touch, startTime, ctx));

               /* ------------------------------------------------------ */
               /* If dynamic method switch, accumulate inner method time */
               /* ------------------------------------------------------ */
               if (primme->dynamicMethodSwitch > 0) 
                  CostModel.time_in_inner += primme_wTimer() - tstart;

              
            } /* end of else blocksize=0 */

            /* If locking with GD and no preconditioning is running, check    */
            /* practically convergence here. The heuristic needs Vlocked'     */
            /* times r, which is going to be store in Rlocked.                */

            HSCALAR *Rlocked = NULL;
            int ldRlocked = primme->numOrthoConst + numLocked;
            int blockSize0 = blockSize;
            if (primme->locking && !primme->massMatrixMatvec &&
                  !primme->correctionParams.precondition &&
                  primme->correctionParams.maxInnerIterations == 0) {
               CHKERR(
                     Num_malloc_SHprimme(ldRlocked * blockSize, &Rlocked, ctx));
            }

            /* Orthogonalize the corrections with respect to each other and   */
            /* the current basis. If basis hasn't expanded with any new       */
            /* vector, try a random vector                                    */

            for (i=0; i < maxNumRandoms; i++) {
               int basisSizeOut;
               CHKERR(Bortho_block_Sprimme(V, ldV, VtBV, ldVtBV, fVtBV, ldfVtBV,
                     NULL, 0, basisSize, basisSize + blockSize - 1, evecs,
                     ldevecs, primme->numOrthoConst + numLocked, BV, ldBV,
                     i == 0 ? Rlocked : NULL, ldRlocked, primme->nLocal,
                     maxRank, &basisSizeOut, ctx));
               blockSize = basisSizeOut - basisSize;
               if (blockSize > 0 || availableBlockSize <= 0) break;
               PRINTF(5, "Warning: adding random vector "
                         "to the search subspace");
               CHKERR(Num_larnv_Sprimme(2, primme->iseed, primme->nLocal,
                     &V[ldV * basisSize], ctx));
               blockSize = 1;
            }

            /* The method failing in expanding the search subspace may be an  */
            /* indication that a) the orthogonalization error is to large for */
            /* distinguishing the new components on the computed correction,  */
            /* or b) there is no more pairs that satisfy the convergence      */
            /* criterion. If prepare_candidates gave a candidate, we guess    */
            /* we are in case a), and otherwise we are in case b).            */

            if (i >= maxNumRandoms) {
               if (availableBlockSize > 0 && blockSize0 <= 0 && reset == 0) {
                  wholeSpace = 1;
               } else {
                  reset = 2;
               }
               blockSize = 0;
               if (primme->locking && !primme->massMatrixMatvec &&
                     !primme->correctionParams.precondition &&
                     primme->correctionParams.maxInnerIterations == 0) {
                  CHKERR(Num_free_SHprimme(Rlocked, ctx));
               }
               break;
            }

            // If locking with GD and no preconditioning is running, check
            // practically convergence here. A pair in the block is marked as
            // practically converged if the next expression is small enough to
            // pass the convergence criterion:
            //    sqrt(|r|^2 - |Rlocked|^2*(1 + |Xx|^2)),
            // where Rlocked is V_locked'*r and Xx is V_locked'*x.
            // The term Xx is not totally clear. It appears in the perturbation
            // analysis of the eigenproblem with V_locked and x. And it is
            // necessary to pass the next test for complex:
            // tests/testi-100-LOBPCG_OrthoBasis_Window-100-primme_closest_abs-primme_proj_refined.F

            if (primme->locking && !primme->massMatrixMatvec &&
                  !primme->correctionParams.precondition &&
                  primme->correctionParams.maxInnerIterations == 0) {
               if (numLocked > 0) {
                  for (i = 0; i < blockSize0 && numConverged < primme->numEvals;
                        i++) {
                     HREAL normXx = 0.0;
                     if (primme->orth == primme_orth_explicit_I) {
                        HSCALAR *Xx = NULL;
                        CHKERR(Num_malloc_SHprimme(numLocked, &Xx, ctx));
                        CHKERR(Num_zero_matrix_SHprimme(
                              Xx, numLocked, 1, numLocked, ctx));
                        CHKERR(Num_gemv_SHprimme("N", numLocked, basisSize, 1.0,
                              &VtBV[ldVtBV * numLocked], ldVtBV,
                              &hVecs[iev[i] * basisSize], 1, 0.0, Xx, 1, ctx));
                        normXx = ABS(
                              Num_dot_SHprimme(numLocked, Xx, 1, Xx, 1, ctx));
                        CHKERR(Num_free_SHprimme(Xx, ctx));
                     }
                     HREAL normRlockedi = ABS(
                           Num_dot_SHprimme(ldRlocked, &Rlocked[ldRlocked * i],
                                 1, &Rlocked[ldRlocked * i], 1, ctx));
                     HREAL newBlockNorm =
                           sqrt(max(blockNorms[i] * blockNorms[i] -
                                          normRlockedi * (1. + normXx),
                                 0.0));
                     CHKERR(check_convergence_Sprimme(&V[(basisSize + i) * ldV],
                           ldV, 1 /* given X */, NULL, 0, 0 /* not given R */,
                           evecs, numLocked, ldevecs, Bevecs, ldBevecs, VtBV,
                           ldVtBV, 0, 1, &flags[iev[i]], &newBlockNorm,
                           &hVals[iev[i]], &reset,
                           -1 /* don't check practically convergence */, ctx));
                     basisNorms[iev[i]] = newBlockNorm;
                     if (flags[iev[i]] == CONVERGED) {
                        flags[iev[i]] = PRACTICALLY_CONVERGED;
                        numConverged++;
                        /* Report a pair was soft converged */
                        CHKERR(monitorFun_Sprimme(hVals, basisSize, flags,
                              &iev[i], 1, basisNorms, numConverged, NULL, 0,
                              NULL, NULL, -1, -1.0, NULL, 0.0,
                              primme_event_converged, startTime, ctx));
                     }
                  }
               }
               CHKERR(Num_free_SHprimme(Rlocked, ctx));

               if (numConverged > numLocked &&
                     primme->target != primme_smallest &&
                     primme->target != primme_largest &&
                     (primme->projectionParams.projection == primme_proj_RR ||
                           primme->target == primme_closest_geq ||
                           primme->target == primme_closest_leq)) {

                  break;
               }
            }

            /* Compute W = A*V for the orthogonalized corrections */

            CHKERR(matrixMatvec_Sprimme(V, primme->nLocal, ldV, W, ldW,
                     basisSize, blockSize, ctx));

            if (numQR) {
               int nQ = basisSize;
               CHKERR(update_Q_Sprimme(BV ? BV : V, primme->nLocal, ldBV, W,
                     ldW, Q, ldQ, R, primme->maxBasisSize, QtQ, ldQtQ, fQtQ,
                     ldfQtQ, primme->targetShifts[targetShiftIndex], basisSize,
                     blockSize, &nQ, ctx));

               if (basisSize + blockSize != nQ) {
                  blockSize = 0;
                  reset = 1;
                  break;
               }
            }

            /* Extend H by blockSize columns and rows and solve the */
            /* eigenproblem for the new H.                          */

            if (H)
               CHKERR(update_projection_Sprimme(V, ldV, W, ldW, H,
                     primme->maxBasisSize, primme->nLocal, basisSize, blockSize,
                     KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx));

            if (QtV)
              CHKERR(update_projection_Sprimme(
                  Q, ldQ, V, ldV, QtV, primme->maxBasisSize, primme->nLocal,
                  basisSize, blockSize, 0 /*unsymmetric*/, ctx));

            /* Copy hVecs into prevhVecs */

            CHKERR(Num_copy_matrix_SHprimme(hVecs, basisSize, basisSize,
                  basisSize, prevhVecs, primme->maxBasisSize, ctx));
            CHKERR(Num_zero_matrix_SHprimme(&prevhVecs[basisSize],
                  primme->maxBasisSize - basisSize, basisSize,
                  primme->maxBasisSize, ctx));
            nprevhVecs = basisSize;

            basisSize += blockSize;
            blockSize = 0;

            CHKERR(solve_H_SHprimme(H, basisSize, primme->maxBasisSize,
                  VtBV
                        ? &VtBV[(primme->numOrthoConst + numLocked) * ldVtBV +
                                          primme->numOrthoConst + numLocked]
                        : NULL,
                  ldVtBV, R, primme->maxBasisSize, QtV, primme->maxBasisSize,
                  QtQ, ldQtQ, hU, basisSize, hVecs, basisSize, hVals, hSVals,
                  numConverged, ctx));

            numArbitraryVecs = 0;
            candidates_prepared = 0;

            /* If QR decomposition accumulates so much error, force it to     */
            /* reset by setting targetShiftIndex to -1. We use the next       */
            /* heuristic. Note that if (s_0, u_0, y_0) is the smallest triplet*/
            /* of R, (A-tau*I)*V = Q*R, and l_0 is the Rayleigh quotient of   */
            /* V*y_0, then                                                    */
            /*    |l_0-tau| = |y_0'*V'*(A-tau*I)*V*y_0| = |y_0'*V'*Q*R*y_0| = */
            /*       = |y_0'*V'*Q*u_0*s_0| <= s_0.                            */
            /* So when |l_0-tau|-machEps*|A| > s_0, we consider to reset the  */
            /* QR factorization. machEps*|A| is the error computing l_0.      */
            /* NOTE: rarely observed |l_0-tau|-machEps*|A| slightly greater   */
            /* than s_0 after resetting. The condition restartsSinceReset > 0 */
            /* avoids infinite loop in those cases.                           */

            double eps_orth;
            CHKERR(machineEpsOrth_Sprimme(&eps_orth, ctx));
            if (primme->projectionParams.projection == primme_proj_refined &&
                  basisSize > 0 && restartsSinceReset > 1 &&
                  targetShiftIndex >= 0 &&
                  ABS((HSCALAR)primme->targetShifts[targetShiftIndex] -
                        hVals[0]) -
                              max(primme->aNorm,
                                    primme->stats.estimateLargestSVal) *
                                    eps_orth >
                        hSVals[0]) {

               reset = 2;
               PRINTF(5, "Detected some errors in QR");

               break;
            }

           /* --------------------------------------------------------------- */
         } /* while (basisSize<maxBasisSize && basisSize<n-orthoConst-numLocked)
            * --------------------------------------------------------------- */

         /* If wholeSpace, reset if the accumulated errors on the residual    */
         /* vectors are too large                                             */

         if (basisSize >= primme->n - primme->numOrthoConst - numLocked) {
            if (primme->stats.maxConvTol < primme->stats.estimateResidualError)
               reset = 1;
         }

         if (reset > 0) break;

         /* ----------------------------------------------------------------- */
         /* Restart basis will need the final coefficient vectors in hVecs    */
         /* to lock out converged vectors and to compute X and R for the next */
         /* iteration. The function prepare_vecs will make sure that hVecs    */
         /* has proper coefficient vectors.                                   */
         /* ----------------------------------------------------------------- */

         if (!candidates_prepared) {
            /* -------------------------------------------------------------- */
            /* TODO: merge all this logic in the regular main loop. After all */
            /* restarting is a regular step that also shrinks the basis.      */
            /* -------------------------------------------------------------- */

            int j,k,l,m;
            double *dummySmallestResNorm, dummyZero = 0.0;

            if (blockSize > 0) {
               availableBlockSize = blockSize;
               maxRecentlyConverged = 0;
            }

            /* When there are more than one target shift,                     */
            /* limit blockSize and the converged values to one.               */

            else if (primme->numTargetShifts > numConverged+1) {
               if (primme->locking) {
                  maxRecentlyConverged =
                     max(min(primme->numEvals, numLocked+1) - numConverged, 0);
               }
               else {
                  maxRecentlyConverged =
                     max(min(primme->numEvals, numConverged+1) - numConverged, 0);
               }
               availableBlockSize = maxRecentlyConverged;
            }

            else {
               maxRecentlyConverged = max(0, primme->numEvals-numConverged);

               /* Limit blockSize to vacant vectors in the basis */

               availableBlockSize = max(0, min(primme->maxBlockSize, primme->maxBasisSize-(numConverged-numLocked)));

               /* Limit blockSize to remaining values to converge plus one */

               availableBlockSize = min(availableBlockSize, maxRecentlyConverged+1);
            }

            /* Limit basisSize to the matrix dimension */

            availableBlockSize =
                  max(0, min(availableBlockSize,
                               primme->n - numLocked - primme->numOrthoConst));

            // Calling prepare_candidates before restarting is compulsory when:
            // - the order of the Ritz pairs depends on the residual vector
            //   norms (e.g. when using refine extraction, or target being
            //   closest to a point but on the left/right); or
            // - the search subspace is exhausted.
            // Otherwise we avoid calling prepare_candidates when the method
            // restarts all the time, like in LOBPCG.

            if (availableBlockSize <= 0 ||
                  primme->minRestartSize +
                              primme->restartingParams.maxPrevRetain +
                              availableBlockSize <
                        primme->maxBasisSize ||
                  primme->numOrthoConst + numLocked + basisSize >= primme->n ||
                  primme->projectionParams.projection != primme_proj_RR) {

               // NOTE: setting smallestResNorm to zero may overpass the inner
               // product condition (ip) and result in a smaller
               // numArbitraryVecs before restarting. This helps or is needed to
               // pass testi-100-LOBPCG_OrthoBasis_Window-100-primme_closest_abs
               // -primme_proj_refined.F

               if (primme->target == primme_closest_abs ||
                     primme->target == primme_largest_abs) {
                  dummySmallestResNorm = &dummyZero;
               } else {
                  dummySmallestResNorm = &smallestResNorm;
               }

               CHKERR(prepare_candidates(V, ldV, W, ldW, BV, ldBV,
                     primme->nLocal, H, primme->maxBasisSize, basisSize, NULL,
                     NULL, NULL,
                     0 /* Not compute approx vectors and residuals */, hVecs,
                     basisSize, hVals, hSVals, flags, maxRecentlyConverged,
                     blockNorms, blockSize, availableBlockSize, evecs,
                     numLocked, ldevecs, Bevecs, ldBevecs, evals, resNorms,
                     targetShiftIndex, iev, &blockSize, &recentlyConverged,
                     &numArbitraryVecs, dummySmallestResNorm, hVecsRot,
                     primme->maxBasisSize, numConverged, basisNorms, &reset,
                     VtBV, ldVtBV, prevhVecs, nprevhVecs, primme->maxBasisSize,
                     0, map, startTime, ctx));
               assert(recentlyConverged >= 0);

               // When QR is computed and there are more than one target shift
               // and some eigenpair has converged values to one, don't provide
               // candidate for the next iteration, because it may be the
               // closest to a different target.

               if (numQR && numConverged + recentlyConverged > numLocked &&
                     primme->numTargetShifts > numLocked + 1) {
                  blockSize = 0;
               }

               // Updated the number of converged pairs. 
               // Intentionally we include the pairs flagged SKIP_UNTIL_RESTART

               for (i = 0, numConverged = numLocked; i < basisSize; i++) {
                  if (flags[i] != UNCONVERGED &&
                        numConverged < primme->numEvals &&
                        (i < primme->numEvals - numLocked
                              /* Refined and prepare_vecs may not completely */
                              // order pairs considering closest_leq/geq; so we
                              // find converged pairs beyond the first remaining
                              // pairs to converge.
                              || primme->target == primme_closest_geq ||
                              primme->target == primme_closest_leq)) {

                     numConverged++;
                  }
               }

               /* Move the converged pairs and the ones in iev at the beginning
                */

               int *iwork; /* new permutation of pairs */
               CHKERR(Num_malloc_iprimme(basisSize, &iwork, ctx));
               for (i = k = l = m = 0; i < basisSize; i++) {
                  int iIsInIev = 0;
                  for (j = 0; j < blockSize; j++)
                     if (iev[j] == i) iIsInIev = 1;
                  if ((flags[i] != UNCONVERGED &&
                            m++ < numConverged - numLocked) ||
                        iIsInIev) {
                     iwork[k++] = i;
                  } else {
                     iwork[numConverged - numLocked + blockSize + l++] = i;
                  }
               }
               CHKERR(KIND(permute_vecs_RHprimme, permute_vecs_SHprimme)(
                     hVals, 1, basisSize, 1, iwork, ctx));
               CHKERR(permute_vecs_SHprimme(
                     hVecs, basisSize, basisSize, basisSize, iwork, ctx));
               CHKERR(permute_vecs_iprimme(flags, basisSize, iwork, ctx));
               if (hVecsRot) {
                  Num_zero_matrix_SHprimme(
                        &hVecsRot[numArbitraryVecs * primme->maxBasisSize],
                        primme->maxBasisSize, basisSize - numArbitraryVecs,
                        primme->maxBasisSize, ctx);
                  for (i = numArbitraryVecs; i < basisSize; i++)
                     hVecsRot[primme->maxBasisSize * i + i] = 1.0;
                  CHKERR(permute_vecs_SHprimme(hVecsRot, basisSize, basisSize,
                        primme->maxBasisSize, iwork, ctx));
                  for (i = j = 0; i < basisSize; i++)
                     if (iwork[i] != i) j = i + 1;
                  numArbitraryVecs = max(numArbitraryVecs, j);
               }
               CHKERR(Num_free_iprimme(iwork, ctx));
            } else {
               blockSize = availableBlockSize;
               i = 0;
               for (i = 0; i < blockSize; i++) iev[i] = i;
               if (ctx.procID == 0) {
                  CHKERR(map_vecs_SHprimme(prevhVecs, basisSize, nprevhVecs,
                        primme->maxBasisSize, hVecs, 0, basisSize, basisSize,
                        map, ctx));
               }
               CHKERR(broadcast_iprimme(map, basisSize, ctx));
            }
         }

         if (reset > 0) break;

         /* Permute prevhVecs respect to map */

         CHKERR(permute_vecs_SHprimme(prevhVecs, basisSize, nprevhVecs,
               primme->maxBasisSize, map, ctx));

         /* ------------------ */
         /* Restart the basis  */
         /* ------------------ */

         assert(ldV == ldW); /* this function assumes ldV == ldW */
         CHKERR(restart_Sprimme(V, W, BV, primme->nLocal, basisSize, ldV, hVals,
               hSVals, flags, iev, &blockSize, blockNorms, evecs, ldevecs,
               Bevecs, ldBevecs, perm, evals, resNorms, evecsHat,
               primme->nLocal, M, maxEvecsSize, Mfact, 0, ipivot, &numConverged,
               &numLocked, lockedFlags, &numConvergedStored, prevhVecs,
               nprevhVecs, primme->maxBasisSize, numGuesses, prevRitzVals,
               &numPrevRitzVals, H, primme->maxBasisSize, VtBV, ldVtBV, fVtBV,
               ldfVtBV, Q, ldQ, R, primme->maxBasisSize, QtV,
               primme->maxBasisSize, QtQ, ldQtQ, fQtQ, ldfQtQ, hU, basisSize, 0,
               hVecs, basisSize, 0, &basisSize, &targetShiftIndex,
               &numArbitraryVecs, hVecsRot, primme->maxBasisSize,
               &restartsSinceReset, startTime, ctx));
         restartsSinceReset++;

         /* If there are any initial guesses remaining, then copy it */
         /* into the basis.                                          */

         if (numGuesses > 0) {
            /* Try to keep minRestartSize guesses in the search subspace */

            int numNew = max(0, min(
                     primme->minRestartSize + numConverged 
                     - (nextGuess - primme->numOrthoConst), numGuesses));

            /* Don't make the resulting basis size larger than maxBasisSize */

            numNew = max(0,
                  min(basisSize+numNew, primme->maxBasisSize) - basisSize);

            /* Don't increase basis size beyond matrix dimension */

            numNew = max(0,
                  min(basisSize+numNew+primme->numOrthoConst+numLocked,
                     primme->n)
                  - primme->numOrthoConst - numLocked - basisSize);

            CHKERR(Num_copy_matrix_Sprimme(&evecs[nextGuess * ldevecs],
                  primme->nLocal, numNew, ldevecs, &V[basisSize * ldV], ldV,
                  ctx));

            nextGuess += numNew;
            numGuesses -= numNew;

            int basisSizeOut;
            CHKERR(Bortho_block_Sprimme(V, ldV, VtBV, ldVtBV, fVtBV, ldfVtBV,
                  NULL, 0, basisSize, basisSize + numNew - 1, evecs, ldevecs,
                  numLocked + primme->numOrthoConst, BV, ldBV, NULL, 0,
                  primme->nLocal, maxRank, &basisSizeOut, ctx));
            numNew = basisSizeOut - basisSize;

            /* Compute W = A*V for the orthogonalized corrections */

            CHKERR(matrixMatvec_Sprimme(V, primme->nLocal, ldV, W, ldW,
                     basisSize, numNew, ctx));

            if (numQR) {
               int nQ = basisSize;
               CHKERR(update_Q_Sprimme(V, primme->nLocal, ldV, W, ldW, Q, ldQ,
                     R, primme->maxBasisSize, QtQ, ldQtQ, fQtQ, ldfQtQ,
                     primme->targetShifts[targetShiftIndex], basisSize, numNew,
                     &nQ, ctx));
               CHKERRM(basisSize + numNew != nQ, -1,
                     "Not supported deficient QR");
            }

            /* Extend H by numNew columns and rows and solve the */
            /* eigenproblem for the new H.                       */

            if (H)
               CHKERR(update_projection_Sprimme(V, ldV, W, ldW, H,
                     primme->maxBasisSize, primme->nLocal, basisSize, numNew,
                     KIND(1 /*symmetric*/, 0 /*unsymmetric */), ctx));
            if (QtV) CHKERR(update_projection_Sprimme(Q, ldQ, V, ldV, QtV,
                     primme->maxBasisSize, primme->nLocal, basisSize, numNew,
                     0/*unsymmetric*/, ctx));
            basisSize += numNew;
            CHKERR(solve_H_SHprimme(H, basisSize, primme->maxBasisSize,
                  VtBV
                        ? &VtBV[(primme->numOrthoConst + numLocked) * ldVtBV +
                                          primme->numOrthoConst + numLocked]
                        : NULL,
                  ldVtBV, R, primme->maxBasisSize, QtV, primme->maxBasisSize,
                  QtQ, ldQtQ, hU, basisSize, hVecs, basisSize, hVals,
                  hSVals, numConverged, ctx));
         }
 
         primme->stats.numRestarts++;

         primme->initSize = numConverged;

         /* ------------------------------------------------------------- */
         /* If dynamic method switching == 1, update model parameters and */
         /* evaluate whether to switch from GD+k to JDQMR. This is after  */
         /* restart. GD+k is also evaluated if a pair converges.          */
         /* ------------------------------------------------------------- */
         if (primme->dynamicMethodSwitch == 1 ) {
            tstart = primme_wTimer();
            CostModel.MV = primme->stats.timeMatvec/primme->stats.numMatvecs;
            update_statistics(&CostModel, primme, tstart, 0, 1,
               numConverged, blockNorms[0]); 
            CHKERR(switch_from_GDpk(&CostModel, ctx));
         } /* ---------------------------------------------------------- */

         /* Reset map */

         for (i = 0; i < primme->maxBasisSize; i++) map[i] = i;

         /* ----------------------------------------------------------- */
      } /* while ((numConverged < primme->numEvals)  (restarting loop)
         * ----------------------------------------------------------- */

      if (reset > 0) continue;

      /* ------------------------------------------------------------ */
      /* If locking is enabled, check to make sure the required       */
      /* number of eigenvalues have been computed, else make sure the */
      /* residual norms of the converged Ritz vectors have remained   */
      /* converged by calling verify_norms.                           */
      /* ------------------------------------------------------------ */

      if (primme->locking) {

         // If maximum iterations/matvecs reached, return as many candidates as
         // possible as converged

         CHKERR(copy_back_candidates(V, ldV, W, ldW, BV, ldBV, primme->nLocal,
               H, primme->maxBasisSize, basisSize, hVecs, primme->maxBasisSize,
               hVals, evecs, ldevecs, evals, resNorms, targetShiftIndex,
               numConverged,
               VtBV ? &VtBV[(primme->numOrthoConst + numLocked) * ldVtBV +
                            primme->numOrthoConst + numLocked]
                    : NULL,
               ldVtBV, numRet, ctx));

         /* if dynamic method, give method recommendation for future runs */
         if (primme->dynamicMethodSwitch > 0 ) {
            if (CostModel.accum_jdq_gdk < 0.96) 
               primme->dynamicMethodSwitch = -2;  /* Use JDQMR_ETol */
            else if (CostModel.accum_jdq_gdk > 1.04) 
               primme->dynamicMethodSwitch = -1;  /* Use GD+k */
            else
               primme->dynamicMethodSwitch = -3;  /* Close call. Use dynamic */
         }

         /* Return flag showing if there has been a locking problem */
         primme->stats.lockingIssue = LockingProblem;

         /* If all of the target eigenvalues have been computed, */
         /* then return success, else return with a failure.     */
 
         if (numConverged == primme->numEvals || wholeSpace) {
            *ret = 0;
         }
         else {
            *ret = PRIMME_MAIN_ITER_FAILURE;
         }
         goto clean;

      }
      else {      /* no locking. Verify that everything is converged  */

         /* Determine if the maximum number of matvecs or iterations has been
          * reached */

         if (primme->stats.numMatvecs < primme->maxMatvecs &&
               (primme->maxOuterIterations == 0 ||
                     primme->stats.numOuterIterations <
                           primme->maxOuterIterations))
            restartLimitReached = 0;
         else
            restartLimitReached = 1;

         /* ---------------------------------------------------------- */
         /* The norms of the converged Ritz vectors must be recomputed */
         /* before return.  This is because the restarting may cause   */
         /* some converged Ritz vectors to become slightly unconverged.*/
         /* If some have become unconverged, then further iterations   */
         /* must be performed to force these approximations back to a  */
         /* converged state.                                           */
         /* ---------------------------------------------------------- */

         CHKERR(verify_norms(V, ldV, W, ldW, BV, ldBV, hVals,
               restartLimitReached ? primme->numEvals : numConverged, resNorms,
               flags, &numConverged, ctx));

         /* ---------------------------------------------------------- */
         /* If the convergence limit is reached or the target vectors  */
         /* have remained converged, then copy the current Ritz values */
         /* and vectors to the output arrays and return, else continue */
         /* iterating.                                                 */
         /* ---------------------------------------------------------- */

         if (restartLimitReached || numConverged >= primme->numEvals ||
               wholeSpace) {
            for (i=0; i < primme->numEvals; i++) {
               evals[i] = hVals[i];
               perm[i] = i;
            }

            CHKERR(Num_copy_matrix_Sprimme(V, primme->nLocal, primme->numEvals,
                  ldV, &evecs[ldevecs * primme->numOrthoConst], ldevecs, ctx));
            *numRet = primme->numEvals;

            /* The target values all remained converged, then return */
            /* successfully, else return with a failure code.        */
            /* Return also the number of actually converged pairs    */
 
            primme->initSize = numConverged;

            /* if dynamic method, give method recommendation for future runs */
            if (primme->dynamicMethodSwitch > 0 ) {
               if (CostModel.accum_jdq_gdk < 0.96) 
                  primme->dynamicMethodSwitch = -2;  /* Use JDQMR_ETol */
               else if (CostModel.accum_jdq_gdk > 1.04) 
                  primme->dynamicMethodSwitch = -1;  /* Use GD+k */
               else
                  primme->dynamicMethodSwitch = -3;  /* Close call.Use dynamic*/
            }

            if (numConverged >= primme->numEvals) {
               *ret = 0;
            }
            else {
               *ret = PRIMME_MAIN_ITER_FAILURE;
            }
            goto clean;

         } else {
            /* ------------------------------------------------------------ */
            /* Reorthogonalize the basis, recompute W=AV, and continue the  */
            /* outer while loop, resolving the epairs. Slow, but robust!    */
            /* ------------------------------------------------------------ */
            CHKERR(Bortho_block_Sprimme(V, ldV, VtBV, ldVtBV, fVtBV, ldfVtBV,
                  NULL, 0, 0, basisSize - 1, evecs, ldevecs,
                  primme->numOrthoConst, BV, ldBV, NULL, 0, primme->nLocal,
                  maxRank, &basisSize, ctx));
            CHKERR(matrixMatvec_Sprimme(V, primme->nLocal, ldV, W, ldW, 0,
                     basisSize, ctx));

            PRINTF(2,
                  "Verifying before return: Some vectors are unconverged");

            restartsSinceReset = 0;
            reset = 0;
            primme->stats.estimateResidualError = 0.0;
            numConverged = 0;

         } 

        /* ------------------------------------------------------------ */
      } /* End of non locking
         * ------------------------------------------------------------ */

     /* -------------------------------------------------------------- */
   } /* while (!converged)  Outer verification loop
      * -------------------------------------------------------------- */

clean:
   if (primme->aNorm <= 0.0L) {
      primme->aNorm =
            primme->stats.estimateLargestSVal / primme->stats.estimateInvBNorm;
   }

   /*----------------------------------------------------------------------*/
   /* If locking is engaged, the converged Ritz vectors are stored in the  */
   /* order they converged.  They must then be permuted so that they       */
   /* correspond to the sorted Ritz values in evals.                       */
   /*----------------------------------------------------------------------*/

   CHKERR(permute_vecs_Sprimme(&evecs[primme->numOrthoConst * primme->ldevecs],
            primme->nLocal, primme->initSize, primme->ldevecs, perm,
            ctx));

   if (primme->massMatrixMatvec) {
      CHKERR(Num_free_Sprimme(BV, ctx));
   }
   CHKERR(Num_free_Sprimme(V, ctx));
   CHKERR(Num_free_Sprimme(W, ctx));
   if (numQR > 0) {
      CHKERR(Num_free_Sprimme(Q, ctx));
      CHKERR(Num_free_SHprimme(R, ctx));
      CHKERR(Num_free_SHprimme(hU, ctx));
   }
   if (primme->projectionParams.projection == primme_proj_harmonic) {
      CHKERR(Num_free_SHprimme(QtV, ctx));
   }
   CHKERR(Num_free_SHprimme(H, ctx));
   CHKERR(Num_free_SHprimme(hVecs, ctx));
   if (primme->orth == primme_orth_explicit_I) {
      CHKERR(Num_free_SHprimme(VtBV, ctx));
      CHKERR(Num_free_SHprimme(fVtBV, ctx));
      CHKERR(Num_free_SHprimme(QtQ, ctx));
      CHKERR(Num_free_SHprimme(fQtQ, ctx));
   }
   if (primme->projectionParams.projection == primme_proj_refined ||
         primme->projectionParams.projection == primme_proj_harmonic) {
      CHKERR(Num_free_SHprimme(hVecsRot, ctx));
   }

   if (primme->correctionParams.precondition &&
         primme->correctionParams.maxInnerIterations != 0 &&
         primme->correctionParams.projectors.RightQ &&
         primme->correctionParams.projectors.SkewQ) {
      CHKERR(Num_free_Sprimme(evecsHat, ctx));
      CHKERR(Num_free_SHprimme(M, ctx));
      CHKERR(Num_free_SHprimme(Mfact, ctx));
   }
   if (primme->massMatrixMatvec && primme->locking) {
      CHKERR(Num_free_Sprimme(Bevecs, ctx));
   }
   CHKERR(Num_free_SHprimme(prevhVecs, ctx));

   CHKERR(KIND(Num_free_RHprimme, Num_free_SHprimme)(hVals, ctx));
   if (numQR > 0) { CHKERR(Num_free_RHprimme(hSVals, ctx)); }
   CHKERR(KIND(Num_free_RHprimme, Num_free_SHprimme)(prevRitzVals, ctx));
   CHKERR(Num_free_RHprimme(blockNorms, ctx));
   CHKERR(Num_free_RHprimme(basisNorms, ctx));

   /* Integer objects */

   CHKERR(Num_free_iprimme(perm, ctx));
   if (primme->locking) { CHKERR(Num_free_iprimme(lockedFlags, ctx)); }
   CHKERR(Num_free_iprimme(flags, ctx));
   CHKERR(Num_free_iprimme(map, ctx));
   CHKERR(Num_free_iprimme(iev, ctx));
   CHKERR(Num_free_iprimme(ipivot, ctx));

   return 0;
}

/******************************************************************************
          Some basic functions within the scope of main_iter
*******************************************************************************/

/*******************************************************************************
 * Subroutine prepare_candidates - This subroutine puts into the block the first
 *    unconverged Ritz pairs, up to maxBlockSize. If needed, compute residuals
 *    and rearrange the coefficient vectors in hVecs.
 * 
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V              The orthonormal basis
 * W              A*V
 * BV             B*V
 * nLocal         Local length of vectors in the basis
 * basisSize      Size of the basis V and W
 * ldV            The leading dimension of V and X
 * ldW            The leading dimension of W and R
 * hVecs          The projected vectors
 * ldhVecs        The leading dimension of hVecs
 * hVals          The Ritz values
 * maxBasisSize   maximum allowed size of the basis
 * numSoftLocked  Number of vectors that have converged (not updated here)
 * remainedEvals  Remained number of eigenpairs that the user wants computed
 * blockNormsSize Number of already computed residuals
 * maxBlockSize   maximum allowed size of the block
 * evecs          Converged eigenvectors
 * evecsSize      The size of evecs
 * numLocked      The number of vectors currently locked (if locking)
 * numConverged   Number of converged pairs (soft+hard locked)
 * prevhVecs      hVecs from previous iteration
 * nprevhVecs     number of columns in prevhVecs
 * ldprevhVecs    The leading dimension in prevhVecs
 * practConvChecking  flag to disable or enforce checking for practically convergence
 * primme         Structure containing various solver parameters
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * X             The eigenvectors put in the block
 * R             The residual vectors put in the block
 * flags         Array indicating which eigenvectors have converged     
 * map           The indices of closest vectors in previous iteration
 * iev           indicates which eigenvalue each block vector corresponds to
 * blockNorms    Residual norms of the Ritz vectors being computed during the
 *               current iteration
 * blockSize     Dimension of the block
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * recentlyConverged  Number of pairs converged
 * reset         flag to reset V and W in the next restart
 * 
 ******************************************************************************/

STATIC int prepare_candidates(SCALAR *V, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, SCALAR *BV, PRIMME_INT ldBV, PRIMME_INT nLocal,
      HSCALAR *H, int ldH, int basisSize, SCALAR *X, SCALAR *R, SCALAR *BX,
      int computeXR, HSCALAR *hVecs, int ldhVecs, HEVAL *hVals, HREAL *hSVals,
      int *flags, int remainedEvals, HREAL *blockNorms, int blockNormsSize,
      int maxBlockSize, SCALAR *evecs, int numLocked, PRIMME_INT ldevecs,
      SCALAR *Bevecs, PRIMME_INT ldBevecs, HEVAL *evals, HREAL *resNorms,
      int targetShiftIndex, int *iev, int *blockSize, int *recentlyConverged,
      int *numArbitraryVecs, double *smallestResNorm, HSCALAR *hVecsRot,
      int ldhVecsRot, int numConverged, HREAL *basisNorms, int *reset,
      HSCALAR *VtBV, int ldVtBV, HSCALAR *prevhVecs, int nprevhVecs,
      int ldprevhVecs, int practConvChecking, int *map, double startTime,
      primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, blki;            /* loop variables */
   HEVAL *hValsBlock;      /* contiguous copy of the hVals to be tested */
   HSCALAR *hVecsBlock;    /* contiguous copy of the hVecs columns to be tested */     
   int *flagsBlock;        /* contiguous copy of the flags to be tested */
   HEVAL *hValsBlock0;     /* workspace for hValsBlock */
   HSCALAR *hVecsBlock0;   /* workspace for hVecsBlock */
   HREAL *XNorms=NULL;     /* 2-norm of V*hVecs */
   int lasti;              /* last tested pair */

   *blockSize = 0;
   CHKERR(KIND(Num_malloc_RHprimme, Num_malloc_SHprimme)(
         maxBlockSize, &hValsBlock0, ctx));
   CHKERR(Num_malloc_SHprimme(ldhVecs*maxBlockSize, &hVecsBlock0, ctx));
   CHKERR(Num_malloc_iprimme(maxBlockSize, &flagsBlock, ctx));
   if (primme->massMatrixMatvec) {
      CHKERR(Num_malloc_RHprimme(maxBlockSize, &XNorms, ctx));
   }
   lasti = -1;

   /* Pack hVals for already computed residual pairs */

   hValsBlock = KIND(Num_compact_vecs_RHprimme, Num_compact_vecs_SHprimme)(
         hVals, 1, blockNormsSize, 1, &iev[*blockSize], hValsBlock0, 1,
         1 /* avoid copy */, ctx);

   /* If some residual norms have already been computed, set the minimum   */
   /* of them as the smallest residual norm. If not, use the value from    */
   /* previous iteration.                                                  */

   if (blockNormsSize > 0) {
      for (*smallestResNorm = HUGE_VAL, i=0; i < blockNormsSize; i++) {
         *smallestResNorm = min(*smallestResNorm, blockNorms[i]);
      }
   }

   /* Assign to each eigenpair a pair from previous iteration that is close in
    * angle */

   if (ctx.procID == 0) {
      CHKERR(map_vecs_SHprimme(prevhVecs, basisSize, nprevhVecs, ldprevhVecs,
            hVecs, 0, basisSize, ldhVecs, map, ctx));
   }
   CHKERR(broadcast_iprimme(map, basisSize, ctx));

   /* Reorder the flags from previous iteration following map */

   CHKERR(permute_vecs_iprimme(flags, basisSize, map, ctx));

   *recentlyConverged = 0;
   while (1) {
      /* Recompute flags in iev(*blockSize:*blockSize+blockNormsize) */
      for (i = *blockSize; i < *blockSize+blockNormsSize; i++)
         flagsBlock[i-*blockSize] = flags[iev[i]];
      CHKERR(check_convergence_Sprimme(X ? &X[(*blockSize) * ldV] : NULL, ldV,
            computeXR, R ? &R[(*blockSize) * ldW] : NULL, ldW, computeXR, evecs,
            numLocked, ldevecs, Bevecs, ldBevecs, VtBV, ldVtBV, 0,
            blockNormsSize, flagsBlock,
            blockNorms ? &blockNorms[*blockSize] : NULL, hValsBlock, reset,
            practConvChecking, ctx));

      /* Compact blockNorms, X and R for the unconverged pairs in    */
      /* iev(*blockSize:*blockSize+blockNormsize). Do the proper     */
      /* actions for converged pairs.                                */

      for (blki=*blockSize, i=0; i < blockNormsSize && *blockSize < maxBlockSize; i++, blki++) {
         /* Write back flags and residual norms */
         flags[iev[blki]] = flagsBlock[i];
         basisNorms[iev[blki]] = blockNorms[blki];

         /* Ignore some cases */
#ifdef USE_HERMITIAN
         double targetShift = primme->targetShifts
                                    ? primme->targetShifts[targetShiftIndex]
                                    : 0.0;
         if ((primme->target == primme_closest_leq
                  && hVals[iev[blki]]-blockNorms[blki] > targetShift) ||
               (primme->target == primme_closest_geq
                && hVals[iev[blki]]+blockNorms[blki] < targetShift)) {
         } else
#endif
               if (flagsBlock[i] != UNCONVERGED &&
                     *recentlyConverged < remainedEvals &&
                     (iev[blki] < primme->numEvals - numLocked
                           /* Refined and prepare_vecs may not completely    */
                           /* order pairs considering closest_leq/geq; so we */
                           /* find converged pairs beyond the first remaining*/
                           /* pairs to converge.                             */
                           || primme->target == primme_closest_geq ||
                           primme->target == primme_closest_leq)) {

            /* Write the current Ritz value in evals and the residual in resNorms;  */
            /* it will be checked by restart routine later.                         */
            /* Also print the converged eigenvalue.                                 */

            if (!primme->locking) {
               evals[iev[blki]] = hVals[iev[blki]];
               resNorms[iev[blki]] = blockNorms[blki];
               if (flagsBlock[i] == CONVERGED)
                  primme->stats.maxConvTol =
                        max(primme->stats.maxConvTol, blockNorms[blki]);
            }

            /* Count the new solution */
            (*recentlyConverged)++;

            /* Reset smallestResNorm if some pair converged */
            if (*blockSize == 0) {
               *smallestResNorm = HUGE_VAL;
            }

            /* Update maxBlockSize */
            maxBlockSize =
                min(maxBlockSize,
                    primme->numEvals + 1 - (*recentlyConverged) - numConverged);

            /* Report a pair was soft converged */
            CHKERR(monitorFun_Sprimme(hVals, basisSize, flags, &iev[blki], 1,
                  basisNorms, numConverged + *recentlyConverged, NULL, 0, NULL,
                  NULL, -1, -1.0, NULL, 0.0, primme_event_converged, startTime,
                  ctx));
         }
         else if (flagsBlock[i] == UNCONVERGED) {
            /* Update the smallest residual norm */
            if (*blockSize == 0) {
               *smallestResNorm = HUGE_VAL;
            }
            *smallestResNorm = min(*smallestResNorm, blockNorms[blki]);

            blockNorms[*blockSize] = blockNorms[blki];
            iev[*blockSize] = iev[blki];
            if (computeXR) {
               CHKERR(Num_copy_matrix_Sprimme(&X[blki * ldV], nLocal, 1, ldV,
                     &X[(*blockSize) * ldV], ldV, ctx));
               CHKERR(Num_copy_matrix_Sprimme(&R[blki * ldW], nLocal, 1, ldW,
                     &R[(*blockSize) * ldW], ldW, ctx));
               if (BX) {
                  CHKERR(Num_copy_matrix_Sprimme(&BX[blki * ldBV], nLocal, 1,
                        ldBV, &BX[(*blockSize) * ldBV], ldBV, ctx));
               }
            }
            (*blockSize)++;
         }

         lasti = iev[blki];
      }

      /* Generate well conditioned coefficient vectors; start from the last   */
      /* position visited (variable i)                                        */

      blki = *blockSize;
      CHKERR(prepare_vecs_SHprimme(basisSize, lasti + 1, maxBlockSize - blki, H,
            ldH, hVals, hSVals, hVecs, ldhVecs, targetShiftIndex,
            numArbitraryVecs, *smallestResNorm, flags, 1, hVecsRot, ldhVecsRot,
            ctx));

      /* Find next candidates, starting from iev(*blockSize)+1 */

      for (i=lasti+1; i<basisSize && blki < maxBlockSize; i++) {
         if (flags[i] == UNCONVERGED) iev[blki++] = i;
      }

      /* If no new candidates or all required solutions converged yet, go out */

      if (blki == *blockSize || *recentlyConverged >= remainedEvals) break;
      blockNormsSize = blki - *blockSize;

      /* Pack hVals & hVecs */

      hValsBlock = KIND(Num_compact_vecs_RHprimme, Num_compact_vecs_SHprimme)(
            hVals, 1, blockNormsSize, 1, &iev[*blockSize], hValsBlock0, 1,
            1 /* avoid copy */, ctx);
      hVecsBlock = Num_compact_vecs_SHprimme(
          hVecs, basisSize, blockNormsSize, ldhVecs, &iev[*blockSize],
          hVecsBlock0, ldhVecs, 1 /* avoid copy */, ctx);

      /* Compute X, R and residual norms for the next candidates */
      /* X(basisSize:) = V*hVecs(*blockSize:*blockSize+blockNormsize) */
      /* BX(basisSize:) = BV*hVecs(*blockSize:*blockSize+blockNormsize) */
      /* R(basisSize:) = W*hVecs(*blockSize:*blockSize+blockNormsize) -
       * BX(basisSize:)*diag(hVals) */
      /* blockNorms(basisSize:) = norms(R(basisSize:)) */

      assert(ldV == ldW); /* This functions only works in this way */
      CHKERR(Num_update_VWXR_Sprimme(V, W, BV, nLocal, basisSize, ldV,
               hVecsBlock, basisSize, ldhVecs, hValsBlock,
               X?&X[(*blockSize)*ldV]:NULL, 0, computeXR?blockNormsSize:0, ldV,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               R?&R[(*blockSize)*ldV]:NULL, 0, computeXR?blockNormsSize:0, ldV, computeXR?&blockNorms[*blockSize]:NULL,
               BX?&BX[(*blockSize)*ldV]:NULL, 0, BX?blockNormsSize:0, ldV,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               &blockNorms[*blockSize], 0, !computeXR?blockNormsSize:0,
               NULL, 0, 0,
               NULL, 0, 0,
               XNorms?&XNorms[*blockSize]:NULL, 0, primme->massMatrixMatvec?blockNormsSize:0,
               ctx));

      /* Don't trust residual norm smaller than the error in the residual norm */

      for (i = *blockSize; i < blockNormsSize; i++) {
         blockNorms[i] = max(blockNorms[i], primme->stats.estimateResidualError);
      }


      /* Update estimation of |inv(B)| */

      if (XNorms) {
         for (i=0; i<blockNormsSize; i++) {
            primme->stats.estimateBNorm = max(
                  primme->stats.estimateBNorm, 1.0 / XNorms[*blockSize + i]);
            primme->stats.estimateInvBNorm = max(
                  primme->stats.estimateInvBNorm, XNorms[*blockSize + i]);
         }
      }
   }

   CHKERR(KIND(Num_free_RHprimme, Num_free_SHprimme)(hValsBlock0, ctx));
   CHKERR(Num_free_SHprimme(hVecsBlock0, ctx));
   CHKERR(Num_free_iprimme(flagsBlock, ctx));
   CHKERR(Num_free_RHprimme(XNorms, ctx));

   return 0;
}

/*******************************************************************************
 * Subroutine copy_back_candidates - This subroutine copy into evecs, evals
 *    and resNorms unconverged Ritz pairs, up to numEvals. The uninitalized
 *    returned pairs are marked setting its resNorms to -1.
 * 
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V              The orthonormal basis
 * W              A*V
 * BV             B*V
 * nLocal         Local length of vectors in the basis
 * basisSize      Size of the basis V and W
 * ldV            The leading dimension of V and X
 * ldW            The leading dimension of W and R
 * hVecs          The projected vectors
 * ldhVecs        The leading dimension of hVecs
 * hVals          The Ritz values
 * primme         Structure containing various solver parameters
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * evecs          Converged eigenvectors
 * evals          Converged values
 * resNorms       Converged residual norms
 * numRet         number of returned pairs
 *
 * RETURN
 * ------
 * error code
 * 
 ******************************************************************************/

STATIC int copy_back_candidates(SCALAR *V, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, SCALAR *BV, PRIMME_INT ldBV, PRIMME_INT nLocal,
      HSCALAR *H, int ldH, int basisSize, HSCALAR *hVecs, int ldhVecs,
      HEVAL *hVals, SCALAR *evecs, PRIMME_INT ldevecs, HEVAL *evals,
      HREAL *resNorms, int targetShiftIndex, int numConverged, HSCALAR *VtBV,
      int ldVtBV, int *numRet, primme_context ctx) {

   (void)ldBV;
#ifndef USE_HERMITIAN
   (void)targetShiftIndex;
#endif

   primme_params *primme = ctx.primme;

   // Quick exit

   if (numConverged >= primme->numEvals || basisSize <= 0) {
      return 0;
   }

   evecs = &evecs[ldevecs * primme->numOrthoConst];

   // To avoid calling prepare_vecs Use Rayleigh-Ritz instead */

   if (primme->projectionParams.projection != primme_proj_RR &&
         primme->numTargetShifts > numConverged) {
      CHKERR(solve_H_RR_SHprimme(H, ldH, VtBV, ldVtBV, hVecs, ldhVecs, hVals,
            basisSize, numConverged, ctx));
   }

   int i = 0;
   while (i < basisSize && numConverged < primme->numEvals) {

      int blockSize =
            max(0, min(primme->numEvals - numConverged, basisSize - i));

      /* Append the new blockSize vectors to evecs and compute residual norms */
      /* evecs(numConverged+(1:blockSize)) = V*hVecs(i+(1:blockSize)) */
      /* blockNorms(1:blockSize) = residual norms */

      assert(ldV == ldW &&
             (BV == NULL ||
                   ldBV == ldV)); /* This functions only works in this way */
      CHKERR(Num_update_VWXR_Sprimme(V, W, BV, nLocal, basisSize, ldV,
               &hVecs[ldhVecs*i], basisSize, ldhVecs, &hVals[i],
               &evecs[ldevecs * numConverged], 0, blockSize, ldevecs,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0, NULL,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               &resNorms[numConverged], 0, blockSize,
               NULL, 0, 0,
               NULL, 0, 0,
               NULL, 0, 0,
               ctx));

      int numConverged0 = numConverged;
      int blki;
      for (blki = 0; blki < blockSize; blki++, i++) {
#ifdef USE_HERMITIAN
         double targetShift = primme->targetShifts
            ? primme->targetShifts[targetShiftIndex]
                                    : 0.0;
         if ((primme->target == primme_closest_leq &&
                   hVals[i] - resNorms[numConverged0 + blki] > targetShift) ||
               (primme->target == primme_closest_geq &&
                     hVals[i] + resNorms[numConverged0 + blki] < targetShift)) {
            /* Ignore some cases */
         } else
#endif
         {
            evals[numConverged] = hVals[i];
            resNorms[numConverged] = resNorms[numConverged0 + blki];
            CHKERR(Num_copy_matrix_Sprimme(
                  &evecs[(numConverged0 + blki) * ldevecs], nLocal, 1, ldevecs,
                  &evecs[numConverged * ldevecs], ldevecs, ctx));
            numConverged++;
         }
      }
   }

   for (i = numConverged; i < primme->numEvals; i++) resNorms[i] = -1;
   *numRet = numConverged;

   return 0;
}

/*******************************************************************************
 * Function verify_norms - This subroutine computes the residual norms of the 
 *    target eigenvectors before the Davidson-type main iteration terminates. 
 *    This is done to insure that the eigenvectors have remained converged. 
 *    If any have become unconverged, this routine performs restarting before 
 *    the eigensolver iteration is once again performed.  
 *
 * Note: This routine assumes it is called immediately after a call to the 
 *       restart subroutine.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V            The orthonormal basis
 *
 * W            A*V
 *
 * hVals        The eigenvalues of V'*A*V
 *
 * basisSize    Size of the basis V
 *
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * resNorms      Residual norms for each eigenpair
 *
 * flag          Indicates which Ritz pairs have converged
 *
 * numConverged  Return the number of consecutive converged pairs
 *
 ******************************************************************************/

STATIC int verify_norms(SCALAR *V, PRIMME_INT ldV, SCALAR *W, PRIMME_INT ldW,
      SCALAR *BV, PRIMME_INT ldBV, HEVAL *hVals, int basisSize, HREAL *resNorms,
      int *flags, int *numConverged, primme_context ctx) {

   int i;         /* Loop variable                                     */

   /* Compute the residual vectors */

   for (i=0; i < basisSize; i++) {
      Num_axpy_Sprimme(ctx.primme->nLocal, -hVals[i],
            BV ? &BV[ldBV * i] : &V[ldV * i], 1, &W[ldW * i], 1, ctx);
      resNorms[i] = REAL_PART(Num_dot_Sprimme(
            ctx.primme->nLocal, &W[ldW * i], 1, &W[ldW * i], 1, ctx));
   }
      
   CHKERR(globalSum_RHprimme(resNorms, basisSize, ctx));
   for (i=0; i < basisSize; i++)
      resNorms[i] = sqrt(resNorms[i]);

   /* Check for convergence of the residual norms. */

   CHKERR(check_convergence_Sprimme(V, ldV, 1 /* given X */, W, ldW,
         1 /* R given */, NULL, 0, 0, NULL, 0, NULL, 0, 0, basisSize, flags,
         resNorms, hVals, NULL, 0, ctx));

   /* Return the number of consecutive pairs converged  */

   for (i = 0; i < basisSize && flags[i] != UNCONVERGED; i++);
   *numConverged = i;

   return 0;
}

/******************************************************************************
           Dynamic Method Switching uses the following functions 
    ---------------------------------------------------------------------
     If primme->dynamicMethodSwitch > 0 find which of GD+k,JDQMR is best
     NOTE: JDQMR requires additional memory for inner iteration, so we 
       assume either the user has called primme_set_method(DYNAMIC) 
       or has allocated the appropriate space.
*******************************************************************************/

/******************************************************************************
 * Function switch_from_JDQMR - 
 *    If primme->dynamicMethodSwitch=2,4 try to switch from JDQMR_ETol to GD+k.
 *    Assumes that the CostModel has been updated through runtime measurements.
 *    Based on this CostModel, the switch occurs only if
 *
 *                          expected_JDQMR_ETol_time
 *                ratio =  --------------------------  > 1.05
 *                             expected_GD+k_time
 *
 *    There are two cases determining when this function is called.
 *
 *    primme->dynamicMethodSwitch = 2 
 *        numEvals < 5 (few eigenvalues). We must dynamically decide the 
 *        best method before an eigenvalue converges. Because a very slow
 *        and expensive inner iteration may take too many inner steps, we 
 *        must re-evaluate the ratio at every outer step, before the call 
 *        to solve the correction equation. 
 *        After a switch, dynamicMethodSwitch is 1.
 *        
 *    primme->dynamicMethodSwitch = 4 
 *        numEvals > 4 (many eigenvalues). We can afford to check how both 
 *        methods converge to an eigenvalue, and then measure the statistics.
 *        In this case we re-evaluate the ratio every time one or more 
 *        eigenvalues converge, just before the next correction equation.
 *        After a switch, dynamicMethodSwitch is 3.
 *
 * INPUT/OUTPUT
 * ------------
 * model        The CostModel that contains all the model relevant parameters
 *              (accum_jdq_gdk, accum_jdq, accum_gdk updated)
 *
 * primme       The solver parameters (dynamicMethodSwitch, maxInnerIterations
 *              changed if there is a switch)
 *
 ******************************************************************************/

STATIC int switch_from_JDQMR(void *model_, primme_context ctx) {

   primme_CostModel *model = (primme_CostModel *)model_;
   primme_params *primme = ctx.primme;
   int switchto=0;
   HREAL est_slowdown, est_ratio_MV_outer, ratio; 

   /* ----------------------------------------------------------------- */
   /* Asymptotic evaluation of the JDQMR versus GD+k for small numEvals */
   /* ----------------------------------------------------------------- */
   if (primme->dynamicMethodSwitch == 2) {

     /* For numEvals<4, (dyn=2), after we first estimate timings, decide if 
      * must always use GD+k (eg., because the operator is very expensive)
      * Use a best case scenario for JDQMR (small slowdown and many inner its)*/
      est_slowdown       = 1.1;     /* a small enough slowdown */
      est_ratio_MV_outer = 1000;    /* practically all time is in inner its */
      ratio = ratio_JDQMR_GDpk(model, 0, est_slowdown, est_ratio_MV_outer);

      /* Average the ratio among the processes */
      CHKERR(globalSum_RHprimme(&ratio, 1, ctx));
      ratio /= (HREAL)primme->numProcs;

      if (ratio > 1.05) { 
         /* Always use GD+k. No further model updates */
         primme->dynamicMethodSwitch = -1;
         primme->correctionParams.maxInnerIterations = 0;
         primme->correctionParams.projectors.RightX  = 1;
         PRINTF(3, "Ratio: %e Switching permanently to GD+k", ratio);
         return 0;
      }
   }

   /* Select method to switch to if needed: 2->1 and 4->3 */
   switch (primme->dynamicMethodSwitch) {
      case 2: switchto = 1; break;
      case 4: switchto = 3; 
   }

   /* ----------------------------------------------------------------- *
    * Compute the ratio of expected times JDQMR/GD+k. To switch to GD+k, the 
    * ratio must be > 1.05. Update accum_jdq_gdk for recommendation to user
    * ----------------------------------------------------------------- */

   ratio = ratio_JDQMR_GDpk(model, 0, model->JDQMR_slowdown, 
                                  model->ratio_MV_outer);

   /* Average the ratio among the processes */
   CHKERR(globalSum_RHprimme(&ratio, 1, ctx));
   ratio /= (HREAL)primme->numProcs;
   
   if (ratio > 1.05) {
      primme->dynamicMethodSwitch = switchto; 
      primme->correctionParams.maxInnerIterations = 0;
      primme->correctionParams.projectors.RightX  = 1;
   }

   model->accum_jdq += model->gdk_plus_MV_PR*ratio;
   model->accum_gdk += model->gdk_plus_MV_PR;
   model->accum_jdq_gdk = model->accum_jdq/model->accum_gdk;

   switch (primme->correctionParams.maxInnerIterations) {
   case 0: PRINTF(3, "Ratio: %e JDQMR switched to GD+k", ratio); break;
   case -1: PRINTF(3, "Ratio: %e Continue with JDQMR", ratio);
   }

   return 0;

} /* end of switch_fromJDQMR() */

/******************************************************************************
 * Function switch_from_GDpk - 
 *    If primme->dynamicMethodSwitch=1,3 try to switch from GD+k to JDQMR_ETol.
 *    Assumes that the CostModel has been updated through runtime measurements.
 *    If no JDQMR measurements exist (1st time), switch to it unconditionally.
 *    Otherwise, based on the CostModel, the switch occurs only if
 *
 *                          expected_JDQMR_ETol_time
 *                ratio =  --------------------------  < 0.95
 *                             expected_GD+k_time
 *
 *    There are two cases determining when this function is called.
 *
 *    primme->dynamicMethodSwitch = 1
 *        numEvals < 5 (few eigenvalues). GD+k does not have an inner iteration
 *        so it is not statistically meaningful to check every outer step. 
 *        For few eigenvalues, we re-evaluate the ratio immediately after
 *        restart of the method. This allows also the cost of restart to 
 *        be included in the measurements.
 *        After a switch, dynamicMethodSwitch is 2.
 *        
 *    primme->dynamicMethodSwitch = 3 
 *        numEvals > 4 (many eigenvalues). As with JDQMR, we can check how both 
 *        methods converge to an eigenvalue, and then measure the statistics.
 *        In this case, we re-evaluate the ratio every time one or more 
 *        eigenvalues converge, just before the next preconditioner application
 *        After a switch, dynamicMethodSwitch is 4.
 *
 * INPUT/OUTPUT
 * ------------
 * model        The CostModel that contains all the model relevant parameters
 *              (accum_jdq_gdk, accum_jdq, accum_gdk updated)
 *
 * primme       The solver parameters (dynamicMethodSwitch, maxInnerIterations
 *              changed if there is a switch)
 *
 ******************************************************************************/
STATIC int switch_from_GDpk(void *model_, primme_context ctx) {

   primme_CostModel *model = (primme_CostModel *)model_;
   primme_params *primme = ctx.primme;
   int switchto=0;
   HREAL ratio;

   /* if no restart has occurred (only possible under dyn=3) current timings */
   /* do not include restart costs. Remain with GD+k until a restart occurs. */
   /* If search space is going to saturate, just use GD for ever.            */
   if (primme->stats.numRestarts == 0 ||
         primme->maxBasisSize + (primme->locking ? primme->numEvals : 0) >=
               primme->n) {
      return 0;
   }

   /* Select method to switch to if needed: 1->2 and 3->4 */
   switch (primme->dynamicMethodSwitch) {
      case 1: switchto = 2; break;
      case 3: switchto = 4; 
   }

   /* If JDQMR never run before, switch to it to get first measurements */
   if (model->qmr_only == 0.0) {
      primme->dynamicMethodSwitch = switchto;
      primme->correctionParams.maxInnerIterations = -1;
      primme->correctionParams.projectors.RightX  = 0;
      PRINTF(3, "Ratio: N/A  GD+k switched to JDQMR (first time)");
      return 0;
   }

   /* ------------------------------------------------------------------- *
    * Compute the ratio of expected times JDQMR/GD+k. To switch to JDQMR, the
    * ratio must be < 0.95. Update accum_jdq_gdk for recommendation to user
    * ------------------------------------------------------------------- */

   ratio = ratio_JDQMR_GDpk(model, 0, model->JDQMR_slowdown, 
                                  model->ratio_MV_outer);

   /* Average the ratio among the processes */
   CHKERR(globalSum_RHprimme(&ratio, 1, ctx));
   ratio /= (HREAL)primme->numProcs;

   if (ratio < 0.95) {
      primme->dynamicMethodSwitch = switchto;
      primme->correctionParams.maxInnerIterations = -1;
      primme->correctionParams.projectors.RightX  = 0;
   } 

   model->accum_jdq += model->gdk_plus_MV_PR*ratio;
   model->accum_gdk += model->gdk_plus_MV_PR;
   model->accum_jdq_gdk = model->accum_jdq/model->accum_gdk;

   switch (primme->correctionParams.maxInnerIterations) {
   case 0: PRINTF(3, "Ratio: %e Continue with GD+k", ratio); break;
   case -1: PRINTF(3, "Ratio: %e GD+k switched to JDQMR", ratio);
   }

   return 0;
}

/******************************************************************************
 * Function update_statistics - 
 *
 *    Performs runtime measurements and updates the cost model that describes
 *    the average cost of Matrix-vector and Preconditioning operations, 
 *    the average cost of running 1 full iteration of GD+k and of JDQMR, 
 *    the number of inner/outer iterations since last update,
 *    the current convergence rate measured for each of the two methods,
 *    and based on these rates, the expected slowdown of JDQMR over GD+k
 *    in terms of matrix-vector operations.
 *    Times are averaged with one previous measurement, and convergence
 *    rates are averaged over a window which is reset over 10 converged pairs.
 *
 *    The function is called right before switch_from_JDQMR/switch_from_GDpk.
 *    Depending on the current method running, this is at two points:
 *       If some eigenvalues just converged, called before solve_correction()
 *       If primme->dynamicMethodSwitch = 2, called before solve_correction()
 *       If primme->dynamicMethodSwitch = 1, called after restart()
 *
 *    The Algorithm
 *    -------------
 *    The dynamic switching algorithm starts with GD+k. After the 1st restart
 *    (for dyn=1) or after an eigenvalue converges (dyn=3), we collect the 
 *    first measurements for GD+k and switch to JDQMR_ETol. We collect the 
 *    first measurements for JDQMR after 1 outer step (dyn=2) or after an
 *    eigenvalue converges (dyn=4). From that time on, the method is chosen
 *    dynamically by switch_from_JDQMR/switch_from_GDpk using the ratio.
 *
 *    Cost/iteration breakdown
 *    ------------------------
 *    kout: # of outer iters. kinn: # of inner QMR iters. nMV = # of Matvecs
 *    All since last call to update_statistics.
 *    1 outer step: costJDQMR = costGD_outer_method + costQMR_Iter*kinn + mv+pr
 *    
 *        <---------1 step JDQMR----------->
 *       (GDout)(---------QMR--------------)
 *        gd mv  pr q+mv+pr .... q+mv+pr mv  
 *                  <-------kinn------->      kinn = nMV/kout - 2
 *        (-----)(------time_in_inner------)  
 *
 *    The model
 *    ---------
 *    time_in_inner = kout (pr+kinn*(q+mv+pr)+mv) 
 *                  = kout (pr+mv) + nMV*(q+mv+pr) - 2kout(q+mv+pr)
 *    time_in_outer = elapsed_time-time_in_inner = kout*(gd+mv)
 *    JDQMR_time = time_in_inner+time_in_outer = 
 *               = kout(gd+mv) + kout(pr+mv) + nMV(q+mv+pr) -2kout(q+mv+pr)
 *               = kout (gd -2q -pr) + nMV(q+mv+pr)
 *    GDpk_time  = gdOuterIters*(gd+mv+pr)
 *
 *    Letting slowdown = (nMV for JDQMR)/(gdOuterIters) we have the ratio:
 *
 *          JDQMR_time     q+mv+pr + kout/nMV (gd-2q-pr)
 *          ----------- = ------------------------------- slowdown
 *           GDpk_time               gd+mv+pr
 *
 *    Because the QMR of JDQMR "wastes" one MV per outer step, and because
 *    its number of outer steps cannot be more than those of GD+k
 *           (kinn+2)/(kinn+1) < slowdown < kinn+2
 *    However, in practice 1.1 < slowdown < 2.5. 
 *
 *    For more details see papers [1] and [2] (listed in primme_S.c)
 *
 *
 * INPUT
 * -----
 * primme           Structure containing the solver parameters
 * current_time     Time stamp taken before update_statistics is called
 * recentConv       Number of converged pairs since last update_statistics
 * calledAtRestart  True if update_statistics is called at restart by dyn=1
 * numConverged     Total number of converged pairs
 * currentResNorm   Residual norm of the next unconverged epair
 *
 * INPUT/OUTPUT
 * ------------
 * model            The model parameters updated
 *
 * Return value
 * ------------
 *     0    Not enough iterations to update model. Continue with current method
 *     1    Model updated. Proceed with relative evaluation of the methods
 *
 ******************************************************************************/
STATIC int update_statistics(void *model_, primme_params *primme,
      double current_time, int recentConv, int calledAtRestart,
      int numConverged, double currentResNorm) {

   primme_CostModel *model = (primme_CostModel *)model_;
   double low_res, elapsed_time, time_in_outer, kinn;
   int kout, nMV;

   /* ------------------------------------------------------- */
   /* Time in outer and inner iteration since last update     */
   /* ------------------------------------------------------- */
   elapsed_time = current_time - model->timer_0;
   time_in_outer = elapsed_time - model->time_in_inner;

   /* ------------------------------------------------------- */
   /* Find # of outer, MV, inner iterations since last update */
   /* ------------------------------------------------------- */
   kout = primme->stats.numOuterIterations - model->numIt_0;
   nMV  = primme->stats.numMatvecs - model->numMV_0;
   if (calledAtRestart) kout++; /* Current outer iteration is complete,  */
                                /*   but outerIterations not incremented yet */
   if (kout == 0) return 0;     /* No outer iterations. No update or evaluation
                                 *   Continue with current method */
   kinn = ((double) nMV)/kout - 2;

   if (primme->correctionParams.maxInnerIterations == -1 && 
       (kinn < 1.0 && model->qmr_only == 0.0L) )  /* No inner iters yet, and */
      return 0;              /* no previous QMR timings. Continue with JDQMR */

   /* --------------------------------------------------------------- */
   /* After one or more pairs converged, currentResNorm corresponds   */
   /* to the next unconverged pair. To measure the residual reduction */
   /* during the previous step, we must use the convergence tolerance */
   /* Also, update how many evals each method found since last reset  */
   /* --------------------------------------------------------------- */
   if (recentConv > 0) {
      /* Use tolerance as the largest residual norm of the converged pairs */
      low_res = primme->stats.maxConvTol;
      /* Update num of evals found */
      if (primme->correctionParams.maxInnerIterations == -1)
          model->nevals_by_jdq += recentConv;
      else
          model->nevals_by_gdk += recentConv;
   }
   else 
      /* For dyn=1 at restart, and dyn=2 at every step. Use current residual */
      low_res = currentResNorm;

   /* ------------------------------------------------------- */
   /* Update model timings and parameters                     */
   /* ------------------------------------------------------- */

   /* update outer iteration time for both GD+k,JDQMR.Average last two updates*/
   if (model->gdk_plus_MV == 0.0L) 
      model->gdk_plus_MV = time_in_outer/kout;
   else 
      model->gdk_plus_MV = (model->gdk_plus_MV + time_in_outer/kout)/2.0L;

   /* ---------------------------------------------------------------- *
    * Reset the conv rate averaging window every 10 converged pairs. 
    * See also Notes 2 and 3 below.
    * Note 1: For large numEvals, we should average the convergence 
    * rate over only the last few converged pairs. To avoid a more 
    * expensive moving window approach, we reset the window when >= 10 
    * additional pairs converge. To avoid a complete reset, we consider 
    * the current average rate to be the new rate for the "last" pair.
    * By scaling down the sums: sum_logResReductions and sum_MV, this 
    * rate does not dominate the averaging of subsequent measurements.
    * ---------------------------------------------------------------- */

   if (numConverged/10 >= model->nextReset) {
      model->gdk_sum_logResReductions /= model->nevals_by_gdk;
      model->gdk_sum_MV /= model->nevals_by_gdk;
      model->jdq_sum_logResReductions /= model->nevals_by_jdq;
      model->jdq_sum_MV /= model->nevals_by_jdq;
      model->nextReset = numConverged/10+1;
      model->nevals_by_gdk = 1;
      model->nevals_by_jdq = 1;
   }

   switch (primme->dynamicMethodSwitch) {

      case 1: case 3: /* Currently running GD+k */
        /* Update Precondition times */
        if (model->PR == 0.0L) 
           model->PR          = model->time_in_inner/kout;
        else 
           model->PR          = (model->PR + model->time_in_inner/kout)/2.0L;
        model->gdk_plus_MV_PR = model->gdk_plus_MV + model->PR;
        model->MV_PR          = model->MV + model->PR;

        /* update convergence rate.
         * ---------------------------------------------------------------- *
         * Note 2: This is NOT a geometric average of piecemeal rates. 
         * This is the actual geometric average of all rates per MV, or 
         * equivalently the total rate as TotalResidualReductions over
         * the corresponding nMVs. If the measurement intervals (in nMVs)
         * are identical, the two are equivalent. 
         * Note 3: In dyn=1,2, we do not record residual norm increases.
         * This overestimates convergence rates a little, but otherwise, a 
         * switch would leave the current method with a bad rate estimate.
         * ---------------------------------------------------------------- */

        if (low_res <= model->resid_0) 
           model->gdk_sum_logResReductions += log(low_res/model->resid_0);
        model->gdk_sum_MV += nMV;
        model->gdk_conv_rate = exp( model->gdk_sum_logResReductions
                                   /model->gdk_sum_MV );
        break;

      case 2: case 4: /* Currently running JDQMR */
        /* Basic timings for QMR iteration (average of last two updates) */
        if (model->qmr_plus_MV_PR == 0.0L) {
           model->qmr_plus_MV_PR=(model->time_in_inner/kout- model->MV_PR)/kinn;
           model->ratio_MV_outer = ((double) nMV)/kout;
        }
        else {
           if (kinn != 0.0) model->qmr_plus_MV_PR = (model->qmr_plus_MV_PR  + 
              (model->time_in_inner/kout - model->MV_PR)/kinn )/2.0L;
           model->ratio_MV_outer =(model->ratio_MV_outer+((double) nMV)/kout)/2;
        }
        model->qmr_only = model->qmr_plus_MV_PR - model->MV_PR;

        /* Update the cost of a hypothetical GD+k, as measured outer + PR */
        model->gdk_plus_MV_PR = model->gdk_plus_MV + model->PR;

        /* update convergence rate */
        if (low_res <= model->resid_0) 
           model->jdq_sum_logResReductions += log(low_res/model->resid_0);
        model->jdq_sum_MV += nMV;
        model->jdq_conv_rate = exp( model->jdq_sum_logResReductions
                                   /model->jdq_sum_MV);
        break;
   }
   update_slowdown(model);

   /* ------------------------------------------------------- */
   /* Reset counters to measure statistics at the next update */
   /* ------------------------------------------------------- */
   model->numIt_0 = primme->stats.numOuterIterations;
   if (calledAtRestart) model->numIt_0++; 
   model->numMV_0 = primme->stats.numMatvecs;
   model->timer_0 = current_time;      
   model->time_in_inner = 0.0;
   model->resid_0 = currentResNorm;

   return 1;
}

/******************************************************************************
 * Function ratio_JDQMR_GDpk -
 *    Using model parameters, computes the ratio of expected times:
 *
 *          JDQMR_time     q+mv+pr + kout/nMV (gd-2q-pr)
 *          ----------- = ------------------------------- slowdown
 *           GDpk_time               gd+mv+pr
 *
 ******************************************************************************/
STATIC double ratio_JDQMR_GDpk(void *model_, int numLocked,
   double estimate_slowdown, double estimate_ratio_MV_outer) {
   
   primme_CostModel *model = (primme_CostModel *)model_;
   return estimate_slowdown* 
     ( model->qmr_plus_MV_PR + model->project_locked*numLocked + 
       (model->gdk_plus_MV - model->qmr_only - model->qmr_plus_MV_PR  
        + (model->reortho_locked - model->project_locked)*numLocked  
       )/estimate_ratio_MV_outer
     ) / (model->gdk_plus_MV_PR + model->reortho_locked*numLocked);
}

/******************************************************************************
 * Function update_slowdown -
 *    Given the model measurements for convergence rates, computes the slowdown
 *           log(GDpk_conv_rate)/log(JDQMR_conv_rate)
 *    subject to the bounds 
 *    max(1.1, (kinn+2)/(kinn+1)) < slowdown < min(2.5, kinn+2)
 *
 ******************************************************************************/
STATIC void update_slowdown(void *model_) {
   primme_CostModel *model = (primme_CostModel *)model_;
  double slowdown;

  if (model->gdk_conv_rate < 1.0) {
     if (model->jdq_conv_rate < 1.0) 
        slowdown = log(model->gdk_conv_rate)/log(model->jdq_conv_rate);
     else if (model->jdq_conv_rate == 1.0)
        slowdown = 2.5;
     else
        slowdown = -log(model->gdk_conv_rate)/log(model->jdq_conv_rate);
  }
  else if (model->gdk_conv_rate == 1.0) 
        slowdown = 1.1;
  else { /* gdk > 1 */
     if (model->jdq_conv_rate < 1.0) 
        slowdown = log(model->gdk_conv_rate)/log(model->jdq_conv_rate);
     else if (model->jdq_conv_rate == 1.0)
        slowdown = 1.1;
     else  /* both gdk, jdq > 1 */
        slowdown = log(model->jdq_conv_rate)/log(model->gdk_conv_rate);
  }

  /* Slowdown cannot be more than the matvecs per outer iteration */
  /* nor less than MV per outer iteration/(MV per outer iteration-1) */
  slowdown = max(model->ratio_MV_outer/(model->ratio_MV_outer-1.0),
                                    min(slowdown, model->ratio_MV_outer));
  /* Slowdown almost always in [1.1, 2.5] */
  model->JDQMR_slowdown = max(1.1, min(slowdown, 2.5));
}

/******************************************************************************
 * Function initializeModel - Initializes model members
 ******************************************************************************/
STATIC void initializeModel(void *model_, primme_params *primme) {
   primme_CostModel *model = (primme_CostModel *)model_;
   model->MV_PR          = 0.0L;
   model->MV             = 0.0L;
   model->PR             = 0.0L;
   model->qmr_only       = 0.0L;
   model->qmr_plus_MV_PR = 0.0L;
   model->gdk_plus_MV_PR = 0.0L;
   model->gdk_plus_MV    = 0.0L;
   model->project_locked = 0.0L;
   model->reortho_locked = 0.0L;

   model->gdk_conv_rate  = 0.0001L;
   model->jdq_conv_rate  = 0.0001L;
   model->JDQMR_slowdown = 1.5L;
   model->ratio_MV_outer = 0.0L;

   model->nextReset      = 1;
   model->gdk_sum_logResReductions = 0.0;
   model->jdq_sum_logResReductions = 0.0;
   model->gdk_sum_MV     = 0.0;
   model->jdq_sum_MV     = 0.0;
   model->nevals_by_gdk  = 0;
   model->nevals_by_jdq  = 0;

   model->numMV_0 = primme->stats.numMatvecs;
   model->numIt_0 = primme->stats.numOuterIterations+1;
   model->timer_0 = primme_wTimer();
   model->time_in_inner  = 0.0L;
   model->resid_0        = -1.0L;

   model->accum_jdq      = 0.0L;
   model->accum_gdk      = 0.0L;
   model->accum_jdq_gdk  = 1.0L;
}

#if 0
/******************************************************************************
 *
 * Function to display the model parameters --- For debugging purposes only
 *
 ******************************************************************************/
STATIC void displayModel(primme_CostModel *model){
   fprintf(stdout," MV %e\n", model->MV);
   fprintf(stdout," PR %e\n", model->PR);
   fprintf(stdout," MVPR %e\n", model->MV_PR);
   fprintf(stdout," QMR %e\n", model->qmr_only);
   fprintf(stdout," QMR+MVPR %e\n", model->qmr_plus_MV_PR);
   fprintf(stdout," GD+MVPR %e\n", model->gdk_plus_MV_PR);
   fprintf(stdout," GD+MV %e\n\n", model->gdk_plus_MV);
   fprintf(stdout," project %e\n", model->project_locked);
   fprintf(stdout," reortho %e\n", model->reortho_locked);
   fprintf(stdout," GD convrate %e\n", model->gdk_conv_rate);
   fprintf(stdout," JD convrate %e\n", model->jdq_conv_rate);
   fprintf(stdout," slowdown %e\n", model->JDQMR_slowdown);
   fprintf(stdout," outer/Totalmv %e\n", model->ratio_MV_outer);
   
   fprintf(stdout," gd sum_logResRed %e\n",model->gdk_sum_logResReductions);
   fprintf(stdout," jd sum_logResRed %e\n",model->jdq_sum_logResReductions);
   fprintf(stdout," gd sum_MV %e\n", model->gdk_sum_MV);
   fprintf(stdout," jd sum_MV %e\n", model->jdq_sum_MV);
   fprintf(stdout," gd nevals %d\n", model->nevals_by_gdk);
   fprintf(stdout," jd nevals %d\n", model->nevals_by_jdq);

   fprintf(stdout," Accumul jd %e\n", model->accum_jdq);
   fprintf(stdout," Accumul gd %e\n", model->accum_gdk);
   fprintf(stdout," Accumul ratio %e\n", model->accum_jdq_gdk);

   fprintf(stdout,"0: MV %d IT %d timer %e timInn %e res %e\n", model->numMV_0,
      model->numIt_0,model->timer_0,model->time_in_inner, model->resid_0);
   fprintf(stdout," ------------------------------\n");
}
#endif

#endif /* SUPPORTED_TYPE */
