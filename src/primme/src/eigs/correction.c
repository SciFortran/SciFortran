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
 * File: correction.c
 *
 * Purpose - Computes the correction corresponding to each block Ritz vectors.
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/correction.c"
#endif

#include "numerical.h"
#include "template_normal.h"
#include "common_eigs.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "correction.h"
#include "inner_solve.h"
#include "auxiliary_eigs.h"
#endif

#ifdef SUPPORTED_TYPE

/*******************************************************************************
 * Subroutine solve_correction - This routine solves the correction equation
 *    for each Ritz vector and residual in the block.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * W              The last blockSize vectors of W contain the residuals
 *
 * evecs          The converged Ritz vectors.  Array is of dimension numLocked.
 *
 * evecsHat       K^{-1}evecs given a preconditioner K. 
 *                (accessed only if skew projector is requested)
 *
 * Mfact          The factorization of the matrix evecs'*evecsHat
 * ipivot         The pivots for the Mfact factorization
 *
 * lockedEvals    The locked eigenvalues
 *
 * numLocked      The number of locked eigenvalues (zero if not locking)
 *
 * numConvergedStored  Number of converged eigenvectors in V that have been
 *                copied in evecs when no locking is employed, to accommodate 
 *                for a requested skew projection with evecs. 
 *                   
 * ritzVals       Array of size basisSize. The Ritz values corresponding to 
 *                the current basis V
 *
 * flags          Array indicating which eigenvectors have converged
 *
 * basisSize      The size of the basis V
 *
 * iev            Array of size blockSize.  The index of each block Ritz value
 *                in ritzVals (and each Ritz vector in V)
 *
 * blockSize      The current block size
 *
 * machEps        machine precision 
 *
 * rwork          Real workspace of size          
 *                3*maxEvecsSize + 2*primme->maxBlockSize 
 *                + (primme->numEvals+primme->maxBasisSize)
 *                        *----------------------------------------------------*
 *                        | The following are optional and mutually exclusive: |
 *                        *------------------------------+                     |
 *                + 4*primme->nLocal + primme->nLocal    | For QMR work and sol|
 *                + primme->nLocal*primme->maxBlockSize  | OLSEN for KinvBx    |
 *                                                       *---------------------*
 *
 * rworkSize      the size of rwork. If less than needed, func returns needed.
 *
 * iwork          Integer workspace of size maxBlockSize 
 * 
 * primme         Structure containing various solver parameters 
 * 
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * V              The orthonormal basis.  The last blockSize vectors of V
 *                contain the Ritz vectors
 *
 * blockNorms     On input, residual norms of the Ritz vectors computed 
 *                during the current outer iteration.
 *                On output, the approximate residual norms of the corrected
 *                Ritz vectors (through JDQMR only).
 *
 *
 * prevRitzVals   Array of size numPrevRitzVals.  The Ritz values from the
 *                previous iteration and converged Ritz values
 *
 * numPrevRitzVals  The of size prevRitzVals updated every outer step
 *
 * touch            Parameter used in inner solve stopping criteria
 *
 * Return Value
 * ------------
 * int  Error code: 0 upon success, nonzero otherwise
 *                 -1 innner solver failure
 *                 >0 the needed rworkSize, if the given was not enough
 *
 ******************************************************************************/
 
TEMPLATE_PLEASE
int solve_correction_Sprimme(SCALAR *V, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, SCALAR *BV, PRIMME_INT ldBV, SCALAR *evecs,
      PRIMME_INT ldevecs, SCALAR *Bevecs, PRIMME_INT ldBevecs, SCALAR *evecsHat,
      PRIMME_INT ldevecsHat, HSCALAR *Mfact, int *ipivot, HEVAL *lockedEvals,
      int numLocked, int numConvergedStored, HEVAL *ritzVals,
      HEVAL *prevRitzVals, int *numPrevRitzVals, int *flags, int basisSize,
      HREAL *blockNorms, int *iev, int blockSize, int *touch, double startTime,
      primme_context ctx) {

   KIND(, (void)lockedEvals);
   KIND(, (void)flags);
   KIND(, (void)evecs);
   KIND(, (void)ldevecs);
   KIND(, (void)Bevecs);
   KIND(, (void)ldBevecs);
   KIND(, (void)evecsHat);
   KIND(, (void)ldevecsHat);
   KIND(, (void)Mfact);
   KIND(, (void)ipivot);
   KIND(, (void)numConvergedStored);
   KIND(, (void)touch);
   KIND(, (void)startTime);

   primme_params *primme = ctx.primme;
   int blockIndex; /* Loop index.  Ranges from 0..blockSize-1.       */
   int *ilev;      /* Array of size blockSize.  Maps the target Ritz */
                   /* values to their positions in the sortedEvals   */
                   /* array.                                         */

   HEVAL *sortedRitzVals; /* Sorted array of current and converged Ritz     */
                          /* values.  Size of array is numLocked+basisSize. */
   KIND(double, PRIMME_COMPLEX_DOUBLE) *
         blockOfShifts;   /* Shifts for (A-shiftI) or (if needed) (K-shiftI)*/
   HREAL *approxOlsenEps; /* Shifts for approximate Olsen implementation    */

   CHKERR(KIND(Num_malloc_dprimme, Num_malloc_zprimme)(
         blockSize, &blockOfShifts, ctx));
   CHKERR(Num_malloc_RHprimme(blockSize, &approxOlsenEps, ctx));

   /*------------------------------------------------------------*/
   /*  Figuring out preconditioning shifts  (robust, Olsen, etc) */
   /*------------------------------------------------------------*/
   /* blockOfShifts will contain the preconditioning shifts:               */
   /* either Ritz values or robustShifts computed below. These shifts      */
   /* will be used in the correction equations or in inverting (K-sigma I) */
   /* approxOlsenEps will contain error approximations for eigenavalues    */
   /* to be used for Olsen's method (when innerIterations =0).             */

#ifdef USE_HERMITIAN
   if (primme->locking && (primme->target == primme_smallest ||
                                primme->target == primme_largest)) {
      /* Combine the sorted list of locked Ritz values with the sorted  */
      /* list of current Ritz values, ritzVals.  The merging of the two */
      /* lists lockedEvals and ritzVals is stored in sortedRitzVals.    */

      CHKERR(Num_malloc_RHprimme(numLocked + basisSize, &sortedRitzVals, ctx));
      CHKERR(Num_malloc_iprimme(blockSize, &ilev, ctx));
      mergeSort(lockedEvals, numLocked, ritzVals, flags, basisSize,
            sortedRitzVals, ilev, blockSize, primme);
   } else
#endif /* USE_HERMITIAN */
   {
      /* In the case of soft-locking or when we look for interior ones  */
      /* the sorted evals are simply the ritzVals, targeted as iev      */

      sortedRitzVals = ritzVals;
      ilev = iev;
   }

   // For interior pairs we an approach based on the target shift given by the
   // user. For Rayleigh-Ritz, the Ritz value does not get monotonically closer
   // to the exact value, and it is not trusted as a good shift until its
   // residual vector norm is small. We capture this hint in the following
   // heuristic: take the closest point in the interval Ritz value +- residual
   // norm to the target as the shift. For refined extraction,
   //     |Ritz value - target| <= the Ritz singular value.
   // So we use the Ritz value as the shift as soon as the Ritz value is closer
   // to the Ritz singular value than to the target.

   if (primme->target != primme_smallest && primme->target != primme_largest) {

      for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {

         double targetShift =
               primme->numTargetShifts > 0
                     ? primme->targetShifts[min(
                             primme->numTargetShifts - 1, numLocked)]
                     : 0.0;
         int sortedIndex = ilev[blockIndex];
         if (EVAL_ABS(sortedRitzVals[sortedIndex] - (HEVAL)targetShift) <
               blockNorms[blockIndex] * sqrt(primme->stats.estimateInvBNorm)) {
            blockOfShifts[blockIndex] = targetShift;
         } else if (primme->projectionParams.projection ==
                    primme_proj_refined) {
            blockOfShifts[blockIndex] = sortedRitzVals[sortedIndex];
         } else {
            blockOfShifts[blockIndex] =
                  sortedRitzVals[sortedIndex] +
                  ((HEVAL)(blockNorms[blockIndex] *
                           sqrt(primme->stats.estimateInvBNorm))) *
                        ((HEVAL)targetShift - sortedRitzVals[sortedIndex]) /
                        EVAL_ABS(
                              (HEVAL)targetShift - sortedRitzVals[sortedIndex]);
         }

         if (sortedIndex < *numPrevRitzVals) {
            approxOlsenEps[blockIndex] = EVAL_ABS(
                  prevRitzVals[sortedIndex] - sortedRitzVals[sortedIndex]);
         }  
         else {
            approxOlsenEps[blockIndex] =
                  blockNorms[blockIndex] * sqrt(primme->stats.estimateInvBNorm);
         }  
      } /* for loop */

      /* Remember the previous ritz values*/
      *numPrevRitzVals = basisSize;
      CHKERR(KIND(Num_copy_RHprimme, Num_copy_SHprimme)(
            *numPrevRitzVals, sortedRitzVals, 1, prevRitzVals, 1, ctx));

   } /* user provided shifts */
#ifdef USE_HERMITIAN
   else {
      /*-----------------------------------------------------------------*/
      /* else it is primme_smallest or primme_largest                    */
      /*-----------------------------------------------------------------*/

      if (primme->correctionParams.robustShifts) { 
         /*----------------------------------------------------*/
         /* Subtract/add a robust shift from/to the Ritz value */
         /*----------------------------------------------------*/

         /* Find the robust shift for each block vector */
         for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {
   
            int sortedIndex = ilev[blockIndex];
            HREAL eval = sortedRitzVals[sortedIndex];

            HREAL robustShift = computeRobustShift(blockIndex,
                  blockNorms[blockIndex], prevRitzVals, *numPrevRitzVals,
                  sortedRitzVals, &approxOlsenEps[blockIndex],
                  numLocked + basisSize, ilev, primme);

            /* Subtract/add the shift if looking for the smallest/largest  */
            /* eigenvalues, Do not go beyond the previous computed eigval  */
       
            if (primme->target == primme_smallest) {
               blockOfShifts[blockIndex] = eval - robustShift;
               if (sortedIndex > 0) blockOfShifts[blockIndex] = 
                  max(blockOfShifts[blockIndex], sortedRitzVals[sortedIndex-1]);
            }
            else {
               blockOfShifts[blockIndex] = eval + robustShift;
               if (sortedIndex > 0) blockOfShifts[blockIndex] = 
                  min(blockOfShifts[blockIndex], sortedRitzVals[sortedIndex-1]);
            } /* robust shifting */
   
         }  /* for loop */
   
      }  /* endif robust shifts */
      else {
         /*--------------------------------------------------------------*/
         /* Otherwise, the shifts for both preconditioner and correction */
         /* equation should be just the Ritz values. For Olsen's method, */
         /* the shifts for r-eps*x, are chosen as the difference in Ritz */
         /* value between successive iterations.                         */
         /*--------------------------------------------------------------*/
   
         for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {
            int ritzIndex   =  iev[blockIndex];
            int sortedIndex = ilev[blockIndex];
            blockOfShifts[blockIndex] = ritzVals[ritzIndex];
            if (sortedIndex < *numPrevRitzVals) {
               approxOlsenEps[blockIndex] = fabs(
                     prevRitzVals[sortedIndex] - sortedRitzVals[sortedIndex]);
            }
            else {
               approxOlsenEps[blockIndex] =
                     blockNorms[blockIndex] *
                     sqrt(primme->stats.estimateInvBNorm);
            }
         } /* for loop */
      } /* else no robust shifts */

      /* Remember the previous ritz values*/
      *numPrevRitzVals = numLocked+basisSize;
      Num_copy_RHprimme(
            *numPrevRitzVals, sortedRitzVals, 1, prevRitzVals, 1, ctx);

   } /* else primme_smallest or primme_largest */
#endif /* USE_HERMITIAN */

   /* Equip the primme struct with the blockOfShifts, in case the user */
   /* wants to precondition (K-sigma_i I)^{-1} with a different shift  */
   /* for each vector                                                  */

   primme->ShiftsForPreconditioner = (double *)blockOfShifts;

   /*------------------------------------------------------------ */
   /*  Generalized Davidson variants -- No inner iterations       */
   /*------------------------------------------------------------ */
   if (primme->correctionParams.maxInnerIterations == 0) {
      /* This is Generalized Davidson or approximate Olsen's method. */
      /* Perform block preconditioning (with or without projections) */

      SCALAR *r = &W[ldW*basisSize];    /* All the block residuals    */
      SCALAR *x = &V[ldV*basisSize];    /* All the block Ritz vectors */
      SCALAR *Bx = &BV[ldBV*basisSize]; /* B*x                        */
      
      if ( primme->correctionParams.projectors.RightX &&
           primme->correctionParams.projectors.SkewX    ) {    
           /* Compute exact Olsen's projected preconditioner. This is */
          /* expensive and rarely improves anything! Included for completeness*/
          
          Olsen_preconditioner_block(r, ldW, x, ldV, Bx, ldBV, blockSize, ctx);
      }
      else {
         // In general the Olsen's approximation is not useful without
         // preconditioning. However expanding the basis with r - approxOlsen*x
         // helps in checking the practical convergence. Check main_iter for
         // more details.

         if (primme->correctionParams.projectors.RightX &&
               ((primme->correctionParams.precondition &&
                      primme->applyPreconditioner) ||
                     Bx != x ||
                     (primme->locking &&
                           primme->orth == primme_orth_implicit_I))) {
#ifdef USE_HERMITIAN
            /*Compute a cheap approximation to OLSENS, where (x'Kinvr)/xKinvx */
            /*is approximated by e: Kinvr-e*KinvBx=Kinv(r-e*x)=Kinv(I-ct*x*x')r*/

            for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {
               /* Compute r_i = r_i - err_i * Bx_i */
               Num_axpy_Sprimme(primme->nLocal, -approxOlsenEps[blockIndex],
                     &Bx[ldBV * blockIndex], 1, &r[ldW * blockIndex], 1, ctx);
            }
#else
            CHKERR(PRIMME_FUNCTION_UNAVAILABLE);
#endif
         }

         /* GD: compute K^{-1}r , or approx.Olsen: K^{-1}(r-eBx) */

         CHKERR(applyPreconditioner_Sprimme(r, primme->nLocal, ldW,
                  x, ldV, blockSize, ctx));
      }
   }
   /* ------------------------------------------------------------ */
   /*  JDQMR --- JD inner-outer variants                           */
   /* ------------------------------------------------------------ */
   else {  /* maxInnerIterations > 0  We perform inner-outer JDQMR */
#ifdef USE_HERMITIAN
      int touch0 = *touch;

      SCALAR *sol, *KinvBx;
      PRIMME_INT ldsol = primme->ldOPs, ldKinvBx = primme->ldOPs;
      CHKERR(Num_malloc_Sprimme(ldsol * blockSize, &sol, ctx));
      CHKERR(Num_malloc_Sprimme(ldKinvBx * blockSize, &KinvBx, ctx));


      HSCALAR *xKinvBx;  /* Stores x'*K^{-1}Bx if needed    */
      CHKERR(Num_malloc_SHprimme(blockSize, &xKinvBx, ctx));

      SCALAR *r = &W[ldW * basisSize];
      SCALAR *x = &V[ldV * basisSize];
      SCALAR *Bx = &BV[ldBV * basisSize];

      /* Set up the left/right/skew projectors for JDQMR.        */
      /* The pointers Lprojector, Rprojector(Q/X) point to the   */
      /* appropriate arrays for use in the projection step       */

      int sizeLprojectorQ; /* Sizes of the various left/right projectors     */
      int sizeLprojectorX; /* These will be 0/1/or numOrthConstr+numLocked   */
      int sizeRprojectorQ; /* or numOrthConstr+numConvergedStored w/o locking*/
      int sizeRprojectorX;
      SCALAR *LprojectorQ;      /* Q pointer for (I-Q*Q'). Usually points to evecs*/
      SCALAR *LprojectorX;      /* Usually points to x */
      SCALAR *LprojectorBQ;     /* BQ pointer for (I-BQ*Q'). Usually points to Bevecs*/
      SCALAR *LprojectorBX;     /* Usually points to Bx */
      SCALAR *RprojectorQ;      /* May point to evecs/evecsHat depending on skewQ */
      SCALAR *RprojectorX;      /* May point to x/Kinvx depending on skewX        */
      PRIMME_INT ldLprojectorQ; /* The leading dimension of LprojectorQ    */
      PRIMME_INT ldLprojectorX; /* The leading dimension of LprojectorX    */
      PRIMME_INT ldLprojectorBQ;/* The leading dimension of LprojectorBQ   */
      PRIMME_INT ldLprojectorBX;/* The leading dimension of LprojectorBX   */
      PRIMME_INT ldRprojectorQ; /* The leading dimension of RprojectorQ    */
      PRIMME_INT ldRprojectorX; /* The leading dimension of RprojectorX    */

      CHKERR(setup_JD_projectors(x, ldV, Bx, ldBV, evecs, ldevecs, Bevecs,
            ldBevecs, evecsHat, ldevecsHat, KinvBx, ldKinvBx, xKinvBx,
            &LprojectorQ, &ldLprojectorQ, &LprojectorX, &ldLprojectorX,
            &LprojectorBQ, &ldLprojectorBQ, &LprojectorBX, &ldLprojectorBX,
            &RprojectorQ, &ldRprojectorQ, &RprojectorX, &ldRprojectorX,
            &sizeLprojectorQ, &sizeLprojectorX, &sizeRprojectorQ,
            &sizeRprojectorX, numLocked, numConvergedStored, blockSize, ctx));

      /* Map the index of the block vector to its corresponding eigenvalue */
      /* index, and the shift for the correction equation. Also make the   */
      /* shift available to PRIMME, in case (K-shift B)^-1 is needed       */

      HEVAL *blockRitzVals;
      CHKERR(KIND(Num_malloc_RHprimme, Num_malloc_SHprimme)(
            blockSize, &blockRitzVals, ctx));
      KIND(Num_compact_vecs_RHprimme, Num_compact_vecs_SHprimme)
      (ritzVals, 1, blockSize, 1, iev, blockRitzVals, 1, 0 /* force copy */,
            ctx);
      primme->ShiftsForPreconditioner = (double *)blockOfShifts;

      /* Pass the original value of touch and update touch as the maximum  */
      /* value that takes for all inner_solve calls                        */
      int touch1 = touch0;

      CHKERR(inner_solve_Sprimme(blockSize, x, ldV, Bx, ldBV, r, ldW,
            blockNorms, evecs, ldevecs, Mfact, ipivot, xKinvBx, LprojectorQ,
            ldLprojectorQ, LprojectorX, ldLprojectorX, LprojectorBQ,
            ldLprojectorBQ, LprojectorBX, ldLprojectorBX, RprojectorQ,
            ldRprojectorQ, RprojectorX, ldRprojectorX, sizeLprojectorQ,
            sizeLprojectorX, sizeRprojectorQ, sizeRprojectorX, sol, ldsol,
            blockRitzVals, blockOfShifts, &touch1, startTime, ctx));
      *touch = max(*touch, touch1);

      Num_copy_matrix_Sprimme(sol, primme->nLocal, blockSize, ldsol,
            &V[ldV * basisSize], ldV, ctx);


      CHKERR(Num_free_SHprimme(xKinvBx, ctx));
      CHKERR(Num_free_Sprimme(sol, ctx));
      CHKERR(Num_free_Sprimme(KinvBx, ctx));
      CHKERR(KIND(Num_free_RHprimme, Num_free_SHprimme)(blockRitzVals, ctx));
#else
      CHKERR(PRIMME_FUNCTION_UNAVAILABLE);
#endif /* USE_HERMITIAN */
   } /* JDqmr variants */

   if (sortedRitzVals != ritzVals)
      CHKERR(KIND(Num_free_RHprimme, Num_free_SHprimme)(sortedRitzVals, ctx));
   CHKERR(KIND(Num_free_dprimme, Num_free_zprimme)(blockOfShifts, ctx));
   CHKERR(Num_free_RHprimme(approxOlsenEps, ctx));
   if (ilev != iev) CHKERR(Num_free_iprimme(ilev, ctx));

   return 0;

}
      

/*******************************************************************************
 * Subroutine computeRobustShift - This function computes the robust shift
 *    to be used in the correction equation.  The standard shift is the current
 *    Ritz value.  The calling routine eventually adds the robust shift, 
 *    epsilon, to the standard shift to accelerate the convergence of the outer 
 *    method.  
 *
 *    If residual vectors norms are compute with approximate eigenvectors with
 *    B-norm 1, then the bounds used to estimate the distance between the exact
 *    and the approximate eigenvalue is as follows [1]:
 *
 *    |val - appr. val| <= |r|_inv(B) <= sqrt(|inv(B)|) * |r|_2
 *
 *    |val - appr. val| <= (|r|_inv(B))^2 / gap <= |inv(B)| * (|r|_2)^2 / gap
 *
 *    [1] http://www.netlib.org/utk/people/JackDongarra/etemplates/node181.html
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * blockIndex    Index of the block vector for which the correction will be 
 *               solved.
 *
 * resNorm       The residual norm of the current Ritz vector
 *
 * prevRitzVals  Sorted array of locked Ritz values and Ritz values from the 
 *               previous outer iteration.  Array size is numPrevRitzVals.
 *
 * numPrevRitzVals The number of values in prevRitzVals
 * 
 * sortedRitzVals  Sorted array of locked Ritz values and current Ritz values.
 *                 The array is of size numLocked+basisSize
 *
 * approxOlsenShift The shift to be used to modify r-shift x, for approx Olsen's
 *
 * numSorted     The size of sortedRitzVals.
 *
 * ilev          Array of size blockSize that maps targeted Ritz values to
 *               their position with the array sortedRitzVals.
 *
 * Return Value
 * ------------
 * double  The value to be added to the current shift.  The calling routine
 *         will compute shift := shift +/- epsilon.
 ******************************************************************************/

STATIC HREAL computeRobustShift(int blockIndex, double resNorm, 
   HREAL *prevRitzVals, int numPrevRitzVals, HREAL *sortedRitzVals, 
   HREAL *approxOlsenShift, int numSorted, int *ilev, primme_params *primme) {

   int sortedIndex;                 /* Index of the current Ritz value in */
                                    /* sortedRitzVals.                    */
   HREAL gap, lowerGap, upperGap;  /* Gaps between the current and       */
                                    /* neighboring Ritz Values.           */ 
   HREAL delta;                    /* The difference between the current */
                                    /* Ritz value and its value in the    */
                                    /* previous iteration.                */
   HREAL epsilon;                  /* The return value.  The calling     */
                                    /* routine will modify the shift by   */
                                    /* this amount.                       */    

   /* If this is the 1st outer iteration, there are no previous */
   /* Ritz values. Return the residual norm                     */

   if (primme->stats.numOuterIterations <= 1) {
      *approxOlsenShift = resNorm * sqrt(primme->stats.estimateInvBNorm);
      return resNorm * sqrt(primme->stats.estimateInvBNorm);
   }

   sortedIndex = ilev[blockIndex];

   /* Compute the gap when the first eigenvalue with respect to the */
   /* current basis is to be computed.                              */

   if (sortedIndex == 0 && numSorted >= 2) {
      lowerGap = HOST_MACHINE_MAX;
      gap = fabs(sortedRitzVals[1] - sortedRitzVals[0]);
   }
   else if (sortedIndex > 0 && numSorted >= 2 && sortedIndex+1 < numSorted) {

      /* Take the smaller of the two gaps if an interior eigenvalue is */
      /* targeted.                                                     */

      lowerGap = fabs(sortedRitzVals[sortedIndex] - 
                      sortedRitzVals[sortedIndex-1]);
      upperGap = fabs(sortedRitzVals[sortedIndex+1] - 
                      sortedRitzVals[sortedIndex]);
      gap = min(lowerGap, upperGap);
   }
   else {
      lowerGap = fabs(sortedRitzVals[sortedIndex] -
                      sortedRitzVals[sortedIndex-1]);
      gap = lowerGap;
   }
   
   /* Compute the change in a Ritz value between successive iterations */

   if (sortedIndex < numPrevRitzVals) {
      delta = fabs(prevRitzVals[sortedIndex] - sortedRitzVals[sortedIndex]);
   }
   else {
      delta = HOST_MACHINE_MAX;
   }

   /* Compute epsilon according to the theory of Davis and Kahan described */
   /* in The Symmetric Eigenvalue Problem by B.N. Parlett.                 */

   if (gap > resNorm) {
      epsilon = min(
            delta, min(resNorm * resNorm * primme->stats.estimateInvBNorm / gap,
                         lowerGap));
   }
   else {
      epsilon = min(resNorm * sqrt(primme->stats.estimateInvBNorm), lowerGap);
   }

   *approxOlsenShift = min(delta, epsilon);

   /* If the above is too large a shift set it to a milder shift */
   /* epsilon = min(delta, epsilon); */

   return epsilon;

}


/*******************************************************************************
 * Subroutine mergeSort -- This routine merges the sorted lockedEvals array and 
 *   the sorted ritzVals array into a single sorted array, sortedRitzVals.
 *   It is only called for extreme (largest, smallest) eigenvalues, not for 
 *   interior ones. Thus it does not cover the interior case. 
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * lockedEvals  Array of size numLocked.  Ritz values that have been locked.
 *
 * numLocked    The size of array lockedEvals.
 *
 * ritzVals     The current Ritz values.  They are the eigenvalues of the
 *              projection V'*A*V.
 *
 * flags        Array of flags indicating which Ritz values in ritzVals are
 *              converged.
 *
 * basisSize    The current size of the basis V.
 *
 * ilev         Maps the blockSize targeted Ritz values to their position
 *              within sortedRitzVals.
 * 
 * blockSize    The number of Ritz values targeted during the current outer
 *              iteration.
 *
 * primme       A structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * sortedRitzVals  The array lockedEvals and ritzVals are merged into.
 *
 ******************************************************************************/

STATIC void mergeSort(HREAL *lockedEvals, int numLocked, HREAL *ritzVals, 
   int *flags, int basisSize, HREAL *sortedRitzVals, int *ilev, int blockSize,
   primme_params *primme) {
   
   int count;         /* The number of Ritz values merged        */
   int eval, ritzVal; /* The index of the next locked Ritz value */
                      /*   and current Ritz values to be merged  */
   int blockIndex;    /* Counter used to index ilev              */

   count = 0;
   eval = 0;
   ritzVal = 0;
   blockIndex = 0;
   
   /* Continue merging until all elements have been merged */

   while (count < numLocked+basisSize) {

      /* I. The list of locked eigenvalues has been completely merged, or  */
      /* if ritzVals[ritzVal] is greater/less than lockedEvals[eval], then */
      /* merge ritzVals[ritzVal] into the list of sorted Ritz values.      */
      /*                                                                   */
      /* II. If the list of current Ritz values has been completely merged,*/
      /* or if lockedEvals[eval] greater/less than ritzVals[ritzVal], then */
      /* merge lockedEvals[eval] into the list of sorted Ritz values.      */

      if (eval >= numLocked ||
            (ritzVal < basisSize &&
                  ((primme->target == primme_largest &&
                         ritzVals[ritzVal] >= lockedEvals[eval]) ||
                        (primme->target == primme_smallest &&
                              ritzVals[ritzVal] <= lockedEvals[eval])))) {
         sortedRitzVals[count] = ritzVals[ritzVal];

         /* If the Ritz value just merged is unconverged and there is room */
         /* enough in the block to target it, then keep track of its index */
         /* in the sorted list.                                            */

         if (blockIndex < blockSize && flags[ritzVal] == UNCONVERGED) {
            ilev[blockIndex] = count;
            blockIndex++;
         }

         ritzVal++;
      }
      else if (ritzVal >= basisSize || (primme->target == primme_largest &&
                                       lockedEvals[eval] >= ritzVals[ritzVal])
                                   || (primme->target == primme_smallest &&
                                       lockedEvals[eval] <= ritzVals[ritzVal])
              )
      {
         sortedRitzVals[count] = lockedEvals[eval];
         eval++;
      } 

      count++;
   }

}

/*******************************************************************************
 * Subroutine Olsen_preconditioner_block - This subroutine applies the projected
 *    preconditioner to a block of blockSize vectors r by computing:
 *       (I - (K^{-1} B x_i x_i^T) / (x_i^T K^{-1} B x_i) ) K^{-1}r_i
 *
 * Input Parameters
 * ----------------
 * r          The vectors the preconditioner and projection will be applied to.
 *
 * blockSize  The number of vectors in r, x
 *
 * primme       Structure containing various solver parameters
 *
 * Output parameters
 * -----------------
 * x The current Ritz vector used in the projection. It also carries the result.
 *
 ******************************************************************************/

STATIC int Olsen_preconditioner_block(SCALAR *r, PRIMME_INT ldr, SCALAR *x,
      PRIMME_INT ldx, SCALAR *Bx, PRIMME_INT ldBx, int blockSize,
      primme_context ctx) {

   primme_params *primme = ctx.primme;

   /* Allocate workspaces */

   SCALAR *KinvBxr;     /* Contain K^{-1} * [Bx r] */
   PRIMME_INT ldKinvBxr = primme->ldOPs;
   CHKERR(Num_malloc_Sprimme(ldKinvBxr * blockSize * 2, &KinvBxr, ctx));
   SCALAR *KinvBx = KinvBxr, *Kinvr = &KinvBxr[ldKinvBxr * blockSize];
   HSCALAR *xKinvBx, *xKinvr;
   CHKERR(Num_malloc_SHprimme(blockSize * 2, &xKinvBx, ctx));
   xKinvr = xKinvBx + blockSize;

   /*------------------------------------------------------------------ */
   /* Compute K^{-1}Bx and K^{-1}r                                      */
   /*------------------------------------------------------------------ */

   CHKERR(applyPreconditioner_Sprimme(
         Bx, primme->nLocal, ldBx, KinvBx, ldKinvBxr, blockSize, ctx));

   CHKERR(applyPreconditioner_Sprimme(
         r, primme->nLocal, ldr, Kinvr, ldKinvBxr, blockSize, ctx));

   /*---------------------------------------------------------------------- */
   /* Compute local x'K^{-1}Bx and x'K^{-1}r for each vector                */
   /*---------------------------------------------------------------------- */

   CHKERR(Num_dist_dots_Sprimme(x, ldx, KinvBx, ldKinvBxr, primme->nLocal,
         blockSize, xKinvBx, ctx));

   CHKERR(Num_dist_dots_Sprimme(x, ldx, Kinvr, ldKinvBxr, primme->nLocal,
         blockSize, xKinvr, ctx));

   /*------------------------------------------------------------------*/
   /* Compute K^(-1)r  - ( xKinvr/xKinvBx ) K^(-1)r for each vector     */
   /*------------------------------------------------------------------*/

   int blockIndex;
   for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {
      CHKERR(Num_copy_matrix_Sprimme(&Kinvr[ldKinvBxr * blockIndex],
            primme->nLocal, 1, ldKinvBxr, &x[ldx * blockIndex], ldx, ctx));

      if (ABS(xKinvBx[blockIndex]) > 0.0) {
         HSCALAR alpha;
         alpha = -xKinvr[blockIndex] / xKinvBx[blockIndex];
         Num_axpy_Sprimme(primme->nLocal, alpha,
               &KinvBx[ldKinvBxr * blockIndex], 1, &x[ldx * blockIndex], 1,
               ctx);
      }
   }

   CHKERR(Num_free_Sprimme(KinvBxr, ctx));
   CHKERR(Num_free_SHprimme(xKinvBx, ctx));

   return 0;

}

/*******************************************************************************
 *   subroutine setup_JD_projectors()
 *
 *   Sets up whether to use any of the following projections in the 
 *   correction equation:
 *
 *  INPUT
 *  -----
 *   x                The Ritz vector
 *   Bx               B*x
 *   evecs            Converged locked eigenvectors (denoted as Q herein)
 *   Bevecs           B*evecs
 *   evecsHat         K^{-1}*B*evecs
 *   numLocked        Number of locked eigenvectors (if locking)
 *   numConverged     Number of converged e-vectors copied in evecs (no locking)
 *   primme           The main data structures that contains the choices for
 *
 *       primme->LeftQ  : evecs in the left projector
 *       primme->LeftX  : x in the left projector
 *       primme->RightQ : evecs in the right projector
 *       primme->RightX : x in the right projector
 *       primme->SkewQ  : evecs in the right one, should be skewed
 *       primme->SkewX  : x in the right one, should be skewed
 *                   
 *  OUTPUT
 *  ------
 *  *KinvBx           The result of K^{-1}Bx (if needed, otherwise NULL)
 * **LprojectorQ      Pointer to the left projector for Q (could be NULL)
 * **LprojectorX      Pointer to the left projector for x (could be NULL)
 * **LprojectorBQ     Pointer to the left projector for BQ (could be NULL)
 * **LprojectorBX     Pointer to the left projector for Bx (could be NULL)
 * **RprojectorQ      Pointer to the right projector for Q (could be NULL)
 * **RprojectorX      Pointer to the right projector for X (could be NULL)
 *   sizeLprojectorQ  Size of the left projector Q(numConverged/numLocked)/+1 or 0
 *   sizeLprojectorX  Size of the left projector X (blockSize or 0)
 *   sizeRprojectorQ  Size of the Q right projector (numConverged/numLocked or 0)
 *   sizeRprojectorX  Size of the X right projector (blockSize or 0)
 *
 * ============================================================================
 * Functionality:
 *
 *   Q = locked vectors
 *   x = Ritz vector
 *   K = preconditioner
 *   r = a vector to apply the 
 *
 *   (I-BQQ')(I-Bxx)(A-shift B)(I-KBx(x'KBx)^(-1)x')(I-KBQ(Q'KBQ)^(-1)Q') K 
 *    -----  -----               -------------------  -------------------
 *     Lq     Lx                   Rx  Sx              Rq Sq
 *
 * Control parameters from primme (shortened here, see above for full names):
 *
 *     Lq = 1/0 whether to apply the left projector I-QQ'
 *     Lx = 1/0 whether to apply the left projector I-xx'
 *     Rx = 1/0 whether to apply the right projector with x (skew or orthogonal)
 *     Rq = 1/0 whether to apply the right projector with Q (skew or orthogonal)
 *     if (Rx == 1)
 *        Sx = 1 applies the right skew projector (I-KBx(x'KBx)^(-1)x')
 *        Sx = 0 applies the right orthogonal projector (I-Bxx')
 *     if (Rq == 1)
 *        Sq = 1 applies the right skew projector (I-KQ(Q'K Q)^(-1)Q') 
 *        Sq = 0 applies the right orthogonal projector (I-BQQ')
 *
 * Examples
 * Lq,Lx,Rx,Sx,Rq,Sq
 *  1  1  1  1  1  1  (the above equation)                   full, expensive JD
 *  0  0  1  1  1  1  (I-KBx(x'KBx)^(-1)x')(I-KBQ(Q'KBQ)^(-1)Q')        classic JD
 *  0  1  1  1  0  0  (I-Bxx')(A-shift B)(I-KBx(x'KBx)^(-1)x')     JD w/o deflaton
 *  0  0  0  0  0  0  (A-shift B) K                                  classic GD
 *  1  1  1  1  0  0  (I-BQQ')(I-Bxx')(A-shift B)(I-KBx(x'KBx)^(-1)x')  RECOMMENDED
 *  1  1  0  0  0  0  (I-BQQ')(I-Bxx')(A-shift B)                         similar
 *      Others...
 *                    Researchers can experiment with other projection schemes,
 *                    although our experience says they are rarely beneficial
 *
 * The left orthogonal projector for x and Q can be performed as one block
 * containing [Q x]. However, the right projections (if either is skew)
 * are performed separately for Q and x. There are memory reasons for 
 * doing so, but also we do not have to factor (Q'KBQ) at every outer step;
 * only when an eval converges. 
 *
 ******************************************************************************/

STATIC int setup_JD_projectors(SCALAR *x, PRIMME_INT ldx, SCALAR *Bx,
      PRIMME_INT ldBx, SCALAR *evecs, PRIMME_INT ldevecs, SCALAR *Bevecs,
      PRIMME_INT ldBevecs, SCALAR *evecsHat, PRIMME_INT ldevecsHat,
      SCALAR *KinvBx, PRIMME_INT ldKinvBx, HSCALAR *xKinvBx,
      SCALAR **LprojectorQ, PRIMME_INT *ldLprojectorQ, SCALAR **LprojectorX,
      PRIMME_INT *ldLprojectorX, SCALAR **LprojectorBQ,
      PRIMME_INT *ldLprojectorBQ, SCALAR **LprojectorBX,
      PRIMME_INT *ldLprojectorBX, SCALAR **RprojectorQ,
      PRIMME_INT *ldRprojectorQ, SCALAR **RprojectorX,
      PRIMME_INT *ldRprojectorX, int *sizeLprojectorQ, int *sizeLprojectorX,
      int *sizeRprojectorQ, int *sizeRprojectorX, int numLocked,
      int numConverged, int blockSize, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int n, sizeEvecs;

   *sizeLprojectorQ = 0;
   *sizeLprojectorX = 0;
   *sizeRprojectorQ = 0;
   *sizeRprojectorX = 0;
   *ldLprojectorQ = 0;
   *ldLprojectorX = 0;
   *ldLprojectorBQ = 0;
   *ldLprojectorBX = 0;
   *ldRprojectorQ = 0;
   *ldRprojectorX = 0;
   *LprojectorQ = NULL;
   *LprojectorX = NULL;
   *LprojectorBQ = NULL;
   *LprojectorBX = NULL;
   *RprojectorQ = NULL;
   *RprojectorX = NULL;

   n = primme->nLocal;
   if (primme->locking) 
      sizeEvecs = primme->numOrthoConst+numLocked;
   else
      sizeEvecs = primme->numOrthoConst+numConverged;
   
   /* ---------------------------------------------------------------------------- */
   /* Set up the left projector arrays. Make x adjacent to Q and Bx adjacent to BQ */
   /* ---------------------------------------------------------------------------- */
   
   if (primme->correctionParams.projectors.LeftQ) {
   
      *sizeLprojectorQ = sizeEvecs;
      *LprojectorQ = evecs;
      *ldLprojectorQ = ldevecs;
      *LprojectorBQ = Bevecs;
      *ldLprojectorBQ = ldBevecs;
      if (primme->correctionParams.projectors.LeftX) {
         if (blockSize <= 1) {
            CHKERR(Num_copy_matrix_Sprimme(x, n, blockSize, ldx,
                  &evecs[sizeEvecs * ldevecs], ldevecs, ctx));
            if (Bx != x) {
               CHKERR(Num_copy_matrix_Sprimme(Bx, n, blockSize, ldBx,
                     &Bevecs[sizeEvecs * ldBevecs], ldBevecs, ctx));
            }
            *sizeLprojectorQ += blockSize;
         }
         else {
            *sizeLprojectorX = blockSize;
            *LprojectorX = x;
            *ldLprojectorX = ldx;
            *LprojectorBX = Bx;
            *ldLprojectorBX = ldBx;
         }
      }
   }
   else {
      if (primme->correctionParams.projectors.LeftX) {
         *LprojectorX = x;
         *ldLprojectorX = ldx;
         *LprojectorBX = Bx;
         *ldLprojectorBX = ldBx;
         *sizeLprojectorX = blockSize;
      }
   }
      
   /* --------------------------------------------------------*/
   /* Set up the right projector arrays. Q and x separately   */
   /* --------------------------------------------------------*/
   
   /* ------------*/
   /* First for Q */
   /* ------------*/
   if (primme->correctionParams.projectors.RightQ) {
   
      if (primme->correctionParams.precondition    &&
          primme->correctionParams.projectors.SkewQ) {
         *RprojectorQ = evecsHat;       /* Use the K^(-1)*B*evecs array */
         *ldRprojectorQ = ldevecsHat;
      }                               
      else {  /* Right Q but not SkewQ */
         *RprojectorQ = Bevecs;          /* Use just the Bevecs array. */
         *ldRprojectorQ = ldBevecs;
      }
      *sizeRprojectorQ = sizeEvecs;
   
   } 
   else { /* if no RightQ projector */
   
      *RprojectorQ = NULL;       
      *sizeRprojectorQ = 0;
   }
   
   /* ------------*/
   /* Then for x  */
   /* ------------*/
   if (primme->correctionParams.projectors.RightX) {
   
      if (primme->correctionParams.precondition   &&
          primme->correctionParams.projectors.SkewX) {
         CHKERR(applyPreconditioner_Sprimme(
               Bx, primme->nLocal, ldBx, KinvBx, ldKinvBx, blockSize, ctx));
         *RprojectorX  = KinvBx;
         *ldRprojectorX  = ldKinvBx;
         CHKERR(Num_dist_dots_Sprimme(x, ldx, KinvBx, ldKinvBx, primme->nLocal,
               blockSize, xKinvBx, ctx));
      }      
      else {
         *RprojectorX = Bx;
         *ldRprojectorX  = ldBx;
         int i;
         for (i=0; i<blockSize; i++) xKinvBx[i] = 1.0;
      }
      *sizeRprojectorX = blockSize;
   }
   else {
      *RprojectorX = NULL;
      *sizeRprojectorX = 0;
   }

   return 0;

} /* setup_JD_projectors */

#endif /* SUPPORTED_TYPE */
