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
 * File: primme_c.c
 *
 * Purpose - Real, SCALAR precision front end to the multimethod eigensolver
 *
 * For the moment please cite the following two papers: 
 *
 *  A. Stathopoulos, Nearly optimal preconditioned methods for hermitian
 *    eigenproblems under limited memory. Part I: Seeking one eigenvalue,
 *    Tech Report: WM-CS-2005-03, July, 2005. To appear in SISC.
 *  A. Stathopoulos and J. R. McCombs, Nearly optimal preconditioned methods
 *    for hermitian eigenproblems under limited memory. Part II: Seeking many
 *    eigenvalues, Tech Report: WM-CS-2006-02, June, 2006.
 *
 * Additional information on the algorithms appears in:
 *
 *  J. R. McCombs and A. Stathopoulos, Iterative Validation of Eigensolvers:
 *    A Scheme for Improving the Reliability of Hermitian Eigenvalue Solvers
 *    Tech Report: WM-CS-2005-02, July, 2005, to appear in SISC.
 *  A. Stathopoulos, Locking issues for finding a large number of eigenvectors 
 *    of hermitian matrices, Tech Report: WM-CS-2005-03, July, 2005.
 *
 *   Some papers to be listed here. For the moment contact the author:
 *
 *                       andreas@cs.wm.edu
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/primme_c.c"
#endif

#include "numerical.h"
#include "template_normal.h"
#include "common_eigs.h"
#include "primme_interface.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "primme_c.h"
#include "main_iter.h"
#include "auxiliary_eigs.h"
#endif

/*******************************************************************************
 * Subroutine Xprimme - This routine is a front end used to perform 
 *    error checking on the input parameters, perform validation, 
 *    and make the call to main_iter. 
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * evals  Contains the converged Ritz values upon return.   Should be of size 
 *        primme->numEvals.
 * 
 * evecs  The local portions of the converged Ritz vectors.  The dimension of
 *        the array is at least primme->nLocal*primme->numEvals
 *
 * resNorms  The residual norms of the converged Ritz vectors.  Should be of 
 *           size primme->numEvals
 *  
 * primme  Structure containing various solver parameters and statistics
 *         See readme.txt for INPUT/OUTPUT variables in primme
 *
 * Return Value
 * ------------
 *  0 - Success
 *  1 - Reporting only memory space required
 * -1 - Failure to allocate workspace
 * -2 - Malloc failure in allocating a permutation integer array
 * -3 - main_iter encountered a problem
 * -4 ...-32 - Invalid input (parameters or primme struct) returned 
 *             by check_input()
 *
 ******************************************************************************/

int Xprimme(XEVAL *evals, XSCALAR *evecs, XREAL *resNorms,
            primme_params *primme) {

   return Xprimme_aux_Sprimme((void *)evals, (void *)evecs, (void *)resNorms, primme,
         PRIMME_OP_SCALAR);
}

// Definition for *hsprimme, *ksprimme, and *kcprimme

#if defined(USE_HALF) || defined(USE_HALFCOMPLEX) ||                           \
      defined(USE_HALF_MAGMA) || defined(USE_HALFCOMPLEX_MAGMA) ||             \
      defined(USE_HALF_CUBLAS) || defined(USE_HALFCOMPLEX_CUBLAS)

#  ifdef USE_HERMITIAN
      // Expand the terms {,magma_}{hs,ks}primme
#     define Xsprimme  CONCAT(CONCAT(STEM,USE_ARITH(hs,ks)),primme)
#  else
      // Expand the terms {,magma_}kcprimme_normal
#     define Xsprimme  WITH_KIND(CONCAT(CONCAT(STEM,kc),primme))
#  endif

int Xsprimme(KIND(float, PRIMME_COMPLEX_FLOAT) * evals, XSCALAR *evecs,
      float *resNorms, primme_params *primme) {

   return Xprimme_aux_Sprimme((void *)evals, (void *)evecs, (void *)resNorms, primme,
         primme_op_float);
}

#  undef Xsprimme
#endif

/*******************************************************************************
 * Subroutine Xprimme_aux - set defaults depending on the callee's type, and
 *    call wrapper_Sprimme with type set in internalPrecision. 
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * evals  Contains the converged Ritz values upon return.   Should be of size 
 *        primme->numEvals.
 * 
 * evecs  The local portions of the converged Ritz vectors.  The dimension of
 *        the array is at least primme->nLocal*primme->numEvals
 *
 * resNorms  The residual norms of the converged Ritz vectors.  Should be of 
 *           size primme->numEvals
 *  
 * primme  the PRIMME parameters
 *
 * evals_resNorms_type The type of the arrays evals and resNorsm.
 *
 * Return Value
 * ------------
 * return  error code
 ******************************************************************************/

TEMPLATE_PLEASE
int Xprimme_aux_Sprimme(void *evals, void *evecs, void *resNorms,
            primme_params *primme, primme_op_datatype evals_resNorms_type) {

#ifdef SUPPORTED_TYPE

   /* Generate context */

   primme_context ctx = primme_get_context(primme);

   /* Set the current type as the default type for user's operators */

   if (primme->matrixMatvec && primme->matrixMatvec_type == primme_op_default)
      primme->matrixMatvec_type = PRIMME_OP_SCALAR;
   if (primme->massMatrixMatvec && primme->massMatrixMatvec_type == primme_op_default)
      primme->massMatrixMatvec_type = PRIMME_OP_SCALAR;
   if (primme->applyPreconditioner && primme->applyPreconditioner_type == primme_op_default)
      primme->applyPreconditioner_type = PRIMME_OP_SCALAR;
   if (primme->globalSumReal && primme->globalSumReal_type == primme_op_default)
      primme->globalSumReal_type = PRIMME_OP_SCALAR;
   if (primme->broadcastReal && primme->broadcastReal_type == primme_op_default)
      primme->broadcastReal_type = PRIMME_OP_SCALAR;
   if (primme->convTestFun && primme->convTestFun_type == primme_op_default)
      primme->convTestFun_type = PRIMME_OP_SCALAR;
   if (primme->monitorFun && primme->monitorFun_type == primme_op_default)
      primme->monitorFun_type = PRIMME_OP_SCALAR;

   /* Number of returned eigenpairs */

   int outInitSize = 0;

   /* call primme for the internal working precision */

   int ret;
   primme_op_datatype t = primme->internalPrecision;
   if (t == primme_op_default) t = PRIMME_OP_SCALAR;
   switch (t) {
#  ifdef SUPPORTED_HALF_TYPE
   case primme_op_half:
      CHKERRVAL(wrapper_Shprimme(evals, evecs, resNorms, evals_resNorms_type,
                      PRIMME_OP_SCALAR, &outInitSize, ctx),
            &ret);
      break;
#  endif
#  ifndef PRIMME_WITHOUT_FLOAT
   case primme_op_float:
      CHKERRVAL(wrapper_Ssprimme(evals, evecs, resNorms, evals_resNorms_type,
                      PRIMME_OP_SCALAR, &outInitSize, ctx),
            &ret);
      break;
#  endif
   case primme_op_double:
      CHKERRVAL(wrapper_Sdprimme(evals, evecs, resNorms, evals_resNorms_type,
                      PRIMME_OP_SCALAR, &outInitSize, ctx),
            &ret);
      break;
#  ifdef PRIMME_WITH_NATIVE_QUAD
   case primme_op_quad:
      CHKERRVAL(wrapper_Sqprimme(evals, evecs, resNorms, evals_resNorms_type,
                      PRIMME_OP_SCALAR, &outInitSize, ctx),
            &ret);
      break;
#  endif
   default: ret = PRIMME_FUNCTION_UNAVAILABLE;
   }

   /* Free context */

   primme_free_context(ctx);

   /* Set the number of returned eigenpairs */

   primme->initSize = outInitSize;

   return ret;
#else

   (void)evals;
   (void)evecs;
   (void)resNorms;
   (void)evals_resNorms_type;

   primme->initSize = 0;
   return PRIMME_FUNCTION_UNAVAILABLE;
#endif /* SUPPORTED_TYPE */
}


#ifdef SUPPORTED_TYPE

/*******************************************************************************
 * Subroutine wrapper_Sprimme - Perform error checking on the input parameters,
 *    and make the call to main_iter. 
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * evals  Contains the converged Ritz values upon return.   Should be of size 
 *        primme->numEvals.
 * 
 * evecs  The local portions of the converged Ritz vectors.  The dimension of
 *        the array is at least primme->nLocal*primme->numEvals
 *
 * resNorms  The residual norms of the converged Ritz vectors.  Should be of 
 *           size primme->numEvals
 *
 * evals_resNorms_type The type of the arrays evals and resNorsm.
 *
 * evecs_type The type of the array evecs
 *
 * outInitSize The number of columns returned back.
 * 
 * ctx    primme context
 *
 * Return Value
 * ------------
 * return  error code
 ******************************************************************************/


TEMPLATE_PLEASE
int wrapper_Sprimme(void *evals, void *evecs, void *resNorms,
      primme_op_datatype evals_resNorms_type, primme_op_datatype evecs_type,
      int *outInitSize, primme_context ctx) {

   primme_params *primme = ctx.primme;

   /* In case of error, return initSize = 0 */

   *outInitSize = 0;

   /* zero out the timer */

   double t0 = primme_wTimer();

   /* Set some defaults for sequential programs */

   if (primme->numProcs <= 1 && evals != NULL && evecs != NULL &&
         resNorms != NULL) {
      primme->nLocal = primme->n;
      primme->procID = 0;
   }

   /* Set some defaults  */

   primme_set_defaults(primme);

   /* Observed orthogonality issues finding the largest/smallest values in  */
   /* single precision. Computing V'*B*V and solving the projected problem  */
   /* V'AVx = V'BVxl mitigates the problem.                                 */
   /* Also if maxBlockSize > 1, the code uses the block orthogonalization   */
   /* instead of Gram-Schmidt. But the current block orthogonalization does */
   /* not produce a machine precision orthonormal basis, and we deal with   */
   /* this by computing V'*B*V also in this case.                           */


   if (primme->orth == primme_orth_default) {
      if (PRIMME_OP_SCALAR <= primme_op_float || primme->maxBlockSize > 1) {
         primme->orth = primme_orth_explicit_I;
      } else {
         primme->orth = primme_orth_implicit_I;
      }
   }

   /* If we are free to choose the leading dimension of V and W, use    */
   /* a multiple of PRIMME_BLOCK_SIZE. This may improve the performance */
   /* of Num_update_VWXR_Sprimme and Bortho_block_Sprimme.              */

   if (primme->ldOPs == -1) {
      if (PRIMME_BLOCK_SIZE < INT_MAX) {
         primme->ldOPs = min(((primme->nLocal + PRIMME_BLOCK_SIZE - 1)
                  /PRIMME_BLOCK_SIZE)*PRIMME_BLOCK_SIZE, primme->nLocal);
      } else {
         primme->ldOPs = primme->nLocal;
      }
   }

   /* Deprecated input:                                              */

   if (evals == NULL && evecs == NULL && resNorms == NULL)
      return 0;

   /* Reset random number seed if inappropriate for DLARENV */
   /* Yields unique quadruples per proc if procID < 4096^3  */

   if (primme->iseed[0]<0 || primme->iseed[0]>4095) primme->iseed[0] = 
      primme->procID % 4096;
   if (primme->iseed[1]<0 || primme->iseed[1]>4095) primme->iseed[1] = 
      (int)(primme->procID/4096+1) % 4096;
   if (primme->iseed[2]<0 || primme->iseed[2]>4095) primme->iseed[2] = 
      (int)((primme->procID/4096)/4096+2) % 4096;
   if (primme->iseed[3]<0 || primme->iseed[3]>4095) primme->iseed[3] = 
      (2*(int)(((primme->procID/4096)/4096)/4096)+1) % 4096;

   /* Set default convTetFun  */

   if (!primme->convTestFun) {
      primme->convTestFun = convTestFunAbsolute;
      primme->convTestFun_type = PRIMME_OP_SCALAR;
      if (primme->eps == 0.0) {
         primme->eps = MACHINE_EPSILON * 1e4;
         /* The default value of eps is too much for half precision */
         if (primme->eps >= 1.0) primme->eps = 0.1;
      }
   }

   /* Set default monitor     */

   if (!primme->monitorFun) {
      primme->monitorFun = default_monitor;
      primme->monitorFun_type = PRIMME_OP_SCALAR;
   }

   /* Check primme input data for bounds, correct values etc. */

   CHKERR(coordinated_exit(check_params_coherence(ctx), ctx));
   CHKERR(check_input(evals, evecs, resNorms, primme))
       
   /* Cast evals, evecs and resNorms to working precision */

   HEVAL *evals0;
   HREAL *resNorms0;
   SCALAR *evecs0;
   CHKERR(KIND(Num_matrix_astype_RHprimme, Num_matrix_astype_SHprimme)(evals, 1,
         primme->numEvals, 1, evals_resNorms_type, (void **)&evals0, NULL,
         PRIMME_OP_HREAL, 1 /* alloc */, 0 /* not copy */, ctx));
   PRIMME_INT ldevecs0;
   CHKERR(Num_matrix_astype_Sprimme(evecs, primme->nLocal,
         primme->numOrthoConst + max(primme->numEvals, primme->initSize),
         primme->ldevecs, evecs_type, (void **)&evecs0, &ldevecs0,
         PRIMME_OP_SCALAR, 1 /* alloc */,
         primme->numOrthoConst + primme->initSize > 0 ? 1 : 0 /* copy? */,
         ctx));
   CHKERR(Num_matrix_astype_RHprimme(resNorms, 1, primme->numEvals,
         1, evals_resNorms_type, (void **)&resNorms0, NULL, PRIMME_OP_HREAL,
         1 /* alloc */, 0 /* not copy */, ctx));

   /* Call the solver */

   int ret, numRet;
   CHKERR(coordinated_exit(main_iter_Sprimme(evals0, evecs0, ldevecs0,
                                 resNorms0, t0, &ret, &numRet, ctx),
         ctx));

   /* Copy back evals, evecs and resNorms */

   CHKERR(KIND(Num_matrix_astype_RHprimme, Num_matrix_astype_SHprimme)(evals0,
         1, numRet, 1, PRIMME_OP_HREAL, (void **)&evals, NULL,
         evals_resNorms_type, -1 /* destroy */, 1 /* copy */, ctx));
   CHKERR(Num_copy_matrix_astype_Sprimme(evecs0, 0, primme->numOrthoConst,
         primme->nLocal, numRet, ldevecs0, PRIMME_OP_SCALAR, evecs, 0,
         primme->numOrthoConst, primme->ldevecs, evecs_type, ctx));
   if (evecs != evecs0) {
      CHKERR(Num_free_Sprimme(evecs0, ctx));
   }
   CHKERR(Num_matrix_astype_RHprimme(resNorms0, 1, numRet, 1, PRIMME_OP_HREAL,
         (void **)&resNorms, NULL, evals_resNorms_type, -1 /* destroy */,
         1 /* copy */, ctx));

   /* If no error, return initSize */

   *outInitSize = primme->initSize;

   primme->stats.elapsedTime = primme_wTimer() - t0;
   return ret;
}


/******************************************************************************
 * Subroutine check_input - checks the value of the input arrays, evals,
 *    evecs, and resNorms and the values of primme_params.
 *
 * INPUT
 * -----
 *  evals, evecs, resNorms   Output arrays for primme
 *  primme                   the main structure of parameters 
 *
 * return value -   0    If input parameters in primme are appropriate
 *                 <0    Inappropriate input parameters were found
 *
 ******************************************************************************/
STATIC int check_input(
      void *evals, void *evecs, void *resNorms, primme_params *primme) {
   int ret;
   ret = 0;

   if (primme == NULL)
      ret = -4;
   else if (primme->n < 0 || primme->nLocal < 0 || primme->nLocal > primme->n) 
      ret = -5;
   else if (primme->numProcs < 1)
      ret = -6;
   else if (primme->matrixMatvec == NULL) 
      ret = -7;
   else if (primme->applyPreconditioner == NULL && 
            primme->correctionParams.precondition > 0 ) 
      ret = -8;
   else if (primme->numEvals > primme->n)
      ret = -10;
   else if (primme->numEvals < 0)
      ret = -11;
   else if (primme->convTestFun != NULL && fabs(primme->eps) != 0.0L &&
            primme->eps < MACHINE_EPSILON)
      ret = -12;
   else if ( primme->target != primme_smallest  &&
             primme->target != primme_largest  &&
             primme->target != primme_largest_abs  &&
             primme->target != primme_closest_geq  &&
             primme->target != primme_closest_leq  &&
             primme->target != primme_closest_abs    )
      ret = -13;
   else if (primme->numOrthoConst < 0 || primme->numOrthoConst > primme->n)
      ret = -16;
   else if (primme->maxBasisSize < 2 && primme->n > 2) 
      ret = -17;
   else if (primme->minRestartSize < 0 || (primme->minRestartSize == 0
                                    && primme->n > 2 && primme->numEvals > 0))
      ret = -18;
   else if (primme->maxBlockSize < 0
             || (primme->maxBlockSize == 0 && primme->numEvals > 0)) 
      ret = -19;
   else if (primme->restartingParams.maxPrevRetain < 0)
      ret = -20;
   else if (primme->initSize < 0) 
      ret = -22;
   else if (primme->locking == 0 && primme->initSize > primme->maxBasisSize)
      ret = -23;
   else if (primme->locking > 0 && primme->initSize > primme->numEvals)
      ret = -24;
   else if (primme->minRestartSize + primme->restartingParams.maxPrevRetain 
                   >= primme->maxBasisSize && primme->n > primme->maxBasisSize)
      ret = -25;
   else if (primme->minRestartSize > primme->n && primme->n > 2)
      ret = -26;
   else if (primme->printLevel < 0 || primme->printLevel > 5)
      ret = -27; 
   else if (primme->correctionParams.convTest != primme_full_LTolerance &&
            primme->correctionParams.convTest != primme_decreasing_LTolerance &&
            primme->correctionParams.convTest != primme_adaptive_ETolerance &&
            primme->correctionParams.convTest != primme_adaptive )
      ret = -28;
   else if (primme->correctionParams.convTest == primme_decreasing_LTolerance &&
            primme->correctionParams.relTolBase <= 1.0L ) 
      ret = -29;
   else if (evals == NULL)
      ret = -30;
   else if (evecs == NULL || Num_check_pointer_Sprimme(evecs))
      ret = -31;
   else if (resNorms == NULL)
      ret = -32;
   else if (primme->locking == 0 && primme->minRestartSize < primme->numEvals &&
            primme->n > 2)
      ret = -33;
   else if (primme->ldevecs < primme->nLocal)
      ret = -34;
   else if (primme->ldOPs != 0 && primme->ldOPs < primme->nLocal)
      ret = -35;
   else if (primme->locking == 0
         && (primme->target == primme_closest_leq
            || primme->target == primme_closest_geq))
      ret = -38;
   else if (primme->massMatrixMatvec &&
            primme->projectionParams.projection != primme_proj_RR)
      ret = -39;
   /* Please keep this if instruction at the end */
   else if ( primme->target == primme_largest_abs ||
             primme->target == primme_closest_geq ||
             primme->target == primme_closest_leq ||
             primme->target == primme_closest_abs   ) {
      if (primme->numTargetShifts <= 0) {
         ret = -14;
      }
      else if (primme->targetShifts == NULL ) {
         ret = -15;
      }
   }

   return ret;
}

/*******************************************************************************
 * Subroutine convTestFunAbsolute - This routine implements primme_params.
 *    convTestFun and return an approximate eigenpair converged when           
 *    resNorm < eps*|A| for standard problems, and
 *    resNorm < (|A| + max(|\lambda|)*|B|)*eps for generalized problems
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * evec         The approximate eigenvector
 * eval         The approximate eigenvalue 
 * rNorm        The norm of the residual vector
 * primme       Structure containing various solver parameters
 *
 * OUTPUT PARAMETERS
 * ----------------------------------
 * isConv      if it isn't zero the approximate pair is marked as converged
 ******************************************************************************/

STATIC void convTestFunAbsolute(double *eval, void *evec, double *rNorm,
      int *isConv, primme_params *primme, int *ierr) {

   (void)eval; /* unused parameter */
   (void)evec; /* unused parameter */

   if (primme->massMatrixMatvec == NULL) {
      *isConv = *rNorm < max(primme->eps, MACHINE_EPSILON * 2) *
                               problemNorm_Sprimme(0, primme);
   }
   else {
      *isConv = *rNorm < max(primme->eps, MACHINE_EPSILON) *
                               problemNorm_Sprimme(0, primme);
   }
   *ierr = 0;
}

/*******************************************************************************
 * Subroutine default_monitor - report iterations, #MV, residual norm,
 *    eigenvalues, etc. at every inner/outer iteration and when some pair
 *    converges.       
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * basisEvals   The approximate eigenvalues of the basis
 * basisSize    The size of the basis
 * basisFlags   The state of every approximate pair of the basis (see conv_flags)
 * iblock       Indices of the approximate pairs in the block
 * blockSize    The size of the block
 * basisNorms   The approximate residual norms of the pairs of the basis
 * numConverged The number of pairs converged in the basis and the locked pairs
 *              (this value isn't monotonic!)
 * lockedEvals  The locked eigenvalues
 * numLocked    The number of pairs locked
 * lockedFlags  The state of each locked eigenpair (see conv_flags)
 * lockedNorms  The residual norms of the locked pairs
 * inner_its    The number of performed QMR iterations in the current correction equation
 * LSRes        The residual norm of the linear system at the current QMR iteration
 * event        The event reported
 * primme       Structure containing various solver parameters and statistics
 *
 * OUTPUT
 * ------
 * err          Error code
 * 
 ******************************************************************************/

STATIC void default_monitor(void *basisEvals_, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms_, int *numConverged,
      void *lockedEvals_, int *numLocked, int *lockedFlags, void *lockedNorms_,
      int *inner_its, void *LSRes_, const char *msg, double *time,
      primme_event *event, primme_params *primme, int *err) {

   (void)basisSize;
   (void)basisFlags;

   XEVAL *basisEvals = (XEVAL *)basisEvals_,
         *lockedEvals = (XEVAL *)lockedEvals_;
   XREAL *basisNorms = (XREAL *)basisNorms_,
         *lockedNorms = (XREAL *)lockedNorms_, *LSRes = (XREAL *)LSRes_;

   assert(event != NULL && primme != NULL);

   /* Only print report if this is proc zero or it is profiling */
   if (primme->outputFile &&
         (primme->procID == 0 || *event == primme_event_profile)) {
      switch(*event) {
      case primme_event_outer_iteration:
         assert(basisSize && (!*basisSize || (basisEvals && basisFlags)) &&
                blockSize && (!*blockSize || (iblock && basisNorms)) &&
                numConverged);
         assert(!primme->locking ||
                (numLocked && (!*numLocked || (lockedEvals && lockedFlags &&
                                                    lockedNorms))));
         if (primme->printLevel >= 3) {
            int i;  /* Loop variable */
            int found;  /* Reported eigenpairs found */

            if (primme->locking) 
               found = *numLocked;
            else 
               found = *numConverged;

            for (i=0; i < *blockSize; i++) {
               fprintf(primme->outputFile,
                     "OUT %" PRIMME_INT_P " conv %d blk %d MV %" PRIMME_INT_P
                     " Sec %E EV %13E " KIND(, "%13E i ") "|r| %.3E\n",
                     primme->stats.numOuterIterations, found, i,
                     primme->stats.numMatvecs, primme->stats.elapsedTime,
                     (double)EVAL_REAL_PART(basisEvals[iblock[i]]),
#ifndef USE_HERMITIAN
                     (double)EVAL_IMAGINARY_PART(basisEvals[iblock[i]]),
#endif

                     (double)basisNorms[iblock[i]]);
            }
         }
         break;
      case primme_event_inner_iteration:
         assert(basisSize && iblock && basisNorms && inner_its && LSRes);
         (void)inner_its;
         if (primme->printLevel >= 4) {
            fprintf(primme->outputFile,
                  "INN MV %" PRIMME_INT_P " Sec %e Eval %13E " KIND(
                        , "%13E i ") "Lin|r| %.3e EV|r| %.3e\n",
                  primme->stats.numMatvecs, primme->stats.elapsedTime,
                  (double)EVAL_REAL_PART(basisEvals[iblock[0]]),
#ifndef USE_HERMITIAN
                  (double)EVAL_IMAGINARY_PART(basisEvals[iblock[0]]),
#endif
                  (double)*LSRes, (double)basisNorms[iblock[0]]);
         }
        break;
      case primme_event_converged:
         assert(numConverged && iblock && basisEvals && basisNorms);
         if ((!primme->locking && primme->printLevel >= 2)
               || (primme->locking && primme->printLevel >= 5))
            fprintf(primme->outputFile,
                  "#Converged %d eval[ %d ]= %13E " KIND(,
                        "%13E i ") "norm %e Mvecs %" PRIMME_INT_P " Time %g\n",
                  *numConverged, iblock[0],
                  (double)EVAL_REAL_PART(basisEvals[iblock[0]]),
#ifndef USE_HERMITIAN
                  (double)EVAL_IMAGINARY_PART(basisEvals[iblock[0]]),
#endif
                  (double)basisNorms[iblock[0]], primme->stats.numMatvecs,
                  primme->stats.elapsedTime);
         break;
      case primme_event_locked:
         assert(numLocked && lockedEvals && lockedNorms && lockedFlags);
         if (primme->printLevel >= 2) {
            fprintf(primme->outputFile,
                  "Lock epair[ %d ]= %13E " KIND(
                        , "%13E i ") "norm %.4e Mvecs %" PRIMME_INT_P
                                     " Time %.4e Flag %d\n",
                  *numLocked - 1,
                  (double)EVAL_REAL_PART(lockedEvals[*numLocked - 1]),
#ifndef USE_HERMITIAN
                  (double)EVAL_IMAGINARY_PART(lockedEvals[*numLocked - 1]),
#endif
                  (double)lockedNorms[*numLocked - 1], primme->stats.numMatvecs,
                  primme->stats.elapsedTime, lockedFlags[*numLocked - 1]);
         }
         break;
      case primme_event_message:
         assert(msg != NULL);
         if (primme->printLevel >= 2) { 
            fprintf(primme->outputFile, 
                  "%s\n", msg);
         }
         break;
      case primme_event_profile:
         assert(msg != NULL && time != NULL);
         if (primme->printLevel >= 3 && *time < 0.0) { 
            fprintf(primme->outputFile, "entering in %s proc %d\n", msg, primme->procID);
         }
         if (primme->printLevel >= 2 && *time >= 0.0) { 
            fprintf(primme->outputFile, "time %g for %s proc %d\n", *time, msg, primme->procID);
         }
         break;
      default:
         break;
      }
      fflush(primme->outputFile);
   }
   *err = 0;
}

/******************************************************************************
 * check_params_coherence - check that all processes has the same values in
 *    critical parameters.
 *
 * INPUT
 * -----
 *  primme         the main structure of parameters 
 *
 * RETURN:
 *    error code
 ******************************************************************************/

STATIC int check_params_coherence(primme_context ctx) {
   primme_params *primme = ctx.primme;

   /* Check number of procs and procs with id zero */

   HREAL aux[2] = {(HREAL)1.0, (HREAL)(ctx.procID == 0 ? 1.0 : 0.0)};
   CHKERR(globalSum_RHprimme(aux, 2, ctx));
   CHKERRM((aux[0] > 1) != (ctx.numProcs > 1),
         -1, "numProcs does not match the actual number of processes");
   CHKERRM(aux[1] > 1, -1, "There is not a single process with ID zero");

   /* Check broadcast */

   HREAL val = 123, val0 = val;
   CHKERR(broadcast_RHprimme(&val, 1, ctx));
   CHKERRM(fabs(val - val0) > val0 * MACHINE_EPSILON * 1.3, -1,
         "broadcast function does not work properly");

   /* Check that all processes has the same value for the next params */

   PARALLEL_CHECK(primme->n);
   PARALLEL_CHECK(primme->numEvals);
   PARALLEL_CHECK(primme->target);
   PARALLEL_CHECK(primme->numTargetShifts);
   PARALLEL_CHECK(primme->dynamicMethodSwitch);
   PARALLEL_CHECK(primme->locking);
   PARALLEL_CHECK(primme->initSize);
   PARALLEL_CHECK(primme->numOrthoConst);
   PARALLEL_CHECK(primme->maxBasisSize);
   PARALLEL_CHECK(primme->minRestartSize);
   PARALLEL_CHECK(primme->maxBlockSize);
   PARALLEL_CHECK(primme->maxMatvecs);
   PARALLEL_CHECK(primme->maxOuterIterations);
   PARALLEL_CHECK(primme->aNorm);
   PARALLEL_CHECK(primme->BNorm);
   PARALLEL_CHECK(primme->invBNorm);
   PARALLEL_CHECK(primme->eps);
   PARALLEL_CHECK(primme->orth);
   PARALLEL_CHECK(primme->initBasisMode);
   PARALLEL_CHECK(primme->projectionParams.projection);
   PARALLEL_CHECK(primme->restartingParams.maxPrevRetain);
   PARALLEL_CHECK(primme->correctionParams.precondition);
   PARALLEL_CHECK(primme->correctionParams.robustShifts);
   PARALLEL_CHECK(primme->correctionParams.maxInnerIterations);
   PARALLEL_CHECK(primme->correctionParams.projectors.LeftQ);
   PARALLEL_CHECK(primme->correctionParams.projectors.LeftX);
   PARALLEL_CHECK(primme->correctionParams.projectors.RightQ);
   PARALLEL_CHECK(primme->correctionParams.projectors.RightX);
   PARALLEL_CHECK(primme->correctionParams.projectors.SkewQ);
   PARALLEL_CHECK(primme->correctionParams.projectors.SkewX);
   PARALLEL_CHECK(primme->correctionParams.convTest);
   PARALLEL_CHECK(primme->correctionParams.relTolBase);

   return 0;
}

/******************************************************************************
 * coordinated_exit - make sure that if main_iter returns error in some process,
 *    then all processes return an error.
 *
 * INPUT
 * -----
 *  primme         the main structure of parameters 
 *
 * RETURN:
 *    error code
 ******************************************************************************/
STATIC int coordinated_exit(int ret, primme_context ctx) {

   primme_params *primme = ctx.primme;

   if (ret != PRIMME_PARALLEL_FAILURE && primme->globalSumReal) {
      HREAL pret = (HREAL)(ret != 0 ? 1 : 0);
      int count = 1, ierr = 0;
      CHKERRM(
            (primme->globalSumReal(&pret, &pret, &count, primme, &ierr), ierr),
            PRIMME_USER_FAILURE, "Error returned by 'globalSumReal' %d", ierr);
      if (pret > 0.0) return ret ? ret : PRIMME_PARALLEL_FAILURE;
   }

   return ret;
}

#endif /* SUPPORTED_TYPE */
