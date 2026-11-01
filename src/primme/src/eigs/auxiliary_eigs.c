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
 * File: auxiliary_eigs.c
 *
 * Purpose - Miscellanea functions used by PRIMME EIGS not dependent on the
 *           Hermiticity/normality of the operator
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/auxiliary_eigs.c"
#endif

#include <string.h> /* memset */
#include "common_eigs.h"
#include "numerical.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "auxiliary_eigs.h"
#endif

#ifdef SUPPORTED_TYPE

#ifdef USE_DOUBLE

/******************************************************************************
 * Function monitor_report - pass to the monitor the reports
 *
 * PARAMETERS
 * ---------------------------
 * fun      function name or message to report
 * time     time spent on the call
 * ctx      primme context
 *
 ******************************************************************************/

static int monitor_report(const char *fun, double time, primme_context ctx) {
   if (ctx.primme && ctx.primme->monitorFun) {
      int err = 0;
      primme_event event =
            (time >= -.5 ? primme_event_profile : primme_event_message);

#ifdef PRIMME_PROFILE
      /* Avoid profiling this function. It will turn out in a recursive call */
      ctx.path = NULL;
#endif

      CHKERRM((ctx.primme->monitorFun(NULL, NULL, NULL, NULL, NULL,
                     NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, fun, &time, &event, ctx.primme, &err),
                    err),
            -1, "Error returned by monitorFun: %d", err);
   }
   return 0;
}

/******************************************************************************
 * Function primme_get_context - return a context from the primme_params
 *
 * PARAMETERS
 * ---------------------------
 * primme      primme_params struct
 *
 ******************************************************************************/

TEMPLATE_PLEASE
primme_context primme_get_context(primme_params *primme) {
   primme_context ctx;
   memset(&ctx, 0, sizeof(primme_context));
   if (primme) {
      ctx.primme = primme;
      ctx.printLevel = primme->printLevel;
      ctx.outputFile = primme->outputFile;
      ctx.numProcs = primme->numProcs;
      ctx.procID = primme->procID;
      ctx.mpicomm = primme->commInfo;
      ctx.globalSum = globalSum_Tprimme;
      ctx.bcast = broadcast_Tprimme; 
      ctx.queue = primme->queue;
      ctx.report = monitor_report;
#ifdef PRIMME_PROFILE
      if (primme->profile) {
         /* Compile regex. If there is no errors, set path to a nonzero       */
         /* value. Set ctx.report to the function that will channel the       */
         /* reports to the monitor. Report errors if they are.                */

         int ierr = regcomp(&ctx.profile, primme->profile, REG_NOSUB);
         if (ierr || MALLOC_PRIMME(1, &ctx.timeoff)) {
            char errmsg[100];
            regerror(ierr, &ctx.profile, errmsg, 100);
            if (ctx.report && ierr != 0) ctx.report(errmsg, -1, ctx);
            regfree(&ctx.profile);
            ctx.path = NULL;
         } else {
            *ctx.timeoff = 0.0;
            ctx.path = "";
         }
      } else {
         ctx.path = NULL;
      }
#endif
   }

   /* All error routines assume that there is frame. We push one here */

   Mem_push_frame(&ctx);

   return ctx;
} 

/******************************************************************************
 * Function primme_free_context - free memory associated to the context
 *
 * PARAMETERS
 * ---------------------------
 * ctx         context
 *
 ******************************************************************************/

TEMPLATE_PLEASE
void primme_free_context(primme_context ctx) {

   /* Pop frame pushed in primme_get_context */

   Mem_pop_frame(&ctx);

   /* Free profiler */

#ifdef PRIMME_PROFILE
   if (ctx.path) regfree(&ctx.profile);
   if (ctx.timeoff) free(ctx.timeoff);
#endif
}

#endif /* USE_DOUBLE */

/*******************************************************************************
 * Subroutine matrixMatvec_ - Computes A*V(:,nv+1) through A*V(:,nv+blksze)
 *           where V(:,nv+1:nv+blksze) are the new correction vectors.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V          The orthonormal basis
 * nLocal     Number of rows of each vector stored on this node
 * ldV        The leading dimension of V
 * ldW        The leading dimension of W
 * basisSize  Number of vectors in V
 * blockSize  The current block size
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * W          A*V
 ******************************************************************************/

TEMPLATE_PLEASE
int matrixMatvec_Sprimme(SCALAR *V, PRIMME_INT nLocal, PRIMME_INT ldV,
      SCALAR *W, PRIMME_INT ldW, int basisSize, int blockSize,
      primme_context ctx) {

   primme_params *primme = ctx.primme;

   if (blockSize <= 0)
      return 0;

   assert(ldV >= nLocal && ldW >= nLocal);
   assert(primme->ldOPs == 0 || primme->ldOPs >= nLocal);

   double t0 = primme_wTimer();

   /* Cast V and W */

   SCALAR *Vb = &V[ldV * basisSize], *Wb = &W[ldW * basisSize];
   void *V0, *W0;
   PRIMME_INT ldV0, ldW0;
   CHKERR(Num_matrix_astype_Sprimme(Vb, nLocal, blockSize, ldV,
         PRIMME_OP_SCALAR, &V0, &ldV0, primme->matrixMatvec_type, 1 /* alloc */,
         1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_Sprimme(Wb, nLocal, blockSize, ldW,
         PRIMME_OP_SCALAR, &W0, &ldW0, primme->matrixMatvec_type, 1 /* alloc */,
         0 /* no copy */, ctx));

   /* W(:,c) = A*V(:,c) for c = basisSize:basisSize+blockSize-1 */

   int ierr = 0;
   CHKERRM(
         (primme->matrixMatvec(V0, &ldV0, W0, &ldW0, &blockSize, primme, &ierr),
               ierr),
         PRIMME_USER_FAILURE, "Error returned by 'matrixMatvec' %d", ierr);

   /* Copy back W */

   CHKERR(Num_matrix_astype_Sprimme(W0, nLocal, blockSize, ldW0,
         primme->matrixMatvec_type, (void **)&Wb, &ldW,
         PRIMME_OP_SCALAR, 0 /* not alloc */, 1 /* copy */, ctx));

   if (Vb != V0) CHKERR(Num_free_Sprimme((SCALAR*)V0, ctx));
   if (Wb != W0) CHKERR(Num_free_Sprimme((SCALAR*)W0, ctx));

   primme->stats.timeMatvec += primme_wTimer() - t0;
   primme->stats.numMatvecs += blockSize;

   return 0;
}

/*******************************************************************************
 * Subroutine massMatrixMatvec - Computes B*V(:,nv+1:nv+blksze)
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V          The orthonormal basis
 * nLocal     Number of rows of each vector stored on this node
 * ldV        The leading dimension of V
 * ldW        The leading dimension of W
 * basisSize  Number of vectors in V
 * blockSize  The current block size
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * BV          B*V
 ******************************************************************************/

TEMPLATE_PLEASE
int massMatrixMatvec_Sprimme(SCALAR *V, PRIMME_INT nLocal, PRIMME_INT ldV,
      SCALAR *BV, PRIMME_INT ldBV, int basisSize, int blockSize,
      primme_context ctx) {

   primme_params *primme = ctx.primme;

   if (blockSize <= 0)
      return 0;

   assert(ldV >= nLocal && ldBV >= nLocal);
   assert(primme->ldOPs == 0 || primme->ldOPs >= nLocal);

   double t0 = primme_wTimer();

   /* Cast V and BV */

   SCALAR *Vb = &V[ldV * basisSize], *BVb = &BV[ldBV * basisSize];
   void *V0, *BV0;
   PRIMME_INT ldV0, ldBV0;
   CHKERR(Num_matrix_astype_Sprimme(Vb, nLocal, blockSize, ldV,
         PRIMME_OP_SCALAR, &V0, &ldV0, primme->massMatrixMatvec_type,
         1 /* alloc */, 1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_Sprimme(BVb, nLocal, blockSize, ldBV,
         PRIMME_OP_SCALAR, &BV0, &ldBV0, primme->massMatrixMatvec_type,
         1 /* alloc */, 0 /* no copy */, ctx));

   /* BV(:,c) = B*V(:,c) for c = basisSize:basisSize+blockSize-1 */

   int ierr = 0;
   CHKERRM((primme->massMatrixMatvec(
                  V0, &ldV0, BV0, &ldBV0, &blockSize, primme, &ierr),
                 ierr),
         PRIMME_USER_FAILURE, "Error returned by 'massMatrixMatvec' %d", ierr);

   /* Copy back BV */

   CHKERR(Num_matrix_astype_Sprimme(BV0, nLocal, blockSize, ldBV0,
         primme->matrixMatvec_type, (void **)&BVb, &ldBV, PRIMME_OP_SCALAR,
         0 /* not alloc */, 1 /* copy */, ctx));

   if (Vb != V0) CHKERR(Num_free_Sprimme((SCALAR*)V0, ctx));
   if (BVb != BV0) CHKERR(Num_free_Sprimme((SCALAR*)BV0, ctx));

   primme->stats.timeMatvec += primme_wTimer() - t0;
   primme->stats.numMatvecs += blockSize;

   return 0;
}


/*******************************************************************************
 * Subroutine applyPreconditioner - apply preconditioner to V
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V          The orthonormal basis
 * nLocal     Number of rows of each vector stored on this node
 * ldV        The leading dimension of V
 * ldW        The leading dimension of W
 * blockSize  The number of columns of V and W.
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * W          M*V
 ******************************************************************************/

TEMPLATE_PLEASE
int applyPreconditioner_Sprimme(SCALAR *V, PRIMME_INT nLocal, PRIMME_INT ldV,
      SCALAR *W, PRIMME_INT ldW, int blockSize, primme_context ctx) {

   primme_params *primme = ctx.primme;

   if (blockSize <= 0) return 0;
   assert(primme->nLocal == nLocal);

   double t0 = primme_wTimer();

   if (primme->correctionParams.precondition) {

      /* Cast V and W */

      void *V0, *W0;
      PRIMME_INT ldV0, ldW0;
      CHKERR(Num_matrix_astype_Sprimme(V, nLocal, blockSize, ldV,
            PRIMME_OP_SCALAR, &V0, &ldV0, primme->applyPreconditioner_type,
            1 /* alloc */, 1 /* copy */, ctx));
      CHKERR(Num_matrix_astype_Sprimme(W, nLocal, blockSize, ldW,
            PRIMME_OP_SCALAR, &W0, &ldW0, primme->applyPreconditioner_type,
            1 /* alloc */, 0 /* no copy */, ctx));

      /* Call user function */

      int ierr=0;
      CHKERRM((primme->applyPreconditioner(
                     V0, &ldV0, W0, &ldW0, &blockSize, primme, &ierr),
                    ierr),
            -1, "Error returned by 'applyPreconditioner' %d", ierr);
      primme->stats.numPreconds += blockSize;

      /* Copy back W and destroy cast matrices */

      if (V != V0) CHKERR(Num_free_Sprimme((SCALAR*)V0, ctx));
      CHKERR(Num_matrix_astype_Sprimme(W0, nLocal, blockSize, ldW0,
            primme->applyPreconditioner_type, (void **)&W, &ldW,
            PRIMME_OP_SCALAR, -1 /* destroy */, 1 /* copy */, ctx));

   }
   else {
      Num_copy_matrix_Sprimme(V, nLocal, blockSize, ldV, W, ldW, ctx);
   }

   primme->stats.timePrecond += primme_wTimer() - t0;

   return 0;
}

#ifdef USE_HOST

TEMPLATE_PLEASE
int globalSum_Sprimme(SCALAR *buffer, int count, primme_context ctx) {

#ifdef USE_COMPLEX
   count *= 2;
#endif

   return globalSum_Tprimme(buffer, PRIMME_OP_SCALAR, count, ctx);
}

TEMPLATE_PLEASE
int broadcast_Sprimme(SCALAR *buffer, int count, primme_context ctx) {

#ifdef USE_COMPLEX
   count *= 2;
#endif

   return broadcast_Tprimme(buffer, PRIMME_OP_SCALAR, count, ctx);
}

#ifdef USE_DOUBLE

TEMPLATE_PLEASE
int globalSum_Tprimme(
      void *buffer, primme_op_datatype buffert, int count, primme_context ctx) {

   primme_params *primme = ctx.primme;

   /* Quick exit */

   if (!primme || primme->numProcs == 1 || !primme->globalSumReal) {
      return 0;
   }


   double t0 = primme_wTimer();

   /* Cast buffer */

   void *buffer0 = NULL;
   CHKERR(Num_matrix_astype_Rprimme(buffer, 1, count, 1, buffert, &buffer0,
         NULL, primme->globalSumReal_type, 1 /* alloc */, 1 /* copy */, ctx));

   int ierr = 0;
   CHKERRM(
         (primme->globalSumReal(buffer0, buffer0, &count, primme, &ierr), ierr),
         PRIMME_USER_FAILURE, "Error returned by 'globalSumReal' %d", ierr);

   /* Copy back buffer0 */

   CHKERR(Num_matrix_astype_Rprimme(buffer0, 1, count, 1,
         primme->globalSumReal_type, (void **)&buffer, NULL, buffert,
         -1 /* dealloc */, 1 /* copy */, ctx));

   primme->stats.numGlobalSum++;
   primme->stats.timeGlobalSum += primme_wTimer() - t0;
   primme->stats.volumeGlobalSum += count;

   return 0;
}

TEMPLATE_PLEASE
int broadcast_Tprimme(
      void *buffer, primme_op_datatype buffert, int count, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int ierr;

   /* Quick exit */

   if (!primme || primme->numProcs == 1) {
      return 0;
   }

   double t0 = primme_wTimer();

   if (primme && primme->broadcastReal) {
      /* Cast buffer */

      void *buffer0 = NULL;
      CHKERR(Num_matrix_astype_dprimme(buffer, 1, count, 1, buffert,
            (void **)&buffer0, NULL, primme->broadcastReal_type, 1 /* alloc */,
            1 /* copy */, ctx));

      CHKERRM((primme->broadcastReal(buffer0, &count, primme, &ierr), ierr),
            PRIMME_USER_FAILURE, "Error returned by 'broadcastReal' %d", ierr);

      /* Copy back buffer0 */

      CHKERR(Num_matrix_astype_Sprimme(buffer0, 1, count, 1,
            primme->broadcastReal_type, (void **)&buffer, NULL, buffert,
            -1 /* dealloc */, 1 /* copy */, ctx));

   } else {
      if (primme->procID != 0) {
         CHKERR(Num_zero_matrix_Tprimme(buffer, buffert, 1, count, 1, ctx));
      }
      CHKERR(globalSum_Tprimme(buffer, buffert, count, ctx));
   }

   primme->stats.numBroadcast++;
   primme->stats.timeBroadcast += primme_wTimer() - t0;
   primme->stats.volumeBroadcast += count;

   return 0;
}

TEMPLATE_PLEASE
int broadcast_iprimme(int *buffer, int count, primme_context ctx) {

   CHKERR(broadcast_Tprimme(buffer, primme_op_int, count, ctx));
   return 0;
}

#endif /* USE_DOUBLE */

#endif /* USE_HOST */

/*******************************************************************************
 * Subroutine machineEpsMatrix - return the machine epsilon considering the
 *    precision of the matrixMatvec and massMatrixMatvec, and the working
 *    precision.
 * 
 * INPUT/OUTPUT PARAMETERS
 * -----------------------
 * eps      return machine epsilon
 *
 * ctx      primme context
 * 
 * OUTPUT
 * ------
 * return   error code
 ******************************************************************************/

TEMPLATE_PLEASE
int machineEpsMatrix_Sprimme(double *eps, primme_context ctx) {
   primme_params *primme = ctx.primme;

   double eps_matvec = MACHINE_EPSILON, eps_massmatvec = MACHINE_EPSILON;
   CHKERR(Num_machine_epsilon_Sprimme(primme->matrixMatvec_type, &eps_matvec));
   if (primme && primme->massMatrixMatvec) {
      CHKERR(Num_machine_epsilon_Sprimme(
            primme->massMatrixMatvec_type, &eps_massmatvec));
   }

   *eps = max(max(MACHINE_EPSILON, eps_matvec), eps_massmatvec);

   return 0;
}

/*******************************************************************************
 * Subroutine machineEpsOrth - return the machine epsilon considering the
 *    precision of the globalSumReal and broadcastReal, and the working
 *    precision.
 * 
 * INPUT/OUTPUT PARAMETERS
 * -----------------------
 * eps      return machine epsilon
 *
 * ctx      primme context
 * 
 * OUTPUT
 * ------
 * return   error code
 ******************************************************************************/

TEMPLATE_PLEASE
int machineEpsOrth_Sprimme(double *eps, primme_context ctx) {
   primme_params *primme = ctx.primme;

   double eps_globalsum = MACHINE_EPSILON, eps_broadcast = MACHINE_EPSILON;
   if (primme && primme->numProcs > 1) {
      CHKERR(Num_machine_epsilon_Sprimme(
            primme->globalSumReal_type, &eps_globalsum));
      if (primme->broadcastReal) {
         CHKERR(Num_machine_epsilon_Sprimme(
               primme->broadcastReal_type, &eps_broadcast));
      }
   }

   *eps = max(max(MACHINE_EPSILON, eps_globalsum), eps_broadcast);

   return 0;
}

/*******************************************************************************
 * Subroutine problemNorm - return an estimation of |B\A|
 * 
 * INPUT PARAMETERS
 * ----------------
 * overrideUserEstimations    if nonzero, use estimations of |A| and |inv(B)| if
 *                            they are larger than primme.aNorm and primme.invBNorm 
 * 
 * OUTPUT
 * ------
 * return                     estimation of |B\A|
 ******************************************************************************/

TEMPLATE_PLEASE
HREAL problemNorm_Sprimme(
      int overrideUserEstimations, struct primme_params *primme) {

   if (!overrideUserEstimations) {
      if (!primme->massMatrixMatvec) {
         return primme->aNorm > 0.0 ? primme->aNorm
                                    : primme->stats.estimateLargestSVal;
      } else {
         return primme->aNorm > 0.0 && primme->invBNorm > 0.0
                      ? primme->aNorm * primme->invBNorm
                      : primme->stats.estimateLargestSVal;
      }
   }
   else {
      if (!primme->massMatrixMatvec) {
         return max(primme->aNorm > 0.0 ? primme->aNorm : 0.0,
               primme->stats.estimateLargestSVal);
      } else {
         return max(primme->aNorm > 0.0 && primme->invBNorm > 0.0
                          ? primme->aNorm * primme->invBNorm
                          : 0.0,
               primme->stats.estimateLargestSVal);
      }
   }
}

/*******************************************************************************
 * Subroutine deltaEig - return an estimation of the minimum distance that
 * that two distinct eigenvalues can have. We estimate that as the smallest
 * e so that the eigenpair (\lambda+e,x) has residual vector norm smaller than 
 * (|Ax| + max(|\lambda|)|Bx|)*\epsilon := |(A,B)|*\epsilon.
 *
 * If (\lambda,x) is an exact eigenpair, then |e| is constrain as
 *
 *    |Ax - (\lambda+e)Bx| = |e|*|Bx| <= |(A,B)|*\epsilon
 *
 * If x is B-normal, then x'*B*x = x'*chol(B)*chol(B)'*x = 1 and
 * |chol(B)'*x| = 1. Therefore
 *
 *    |Bx| = |chol(B)*chol(B)'x| >= minsval(chol(B)) * |chol(B)'*x|
 *                               >= sqrt(minsval(B))
 *
 * Therefore, |e| <= sqrt(|inv(B)|) * |(A,B)| * \epsilon
 * 
 * INPUT PARAMETERS
 * ----------------
 * overrideUserEstimations    if nonzero, use estimations of |A| and |B| if
 *                            they are larger than primme.aNorm and primme.BNorm 
 * 
 * OUTPUT
 * ------
 * return                     estimation of the minimum distance
 ******************************************************************************/

TEMPLATE_PLEASE
HREAL deltaEig_Sprimme(
      int overrideUserEstimations, primme_context ctx) {

   primme_params *primme = ctx.primme;
   HREAL BNorm;

   if (overrideUserEstimations) {
      BNorm = max(primme->BNorm, primme->stats.estimateBNorm);
   }
   else {
      BNorm =
            (primme->BNorm > 0.0 ? primme->BNorm : primme->stats.estimateBNorm);
   }

   double eps_matrix;
   CHKERR(machineEpsMatrix_Sprimme(&eps_matrix, ctx));
   return problemNorm_Sprimme(overrideUserEstimations, primme) /
          sqrt(BNorm) * eps_matrix;
}

/*******************************************************************************
 * Function dist_dots - Computes several dot products in parallel
 *
 * Input Parameters
 * ----------------
 * x, y     Operands of the dot product operations
 *
 * ldx, ldy Leading dimension of x and y
 *
 * m        Length of the vectors x and y
 *
 * n        Number of columns in x and y
 *
 * ctx      Structure containing various solver parameters
 *
 * result   The inner products
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_dist_dots_Sprimme(SCALAR *x, PRIMME_INT ldx, SCALAR *y, PRIMME_INT ldy,
      PRIMME_INT m, int n, HSCALAR *result, primme_context ctx) {

   int i;
   for (i=0; i<n; i++) {
      result[i] = Num_dot_Sprimme(m, &x[ldx * i], 1, &y[ldy * i], 1, ctx);
   }
   CHKERR(globalSum_SHprimme(result, n, ctx));

   return 0;
}

/*******************************************************************************
 * Function dist_dots_real - Computes several dot products in parallel.
 *    Returns only the real part.
 *
 * Input Parameters
 * ----------------
 * x, y     Operands of the dot product operations
 *
 * ldx, ldy Leading dimension of x and y
 *
 * m        Length of the vectors x and y
 *
 * n        Number of columns in x and y
 *
 * ctx      Structure containing various solver parameters
 *
 * result   The inner products
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_dist_dots_real_Sprimme(SCALAR *x, PRIMME_INT ldx, SCALAR *y,
      PRIMME_INT ldy, PRIMME_INT m, int n, HREAL *result, primme_context ctx) {

   int i;
   for (i=0; i<n; i++) {
      result[i] =
            REAL_PART(Num_dot_Sprimme(m, &x[ldx * i], 1, &y[ldy * i], 1, ctx));
   }
   CHKERR(globalSum_RHprimme(result, n, ctx));

   return 0;
}

 
#endif /* SUPPORTED_TYPE */
