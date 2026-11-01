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
 * File: primme_svds_c.c
 *
 * Purpose - front end to svd problems. 
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../svds/primme_svds_c.c"
#endif

#include <string.h>  
#include "numerical.h"
#include "primme_interface.h"
#include "primme_svds_interface.h"
#include "../eigs/common_eigs.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "primme_svds_c.h"
#endif


#ifdef SUPPORTED_TYPE

#define UPDATE_STATS(PRIMME_SVDS_STATS, OP, PRIMME_STATS) {\
   (PRIMME_SVDS_STATS).numOuterIterations OP  (PRIMME_STATS).numOuterIterations;\
   (PRIMME_SVDS_STATS).numRestarts        OP  (PRIMME_STATS).numRestarts       ;\
   /* NOTE: for the augmented and for normal equations, every matvec for the  */\
   /* eigensolver involves the direct and the transpose matrix-vector product */\
   (PRIMME_SVDS_STATS).numMatvecs         OP  (PRIMME_STATS).numMatvecs*2      ;\
   (PRIMME_SVDS_STATS).numPreconds        OP  (PRIMME_STATS).numPreconds       ;\
   (PRIMME_SVDS_STATS).numGlobalSum       OP  (PRIMME_STATS).numGlobalSum      ;\
   (PRIMME_SVDS_STATS).numBroadcast       OP  (PRIMME_STATS).numBroadcast      ;\
   (PRIMME_SVDS_STATS).volumeGlobalSum    OP  (PRIMME_STATS).volumeGlobalSum   ;\
   (PRIMME_SVDS_STATS).volumeBroadcast    OP  (PRIMME_STATS).volumeBroadcast   ;\
   (PRIMME_SVDS_STATS).numOrthoInnerProds OP  (PRIMME_STATS).numOrthoInnerProds;\
   (PRIMME_SVDS_STATS).elapsedTime        OP  (PRIMME_STATS).elapsedTime       ;\
   (PRIMME_SVDS_STATS).timeMatvec         OP  (PRIMME_STATS).timeMatvec        ;\
   (PRIMME_SVDS_STATS).timePrecond        OP  (PRIMME_STATS).timePrecond       ;\
   (PRIMME_SVDS_STATS).timeOrtho          OP  (PRIMME_STATS).timeOrtho         ;\
   (PRIMME_SVDS_STATS).timeGlobalSum      OP  (PRIMME_STATS).timeGlobalSum     ;\
   (PRIMME_SVDS_STATS).timeBroadcast      OP  (PRIMME_STATS).timeBroadcast     ;\
   (PRIMME_SVDS_STATS).lockingIssue       OP  (PRIMME_STATS).lockingIssue      ;\
}

#endif /* SUPPORTED_TYPE */

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
   if (ctx.primme_svds && ctx.primme_svds->monitorFun) {
      int err = 0;
      primme_event event =
            (time >= -.5 ? primme_event_profile : primme_event_message);

#ifdef PRIMME_PROFILE
      /* Avoid profiling this function. It will turn out in a recursive call */
      ctx.path = NULL;
#endif

      CHKERRM((ctx.primme_svds->monitorFun(NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun, &time,
                     &event, NULL, ctx.primme_svds, &err),
                    err),
            PRIMME_USER_FAILURE, "Error code returned by 'monitorFun' %d", err);
   }
   return 0;
}


/******************************************************************************
 * Function primme_svds_get_context - return a context from the primme_svds_params
 *
 * PARAMETERS
 * ---------------------------
 * primme_svds      primme_svds_params struct
 *
 ******************************************************************************/

static primme_context primme_svds_get_context(primme_svds_params *primme_svds) {
   primme_context ctx;
   memset(&ctx, 0, sizeof(primme_context));
   if (primme_svds) {
      ctx.primme_svds = primme_svds;
      ctx.printLevel = primme_svds->printLevel;
      ctx.outputFile = primme_svds->outputFile;
      ctx.numProcs = primme_svds->numProcs;
      ctx.procID = primme_svds->procID;
      ctx.mpicomm = primme_svds->commInfo;
      ctx.queue = primme_svds->queue;
      ctx.report = monitor_report;
#ifdef PRIMME_PROFILE
      if (primme_svds->profile) {
         /* Compile regex. If there is no errors, set path to a nonzero       */
         /* value. Set ctx.report to the function that will channel the       */
         /* reports to the monitor. Report errors if they are.                */

         int ierr = regcomp(&ctx.profile, primme_svds->profile, REG_NOSUB);
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
 * Function primme_svds_free_context - free memory associated to the context
 *
 * PARAMETERS
 * ---------------------------
 * ctx   context
 *
 ******************************************************************************/

static void primme_svds_free_context(primme_context ctx) {

   /* Pop frame pushed in primme_get_context */

   Mem_pop_frame(&ctx);

   /* Free profiler */

#ifdef PRIMME_PROFILE
   if (ctx.path) regfree(&ctx.profile);
   if (ctx.timeoff) free(ctx.timeoff);
#endif
}

#endif /*USE_DOUBLE */

/*******************************************************************************
 * Subroutine Sprimme_svds - This routine is a front end used to perform 
 *    error checking on the input parameters, perform validation, 
 *    and make the call to main_iter. 
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * svals  Contains the converged Ritz values upon return.   Should be of size 
 *        primme_svds->numSvals.
 * 
 * svecs  The local portions of the converged Ritz vectors.  The dimension of
 *        the array is at least
 *        (primme_svds->mLocal+primme_svds->nLocal)*primme_svds->numSvals
 *
 * resNorms  The residual norms of the converged Ritz vectors.  Should be of 
 *           size primme_svds->numSvals
 *  
 * primme  Structure containing various solver parameters and statistics
 *         See readme.txt for INPUT/OUTPUT variables in primme
 *
 * Return Value
 * ------------
 *  0 - Success
 * -1 - Failure to allocate workspace
 * -2 - Malloc failure in allocating a permutation integer array
 * -3 - main_iter encountered a problem
 * -4 ...-32 - Invalid input (parameters or primme struct) returned 
 *             by check_input()
 * -100...-199 - PRIMME error code from first stage
 * -200...-299 - PRIMME error code from second stage
 *
 ******************************************************************************/

int Sprimme_svds(XREAL *svals, XSCALAR *svecs, XREAL *resNorms, 
      primme_svds_params *primme_svds) {

   return Xprimme_svds_aux(
         svals, svecs, resNorms, primme_svds, PRIMME_OP_SCALAR);
}

// Definition for *hsprimme_svds and *ksprimme_primme

#if defined(USE_HALF) || defined(USE_HALFCOMPLEX) ||                           \
      defined(USE_HALF_MAGMA) || defined(USE_HALFCOMPLEX_MAGMA) ||             \
      defined(USE_HALF_CUBLAS) || defined(USE_HALFCOMPLEX_CUBLAS)

// Expand the terms {,magma_}{hs,ks}primme
#  define Xsprimme_svds CONCAT(CONCAT(STEM, USE_ARITH(hs, ks)), primme_svds)

int Xsprimme_svds(float *svals, XSCALAR *svecs, float *resNorms,
      primme_svds_params *primme_svds) {

   return Xprimme_svds_aux(
         svals, svecs, resNorms, primme_svds, primme_op_float);
}

#  undef Xsprimme_svds
#endif

/*******************************************************************************
 * Subroutine Sprimme_svds - set defaults depending on the callee's type, and
 *    call wrapper_Sprimme with type set in internalPrecision. 
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * svals  Contains the converged Ritz values upon return.   Should be of size 
 *        primme_svds->numSvals.
 * 
 * svecs  The local portions of the converged Ritz vectors.  The dimension of
 *        the array is at least
 *        (primme_svds->mLocal+primme_svds->nLocal)*primme_svds->numSvals
 *
 * resNorms  The residual norms of the converged Ritz vectors.  Should be of 
 *           size primme_svds->numSvals
 *  
 * primme  Structure containing various solver parameters and statistics
 *         See readme.txt for INPUT/OUTPUT variables in primme
 *
 * svals_resNorms_type The type of the given arrays svals and resNorsm.
 *
 * Return Value
 * ------------
 * return error code
 ******************************************************************************/

STATIC int Xprimme_svds_aux(void *svals, XSCALAR *svecs, void *resNorms, 
      primme_svds_params *primme_svds, primme_op_datatype svals_resNorms_type) {

#ifdef SUPPORTED_TYPE

   /* Generate context */

   primme_context ctx = primme_svds_get_context(primme_svds);

   /* Set the current type as the default type for user's operators */

   if (primme_svds->matrixMatvec && primme_svds->matrixMatvec_type == primme_op_default)
      primme_svds->matrixMatvec_type = PRIMME_OP_SCALAR;
   if (primme_svds->applyPreconditioner && primme_svds->applyPreconditioner_type == primme_op_default)
      primme_svds->applyPreconditioner_type = PRIMME_OP_SCALAR;
   if (primme_svds->globalSumReal && primme_svds->globalSumReal_type == primme_op_default)
      primme_svds->globalSumReal_type = PRIMME_OP_SCALAR;
   if (primme_svds->broadcastReal && primme_svds->broadcastReal_type == primme_op_default)
      primme_svds->broadcastReal_type = PRIMME_OP_SCALAR;
   if (primme_svds->convTestFun && primme_svds->convTestFun_type == primme_op_default)
      primme_svds->convTestFun_type = PRIMME_OP_SCALAR;
   if (primme_svds->monitorFun && primme_svds->monitorFun_type == primme_op_default)
      primme_svds->monitorFun_type = PRIMME_OP_SCALAR;

   /* Number of returned singular triplets */

   int outInitSize = 0;

   /* call primme for the internal working precision */

   int ret;
   primme_op_datatype t = primme_svds->internalPrecision;
   if (t == primme_op_default) t = PRIMME_OP_SCALAR;
   switch (t) {
#  ifdef SUPPORTED_HALF_TYPE
   case primme_op_half:
      CHKERRVAL(wrapper_svds_Shprimme(svals, svecs, resNorms,
                      svals_resNorms_type, PRIMME_OP_SCALAR, &outInitSize, ctx),
            &ret);
      break;
#  endif
#  ifndef PRIMME_WITHOUT_FLOAT
   case primme_op_float:
      CHKERRVAL(wrapper_svds_Ssprimme(svals, svecs, resNorms,
                      svals_resNorms_type, PRIMME_OP_SCALAR, &outInitSize, ctx),
            &ret);
      break;
#  endif
   case primme_op_double:
      CHKERRVAL(wrapper_svds_Sdprimme(svals, svecs, resNorms,
                      svals_resNorms_type, PRIMME_OP_SCALAR, &outInitSize, ctx),
            &ret);
      break;
#  ifdef PRIMME_WITH_NATIVE_QUAD
   case primme_op_quad:
      CHKERRVAL(wrapper_svds_Sqprimme(svals, svecs, resNorms,
                      svals_resNorms_type, PRIMME_OP_SCALAR, &outInitSize, ctx),
            &ret);
      break;
#  endif
   default: ret = PRIMME_FUNCTION_UNAVAILABLE;
   }

   /* Free context */

   primme_svds_free_context(ctx);

   /* Set the number of returned triplets */

   primme_svds->initSize = outInitSize;

   return ret;
#else

   (void)svals;
   (void)svecs;
   (void)resNorms;
   (void)svals_resNorms_type;

   primme_svds->initSize = 0;
   return PRIMME_FUNCTION_UNAVAILABLE;
#endif /* SUPPORTED_TYPE */
}


#ifdef SUPPORTED_TYPE

/*******************************************************************************
 * Subroutine Sprimme_svds - Perform error checking on the input parameters,
 *    configure the eigensolvers, and make the solver calls. 
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * svals  Contains the converged Ritz values upon return.   Should be of size 
 *        primme_svds->numSvals.
 * 
 * svecs  The local portions of the converged Ritz vectors.  The dimension of
 *        the array is at least
 *        (primme_svds->mLocal+primme_svds->nLocal)*primme_svds->numSvals
 *
 * resNorms  The residual norms of the converged Ritz vectors.  Should be of 
 *           size primme_svds->numSvals
 *  
 * svals_resNorms_type The type of the arrays svals and resNorsm.
 *
 * svecs_type The type of the array svecs
 *
 * outInitSize The number of columns returned back.
 *
 * ctx    primme context
 *
 * Return Value
 * ------------
 * return error code
 ******************************************************************************/

TEMPLATE_PLEASE
int wrapper_svds_Sprimme(void *svals_, void *svecs_, void *resNorms_,
      primme_op_datatype svals_resNorms_type, primme_op_datatype svecs_type,
      int *outInitSize, primme_context ctx) {

   primme_svds_params *primme_svds = ctx.primme_svds;

   /* In case of error, return initSize = 0 */

   *outInitSize = 0;

   /* ------------------------------------ */
   /* Set defaults for sequential programs */
   /* ------------------------------------ */
   if (primme_svds->numProcs <= 1 && svals_ != NULL && svecs_ != NULL &&
         resNorms_ != NULL) {
      primme_svds->mLocal = primme_svds->m;
      primme_svds->nLocal = primme_svds->n;
      primme_svds->procID = 0;
      primme_svds->numProcs = 1;
   }

   /* ------------------ */
   /* Set some defaults  */
   /* ------------------ */
   primme_svds_set_defaults(primme_svds);

   /* Deprecated input  */

   if (svals_ == NULL && svecs_ == NULL && resNorms_ == NULL)
      return 0;

   /* ----------------------------------------------------------- */
   /* Primme_svds_initialize must be called by users unless users */  
   /* specify all parameters in primme_svds structure. Check if   */
   /* primme_svds inputs are good for bounds, correct values etc. */
   /* ----------------------------------------------------------- */
   CHKERR(primme_svds_check_input(svals_, svecs_, resNorms_, primme_svds)); 

   /* Set default convTetFun  */

   if (!primme_svds->convTestFun) {
      primme_svds->convTestFun = default_convTestFun;
      primme_svds->convTestFun_type = PRIMME_OP_SCALAR;
      if (primme_svds->eps == 0.0) {
         primme_svds->eps = MACHINE_EPSILON * 1e4;
      }
   }

   /* Set default monitor     */

   if (!primme_svds->monitorFun) {
      primme_svds->monitorFun = default_monitor_svds;
      primme_svds->monitorFun_type = PRIMME_OP_SCALAR;
   }

   /* ----------------------- */
   /* Reset stats             */
   /* ----------------------- */

   primme_svds->stats.numOuterIterations            = 0; 
   primme_svds->stats.numRestarts                   = 0;
   primme_svds->stats.numMatvecs                    = 0;
   primme_svds->stats.numPreconds                   = 0;
   primme_svds->stats.numGlobalSum                  = 0;
   primme_svds->stats.volumeGlobalSum               = 0;
   primme_svds->stats.numOrthoInnerProds            = 0.0;
   primme_svds->stats.elapsedTime                   = 0.0;
   primme_svds->stats.timeMatvec                    = 0.0;
   primme_svds->stats.timePrecond                   = 0.0;
   primme_svds->stats.timeOrtho                     = 0.0;
   primme_svds->stats.timeGlobalSum                 = 0.0;
   primme_svds->stats.lockingIssue                  = 0;

   /* Cast svals, svecs and resNorms to working precision */

   PRIMME_INT ldsvecs = primme_svds->mLocal + primme_svds->nLocal;
   HREAL *svals, *resNorms;
   SCALAR *svecs;
   CHKERR(Num_matrix_astype_RHprimme(svals_, 1, primme_svds->numSvals, 1,
         svals_resNorms_type, (void **)&svals, NULL, PRIMME_OP_HREAL,
         1 /* alloc */, 0 /* not copy */, ctx));
   CHKERR(Num_matrix_astype_Sprimme(svecs_,
         primme_svds->mLocal + primme_svds->nLocal,
         primme_svds->numOrthoConst +
               max(primme_svds->numSvals, primme_svds->initSize),
         ldsvecs, svecs_type, (void **)&svecs, NULL, PRIMME_OP_SCALAR,
         1 /* alloc */,
         primme_svds->numOrthoConst + primme_svds->initSize > 0 ? 1
                                                                : 0 /* copy? */,
         ctx));
   CHKERR(Num_matrix_astype_RHprimme(resNorms_, 1, primme_svds->numSvals,
         1, svals_resNorms_type, (void **)&resNorms, NULL, PRIMME_OP_HREAL,
         1 /* alloc */, 0 /* not copy */, ctx));

   /* --------------- */
   /* Execute stage 1 */
   /* --------------- */

   int ret, allocatedTargetShifts;

   SCALAR *svecs0;
   CHKERR(copy_last_params_from_svds(0, NULL, svecs,
            NULL, &allocatedTargetShifts, &svecs0, ctx));

   ret = Xprimme_aux_Sprimme(svals, (XSCALAR *)svecs0, resNorms,
         &primme_svds->primme, PRIMME_OP_HREAL);

   CHKERR(copy_last_params_to_svds(
            0, svals, svecs, resNorms, allocatedTargetShifts, ctx));

   if (ret != 0) ret = ret - 100;

   /* --------------- */
   /* Execute stage 2 */
   /* --------------- */

   if (primme_svds->methodStage2 != primme_svds_op_none && ret == 0) {
      CHKERR(copy_last_params_from_svds(
            1, svals, svecs, resNorms, &allocatedTargetShifts, &svecs0, ctx));

      /* The value numSvals-primme->numEvals indicates how many svals */
      /* are already converged. So shift svals and resnorms that much */
      int nconv = primme_svds->numSvals - primme_svds->primmeStage2.numEvals;

      ret = Xprimme_aux_Sprimme(svals + nconv, (XSCALAR *)svecs0, resNorms + nconv,
            &primme_svds->primmeStage2, PRIMME_OP_HREAL);

      CHKERR(copy_last_params_to_svds(
            1, svals, svecs, resNorms, allocatedTargetShifts, ctx));

      if (ret != 0)  ret = ret - 200;
   }

   /* Copy back evals, evecs and resNorms */

   CHKERR(Num_matrix_astype_RHprimme(svals,
         1, primme_svds->initSize, 1, PRIMME_OP_HREAL, (void **)&svals_, NULL,
         svals_resNorms_type, -1 /* destroy */, 1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_Sprimme(svecs,
         primme_svds->mLocal + primme_svds->nLocal,
         primme_svds->numOrthoConst + primme_svds->initSize, ldsvecs,
         PRIMME_OP_SCALAR, (void **)&svecs_, &ldsvecs, svecs_type,
         -1 /* destroy */, 1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_RHprimme(resNorms, 1, primme_svds->initSize, 1,
         PRIMME_OP_HREAL, (void **)&resNorms_, NULL, svals_resNorms_type,
         -1 /* destroy */, 1 /* copy */, ctx));

   /* If no error, return initSize */

   *outInitSize = primme_svds->initSize;

   return ret;
}

STATIC int comp_double(const void *a, const void *b)
{
   return *(double*)a <= *(double*)b ? -1 : 1;
}

STATIC int copy_last_params_from_svds(int stage, HREAL *svals, SCALAR *svecs,
      HREAL *rnorms, int *allocatedTargetShifts,
      SCALAR **out_svecs, primme_context ctx) {

   primme_svds_params *primme_svds = ctx.primme_svds;
   primme_params *primme;
   primme_svds_operator method;
   SCALAR *aux;
   *out_svecs = svecs;
   int n, nMax, i;

   primme = stage == 0 ? &primme_svds->primme : &primme_svds->primmeStage2;
   method = stage == 0 ? primme_svds->method : primme_svds->methodStage2;

   *allocatedTargetShifts = 0;

   if (method == primme_svds_op_none) {
      primme->maxMatvecs = 0;
      return 0;
   }

   if (!primme->matrixMatvec) {
      switch (primme_svds->matrixMatvec_type) {
#ifdef SUPPORTED_HALF_TYPE
      case primme_op_half:   primme->matrixMatvec = matrixMatvec_eigs_Shprimme; break;
#endif
#ifndef PRIMME_WITHOUT_FLOAT
      case primme_op_float:  primme->matrixMatvec = matrixMatvec_eigs_Ssprimme; break;
#endif
      case primme_op_double: primme->matrixMatvec = matrixMatvec_eigs_Sdprimme; break;
#ifdef PRIMME_WITH_NATIVE_QUAD
      case primme_op_quad:   primme->matrixMatvec = matrixMatvec_eigs_Sqprimme; break;
#endif
      default: return PRIMME_FUNCTION_UNAVAILABLE;
      }
      primme->matrixMatvec_type = primme_svds->matrixMatvec_type;
      primme->matrix = primme_svds;
   }
   if (primme_svds->applyPreconditioner && !primme->applyPreconditioner) {
      switch (primme_svds->applyPreconditioner_type) {
#ifdef SUPPORTED_HALF_TYPE
      case primme_op_half:   primme->applyPreconditioner = applyPreconditioner_eigs_Shprimme; break;
#endif
#ifndef PRIMME_WITHOUT_FLOAT
      case primme_op_float:  primme->applyPreconditioner = applyPreconditioner_eigs_Ssprimme; break;
#endif
      case primme_op_double: primme->applyPreconditioner = applyPreconditioner_eigs_Sdprimme; break;
#ifdef PRIMME_WITH_NATIVE_QUAD
      case primme_op_quad:   primme->applyPreconditioner = applyPreconditioner_eigs_Sqprimme; break;
#endif
      default: return PRIMME_FUNCTION_UNAVAILABLE;
      }
      primme->applyPreconditioner_type = primme_svds->applyPreconditioner_type;
      primme->preconditioner = primme_svds;
   }

   if (primme_svds->aNorm > 0.0) {
      switch (method) {
         case primme_svds_op_AtA:
         case primme_svds_op_AAt:
            primme->aNorm = primme_svds->aNorm * primme_svds->aNorm;
            break;
         case primme_svds_op_augmented:
            primme->aNorm = primme_svds->aNorm;
            break;
         case primme_svds_op_none:
            break;
      }
   }

   switch (method) {
      case primme_svds_op_AtA:
      case primme_svds_op_AAt:
         primme->convTestFun = convTestFunATA;
         break;
      case primme_svds_op_augmented:
         primme->convTestFun = convTestFunAug;
         break;
      case primme_svds_op_none:
         break;
   }
   primme->convTestFun_type = PRIMME_OP_SCALAR;

   /* Set properly initial vectors. Now svecs = [Uc U0 Vc V0], where          */
   /* Uc, m x numOrthoConst, left constrain vectors;                          */
   /* U0, m x initSize, left initial vectors;                                 */
   /* Vc, n x numOrthoConst, right constrain vectors;                         */
   /* V0, n x numOrthoConst, right initial vectors.                           */

   primme->initSize = primme_svds->initSize;
   primme->numOrthoConst = primme_svds->numOrthoConst;
   n = primme_svds->initSize + primme_svds->numOrthoConst;
   nMax = max(primme_svds->initSize, primme_svds->numSvals) +
      primme_svds->numOrthoConst;
   switch (method) {
      case primme_svds_op_AtA:
      case primme_svds_op_AAt:
         /* Move Vc V0 to the rightmost position in svecs (aux).
            If using AA', only move Vc */
         aux = &svecs[nMax * primme_svds->mLocal];
         Num_copy_matrix_Sprimme(
               &svecs[primme_svds->mLocal * n], primme_svds->nLocal,
               method == primme_svds_op_AtA ? n : primme_svds->numOrthoConst,
               primme_svds->nLocal, aux, primme_svds->nLocal, ctx);
         if (method == primme_svds_op_AtA)
            *out_svecs = aux;
         break;
      case primme_svds_op_augmented:
         /* Shuffle svecs so that svecs = [V; U] */
         assert(primme->nLocal == primme_svds->mLocal + primme_svds->nLocal);
         CHKERR(Num_malloc_Sprimme(primme->nLocal * n, &aux, ctx));
         Num_copy_Sprimme(primme->nLocal * n, svecs, 1, aux, 1, ctx);
         Num_copy_matrix_Sprimme(&aux[primme_svds->mLocal * n], primme_svds->nLocal,
               n, primme_svds->nLocal, svecs, primme->nLocal, ctx);
         Num_copy_matrix_Sprimme(aux, primme_svds->mLocal, n, primme_svds->mLocal,
               &svecs[primme_svds->nLocal], primme->nLocal, ctx);
         CHKERR(Num_free_Sprimme(aux, ctx));

         /* Normalize the orthogonal constrains */
         Num_scal_Sprimme(primme->nLocal * primme_svds->numOrthoConst, 1. / sqrt(2.),
               svecs, 1, ctx);
         break;
      case primme_svds_op_none:
         break;
   }
   primme->iseed[0] = primme_svds->iseed[0];
   primme->iseed[1] = primme_svds->iseed[1];
   primme->iseed[2] = primme_svds->iseed[2];
   primme->iseed[3] = primme_svds->iseed[3];
   if (stage == 0) {
      primme->maxMatvecs = primme_svds->maxMatvecs / 2;
   } else {
      primme->maxMatvecs =
         primme_svds->maxMatvecs / 2 - primme_svds->primme.stats.numMatvecs;
   }

   if ((stage == 0 && primme_svds->numTargetShifts > 0) ||
         (stage == 1 && primme->targetShifts == NULL &&
          primme_svds->target == primme_svds_closest_abs)) {
      primme->numTargetShifts = primme_svds->numTargetShifts;
      if (stage == 0 &&
            (method == primme_svds_op_AtA || method == primme_svds_op_AAt)) {
         *allocatedTargetShifts = 1;
         CHKERR(Num_malloc_dprimme(primme_svds->numSvals, &primme->targetShifts,
                  ctx));
         /* Previous allocation is going to be freed at                       */
         /* copy_last_params_from_svds. To avoid complaining by the memory    */
         /* manager, ask to keep the frame.                                   */
         Mem_keep_frame(ctx);
         for (i = 0; i < primme->numTargetShifts; i++) {
            primme->targetShifts[i] =
               primme_svds->targetShifts[i] * primme_svds->targetShifts[i];
         }
      } else {
         primme->targetShifts = primme_svds->targetShifts;
      }
   } else if (stage == 1 && primme->targetShifts == NULL &&
         primme_svds->target == primme_svds_smallest) {

      assert(method == primme_svds_op_augmented);
      *allocatedTargetShifts = 1;
      CHKERR(Num_malloc_dprimme(
               primme_svds->numSvals, &primme->targetShifts, ctx));

      /* Previous allocation is going to be freed at                       */
      /* copy_last_params_from_svds. To avoid complaining by the memory    */
      /* manager, ask to keep the frame.                                   */
      Mem_keep_frame(ctx);

      /* primme was configured to find the closest but greater values than */
      /* some shift. The eigensolver is not able to distinguish eigenvalues*/
      /* separated by less than machEps*|A|. The augmented matrix has      */
      /* |m-n| eigenpairs with value zero that don't correspond to         */
      /* singular triplets of A. To avoid to return incorrect triplets set */
      /* shifts not smaller than machEps*|A|.                              */
      /* If d^2 and d'^2 are the exact and the approximate eigenvalues     */
      /* from normal equations respectively, and assuming that d' >= d,    */
      /* then d can be lower bounded as:                                   */
      /*    d'^2 - d^2 <= |r_AtA| -> sqrt(d'^2-|r_AtA|) <= d               */
      /*                             sqrt(d'^2-|r|*d')  <= d               */

      double min_val = primme_svds->aNorm * MACHINE_EPSILON;
      for (i = 0; i < primme_svds->initSize; i++) {
         primme->targetShifts[i] =
            max(sqrt(fabs(max(svals[i] - rnorms[i], 0.0) * svals[i])), min_val);
      }
      for (; i < primme_svds->numSvals; i++) {
         primme->targetShifts[i] = min_val;
      }

      /* Sort the shifts in ascending order */

      qsort(primme->targetShifts, primme_svds->numSvals, sizeof(double),
            comp_double);
      primme->numTargetShifts = primme_svds->numSvals;

   } else if (method == primme_svds_op_augmented &&
         primme_svds->target == primme_svds_smallest &&
         primme->targetShifts == NULL) {

      CHKERR(Num_malloc_dprimme(1, &primme->targetShifts, ctx));

      /* Previous allocation is going to be freed at                       */
      /* copy_last_params_from_svds. To avoid complaining by the memory    */
      /* manager, ask to keep the frame.                                   */
      Mem_keep_frame(ctx);

      *allocatedTargetShifts = 1;
      primme->targetShifts[0] = 0.0;
      primme->numTargetShifts = 1;
   }

   /* Set an initial guess [x; A'x] or [Ax; x] if there is no initial guess   */
   /* and augmented matrix will be used                                       */

   if (method == primme_svds_op_augmented && primme->initSize <= 0) {
      HREAL norms2_[2], norms2[2];
      SCALAR *svecs0 = &svecs[primme->numOrthoConst*primme->nLocal];
      if (primme_svds->m >= primme_svds->n) {
         Num_larnv_Sprimme(2, primme->iseed, primme_svds->mLocal,
               &svecs0[primme_svds->nLocal], ctx);
         CHKERR(matrixMatvecSVDS_Sprimme(
                     &svecs0[primme_svds->nLocal], primme_svds->mLocal, svecs0,
                     primme_svds->nLocal, 0, 1, 1 /* trans */, ctx));
      } else {
         Num_larnv_Sprimme(2, primme->iseed, primme_svds->nLocal, svecs0, ctx);
         CHKERR(matrixMatvecSVDS_Sprimme(
                     svecs0, primme_svds->nLocal, &svecs0[primme_svds->nLocal],
                     primme_svds->mLocal, 0, 1, 0 /* no trans */, ctx));
      }
      norms2_[0] = REAL_PART(
            Num_dot_Sprimme(primme_svds->nLocal, svecs0, 1, svecs0, 1, ctx));
      norms2_[1] = REAL_PART(
            Num_dot_Sprimme(primme_svds->mLocal, &svecs0[primme_svds->nLocal], 1,
               &svecs0[primme_svds->nLocal], 1, ctx));
      CHKERR(globalSum_Rprimme_svds(norms2_, norms2, 2, ctx));
      Num_scal_Sprimme(primme_svds->nLocal, 1.0 / sqrt(norms2[0]), svecs0, 1, ctx);
      Num_scal_Sprimme(primme_svds->mLocal, 1.0 / sqrt(norms2[1]),
            &svecs0[primme_svds->nLocal], 1, ctx);
      primme->initSize = 1;
      if (rnorms)
         rnorms[0] = HOST_MACHINE_MAX;
      primme->initBasisMode = primme_init_user;
   }

   /* If second stage, set as numOrthoConst the first ones that pass */
   /* the convergence criterion.                                     */

   if (stage == 1) {
      assert(method == primme_svds_op_augmented);
      int *flags;
      CHKERR(Num_malloc_iprimme(primme->initSize, &flags, ctx));

      for (i = 0; primme->initSize > 0; i++) {
         /* NOTE: convTestFun at this stage expects the residual norm for the */
         /*       the augmented problem; this is why the residual norm is     */
         /*       divided by sqrt(2).                                         */
         int isConv = 0;
         CHKERR(convTestFunSVDS_Sprimme(svals[i],
               &svecs[primme->nLocal * primme->numOrthoConst +
                      primme_svds->nLocal],
               1 /* vector given */,
               &svecs[primme->nLocal * primme->numOrthoConst],
               1 /* vector given */, rnorms[i], method, &isConv, ctx));
         if (!isConv) break;

         /* Report a triplet is locked */
         int numLocked = i + 1;
         flags[i] = CONVERGED;
         CHKERR(monitorFunSVDS_Sprimme(NULL, 0, NULL, NULL, 0, NULL, 0, svals,
               numLocked, flags, rnorms, 0, 0.0, NULL, 0.0, primme_event_locked,
               0 /* stage 0 */, -1.0 /* don't update elapsedTime */, ctx));

         primme->numOrthoConst++;
         primme->initSize--;
         primme->numEvals--;
      }

      CHKERR(Num_free_iprimme(flags, ctx));
   }

   /* Set locking */

   if (primme_svds->locking >= 0) {
      primme->locking = primme_svds->locking;
   }

   /* Set monitor */

   if (primme->monitorFun == NULL) {
      if (primme_svds->methodStage2 == primme_svds_op_none) {
         primme->monitorFun = monitor_single_stage;
      } else if (stage == 0) {
         primme->monitorFun = monitor_stage1;
      } else {
         primme->monitorFun = monitor_stage2;
      }
      primme->monitorFun_type = PRIMME_OP_SCALAR;
   }

   /* Copy queue */
   primme->queue = primme_svds->queue;

   /* Copy profile */
   primme->profile = primme_svds->profile;

   return 0;
}

STATIC int copy_last_params_to_svds(int stage, HREAL *svals, SCALAR *svecs,
      HREAL *rnorms, int allocatedTargetShifts,
      primme_context ctx) {

   primme_svds_params *primme_svds = ctx.primme_svds;
   primme_params *primme;
   primme_svds_operator method;
   int n, nMax, i;

   primme = stage == 0 ? &primme_svds->primme : &primme_svds->primmeStage2;
   method = stage == 0 ? primme_svds->method : primme_svds->methodStage2;

   if (method == primme_svds_op_none) {
      primme->maxMatvecs = 0;
      return 0;
   }

   /* Pass back the converged vectors in first stage to regular vectors */

   if (stage == 1) {
      int nconv = primme_svds->numSvals - primme->numEvals;
      primme->initSize += nconv;
      primme->numOrthoConst -= nconv;
      primme->numEvals += nconv;
   }

   /* Record performance measurements */
   UPDATE_STATS(primme_svds->stats, +=, primme->stats);

   if (primme->aNorm > 0.0) {
      switch (method) {
         case primme_svds_op_AtA:
         case primme_svds_op_AAt:
            primme_svds->aNorm = sqrt(primme->aNorm);
            break;
         case primme_svds_op_augmented:
            primme_svds->aNorm = primme->aNorm;
            break;
         case primme_svds_op_none:
            break;
      }
   }

   if (method == primme_svds_op_AtA || method == primme_svds_op_AAt) {
      for (i = 0; i < primme->initSize; i++) {
         svals[i] = sqrt(max(0.0, svals[i]));
      }
   }

   /* Set svecs = [Uc U Vc V] */
   nMax = max(primme_svds->initSize, primme_svds->numSvals) +
      primme_svds->numOrthoConst;
   primme_svds->initSize = primme->initSize;
   n = primme_svds->initSize + primme_svds->numOrthoConst;
   switch (method) {
      case primme_svds_op_AtA:
         /* Transform svecs to [Uc A*V/Sigma Vc V] */
         CHKERR(matrixMatvecSVDS_Sprimme(
               &svecs[primme_svds->mLocal * nMax +
                      primme->nLocal * primme_svds->numOrthoConst],
               primme_svds->nLocal,
               &svecs[primme_svds->mLocal * primme_svds->numOrthoConst],
               primme_svds->mLocal, 0, primme_svds->initSize, 0 /* no trans */,
               ctx));
         CHKERR(Num_scalInv_Smatrix(
               &svecs[primme_svds->mLocal * primme_svds->numOrthoConst],
               primme_svds->mLocal, primme_svds->initSize, primme_svds->mLocal,
               svals, ctx));
         Num_copy_matrix_Sprimme(&svecs[primme_svds->mLocal * nMax],
               primme_svds->nLocal, n, primme_svds->nLocal,
               &svecs[primme_svds->mLocal * n],
               primme_svds->nLocal, ctx);
         break;
      case primme_svds_op_AAt:
         /* Transform svecs to [Uc U Vc A'*U/Sigma] */
         Num_copy_matrix_Sprimme(
               &svecs[primme_svds->mLocal * nMax], primme_svds->nLocal,
               primme_svds->numOrthoConst, primme_svds->nLocal,
               &svecs[primme_svds->mLocal * n], primme_svds->nLocal, ctx);
         CHKERR(matrixMatvecSVDS_Sprimme(
               &svecs[primme_svds->mLocal * primme_svds->numOrthoConst],
               primme_svds->mLocal,
               &svecs[primme_svds->mLocal * n +
                      primme->nLocal * primme_svds->numOrthoConst],
               primme_svds->nLocal, 0, primme_svds->initSize, 1 /* trans */,
               ctx));
         CHKERR(Num_scalInv_Smatrix(
               &svecs[primme_svds->mLocal * n +
                      primme->nLocal * primme_svds->numOrthoConst],
               primme_svds->nLocal, primme_svds->initSize, primme_svds->nLocal,
               svals, ctx));
         break;
      case primme_svds_op_augmented:
         assert(primme->nLocal == primme_svds->mLocal + primme_svds->nLocal);

         /* Normalize back the orthogonal constrains */
         Num_scal_Sprimme(primme->nLocal * primme_svds->numOrthoConst, sqrt(2.),
               svecs, 1, ctx);

         /* Shuffle svecs from [Vc V; Uc U] to [Uc U Vc V] */
         SCALAR *aux;
         CHKERR(Num_malloc_Sprimme(primme->nLocal * n, &aux, ctx));
         Num_copy_Sprimme(primme->nLocal * n, svecs, 1, aux, 1, ctx);
         Num_copy_matrix_Sprimme(aux, primme_svds->nLocal, n, primme->nLocal,
               &svecs[primme_svds->mLocal * n],
               primme_svds->nLocal, ctx);
         Num_copy_matrix_Sprimme(&aux[primme_svds->nLocal], primme_svds->mLocal, n,
               primme->nLocal, svecs, primme_svds->mLocal, ctx);
         CHKERR(Num_free_Sprimme(aux, ctx));

         /* Normalize every column in U and V */
         HREAL *norms2;
         CHKERR(Num_malloc_RHprimme(2 * n, &norms2, ctx));
         for (i = 0; i < n; i++) {
            norms2[i] = REAL_PART(Num_dot_Sprimme(
                     primme_svds->mLocal, &svecs[primme_svds->mLocal * i], 1,
                     &svecs[primme_svds->mLocal * i], 1, ctx));
         }
         for (i = 0; i < n; i++) {
            norms2[n + i] = REAL_PART(Num_dot_Sprimme(
                     primme_svds->nLocal,
                     &svecs[primme_svds->mLocal * n + primme_svds->nLocal * i], 1,
                     &svecs[primme_svds->mLocal * n + primme_svds->nLocal * i], 1, ctx));
         }
         CHKERR(globalSum_Rprimme_svds(norms2, norms2, 2 * n, ctx));
         for (i = 0; i < n; i++) {
            Num_scal_Sprimme(primme_svds->mLocal, 1.0 / sqrt(norms2[i]),
                  &svecs[primme_svds->mLocal * i], 1, ctx);
         }
         for (i = 0; i < n; i++) {
            Num_scal_Sprimme(
                  primme_svds->nLocal, 1.0 / sqrt(norms2[n + i]),
                  &svecs[primme_svds->mLocal * n + primme_svds->nLocal * i], 1, ctx);
         }
         CHKERR(Num_free_RHprimme(norms2, ctx));
         break;
      case primme_svds_op_none:
         break;
   }

   primme_svds->iseed[0] = primme->iseed[0];
   primme_svds->iseed[1] = primme->iseed[1];
   primme_svds->iseed[2] = primme->iseed[2];
   primme_svds->iseed[3] = primme->iseed[3];

   if (allocatedTargetShifts) {
      CHKERR(Num_free_dprimme(primme->targetShifts, ctx));
      primme->targetShifts = NULL;
   }

   /* Update residual norms. For normal equations we have to divide by the    */
   /* the singular value. For the augmented, nothing is required because the  */
   /* user defined stopping criterion computes the actual residual vector     */
   /* norm and update the value.                                              */

   switch (method) {
      case primme_svds_op_AtA:
      case primme_svds_op_AAt:
         for (i = 0; i < primme_svds->initSize; i++) {
            rnorms[i] = min(rnorms[i] / svals[i], primme_svds->aNorm);
         }
         break;
      case primme_svds_op_augmented:
         for (i = 0; i < primme_svds->initSize; i++) {
            rnorms[i] *= sqrt(2.0);
         }
         break;
      case primme_svds_op_none:
         break;
   }


   return 0;
}

/******************************************************************************
 *
 * int primme_svds_check_input(double *svals, SCALAR *svecs, double *resNorms, 
 *                        primme_svds_params *primme_svds) 
 *
 * INPUT
 * -----
 *  svals, svecs, resNorms   Output arrays for primme
 *  primme_svds              the main structure of parameters 
 *
 * return value -   0    If input parameters in primme are appropriate
 *              -4..-19  Inappropriate input parameters were found
 *
 ******************************************************************************/
STATIC int primme_svds_check_input(void *svals, void *svecs, void *resNorms, 
      primme_svds_params *primme_svds) {
   int ret;
   ret = 0;

   if (primme_svds == NULL)
      ret = -4;
   else if (primme_svds->n < 0 || primme_svds->m < 0 || primme_svds->nLocal < 0
         || primme_svds->mLocal < 0 || primme_svds->nLocal > primme_svds->n
         || primme_svds->mLocal > primme_svds->m) 
      ret = -5;
   else if (primme_svds->numProcs < 1)
      ret = -6;
   else if (primme_svds->matrixMatvec == NULL) 
      ret = -7;
   else if (primme_svds->applyPreconditioner == NULL && 
         primme_svds->precondition == 1) 
      ret = -8;
   else if (primme_svds->numProcs >1 && primme_svds->globalSumReal == NULL)
      ret = -9;
   else if (primme_svds->numSvals > min(primme_svds->n, primme_svds->m))
      ret = -10;
   else if (primme_svds->numSvals < 1)
      ret = -11;
   else if ( primme_svds->target != primme_svds_smallest  &&
         primme_svds->target != primme_svds_largest   &&
         primme_svds->target != primme_svds_closest_abs)
      ret = -13;
   else if ( primme_svds->method != primme_svds_op_AtA &&
         primme_svds->method != primme_svds_op_AAt &&
         primme_svds->method != primme_svds_op_augmented)
      ret = -14;
   else if ( (primme_svds->method == primme_svds_op_augmented &&
            primme_svds->methodStage2 != primme_svds_op_none) ||
         (primme_svds->method != primme_svds_op_augmented &&
          primme_svds->methodStage2 != primme_svds_op_augmented &&
          primme_svds->methodStage2 != primme_svds_op_none))
      ret = -15;
   else if (primme_svds->printLevel < 0 || primme_svds->printLevel > 5)
      ret = -16; 
   else if (svals == NULL)
      ret = -17;
   else if (svecs == NULL)
      ret = -18;
   else if (resNorms == NULL)
      ret = -19;

   return ret;
   /***************************************************************************/
} /* end of check_input
   ***************************************************************************/

/*******************************************************************************
 * Subroutine matrixMatvec_ - Computes A*V(:,nv+1:,nv+blksze), or A'*...,
 *    where A is the operator defined by primme_svds->matrixMatvec
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V          The orthonormal basis
 * ldV        The leading dimension of V
 * ldW        The leading dimension of W
 * basisSize  Number of vectors in V
 * blockSize  The current block size
 * conj       A*V for conj==0, and A'*V otherwise
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * W          A*V or A'*V
 ******************************************************************************/

TEMPLATE_PLEASE
int matrixMatvecSVDS_Sprimme(SCALAR *V, PRIMME_INT ldV, SCALAR *W, PRIMME_INT ldW,
      int basisSize, int blockSize, int conj, primme_context ctx) {

   primme_svds_params *primme_svds = ctx.primme_svds;

   if (blockSize <= 0)
      return 0;

   PRIMME_INT mV, mW;
   if (conj) {
      mV = primme_svds->mLocal;
      mW = primme_svds->nLocal;
   } else {
      mV = primme_svds->nLocal;
      mW = primme_svds->mLocal;
   }
   assert(ldV >= mV && ldW >= mW);

   double t0 = primme_wTimer();

   /* Cast V and W */

   SCALAR *Vb = &V[ldV * basisSize], *Wb = &W[ldW * basisSize];
   void *V0, *W0;
   PRIMME_INT ldV0, ldW0;
   CHKERR(Num_matrix_astype_Sprimme(Vb, mV, blockSize, ldV, PRIMME_OP_SCALAR,
         &V0, &ldV0, primme_svds->matrixMatvec_type, 1 /* alloc */,
         1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_Sprimme(Wb, mW, blockSize, ldW, PRIMME_OP_SCALAR,
         &W0, &ldW0, primme_svds->matrixMatvec_type, 1 /* alloc */,
         0 /* no copy */, ctx));

   /* W(:,c) = A*V(:,c) for c = basisSize:basisSize+blockSize-1 */

   int ierr = 0;
   CHKERRM((primme_svds->matrixMatvec(
                  V0, &ldV0, W0, &ldW0, &blockSize, &conj, primme_svds, &ierr),
                 ierr),
         PRIMME_USER_FAILURE, "Error returned by 'matrixMatvec' %d", ierr);

   /* Copy back W */

   CHKERR(Num_matrix_astype_Sprimme(W0, mW, blockSize, ldW0,
         primme_svds->matrixMatvec_type, (void **)&Wb, &ldW,
         PRIMME_OP_SCALAR, 0 /* not alloc */, 1 /* copy */, ctx));

   if (Vb != V0) CHKERR(Num_free_Sprimme((SCALAR*)V0, ctx));
   if (Wb != W0) CHKERR(Num_free_Sprimme((SCALAR*)W0, ctx));

   primme_svds->stats.timeMatvec += primme_wTimer() - t0;
   primme_svds->stats.numMatvecs += blockSize;

   return 0;
}

/*******************************************************************************
 * Subroutine convTestFun - wrapper around primme.convTestFun; evaluate if
 *    the approximate eigenpair eval, evec with given residual norm is
 *    considered as converged.
 *
 * INPUT PARAMETERS
 * ----------------
 * sval       the eigenvalue
 * svec       the eigenvector
 * givenSvec  whether eigenvector is provided
 * rNorm      the residual vector norm
 * 
 * OUTPUT
 * ------
 * isconv   if non-zero, the triplet is considered converged.
 ******************************************************************************/

TEMPLATE_PLEASE
int convTestFunSVDS_Sprimme(HREAL sval, SCALAR *leftsvec, int givenLeftSvec,
      SCALAR *rightsvec, int givenRightSvec, HREAL rNorm, int method,
      int *isconv, primme_context ctx) {

   primme_svds_params *primme_svds = ctx.primme_svds;

   /* Cast sval and rNorm */

   double svald = sval, rNormd = rNorm;

   /* Cast svecs if given */

   void *leftsvec0 = NULL, *rightsvec0 = NULL;
   if (leftsvec && givenLeftSvec) {
      CHKERR(Num_matrix_astype_Sprimme(leftsvec, primme_svds->mLocal, 1,
            primme_svds->mLocal, PRIMME_OP_SCALAR, &leftsvec0, NULL,
            primme_svds->convTestFun_type, 1 /* alloc */, 1 /* copy */, ctx));
   }
   if (rightsvec && givenRightSvec) {
      CHKERR(Num_matrix_astype_Sprimme(rightsvec, primme_svds->nLocal, 1,
            primme_svds->nLocal, PRIMME_OP_SCALAR, &rightsvec0, NULL,
            primme_svds->convTestFun_type, 1 /* alloc */, 1 /* copy */, ctx));
   }

   // If a svec is going to be passed to convTestFun, but mLocal or nLocal is 0
   // then fake the svec with a nonzero pointer in order to not be mistaken   
   // by not passing a vector.                                                
   SCALAR dummy;

   if (primme_svds->mLocal == 0 && givenLeftSvec) leftsvec0 = &dummy;
   if (primme_svds->nLocal == 0 && givenRightSvec) rightsvec0 = &dummy;

   int ierr=0;
   CHKERRM((primme_svds->convTestFun(&svald, givenLeftSvec ? leftsvec : NULL,
                  givenRightSvec ? rightsvec : NULL, &rNormd, &method, isconv,
                  primme_svds, &ierr),
                 ierr),
         -1, "Error returned by 'convTestFun' %d", ierr);

   if (primme_svds->mLocal > 0 && leftsvec && givenLeftSvec &&
         leftsvec != (SCALAR *)leftsvec0) {
      CHKERR(Num_free_Sprimme((SCALAR*)leftsvec0, ctx));
   }
   if (primme_svds->nLocal > 0 && rightsvec && givenRightSvec &&
         rightsvec != (SCALAR *)rightsvec0) {
      CHKERR(Num_free_Sprimme((SCALAR*)rightsvec0, ctx));
   }

   return 0;
}

TEMPLATE_PLEASE
int monitorFunSVDS_Sprimme(HREAL *basisSvals, int basisSize, int *basisFlags,
      int *iblock, int blockSize, HREAL *basisNorms, int numConverged,
      HREAL *lockedSvals, int numLocked, int *lockedFlags, HREAL *lockedNorms,
      int inner_its, HREAL LSRes, const char *msg, double time,
      primme_event event, int stage, double startTime, primme_context ctx) {

   /* Quick exit */

   primme_svds_params *primme_svds = ctx.primme_svds;
   if (!primme_svds->monitorFun) return 0;

   /* Cast basisSvals, basisNorms, lockedSvals, lockedNorms and LSRes */
   /* to monitorFun_type                                              */

   void *basisSvals0, *basisNorms0, *lockedSvals0, *lockedNorms0, *LSRes0;
   CHKERR(Num_matrix_astype_RHprimme(basisSvals, 1, basisSize, 1,
         PRIMME_OP_HREAL, (void **)&basisSvals0, NULL,
         primme_svds->monitorFun_type, 1 /* alloc */, 1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_RHprimme(basisNorms, 1, basisSize, 1,
         PRIMME_OP_HREAL, (void **)&basisNorms0, NULL,
         primme_svds->monitorFun_type, 1 /* alloc */, 1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_RHprimme(lockedSvals, 1, numLocked, 1,
         PRIMME_OP_HREAL, (void **)&lockedSvals0, NULL,
         primme_svds->monitorFun_type, 1 /* alloc */, 1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_RHprimme(lockedNorms, 1, numLocked, 1,
         PRIMME_OP_HREAL, (void **)&lockedNorms0, NULL,
         primme_svds->monitorFun_type, 1 /* alloc */, 1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_RHprimme(&LSRes, 1, 1, 1, PRIMME_OP_HREAL,
         (void **)&LSRes0, NULL, primme_svds->monitorFun_type, 1 /* alloc */,
         1 /* copy */, ctx));

   /* Call the user-defined functions */

   if (startTime > 0.0) {
      primme_svds->stats.elapsedTime = primme_wTimer() - startTime;
   }

   int err = 0;
   CHKERRM(
         (primme_svds->monitorFun(basisSvals0, &basisSize, basisFlags, iblock,
                &blockSize, basisNorms0, &numConverged, lockedSvals0,
                &numLocked, lockedFlags, lockedNorms0,
                inner_its >= 0 ? &inner_its : NULL, LSRes >= 0 ? LSRes0 : NULL,
                msg, &time, &event, &stage, primme_svds, &err),
               err),
         -1, "Error returned by monitorFun: %d", err);

   if (basisSvals != (HREAL *)basisSvals0)
      CHKERR(Num_free_RHprimme((HREAL *)basisSvals0, ctx));
   if (basisNorms != (HREAL *)basisNorms0)
      CHKERR(Num_free_RHprimme((HREAL *)basisNorms0, ctx));
   if (lockedSvals != (HREAL *)lockedSvals0)
      CHKERR(Num_free_RHprimme((HREAL *)lockedSvals0, ctx));
   if (lockedNorms != (HREAL *)lockedNorms0)
      CHKERR(Num_free_RHprimme((HREAL *)lockedNorms0, ctx));
   if (&LSRes != (HREAL *)LSRes0)
      CHKERR(Num_free_RHprimme((HREAL *)LSRes0, ctx));

   return 0;
}


/*******************************************************************************
 * Subroutine matrixMatvec_eigs - Computes y = B*x where B can be A'*A or A*A' or
 *    the augmented operator, A is the operator defined by
 *    primme_svds->matrixMatvec
 *
 * INPUT/OUPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * x          The input vectors
 * ldx        The leading dimension of x
 * y          The output vectors
 * ldy        The leading dimension of y
 * blockSize  Number of vectors in x and y
 * primme     The primme params structure of the eigensolver working on
 *            the equivalent eigenvalue problem to the original singular value
 *            problem.
 * ierr       Error code. 0 is ok; otherwise, error
 * 
 ******************************************************************************/

TEMPLATE_PLEASE
void matrixMatvec_eigs_Sprimme(void *x_, PRIMME_INT *ldx, void *y_,
      PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {

   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;
   primme_context ctx = primme_svds_get_context(primme_svds);
   int trans = 1, notrans = 0;
   SCALAR *x = (SCALAR*)x_, *y = (SCALAR*)y_, *aux;
   primme_svds_operator method = &primme_svds->primme == primme ?
      primme_svds->method : primme_svds->methodStage2;
   int i, bs;
   *ierr = 0;

   switch(method) {
      case primme_svds_op_AtA:
         CHKERRA(Num_malloc_Sprimme(primme_svds->mLocal *
                                          min(primme->maxBlockSize, *blockSize),
                       &aux, ctx),
               *ierr = 1);
         for (i=0, bs=min((*blockSize-i), primme->maxBlockSize); bs>0;
               i+= bs, bs=min((*blockSize-i), primme->maxBlockSize))
         {
            primme_svds->matrixMatvec(&x[*ldx * i], ldx, aux,
                  &primme_svds->mLocal, &bs, &notrans, primme_svds, ierr);
            CHKERRA(*ierr, /* do nothing */);
            primme_svds->matrixMatvec(aux, &primme_svds->mLocal,
                  &y[*ldy*i], ldy, &bs, &trans, primme_svds, ierr);
            CHKERRA(*ierr, /* do nothing */);
         }
         CHKERRA(Num_free_Sprimme(aux, ctx), *ierr = 1);
         break;
      case primme_svds_op_AAt:
         CHKERRA(Num_malloc_Sprimme(primme_svds->nLocal *
                                          min(primme->maxBlockSize, *blockSize),
                       &aux, ctx),
               *ierr = 1);
         for (i=0, bs=min((*blockSize-i), primme->maxBlockSize); bs>0;
               i+= bs, bs=min((*blockSize-i), primme->maxBlockSize))
         {
            primme_svds->matrixMatvec(&x[*ldx*i], ldx, aux,
                  &primme_svds->nLocal, &bs, &trans, primme_svds, ierr);
            CHKERRA(*ierr, /* do nothing */);
            primme_svds->matrixMatvec(aux, &primme_svds->nLocal,
                  &y[*ldy*i], ldy, &bs, &notrans, primme_svds, ierr);
            CHKERRA(*ierr, /* do nothing */);
         }
         CHKERRA(Num_free_Sprimme(aux, ctx), *ierr = 1);
         break;
      case primme_svds_op_augmented:
         primme_svds->matrixMatvec(&x[primme_svds->nLocal], ldx, y, ldy, blockSize,
               &trans, primme_svds, ierr);
         CHKERRA(*ierr, /* do nothing */);
         primme_svds->matrixMatvec(x, ldx, &y[primme_svds->nLocal],
               ldy, blockSize, &notrans, primme_svds, ierr);
         CHKERRA(*ierr, /* do nothing */);
         break;
      case primme_svds_op_none:
         break;
   }

   primme_svds_free_context(ctx);
}


/*******************************************************************************
 * Subroutine applyPreconditioner_eigs - Computes y = M*x where M is the
 *    preconditioner defined by primme_svds->applyPreconditioner
 *
 * INPUT/OUPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * x          The input vectors
 * ldx        The leading dimension of x
 * y          The output vectors
 * ldy        The leading dimension of y
 * blockSize  Number of vectors in x and y
 * primme     The primme params structure of the eigensolver working on
 *            the equivalent eigenvalue problem to the original singular value
 *            problem.
 * ierr       Error code. 0 is ok; otherwise, error
 * 
 ******************************************************************************/

TEMPLATE_PLEASE
void applyPreconditioner_eigs_Sprimme(void *x, PRIMME_INT *ldx, void *y,
      PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {

   primme_svds_params *primme_svds =
         (primme_svds_params *)primme->preconditioner;
   int method =
         (int)(&primme_svds->primme == primme ? primme_svds->method
                                              : primme_svds->methodStage2);

   primme_svds->applyPreconditioner(
         x, ldx, y, ldy, blockSize, &method, primme_svds, ierr);
}

STATIC int Num_scalInv_Smatrix(SCALAR *x, PRIMME_INT m, int n, PRIMME_INT ldx,
      HREAL *factors, primme_context ctx) {

   int i;
   HREAL norm, norm0, factor;

   assert(ldx >= m);
   for (i=0; i<n; i++) {
      if (factors[i] > 0.0 && 1.0L/factors[i] < HOST_MACHINE_MAX) {
         factor = factors[i];
      }
      else {
         norm0 = REAL_PART(Num_dot_Sprimme(m, &x[i*ldx], 1, &x[i*ldx], 1, ctx));
         CHKERR(globalSum_Rprimme_svds(&norm0, &norm, 1, ctx));
         factor = sqrt(norm);
      }
      Num_scal_Sprimme(m, 1.0/factor, &x[i*ldx], 1, ctx);
   }

   return 0;
}

STATIC int globalSum_Rprimme_svds(
      HREAL *sendBuf, HREAL *recvBuf, int count, primme_context ctx) {

   primme_svds_params *primme_svds = ctx.primme_svds;

   if (primme_svds && primme_svds->globalSumReal) {
      double t0 = primme_wTimer();

      /* Cast sendBuf and recvBuf */

      void *sendBuf0, *recvBuf0;
      CHKERR(Num_matrix_astype_Rprimme(sendBuf, 1, count, 1,
            primme_op_default, &sendBuf0, NULL, primme_svds->globalSumReal_type,
            1 /* alloc */, 1 /* copy */, ctx));
      if (sendBuf == recvBuf) {
         recvBuf0 = sendBuf0;
      } else {
         CHKERR(Num_matrix_astype_Rprimme(recvBuf, 1, count, 1,
            primme_op_default, &recvBuf0, NULL, primme_svds->globalSumReal_type,
            1 /* alloc */, 0 /* no copy */, ctx));
      }

      int ierr;
      CHKERRM((primme_svds->globalSumReal(
                     sendBuf0, recvBuf0, &count, primme_svds, &ierr),
                    ierr),
            PRIMME_USER_FAILURE, "Error returned by 'globalSumReal' %d", ierr);

      /* Copy back recvBuf0 */

      CHKERR(Num_matrix_astype_Rprimme(recvBuf0, 1, count, 1,
            primme_svds->globalSumReal_type, (void **)&recvBuf, NULL,
            primme_op_default, 0 /* no alloc */, 1 /* copy */, ctx));

      if (sendBuf != sendBuf0)
         CHKERR(Num_free_Sprimme((SCALAR *)sendBuf0, ctx));
      if (sendBuf != recvBuf && recvBuf != recvBuf0)
         CHKERR(Num_free_Sprimme((SCALAR *)recvBuf0, ctx));

      primme_svds->stats.numGlobalSum++;
      primme_svds->stats.timeGlobalSum += primme_wTimer() - t0;
      primme_svds->stats.volumeGlobalSum += count;
   }
   else {
      Num_copy_RHprimme(count, sendBuf, 1, recvBuf, 1, ctx);
   }

   return 0;
}

/*******************************************************************************
 * Subroutine compute_resNorm - This routine computes the residual norm of a
 *    given triplet (u,s,v):
 *
 *    sqrt(||A*v - s*u||^2 + ||A'*u - s*v||^2)
 *
 * NOTE:
 *    - The given u and v may not have norm one.
 *    - The computation requires two matvecs.
 * 
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * leftsvec     The approximate left singular vector
 * rightsvec    The approximate right singular vector
 * primme_svds  Structure containing various solver parameters
 *
 * OUTPUT PARAMETERS
 * ----------------------------------
 * rNorm        The norm of the residual vector
 * ierr         Error code
 ******************************************************************************/

STATIC int compute_resNorm(SCALAR *leftsvec, SCALAR *rightsvec, HREAL *rNorm,
      primme_context ctx) {

   primme_svds_params *primme_svds = ctx.primme_svds;
   PRIMME_INT nLocal = primme_svds->mLocal+primme_svds->nLocal;
   SCALAR *Atu;
   CHKERR(Num_malloc_Sprimme(nLocal, &Atu, ctx));
   SCALAR *Av = &Atu[primme_svds->nLocal];

   /* Av = A * v; Atu = A'u */

   CHKERR(matrixMatvecSVDS_Sprimme(leftsvec, primme_svds->mLocal, Atu,
         primme_svds->nLocal, 0, 1, 1 /* trans */, ctx));
   CHKERR(matrixMatvecSVDS_Sprimme(rightsvec, primme_svds->nLocal, Av,
         primme_svds->mLocal, 0, 1, 0 /* notrans */, ctx));

   /* ip[0] = ||v|| */
   /* ip[1] = ||u|| */
   /* ip[2] = u'*A*v = u'*Av */

   HREAL ip[3];
   ip[0] = REAL_PART(Num_dot_Sprimme(primme_svds->nLocal, rightsvec, 1,
            rightsvec, 1, ctx));
   ip[1] = REAL_PART(Num_dot_Sprimme(primme_svds->mLocal, leftsvec, 1,
            leftsvec, 1, ctx));
   ip[2] = REAL_PART(
         Num_dot_Sprimme(primme_svds->mLocal, leftsvec, 1, Av, 1, ctx));
   CHKERR(globalSum_Rprimme_svds(ip, ip, 3, ctx));

   ip[0] = sqrt(ip[0]);
   ip[1] = sqrt(ip[1]);
   HREAL sval = ip[2]/ip[0]/ip[1];

   /* If u'*A*v is negative, set rNorm as a large number */

   if (sval < -0.0) {
      *rNorm = HOST_MACHINE_MAX;
      CHKERR(Num_free_Sprimme(Atu, ctx));
      return 0;
   }

   /* Atu = A'*u/||u|| - sval*v/||v|| */

   Num_scal_Sprimme(primme_svds->nLocal, 1.0/ip[1], Atu, 1, ctx);
   Num_axpy_Sprimme(
         primme_svds->nLocal, -sval / ip[0], rightsvec, 1, Atu, 1, ctx);

   /* Av = A*v/||v|| - sval*u/||u|| */

   Num_scal_Sprimme(primme_svds->mLocal, 1.0/ip[0], Av, 1, ctx);
   Num_axpy_Sprimme(primme_svds->mLocal, -sval/ip[1], leftsvec, 1, Av, 1, ctx);

   /* resNorm = sqrt(||A*v - s*u||^2 + ||A'*u - s*v||^2) = norm([Atu; Av]) */

   HREAL normr0;
   normr0 = REAL_PART(Num_dot_Sprimme(nLocal, Atu, 1, Atu, 1, ctx));
   CHKERR(globalSum_Rprimme_svds(&normr0, rNorm, 1, ctx));
   *rNorm = sqrt(*rNorm);

   CHKERR(Num_free_Sprimme(Atu, ctx));
   return 0;
}

/*******************************************************************************
 * Subroutine default_convTestFun - This routine implements primme_params.
 *    convTestFun and returns an approximate triplet converged when           
 *    resNorm < eps * ||A||.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * sval         The approximate singular value 
 * leftsvec     The approximate left singular vector
 * rightsvec    The approximate right singular vector
 * rNorm        The norm of the residual vector
 * primme_svds  Structure containing various solver parameters
 *
 * OUTPUT PARAMETERS
 * ----------------------------------
 * isConv      if it isn't zero the approximate pair is marked as converged
 * ierr        error code
 ******************************************************************************/

STATIC void default_convTestFun(double *sval, void *leftsvec_, void *rightsvec_,
      double *rNorm, int *method, int *isConv, primme_svds_params *primme_svds,
      int *ierr) {

   (void)sval; /* unused parameter */
   const double aNorm = primme_svds->aNorm;
   SCALAR *leftsvec = (SCALAR*)leftsvec_, *rightsvec = (SCALAR*)rightsvec_;

   *isConv = *rNorm < max(primme_svds->eps, MACHINE_EPSILON * 3.16) * aNorm;

   /* If solving the augmented problem, the reported residual norm is an      */
   /* approximation. Recheck the convergence criterion with the actual        */
   /* residual norm when the convergence criterion is passed and the residual */
   /* vector norm is from the augmented problem.                              */

   if (*isConv && *method == primme_svds_op_augmented && leftsvec &&
         rightsvec) {

      HREAL rnorm;
      primme_context ctx = primme_svds_get_context(primme_svds);
      CHKERRA(compute_resNorm(leftsvec, rightsvec, &rnorm, ctx), *ierr = 1);
      primme_svds_free_context(ctx);

      *isConv = rnorm < max(primme_svds->eps, MACHINE_EPSILON * 3.16) * aNorm;
   }

   *ierr = 0;
}

/*******************************************************************************
 * Subroutine convTestFunATA - This routine implements primme_params.
 *    convTestFun and calls primme_svds.convTestFun when solving normal
 *    equations.
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

STATIC void convTestFunATA(double *eval, void *evec, double *rNorm, int *isConv,
      primme_params *primme, int *ierr) {

   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;
   primme_svds_operator method = &primme_svds->primme == primme ?
      primme_svds->method : primme_svds->methodStage2;
   assert(method == primme_svds_op_AtA || method == primme_svds_op_AAt);
   double aNorm = (primme->aNorm > 0.0) ?
      primme->aNorm : primme->stats.estimateLargestSVal;
   primme_context ctx = primme_svds_get_context(primme_svds);

   /* Check machine precision limit */

   double maxaNorm = max(primme->aNorm, primme->stats.estimateLargestSVal);
   double eps_matvec = MACHINE_EPSILON;
   CHKERRA(Num_machine_epsilon_Sprimme(
                 primme_svds->matrixMatvec_type, &eps_matvec),
         *ierr = 1);
   if (rNorm && *rNorm < eps_matvec * maxaNorm * 3.16) {
      *isConv = 1;
      *ierr = 0;
      return;
   }

   /* Update primme_svds->aNorm */

   double oldaNorm = primme_svds->aNorm;
   if (primme_svds->aNorm <= 0.0)
      primme_svds->aNorm = sqrt(aNorm);

   /* Call the callback */

   double sval = eval ? sqrt(fabs(*eval)) : 0.0;
   double srNorm = (rNorm && eval) ? *rNorm / sval : 0.0;
   CHKERRA(convTestFunSVDS_Sprimme(sval,
                 (method == primme_svds_op_AAt && evec) ? (SCALAR *)evec : NULL,
                 method == primme_svds_op_AAt,
                 (method == primme_svds_op_AtA && evec) ? (SCALAR *)evec : NULL,
                 method == primme_svds_op_AtA, srNorm, method, isConv, ctx),
         *ierr = 1);

   /* Restore aNorm */

   primme_svds->aNorm = oldaNorm;
   primme_svds_free_context(ctx);
}


/*******************************************************************************
 * Subroutine convTestFunAug - This routine implements primme_params.
 *    convTestFun and calls primme_svds.convTestFun when solving augmented
 *    problem.
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

STATIC void convTestFunAug(double *eval, void *evec, double *rNorm, int *isConv,
      primme_params *primme, int *ierr) {

   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;
   primme_svds_operator method = &primme_svds->primme == primme ?
      primme_svds->method : primme_svds->methodStage2;
   assert(method == primme_svds_op_augmented);
   double aNorm = (primme->aNorm > 0.0) ?
      primme->aNorm : primme->stats.estimateLargestSVal;
   primme_context ctx = primme_svds_get_context(primme_svds);

   /* NOTE: Don't check machine precision limit of the residual norm.      */
   /* Regardless of how small the residual is, we don't want to mark as    */
   /* converged pairs that approximate the null space of the augmented     */
   /* problem.                                                             */

   /* Update primme_svds->aNorm */

   double oldaNorm = primme_svds->aNorm;
   if (primme_svds->aNorm <= 0.0)
      primme_svds->aNorm = aNorm;

   /* Call the callback */

   double sval = eval ? fabs(*eval) : 0.0;
   double srNorm = rNorm ? *rNorm * sqrt(2.0) : 0.0;
   CHKERRA(convTestFunSVDS_Sprimme(sval,
                 evec ? &((SCALAR *)evec)[primme_svds->nLocal] : NULL,
                 1 /* vector given */, (SCALAR *)evec, 1 /* vector given */,
                 srNorm, method, isConv, ctx),
         *ierr = 1);

   /* Restore aNorm */

   primme_svds->aNorm = oldaNorm;
   primme_svds_free_context(ctx);
}


/*******************************************************************************
 * Subroutine default_monitor_svds - report iterations, #MV, residual norm,
 *    singular values, etc. at every inner/outer iteration and when some triplet
 *    converges.       
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * basisSvals   The approximate singular values of the basis
 * basisSize    The size of the basis
 * basisFlags   The state of every approximate triplet of the basis (see conv_flags)
 * iblock       Indices of the approximate triplet in the block
 * blockSize    The size of the block
 * basisNorms   The approximate residual norms of the triplet of the basis
 * numConverged The number of triplets converged in the basis and the locked triplets
 *              (this value isn't monotonic!)
 * lockedSvals  The locked singular values
 * numLocked    The number of triplets locked
 * lockedFlags  The state of each locked triplet (see conv_flags)
 * lockedNorms  The residual norms of the locked pairs
 * inner_its    The number of performed QMR iterations in the current correction equation
 * LSRes        The residual norm of the linear system at the current QMR iteration
 * event        The event reported
 * stage        0 for first stage and 1 for second stage
 * primme_svds  Structure containing various solver parameters and statistics
 *
 * OUTPUT
 * ------
 * err          Error code
 * 
 ******************************************************************************/

STATIC void default_monitor_svds(void *basisSvals_, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms_, int *numConverged,
      void *lockedSvals_, int *numLocked, int *lockedFlags, void *lockedNorms_,
      int *inner_its, void *LSRes_, const char *msg, double *time,
      primme_event *event, int *stage, primme_svds_params *primme_svds,
      int *err) {

   (void)basisSize;
   (void)basisFlags;

   XREAL *basisSvals = (XREAL*)basisSvals_, *basisNorms = (XREAL*)basisNorms_,
         *lockedSvals = (XREAL*)lockedSvals_, *lockedNorms = (XREAL*)lockedNorms_,
         *LSRes = (XREAL*)LSRes_;
   assert(event != NULL && primme_svds != NULL &&
          (stage != NULL || *event == primme_event_message));

   /* Only print report if this is proc zero or it is profiling */
   if (primme_svds->outputFile &&
         (primme_svds->procID == 0 || *event == primme_event_profile)) {
      switch(*event) {
         case primme_event_outer_iteration:
            assert(basisSvals && basisSize && basisFlags && iblock && blockSize
                  && basisNorms && numConverged);
            if (primme_svds->printLevel >= 3) {
               int i;  /* Loop variable */
               for (i=0; i < *blockSize; i++) {
                  fprintf(primme_svds->outputFile, 
                        "OUT %" PRIMME_INT_P " conv %d blk %d MV %" PRIMME_INT_P " Sec %E SV %13E |r| %.3E stage %d\n",
                        primme_svds->stats.numOuterIterations, *numConverged, i,
                        primme_svds->stats.numMatvecs,
                        primme_svds->stats.elapsedTime, (double)basisSvals[iblock[i]],
                        (double)basisNorms[iblock[i]], *stage+1);
               }
            }
            break;
         case primme_event_inner_iteration:
            assert(basisSize && iblock && basisNorms && inner_its && LSRes);
            (void)inner_its;
            if (primme_svds->printLevel >= 4) {
               fprintf(primme_svds->outputFile,
                     "INN MV %" PRIMME_INT_P " Sec %e Sval %e Lin|r| %.3e SV|r| %.3e stage %d\n",
                     primme_svds->stats.numMatvecs, primme_svds->stats.elapsedTime,
                     (double)basisSvals[iblock[0]], (double)*LSRes,
                     (double)basisNorms[iblock[0]], *stage+1);
            }
            break;
         case primme_event_restart:
            break;
         case primme_event_reset:
            break;
         case primme_event_converged:
            if ((*stage == 0 && primme_svds->printLevel >= 2)
                  || (primme_svds->printLevel >= 5))
               fprintf(primme_svds->outputFile, 
                     "#Converged %d sval[ %d ]= %e norm %e Mvecs %" PRIMME_INT_P " Time %g stage %d\n",
                     *numConverged, iblock[0], (double)basisSvals[iblock[0]],
                     (double)basisNorms[iblock[0]], primme_svds->stats.numMatvecs,
                     primme_svds->stats.elapsedTime, *stage+1);
            break;
         case primme_event_locked:
            if (primme_svds->printLevel >= 2) { 
               fprintf(primme_svds->outputFile, 
                     "Lock striplet[ %d ]= %e norm %.4e Mvecs %" PRIMME_INT_P " Time %.4e Flag %d stage %d\n",
                     *numLocked-1, (double)lockedSvals[*numLocked-1],
                     (double)lockedNorms[*numLocked-1], primme_svds->stats.numMatvecs,
                     primme_svds->stats.elapsedTime, lockedFlags[*numLocked-1],
                     *stage+1);
            }
            break;
         case primme_event_message:
            assert(msg != NULL);
            if (primme_svds->printLevel >= 2) { 
               fprintf(primme_svds->outputFile, 
                     "%s\n", msg);
            }
            break;
         case primme_event_profile:
            assert(msg != NULL && time != NULL);
            if (primme_svds->printLevel >= 3 && *time < 0.0) { 
               fprintf(primme_svds->outputFile, "entering in %s proc %d\n", msg, primme_svds->procID);
            }
            if (primme_svds->printLevel >= 2 && *time >= 0.0) { 
               fprintf(primme_svds->outputFile, "time for %s : %g proc %d\n", msg, *time, primme_svds->procID);
            }
            break;
      }
      fflush(primme_svds->outputFile);
   }
   *err = 0;
}


/*******************************************************************************
 * Subroutine monitor_single_stage - report iterations, #MV, residual norm,
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

STATIC void monitor_single_stage(void *basisEvals_, int *basisSize,
      int *basisFlags, int *iblock, int *blockSize, void *basisNorms_,
      int *numConverged, void *lockedEvals_, int *numLocked, int *lockedFlags,
      void *lockedNorms_, int *inner_its, void *LSRes_, const char *msg,
      double *time, primme_event *event, primme_params *primme, int *err) {

   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;
   primme_context ctx = primme_svds_get_context(primme_svds);

   int i;
   XREAL *basisEvals = (XREAL*)basisEvals_, *basisNorms = (XREAL*)basisNorms_,
         *lockedEvals = (XREAL*)lockedEvals_, *lockedNorms = (XREAL*)lockedNorms_,
         *LSRes = (XREAL*)LSRes_;
   assert(event != NULL && primme != NULL);

   HREAL *basisSvals, *basisSVNorms, *lockedSvals, *lockedSVNorms;

   CHKERRA(Num_malloc_RHprimme(
                 basisEvals && basisSize ? *basisSize : 0, &basisSvals, ctx),
         *err = 1);
   CHKERRA(Num_malloc_RHprimme(
                 basisEvals && basisSize ? *basisSize : 0, &basisSVNorms, ctx),
         *err = 1);
   CHKERRA(Num_malloc_RHprimme(
                 lockedEvals && numLocked ? *numLocked : 0, &lockedSvals, ctx),
         *err = 1);
   CHKERRA(Num_malloc_RHprimme(lockedEvals && numLocked ? *numLocked : 0,
                 &lockedSVNorms, ctx),
         *err = 1);

   if (primme_svds->method == primme_svds_op_AtA
         || primme_svds->method == primme_svds_op_AAt) {
      /* sval = sqrt(abs(eval)) and SVrnorm = rnorm/sval */

      if (basisEvals && basisSize) {
         for (i = 0; i < *basisSize; i++) {
            basisSvals[i] = (HREAL)sqrt(fabs((HREAL)basisEvals[i]));
            basisSVNorms[i] =
                  (basisSvals[i] > 0.0 ? (XREAL)((HREAL)basisNorms[i] /
                                                 (HREAL)basisSvals[i])
                                       : basisNorms[i]);
         }
      }

      if (lockedEvals && numLocked) {
         for (i = 0; i < *numLocked; i++) {
            lockedSvals[i] = (HREAL)sqrt(fabs((HREAL)lockedEvals[i]));
            lockedSVNorms[i] =
                  (lockedSvals[i] > 0.0 ? (XREAL)((HREAL)lockedNorms[i] /
                                                  (HREAL)lockedSvals[i])
                                        : lockedNorms[i]);
         }
      }
   }
   else if (primme_svds->method == primme_svds_op_augmented) {
      /* SVrnorm = rnorm/sqrt(2) */

      if (basisEvals && basisNorms && basisSize) {
         for (i = 0; i < *basisSize; i++) {
            basisSVNorms[i] = basisNorms[i] / sqrt(2.0);
         }
      }

      if (lockedEvals && numLocked) {
         for (i = 0; i < *numLocked; i++) {
            lockedSVNorms[i] = lockedNorms[i] / sqrt(2.0);
         }
      }
   }

   /* Prefix msg with Sprimme0 if the event is a profile */

   char *new_msg = NULL;
   if (*event == primme_event_profile && msg) {
      int len = 12 + strlen(msg);
      if (MALLOC_PRIMME(len, &new_msg) == 0) {
         snprintf(new_msg, len, "~Sprimme0%s", msg);
         msg = new_msg;
      }
   }

   /* Record performance measurements */ 

   primme_svds_stats stats = primme_svds->stats;
   UPDATE_STATS(primme_svds->stats, +=, primme->stats);

   /* Call the user function report */

   CHKERRA(monitorFunSVDS_Sprimme(basisSvals, basisSize ? *basisSize : 0,
                 basisFlags, iblock, blockSize ? *blockSize : 0, basisSVNorms,
                 numConverged ? *numConverged : 0, lockedSvals,
                 numLocked ? *numLocked : 0, lockedFlags, lockedSVNorms,
                 inner_its ? *inner_its : 0, LSRes ? *LSRes : 0.0, msg,
                 time ? *time : 0.0, *event, 0 /* stage 0 */,
                 -1.0 /* don't update elapsedTime */, ctx),
         *err = 1);
   primme_svds->stats = stats; /* restore original values */

   CHKERRA(Num_free_RHprimme(basisSvals, ctx), *err = 1);
   CHKERRA(Num_free_RHprimme(basisSVNorms, ctx), *err = 1);
   CHKERRA(Num_free_RHprimme(lockedSvals, ctx), *err = 1);
   CHKERRA(Num_free_RHprimme(lockedSVNorms, ctx), *err = 1);
   if (new_msg) free(new_msg);
   primme_svds_free_context(ctx);
}


/*******************************************************************************
 * Subroutine monitor_stage1 - translate monitored information from eigenvalues
 *    to singular values and call the monitor in primme_svds. Notice that
 *    because there is a second stage, the locked pairs at this stage are
 *    reported as converged. 
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

STATIC void monitor_stage1(void *basisEvals_, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms_, int *numConverged,
      void *lockedEvals_, int *numLocked, int *lockedFlags, void *lockedNorms_,
      int *inner_its, void *LSRes_, const char *msg, double *time,
      primme_event *event, primme_params *primme, int *err) {

   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;
   primme_context ctx = primme_svds_get_context(primme_svds);

   XREAL *basisEvals = (XREAL*)basisEvals_, *basisNorms = (XREAL*)basisNorms_,
         *lockedEvals = (XREAL*)lockedEvals_, *lockedNorms = (XREAL*)lockedNorms_,
         *LSRes = (XREAL*)LSRes_;
   assert(event != NULL && primme != NULL);

   /* Ignore the converged events if locking is active and printLevel <= 4 */

   if (*event == primme_event_converged && primme->locking
         && primme->printLevel <= 4) {
      *err = 0;
      primme_svds_free_context(ctx);
      return;
   }

   /* Show locked pairs as converged pairs of the basis */

   int numLocked0 = lockedEvals&&numLocked?*numLocked:0;
   int basisSize0 = (basisEvals&&basisSize?*basisSize:0) + numLocked0;
   HREAL *basisSvals = NULL, *basisSVNorms = NULL;
   int *basisSVFlags = NULL, *iblockSV = NULL;
   CHKERRA(Num_malloc_RHprimme(basisSize0, &basisSvals, ctx), *err = 1);
   CHKERRA(Num_malloc_RHprimme(basisSize0, &basisSVNorms, ctx), *err = 1);
   CHKERRA(Num_malloc_iprimme(basisSize0, &basisSVFlags, ctx), *err = 1);
   CHKERRA(Num_malloc_iprimme(
                 blockSize && *blockSize > 0 ? *blockSize : 1, &iblockSV, ctx),
         *err = 1);
   int numConvergedSV = (numConverged?*numConverged:numLocked0);

   assert(primme_svds->method == primme_svds_op_AtA
         || primme_svds->method == primme_svds_op_AAt);

   /* sval = sqrt(abs(eval)) and SVrnorm = rnorm/sval */

   int i, j = 0;
   if (lockedEvals && numLocked) {
      for (i = 0; i < *numLocked; i++, j++) {
         basisSvals[j] = (HREAL)sqrt(fabs((HREAL)lockedEvals[i]));
         basisSVNorms[j] = (basisSvals[i] > 0.0 ? lockedNorms[i] / basisSvals[i]
                                                : lockedNorms[i]);
         basisSVFlags[j] = lockedFlags[i];
      }
   }

   if (basisEvals && basisSize) {
      for (i = 0; i < *basisSize; i++, j++) {
         basisSvals[j] = (HREAL)sqrt(fabs((HREAL)basisEvals[i]));
         basisSVNorms[j] = (basisSvals[i] > 0.0 ? basisNorms[i] / basisSvals[i]
                                                : basisNorms[i]);
         basisSVFlags[j] = (basisFlags ? basisFlags[i] : UNCONVERGED);
      }
   }

   if (iblock && blockSize) for (i=0; i<*blockSize; i++) {
      iblockSV[i] = iblock[i] + numLocked0;
   }

   primme_event eventSV = *event;
   if (*event == primme_event_locked) {
      eventSV = primme_event_converged;
      iblockSV[0] = *numLocked-1;
   }

   /* Prefix msg with Sprimme0 if the event is a profile */

   char *new_msg = NULL;
   if (*event == primme_event_profile && msg) {
      int len = 12 + strlen(msg);
      if (MALLOC_PRIMME(len, &new_msg) == 0) {
         snprintf(new_msg, len, "~Sprimme0%s", msg);
         msg = new_msg;
      }
   }

   /* Record performance measurements */ 

   primme_svds_stats stats = primme_svds->stats;
   UPDATE_STATS(primme_svds->stats, +=, primme->stats);

   /* Call the user function report */

   CHKERRA(
         monitorFunSVDS_Sprimme(basisSvals, basisSize0, basisSVFlags, iblockSV,
               blockSize ? *blockSize : 0, basisSVNorms, numConvergedSV, NULL,
               0, NULL, NULL, inner_its ? *inner_its : 0, LSRes ? *LSRes : 0.0,
               msg, time ? *time : 0.0, eventSV, 0 /* stage 0 */,
               -1.0 /* don't update elapsedTime */, ctx),
         *err = 1);
   primme_svds->stats = stats; /* restore original values */

   CHKERRA(Num_free_RHprimme(basisSvals, ctx), *err = 1);
   CHKERRA(Num_free_RHprimme(basisSVNorms, ctx), *err = 1);
   CHKERRA(Num_free_iprimme(basisSVFlags, ctx), *err = 1);
   CHKERRA(Num_free_iprimme(iblockSV, ctx), *err = 1);
   if (new_msg) free(new_msg);
   primme_svds_free_context(ctx);
}

/*******************************************************************************
 * Subroutine monitor_stage2 - report iterations, #MV, residual norm,
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

STATIC void monitor_stage2(void *basisEvals_, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms_, int *numConverged,
      void *lockedEvals_, int *numLocked, int *lockedFlags, void *lockedNorms_,
      int *inner_its, void *LSRes_, const char *msg, double *time,
      primme_event *event, primme_params *primme, int *err) {

   XREAL *basisEvals = (XREAL*)basisEvals_, *basisNorms = (XREAL*)basisNorms_,
         *lockedEvals = (XREAL*)lockedEvals_, *lockedNorms = (XREAL*)lockedNorms_,
         *LSRes = (XREAL*)LSRes_;
   assert(event != NULL && primme != NULL);
   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;
   primme_context ctx = primme_svds_get_context(primme_svds);
   int i;

   /* Included the converged triplets after the first stage as locked */

   int numLockedExtra = lockedEvals&&numLocked ?
      primme_svds->numSvals - primme->numEvals : 0;
   int numLockedSV = (lockedEvals&&numLocked?*numLocked:0) + numLockedExtra;
   int basisSize0 = (basisEvals&&basisSize?*basisSize:0);
   HREAL *basisSvals, *basisSVNorms, *lockedSvals, *lockedSVNorms;
   int *lockedSVFlags;
   CHKERRA(Num_malloc_RHprimme(basisSize0, &basisSvals, ctx), *err = 1);
   CHKERRA(Num_malloc_RHprimme(basisSize0, &basisSVNorms, ctx), *err = 1);
   CHKERRA(Num_malloc_RHprimme(numLockedSV, &lockedSvals, ctx), *err = 1);
   CHKERRA(Num_malloc_RHprimme(numLockedSV, &lockedSVNorms, ctx), *err = 1);
   CHKERRA(Num_malloc_iprimme(numLockedSV, &lockedSVFlags, ctx), *err = 1);

   /* SVrnorm = rnorm/sqrt(2) */

   for (i=0; i<basisSize0; i++) {
      basisSvals[i] = basisEvals[i];
      basisSVNorms[i] = basisNorms[i]/sqrt(2.0);
   }

   // NOTE: Always lockedEvals and lockedNorms are either NULL or a position
   // of the vectors svals and resNorms in wrapper_svds_Sprimme. So we can
   // access to the singular values and residual norms at stage 0.

   lockedEvals -= numLockedExtra;
   lockedNorms -= numLockedExtra;

   for (i=0; i<numLockedExtra; i++) {
      lockedSvals[i] = lockedEvals[i];
      lockedSVNorms[i] = lockedNorms[i];
      lockedSVFlags[i] = CONVERGED;
   }

   for (i=numLockedExtra; i<numLockedSV; i++) {
      lockedSvals[i] = lockedEvals[i];
      lockedSVNorms[i] = lockedNorms[i] / sqrt(2.0);
      lockedSVFlags[i] = lockedFlags[i - numLockedExtra];
   }

   /* Prefix msg with Sprimme1 if the event is a profile */

   char *new_msg = NULL;
   if (*event == primme_event_profile && msg) {
      int len = 12 + strlen(msg);
      if (MALLOC_PRIMME(len, &new_msg) == 0) {
         snprintf(new_msg, len, "~Sprimme1%s", msg);
         msg = new_msg;
      }
   }

   /* Record performance measurements */ 

   primme_svds_stats stats = primme_svds->stats;
   UPDATE_STATS(primme_svds->stats, +=, primme->stats);

   /* Call the user function report */

   CHKERRA(monitorFunSVDS_Sprimme(basisSvals, basisSize0,
                 basisFlags, iblock, blockSize ? *blockSize : 0, basisSVNorms,
                 numConverged ? *numConverged : 0, lockedSvals, numLockedSV,
                 lockedSVFlags, lockedSVNorms, inner_its ? *inner_its : 0,
                 LSRes ? *LSRes : 0.0, msg, time ? *time : 0.0, *event,
                 1 /* stage 2 */, -1.0 /* don't update elapsedTime */, ctx),
         *err = 1);
   primme_svds->stats = stats; /* restore original values */

   CHKERRA(Num_free_RHprimme(basisSvals, ctx), *err = 1);
   CHKERRA(Num_free_RHprimme(basisSVNorms, ctx), *err = 1);
   CHKERRA(Num_free_RHprimme(lockedSvals, ctx), *err = 1);
   CHKERRA(Num_free_RHprimme(lockedSVNorms, ctx), *err = 1);
   CHKERRA(Num_free_iprimme(lockedSVFlags, ctx), *err = 1);
   if (new_msg) free(new_msg);
   primme_svds_free_context(ctx);
}

#endif /* SUPPORTED_TYPE */
