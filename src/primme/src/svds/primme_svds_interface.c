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
 **********************************************************************
 * File: primme_svds_interface.c
 *
 * Purpose - Contains interface functions to PRIMME SVDS named primme_svds_*.
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../svds/primme_svds_interface.c"
#endif


#include <stdlib.h>   /* mallocs, free */
#include <stdio.h>    
#include <math.h>    
#include <string.h>    
#include <limits.h>    
#include "numerical.h"
#include "primme_svds_interface.h"
#include "primme_interface.h"

/* Only define these functions ones */
#ifdef USE_DOUBLE
#include "notemplate.h"

typedef union {
   PRIMME_INT int_v;
   void (*matFunc_v) (void*,PRIMME_INT*,void*,PRIMME_INT*,int*,int*,struct primme_svds_params*,int*);
   void (*globalSumRealFunc_v) (void *,void *,int *,struct primme_svds_params*,int*);
   void (*broadcastRealFunc_v) (void *,int *,struct primme_svds_params*,int*);
   void (*convTestFun_v)(double *sval, void *leftsvec, void *rightsvec,
         double *rNorm, int *method, int *isconv,
         struct primme_svds_params *primme, int *ierr);
   void (*monitorFun_v)(void *basisSvals, int *basisSize, int *basisFlags,
         int *iblock, int *blockSize, void *basisNorms, int *numConverged,
         void *lockedSvals, int *numLocked, int *lockedFlags,
         void *lockedNorms, int *inner_its, void *LSRes, const char *msg,
         double *time, primme_event *event, int *stage,
         struct primme_svds_params *primme_svds, int *err);
} svds_value_t;
typedef void *ptr_v;
typedef const char *str_v;

static void copy_params_from_svds(primme_svds_params *primme_svds, int stage);
static void globalSumRealSvds(void *sendBuf, void *recvBuf, int *count, 
                         primme_params *primme, int *ierr);
static void broadcastRealSvds(
      void *buffer, int *count, primme_params *primme, int *ierr);

/*****************************************************************************
 * Initialize handles also the allocation of primme_svds structure 
 *****************************************************************************/
primme_svds_params * primme_svds_params_create(void) {

   primme_svds_params *primme_svds = NULL;
   if (MALLOC_PRIMME(1, &primme_svds) == 0)
      primme_svds_initialize(primme_svds);
    return primme_svds;
}

/*****************************************************************************
 *  * Free the internally allocated work arrays of the primme_svds structure 
 *   *****************************************************************************/
int primme_svds_params_destroy(primme_svds_params *primme_svds) {
    free(primme_svds);
    return 0;
}

/*******************************************************************************
 * Subroutine primme_svds_initialize - Set primme_svds_params members to default
 *    values.
 * 
 * INPUT/OUTPUT PARAMETERS
 * ----------------------------------
 * primme_svds  Structure containing various solver parameters and statistics
 *
 ******************************************************************************/

void primme_svds_initialize(primme_svds_params *primme_svds) {

   /* Essential parameters */
   primme_svds->m                       = 0;
   primme_svds->n                       = 0;
   primme_svds->numSvals                = 1;
   primme_svds->target                  = primme_svds_largest;
   primme_svds->method                  = primme_svds_op_none;
   primme_svds->methodStage2            = primme_svds_op_none;

   /* Shifts for primme_svds_augmented method */
   primme_svds->numTargetShifts         = 0;
   primme_svds->targetShifts            = NULL;

   /* Parallel computing parameters */
   primme_svds->numProcs                = 1;
   primme_svds->procID                  = 0;
   primme_svds->mLocal                  = -1;
   primme_svds->nLocal                  = -1;
   primme_svds->commInfo                = NULL;
   primme_svds->globalSumReal           = NULL;
   primme_svds->globalSumReal_type      = primme_op_default;
   primme_svds->broadcastReal           = NULL;
   primme_svds->broadcastReal_type      = primme_op_default;
   primme_svds->internalPrecision       = primme_op_default;

   /* Use these pointers to provide matrix/preconditioner */
   primme_svds->matrix                  = NULL;
   primme_svds->preconditioner          = NULL;

   /* Matvec and preconditioner */
   primme_svds->matrixMatvec            = NULL;
   primme_svds->matrixMatvec_type       = primme_op_default;
   primme_svds->applyPreconditioner     = NULL;
   primme_svds->applyPreconditioner_type= primme_op_default;

   /* Other important parameters users may set */
   primme_svds->aNorm                   = 0.0L;
   primme_svds->eps                     = 0.0;
   primme_svds->precondition            = -1;
   primme_svds->initSize                = 0;
   primme_svds->maxBasisSize            = 0;
   primme_svds->maxBlockSize            = 0;
   primme_svds->maxMatvecs              = INT_MAX;
   primme_svds->printLevel              = 1;
   primme_svds->outputFile              = stdout;
   primme_svds->locking                 = -1;
   primme_svds->numOrthoConst           = 0;

   /* Reporting performance */
   primme_svds->stats.numOuterIterations            = 0; 
   primme_svds->stats.numRestarts                   = 0;
   primme_svds->stats.numMatvecs                    = 0;
   primme_svds->stats.numPreconds                   = 0;
   primme_svds->stats.numGlobalSum                  = 0;
   primme_svds->stats.volumeGlobalSum               = 0;
   primme_svds->stats.numBroadcast                  = 0;
   primme_svds->stats.volumeBroadcast               = 0;
   primme_svds->stats.numOrthoInnerProds            = 0.0;
   primme_svds->stats.elapsedTime                   = 0.0;
   primme_svds->stats.timeMatvec                    = 0.0;
   primme_svds->stats.timePrecond                   = 0.0;
   primme_svds->stats.timeOrtho                     = 0.0;
   primme_svds->stats.timeGlobalSum                 = 0.0;
   primme_svds->stats.timeBroadcast                 = 0.0;

   /* Internally used variables */
   primme_svds->iseed[0] = -1;   /* To set iseed, we first need procID           */ 
   primme_svds->iseed[1] = -1;   /* Thus we set all iseeds to -1                 */
   primme_svds->iseed[2] = -1;   /* Unless users provide their own iseeds,       */
   primme_svds->iseed[3] = -1;   /* PRIMME will set thse later uniquely per proc */
   primme_svds->convTestFun             = NULL;
   primme_svds->convTestFun_type        = primme_op_default;
   primme_svds->convtest                = NULL;
   primme_svds->monitorFun              = NULL;
   primme_svds->monitorFun_type         = primme_op_default;
   primme_svds->monitor                 = NULL;
   primme_svds->queue                   = NULL;
   primme_svds->profile                 = NULL;

   primme_initialize(&primme_svds->primme);
   primme_initialize(&primme_svds->primmeStage2);
 
}

/******************************************************************************
 * Subroutine primme_svds_set_method - Set the approach to solve the singular
 *    value problem.
 *
 * INPUT PARAMETERS
 * ----------------
 *    method   singular value method, one of:
 *
 *       primme_svds_default, currently set as primme_svds_hybrid.
 *       primme_svds_normalequations, compute the eigenvectors of A'*A or A*A'.
 *       primme_svds_augmented|, compute the eigenvectors of the augmented
 *          matrix, [zeros() A'; A zeros()].
 *       primme_svds_hybrid, start with primme_svds_normalequations; use the
 *       resulting approximate singular vectors as initial vectors for
 *       primme_svds_augmented if the required accuracy was not achieved.
 *
 *    methodStage1: preset method to compute the eigenpairs at the first stage.
 *
 *    methodStage2: preset method to compute the eigenpairs with the second
 *       stage of primme_svds_hybrid.
 *
 * INPUT/OUTPUT PARAMATERS
 * -----------------------
 *    primme_svds   parameters structure 
 *
 ******************************************************************************/

int primme_svds_set_method(primme_svds_preset_method method,
      primme_preset_method methodStage1, primme_preset_method methodStage2,
      primme_svds_params *primme_svds) {

   /* Set method and methodStage2 in primme_svds_params */
   switch(method) {
   case primme_svds_default:
   case primme_svds_hybrid:
      primme_svds->method = (primme_svds->n <= primme_svds->m ?
            primme_svds_op_AtA : primme_svds_op_AAt);
      primme_svds->methodStage2 = primme_svds_op_augmented;
      break;
   case primme_svds_normalequations:
      primme_svds->method = (primme_svds->n <= primme_svds->m ?
            primme_svds_op_AtA : primme_svds_op_AAt);
      primme_svds->methodStage2 = primme_svds_op_none;
      break;
   case primme_svds_augmented:
      primme_svds->method = primme_svds_op_augmented;
      primme_svds->methodStage2 = primme_svds_op_none;
      break;
   }

   /* Setup underneath eigensolvers based on primme_svds configuration */
   primme_svds_set_defaults(primme_svds);

   /* Set method for the first stage */
   primme_set_method(methodStage1, &primme_svds->primme);

   /* Set method for the second state */
   if (methodStage2 == PRIMME_DEFAULT_METHOD
         && primme_svds->target != primme_svds_largest)
      methodStage2 = PRIMME_JDQMR;
   if (primme_svds->methodStage2 != primme_svds_op_none)
      primme_set_method(methodStage2, &primme_svds->primmeStage2);

   return 0;
}

/******************************************************************************
 * Subroutine primme_svds_set_defaults - Set valid values for options that still
 *    has the initial invalid value set in primme_svds_initialize.
 *
 * INPUT/OUTPUT PARAMETERS
 * ----------------------------------
 *    primme_svds   parameters structure 
 *
 ******************************************************************************/

void primme_svds_set_defaults(primme_svds_params *primme_svds) {

   /* Set svds method if none set */
   if (primme_svds->method == primme_svds_op_none) {
      primme_svds_set_method(primme_svds_default, PRIMME_DEFAULT_METHOD,
            PRIMME_DEFAULT_METHOD, primme_svds);
   }

   /* Copy values set in primme_svds to the first stage underneath eigensolver */
   copy_params_from_svds(primme_svds, 0);

   if (primme_svds->methodStage2 != primme_svds_op_none) {
      /* Copy values set in primme_svds to the second stage underneath eigensolver */
      copy_params_from_svds(primme_svds, 1);
   }
}

/*******************************************************************************
 * Subroutine copy_params_from_svds - Transfer some options in primme_svds to
 *    the primme_params for the given stage.
 * 
 * INPUT/OUTPUT PARAMETERS
 * ----------------------------------
 * stage        First stage (0) or second stage (1)
 * primme_svds  Structure containing various solver parameters and statistics
 *
 ******************************************************************************/

static void copy_params_from_svds(primme_svds_params *primme_svds, int stage) {
   primme_params *primme;
   primme_svds_operator method;

   primme = stage == 0 ? &primme_svds->primme : &primme_svds->primmeStage2;
   method = stage == 0 ? primme_svds->method : primme_svds->methodStage2;

   if (method == primme_svds_op_none) {
      primme->maxMatvecs = 1;
      return;
   }

   /* -----------------------------------------------*/
   /* Set important parameters for primme structure  */
   /* -----------------------------------------------*/
   primme->numEvals = primme_svds->numSvals;
   if (primme_svds->aNorm > 0.0) {
      switch(method) {
      case primme_svds_op_AtA:
      case primme_svds_op_AAt:
         primme->aNorm = primme_svds->aNorm*primme_svds->aNorm;
         break;
      case primme_svds_op_augmented:
         primme->aNorm = primme_svds->aNorm*sqrt(2.0);
         break;
      case primme_svds_op_none:
         break;
      }
   }
   primme->eps = primme_svds->eps;
   primme->initSize = primme_svds->initSize;
   if (primme_svds->maxBasisSize > 0)
      primme->maxBasisSize = primme_svds->maxBasisSize;
   if (primme_svds->maxBlockSize > 0)
      primme->maxBlockSize = primme_svds->maxBlockSize;
   primme->maxMatvecs = primme_svds->maxMatvecs;
   primme->printLevel = primme_svds->printLevel;
   primme->outputFile = primme_svds->outputFile;
   primme->numOrthoConst = primme_svds->numOrthoConst;

   /* ---------------------------------------------- */
   /* Set some parameters only for parallel programs */
   /* ---------------------------------------------- */
   if (primme_svds->numProcs > 1) {
      primme->procID = primme_svds->procID;
      primme->numProcs = primme_svds->numProcs;
      primme->commInfo = primme_svds->commInfo;
   }
   if (primme_svds->globalSumReal != NULL) {
      primme->globalSumReal = globalSumRealSvds;
   }
   if (primme_svds->broadcastReal != NULL) {
      primme->broadcastReal = broadcastRealSvds;
   }

   switch(method) {
   case primme_svds_op_AtA:
      primme->n = primme_svds->n;
      if (primme->nLocal == -1 && primme_svds->nLocal != -1)
         primme->nLocal = primme_svds->nLocal;
      break;
   case primme_svds_op_AAt:
      primme->n = primme_svds->m;
      if (primme->nLocal == -1 && primme_svds->mLocal != -1)
         primme->nLocal = primme_svds->mLocal;
      break;
   case primme_svds_op_augmented:
      primme->n = primme_svds->m + primme_svds->n;
      if (primme->nLocal == -1 && primme_svds->mLocal != -1
            && primme_svds->nLocal != -1)
         primme->nLocal = primme_svds->mLocal + primme_svds->nLocal;
      break;
   case primme_svds_op_none:
      break;
   }

   switch (primme_svds->target) {
   case primme_svds_largest:
      primme->target = primme_largest;
      break;
   case primme_svds_smallest:
      primme->target = (method == primme_svds_op_augmented) ?
         primme_closest_geq : primme_smallest;
      break;
   case primme_svds_closest_abs:
      primme->target = primme_closest_abs;
      primme->numTargetShifts = primme_svds->numTargetShifts;
      break;
   }

   if (stage == 1 && primme->initBasisMode == primme_init_default) {
      primme->initBasisMode = primme_init_user;
   }

   if (((method == primme_svds_op_augmented && 
         primme_svds->target != primme_svds_largest) ||
        primme_svds->target == primme_svds_closest_abs) &&
         primme->projectionParams.projection == primme_proj_default) {
      /* NOTE: refined extraction seems to work better than RR */
      primme->projectionParams.projection = primme_proj_refined;
   }

   if (primme_svds->locking >= 0) {
      primme->locking = primme_svds->locking;
   }

   if (primme_svds->precondition >= 0) {
      primme->correctionParams.precondition = primme_svds->precondition;
   }
   else if (primme->correctionParams.precondition < 0) {
      primme->correctionParams.precondition = primme_svds->applyPreconditioner ? 1 : 0;
   }

}

/*******************************************************************************
 * Subroutine primme_svds_display_params - Displays the current configuration of
 *    primme svds data structure.
 * 
 * INPUT PARAMETERS
 * ----------------------------------
 * primme_svds  Structure containing various solver parameters and statistics
 *
 ******************************************************************************/

void primme_svds_display_params(primme_svds_params primme_svds) {

   int i;
   FILE *outputFile = primme_svds.outputFile;

#define PRINT(P,L) fprintf(outputFile, "primme_svds." #P " = " #L "\n", primme_svds. P);
#define PRINTIF(P,V) if (primme_svds. P == V) fprintf(outputFile, "primme_svds." #P " = " #V "\n");
#define PRINTParams(P,S,L) fprintf(outputFile, "primme_svds." #P "." #S " = " #L "\n", \
                                    primme_svds. P ## Params. S);
#define PRINTParamsIF(P,S,V) if (primme_svds. P ## Params. S == V) \
                                 fprintf(outputFile, "primme_svds." #P "." #S " = " #V "\n");
#define PRINT_PRIMME_INT(P) fprintf(outputFile, "primme_svds." #P " = %" PRIMME_INT_P "\n", primme_svds. P);
 
fprintf(outputFile, "// ---------------------------------------------------\n"
                    "//            primme_svds configuration               \n"
                    "// ---------------------------------------------------\n");

   PRINT_PRIMME_INT(m);
   PRINT_PRIMME_INT(n);
   PRINT_PRIMME_INT(mLocal);
   PRINT_PRIMME_INT(nLocal);
   PRINT(numProcs, %d);
   PRINT(procID, %d);

   fprintf(outputFile, "\n// Output and reporting\n");
   PRINT(printLevel, %d);

   fprintf(outputFile, "\n// Solver parameters\n");
   PRINT(numSvals, %d);
   PRINT(aNorm, %e);
   PRINT(eps, %e);
   PRINT(maxBasisSize, %d);
   PRINT(maxBlockSize, %d);
   PRINT_PRIMME_INT(maxMatvecs);

   PRINTIF(target, primme_svds_smallest);
   PRINTIF(target, primme_svds_largest);
   PRINTIF(target, primme_svds_closest_abs);

   PRINT(numTargetShifts, %d);
   if (primme_svds.numTargetShifts > 0) {
      fprintf(outputFile, "primme_svds.targetShifts =");
      for (i=0; i<primme_svds.numTargetShifts;i++) {
         fprintf(outputFile, " %e",primme_svds.targetShifts[i]);
      }
      fprintf(outputFile, "\n");
   }

   PRINT(locking, %d);
   PRINT(initSize, %d);
   PRINT(numOrthoConst, %d);
   fprintf(outputFile, "primme_svds.iseed =");
   for (i=0; i<4;i++) {
      fprintf(outputFile, " %" PRIMME_INT_P, primme_svds.iseed[i]);
   }
   fprintf(outputFile, "\n");

   PRINT(precondition, %d);

   PRINTIF(method, primme_svds_op_none);
   PRINTIF(method, primme_svds_op_AtA);
   PRINTIF(method, primme_svds_op_AAt);
   PRINTIF(method, primme_svds_op_augmented);

   PRINTIF(methodStage2, primme_svds_op_none);
   PRINTIF(methodStage2, primme_svds_op_AtA);
   PRINTIF(methodStage2, primme_svds_op_AAt);
   PRINTIF(methodStage2, primme_svds_op_augmented);

   PRINTIF(internalPrecision, primme_op_half);
   PRINTIF(internalPrecision, primme_op_float);
   PRINTIF(internalPrecision, primme_op_double);
   PRINTIF(internalPrecision, primme_op_quad);


   if (primme_svds.method != primme_svds_op_none) {
      fprintf(outputFile, "\n"
                          "// ---------------------------------------------------\n"
                          "//            1st stage primme configuration          \n"
                          "// ---------------------------------------------------\n");
      primme_svds.primme.outputFile = outputFile;
      primme_display_params_prefix("primme", primme_svds.primme);
   }

   if (primme_svds.methodStage2 != primme_svds_op_none) {
      fprintf(outputFile, "\n"
                          "// ---------------------------------------------------\n"
                          "//            2st stage primme configuration          \n"
                          "// ---------------------------------------------------\n");
      primme_svds.primmeStage2.outputFile = outputFile;
      primme_display_params_prefix("primmeStage2", primme_svds.primmeStage2);
   }
   fflush(outputFile);

}


/*******************************************************************************
 * Subroutine primme_svds_free - Free memory resources allocated by Sprimme_svds
 * 
 * INPUT/OUTPUT PARAMETERS
 * ----------------------------------
 * primme_svds  Structure containing various solver parameters and statistics
 *
 ******************************************************************************/

void primme_svds_free(primme_svds_params *primme) {
   (void)primme;
   /* No function */    
}

/*******************************************************************************
 * Subroutine globalSumRealSvds - implementation of primme_params' globalSumReal
 *    that uses the callback defined in primme_svds_params.
 * 
 ******************************************************************************/

static void globalSumRealSvds(void *sendBuf, void *recvBuf, int *count, 
                         primme_params *primme, int *ierr) {
   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;
   primme_svds->globalSumReal(sendBuf, recvBuf, count, primme_svds, ierr);
}

/*******************************************************************************
 * Subroutine broadcastRealSvds - implementation of primme_params' broadcastReal
 *    that uses the callback defined in primme_svds_params.
 * 
 ******************************************************************************/

static void broadcastRealSvds(
      void *buffer, int *count, primme_params *primme, int *ierr) {
   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;
   primme_svds->broadcastReal(buffer, count, primme_svds, ierr);
}

/*******************************************************************************
 * Subroutine primme_svds_get_member - get the value of a parameter in
 *    primme_svds_params
 * 
 * INPUT PARAMETERS
 * ----------------
 * primme_svds  Structure containing various solver parameters and statistics
 * label   reference to the parameter
 *
 * OUTPUT PARAMETERS
 * -----------------
 * value   value of the parameter
 *
 * RETURN
 * ------
 * error code  zero if ok
 *
 ******************************************************************************/

int primme_svds_get_member(primme_svds_params *primme_svds,
      primme_svds_params_label label, void *value) {

   int i;
   svds_value_t *v = (svds_value_t*)value;

   switch(label) {
      case PRIMME_SVDS_primme :
         *(ptr_v*)value = &primme_svds->primme;
         break;
      case PRIMME_SVDS_primmeStage2 :
         *(ptr_v*)value = &primme_svds->primmeStage2;
         break;
      case PRIMME_SVDS_m :
         *(PRIMME_INT*)value = primme_svds->m;
         break;
      case PRIMME_SVDS_n :
         *(PRIMME_INT*)value = primme_svds->n;
         break;
      case PRIMME_SVDS_matrixMatvec :
         v->matFunc_v = primme_svds->matrixMatvec;
         break;
      case PRIMME_SVDS_matrixMatvec_type:
         *(PRIMME_INT*)value = primme_svds->matrixMatvec_type;
         break;
      case PRIMME_SVDS_applyPreconditioner :
         v->matFunc_v = primme_svds->applyPreconditioner;
         break;
      case PRIMME_SVDS_applyPreconditioner_type:
         *(PRIMME_INT*)value = primme_svds->applyPreconditioner_type;
         break;
      case PRIMME_SVDS_numProcs :
         *(PRIMME_INT*)value = primme_svds->numProcs;
         break;
      case PRIMME_SVDS_procID :
         *(PRIMME_INT*)value = primme_svds->procID;
         break;
      case PRIMME_SVDS_mLocal :
         *(PRIMME_INT*)value = primme_svds->mLocal;
         break;
      case PRIMME_SVDS_nLocal :
         *(PRIMME_INT*)value = primme_svds->nLocal;
         break;
      case PRIMME_SVDS_commInfo :
         *(ptr_v*)value = primme_svds->commInfo;
         break;
      case PRIMME_SVDS_globalSumReal :
         v->globalSumRealFunc_v = primme_svds->globalSumReal;
         break;
      case PRIMME_SVDS_globalSumReal_type:
         *(PRIMME_INT*)value = primme_svds->globalSumReal_type;
         break;
      case PRIMME_SVDS_broadcastReal :
         v->broadcastRealFunc_v = primme_svds->broadcastReal;
         break;
      case PRIMME_SVDS_broadcastReal_type:
         *(PRIMME_INT*)value = primme_svds->broadcastReal_type;
         break;
      case PRIMME_SVDS_internalPrecision:
         *(PRIMME_INT*)value = primme_svds->internalPrecision;
         break;
      case PRIMME_SVDS_numSvals :
         *(PRIMME_INT*)value = primme_svds->numSvals;
         break;
      case PRIMME_SVDS_target :
         *(PRIMME_INT*)value = primme_svds->target;
         break;
      case PRIMME_SVDS_numTargetShifts :
         *(PRIMME_INT*)value = primme_svds->numTargetShifts;
         break;
      case PRIMME_SVDS_targetShifts :
         *(ptr_v*)value = primme_svds->targetShifts;
         break;
      case PRIMME_SVDS_method :
         *(PRIMME_INT*)value = primme_svds->method;
         break;
      case PRIMME_SVDS_methodStage2 :
         *(PRIMME_INT*)value = primme_svds->methodStage2;
         break;
      case PRIMME_SVDS_matrix :
         *(ptr_v*)value = primme_svds->matrix;
         break;
      case PRIMME_SVDS_preconditioner :
         *(ptr_v*)value = primme_svds->preconditioner;
         break;
      case PRIMME_SVDS_locking :
         *(PRIMME_INT*)value = primme_svds->locking;
         break;
      case PRIMME_SVDS_numOrthoConst :
         *(PRIMME_INT*)value = primme_svds->numOrthoConst;
         break;
      case PRIMME_SVDS_aNorm :
         *(double*)value = primme_svds->aNorm;
         break;
      case PRIMME_SVDS_eps :
         *(double*)value = primme_svds->eps;
         break;
      case PRIMME_SVDS_precondition :
         *(PRIMME_INT*)value = primme_svds->precondition;
         break;
      case PRIMME_SVDS_initSize :
         *(PRIMME_INT*)value = primme_svds->initSize;
         break;
      case PRIMME_SVDS_maxBasisSize :
         *(PRIMME_INT*)value = primme_svds->maxBasisSize;
         break;
      case PRIMME_SVDS_maxBlockSize :
         *(PRIMME_INT*)value = primme_svds->maxBlockSize;
         break;
      case PRIMME_SVDS_maxMatvecs :
         *(PRIMME_INT*)value = primme_svds->maxMatvecs;
         break;
      case PRIMME_SVDS_iseed :
         for (i=0; i<4; i++) {
            ((PRIMME_INT*)value)[i] = primme_svds->iseed[i];
         }
         break;
      case PRIMME_SVDS_printLevel :
         *(PRIMME_INT*)value = primme_svds->printLevel;
         break;
      case PRIMME_SVDS_outputFile :
         *(FILE**)value = primme_svds->outputFile;
         break;
      case PRIMME_SVDS_stats_numOuterIterations :
         *(PRIMME_INT*)value = primme_svds->stats.numOuterIterations;
         break;
      case PRIMME_SVDS_stats_numRestarts :
         *(PRIMME_INT*)value = primme_svds->stats.numRestarts;
         break;
      case PRIMME_SVDS_stats_numMatvecs :
         *(PRIMME_INT*)value = primme_svds->stats.numMatvecs;
         break;
      case PRIMME_SVDS_stats_numPreconds :
         *(PRIMME_INT*)value = primme_svds->stats.numPreconds;
         break;
      case PRIMME_SVDS_stats_numGlobalSum:
         *(PRIMME_INT*)value = primme_svds->stats.numGlobalSum;
         break;
      case PRIMME_SVDS_stats_volumeGlobalSum:
         *(PRIMME_INT*)value = primme_svds->stats.volumeGlobalSum;
         break;
      case PRIMME_SVDS_stats_numBroadcast:
         *(PRIMME_INT*)value = primme_svds->stats.numBroadcast;
         break;
      case PRIMME_SVDS_stats_volumeBroadcast:
         *(PRIMME_INT*)value = primme_svds->stats.volumeBroadcast;
         break;
      case PRIMME_SVDS_stats_numOrthoInnerProds:
         *(double*)value = primme_svds->stats.numOrthoInnerProds;
         break;
      case PRIMME_SVDS_stats_elapsedTime :
         *(double*)value = primme_svds->stats.elapsedTime;
         break;
      case PRIMME_SVDS_stats_timeMatvec:
         *(double*)value = primme_svds->stats.timeMatvec;
         break;
      case PRIMME_SVDS_stats_timePrecond:
         *(double*)value = primme_svds->stats.timePrecond;
         break;
      case PRIMME_SVDS_stats_timeOrtho:
         *(double*)value = primme_svds->stats.timeOrtho;
         break;
      case PRIMME_SVDS_stats_timeGlobalSum:
         *(double*)value = primme_svds->stats.timeGlobalSum;
         break;
      case PRIMME_SVDS_stats_timeBroadcast:
         *(double*)value = primme_svds->stats.timeBroadcast;
         break;
      case PRIMME_SVDS_stats_lockingIssue:
         *(PRIMME_INT*)value = primme_svds->stats.lockingIssue;
         break;
      case PRIMME_SVDS_convTestFun:
         v->convTestFun_v = primme_svds->convTestFun;
         break;
      case PRIMME_SVDS_convTestFun_type:
         *(PRIMME_INT*)value = primme_svds->convTestFun_type;
         break;
      case PRIMME_SVDS_convtest:
         *(ptr_v*)value = primme_svds->convtest;
         break;
      case PRIMME_SVDS_monitorFun:
         v->monitorFun_v = primme_svds->monitorFun;
         break;
      case PRIMME_SVDS_monitorFun_type:
         *(PRIMME_INT*)value = primme_svds->monitorFun_type;
         break;
      case PRIMME_SVDS_monitor:
         *(ptr_v*)value = primme_svds->monitor;
         break;
      case PRIMME_SVDS_queue:
         *(ptr_v*)value = primme_svds->queue;
         break;
      case PRIMME_SVDS_profile:
         *(str_v*)value = primme_svds->profile;
         break;
      default:
         return 1;
   }
   return 0;
}

/*******************************************************************************
 * Subroutine primme_svds_set_member - set the value to a parameter in
 *    primme_svds_params
 * 
 * INPUT PARAMETERS
 * ----------------
 * label   reference to the parameter
 * value   value of the parameter
 *
 * INPUT/OUTPUT PARAMETERS
 * -----------------
 * primme_svds  Structure containing various solver parameters and statistics
 *
 * RETURN
 * ------
 * error code  zero if ok
 *
 ******************************************************************************/

int primme_svds_set_member(primme_svds_params *primme_svds,
      primme_svds_params_label label, void *value) {
   int i;
   // Workaround to avoid warnings assigning void* pointers to function pointers
   svds_value_t v = *(svds_value_t*)&value;
   switch(label) {
      case PRIMME_SVDS_primme :
         return 1;
         break;
      case PRIMME_SVDS_primmeStage2 :
         return 1;
         break;
      case PRIMME_SVDS_m :
         primme_svds->m = *(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_n :
         primme_svds->n = *(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_matrixMatvec :
         primme_svds->matrixMatvec = v.matFunc_v;
         break;
      case PRIMME_SVDS_matrixMatvec_type:
         primme_svds->matrixMatvec_type = (primme_op_datatype)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_applyPreconditioner :
         primme_svds->applyPreconditioner = v.matFunc_v;
         break;
      case PRIMME_SVDS_applyPreconditioner_type:
         primme_svds->applyPreconditioner_type = (primme_op_datatype)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_numProcs :
         if (*(PRIMME_INT*)value > INT_MAX) return 1; else
         primme_svds->numProcs = (int)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_procID :
         if (*(PRIMME_INT*)value > INT_MAX) return 1; else 
         primme_svds->procID = (int)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_mLocal :
         primme_svds->mLocal = *(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_nLocal :
         primme_svds->nLocal = *(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_commInfo :
         primme_svds->commInfo = (ptr_v)value;
         break;
      case PRIMME_SVDS_globalSumReal :
         primme_svds->globalSumReal = v.globalSumRealFunc_v;
         break;
      case PRIMME_SVDS_globalSumReal_type:
         primme_svds->globalSumReal_type = (primme_op_datatype)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_broadcastReal :
         primme_svds->broadcastReal = v.broadcastRealFunc_v;
         break;
      case PRIMME_SVDS_internalPrecision:
         primme_svds->internalPrecision = (primme_op_datatype)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_broadcastReal_type:
         primme_svds->broadcastReal_type = (primme_op_datatype)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_numSvals :
         if (*(PRIMME_INT*)value > INT_MAX) return 1; else 
         primme_svds->numSvals = (int)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_target :
         primme_svds->target = (primme_svds_target)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_numTargetShifts :
         if (*(PRIMME_INT*)value > INT_MAX) return 1; else 
         primme_svds->numTargetShifts = (int)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_targetShifts :
         primme_svds->targetShifts = (double*)value;
         break;
      case PRIMME_SVDS_method :
         primme_svds->method = (primme_svds_operator)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_methodStage2 :
         primme_svds->methodStage2 = (primme_svds_operator)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_matrix :
         primme_svds->matrix = (ptr_v)value;
         break;
      case PRIMME_SVDS_preconditioner :
         primme_svds->preconditioner = (ptr_v)value;
         break;
      case PRIMME_SVDS_locking :
         if (*(PRIMME_INT*)value > INT_MAX) return 1; else 
         primme_svds->locking = (int)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_numOrthoConst :
         if (*(PRIMME_INT*)value > INT_MAX) return 1; else 
         primme_svds->numOrthoConst = (int)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_aNorm :
         primme_svds->aNorm = *(double*)value;
         break;
      case PRIMME_SVDS_eps :
         primme_svds->eps = *(double*)value;
         break;
      case PRIMME_SVDS_precondition :
         if (*(PRIMME_INT*)value > INT_MAX) return 1; else 
         primme_svds->precondition = (int)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_initSize :
         if (*(PRIMME_INT*)value > INT_MAX) return 1; else 
         primme_svds->initSize = (int)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_maxBasisSize :
         if (*(PRIMME_INT*)value > INT_MAX) return 1; else 
         primme_svds->maxBasisSize = (int)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_maxBlockSize :
         if (*(PRIMME_INT*)value > INT_MAX) return 1; else 
         primme_svds->maxBlockSize = (int)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_maxMatvecs :
         primme_svds->maxMatvecs = *(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_iseed :
         for (i=0; i<4; i++) {
            primme_svds->iseed[i] = ((PRIMME_INT*)value)[i];
         }
         break;
      case PRIMME_SVDS_printLevel :
         if (*(PRIMME_INT*)value > INT_MAX) return 1; else 
         primme_svds->printLevel = (int)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_outputFile :
         primme_svds->outputFile = (FILE*)value;
         break;
      case PRIMME_SVDS_stats_numOuterIterations :
         primme_svds->stats.numOuterIterations = *(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_stats_numRestarts :
         primme_svds->stats.numRestarts = *(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_stats_numMatvecs :
         primme_svds->stats.numMatvecs = *(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_stats_numPreconds :
         primme_svds->stats.numPreconds = *(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_stats_volumeGlobalSum:
         primme_svds->stats.volumeGlobalSum = *(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_stats_volumeBroadcast:
         primme_svds->stats.volumeBroadcast = *(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_stats_numOrthoInnerProds:
         primme_svds->stats.numOrthoInnerProds = *(double*)value;
         break;
      case PRIMME_SVDS_stats_elapsedTime :
         primme_svds->stats.elapsedTime = *(double*)value;
         break;
      case PRIMME_SVDS_stats_timeMatvec:
         primme_svds->stats.timeMatvec = *(double*)value;
         break;
      case PRIMME_SVDS_stats_timePrecond:
         primme_svds->stats.timePrecond = *(double*)value;
         break;
      case PRIMME_SVDS_stats_timeOrtho:
         primme_svds->stats.timeOrtho = *(double*)value;
         break;
      case PRIMME_SVDS_stats_timeGlobalSum:
         primme_svds->stats.timeGlobalSum = *(double*)value;
         break;
      case PRIMME_SVDS_stats_timeBroadcast:
         primme_svds->stats.timeBroadcast = *(double*)value;
         break;
      case PRIMME_SVDS_convTestFun:
         primme_svds->convTestFun = v.convTestFun_v;
         break;
      case PRIMME_SVDS_convTestFun_type:
         primme_svds->convTestFun_type = (primme_op_datatype)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_convtest:
         primme_svds->convtest = (ptr_v)value;
         break;
      case PRIMME_SVDS_monitorFun:
         primme_svds->monitorFun = v.monitorFun_v;
         break;
      case PRIMME_SVDS_monitorFun_type:
         primme_svds->monitorFun_type = (primme_op_datatype)*(PRIMME_INT*)value;
         break;
      case PRIMME_SVDS_monitor:
         primme_svds->monitor = (ptr_v)value;
         break;
      case PRIMME_SVDS_queue:
         primme_svds->queue = (ptr_v)value;
         break;
      case PRIMME_SVDS_profile:
         primme_svds->profile = (str_v)value;
         break;
      default:
         return 1;
   }
   return 0;
}

/*******************************************************************************
 * Subroutine primme_svds_member_info - return the label value or the label
 *    name, the type and the arity of some member of primme_svds_params.
 * 
 * INPUT/OUTPUT PARAMETERS
 * ----------------
 * label       reference to the parameter by the value in primme_svds_params_label
 * label_name  reference by the string associated to the parameter
 *
 * OUTPUT PARAMETERS
 * -----------------
 * type        kind of value
 * arity       number of elements
 *
 * RETURN
 * ------
 * error code  zero if ok
 *
 ******************************************************************************/

int primme_svds_member_info(primme_svds_params_label *label,
      const char** label_name, primme_type *type, int *arity) {

   /* Quick exit when neither label nor label_name is given */

   if (label == NULL && (label_name == NULL || *label_name == NULL)) {
      return 1;
   }

   /* Get the label from label_name_ and label_name from label_ */
   int set = 0;
#define IF_IS(F)                                                               \
   if ((label_name && *label_name && strcmp(#F, *label_name) == 0) ||          \
         (label && *label == PRIMME_SVDS_##F)) {                               \
      if (label) *label = PRIMME_SVDS_##F;                                     \
      if (label_name) *label_name = #F;                                        \
      set = 1;                                                                 \
   }

   IF_IS(primme);
   IF_IS(primmeStage2);
   IF_IS(m);
   IF_IS(n);
   IF_IS(matrixMatvec);
   IF_IS(matrixMatvec_type);
   IF_IS(applyPreconditioner);
   IF_IS(applyPreconditioner_type);
   IF_IS(numProcs);
   IF_IS(procID);
   IF_IS(mLocal);
   IF_IS(nLocal);
   IF_IS(commInfo);
   IF_IS(globalSumReal);
   IF_IS(globalSumReal_type);
   IF_IS(broadcastReal);
   IF_IS(broadcastReal_type);
   IF_IS(internalPrecision);
   IF_IS(numSvals);
   IF_IS(target);
   IF_IS(numTargetShifts);
   IF_IS(targetShifts);
   IF_IS(method);
   IF_IS(methodStage2);
   IF_IS(matrix);
   IF_IS(preconditioner);
   IF_IS(locking);
   IF_IS(numOrthoConst);
   IF_IS(aNorm);
   IF_IS(eps);
   IF_IS(precondition);
   IF_IS(initSize);
   IF_IS(maxBasisSize);
   IF_IS(maxBlockSize);
   IF_IS(maxMatvecs);
   IF_IS(iseed);
   IF_IS(printLevel);
   IF_IS(outputFile);
   IF_IS(stats_numOuterIterations);
   IF_IS(stats_numRestarts);
   IF_IS(stats_numMatvecs);
   IF_IS(stats_numPreconds);
   IF_IS(stats_numGlobalSum);
   IF_IS(stats_volumeGlobalSum);
   IF_IS(stats_numBroadcast);
   IF_IS(stats_volumeBroadcast);
   IF_IS(stats_numOrthoInnerProds);
   IF_IS(stats_elapsedTime);
   IF_IS(stats_timeMatvec);
   IF_IS(stats_timePrecond);
   IF_IS(stats_timeOrtho);
   IF_IS(stats_timeGlobalSum);
   IF_IS(stats_timeBroadcast);
   IF_IS(stats_lockingIssue);
   IF_IS(convTestFun);
   IF_IS(convTestFun_type);
   IF_IS(convtest);
   IF_IS(monitorFun);
   IF_IS(monitorFun_type);
   IF_IS(monitor);
   IF_IS(queue);
   IF_IS(profile);
#undef IF_IS

   /* Return error if no label was found */
   if (!set) return 1;

   /* Return type and arity */

   switch(*label) {
      /* members with type int */

      case PRIMME_SVDS_matrixMatvec_type: 
      case PRIMME_SVDS_applyPreconditioner_type:
      case PRIMME_SVDS_globalSumReal_type:
      case PRIMME_SVDS_broadcastReal_type:
      case PRIMME_SVDS_internalPrecision:
      case PRIMME_SVDS_m: 
      case PRIMME_SVDS_n:
      case PRIMME_SVDS_numSvals:
      case PRIMME_SVDS_target:
      case PRIMME_SVDS_method:
      case PRIMME_SVDS_methodStage2:
      case PRIMME_SVDS_locking:
      case PRIMME_SVDS_numOrthoConst:
      case PRIMME_SVDS_precondition:
      case PRIMME_SVDS_initSize:
      case PRIMME_SVDS_maxBasisSize:
      case PRIMME_SVDS_maxBlockSize:
      case PRIMME_SVDS_maxMatvecs:
      case PRIMME_SVDS_printLevel:
      case PRIMME_SVDS_stats_numOuterIterations:
      case PRIMME_SVDS_stats_numRestarts:
      case PRIMME_SVDS_stats_numMatvecs:
      case PRIMME_SVDS_stats_numPreconds:
      case PRIMME_SVDS_stats_numGlobalSum:
      case PRIMME_SVDS_stats_volumeGlobalSum:
      case PRIMME_SVDS_stats_numBroadcast:
      case PRIMME_SVDS_stats_volumeBroadcast:
      case PRIMME_SVDS_stats_lockingIssue:
      case PRIMME_SVDS_iseed:
      case PRIMME_SVDS_numProcs: 
      case PRIMME_SVDS_procID: 
      case PRIMME_SVDS_mLocal: 
      case PRIMME_SVDS_nLocal: 
      case PRIMME_SVDS_numTargetShifts:
      case PRIMME_SVDS_convTestFun_type:
      case PRIMME_SVDS_monitorFun_type:
      if (type) *type = primme_int;
      if (arity) *arity = 1;
      break;

      /* members with type double */

      case PRIMME_SVDS_aNorm:
      case PRIMME_SVDS_eps:
      case PRIMME_SVDS_stats_elapsedTime:
      case PRIMME_SVDS_stats_numOrthoInnerProds:
      case PRIMME_SVDS_stats_timeMatvec:
      case PRIMME_SVDS_stats_timePrecond:
      case PRIMME_SVDS_stats_timeOrtho:
      case PRIMME_SVDS_stats_timeGlobalSum:
      case PRIMME_SVDS_stats_timeBroadcast:
      if (type) *type = primme_double;
      if (arity) *arity = 1;
      break;
  
      case PRIMME_SVDS_targetShifts:
      if (type) *type = primme_double;
      if (arity) *arity = 0;
      break;

      /* members with type pointer */

      case PRIMME_SVDS_primme:
      case PRIMME_SVDS_primmeStage2:
      case PRIMME_SVDS_matrixMatvec: 
      case PRIMME_SVDS_applyPreconditioner:
      case PRIMME_SVDS_commInfo:
      case PRIMME_SVDS_globalSumReal:
      case PRIMME_SVDS_broadcastReal:
      case PRIMME_SVDS_matrix:
      case PRIMME_SVDS_preconditioner:
      case PRIMME_SVDS_outputFile:
      case PRIMME_SVDS_convTestFun:
      case PRIMME_SVDS_convtest:
      case PRIMME_SVDS_monitorFun:
      case PRIMME_SVDS_monitor:
      case PRIMME_SVDS_queue:
      if (type) *type = primme_pointer;
      if (arity) *arity = 1;
      break;

      case PRIMME_SVDS_profile:
      if (type) *type = primme_string;
      if (arity) *arity = 1;
      break;

      default: 
      return 1;
   }

   return 0;
}
 
/*******************************************************************************
 * Subroutine primme_svds_constant_info - return the value of a primme svds enum
 *    constant.
 * 
 * INPUT PARAMETERS
 * ----------------
 * label_name  name of the constant
 *
 * OUTPUT PARAMETERS
 * -----------------
 * value       numerical value of the constant
 *
 * RETURN
 * ------
 * error code  zero if ok
 *
 ******************************************************************************/

int primme_svds_constant_info(const char* label_name, int *value) {

#define IF_IS(F) if (strcmp(#F, label_name) == 0) {*value = (int)F; return 0;}

   /* method in primme_svds_set_method */
   
   IF_IS(primme_svds_default);
   IF_IS(primme_svds_hybrid);
   IF_IS(primme_svds_normalequations);
   IF_IS(primme_svds_augmented);
   
   /* enum members for targeting and operator */
   
   IF_IS(primme_svds_largest);
   IF_IS(primme_svds_smallest);
   IF_IS(primme_svds_closest_abs);
   IF_IS(primme_svds_op_none);
   IF_IS(primme_svds_op_AtA);
   IF_IS(primme_svds_op_AAt);
   IF_IS(primme_svds_op_augmented);
#undef IF_IS

   /* try primme constants */

   return primme_constant_info(label_name, value);
}

/*******************************************************************************
 * Subroutine primme_svds_enum_member_info - return the value of a string
 * representing an enum constant, or vice versa.
 *
 * INPUT/OUTPUT PARAMETERS
 * -----------------------
 * label       member to which the constant relates
 * value       (in) if *value >= 0, value to get the associated string,
 *             (in/out) if *value < 0, return the value associated to
 *             value_name.
 * value_name  (in) if *value_name > 0, string for which to seek the value
 *             (in/out) if *value_name == 0, return the associated to value.
 *
 * RETURN
 * ------
 * error code   0: OK
 *             -1: Invalid input
 *             -2: either value or value_name was not found
 *
 ******************************************************************************/

int primme_svds_enum_member_info(
      primme_svds_params_label label, int *value, const char **value_name) {

   if (!value || !value_name || (*value >= 0 && *value_name) ||
         (*value < 0 && !*value_name)) {
      return -1;
   }

#define IF_IS(F)                                                               \
   if (*value == (int)(F) || (*value_name && strcmp(#F, *value_name) == 0)) {  \
      *value = (int)(F);                                                       \
      *value_name = #F;                                                        \
      return 0;                                                                \
   }

   switch(label) {
   // Hack: Check method
   case PRIMME_SVDS_commInfo:
   IF_IS(primme_svds_default);
   IF_IS(primme_svds_hybrid);
   IF_IS(primme_svds_normalequations);
   IF_IS(primme_svds_augmented);
   break;
   
   case PRIMME_SVDS_target:
   IF_IS(primme_svds_largest);
   IF_IS(primme_svds_smallest);
   IF_IS(primme_svds_closest_abs);
   break;

   case PRIMME_SVDS_method:
   case PRIMME_SVDS_methodStage2:
   IF_IS(primme_svds_op_none);
   IF_IS(primme_svds_op_AtA);
   IF_IS(primme_svds_op_AAt);
   IF_IS(primme_svds_op_augmented);
   break;

   default: break;
   }
#undef IF_IS

   /* return error if label not found */

   return -2;
}
 
#endif /* USE_DOUBLE */
