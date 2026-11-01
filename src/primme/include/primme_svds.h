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
 * File: primme_svds.h
 *
 * Purpose - Main header with the PRIMME SVDS C interface functions.
 *
 ******************************************************************************/


#ifndef PRIMME_SVDS_H
#define PRIMME_SVDS_H

#include "primme_eigs.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
   primme_svds_largest,
   primme_svds_smallest,
   primme_svds_closest_abs
} primme_svds_target;

typedef enum {
   primme_svds_default,
   primme_svds_hybrid,
   primme_svds_normalequations, /* At*A or A*At */
   primme_svds_augmented
} primme_svds_preset_method;

typedef enum {
   primme_svds_op_none,
   primme_svds_op_AtA,
   primme_svds_op_AAt,
   primme_svds_op_augmented
} primme_svds_operator;

typedef struct primme_svds_stats {
   PRIMME_INT numOuterIterations;
   PRIMME_INT numRestarts;
   PRIMME_INT numMatvecs;
   PRIMME_INT numPreconds;
   PRIMME_INT numGlobalSum;         /* times called globalSumReal */
   PRIMME_INT numBroadcast;         /* times called broadcastReal */
   PRIMME_INT volumeGlobalSum;      /* number of SCALARs reduced by globalSumReal */
   PRIMME_INT volumeBroadcast;      /* number of SCALARs broadcast by broadcastReal */
   double numOrthoInnerProds;       /* number of inner prods done by Ortho */
   double elapsedTime; 
   double timeMatvec;               /* time expend by matrixMatvec */
   double timePrecond;              /* time expend by applyPreconditioner */
   double timeOrtho;                /* time expend by ortho  */
   double timeGlobalSum;            /* time expend by globalSumReal  */
   double timeBroadcast;            /* time expend by broadcastReal  */
   PRIMME_INT lockingIssue;         /* Some converged with a weak criterion */
} primme_svds_stats;

typedef struct primme_svds_params {
   /**** Low interface: configuration for the eigensolver */
   primme_params primme;  /* Keep it as first field to access primme_svds_params from
                             primme_params */
   primme_params primmeStage2; /* other primme_params, used by hybrid */

   /* Specify the size of the rectangular matrix A */
   PRIMME_INT m; /* number of rows */ 
   PRIMME_INT n; /* number of columns */

   /***** High interface: these values are transferred to primme and primmeStage2 properly */
   void (*matrixMatvec) 
      (void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
       int *transpose, struct primme_svds_params *primme_svds, int *ierr);
   primme_op_datatype matrixMatvec_type;
   void (*applyPreconditioner)
      (void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
       int *transpose, struct primme_svds_params *primme_svds, int *ierr);
   primme_op_datatype applyPreconditioner_type;

   /* Input for the following is only required for parallel programs */
   int numProcs;
   int procID;
   PRIMME_INT mLocal;
   PRIMME_INT nLocal;
   void *commInfo;
   void (*globalSumReal)
      (void *sendBuf, void *recvBuf, int *count,
       struct primme_svds_params *primme_svds, int *ierr);
   primme_op_datatype globalSumReal_type;
   void (*broadcastReal)(void *buffer, int *count,
         struct primme_svds_params *primme_svds, int *ierr);
   primme_op_datatype broadcastReal_type;

   /* Though primme_svds_initialize will assign defaults, most users will set these */
   int numSvals;
   primme_svds_target target;
   int numTargetShifts;    /* For primme_svds_augmented method, user has to */ 
   double *targetShifts;   /* make sure  at least one shift must also be set */
   primme_svds_operator method; /* one of primme_svds_AtA, primme_svds_AAt or primme_svds_augmented */
   primme_svds_operator methodStage2; /* hybrid second stage method; accepts the same values as method */

   /* These pointers may be used for users to provide matrix/preconditioner */
   void *matrix;
   void *preconditioner;

   /* The following will be given default values depending on the method */
   int locking;
   int numOrthoConst;
   double aNorm;
   double eps;

   int precondition;
   int initSize;
   int maxBasisSize;
   int maxBlockSize;
   PRIMME_INT maxMatvecs;
   PRIMME_INT iseed[4];
   int printLevel;
   primme_op_datatype internalPrecision; /* force primme to work in that precision */
   FILE *outputFile;
   struct primme_svds_stats stats;

   void (*convTestFun)(double *sval, void *leftsvec, void *rightsvec,
         double *rNorm, int *method, int *isconv,
         struct primme_svds_params *primme, int *ierr);
   primme_op_datatype convTestFun_type; /* expected type of evec */
   void *convtest;
   void (*monitorFun)(void *basisSvals, int *basisSize, int *basisFlags,
         int *iblock, int *blockSize, void *basisNorms, int *numConverged,
         void *lockedSvals, int *numLocked, int *lockedFlags, void *lockedNorms,
         int *inner_its, void *LSRes, const char *msg, double *time,
         primme_event *event, int *stage,
         struct primme_svds_params *primme_svds, int *err);
   primme_op_datatype monitorFun_type; /* expected type of float-point arrays */
   void *monitor;
   void *queue;   	/* magma device queue (magma_queue_t*) */
   const char *profile; /* regex expression with functions to monitor times */
} primme_svds_params;

typedef enum {
   /* NOTE: you can maintain the column of numbers with g+Ctrl-A in vim */
   PRIMME_SVDS_invalid_label                = 0,
   PRIMME_SVDS_primme                       = 1,
   PRIMME_SVDS_primmeStage2                 = 2,
   PRIMME_SVDS_m                            = 3,
   PRIMME_SVDS_n                            = 4,
   PRIMME_SVDS_matrixMatvec                 = 5,
   PRIMME_SVDS_matrixMatvec_type            = 6,
   PRIMME_SVDS_applyPreconditioner          = 7,
   PRIMME_SVDS_applyPreconditioner_type     = 8,
   PRIMME_SVDS_numProcs                     = 9,
   PRIMME_SVDS_procID                       = 10,
   PRIMME_SVDS_mLocal                       = 11,
   PRIMME_SVDS_nLocal                       = 12,
   PRIMME_SVDS_commInfo                     = 13,
   PRIMME_SVDS_globalSumReal                = 14,
   PRIMME_SVDS_globalSumReal_type           = 15,
   PRIMME_SVDS_broadcastReal                = 16,
   PRIMME_SVDS_broadcastReal_type           = 17,
   PRIMME_SVDS_numSvals                     = 18,
   PRIMME_SVDS_target                       = 19,
   PRIMME_SVDS_numTargetShifts              = 20,
   PRIMME_SVDS_targetShifts                 = 21,
   PRIMME_SVDS_method                       = 22,
   PRIMME_SVDS_methodStage2                 = 23,
   PRIMME_SVDS_matrix                       = 24,
   PRIMME_SVDS_preconditioner               = 25,
   PRIMME_SVDS_locking                      = 26,
   PRIMME_SVDS_numOrthoConst                = 27,
   PRIMME_SVDS_aNorm                        = 28,
   PRIMME_SVDS_eps                          = 29,
   PRIMME_SVDS_precondition                 = 30,
   PRIMME_SVDS_initSize                     = 31,
   PRIMME_SVDS_maxBasisSize                 = 32,
   PRIMME_SVDS_maxBlockSize                 = 33,
   PRIMME_SVDS_maxMatvecs                   = 34,
   PRIMME_SVDS_iseed                        = 35,
   PRIMME_SVDS_printLevel                   = 36,
   PRIMME_SVDS_internalPrecision            = 37,
   PRIMME_SVDS_outputFile                   = 38,
   PRIMME_SVDS_stats_numOuterIterations     = 39,
   PRIMME_SVDS_stats_numRestarts            = 40,
   PRIMME_SVDS_stats_numMatvecs             = 41,
   PRIMME_SVDS_stats_numPreconds            = 42,
   PRIMME_SVDS_stats_numGlobalSum           = 43,
   PRIMME_SVDS_stats_volumeGlobalSum        = 44,
   PRIMME_SVDS_stats_numBroadcast           = 45,
   PRIMME_SVDS_stats_volumeBroadcast        = 46,
   PRIMME_SVDS_stats_numOrthoInnerProds     = 47,
   PRIMME_SVDS_stats_elapsedTime            = 48,
   PRIMME_SVDS_stats_timeMatvec             = 49,
   PRIMME_SVDS_stats_timePrecond            = 50,
   PRIMME_SVDS_stats_timeOrtho              = 51,
   PRIMME_SVDS_stats_timeGlobalSum          = 52,
   PRIMME_SVDS_stats_timeBroadcast          = 53,
   PRIMME_SVDS_stats_lockingIssue           = 54,
   PRIMME_SVDS_convTestFun                  = 55,
   PRIMME_SVDS_convTestFun_type             = 56,
   PRIMME_SVDS_convtest                     = 57,
   PRIMME_SVDS_monitorFun                   = 58,
   PRIMME_SVDS_monitorFun_type              = 59,
   PRIMME_SVDS_monitor                      = 60,
   PRIMME_SVDS_queue                        = 61,
   PRIMME_SVDS_profile                      = 62 
} primme_svds_params_label;

int hprimme_svds(PRIMME_HALF *svals, PRIMME_HALF *svecs, PRIMME_HALF *resNorms,
      primme_svds_params *primme_svds);
int kprimme_svds(PRIMME_HALF *svals, PRIMME_COMPLEX_HALF *svecs, PRIMME_HALF *resNorms,
      primme_svds_params *primme_svds);
int sprimme_svds(float *svals, float *svecs, float *resNorms,
      primme_svds_params *primme_svds);
int cprimme_svds(float *svals, PRIMME_COMPLEX_FLOAT *svecs, float *resNorms,
      primme_svds_params *primme_svds);
int dprimme_svds(double *svals, double *svecs, double *resNorms,
      primme_svds_params *primme_svds);
int zprimme_svds(double *svals, PRIMME_COMPLEX_DOUBLE *svecs, double *resNorms,
      primme_svds_params *primme_svds);
int magma_hprimme_svds(PRIMME_HALF *svals, PRIMME_HALF *svecs, PRIMME_HALF *resNorms,
      primme_svds_params *primme_svds);
int magma_kprimme_svds(PRIMME_HALF *svals, PRIMME_COMPLEX_HALF *svecs, PRIMME_HALF *resNorms,
      primme_svds_params *primme_svds);
int magma_sprimme_svds(float *svals, float *svecs, float *resNorms,
      primme_svds_params *primme_svds);
int magma_cprimme_svds(float *svals, PRIMME_COMPLEX_FLOAT *svecs, float *resNorms,
      primme_svds_params *primme_svds);
int magma_dprimme_svds(double *svals, double *svecs, double *resNorms,
      primme_svds_params *primme_svds);
int magma_zprimme_svds(double *svals, PRIMME_COMPLEX_DOUBLE *svecs, double *resNorms,
      primme_svds_params *primme_svds);
int cublas_hprimme_svds(PRIMME_HALF *svals, PRIMME_HALF *svecs, PRIMME_HALF *resNorms,
      primme_svds_params *primme_svds);
int cublas_kprimme_svds(PRIMME_HALF *svals, PRIMME_COMPLEX_HALF *svecs, PRIMME_HALF *resNorms,
      primme_svds_params *primme_svds);
int cublas_sprimme_svds(float *svals, float *svecs, float *resNorms,
      primme_svds_params *primme_svds);
int cublas_cprimme_svds(float *svals, PRIMME_COMPLEX_FLOAT *svecs, float *resNorms,
      primme_svds_params *primme_svds);
int cublas_dprimme_svds(double *svals, double *svecs, double *resNorms,
      primme_svds_params *primme_svds);
int cublas_zprimme_svds(double *svals, PRIMME_COMPLEX_DOUBLE *svecs, double *resNorms,
      primme_svds_params *primme_svds);


int hsprimme_svds(float *svals, PRIMME_HALF *svecs, float *resNorms,
      primme_svds_params *primme_svds);
int ksprimme_svds(float *svals, PRIMME_COMPLEX_HALF *svecs, float *resNorms,
      primme_svds_params *primme_svds);
int magma_hsprimme_svds(float *svals, PRIMME_HALF *svecs, float *resNorms,
      primme_svds_params *primme_svds);
int magma_ksprimme_svds(float *svals, PRIMME_COMPLEX_HALF *svecs, float *resNorms,
      primme_svds_params *primme_svds);
int cublas_hsprimme_svds(float *svals, PRIMME_HALF *svecs, float *resNorms,
      primme_svds_params *primme_svds);
int cublas_ksprimme_svds(float *svals, PRIMME_COMPLEX_HALF *svecs, float *resNorms,
      primme_svds_params *primme_svds);



primme_svds_params * primme_svds_params_create(void);
int primme_svds_params_destroy(primme_svds_params *primme_svds);
void primme_svds_initialize(primme_svds_params *primme_svds);
int primme_svds_set_method(primme_svds_preset_method method,
      primme_preset_method methodStage1, primme_preset_method methodStage2,
      primme_svds_params *primme_svds);
void primme_svds_display_params(primme_svds_params primme_svds);
void primme_svds_free(primme_svds_params *primme_svds);
int primme_svds_get_member(primme_svds_params *primme_svds,
      primme_svds_params_label label, void *value);
int primme_svds_set_member(primme_svds_params *primme_svds,
      primme_svds_params_label label, void *value);
int primme_svds_member_info(primme_svds_params_label *label,
      const char** label_name, primme_type *type, int *arity);
int primme_svds_constant_info(const char* label_name, int *value);
int primme_svds_enum_member_info(
      primme_svds_params_label label, int *value, const char **value_name);

#ifdef __cplusplus
}
#endif

#endif /* PRIMME_SVDS_H */  
