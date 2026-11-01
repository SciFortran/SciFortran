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
 * File: primme_eigs.h
 * 
 * Purpose - Main header with the PRIMME EIGS C interface functions.
 * 
 ******************************************************************************/

#ifndef PRIMME_EIGS_H
#define PRIMME_EIGS_H

#include <stdio.h>

#include "primme.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
   primme_smallest,        /* leftmost eigenvalues */
   primme_largest,         /* rightmost eigenvalues */
   primme_closest_geq,     /* leftmost but greater than the target shift */
   primme_closest_leq,     /* rightmost but less than the target shift */
   primme_closest_abs,     /* the closest to the target shift */
   primme_largest_abs      /* the farthest to the target shift */
} primme_target;

/* projection methods for extraction */
typedef enum {
   primme_proj_default,
   primme_proj_RR,          /* Rayleigh-Ritz */
   primme_proj_harmonic,    /* Harmonic Rayleigh-Ritz */
   primme_proj_refined      /* refined with fixed target */
} primme_projection;

typedef enum {         /* Initially fill up the search subspace with: */
   primme_init_default,
   primme_init_krylov, /* a) Krylov with the last vector provided by the user or random */
   primme_init_random, /* b) just random vectors */
   primme_init_user    /* c) provided vectors or a single random vector */
} primme_init;

typedef enum {
   primme_full_LTolerance,
   primme_decreasing_LTolerance,
   primme_adaptive_ETolerance,
   primme_adaptive
} primme_convergencetest;


/* Identifies the type of event for which monitor is being called */
typedef enum {
   primme_event_outer_iteration,    /* report at every outer iteration        */
   primme_event_inner_iteration,    /* report at every QMR iteration          */
   primme_event_restart,            /* report at every basis restart          */
   primme_event_reset,              /* event launch if basis reset            */
   primme_event_converged,          /* report new pair marked as converged    */
   primme_event_locked,             /* report new pair marked as locked       */
   primme_event_message,            /* report warning                         */
   primme_event_profile             /* report time from consumed by a function*/          
} primme_event;

/* Orthogonalization variant */
typedef enum {
   primme_orth_default,
   primme_orth_implicit_I,          /* assume for search subspace V, V'*B*V = I */
   primme_orth_explicit_I           /* explicitly compute V'*B*V */
} primme_orth;

/* Datatype of vectors passed on matrixMatvec, applyPreconditioner,           */
/* massMatrixMatvec and globalSumReal                                         */
typedef enum {
   primme_op_default,   /* The same type as the primme call */
   primme_op_half,
   primme_op_float,
   primme_op_double,
   primme_op_quad,
   primme_op_int
} primme_op_datatype;

typedef struct primme_stats {
   PRIMME_INT numOuterIterations;
   PRIMME_INT numRestarts;
   PRIMME_INT numMatvecs;
   PRIMME_INT numPreconds;
   PRIMME_INT numGlobalSum;         /* times called globalSumReal */
   PRIMME_INT numBroadcast;         /* times called broadcastReal */
   PRIMME_INT volumeGlobalSum;      /* number of SCALARs reduced by globalSumReal */
   PRIMME_INT volumeBroadcast;      /* number of SCALARs broadcast by broadcastReal */
   double flopsDense;               /* FLOPS done by Num_update_VWXR_Sprimme */
   double numOrthoInnerProds;       /* number of inner prods done by Ortho */
   double elapsedTime; 
   double timeMatvec;               /* time expend by matrixMatvec */
   double timePrecond;              /* time expend by applyPreconditioner */
   double timeOrtho;                /* time expend by ortho  */
   double timeGlobalSum;            /* time expend by globalSumReal  */
   double timeBroadcast;            /* time expend by broadcastReal  */
   double timeDense;                /* time expend by Num_update_VWXR_Sprimme */
   double estimateMinEVal;          /* the leftmost Ritz value seen */
   double estimateMaxEVal;          /* the rightmost Ritz value seen */
   double estimateLargestSVal;      /* absolute value of the farthest to zero Ritz value seen */
   double estimateBNorm;            /* estimation of norm of B */
   double estimateInvBNorm;         /* estimation of norm of inv(B) */
   double maxConvTol;               /* largest norm residual of a locked eigenpair */
   double estimateResidualError;    /* accumulated error in V and W */
   PRIMME_INT lockingIssue;         /* Some converged with a weak criterion */
} primme_stats;

typedef struct JD_projectors {
   int LeftQ;
   int LeftX;
   int RightQ;
   int RightX;
   int SkewQ;
   int SkewX;
} JD_projectors;

typedef struct projection_params {
   primme_projection projection;
} projection_params;

typedef struct correction_params {
   int precondition;
   int robustShifts;
   int maxInnerIterations;
   struct JD_projectors projectors;
   primme_convergencetest convTest;
   double relTolBase;
} correction_params;


typedef struct restarting_params {
   int maxPrevRetain;
} restarting_params;


/*--------------------------------------------------------------------------*/
typedef struct primme_params {

   /* The user must input at least the following two arguments */
   PRIMME_INT n;
   void (*matrixMatvec)
      ( void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
        struct primme_params *primme, int *ierr);
   primme_op_datatype matrixMatvec_type; /* expected type of x and y */

   /* Preconditioner applied on block of vectors (if available) */
   void (*applyPreconditioner)
      ( void *x, PRIMME_INT *ldx,  void *y, PRIMME_INT *ldy, int *blockSize,
        struct primme_params *primme, int *ierr);
   primme_op_datatype applyPreconditioner_type; /* expected type of x and y */

   /* Matrix times a multivector for mass matrix B for generalized Ax = xBl */
   void (*massMatrixMatvec)
      ( void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
        struct primme_params *primme, int *ierr);
   primme_op_datatype massMatrixMatvec_type; /* expected type of x and y */

   /* input for the following is only required for parallel programs */
   int numProcs;
   int procID;
   PRIMME_INT nLocal;
   void *commInfo;
   void (*globalSumReal)
      (void *sendBuf, void *recvBuf, int *count, struct primme_params *primme,
       int *ierr );
   primme_op_datatype globalSumReal_type; /* expected type of sendBuf and recvBuf */
   void (*broadcastReal)(
         void *buffer, int *count, struct primme_params *primme, int *ierr);
   primme_op_datatype broadcastReal_type; /* expected type of buffer */

   /*Though primme_initialize will assign defaults, most users will set these */
   int numEvals;          
   primme_target target; 
   int numTargetShifts;              /* For targeting interior epairs,      */
   double *targetShifts;             /* at least one shift must also be set */

   /* the following will be given default values depending on the method */
   int dynamicMethodSwitch;
   int locking;
   int initSize;
   int numOrthoConst;
   int maxBasisSize;
   int minRestartSize;
   int maxBlockSize;
   PRIMME_INT maxMatvecs;
   PRIMME_INT maxOuterIterations;
   PRIMME_INT iseed[4];
   double aNorm;
   double BNorm;                 /* Approximate 2-norm of B */
   double invBNorm;              /* Approximate 2-norm of inv(B) */
   double eps;
   primme_orth orth;
   primme_op_datatype internalPrecision; /* force primme to work in that precision */

   int printLevel;
   FILE *outputFile;

   void *matrix;
   void *preconditioner;
   void *massMatrix;
   double *ShiftsForPreconditioner;
   primme_init initBasisMode;
   PRIMME_INT ldevecs;
   PRIMME_INT ldOPs;

   struct projection_params projectionParams; 
   struct restarting_params restartingParams;
   struct correction_params correctionParams;
   struct primme_stats stats;

   void (*convTestFun)(double *eval, void *evec, double *rNorm, int *isconv, 
         struct primme_params *primme, int *ierr);
   primme_op_datatype convTestFun_type; /* expected type of evec */
   void *convtest;
   void (*monitorFun)(void *basisEvals, int *basisSize, int *basisFlags,
         int *iblock, int *blockSize, void *basisNorms, int *numConverged,
         void *lockedEvals, int *numLocked, int *lockedFlags, void *lockedNorms,
         int *inner_its, void *LSRes, const char *msg, double *time,
         primme_event *event, struct primme_params *primme, int *err);
   primme_op_datatype monitorFun_type; /* expected type of float-point arrays */
   void *monitor;
   void *queue;      /* magma device queue (magma_queue_t*) */
   const char *profile; /* regex expression with functions to monitor times */
} primme_params;
/*---------------------------------------------------------------------------*/

typedef enum {
   PRIMME_DEFAULT_METHOD,
   PRIMME_DYNAMIC,
   PRIMME_DEFAULT_MIN_TIME,
   PRIMME_DEFAULT_MIN_MATVECS,
   PRIMME_Arnoldi,
   PRIMME_GD,
   PRIMME_GD_plusK,
   PRIMME_GD_Olsen_plusK,
   PRIMME_JD_Olsen_plusK,
   PRIMME_RQI,
   PRIMME_JDQR,
   PRIMME_JDQMR,
   PRIMME_JDQMR_ETol,
   PRIMME_STEEPEST_DESCENT,
   PRIMME_LOBPCG_OrthoBasis,
   PRIMME_LOBPCG_OrthoBasis_Window
} primme_preset_method;


/* Indicates the type of primme_param's members. Used in primme_member_info */
/* and in primme_svds_member_info                                           */

typedef enum {
   primme_int,
   primme_double,
   primme_pointer,
   primme_string
} primme_type;

typedef enum {
   /* NOTE: you can maintain the column of numbers with g+Ctrl-A in vim */
   PRIMME_invalid_label                          = 0  ,
   PRIMME_n                                      = 1  ,
   PRIMME_matrixMatvec                           = 2  ,
   PRIMME_matrixMatvec_type                      = 3  ,
   PRIMME_applyPreconditioner                    = 4  ,
   PRIMME_applyPreconditioner_type               = 5  ,
   PRIMME_massMatrixMatvec                       = 6  ,
   PRIMME_massMatrixMatvec_type                  = 7  ,
   PRIMME_numProcs                               = 8  ,
   PRIMME_procID                                 = 9  ,
   PRIMME_commInfo                               = 10  ,
   PRIMME_nLocal                                 = 11  ,
   PRIMME_globalSumReal                          = 12  ,
   PRIMME_globalSumReal_type                     = 13  ,
   PRIMME_broadcastReal                          = 14  ,
   PRIMME_broadcastReal_type                     = 15  ,
   PRIMME_numEvals                               = 16  ,
   PRIMME_target                                 = 17  ,
   PRIMME_numTargetShifts                        = 18  ,
   PRIMME_targetShifts                           = 19  ,
   PRIMME_locking                                = 20  ,
   PRIMME_initSize                               = 21  ,
   PRIMME_numOrthoConst                          = 22  ,
   PRIMME_maxBasisSize                           = 23  ,
   PRIMME_minRestartSize                         = 24  ,
   PRIMME_maxBlockSize                           = 25  ,
   PRIMME_maxMatvecs                             = 26  ,
   PRIMME_maxOuterIterations                     = 27  ,
   PRIMME_iseed                                  = 28  ,
   PRIMME_aNorm                                  = 29  ,
   PRIMME_BNorm                                  = 30  ,
   PRIMME_invBNorm                               = 31  ,
   PRIMME_eps                                    = 32  ,
   PRIMME_orth                                   = 33  ,
   PRIMME_internalPrecision                      = 34  ,
   PRIMME_printLevel                             = 35  ,
   PRIMME_outputFile                             = 36  ,
   PRIMME_matrix                                 = 37  ,
   PRIMME_massMatrix                             = 38  ,
   PRIMME_preconditioner                         = 39  ,
   PRIMME_ShiftsForPreconditioner                = 40  ,
   PRIMME_initBasisMode                          = 41  ,
   PRIMME_projectionParams_projection            = 42  ,
   PRIMME_restartingParams_maxPrevRetain         = 43  ,
   PRIMME_correctionParams_precondition          = 44  ,
   PRIMME_correctionParams_robustShifts          = 45  ,
   PRIMME_correctionParams_maxInnerIterations    = 46  ,
   PRIMME_correctionParams_projectors_LeftQ      = 47  ,
   PRIMME_correctionParams_projectors_LeftX      = 48  ,
   PRIMME_correctionParams_projectors_RightQ     = 49  ,
   PRIMME_correctionParams_projectors_RightX     = 50  ,
   PRIMME_correctionParams_projectors_SkewQ      = 51  ,
   PRIMME_correctionParams_projectors_SkewX      = 52  ,
   PRIMME_correctionParams_convTest              = 53  ,
   PRIMME_correctionParams_relTolBase            = 54  ,
   PRIMME_stats_numOuterIterations               = 55  ,
   PRIMME_stats_numRestarts                      = 56  ,
   PRIMME_stats_numMatvecs                       = 57  ,
   PRIMME_stats_numPreconds                      = 58  ,
   PRIMME_stats_numGlobalSum                     = 59  ,
   PRIMME_stats_volumeGlobalSum                  = 60  ,
   PRIMME_stats_numBroadcast                     = 61  ,
   PRIMME_stats_volumeBroadcast                  = 62  ,
   PRIMME_stats_flopsDense                       = 63  ,
   PRIMME_stats_numOrthoInnerProds               = 64  ,
   PRIMME_stats_elapsedTime                      = 65  ,
   PRIMME_stats_timeMatvec                       = 66  ,
   PRIMME_stats_timePrecond                      = 67  ,
   PRIMME_stats_timeOrtho                        = 68  ,
   PRIMME_stats_timeGlobalSum                    = 69  ,
   PRIMME_stats_timeBroadcast                    = 70  ,
   PRIMME_stats_timeDense                        = 71  ,
   PRIMME_stats_estimateMinEVal                  = 72  ,
   PRIMME_stats_estimateMaxEVal                  = 73  ,
   PRIMME_stats_estimateLargestSVal              = 74  ,
   PRIMME_stats_estimateBNorm                    = 75  ,
   PRIMME_stats_estimateInvBNorm                 = 76  ,
   PRIMME_stats_maxConvTol                       = 77  ,
   PRIMME_stats_lockingIssue                     = 78  ,
   PRIMME_dynamicMethodSwitch                    = 79  ,
   PRIMME_convTestFun                            = 80  ,
   PRIMME_convTestFun_type                       = 81  ,
   PRIMME_convtest                               = 82  ,
   PRIMME_ldevecs                                = 83  ,
   PRIMME_ldOPs                                  = 84  ,
   PRIMME_monitorFun                             = 85  ,
   PRIMME_monitorFun_type                        = 86  ,
   PRIMME_monitor                                = 87  ,
   PRIMME_queue                                  = 88  ,
   PRIMME_profile                                = 89  
} primme_params_label;

/* Hermitian operator */

int hprimme(PRIMME_HALF *evals, PRIMME_HALF *evecs, PRIMME_HALF *resNorms, 
      primme_params *primme);
int kprimme(PRIMME_HALF *evals, PRIMME_COMPLEX_HALF *evecs, PRIMME_HALF *resNorms, 
      primme_params *primme);
int sprimme(float *evals, float *evecs, float *resNorms, 
      primme_params *primme);
int cprimme(float *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, 
      primme_params *primme);
int dprimme(double *evals, double *evecs, double *resNorms, 
      primme_params *primme);
int zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, 
      primme_params *primme);
int magma_hprimme(PRIMME_HALF *evals, PRIMME_HALF *evecs, PRIMME_HALF *resNorms, 
      primme_params *primme);
int magma_kprimme(PRIMME_HALF *evals, PRIMME_COMPLEX_HALF *evecs, PRIMME_HALF *resNorms, 
      primme_params *primme);
int magma_sprimme(float *evals, float *evecs, float *resNorms, 
      primme_params *primme);
int magma_cprimme(float *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, 
      primme_params *primme);
int magma_dprimme(double *evals, double *evecs, double *resNorms, 
      primme_params *primme);
int magma_zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, 
      primme_params *primme);
int cublas_hprimme(PRIMME_HALF *evals, PRIMME_HALF *evecs, PRIMME_HALF *resNorms, 
      primme_params *primme);
int cublas_kprimme(PRIMME_HALF *evals, PRIMME_COMPLEX_HALF *evecs, PRIMME_HALF *resNorms, 
      primme_params *primme);
int cublas_sprimme(float *evals, float *evecs, float *resNorms, 
      primme_params *primme);
int cublas_cprimme(float *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, 
      primme_params *primme);
int cublas_dprimme(double *evals, double *evecs, double *resNorms, 
      primme_params *primme);
int cublas_zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, 
      primme_params *primme);


int hsprimme(float *evals, PRIMME_HALF *evecs, float *resNorms, 
      primme_params *primme);
int ksprimme(float *evals, PRIMME_COMPLEX_HALF *evecs, float *resNorms, 
      primme_params *primme);
int magma_hsprimme(float *evals, PRIMME_HALF *evecs, float *resNorms, 
      primme_params *primme);
int magma_ksprimme(float *evals, PRIMME_COMPLEX_HALF *evecs, float *resNorms, 
      primme_params *primme);
int cublas_hsprimme(float *evals, PRIMME_HALF *evecs, float *resNorms, 
      primme_params *primme);
int cublas_ksprimme(float *evals, PRIMME_COMPLEX_HALF *evecs, float *resNorms, 
      primme_params *primme);


/* Normal operator */

int kprimme_normal(PRIMME_COMPLEX_HALF *evals, PRIMME_COMPLEX_HALF *evecs, PRIMME_HALF *resNorms, 
      primme_params *primme);
int cprimme_normal(PRIMME_COMPLEX_FLOAT *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, 
      primme_params *primme);
int zprimme_normal(PRIMME_COMPLEX_DOUBLE *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, 
      primme_params *primme);
int magma_kprimme_normal(PRIMME_COMPLEX_HALF *evals, PRIMME_COMPLEX_HALF *evecs, PRIMME_HALF *resNorms, 
      primme_params *primme);
int magma_cprimme_normal(PRIMME_COMPLEX_FLOAT *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, 
      primme_params *primme);
int magma_zprimme_normal(PRIMME_COMPLEX_DOUBLE *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, 
      primme_params *primme);
int cublas_kprimme_normal(PRIMME_COMPLEX_HALF *evals, PRIMME_COMPLEX_HALF *evecs, PRIMME_HALF *resNorms, 
      primme_params *primme);
int cublas_cprimme_normal(PRIMME_COMPLEX_FLOAT *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, 
      primme_params *primme);
int cublas_zprimme_normal(PRIMME_COMPLEX_DOUBLE *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, 
      primme_params *primme);


int kcprimme_normal(PRIMME_COMPLEX_FLOAT *evals, PRIMME_COMPLEX_HALF *evecs, float *resNorms, 
      primme_params *primme);
int magma_kcprimme_normal(PRIMME_COMPLEX_FLOAT *evals, PRIMME_COMPLEX_HALF *evecs, float *resNorms, 
      primme_params *primme);
int cublas_kcprimme_normal(PRIMME_COMPLEX_FLOAT *evals, PRIMME_COMPLEX_HALF *evecs, float *resNorms, 
      primme_params *primme);

primme_params* primme_params_create(void);
int primme_params_destroy(primme_params *primme);
void primme_initialize(primme_params *primme);
int  primme_set_method(primme_preset_method method, primme_params *params);
void primme_display_params(primme_params primme);
void primme_free(primme_params *primme);
int primme_get_member(primme_params *primme, primme_params_label label,
      void *value);
int primme_set_member(primme_params *primme, primme_params_label label,
      void *value);
int primme_member_info(primme_params_label *label, const char** label_name,
      primme_type *type, int *arity);
int primme_constant_info(const char* label_name, int *value);
int primme_enum_member_info(
      primme_params_label label, int *value, const char **value_name);


#ifdef __cplusplus
}
#endif

#endif /* PRIMME_EIGS_H */
