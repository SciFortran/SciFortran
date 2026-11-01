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
 * File: common.h
 *
 * Purpose - Includes headers and defines macros used on this project.
 *           This should be effectively the first file include on every source
 *           file.
 *
 ******************************************************************************/

#ifndef COMMON_H
#define COMMON_H

/* Including MAGMA headers before C99 complex.h avoid compiler issues */

#ifdef PRIMME_WITH_MAGMA

/* NOTE: MAGMA supports half from version 2.5 and still does not support   */
/*       half complex (because CUBLAS does not support it either)          */

#  include <magma_v2.h>

#  if MAGMA_VERSION_MAJOR >= 2 && MAGMA_VERSION_MINOR >= 5
#     define MAGMA_WITH_HALF
#  endif
#endif /* PRIMME_WITH_MAGMA */


#include <limits.h>    
#include <float.h>
#include <stdint.h>
#include <stdlib.h>   /* malloc, free */
#include <string.h>   /* strlen */
#include <stdio.h>    /* snprintf */
#include "primme.h"

/*****************************************************************************/
/* General                                                                   */
/*****************************************************************************/

/* Macros emitted by ctemplate and used here */
#define CONCAT(a, b) CONCATX(a, b)
#define CONCATX(a, b) a ## b
#define STR(X) STR0(X)
#define STR0(X) #X

/*****************************************************************************/
/* Arithmetic                                                                */
/*****************************************************************************/

/* BLAS/LAPACK types: allowed arithmetic (copy, add, multiply) on cpu */

typedef PRIMME_HALF                        dummy_type_hprimme;
typedef PRIMME_COMPLEX_HALF                dummy_type_kprimme;
typedef float                              dummy_type_sprimme;
typedef PRIMME_COMPLEX_FLOAT               dummy_type_cprimme;
typedef double                             dummy_type_dprimme;
typedef PRIMME_COMPLEX_DOUBLE              dummy_type_zprimme;
typedef PRIMME_QUAD                        dummy_type_qprimme;
typedef PRIMME_COMPLEX_QUAD                dummy_type_wprimme;

/* MAGMA types: not allowed arithmetic on cpu */

typedef struct { PRIMME_HALF a; }          dummy_type_magma_hprimme;
typedef struct { PRIMME_COMPLEX_HALF a; }  dummy_type_magma_kprimme;
typedef struct { float a; }                dummy_type_magma_sprimme;
typedef struct { PRIMME_COMPLEX_FLOAT a; } dummy_type_magma_cprimme;
typedef struct { double a; }               dummy_type_magma_dprimme;
typedef struct { PRIMME_COMPLEX_DOUBLE a; }dummy_type_magma_zprimme;
typedef struct { PRIMME_QUAD a; }          dummy_type_magma_qprimme;
typedef struct { PRIMME_COMPLEX_QUAD a; }  dummy_type_magma_wprimme;

/* CUBLAS types: not allowed arithmetic on cpu */

typedef struct { PRIMME_HALF a; }          dummy_type_cublas_hprimme;
typedef struct { PRIMME_COMPLEX_HALF a; }  dummy_type_cublas_kprimme;
typedef struct { float a; }                dummy_type_cublas_sprimme;
typedef struct { PRIMME_COMPLEX_FLOAT a; } dummy_type_cublas_cprimme;
typedef struct { double a; }               dummy_type_cublas_dprimme;
typedef struct { PRIMME_COMPLEX_DOUBLE a; }dummy_type_cublas_zprimme;
typedef struct { PRIMME_QUAD a; }          dummy_type_cublas_qprimme;
typedef struct { PRIMME_COMPLEX_QUAD a; }  dummy_type_cublas_wprimme;



/**********************************************************************
 *
 * Macro SCALAR - type of the matrices problem and the eigenvectors.
 *    It is a fake type, no arithmetic allowed. Compiler will complain when
 *    assigning a value and add, subtract, multiply and divide.
 *
 * Macro REAL - type of the real part of SCALAR, and the type of the
 *    eigenvalues.  It is a fake type, no arithmetic allowed.
 *
 * Macro XSCALAR - type of the matrices problem and the eigenvectors.
 *    Arithmetic is allowed.
 *
 * Macro XREAL - type of the real part of SCALAR, and the type of the
 *    eigenvalues. Arithmetic is allowd
 *
 * Macro SCALAR_SUF - suffix appended in functions whose prototypes have
 *    SCALAR and REAL arguments and is called from other files.
 *
 * Macro REAL_SUF - suffix appended to call the REAL version of the
 *    function.
 *
 * Macro HOST_SCALAR_SUF and HOST_REAL_SUF - suffix appended to call the
 *    CPU version.
 *
 * Macro HSCALAR and HREAL - cpu versions of the types.
 *
 **********************************************************************/

#define SCALAR_SUF      CONCAT(CONCAT(STEM,ARITH(h,k,s,c,d,z,q,w)),primme)
#define XSCALAR_SUF     CONCAT(CONCAT(HOST_STEM,ARITH(h,k,s,c,d,z,q,w)),primme)
#define HOST_SCALAR_SUF CONCAT(CONCAT(HOST_STEM,ARITH(s,c,s,c,d,z,q,w)),primme)

#define SCALAR CONCAT(dummy_type_, SCALAR_SUF)
#define XSCALAR CONCAT(dummy_type_, XSCALAR_SUF)
#define HSCALAR CONCAT(dummy_type_, HOST_SCALAR_SUF)

#define REAL_SUF      CONCAT(CONCAT(STEM,REAL_ARITH(h,k,s,c,d,z,q,w)),primme)
#define XREAL_SUF     CONCAT(CONCAT(HOST_STEM,REAL_ARITH(h,k,s,c,d,z,q,w)),primme)
#define HOST_REAL_SUF CONCAT(CONCAT(HOST_STEM,REAL_ARITH(s,c,s,c,d,z,q,w)),primme)

#define REAL CONCAT(dummy_type_, REAL_SUF)
#define XREAL CONCAT(dummy_type_, XREAL_SUF)
#define HREAL CONCAT(dummy_type_, HOST_REAL_SUF)

#define MACHINE_EPSILON                                                        \
   ARITH(0.000977, 0.000977, FLT_EPSILON, FLT_EPSILON, DBL_EPSILON,            \
         DBL_EPSILON, 1.92593e-34, 1.92593e-34)

#define MACHINE_MAX                                                            \
   ARITH(65504.0, 65504.0, FLT_MAX, FLT_MAX, DBL_MAX, DBL_MAX, FLT128_MAX,     \
         FLT128_MAX)

#define HOST_MACHINE_MAX                                                       \
   ARITH(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, DBL_MAX, DBL_MAX, FLT128_MAX,     \
         FLT128_MAX)

#define PRIMME_OP_SCALAR                                                       \
   ARITH(primme_op_half, primme_op_half, primme_op_float, primme_op_float,     \
         primme_op_double, primme_op_double, primme_op_quad, primme_op_quad)

#define PRIMME_OP_REAL PRIMME_OP_SCALAR

#define PRIMME_OP_HSCALAR                                                      \
   ARITH(primme_op_float, primme_op_float, primme_op_float, primme_op_float,   \
         primme_op_double, primme_op_double, primme_op_quad, primme_op_quad)

#define PRIMME_OP_HREAL PRIMME_OP_HSCALAR

#ifndef __cplusplus
/* pgcc does not provide tgmath.h and the system's one may not support PGI. */
/* Reported issues with icc version 2021 on Ubuntu:
 * https://www.sourceware.org/bugzilla/show_bug.cgi?id=26621 */
#  if defined(__PGI) || (__ICC >= 2021)
#    include <math.h>
#  else
#    include <tgmath.h>   /* select proper function abs from fabs, cabs... */
#  endif
#endif

/* exp and log macros cause warnings on some systems. PRIMME only uses  */
/* the non-complex version of these functions.                          */

#ifdef pow
#  undef pow
#endif
#ifdef exp
#  undef exp
#endif
#ifdef log
#  undef log
#endif

#ifndef __cplusplus
#  define ISFINITE isfinite
#else
#  define ISFINITE std::isfinite
#endif

/*****************************************************************************/
/* Memory management                                                         */
/*****************************************************************************/

/* Recent versions of PRIMME are using dynamic memory to manage the memory.  */
/* In general the use of dynamic memory simplifies the code by not needing   */
/* to take care of providing enough working space for all functions in       */
/* PRIMME. The small drawback of dynamic memory is to mingle with error      */
/* management. C++ shines in this situation. But we are restricted to C,     */
/* probably, not for good reasons. The goal is to avoid writing specific code*/
/* to free allocated memory in case of an error happening in the body of a   */
/* function. The next macros together with the functions in memman.c         */
/* provide an interface to track memory allocations and free them in case    */
/* of error. These macros and functions are going to be used mostly by       */
/* error management macros, eg CHKERR, and memory allocation functions, eg   */
/* Num_malloc_Sprimme. The error manager is going to assume that every       */
/* allocation is going to be freed at the end of the function. If            */
/* that is not the case, call Mem_keep_frame(ctx). Examples of this are the  */
/* functions Num_malloc_Sprimme.                                             */

#include "memman.h"

/**********************************************************************
 * Macro MEM_PUSH_FRAME - Set a new frame in which new allocations are going
 *    to be registered.
 *
 * NOTE: only one call is allowed in the same block of code (anything between
 *       curly brackets).
 **********************************************************************/

#define MEM_PUSH_FRAME \
   primme_frame __frame = {NULL, 0, ctx.mm}; \
   ctx.mm = &__frame;

/**********************************************************************
 * Macro MEM_POP_FRAME(ERRN) - If *ERRN is 0, it just removes all registered
 *    allocations. Otherwise it also frees all allocations.
 *
 * NOTE: only one call is allowed in the same block of code (anything between
 *       curly brackets).
 **********************************************************************/

#ifndef NDEBUG
#define MEM_POP_FRAME(ERRN) { \
   if (*(ERRN)) {\
      Mem_pop_clean_frame(ctx);\
   } else {\
      Mem_debug_frame(__FILE__ ": " STR(__LINE__), ctx);\
      if (Mem_pop_frame(&ctx)) { \
         Mem_pop_clean_frame(ctx); \
         PRINTFALL(1, "PRIMME: Error popping frame, most likely forgotten call to Mem_keep_frame.");\
         *(ERRN) = -1; \
      }\
   }\
}
#else
#define MEM_POP_FRAME(ERRN) { \
   if (*(ERRN)) {\
      Mem_pop_clean_frame(ctx);\
   } else {\
      if (Mem_pop_frame(&ctx)) { \
         Mem_pop_clean_frame(ctx); \
         PRINTFALL(1, "PRIMME: Error popping frame, most likely forgotten call to Mem_keep_frame.");\
         *(ERRN) = -1; \
      }\
   }\
}
#endif

/*****************************************************************************/
/* Reporting                                                                 */
/*****************************************************************************/

/**********************************************************************
 * Macro PRINTFALLCTX - report some information. It invokes ctx.report if the
 *    level >= CTX.printLevel regardless of the processor's id.
 *
 * INPUT PARAMETERS
 * ----------------
 * L     print level
 * ...   printf arguments
 *
 * EXAMPLE
 * -------
 *   PRINTFALLCTX(ctx, 2, "Alpha=%f is to close to zero", alpha);
 *
 **********************************************************************/

#define PRINTFALLCTX(CTX, L, ...)                                              \
   {                                                                           \
      if ((CTX).report && (L) <= (CTX).printLevel) {                           \
         int len = snprintf(NULL, 0, __VA_ARGS__) + 1;                         \
         char *str = (char *)malloc(len);                                      \
         snprintf(str, len, __VA_ARGS__);                                      \
         (CTX).report(str, -1.0, (CTX));                                       \
         free(str);                                                            \
      }                                                                        \
   }

/**********************************************************************
 * Macro PRINTFALL - report some information. It invokes ctx.report if the
 *    level >= ctx.printLevel regardless of the processor's id.
 *
 * INPUT PARAMETERS
 * ----------------
 * L     print level
 * ...   printf arguments
 *
 * EXAMPLE
 * -------
 *   PRINTFALL(2, "Alpha=%f is to close to zero", alpha);
 *
 **********************************************************************/

#define PRINTFALL(L, ...) PRINTFALLCTX(ctx, (L), __VA_ARGS__)

/**********************************************************************
 * Macro PRINTF - report some information. It invokes ctx.report if the
 *    level >= ctx.printLevel and it is processor 0.
 *
 * INPUT PARAMETERS
 * ----------------
 * L     print level
 * ...   printf arguments
 *
 * EXAMPLE
 * -------
 *   PRINTF(2, "Alpha=%f is to close to zero", alpha);
 *
 **********************************************************************/

#define PRINTF(L, ...) {if (ctx.procID == 0) PRINTFALL((L), __VA_ARGS__);}

/*****************************************************************************/
/* Profiling                                                                 */
/*****************************************************************************/

#include "wtime.h"

#ifdef PRIMME_PROFILE

#include <regex.h>

static inline const char *__compose_function_name(const char *path,
      const char *call, const char *filename, const char *line) {

   int len = strlen(path) + 1 + 45 + strlen(filename) + 1 + strlen(line) + 10;
   char *s = (char*)malloc(len);
   snprintf(s, len-1, "%s~%.40s@%s:%s", path, call, filename, line);
   s[len-1] = 0;
   return s;
}

#define PROFILE_BEGIN(CALL) \
   double ___t0 = 0; \
   const char *old_path = ctx.path; \
   if (ctx.path && ctx.report) { \
      ctx.path = __compose_function_name(ctx.path, CALL, __FILE__, STR(__LINE__)); \
      if (regexec(&ctx.profile, ctx.path, 0, NULL, 0) == 0) { \
            ctx.report(ctx.path, -.5, ctx); \
         ___t0 = primme_wTimer() - *ctx.timeoff; \
      } \
   }

#define PROFILE_END \
   if (ctx.path) { \
      if (___t0 > 0) { \
         double ___t1 = primme_wTimer(); \
         ___t0 = ___t1 - ___t0 - *ctx.timeoff; \
         if (___t0 > 0 && ctx.report) { \
            ctx.report(ctx.path, ___t0, ctx); \
            *ctx.timeoff += primme_wTimer() - ___t1; \
         } \
      } \
      free((void*)ctx.path); \
      ctx.path = old_path; \
   }
 
#else
#define PROFILE_BEGIN(CALL)
#define PROFILE_END
#endif

/*****************************************************************************/
/* Parallel checks                                                           */
/*****************************************************************************/

static inline uint32_t hash_call(const char *str, double value) {
   uint32_t hash = 5381;
   int c;

   union { double d; uint32_t i[2]; } t;
   t.d = value;
   while ((c = *str++)) hash = hash * 33 + c;
   hash = hash * 33 + t.i[0];
   hash = hash * 33 + t.i[1];

   return hash;
}

/**********************************************************************
 * Macro PARALLEL_CHECK(CALL) - Check all processes are executing CALL
 *    and it returns the same value. If a process is executing a different
 *    PARALLEL_CHECK or CALL returns a different value, it is returned
 *    PRIMME_PARALLEL_FAILURE.
 *
 *    NOTE: CALL is only evaluated once. Avoid to put several
 *    PARALLEL_CHECK with same expression CALL on the same line.
 *
 * INPUT PARAMETERS
 * ----------------
 * CALL    value that is check
 *
 **********************************************************************/

#define PARALLEL_CHECK(CALL)                                                   \
   {                                                                           \
      double __value = CALL;                                                   \
      uint32_t __hash_call =                                                   \
                     hash_call(STR(CALL) __FILE__ STR(__LINE__), __value),     \
               __hash_call0 = __hash_call;                                     \
      CHKERR(ctx.bcast(&__hash_call0, primme_op_float, 1, ctx));               \
      float __not_is_equal = (__hash_call != __hash_call0 ? 1 : 0),            \
            __not_is_equal_global = __not_is_equal;                            \
      CHKERR(ctx.globalSum(&__not_is_equal_global, primme_op_float, 1, ctx));  \
      if (__not_is_equal_global != 0 && ctx.procID == 0) {                     \
         PRINTFALL(1,                                                          \
               "PRIMME: Process 0 has value %12g for '" STR(                   \
                     CALL) "' (" __FILE__ ":" STR(__LINE__) ")\n",             \
               __value);                                                       \
      }                                                                        \
      if (__not_is_equal != 0) {                                               \
         PRINTFALL(1,                                                          \
               "PRIMME: Process %d has different value %12g for '" STR(        \
                     CALL) "' (" __FILE__ ":" STR(__LINE__) ")\n",             \
               ctx.procID, __value);                                   \
      }                                                                        \
      if (__not_is_equal_global != 0) CHKERR(PRIMME_PARALLEL_FAILURE);         \
   }


/*****************************************************************************/
/* Error management                                                          */
/*****************************************************************************/

/* The error code is checked with an assert. Then by default the code will   */
/* abort in case of error, unless it is compiled with the macro NDEBUG       */
/* defined. In that case it is printed out the file and the line in the      */
/* source that caused the error in the file descriptor primme->outputFile    */
/* or primme_svds->outputFile. Use CHKERR and CHKERRM in functions that has  */
/* primme_params* primme, and CHKERRS and CHKERRMS in functions that has     */
/* primme_svds_params* primme_svds.                                          */

#include <assert.h>

/**********************************************************************
 * Macro CHKERR - assert that ERRN == 0. Otherwise it is printed out in
 *    ctx.outputFile the file and the line of the caller, and
 *    forced the caller function to return ERRN.
 *
 *    ERRN is only evaluated once.
 *
 * INPUT PARAMETERS
 * ----------------
 * ERRN    Expression that returns an error code
 *
 **********************************************************************/

#define CHKERR(ERRN) { \
   MEM_PUSH_FRAME; \
   PROFILE_BEGIN(STR(ERRN)); \
   int __err = (ERRN); assert(__err==0);\
   PROFILE_END; \
   MEM_POP_FRAME(&__err); \
   if (__err) {\
      PRINTFALL(1, "PRIMME: Error %d in (" __FILE__ ":%d): %s", __err, __LINE__, #ERRN );\
      return __err;\
   }\
}

/**********************************************************************
 * Macro CHKERRM - assert that ERRN == 0. Otherwise it is printed out in
 *    ctx.outputFile the file and the line of the caller, in
 *    addition to an error message passed as arguments. As CHKERR
 *    the caller function is forced to return RETURN in case of error.
 *
 *    ERRN is only evaluated once.
 *
 * INPUT PARAMETERS
 * ----------------
 * ERRN    Expression that returns an error code
 * RETURN  Value that the caller function will return in case of error
 *
 * EXAMPLE
 * -------
 *   CHKERRM((ptr = malloc(n*sizeof(double))) == NULL, -1,
 *        "malloc could not allocate %d doubles\n", n);
 *
 **********************************************************************/

#define CHKERRM(ERRN, RETURN, ...) { \
   MEM_PUSH_FRAME; \
   PROFILE_BEGIN(STR(ERRN)); \
   int __err = (ERRN); assert(__err==0);\
   PROFILE_END; \
   MEM_POP_FRAME(&__err); \
   if (__err) {\
      PRINTFALL(1, "PRIMME: Error %d in (" __FILE__ ":%d): %s", __err, __LINE__, #ERRN );\
      PRINTFALL(1, "PRIMME: " __VA_ARGS__);\
      return (RETURN);\
   }\
}


/**********************************************************************
 * Macro CHKERRA - assert that ERRN == 0. Otherwise it is printed out in
 *    ctx.outputFile the file and the line of the caller, and
 *    forced the caller function to execute ACTION.
 *
 *    ERRN is only evaluated once.
 *
 * INPUT PARAMETERS
 * ----------------
 * ERRN    Expression that returns an error code
 * ACTION  Expression that the caller function will execute in case of error
 *
 **********************************************************************/

#define CHKERRA(ERRN, ACTION) { \
   MEM_PUSH_FRAME; \
   PROFILE_BEGIN(STR(ERRN)); \
   int __err = (ERRN); assert(__err==0);\
   PROFILE_END; \
   MEM_POP_FRAME(&__err); \
   if (__err) {\
      PRINTFALL(1, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", __err, __LINE__, #ERRN );\
      ACTION;\
      return;\
   } \
}

/**********************************************************************
 * Macro CHKERRVAL - set the error code on the second macro argument.
 *
 *    ERRN is only evaluated once.
 *
 * INPUT PARAMETERS
 * ----------------
 * ERRN    Expression that returns an error code
 * RET     Variable where to assign the error code
 *
 **********************************************************************/

#define CHKERRVAL(ERRN, RET) { \
   MEM_PUSH_FRAME; \
   PROFILE_BEGIN(STR(ERRN)); \
   *(RET) = (ERRN);\
   PROFILE_END; \
   MEM_POP_FRAME((RET)); \
   if (*(RET)) {\
      PRINTFALL(1, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", *(RET), __LINE__, #ERRN );\
   } \
}

/*****************************************************************************/
/* Memory, error and params  management                                      */
/*****************************************************************************/


/**********************************************************************
 * Macro MALLOC_PRIMME - malloc NELEM of type **X and write down the
 *    pointer in *X.
 *
 * INPUT PARAMETERS
 * ----------------
 * NELEM   Number of sizeof(**X) to allocate
 * X       Where to store the pointer
 *
 * RETURN VALUE
 * ------------
 * error code
 *
 * EXAMPLE
 * -------
 *   double *values;
 *   CHKERRM(MALLOC_PRIMME(n, &values), -1,
 *        "malloc could not allocate %d doubles\n", n);
 *
 **********************************************************************/

#define MALLOC_PRIMME(NELEM, X)                                                \
   (*((void **)X) = malloc((NELEM) * sizeof(**(X))),                           \
         *(X) == NULL ? PRIMME_MALLOC_FAILURE : 0)

typedef struct primme_context_str {
   /* For PRIMME */
   primme_params *primme;
   primme_svds_params *primme_svds;

   /* For output */
   int printLevel;
   FILE *outputFile;
   int (*report)(const char *fun, double time, struct primme_context_str ctx);

   /* For memory management */
   primme_frame *mm;

   /* for MPI */
   int numProcs;     /* number of processes */
   int procID;       /* process id */
   void *mpicomm;    /* MPI communicator */
   int (*bcast)(void *buffer, primme_op_datatype buffer_type, int count,
         struct primme_context_str ctx); /* broadcast */
   int (*globalSum)(void *buffer, primme_op_datatype buffer_type, int count,
         struct primme_context_str ctx); /* global reduction */

   /* For MAGMA */
   void *queue;      /* magma device queue (magma_queue_t*) */

   #ifdef PRIMME_PROFILE
   /* For profiling */
   regex_t profile;  /* Pattern of the functions to profile */
   const char *path; /* Path of the current function */
   double *timeoff;  /* Accumulated overhead of profiling and reporting */
   #endif
} primme_context;

/*****************************************************************************/
/* Miscellanea                                                               */
/*****************************************************************************/

#define TO_INT(X) ((X) < INT_MAX ? (X) : INT_MAX)

#ifdef F77UNDERSCORE
#define FORTRAN_FUNCTION(X) CONCAT(X,_)
#else
#define FORTRAN_FUNCTION(X) X
#endif

#ifndef max
#  define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#  define min(a, b) ((a) < (b) ? (a) : (b))
#endif

/* Number of rows processed as once for streaming operation */
#define PRIMME_BLOCK_SIZE IMPL(512, INT_MAX)

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif

#define AS_FORTRAN(X) FORTRAN_FUNCTION(CONCAT(X,_f77))

#endif /* COMMON_H */
