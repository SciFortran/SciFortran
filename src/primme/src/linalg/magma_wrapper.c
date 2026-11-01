/*******************************************************************************
 * Copyright (c) 2017, College of William & Mary
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
 * File: magma_wrapper.c
 *
 * Purpose - This file contains mostly C wrapper routines for
 *           calling MAGMA functions
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../linalg/magma_wrapper.c"
#endif

#include "numerical.h"

#ifdef SUPPORTED_TYPE

#ifdef USE_MAGMA

#if defined(USE_HALF_MAGMA)
#  define USE_CUBLAS_HERE
#endif

#ifndef CHECK_TEMPLATE
#  include <cuda_runtime.h>
#  ifdef USE_CUBLAS_HERE
#    include <cublas_v2.h>
#  endif
#endif

/* The next functions are only exposed using C++.                    */
/* These signatures are copied from MAGMA header file magmablas_h.h. */

#if !defined(MAGMA_WRAPPER_PRIVATE_HALF) && defined(MAGMA_WITH_HALF) && !defined(__cplusplus)
#define MAGMA_WRAPPER_PRIVATE_HALF

void
magmablas_slag2h(
    magma_int_t m, magma_int_t n,
    float const * dA, magma_int_t lda,
    magmaHalf* dHA, magma_int_t ldha,
    magma_int_t *info, magma_queue_t queue);

void
magmablas_hlag2s(
    magma_int_t m, magma_int_t n,
    magmaHalf_const_ptr dA, magma_int_t lda,
    float             *dSA, magma_int_t ldsa,
    magma_queue_t queue );

void
magma_hgemm(
    magma_trans_t transA, magma_trans_t transB,
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaHalf alpha,
    magmaHalf_const_ptr dA, magma_int_t ldda,
    magmaHalf_const_ptr dB, magma_int_t lddb,
    magmaHalf beta,
    magmaHalf_ptr       dC, magma_int_t lddc,
    magma_queue_t queue );

#endif

#ifndef MAGMA_WRAPPER_PRIVATE
#define MAGMA_WRAPPER_PRIVATE

#define MAGMA_SCALAR                                                           \
   ARITH(magmaHalf, , float, magmaFloatComplex, double, magmaDoubleComplex, , )

#define XLASET  CONCAT(magmablas_,ARITH( , , slaset, claset, dlaset, zlaset, , ))

#define XCOPY     CONCAT(magma_,ARITH( , , scopy , ccopy , dcopy , zcopy , , ))   
#define XGEMM     CONCAT(magma_,ARITH( , , sgemm , cgemm , dgemm , zgemm , , ))
#define XGEMV     CONCAT(magma_,ARITH( , , sgemv , cgemv , dgemv , zgemv , , ))
#define XAXPY     CONCAT(magma_,ARITH( , , saxpy , caxpy , daxpy , zaxpy , , ))
#define XDOT      CONCAT(magma_,ARITH( , , sdot  , cdotc , ddot  , zdotc , , ))
#define XSCAL     CONCAT(magma_,ARITH( , , sscal , cscal , dscal , zscal , , ))
#define XTRSM     CONCAT(magma_,ARITH( , , strsm , ctrsm , dtrsm , ztrsm , , ))

static int free_fn_dummy (void *p, primme_context ctx) {
   (void)ctx;
   return magma_free(p) == MAGMA_SUCCESS ? 0 : PRIMME_MALLOC_FAILURE;
}

#endif /* MAGMA_WRAPPER_PRIVATE */

/******************************************************************************
 * Function Num_check_pointer - Return no error code if the pointer is on a
 *    device allocation.
 *
 * PARAMETERS
 * ---------------------------
 * x           pointer to check
 * 
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_check_pointer_Sprimme(void *x) {
   struct cudaPointerAttributes ptr_attr;
   if (cudaPointerGetAttributes(&ptr_attr, x) != cudaSuccess) return -1;

#if CUDART_VERSION >= 10000
   if (ptr_attr.type == cudaMemoryTypeHost) return -1;
#else
   if (!ptr_attr.isManaged && ptr_attr.memoryType == cudaMemoryTypeHost)
      return -1;
#endif
   return 0;
}

/******************************************************************************
 * Function Num_malloc - Allocate a vector of scalars
 *
 * PARAMETERS
 * ---------------------------
 * n           The number of elements
 * v           returned pointer
 * 
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_malloc_Sprimme(PRIMME_INT n, SCALAR **x, primme_context ctx) {

   /* Quick exit */

   if (n <= 0) {
      *x = NULL;
      return 0;
   }

   /* Allocate memory */

   if (magma_malloc((void**)x, n * sizeof(SCALAR)) != MAGMA_SUCCESS)
      return PRIMME_MALLOC_FAILURE;

   /* Register the allocation */

   Mem_keep_frame(ctx);
   Mem_register_alloc(*x, free_fn_dummy, ctx);

   return 0;
}

/******************************************************************************
 * Function Num_free - Free allocated a vector of scalars
 *
 * PARAMETERS
 * ---------------------------
 * v           allocated pointer
 * 
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_free_Sprimme(SCALAR *x, primme_context ctx) {

   /* Quick exit */

   if (!x) return 0;

   /* Deregister the allocation */

   Mem_deregister_alloc(x, ctx);

   /* Free pointer */

   return magma_free(x) == MAGMA_SUCCESS ? 0 : PRIMME_MALLOC_FAILURE;
}

/******************************************************************************
 * Function Num_copy_Tmatrix - Copy the matrix x with type given as a parameter
 *    into y
 *
 * PARAMETERS
 * ---------------------------
 * x           The source matrix
 * xt          The type of x
 * m           The number of rows of x
 * n           The number of columns of x
 * ldx         The leading dimension of x
 * y           On output y = x
 * ldy         The leading dimension of y
 *
 * NOTE: x and y *can* overlap
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_copy_Tmatrix_Sprimme(void *x, primme_op_datatype xt, PRIMME_INT m,
      PRIMME_INT n, PRIMME_INT ldx, SCALAR *y, PRIMME_INT ldy,
      primme_context ctx) {
   (void)ctx;

   /* Quick exit */

   if (xt == primme_op_default || xt == PRIMME_OP_SCALAR) {
      CHKERR(Num_copy_matrix_Sprimme((SCALAR*)x, m, n, ldx, y, ldy, ctx));
      return 0;
   }

   if (m == 0 || n == 0) return 0;

   /* In-place casting is not supported */

   if (x == y) return PRIMME_FUNCTION_UNAVAILABLE;

   /* Call the equivalent real version if needed. MAGMA only provides casting */
   /* functions for real.                                                     */

#ifdef USE_COMPLEX
   return Num_copy_Tmatrix_Rprimme(
         x, xt, m * 2, n, ldx * 2, (REAL *)y, ldy * 2, ctx);
#else

   int info;
   void *y0 = y; PRIMME_INT ldy0 = ldy;
   switch (xt) {
#  ifdef MAGMA_WITH_HALF
      case primme_op_half:
         /* If y is float, cast x directly to float with halg2s. Otherwise    */
         /* cast x to float into an auxiliary array, and then cast the last   */
         /* one into y.                                                       */

         if (PRIMME_OP_SCALAR > primme_op_float) {
            CHKERR(Num_malloc_Rsprimme(
                  m * n, (dummy_type_magma_sprimme **)&y0, ctx));
            ldy0 = m;
         }
         magmablas_hlag2s(m, n, (magmaHalf *)x, ldx, (float *)y0, ldy0,
               *(magma_queue_t *)ctx.queue);
         if (y0 != y) {
            CHKERR(Num_copy_Tmatrix_Rprimme(
                  y0, primme_op_float, m, n, ldy0, y, ldy, ctx));
            CHKERR(Num_free_Rsprimme((dummy_type_magma_sprimme *)y0, ctx));
         }
         break;
#  endif /* MAGMA_WITH_HALF */

      case primme_op_float:
#  if defined(MAGMA_WITH_HALF) && defined(USE_HALF_MAGMA)
         /* Cast x into y from float to half */

         magmablas_slag2h(m, n, (float *)x, ldx, (magmaHalf *)y, ldy, &info,
               *(magma_queue_t *)ctx.queue);
         CHKERRM(info != 0, PRIMME_UNEXPECTED_FAILURE,
               "Error in slag2h with info %d", info);

#  elif defined(USE_DOBULE_MAGMA)
         /* Cast x into y from float to double */

         magmablas_slag2d(m, n, (float *)x, ldx, (double *)y, ldy,
               *(magma_queue_t *)ctx.queue, &info);
         CHKERRM(info != 0, PRIMME_UNEXPECTED_FAILURE,
               "Error in slag2d with info %d", info);
#  else
         /* Quadruple precision is not supported */

         return PRIMME_FUNCTION_UNAVAILABLE;
#  endif /* defined(MAGMA_WITH_HALF) && defined(USE_HALF_MAGMA) */
         break;

      case primme_op_double:
         /* Quadruple precision is not supported */

        if (PRIMME_OP_SCALAR > primme_op_double)
            return PRIMME_FUNCTION_UNAVAILABLE;

         /* If y is float, cast x directly to float with dalg2s. Otherwise    */
         /* cast x to float into an auxiliary array, and then cast the last   */
         /* one into y.                                                       */

         if (PRIMME_OP_SCALAR < primme_op_float) {
            CHKERR(Num_malloc_Rsprimme(
                  m * n, (dummy_type_magma_sprimme **)&y0, ctx));
            ldy0 = m;
         }
         magmablas_dlag2s(m, n, (double *)x, ldx, (float *)y0, ldy0,
               *(magma_queue_t *)ctx.queue, &info);
         CHKERRM(info != 0, PRIMME_UNEXPECTED_FAILURE,
               "Error in slag2d with info %d", info);
         if (y0 != y) {
            CHKERR(Num_copy_Tmatrix_Rprimme(
                  y0, primme_op_float, m, n, ldy0, y, ldy, ctx));
            CHKERR(Num_free_Rsprimme((dummy_type_magma_sprimme *)y0, ctx));
         }
         break;
      default: CHKERR(PRIMME_FUNCTION_UNAVAILABLE);
   }

   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   return 0;
#endif /* USE_COMPLEX */
}

/*******************************************************************************
 * Subroutine Num_copy_Sprimme - y(0:n*incy-1:incy) = x(0:n*incx-1:incx)
 ******************************************************************************/
 
TEMPLATE_PLEASE
int Num_copy_Sprimme(PRIMME_INT n, SCALAR *x, int incx, SCALAR *y, int incy,
      primme_context ctx) {

   CHKERR(Num_copy_matrix_Sprimme(x, 1, n, incx, y, incy, ctx));

   return 0;
}

/*******************************************************************************
 * Subroutine Num_gemm_Sprimme - C = op(A)*op(B), with C size m x n
 * NOTE: A, B and C are in device
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_gemm_Sprimme(const char *transa, const char *transb, int m, int n,
      int k, HSCALAR alpha, SCALAR *a, int lda, SCALAR *b, int ldb,
      HSCALAR beta, SCALAR *c, int ldc, primme_context ctx) {

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return 0;

   /* Quick exit */

   if (k == 0 || ABS(alpha) == 0.0) {
      if (ABS(beta) == 0.0) {
         Num_zero_matrix_Sprimme(c, m, n, ldc, ctx);
      }
      else if (beta != (HSCALAR)1.0) {
         int i;
         for (i=0; i<n; i++) {
            Num_scal_Sprimme(m, beta, &c[ldc*i], 1, ctx);
         }
      }
      return 0;
   }
   if (n == 1) {
      PRIMME_INT mA; int nA;
      if (*transa == 'n' || *transa == 'N') mA = m, nA = k;
      else mA = k, nA = m;
      int incb = ((*transb == 'n' || *transb == 'N') ? 1 : ldb);
      return Num_gemv_Sprimme(
            transa, mA, nA, alpha, a, lda, b, incb, beta, c, 1, ctx);
   }

#if defined(USE_HALFCOMPLEX_MAGMA)
   return PRIMME_FUNCTION_UNAVAILABLE;

#elif defined(USE_HALF_MAGMA)
   magma_hgemm(magma_trans_const(*transa), magma_trans_const(*transb), m, n, k,
         *(MAGMA_SCALAR *)&alpha, (MAGMA_SCALAR *)a, lda, (MAGMA_SCALAR *)b,
         ldb, *(MAGMA_SCALAR *)&beta, (MAGMA_SCALAR *)c, ldc,
         *(magma_queue_t *)ctx.queue);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   return 0;

#else
   XGEMM(magma_trans_const(*transa), magma_trans_const(*transb), m, n, k,
         *(MAGMA_SCALAR *)&alpha, (MAGMA_SCALAR *)a, lda, (MAGMA_SCALAR *)b,
         ldb, *(MAGMA_SCALAR *)&beta, (MAGMA_SCALAR *)c, ldc,
         *(magma_queue_t *)ctx.queue);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   return 0;
#endif
}


/*******************************************************************************
 * Subroutine Num_gemm_dhd_Sprimme - C = op(A)*op(B), with C size m x n
 * NOTE: A and C are hosted on device and B is hosted on cpu
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_gemm_dhd_Sprimme(const char *transa, const char *transb, int m, int n,
      int k, HSCALAR alpha, SCALAR *a, int lda, HSCALAR *b, int ldb,
      HSCALAR beta, SCALAR *c, int ldc, primme_context ctx) {

   /* Quick exit */

   if (m == 0 || n == 0) return 0;

#if defined(USE_HALFCOMPLEX_MAGMA)
   (void)transa;
   (void)transb;
   (void)m;
   (void)n;
   (void)k;
   (void)alpha;
   (void)a;
   (void)lda;
   (void)b;
   (void)ldb;
   (void)beta;
   (void)c;
   (void)ldc;
   (void)ctx;
   return PRIMME_FUNCTION_UNAVAILABLE;

#else
   int mb = (*transb == 'N' || *transb == 'n') ? k : n;
   int nb = (*transb == 'N' || *transb == 'n') ? n : k;

   SCALAR *b_dev; /* copy of b on device */
   CHKERR(Num_malloc_Sprimme(mb*nb, &b_dev, ctx));
   if (mb != 0 && nb != 0) {
      magma_setmatrix(mb, nb, sizeof(SCALAR), (MAGMA_SCALAR *)b, ldb,
            (MAGMA_SCALAR *)b_dev, mb, *(magma_queue_t *)ctx.queue);
      CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
            "Unexpected CUDA error!");
   }
   CHKERR(Num_gemm_Sprimme(
         transa, transb, m, n, k, alpha, a, lda, b_dev, mb, beta, c, ldc, ctx));
   CHKERR(Num_free_Sprimme(b_dev, ctx));

   return 0;
#endif /* USE_HALFCOMPLEX_MAGMA */
}

/*******************************************************************************
 * Subroutine Num_gemm_ddh_Sprimme - C = op(A)*op(B), with C size m x n
 * NOTE: A and B are hosted on device and C is hosted on cpu
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_gemm_ddh_Sprimme(const char *transa, const char *transb, int m, int n,
      int k, HSCALAR alpha, SCALAR *a, int lda, SCALAR *b, int ldb,
      HSCALAR beta, HSCALAR *c, int ldc, primme_context ctx) {

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return 0;

#if defined(USE_HALFCOMPLEX_MAGMA)
   (void)transa;
   (void)transb;
   (void)m;
   (void)n;
   (void)k;
   (void)alpha;
   (void)a;
   (void)lda;
   (void)b;
   (void)ldb;
   (void)beta;
   (void)c;
   (void)ldc;
   (void)ctx;
   return PRIMME_FUNCTION_UNAVAILABLE;

#else
   SCALAR *c_dev; /* copy of c on device */
   CHKERR(Num_malloc_Sprimme(m*n, &c_dev, ctx));
   if (ABS(beta) != 0) {
      magma_setmatrix(m, n, sizeof(SCALAR), (MAGMA_SCALAR *)c, ldc,
            (MAGMA_SCALAR *)c_dev, m, *(magma_queue_t *)ctx.queue);
      CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
            "Unexpected CUDA error!");
   }
   else {
      Num_zero_matrix_Sprimme(c_dev, m, n, m, ctx);
   }
   CHKERR(Num_gemm_Sprimme(
         transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c_dev, m, ctx));
   magma_getmatrix(m, n, sizeof(SCALAR), (MAGMA_SCALAR *)c_dev, m,
         (MAGMA_SCALAR *)c, ldc, *(magma_queue_t *)ctx.queue);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   CHKERR(Num_free_Sprimme(c_dev, ctx));

   return 0;
#endif /* USE_HALFCOMPLEX_MAGMA */
}

/*******************************************************************************
 * Subroutine Num_gemv_Sprimme - y = alpha*A*x + beta*y, with A size m x n
 * NOTE: A, x and y are in device
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_gemv_Sprimme(const char *transa, PRIMME_INT m, int n, HSCALAR alpha,
      SCALAR *a, int lda, SCALAR *x, int incx, HSCALAR beta, SCALAR *y,
      int incy, primme_context ctx) {

   /* Zero dimension matrix may cause problems */
   int tA = (*transa != 'n' && *transa != 'N' ? 1 : 0);
   PRIMME_INT mA = tA ? n : m, nA = tA ? m : n;
   if (mA == 0) return 0;

   /* Quick exit */

   if (nA == 0) {
      if (ABS(beta) == 0.0) {
         Num_zero_matrix_Sprimme(y, 1, mA, incy, ctx);
      }
      else {
         Num_scal_Sprimme(mA, beta, y, incy, ctx);
      }
      return 0;
   }

#if defined(USE_HALFCOMPLEX_MAGMA)
   (void)alpha;
   (void)a;
   (void)lda;
   (void)x;
   (void)incx;
   return PRIMME_FUNCTION_UNAVAILABLE;

#elif defined(USE_HALF_MAGMA)
   /* Neither CUBLAS nor MAGMA provide gemv for half precision. So we do an   */
   /* equivalent call to gemm. Consider that x and y are given as row vectors */
   /* with leading dimensions incx and incy respectively. We like to do       */
   /* A * x' = y', but gemm does not allow returning a matrix implicit        */
   /* transposed. We do instead x * A' = y. Conjugacy required for complex!   */

   CHKERR(Num_gemm_Sprimme("N", tA ? "N" : "C", 1, mA, nA, alpha, x, incx, a,
         lda, beta, y, incy, ctx));
   return 0;

#else
   XGEMV(magma_trans_const(*transa), m, n, *(MAGMA_SCALAR *)&alpha,
         (MAGMA_SCALAR *)a, lda, (MAGMA_SCALAR *)x, incx,
         *(MAGMA_SCALAR *)&beta, (MAGMA_SCALAR *)y, incy,
         *(magma_queue_t *)ctx.queue);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   return 0;
#endif
}

/*******************************************************************************
 * Subroutine Num_gemv_ddh_Sprimme - y = alpha*A*x + beta*y, with A size m x n
 * NOTE: A and x are in device and y is in host
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_gemv_ddh_Sprimme(const char *transa, PRIMME_INT m, int n, HSCALAR alpha,
      SCALAR *a, int lda, SCALAR *x, int incx, HSCALAR beta, HSCALAR *y,
      int incy, primme_context ctx) {

 #if defined(USE_HALFCOMPLEX_MAGMA)
   (void)transa;
   (void)m;
   (void)n;
   (void)alpha;
   (void)a;
   (void)lda;
   (void)x;
   (void)incx;
   (void)beta;
   (void)y;
   (void)incy;
   (void)ctx;
   return PRIMME_FUNCTION_UNAVAILABLE;

#else
   int my = (*transa == 'N' || *transa == 'n') ? m : n;
   if (my == 0) return 0;

   SCALAR *y_dev; /* copy of y on device */
   CHKERR(Num_malloc_Sprimme(my, &y_dev, ctx));
   if (ABS(beta) != 0) {
      magma_setvector(my, sizeof(SCALAR), (MAGMA_SCALAR *)y, incy,
            (MAGMA_SCALAR *)y_dev, 1, *(magma_queue_t *)ctx.queue);
   } else {
      Num_zero_matrix_Sprimme(y_dev, my, 1, my, ctx);
   }
   CHKERR(Num_gemv_Sprimme(
         transa, m, n, alpha, a, lda, x, incx, beta, y_dev, 1, ctx));
   magma_getvector(my, sizeof(SCALAR), (MAGMA_SCALAR *)y_dev, 1,
         (MAGMA_SCALAR *)y, incy, *(magma_queue_t *)ctx.queue);
   CHKERR(Num_free_Sprimme(y_dev, ctx));

   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   return 0;
#endif /* USE_HALFCOMPLEX_MAGMA */
}

/*******************************************************************************
 * Subroutine Num_gemv_dhd_Sprimme - y = alpha*A*x + beta*y, with A size m x n
 * NOTE: A and y are in device and x is in host
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_gemv_dhd_Sprimme(const char *transa, PRIMME_INT m, int n, HSCALAR alpha,
      SCALAR *a, int lda, HSCALAR *x, int incx, HSCALAR beta, SCALAR *y,
      int incy, primme_context ctx) {

 #if defined(USE_HALFCOMPLEX_MAGMA)
   (void)transa;
   (void)m;
   (void)n;
   (void)alpha;
   (void)a;
   (void)lda;
   (void)x;
   (void)incx;
   (void)beta;
   (void)y;
   (void)incy;
   (void)ctx;
   return PRIMME_FUNCTION_UNAVAILABLE;

#else
   int mx = (*transa == 'N' || *transa == 'n') ? n : m;
   if (mx == 0) return 0;

   SCALAR *x_dev; /* copy of x on device */
   CHKERR(Num_malloc_Sprimme(mx, &x_dev, ctx));
   magma_setvector(mx, sizeof(SCALAR), (MAGMA_SCALAR *)x, incx,
         (MAGMA_SCALAR *)x_dev, 1, *(magma_queue_t *)ctx.queue);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   CHKERR(Num_gemv_Sprimme(
         transa, m, n, alpha, a, lda, x_dev, 1, beta, y, incy, ctx));
   CHKERR(Num_free_Sprimme(x_dev, ctx));

   return 0;
#endif /* USE_HALFCOMPLEX_MAGMA */
}

/*******************************************************************************
 * Subroutine Num_axpy_Sprimme - y += alpha*x
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_axpy_Sprimme(PRIMME_INT n, HSCALAR alpha, SCALAR *x, int incx, 
   SCALAR *y, int incy, primme_context ctx) {

   /* Zero dimension matrix may cause problems */
   if (n == 0) return 0;

#if defined(USE_HALFCOMPLEX_MAGMA)
   (void)alpha;
   (void)x;
   (void)incx;
   (void)y;
   (void)incy;
   (void)ctx;
   return PRIMME_FUNCTION_UNAVAILABLE;

#elif defined(USE_HALF_MAGMA)
   cublasHandle_t cublas_handle =
         magma_queue_get_cublas_handle(*(magma_queue_t *)ctx.queue);
   assert(PRIMME_OP_HSCALAR == primme_op_float);
   CHKERR(
         cublasAxpyEx(cublas_handle, n, &alpha, CUDA_R_32F, x, CUDA_R_16F, incx,
               y, CUDA_R_16F, incy, CUDA_R_16F) != CUBLAS_STATUS_SUCCESS);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   return 0;

#else
   XAXPY(n, *(MAGMA_SCALAR *)&alpha, (MAGMA_SCALAR *)x, incx,
         (MAGMA_SCALAR *)y, incy, *(magma_queue_t *)ctx.queue);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   return 0;
#endif
}

/*******************************************************************************
 * Subroutine Num_dot_Sprimme - y'*x
 ******************************************************************************/

TEMPLATE_PLEASE
HSCALAR Num_dot_Sprimme(PRIMME_INT n, SCALAR *x, int incx, SCALAR *y, int incy,
      primme_context ctx) {

   /* Zero dimension matrix may cause problems */
   if (n == 0) return 0;

#if defined(USE_HALFCOMPLEX_MAGMA)
   (void)x;
   (void)incx;
   (void)y;
   (void)incy;
   (void)ctx;
   float nan2[2] = {NAN, NAN};
   return *(HSCALAR*)&nan2;

#elif defined(USE_HALF_MAGMA)
   cublasHandle_t cublas_handle =
         magma_queue_get_cublas_handle(*(magma_queue_t *)ctx.queue);
   assert(PRIMME_OP_HSCALAR == primme_op_float);
   HSCALAR result;
   CHKERR(
         cublasDotEx(cublas_handle, n, x, CUDA_R_16F, incx, y, CUDA_R_16F, incy,
               &result, CUDA_R_32F, CUDA_R_16F) != CUBLAS_STATUS_SUCCESS);
   CHKERRM(cudaSuccess != cudaGetLastError(), (HSCALAR)NAN,
         "Unexpected CUDA error!");
   return result;

#else
   MAGMA_SCALAR result = XDOT(n, (MAGMA_SCALAR *)x, incx, (MAGMA_SCALAR *)y,
         incy, *(magma_queue_t *)ctx.queue);
   CHKERRM(cudaSuccess != cudaGetLastError(), (HSCALAR)NAN,
         "Unexpected CUDA error!");
   return (HSCALAR)*(XSCALAR*)&result;
#endif
}

/*******************************************************************************
 * Subroutine Num_scal_Sprimme - x(0:n*incx-1:incx) *= alpha
 ******************************************************************************/
 
TEMPLATE_PLEASE
int Num_scal_Sprimme(PRIMME_INT n, HSCALAR alpha, SCALAR *x, int incx,
      primme_context ctx) {

   /* Zero dimension matrix may cause problems */
   if (n == 0) return 0;

#if defined(USE_HALFCOMPLEX_MAGMA)
   (void)alpha;
   (void)x;
   (void)incx;
   (void)ctx;
   return PRIMME_FUNCTION_UNAVAILABLE;

#elif defined(USE_HALF_MAGMA)
   cublasHandle_t cublas_handle =
         magma_queue_get_cublas_handle(*(magma_queue_t *)ctx.queue);
   assert(PRIMME_OP_HSCALAR == primme_op_float);
   CHKERR(cublasScalEx(cublas_handle, n, &alpha, CUDA_R_32F, x, CUDA_R_16F,
                incx, CUDA_R_16F) != CUBLAS_STATUS_SUCCESS);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   return 0;

#else
   XSCAL(n, *(MAGMA_SCALAR *)&alpha, (MAGMA_SCALAR *)x, incx,
         *(magma_queue_t *)ctx.queue);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   return 0;
#endif
}

/*******************************************************************************
 * Subroutine Num_larnv_Sprimme - x(0:n*incy-1:incy) = rand(0:n-1)
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_larnv_Sprimme(int idist, PRIMME_INT *iseed, PRIMME_INT length,
      SCALAR *x, primme_context ctx) {

   /* Zero dimension matrix may cause problems */
   if (length == 0) return 0;

   HSCALAR *x_host;
   CHKERR(Num_malloc_SHprimme(length, &x_host, ctx));
   CHKERR(Num_larnv_SHprimme(idist, iseed, length, x_host, ctx));
   magma_setvector(length, sizeof(SCALAR), x_host, 1, x, 1,
         *(magma_queue_t *)ctx.queue);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   CHKERR(Num_free_SHprimme(x_host, ctx));

   return 0;
}

/******************************************************************************
 * Function Num_copy_matrix - Copy the matrix x into y
 *
 * PARAMETERS
 * ---------------------------
 * x           The source matrix
 * m           The number of rows of x
 * n           The number of columns of x
 * ldx         The leading dimension of x
 * y           On output y = x
 * ldy         The leading dimension of y
 *
 * NOTE: x and y *can* overlap
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_copy_matrix_Sprimme(SCALAR *x, PRIMME_INT m, PRIMME_INT n,
      PRIMME_INT ldx, SCALAR *y, PRIMME_INT ldy,
      primme_context ctx) {

#if defined(USE_HALFCOMPLEX_MAGMA)
   (void)x;
   (void)m;
   (void)n;
   (void)ldx;
   (void)y;
   (void)ldy;
   (void)ctx;
   return PRIMME_FUNCTION_UNAVAILABLE;

#else
  assert(m == 0 || n == 0 || (ldx >= m && ldy >= m));

   /* Do nothing if x and y are the same matrix */
   if ((x == y && ldx == ldy) || m == 0 || n == 0)
      return 0;

   magma_copymatrix(m, n, sizeof(SCALAR), (MAGMA_SCALAR *)x, ldx,
         (MAGMA_SCALAR *)y, ldy, *(magma_queue_t *)ctx.queue);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");

   return 0;
#endif /* USE_HALFCOMPLEX_MAGMA */
}

/******************************************************************************
 * Function Num_zero_matrix - Zero the matrix
 *
 * PARAMETERS
 * ---------------------------
 * x           The matrix
 * m           The number of rows of x
 * n           The number of columns of x
 * ldx         The leading dimension of x
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_zero_matrix_Sprimme(SCALAR *x, PRIMME_INT m, PRIMME_INT n,
      PRIMME_INT ldx, primme_context ctx) {

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return 0;

   CHKERR(cudaMemset2D(x, sizeof(SCALAR) * ldx, 0, sizeof(SCALAR) * m, n) !=
          cudaSuccess);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");

   return 0;
} 

/*******************************************************************************
 * Subroutine Num_trsm_hd_Sprimme - b = op(A)\b, where A is triangular
 * NOTE: A is in host and b is in device
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_trsm_hd_Sprimme(const char *side, const char *uplo, const char *transa,
      const char *diag, int m, int n, HSCALAR alpha, HSCALAR *a, int lda,
      SCALAR *b, int ldb, primme_context ctx) {

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return 0;

#if defined(USE_HALFCOMPLEX_MAGMA)
   (void)side;
   (void)uplo;
   (void)transa;
   (void)diag;
   (void)m;
   (void)n;
   (void)alpha;
   (void)a;
   (void)lda;
   (void)b;
   (void)ldb;
   (void)ctx;
   return PRIMME_FUNCTION_UNAVAILABLE;

#elif defined(USE_HALF_MAGMA)
   int mA = (*side == 'R' || *side == 'r') ? n : m;

   /* Create an auxiliary matrix mxm and set the identity matrix */

   HSCALAR *ainv;
   CHKERR(Num_malloc_SHprimme(mA * mA, &ainv, ctx));
   CHKERR(Num_zero_matrix_SHprimme(ainv, mA, mA, mA, ctx));
   int i;
   for (i = 0; i < mA; i++) ainv[mA * i + i] = 1.0;

   /* Compute ainv = inv(A) as A\I */

   CHKERR(Num_trsm_SHprimme(
         "L", uplo, transa, diag, mA, mA, 1.0, a, lda, ainv, mA, ctx));

   /* Compute bainv = b * ainv or ainv * b */

   SCALAR *bainv;
   CHKERR(Num_malloc_Sprimme(m * n, &bainv, ctx));
   CHKERR(Num_zero_matrix_Sprimme(bainv, m, n, m, ctx));
   if (*side == 'R' || *side == 'r') {
      CHKERR(Num_gemm_dhd_Sprimme(
            "N", "N", m, n, n, alpha, b, ldb, ainv, n, 0.0, bainv, m, ctx));
   } else {
      return PRIMME_FUNCTION_UNAVAILABLE;
   }

   /* Copy bainv into b */

   CHKERR(Num_copy_matrix_Sprimme(bainv, m, n, m, b, ldb, ctx));

   CHKERR(Num_free_Sprimme(bainv, ctx));
   CHKERR(Num_free_SHprimme(ainv, ctx));
   return 0;

#else
   int mA = (*side == 'R' || *side == 'r') ? n : m;
   SCALAR *a_dev; /* copy of a on device */
   CHKERR(Num_malloc_Sprimme(mA * mA, &a_dev, ctx));
   magma_setmatrix(mA, mA, sizeof(SCALAR), (MAGMA_SCALAR *)a, lda,
         (MAGMA_SCALAR *)a_dev, mA, *(magma_queue_t *)ctx.queue);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");

   XTRSM(magma_side_const(*side), magma_uplo_const(*uplo),
         magma_trans_const(*transa), magma_diag_const(*diag), m, n,
         *(MAGMA_SCALAR *)&alpha, (MAGMA_SCALAR *)a_dev, mA, (MAGMA_SCALAR *)b,
         ldb, *(magma_queue_t *)ctx.queue);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   CHKERR(Num_free_Sprimme(a_dev, ctx));
   return 0;

#endif
}

/*******************************************************************************
 * Subroutine Num_compute_gramm - Computes the upper part of the Gramm matrix
 *    X' * Y if the result is Hermitian, or the full matrix otherwise, and
 *    do H = X' * Y + alpha * H.
 *
 * Input/Output parameters
 * -----------------------
 * X, Y     The input matrices
 *
 * m, n     Number of rows and columns of X and Y
 *
 * ldX,ldY  Leading dimension of X and Y
 *
 * H        Output matrix storing alpha * H + X' * Y
 *
 * ldH      Leading dimension of H
 *
 * isherm   Whether X' * Y is Hermitian
 * 
 * deep     Maximum number of recursive calls
 * 
 * Return
 * ------
 * error code
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_compute_gramm_Sprimme(SCALAR *X, PRIMME_INT m, int n, int ldX,
      SCALAR *Y, PRIMME_INT ldY, HSCALAR alpha, SCALAR *H, int ldH, int isherm,
      int deep, primme_context ctx) {

   /* Quick exit */

   if (m == 0 || n == 0) return 0;

   // orig-> ---------------- 
   //        |     |        |
   //        |lxl  |        |
   //   l -> |-----| LxL    |
   //        |     |        |
   //   L -> |-----+--------|
   //        |     |  |     |
   //        |     |  | lxl |
   //        ----------------
   //    l == n-L -^  ^
   //    L == n-l ----|

   int L = (n + 1) / 2;
   int l = n / 2;

   // Avoid recursion for non-Hermitian matrices, or matrices smaller than 3 or
   // deep <= 0

   if (!isherm || n <= 2 || deep <= 0) l = 0, L = n;

   CHKERR(Num_compute_gramm_Sprimme(X, m, l, ldX, Y, ldY, alpha, H, ldH,
         deep - 1, 1 /* symmetric */, ctx));
   CHKERR(Num_compute_gramm_Sprimme(&X[ldX * L], m, l, ldX, &Y[ldY * L], ldY,
         alpha, &H[ldH * L + L], ldH, deep - 1, 1 /* symmetric */, ctx));
   CHKERR(Num_gemm_Sprimme("C", "N", L, L, m, 1.0, X, ldX, &Y[ldY * l], ldY,
         alpha, &H[ldH * l], ldH, ctx));

   return 0;
}

/*******************************************************************************
 * Subroutine Num_compute_gramm_ddh - Computes the upper part of the Gramm matrix
 *    X' * Y if the result is Hermitian, or the full matrix otherwise, and
 *    do H = X' * Y + alpha * H.
 *
 * Input/Output parameters
 * -----------------------
 * X, Y     The input matrices
 *
 * m, n     Number of rows and columns of X and Y
 *
 * ldX,ldY  Leading dimension of X and Y
 *
 * H        Output matrix storing alpha * H + X' * Y
 *
 * ldH      Leading dimension of H
 *
 * isherm   Whether X' * Y is Hermitian
 * 
 * Return
 * ------
 * error code
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_compute_gramm_ddh_Sprimme(SCALAR *X, PRIMME_INT m, int n, int ldX,
      SCALAR *Y, PRIMME_INT ldY, HSCALAR alpha, HSCALAR *H, int ldH, int isherm,
      primme_context ctx) {

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return 0;

   if (!isherm) {
      return Num_gemm_ddh_Sprimme(
            "C", "N", n, n, m, 1.0, X, ldX, Y, ldY, alpha, H, ldH, ctx);
   }

#if defined(USE_HALFCOMPLEX_MAGMA)
   return PRIMME_FUNCTION_UNAVAILABLE;

#else
   SCALAR *H_dev; /* copy of H on device */
   CHKERR(Num_malloc_Sprimme(n * n, &H_dev, ctx));
   if (ABS(alpha) == 0.0) {
      CHKERR(Num_zero_matrix_Sprimme(H_dev, n, n, n, ctx));
   } else {
      magma_setmatrix(n, n, sizeof(SCALAR), (MAGMA_SCALAR *)H, ldH,
            (MAGMA_SCALAR *)H_dev, ldH, *(magma_queue_t *)ctx.queue);
      CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
            "Unexpected CUDA error!");
   }
   CHKERR(Num_compute_gramm_Sprimme(
         X, m, n, ldX, Y, ldY, alpha, H_dev, n, 1 /* symmetric */, 4, ctx));
   magma_getmatrix(n, n, sizeof(SCALAR), (MAGMA_SCALAR *)H_dev, n,
         (MAGMA_SCALAR *)H, ldH, *(magma_queue_t *)ctx.queue);
   CHKERRM(cudaSuccess != cudaGetLastError(), PRIMME_UNEXPECTED_FAILURE,
         "Unexpected CUDA error!");
   CHKERR(Num_free_Sprimme(H_dev, ctx));

   return 0;
#endif
}


#undef USE_CUBLAS_HERE

#endif /* USE_MAGMA */

#endif /* SUPPORTED_TYPE */
