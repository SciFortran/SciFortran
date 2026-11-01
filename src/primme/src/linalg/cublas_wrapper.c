/*******************************************************************************
 * Copyright (c) 2020, College of William & Mary
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
 *           calling CUBLAS functions
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../linalg/cublas_wrapper.c"
#endif

#include "numerical.h"

#ifdef SUPPORTED_TYPE

#ifdef USE_CUBLAS

#ifndef CUBLAS_WRAPPER_PRIVATE
#define CUBLAS_WRAPPER_PRIVATE

#ifndef CHECK_TEMPLATE
#    ifdef PRIMME_WITH_CUBLAS
#        include <cuda_runtime.h>
#        include <cublas_v2.h>
#        define GPU_SELECT(X, Y) X
#    elif defined(PRIMME_WITH_HIPBLAS)
#        include <hip/hip_runtime.h>
#        include <hipblas.h>
#        define GPU_SELECT(X, Y) Y
#    else
#        error Unsupported platform
#    endif
#else
#    define GPU_SELECT(X, Y) X
#endif


#define GPU_SYMBOL(X) CONCAT(GPU_SELECT(cuda, hip), X)
#define GPUBLAS_SYMBOL(X) CONCAT(GPU_SELECT(cublas, hipblas), X)
#define GPUBLAS_CONST(X) CONCAT(GPU_SELECT(CUBLAS_, HIPBLAS_), X)

#define CUBLAS_SCALAR                                                          \
   ARITH(, , float, GPU_SELECT(cuComplex, hipblasComplex), double,             \
         GPU_SELECT(cuDoubleComplex, hipblasDoubleComplex), , )

#define XGEMV     CONCAT(GPU_SELECT(cublas,hipblas),ARITH( , , Sgemv , Cgemv , Dgemv , Zgemv , , ))
#define XTRSM     CONCAT(GPU_SELECT(cublas,hipblas),ARITH( , , Strsm , Ctrsm , Dtrsm , Ztrsm , , ))

static int free_fn_dummy (void *p, primme_context ctx) {
   (void)ctx;
   return GPU_SYMBOL(Free)(p) == GPU_SYMBOL(Success) ? 0 : PRIMME_MALLOC_FAILURE;
}

typedef GPU_SELECT(cudaDataType_t, hipblasDatatype_t) gpuDataType;
typedef GPU_SELECT(cublasComputeType_t, hipblasDatatype_t) gpuComputeDataType;

static GPUBLAS_SYMBOL(Operation_t) toCublasOperation(const char trans) {
   switch(trans) {
      case 'n': case 'N': return GPUBLAS_CONST(OP_N);
      case 't': case 'T': return GPUBLAS_CONST(OP_T);
      case 'c': case 'C': return GPUBLAS_CONST(OP_C);
   }
   return (GPUBLAS_SYMBOL(Operation_t))-1;
}

#define CHKERRGPUBLAS(ERRN)                                                    \
   {                                                                           \
      GPUBLAS_SYMBOL(Status_t) ret_;                                           \
      CHKERRM(ret_ = (ERRN),                                                   \
            ret_ == GPUBLAS_CONST(STATUS_SUCCESS)                              \
                  ? 0                                                          \
                  : (ret_ == GPUBLAS_CONST(STATUS_ALLOC_FAILED)                \
                                ? PRIMME_MALLOC_FAILURE                        \
                                : PRIMME_UNEXPECTED_FAILURE),                  \
            "Error in " GPU_SELECT("CUBLAS", "HIPBLAS") " function: %d",       \
            (int)ret_);                                                        \
   }

#define CHKERRCUDA(ERRN)                                                       \
   {                                                                           \
      GPU_SYMBOL(Error_t) ret_;                                                \
      CHKERRM(ret_ = (ERRN),                                                   \
            ret_ == GPU_SYMBOL(Success) ? 0 : PRIMME_UNEXPECTED_FAILURE,       \
            "Error in " GPU_SELECT("CUDA", "HIP") " function: %d", (int)ret_); \
   }

#endif /* CUBLAS_WRAPPER_PRIVATE */

STATIC void toCudaDataType(primme_op_datatype xt, void *r) {
   if (xt == primme_op_default) xt = PRIMME_OP_SCALAR;
   gpuDataType t;
   switch(xt) {
#ifndef USE_COMPLEX
   case primme_op_half:    t = GPU_SELECT(CUDA_R_16F, HIPBLAS_R_16F); break;
   case primme_op_float:   t = GPU_SELECT(CUDA_R_32F, HIPBLAS_R_32F); break;
   case primme_op_double:  t = GPU_SELECT(CUDA_R_64F, HIPBLAS_R_64F); break;
#else
   case primme_op_float:   t = GPU_SELECT(CUDA_C_32F, HIPBLAS_C_32F); break;
   case primme_op_double:  t = GPU_SELECT(CUDA_C_64F, HIPBLAS_C_64F); break;
#endif
   default:                t = (gpuDataType)-1;
   }
   *(gpuDataType*)r = t;
}

STATIC void toCublasComputeDataType(primme_op_datatype xt, void *r) {
#ifdef PRIMME_WITH_HIPBLAS
toCudaDataType(xt, r);
#else
   gpuComputeDataType t;
   if (xt == primme_op_default) xt = PRIMME_OP_SCALAR;
   switch(xt) {
   case primme_op_half:    t = CUBLAS_COMPUTE_16F; break;
   case primme_op_float:   t = CUBLAS_COMPUTE_32F; break;
   case primme_op_double:  t = CUBLAS_COMPUTE_64F; break;
   default:                t = (gpuComputeDataType)-1;
   }
   *(gpuComputeDataType*)r = t;
#endif
}

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
   struct GPU_SELECT(cudaPointerAttributes, hipPointerAttribute_t) ptr_attr;
   if (GPU_SYMBOL(PointerGetAttributes)(&ptr_attr, x) != GPU_SYMBOL(Success))
      return -1;

#if defined(PRIMME_WITH_CUBLAS) && CUDART_VERSION >= 10000
   if (ptr_attr.type == GPU_SYMBOL(MemoryTypeHost)) return -1;
#elif defined(PRIMME_WITH_HIPBLAS) || (defined(PRIMME_WITH_CUBLAS) && CUDART_VERSION < 10000)
   if (!ptr_attr.isManaged && ptr_attr.memoryType == GPU_SYMBOL(MemoryTypeHost))
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

   CHKERRCUDA(GPU_SYMBOL(Malloc)((void**)x, n * sizeof(SCALAR)));

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

   CHKERRCUDA(GPU_SYMBOL(Free(x)));

   return 0;
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

   if (x == y) CHKERR(PRIMME_FUNCTION_UNAVAILABLE);

   /* Call the equivalent real version if needed. It is cheaper to call the */
   /* real variant of GEMM                                                  */

#ifdef USE_COMPLEX
   return Num_copy_Tmatrix_Rprimme(
         x, xt, m * 2, n, ldx * 2, (REAL *)y, ldy * 2, ctx);
#else
   GPUBLAS_SYMBOL(Handle_t) gpu_handle = *(GPUBLAS_SYMBOL(Handle_t) *)ctx.queue;
   XSCALAR zero = 0.0, one = 1.0;
   SET_COMPLEX(zero, 0.0);
   SET_COMPLEX(one, 1.0);
   SCALAR *one_dev;
   CHKERR(Num_malloc_Sprimme(1, &one_dev, ctx));
   HSCALAR one_host;
   SET_COMPLEX(one_host, 1.0);
   CHKERR(Num_set_matrix_Sprimme(&one_host, 1, 1, 1, one_dev, 1, ctx));

   /* Because of calling a BLAS function, output matrix should be initialized */
   CHKERR(Num_zero_matrix_Sprimme(y, m, n, ldy, ctx));

   /* Perform conversion by doing one gemm per column */
   gpuDataType tx; toCudaDataType(xt, &tx);
   gpuDataType ty; toCudaDataType(PRIMME_OP_SCALAR, &ty);
   gpuComputeDataType tc; toCublasComputeDataType(PRIMME_OP_SCALAR, &tc);
   CHKERRGPUBLAS(GPUBLAS_SYMBOL(GemmStridedBatchedEx)(gpu_handle,
         toCublasOperation('n'), toCublasOperation('n'), m, 1, 1, &one, x, tx,
         m, ldx, one_dev, ty, 1, 0, &zero, y, ty, m, ldy, n, tc,
         GPUBLAS_CONST(GEMM_DEFAULT)));

   CHKERR(Num_free_Sprimme(one_dev, ctx));

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

/******************************************************************************
 * Function Num_set_matrix_Sprimme - Copy the matrix x from host into y on device
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
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_set_matrix_Sprimme(HSCALAR *x, PRIMME_INT m, PRIMME_INT n,
      PRIMME_INT ldx, SCALAR *y, PRIMME_INT ldy, primme_context ctx) {

   /* Quick exit */

   if (m == 0 || n == 0) return 0;

   XSCALAR *x0;
   PRIMME_INT ldx0;
   CHKERR(Num_matrix_astype_SHprimme(x, m, n, ldx, PRIMME_OP_HSCALAR,
         (void **)&x0, &ldx0, PRIMME_OP_SCALAR, 1 /* alloc */, 1 /* copy */,
         ctx));
   CHKERRCUDA(GPU_SYMBOL(Memcpy2D)(y, sizeof(SCALAR) * ldy, x0, sizeof(SCALAR) * ldx0,
         sizeof(SCALAR) * m, n, GPU_SYMBOL(MemcpyHostToDevice)));
   CHKERRCUDA(GPU_SYMBOL(DeviceSynchronize)());
   if ((XSCALAR *)x != x0) CHKERR(Num_free_SXprimme(x0, ctx));

   return 0;
}

/******************************************************************************
 * Function Num_get_matrix_Sprimme - Copy the matrix x from device into y on host
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
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_get_matrix_Sprimme(SCALAR *x, PRIMME_INT m, PRIMME_INT n,
      PRIMME_INT ldx, HSCALAR *y, PRIMME_INT ldy, primme_context ctx) {

   /* Quick exit */

   if (m == 0 || n == 0) return 0;

   XSCALAR *y0;
   PRIMME_INT ldy0;
   CHKERR(Num_matrix_astype_SHprimme(y, m, n, ldy, PRIMME_OP_HSCALAR,
         (void **)&y0, &ldy0, PRIMME_OP_SCALAR, 1 /* alloc */,
         0 /* don't copy */, ctx));
   CHKERRCUDA(GPU_SYMBOL(Memcpy2D)(y0, sizeof(SCALAR) * ldy0, x, sizeof(SCALAR) * ldx,
         sizeof(SCALAR) * m, n, GPU_SYMBOL(MemcpyDeviceToHost)));
   CHKERRCUDA(GPU_SYMBOL(DeviceSynchronize)());
   CHKERR(Num_matrix_astype_SHprimme(y0, m, n, ldy0, PRIMME_OP_SCALAR,
         (void **)&y, &ldy, PRIMME_OP_HSCALAR, -1 /* destroy */, 1 /* copy */,
         ctx));
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
         CHKERR(Num_zero_matrix_Sprimme(c, m, n, ldc, ctx));
      }
      else if (beta != (HSCALAR)1.0) {
         int i;
         for (i=0; i<n; i++) {
            CHKERR(Num_scal_Sprimme(m, beta, &c[ldc*i], 1, ctx));
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

#if defined(USE_HALFCOMPLEX_CUBLAS)
   return PRIMME_FUNCTION_UNAVAILABLE;

#else
   GPUBLAS_SYMBOL(Handle_t) gpu_handle = *(GPUBLAS_SYMBOL(Handle_t) *)ctx.queue;
   XSCALAR salpha, sbeta;
   SET_COMPLEX(salpha, alpha);
   SET_COMPLEX(sbeta, beta);
   gpuDataType t; toCudaDataType(PRIMME_OP_SCALAR, &t);
   gpuComputeDataType tc; toCublasComputeDataType(PRIMME_OP_SCALAR, &tc);
   CHKERRGPUBLAS(GPUBLAS_SYMBOL(GemmEx)(gpu_handle, toCublasOperation(*transa),
         toCublasOperation(*transb), m, n, k, (const void *)&salpha,
         (const void *)a, t, lda, (const void *)b, t, ldb, (const void *)&sbeta,
         (CUBLAS_SCALAR *)c, t, ldc, tc, GPUBLAS_CONST(GEMM_DEFAULT)));
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

   int mb = (*transb == 'N' || *transb == 'n') ? k : n;
   int nb = (*transb == 'N' || *transb == 'n') ? n : k;

   SCALAR *b_dev; /* copy of b on device */
   CHKERR(Num_malloc_Sprimme(mb*nb, &b_dev, ctx));
   CHKERR(Num_set_matrix_Sprimme(b, mb, nb, ldb, b_dev, mb, ctx));
   CHKERR(Num_gemm_Sprimme(
         transa, transb, m, n, k, alpha, a, lda, b_dev, mb, beta, c, ldc, ctx));
   CHKERR(Num_free_Sprimme(b_dev, ctx));

   return 0;
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

   SCALAR *c_dev; /* copy of c on device */
   CHKERR(Num_malloc_Sprimme(m*n, &c_dev, ctx));
   if (ABS(beta) != 0) {
      CHKERR(Num_set_matrix_Sprimme(c, m, n, ldc, c_dev, m, ctx));
   } else {
      CHKERR(Num_zero_matrix_Sprimme(c_dev, m, n, m, ctx));
   }
   CHKERR(Num_gemm_Sprimme(
         transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c_dev, m, ctx));
   CHKERR(Num_get_matrix_Sprimme(c_dev, m, n, m, c, ldc, ctx));
   CHKERR(Num_free_Sprimme(c_dev, ctx));

   return 0;
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
         CHKERR(Num_zero_matrix_Sprimme(y, 1, mA, incy, ctx));
      }
      else {
         CHKERR(Num_scal_Sprimme(mA, beta, y, incy, ctx));
      }
      return 0;
   }

#if defined(USE_HALFCOMPLEX_CUBLAS)
   (void)alpha;
   (void)a;
   (void)lda;
   (void)x;
   (void)incx;
   return PRIMME_FUNCTION_UNAVAILABLE;

#elif defined(USE_HALF_CUBLAS)
   /* CUBLAS does not provide gemv for half precision. So we do an            */
   /* equivalent call to gemm. Consider that x and y are given as row vectors */
   /* with leading dimensions incx and incy respectively. We like to do       */
   /* A * x' = y', but gemm does not allow returning a matrix implicit        */
   /* transposed. We do instead x * A' = y. Conjugacy required for complex!   */

   CHKERR(Num_gemm_Sprimme("N", tA ? "N" : "C", 1, mA, nA, alpha, x, incx, a,
         lda, beta, y, incy, ctx));
   return 0;

#else
   GPUBLAS_SYMBOL(Handle_t) gpu_handle = *(GPUBLAS_SYMBOL(Handle_t) *)ctx.queue;
   XSCALAR salpha, sbeta;
   SET_COMPLEX(salpha, alpha);
   SET_COMPLEX(sbeta, beta);
   CHKERRGPUBLAS(XGEMV(gpu_handle, toCublasOperation(*transa), m, n,
         (const CUBLAS_SCALAR *)&salpha, (const CUBLAS_SCALAR *)a, lda,
         (const CUBLAS_SCALAR *)x, incx, (const CUBLAS_SCALAR *)&sbeta,
         (CUBLAS_SCALAR *)y, incy));
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

   int my = (*transa == 'N' || *transa == 'n') ? m : n;
   if (my == 0) return 0;

   SCALAR *y_dev; /* copy of y on device */
   CHKERR(Num_malloc_Sprimme(my, &y_dev, ctx));
   if (ABS(beta) != 0) {
      CHKERR(Num_set_matrix_Sprimme(y, 1, my, incy, y_dev, 1, ctx));
   } else {
      CHKERR(Num_zero_matrix_Sprimme(y_dev, my, 1, my, ctx));
   }
   CHKERR(Num_gemv_Sprimme(
         transa, m, n, alpha, a, lda, x, incx, beta, y_dev, 1, ctx));
   CHKERR(Num_get_matrix_Sprimme(y_dev, 1, my, 1, y, incy, ctx));
   CHKERR(Num_free_Sprimme(y_dev, ctx));

   return 0;
}

/*******************************************************************************
 * Subroutine Num_gemv_dhd_Sprimme - y = alpha*A*x + beta*y, with A size m x n
 * NOTE: A and y are in device and x is in host
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_gemv_dhd_Sprimme(const char *transa, PRIMME_INT m, int n, HSCALAR alpha,
      SCALAR *a, int lda, HSCALAR *x, int incx, HSCALAR beta, SCALAR *y,
      int incy, primme_context ctx) {

   int mx = (*transa == 'N' || *transa == 'n') ? n : m;
   if (mx == 0) return 0;

   SCALAR *x_dev; /* copy of x on device */
   CHKERR(Num_malloc_Sprimme(mx, &x_dev, ctx));
   CHKERR(Num_set_matrix_Sprimme(x, 1, mx, incx, x_dev, 1, ctx));
   CHKERR(Num_gemv_Sprimme(
         transa, m, n, alpha, a, lda, x_dev, 1, beta, y, incy, ctx));
   CHKERR(Num_free_Sprimme(x_dev, ctx));

   return 0;
}

/*******************************************************************************
 * Subroutine Num_axpy_Sprimme - y += alpha*x
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_axpy_Sprimme(PRIMME_INT n, HSCALAR alpha, SCALAR *x, int incx, 
   SCALAR *y, int incy, primme_context ctx) {

   /* Zero dimension matrix may cause problems */
   if (n == 0) return 0;

#if defined(USE_HALFCOMPLEX_CUBLAS)
   (void)alpha;
   (void)x;
   (void)incx;
   (void)y;
   (void)incy;
   (void)ctx;
   return PRIMME_FUNCTION_UNAVAILABLE;

#else
   GPUBLAS_SYMBOL(Handle_t) gpu_handle = *(GPUBLAS_SYMBOL(Handle_t) *)ctx.queue;
   XSCALAR salpha;
   SET_COMPLEX(salpha, alpha);
   gpuDataType t; toCudaDataType(PRIMME_OP_SCALAR, &t);
   CHKERRGPUBLAS(GPUBLAS_SYMBOL(AxpyEx)(gpu_handle, n, (const void *)&salpha, t,
         (const void *)x, t, incx, y, t, incy, t));
   return 0;
#endif
}

/*******************************************************************************
 * Subroutine Num_dot_Sprimme - x'*y
 ******************************************************************************/

TEMPLATE_PLEASE
HSCALAR Num_dot_Sprimme(PRIMME_INT n, SCALAR *x, int incx, SCALAR *y, int incy,
      primme_context ctx) {

   /* Zero dimension matrix may cause problems */
   if (n == 0) return 0;

#if defined(USE_HALFCOMPLEX_CUBLAS)
   (void)x;
   (void)incx;
   (void)y;
   (void)incy;
   (void)ctx;
   float nan2[2] = {NAN, NAN};
   return *(HSCALAR*)&nan2;

#else
   GPUBLAS_SYMBOL(Handle_t) gpu_handle = *(GPUBLAS_SYMBOL(Handle_t) *)ctx.queue;
   HSCALAR result;
   gpuDataType t; toCudaDataType(PRIMME_OP_SCALAR, &t);
   gpuDataType tr; toCudaDataType(PRIMME_OP_HSCALAR, &tr);
   CHKERRGPUBLAS(GPUBLAS_SYMBOL(DotcEx)(
         gpu_handle, n, x, t, incx, y, t, incy, &result, tr, t));
   return result;
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

#if defined(USE_HALFCOMPLEX_CUBLAS)
   (void)alpha;
   (void)x;
   (void)incx;
   (void)ctx;
   return PRIMME_FUNCTION_UNAVAILABLE;

#else
   GPUBLAS_SYMBOL(Handle_t) gpu_handle = *(GPUBLAS_SYMBOL(Handle_t) *)ctx.queue;
   XSCALAR salpha;
   SET_COMPLEX(salpha, alpha);
   gpuDataType t; toCudaDataType(PRIMME_OP_SCALAR, &t);
   CHKERRGPUBLAS(
         GPUBLAS_SYMBOL(ScalEx)(gpu_handle, n, &salpha, t, x, t, incx, t));
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
   CHKERR(Num_set_matrix_Sprimme(x_host, length, 1, length, x, length, ctx));
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

   assert(m == 0 || n == 0 || (ldx >= m && ldy >= m));

   /* Do nothing if x and y are the same matrix */
   if ((x == y && ldx == ldy) || m == 0 || n == 0)
      return 0;

   CHKERRCUDA(GPU_SYMBOL(Memcpy2D)(y, sizeof(SCALAR) * ldy, x, sizeof(SCALAR) * ldx,
         sizeof(SCALAR) * m, n, GPU_SYMBOL(MemcpyDeviceToDevice)));

   return 0;
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

   CHKERRCUDA(GPU_SYMBOL(Memset2D)(x, sizeof(SCALAR) * ldx, 0, sizeof(SCALAR) * m, n));

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

#if defined(USE_HALF_CUBLAS) || defined(USE_HALFCOMPLEX_CUBLAS)
   if (!(*side == 'R' || *side == 'r')) return PRIMME_FUNCTION_UNAVAILABLE;

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

   /* Check that ainv is representable in half precision. If not, normalize */

   HREAL *ainv_scale;
   CHKERR(Num_malloc_RHprimme(n, &ainv_scale, ctx));
   for (i = 0; i < mA; i++) {
      HREAL max_val = 0.0;
      int j;
      for (j = 0; j < mA; j++) max_val = max(max_val, ABS(ainv[i * mA + j]));
      if (max_val >= MACHINE_MAX) {
         ainv_scale[i] = MACHINE_MAX / max_val / 3;
      } else {
         ainv_scale[i] = 1.0;
      }
   }
   CHKERR(Num_scale_matrix_SHprimme(
         ainv, mA, mA, mA, ainv_scale, ainv, mA, ctx));

   /* Compute bainv = b * ainv or ainv * b */

   SCALAR *bainv;
   CHKERR(Num_malloc_Sprimme(m * n, &bainv, ctx));
   CHKERR(Num_zero_matrix_Sprimme(bainv, m, n, m, ctx));
      CHKERR(Num_gemm_dhd_Sprimme(
            "N", "N", m, n, n, alpha, b, ldb, ainv, n, 0.0, bainv, m, ctx));

   /* Copy bainv into b and undo the scale */

   for (i = 0; i < mA; i++) ainv_scale[i] = (HREAL)1.0 / ainv_scale[i];
   CHKERR(Num_scale_matrix_Sprimme(bainv, m, n, m, ainv_scale, b, ldb, ctx));

   CHKERR(Num_free_Sprimme(bainv, ctx));
   CHKERR(Num_free_SHprimme(ainv, ctx));
   CHKERR(Num_free_RHprimme(ainv_scale, ctx));
   return 0;

#else
   int mA = (*side == 'R' || *side == 'r') ? n : m;
   SCALAR *a_dev; /* copy of a on device */
   CHKERR(Num_malloc_Sprimme(mA * mA, &a_dev, ctx));
   CHKERR(Num_set_matrix_Sprimme(a, mA, mA, lda, a_dev, mA, ctx));

   GPUBLAS_SYMBOL(Handle_t) gpu_handle = *(GPUBLAS_SYMBOL(Handle_t) *)ctx.queue;
   XSCALAR salpha;
   SET_COMPLEX(salpha, alpha);
   CHKERRGPUBLAS(XTRSM(gpu_handle,
         (*side == 'R' || *side == 'r') ? GPUBLAS_CONST(SIDE_RIGHT)
                                        : GPUBLAS_CONST(SIDE_LEFT),
         (*uplo == 'U' || *uplo == 'u') ? GPUBLAS_CONST(FILL_MODE_UPPER)
                                        : GPUBLAS_CONST(FILL_MODE_LOWER),
         toCublasOperation(*transa),
         (*diag == 'u' || *diag == 'U') ? GPUBLAS_CONST(DIAG_UNIT)
                                        : GPUBLAS_CONST(DIAG_NON_UNIT),
         m, n, (const CUBLAS_SCALAR *)&salpha, (GPU_SELECT(const,) CUBLAS_SCALAR *)a_dev, mA,
         (CUBLAS_SCALAR *)b, ldb));
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

   SCALAR *H_dev; /* copy of H on device */
   CHKERR(Num_malloc_Sprimme(n * n, &H_dev, ctx));
   if (ABS(alpha) == 0.0) {
      CHKERR(Num_zero_matrix_Sprimme(H_dev, n, n, n, ctx));
   } else {
      CHKERR(Num_set_matrix_Sprimme(H, n, n, ldH, H_dev, n, ctx));
   }
   CHKERR(Num_compute_gramm_Sprimme(
         X, m, n, ldX, Y, ldY, alpha, H_dev, n, 1 /* symmetric */, 4, ctx));
   CHKERR(Num_get_matrix_Sprimme(H_dev, n, n, n, H, ldH, ctx));
   CHKERR(Num_free_Sprimme(H_dev, ctx));

   return 0;
}

#endif /* USE_CUBLAS */

#endif /* SUPPORTED_TYPE */
