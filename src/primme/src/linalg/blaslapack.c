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
 * File: blaslapack.c
 *
 * Purpose - This file contains mostly C wrapper routines for
 *           calling various BLAS and LAPACK FORTRAN routines.
 *
 ******************************************************************************/

#define THIS_FILE "../linalg/blaslapack.c"

#include <string.h>   /* memmove */
#include <assert.h>
#include <math.h>
#include "numerical.h"

#ifdef SUPPORTED_TYPE

#ifdef USE_HOST

#include "blaslapack_private.h"

#ifdef USE_DOUBLE
static int free_fn_dummy (void *p, primme_context ctx) {
   (void)ctx;
   free(p);
   return 0;
}

static int to_blas_int(PRIMME_INT i, PRIMME_BLASINT *out) {
   // If the matrix size is too large for the current BLAS, report an error
   if (i > PRIMME_BLASINT_MAX) return PRIMME_UNEXPECTED_FAILURE;

   *out = (PRIMME_BLASINT)i;
   return 0;
}
#endif

/******************************************************************************
 * Function Num_check_pointer - Return no error code if the pointer is on host.
 * 
 * NOTE: Any check is performed for now.
 *
 * PARAMETERS
 * ---------------------------
 * x           pointer to check
 * 
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_check_pointer_Sprimme(void *x) {
   (void)x;
   return 0;
}

/******************************************************************************
 * Function Num_malloc_Sprimme - Allocate a vector of scalars
 *
 * PARAMETERS
 * ---------------------------
 * n           The number of elements
 * v           returned pointer
 * 
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_malloc_Sprimme(PRIMME_INT n, SCALAR **x, primme_context ctx) {
   (void)ctx;

   /* Quick exit */

   if (n <= 0) {
      *x = NULL;
      return 0;
   }

   /* Allocate memory */

   *x = (SCALAR *)malloc(sizeof(SCALAR) * n);
   if (*x == NULL) return PRIMME_MALLOC_FAILURE;

   /* Register the allocation */

   Mem_keep_frame(ctx);
   Mem_register_alloc(*x, free_fn_dummy, ctx);

   return 0;
}

/******************************************************************************
 * Function Num_free_Sprimme - Free allocated a vector of scalars
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

   free(x);

   return 0;
}

#ifdef USE_DOUBLE

/******************************************************************************
 * Function Num_malloc_iprimme - Allocate a vector of integers
 *
 * PARAMETERS
 * ---------------------------
 * n           The number of elements
 * v           returned pointer
 * 
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_malloc_iprimme(PRIMME_INT n, int **x, primme_context ctx) {
   (void)ctx;

   /* Quick exit */

   if (n <= 0) {
      *x = NULL;
      return 0;
   }

   /* Allocate memory */

   *x = (int *)malloc(sizeof(int) * n);
   if (*x == NULL) return PRIMME_MALLOC_FAILURE;

   /* Register the allocation */

   Mem_keep_frame(ctx);
   Mem_register_alloc(*x, free_fn_dummy, ctx);

   return 0;
}

/******************************************************************************
 * Function Num_free_iprimme - Free allocated a vector of integers
 *
 * PARAMETERS
 * ---------------------------
 * v           allocated pointer
 * 
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_free_iprimme(int *x, primme_context ctx) {
   /* Quick exit */

   if (!x) return 0;

   /* Deregister the allocation */

   Mem_deregister_alloc(x, ctx);

   /* Free pointer */

   free(x);

   return 0;
}



/******************************************************************************
 * Function Num_malloc_iblasprimme - Allocate a vector of integers
 *
 * PARAMETERS
 * ---------------------------
 * n           The number of elements
 * v           returned pointer
 * 
 ******************************************************************************/

static int Num_malloc_iblasprimme(
      PRIMME_INT n, PRIMME_BLASINT **x, primme_context ctx) {
   (void)ctx;

   /* Quick exit */

   if (n <= 0) {
      *x = NULL;
      return 0;
   }

   /* Allocate memory */

   *x = (PRIMME_BLASINT *)malloc(sizeof(PRIMME_BLASINT) * n);
   if (*x == NULL) return PRIMME_MALLOC_FAILURE;

   /* Register the allocation */

   Mem_keep_frame(ctx);
   Mem_register_alloc(*x, free_fn_dummy, ctx);

   return 0;
}

/******************************************************************************
 * Function Num_free_iblasprimme - Free allocated a vector of integers
 *
 * PARAMETERS
 * ---------------------------
 * v           allocated pointer
 * 
 ******************************************************************************/

static int Num_free_iblasprimme(PRIMME_BLASINT *x, primme_context ctx) {
   /* Quick exit */

   if (!x) return 0;

   /* Deregister the allocation */

   Mem_deregister_alloc(x, ctx);

   /* Free pointer */

   free(x);

   return 0;
}
#endif /* USE_DOUBLE */



/*******************************************************************************
 * Subroutine Num_copy_Sprimme - y(0:n*incy-1:incy) = x(0:n*incx-1:incx)
 ******************************************************************************/
 
TEMPLATE_PLEASE
int Num_copy_Sprimme(PRIMME_INT n, SCALAR *x, PRIMME_INT incx, SCALAR *y,
      PRIMME_INT incy, primme_context ctx) {

   /* Quick exit */

   if (x == y && incx == incy) return 0;

#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
   (void)ctx;
   PRIMME_BLASINT ln, lincx, lincy;
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(incx, &lincx));
   CHKERR(to_blas_int(incy, &lincy));

   XCOPY(&ln, (const BLASSCALAR *)x, &lincx, (BLASSCALAR *)y, &lincy);
   return 0;
#else
   return Num_copy_matrix_Sprimme(x, 1, n, incx, y, incy, ctx);
#endif
}

/******************************************************************************
 * Function Num_copy_matrix_Tprimme - Copy the matrix x into y. The types of
 *    x and y can be different.
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

#define Num_copy_matrix_Tprimme(x, xcast, m, n, ldx, y, ldy, ctx)              \
   {                                                                           \
      PRIMME_INT i, j;                                                         \
                                                                               \
      assert((m) == 0 || (n) == 0 || ((ldx) >= (m) && (ldy) >= (m)));          \
                                                                               \
      if ((void *)x != (void *)y)                                              \
         for (i = 0; i < (n); i++)                                             \
            for (j = 0; j < (m); j++)                                          \
               (y)[i * (ldy) + j] = xcast(x)[i * (ldx) + j];                   \
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

#ifdef USE_COMPLEX
   return Num_copy_Tmatrix_Rprimme(
         x, xt, m * 2, n, ldx * 2, (REAL *)y, ldy * 2, ctx);
#else

 
#if defined(USE_HALF) || defined(USE_HALFCOMPLEX)
#  define CAST (float)
#else
#  define CAST
#endif

   switch (xt) {
#  ifdef SUPPORTED_HALF_TYPE
      case primme_op_half:   Num_copy_matrix_Tprimme((PRIMME_HALF*)x, CAST, m, n, ldx, (XSCALAR*)y, ldy, ctx); break;
#  endif
      case primme_op_float:  Num_copy_matrix_Tprimme((float*)      x, CAST, m, n, ldx, (XSCALAR*)y, ldy, ctx); break;
      case primme_op_double: Num_copy_matrix_Tprimme((double*)     x, CAST, m, n, ldx, (XSCALAR*)y, ldy, ctx); break;
      case primme_op_quad:   Num_copy_matrix_Tprimme((PRIMME_QUAD*)x, CAST, m, n, ldx, (XSCALAR*)y, ldy, ctx); break;
      case primme_op_int:    Num_copy_matrix_Tprimme((int*)        x, CAST, m, n, ldx, (XSCALAR*)y, ldy, ctx); break;
      default: CHKERR(PRIMME_FUNCTION_UNAVAILABLE);
   }
#undef CAST

   return 0;
#endif /* USE_COMPLEX */
}

#ifdef USE_DOUBLE
TEMPLATE_PLEASE
int Num_copy_Tmatrix_iprimme(void *x, primme_op_datatype xt, PRIMME_INT m,
      PRIMME_INT n, PRIMME_INT ldx, int *y, PRIMME_INT ldy,
      primme_context ctx) {
   (void)ctx;

#if defined(USE_HALF) || defined(USE_HALFCOMPLEX)
#  define CAST (float)
#else
#  define CAST
#endif

   switch (xt) {
#  ifdef SUPPORTED_HALF_TYPE
      case primme_op_half:   Num_copy_matrix_Tprimme((PRIMME_HALF*)x, CAST, m, n, ldx, y, ldy, ctx); break;
#  endif
      case primme_op_float:  Num_copy_matrix_Tprimme((float*)      x, CAST, m, n, ldx, y, ldy, ctx); break;
      case primme_op_double: Num_copy_matrix_Tprimme((double*)     x, CAST, m, n, ldx, y, ldy, ctx); break;
      case primme_op_quad:   Num_copy_matrix_Tprimme((PRIMME_QUAD*)x, CAST, m, n, ldx, y, ldy, ctx); break;
      case primme_op_int:    Num_copy_matrix_Tprimme((int*)        x, CAST, m, n, ldx, y, ldy, ctx); break;
      default: CHKERR(PRIMME_FUNCTION_UNAVAILABLE);
   }
#undef CAST

   return 0;
}
#endif /* USE_DOUBLE */

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
   (void)ctx;

   PRIMME_INT i, j;

   assert(m == 0 || n == 0 || (ldx >= m && ldy >= m));

   /* Do nothing if x and y are the same matrix */
   if (x == y && ldx == ldy)
      return 0;

   /* Do nothing for zero-size matrices */
   if (m <= 0 || n <= 0)
      return 0;

   /* Copy a contiguous memory region */
   if (ldx == ldy && ldx == m) {
      memmove(y, x, sizeof(SCALAR) * m * n);
   }

   /* Copy the matrix some rows back or forward */
   else if (ldx == ldy && (y > x ? y - x : x - y) < ldx) {
      for (i = 0; i < n; i++)
         memmove(&y[i * ldy], &x[i * ldx], sizeof(SCALAR) * m);
   }

   /* Copy the matrix some columns forward */
   else if (ldx == ldy && y > x && y - x > ldx) {
      for (i = n - 1; i >= 0; i--)
         for (j = 0; j < m; j++)
            y[i * ldy + j] = x[i * ldx + j];
   }

   /* Copy the matrix some columns backward, and other cases */
   else {
      /* TODO: assert x and y don't overlap */
      for (i = 0; i < n; i++) {
         for (j = 0; j < m; j++) {
            y[i * ldy + j] = x[i * ldx + j];
         }
      }
   }

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
  (void)ctx;

   PRIMME_INT i,j;

   for (i=0; i<n; i++)
      for (j=0; j<m; j++)
         SET_ZERO(x[i*ldx+j]);

   return 0;
} 


/*******************************************************************************
 * Subroutine Num_gemm_Sprimme - C = op(A)*op(B), with C size m x n
 ******************************************************************************/

#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
TEMPLATE_PLEASE
int Num_gemm_Sprimme(const char *transa, const char *transb, PRIMME_INT m,
      PRIMME_INT n, PRIMME_INT k, HSCALAR alpha, SCALAR *a, PRIMME_INT lda,
      SCALAR *b, PRIMME_INT ldb, HSCALAR beta, SCALAR *c, PRIMME_INT ldc,
      primme_context ctx) {

   PRIMME_BLASINT lm, ln, lk, llda, lldb, lldc;
   CHKERR(to_blas_int(m, &lm));
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(k, &lk));
   CHKERR(to_blas_int(lda, &llda));
   CHKERR(to_blas_int(ldb, &lldb));
   CHKERR(to_blas_int(ldc, &lldc));

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return 0;

   /* Quick exit */
   if (k == 0 || ABS(alpha) == 0.0) {
      if (ABS(beta) == 0.0) {
         Num_zero_matrix_Sprimme(c, m, n, ldc, ctx);
      }
      else if (beta != (HSCALAR)1.0) {
         PRIMME_INT i;
         for (i=0; i<n; i++) {
            Num_scal_Sprimme(m, beta, &c[ldc*i], 1, ctx);
         }
      }
      return 0;
   }
   if (n == 1
#ifdef USE_COMPLEX
         && *transb != 'c' && *transb != 'C'
#endif
   ) {
      PRIMME_INT mA; PRIMME_INT nA;
      if (*transa == 'n' || *transa == 'N') mA = m, nA = k;
      else mA = k, nA = m;
      PRIMME_INT incb = ((*transb == 'n' || *transb == 'N') ? 1 : ldb);
      return Num_gemv_Sprimme(
            transa, mA, nA, alpha, a, lda, b, incb, beta, c, 1, ctx);
   }

   XGEMM(transa, transb, &lm, &ln, &lk, (const BLASSCALAR *)&alpha,
         (const BLASSCALAR *)a, &llda, (const BLASSCALAR *)b, &lldb,
         (const BLASSCALAR *)&beta, (BLASSCALAR *)c,
         &lldc STRING_LEN1 STRING_LEN1);

   return 0;
}
#endif /* (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF) */

TEMPLATE_PLEASE
int Num_gemm_dhd_Sprimme(const char *transa, const char *transb, PRIMME_INT m,
      PRIMME_INT n, PRIMME_INT k, HSCALAR alpha, SCALAR *a, PRIMME_INT lda,
      HSCALAR *b, PRIMME_INT ldb, HSCALAR beta, SCALAR *c, PRIMME_INT ldc,
      primme_context ctx) {

   /* Quick exit */

   if (m == 0 || n == 0 ||
         ((k == 0 || ABS(alpha) == 0.0) && beta == (HSCALAR)1.0))
      return 0;

   /* Cast the matrices a and c to HSCALAR */

   HSCALAR *af = NULL, *cf = NULL;
   PRIMME_INT ldaf, ldcf, ldc0 = ldc;
   PRIMME_INT ma = (*transa == 'N' || *transa == 'n') ? m : k;
   PRIMME_INT na = (*transa == 'N' || *transa == 'n') ? k : m;

   CHKERR(Num_matrix_astype_Sprimme(a, ma, na, lda, PRIMME_OP_SCALAR,
         (void **)&af, &ldaf, PRIMME_OP_HSCALAR, 1 /* alloc */, 1 /* copy */,
         ctx));
   CHKERR(Num_matrix_astype_Sprimme(c, m, n, ldc, PRIMME_OP_SCALAR,
         (void **)&cf, &ldcf, PRIMME_OP_HSCALAR, 1 /* alloc */, 1 /* copy */,
         ctx));

   /* Call the kernel */

   CHKERR(Num_gemm_SHprimme(transa, transb, m, n, k, alpha, af, ldaf, b, ldb,
         beta, cf, ldcf, ctx));

   /* Copy back c and destroy the cast matrices */

   if (a != (SCALAR*)af) CHKERR(Num_free_SHprimme(af, ctx));
   CHKERR(Num_matrix_astype_Sprimme(cf, m, n, ldcf, PRIMME_OP_HSCALAR,
         (void **)&c, &ldc0, PRIMME_OP_SCALAR, -1 /* destroy */, 1 /* copy */,
         ctx));

   return 0;
}

TEMPLATE_PLEASE
int Num_gemm_ddh_Sprimme(const char *transa, const char *transb, PRIMME_INT m,
      PRIMME_INT n, PRIMME_INT k, HSCALAR alpha, SCALAR *a, PRIMME_INT lda,
      SCALAR *b, PRIMME_INT ldb, HSCALAR beta, HSCALAR *c, PRIMME_INT ldc,
      primme_context ctx) {

   /* Quick exit */

   if (m == 0 || n == 0 ||
         ((k == 0 || ABS(alpha) == 0.0) && beta == (HSCALAR)1.0))
      return 0;

   /* If input matrices are going to be cast and the operation is the inner */
   /* product between, that is, C = A'*B, then stream the operation */

   PRIMME_INT K = k;
   if (PRIMME_OP_SCALAR != PRIMME_OP_HSCALAR &&
         (*transa == 'C' || *transa == 'C') &&
         (*transb == 'N' || *transb == 'n') && k > PRIMME_BLOCK_SIZE) {
      K = PRIMME_BLOCK_SIZE;
   }

   /* Cast the matrices a and b to HSCALAR */

   HSCALAR *af = NULL, *bf = NULL;
   PRIMME_INT ldaf, ldbf;
   PRIMME_INT na = (*transa == 'N' || *transa == 'n') ? k : m;
   PRIMME_INT nb = (*transb == 'N' || *transb == 'n') ? n : k;

   PRIMME_INT i;
   for (i=0; i<k; i+=K, K=min(K,k-i)) {
      CHKERR(Num_matrix_astype_Sprimme(&a[i], K, na, lda, PRIMME_OP_SCALAR,
            (void **)&af, &ldaf, PRIMME_OP_HSCALAR,
            i == 0 /* alloc the first time */, 1 /* copy */, ctx));
      CHKERR(Num_matrix_astype_Sprimme(&b[i], K, nb, ldb, PRIMME_OP_SCALAR,
            (void **)&bf, &ldbf, PRIMME_OP_HSCALAR,
            i == 0 /* alloc the first time */, 1 /* copy */, ctx));

      /* Call the kernel */

      CHKERR(Num_gemm_SHprimme(transa, transb, m, n, K, alpha, af, ldaf, bf,
            ldbf, beta, c, ldc, ctx));

      /* Update beta */

      beta = (HSCALAR)1.0;
   }

   if (a != (SCALAR *)af) CHKERR(Num_free_SHprimme(af, ctx));
   if (b != (SCALAR *)bf) CHKERR(Num_free_SHprimme(bf, ctx));

   return 0;
}

/*******************************************************************************
 * Subroutine Num_hemm_Sprimme - C = A*B or B*A where A is Hermitian,
 *    where C size m x n.
 ******************************************************************************/

#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
TEMPLATE_PLEASE
int Num_hemm_Sprimme(const char *side, const char *uplo, PRIMME_INT m,
      PRIMME_INT n, HSCALAR alpha, SCALAR *a, PRIMME_INT lda, SCALAR *b,
      PRIMME_INT ldb, HSCALAR beta, SCALAR *c, PRIMME_INT ldc,
      primme_context ctx) {

   (void)ctx;
   PRIMME_BLASINT lm, ln, llda, lldb, lldc;
   CHKERR(to_blas_int(m, &lm));
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(lda, &llda));
   CHKERR(to_blas_int(ldb, &lldb));
   CHKERR(to_blas_int(ldc, &lldc));

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return 0;

   XHEMM(side, uplo, &lm, &ln, (const BLASSCALAR *)&alpha,
         (const BLASSCALAR *)a, &llda, (const BLASSCALAR *)b, &lldb,
         (const BLASSCALAR *)&beta, (BLASSCALAR *)c,
         &lldc STRING_LEN1 STRING_LEN1);

   return 0;
}
#endif

/*******************************************************************************
 * Subroutine Num_trmm_Sprimme - B = A*B or B*A where A is triangular,
 *    with B size m x n.
 ******************************************************************************/

#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
TEMPLATE_PLEASE
int Num_trmm_Sprimme(const char *side, const char *uplo, const char *transa,
      const char *diag, PRIMME_INT m, PRIMME_INT n, HSCALAR alpha, SCALAR *a,
      PRIMME_INT lda, SCALAR *b, PRIMME_INT ldb, primme_context ctx) {

   (void)ctx;
   PRIMME_BLASINT lm, ln, llda, lldb;
   CHKERR(to_blas_int(m, &lm));
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(lda, &llda));
   CHKERR(to_blas_int(ldb, &lldb));

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return 0;

   XTRMM(side, uplo, transa, diag, &lm, &ln, (const BLASSCALAR *)&alpha,
         (const BLASSCALAR *)a, &llda, (BLASSCALAR *)b,
         &lldb STRING_LEN1 STRING_LEN1 STRING_LEN1 STRING_LEN1);

   return 0;
}
#endif /* (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF) */

/*******************************************************************************
 * Subroutine Num_gemv_Sprimme - y = alpha*A*x + beta*y, with A size m x n
 ******************************************************************************/

#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
TEMPLATE_PLEASE
int Num_gemv_Sprimme(const char *transa, PRIMME_INT m, PRIMME_INT n,
      HSCALAR alpha, SCALAR *a, PRIMME_INT lda, SCALAR *x, PRIMME_INT incx,
      HSCALAR beta, SCALAR *y, PRIMME_INT incy, primme_context ctx) {

   PRIMME_BLASINT lm, ln, llda, lincx, lincy;
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(lda, &llda));
   CHKERR(to_blas_int(incx, &lincx));
   CHKERR(to_blas_int(incy, &lincy));

   /* Zero dimension matrix may cause problems */
   int tA = (*transa != 'n' && *transa != 'N' ? 1 : 0);
   PRIMME_INT mA = tA ? n : m, nA = tA ? m : n;
   if (mA == 0) return 0;

   /* Quick exit */
   if (nA == 0 || ABS(alpha) == 0.0) {
      if (ABS(beta) == 0.0) {
         Num_zero_matrix_Sprimme(y, 1, mA, incy, ctx);
      }
      else {
         Num_scal_Sprimme(mA, beta, y, incy, ctx);
      }
      return 0;
   }
   if (mA == 1) {
      if (ABS(beta) != 0.0)
         y[0] *= beta;
      else
         y[0] = 0.0;
      if (tA) {
         y[0] += Num_dot_Sprimme(nA, a, 1, x, incx, ctx) * alpha;
      } else {
         y[0] += Num_dot_Sprimme(nA, a, lda, x, incx, ctx) * alpha;
      }
      return 0;
   }

   while(m > 0) {
      lm = (PRIMME_BLASINT)min(m, PRIMME_BLASINT_MAX-1);
      XGEMV(transa, &lm, &ln, (const BLASSCALAR *)&alpha, (const BLASSCALAR *)a,
            &llda, (const BLASSCALAR *)x, &lincx, (const BLASSCALAR *)&beta,
            (BLASSCALAR *)y, &lincy STRING_LEN1);
      m -= (PRIMME_INT)lm;
      a += lm;
      if (!tA) {
         y += lm;
      }
      else {
         x += lm;
         beta = 1.0;
      }
   }

   return 0;
}
#endif /* (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF) */

/*******************************************************************************
 * Subroutine Num_gemv_ddh_Sprimme - y = alpha*A*x + beta*y, with A size m x n
 * NOTE: A and x are in device and y is in host
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_gemv_ddh_Sprimme(const char *transa, PRIMME_INT m, PRIMME_INT n,
      HSCALAR alpha, SCALAR *a, PRIMME_INT lda, SCALAR *x, PRIMME_INT incx,
      HSCALAR beta, HSCALAR *y, PRIMME_INT incy, primme_context ctx) {

   /* Cast the matrices a and x to HSCALAR */

   HSCALAR *af = NULL, *xf = NULL;
   PRIMME_INT ldaf, incxf;
   PRIMME_INT mx = (*transa == 'N' || *transa == 'n') ? n : m;

   CHKERR(Num_matrix_astype_Sprimme(a, m, n, lda, PRIMME_OP_SCALAR,
         (void **)&af, &ldaf, PRIMME_OP_HSCALAR, 1 /* alloc */, 1 /* copy */,
         ctx));
   CHKERR(Num_matrix_astype_Sprimme(x, 1, mx, incx, PRIMME_OP_SCALAR,
         (void **)&xf, &incxf, PRIMME_OP_HSCALAR, 1 /* alloc */, 1 /* copy */,
         ctx));

   /* Call the kernel */

   CHKERR(Num_gemv_SHprimme(
         transa, m, n, alpha, af, m, xf, 1, beta, y, incy, ctx));

   if (a != (SCALAR*)af) CHKERR(Num_free_SHprimme(af, ctx));
   if (x != (SCALAR*)xf) CHKERR(Num_free_SHprimme(xf, ctx));

   return 0;
}

/*******************************************************************************
 * Subroutine Num_gemv_dhd_Sprimme - y = alpha*A*x + beta*y, with A size m x n
 * NOTE: A and y are in device and x is in host
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_gemv_dhd_Sprimme(const char *transa, PRIMME_INT m, PRIMME_INT n,
      HSCALAR alpha, SCALAR *a, PRIMME_INT lda, HSCALAR *x, PRIMME_INT incx,
      HSCALAR beta, SCALAR *y, PRIMME_INT incy, primme_context ctx) {

   /* Cast the matrices a and y to HSCALAR */

   HSCALAR *af = NULL, *yf = NULL;
   PRIMME_INT ldaf, incyf, incy0 = incy;
   PRIMME_INT my = (*transa == 'N' || *transa == 'n') ? m : n;

   CHKERR(Num_matrix_astype_Sprimme(a, m, n, lda, PRIMME_OP_SCALAR,
         (void **)&af, &ldaf, PRIMME_OP_HSCALAR, 1 /* alloc */, 1 /* copy */,
         ctx));
   CHKERR(Num_matrix_astype_Sprimme(y, 1, my, incy, PRIMME_OP_SCALAR,
         (void **)&yf, &incyf, PRIMME_OP_HSCALAR, 1 /* alloc */, 1 /* copy */,
         ctx));

   CHKERR(Num_gemv_SHprimme(
         transa, m, n, alpha, af, m, x, incx, beta, yf, incyf, ctx));

   /* Copy back y and destroy cast matrices */

   if (a != (SCALAR*)af) CHKERR(Num_free_SHprimme(af, ctx));
   CHKERR(Num_matrix_astype_Sprimme(yf, 1, my, incyf, PRIMME_OP_HSCALAR,
         (void **)&y, &incy0, PRIMME_OP_SCALAR, -1 /* destroy */, 1 /* copy */,
         ctx));

   return 0;
}

/*******************************************************************************
 * Subroutine Num_axpy_Sprimme - y += alpha*x
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_axpy_Sprimme(PRIMME_INT n, HSCALAR alpha, SCALAR *x, PRIMME_INT incx,
      SCALAR *y, PRIMME_INT incy, primme_context ctx) {

   (void)ctx;
#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
   PRIMME_BLASINT ln, lincx, lincy;
   CHKERR(to_blas_int(incx, &lincx));
   CHKERR(to_blas_int(incy, &lincy));

   while(n > 0) {
      ln = (PRIMME_BLASINT)min(n, PRIMME_BLASINT_MAX-1);
      XAXPY(&ln, (const BLASSCALAR *)&alpha, (const BLASSCALAR *)x, &lincx,
            (BLASSCALAR *)y, &lincy);
      n -= (PRIMME_INT)ln;
      x += ln;
      y += ln;
   }
#else
   PRIMME_INT i;
   for (i = 0; i < n; i++)
      PLUS_EQUAL(y[incy * i], alpha * TO_COMPLEX(x[incx * i]));
#endif

   return 0;
}

/*******************************************************************************
 * Subroutine Num_dot_Sprimme - x'*y
 ******************************************************************************/

TEMPLATE_PLEASE
HSCALAR Num_dot_Sprimme(PRIMME_INT n, SCALAR *x, PRIMME_INT incx, SCALAR *y,
      PRIMME_INT incy, primme_context ctx) {
   (void)ctx;
/* NOTE: vecLib doesn't follow BLAS reference for sdot */
#if defined(USE_COMPLEX) || (defined(USE_FLOAT) && (defined(__APPLE__) || defined(__MACH__))) || (defined(USE_HALF) && !defined(BLASLAPACK_WITH_HALF))
/* ---- Explicit implementation of the zdotc() --- */
   PRIMME_INT i;
   HSCALAR zdotc = 0.0;
   if (n <= 0) return(zdotc);
   if (incx == 1 && incy == 1) {
      for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
         zdotc += CONJ(TO_COMPLEX(x[i])) * TO_COMPLEX(y[i]);
      }
   }
   else {
      for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
         zdotc += CONJ(TO_COMPLEX(x[i*incx])) * TO_COMPLEX(y[i*incy]);
      }
   }
   return zdotc;
/* -- end of explicit implementation of the zdotc() - */
#else
   PRIMME_BLASINT ln, lincx, lincy;
   CHKERR(to_blas_int(incx, &lincx));
   CHKERR(to_blas_int(incy, &lincy));
   HSCALAR r = 0.0;

   while(n > 0) {
      ln = (PRIMME_BLASINT)min(n, PRIMME_BLASINT_MAX-1);
      r += XDOT(&ln, (const BLASSCALAR *)x, &lincx, (BLASSCALAR *)y, &lincy);
      n -= (PRIMME_INT)ln;
      x += ln;
      y += ln;
   }

   return r;
#endif
}

/*******************************************************************************
 * Subroutine Num_larnv_Sprimme - x(0:n*incy-1:incy) = rand(0:n-1)
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_larnv_Sprimme(int idist, PRIMME_INT *iseed, PRIMME_INT length,
      SCALAR *x, primme_context ctx) {

   (void)ctx;

#ifdef USE_COMPLEX
   /* Lapack's R core library doesn't have zlarnv. The functionality is */
   /* replaced by calling the REAL version doubling the length.         */

   assert(idist < 4); /* complex distributions are not supported */
   return Num_larnv_Rprimme(idist, iseed, length*2, (REAL*)x, ctx);

#elif (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
   PRIMME_BLASINT lidist;
   CHKERR(to_blas_int(idist, &lidist));
   PRIMME_BLASINT llength;
   PRIMME_BLASINT temp[4];
   PRIMME_BLASINT *liseed;
   int i;

   if (sizeof(PRIMME_INT) == sizeof(PRIMME_BLASINT)) {
      liseed = (PRIMME_BLASINT*)iseed; /* cast avoid compiler warning */
   } else {
      liseed = temp;
      for(i=0; i<4; i++)
         liseed[i] = (PRIMME_BLASINT)iseed[i];
   }

   while(length > 0) {
      llength = (PRIMME_BLASINT)min(length, PRIMME_BLASINT_MAX-1);
      XLARNV(&lidist, liseed, &llength, x);
      length -= (PRIMME_INT)llength;
      x += llength;
   }

   if (sizeof(PRIMME_INT) != sizeof(PRIMME_BLASINT))
      for(i=0; i<4; i++)
         iseed[i] = (PRIMME_INT)liseed[i];

   return 0;
#else
   /* Call the HSCALAR kernel and copy back the result */
   HSCALAR *xf = NULL;
   CHKERR(Num_malloc_SHprimme(length, &xf, ctx));
   CHKERR(Num_larnv_SHprimme(idist, iseed, length, xf, ctx));
   CHKERR(Num_matrix_astype_SHprimme(xf, 1, length, 1, PRIMME_OP_HSCALAR,
         (void **)&x, NULL, PRIMME_OP_SCALAR, -1 /* destroy */, 1 /* copy */,
         ctx));
   return 0;
#endif
}

/*******************************************************************************
 * Subroutine Num_scal_Sprimme - x(0:n*incx-1:incx) *= alpha
 ******************************************************************************/
 
TEMPLATE_PLEASE
int Num_scal_Sprimme(PRIMME_INT n, HSCALAR alpha, SCALAR *x, PRIMME_INT incx,
      primme_context ctx) {

   (void)ctx;
#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
   PRIMME_BLASINT ln, lincx;
   CHKERR(to_blas_int(incx, &lincx));

   while(n > 0) {
      ln = (PRIMME_BLASINT)min(n, PRIMME_BLASINT_MAX-1);
      XSCAL(&ln, (const BLASSCALAR *)&alpha, (BLASSCALAR *)x, &lincx);
      n -= (PRIMME_INT)ln;
      x += ln;
   }
#else
   PRIMME_INT i;
   for (i = 0; i < n; i++) MULT_EQUAL(x[incx * i], alpha);
#endif

   return 0;
}

/*******************************************************************************
 * Subroutines for dense eigenvalue decomposition
 * NOTE: xheevx is used instead of xheev because xheev is not in ESSL
 ******************************************************************************/
 
#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
TEMPLATE_PLEASE
int Num_heev_Sprimme(const char *jobz, const char *uplo, PRIMME_INT n,
      SCALAR *a, PRIMME_INT lda, REAL *w, primme_context ctx) {
#ifndef USE_XHEEV

   PRIMME_BLASINT ln, llda;
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(lda, &llda));
   PRIMME_BLASINT lldwork = 0;
   PRIMME_BLASINT linfo = 0;
   SCALAR *z;
   REAL abstol=0.0;
#  ifdef USE_COMPLEX
   REAL *rwork;
#  endif
   PRIMME_BLASINT *iwork, *ifail;
   REAL   dummyr=0;
   PRIMME_BLASINT dummyi=0;

   /* Zero dimension matrix may cause problems */
   if (n == 0) return 0;

   /* Allocate arrays */

   CHKERR(Num_malloc_Sprimme(n*n, &z, ctx)); 
#  ifdef USE_COMPLEX
   CHKERR(Num_malloc_Rprimme(7*n, &rwork, ctx));
#  endif
   CHKERR(Num_malloc_iblasprimme(5*n, &iwork, ctx));
   CHKERR(Num_malloc_iblasprimme(n, &ifail, ctx));

   /* Call to know the optimal workspace */

   lldwork = -1;
   SCALAR lwork0 = 0;
   XHEEVX(jobz, "A", uplo, &ln, a, &llda, &dummyr, &dummyr,
         &dummyi, &dummyi, &abstol, &dummyi, w, z, &ln, &lwork0, &lldwork,
#  ifdef USE_COMPLEX
         &dummyr,
#  endif
         iwork, ifail, &linfo);
   lldwork = REAL_PART(lwork0);
   // ATLAS LAPACK and MacOS LAPACK may suggest a wrong value to lwork for n=1
   if (lldwork < 2 * ln) lldwork = 2 * ln;

   if (linfo == 0) {
      SCALAR *work = NULL;
      CHKERR(Num_malloc_Sprimme(lldwork, &work, ctx));
      XHEEVX(jobz, "A", uplo, &ln, a, &llda, &dummyr, &dummyr,
            &dummyi, &dummyi, &abstol, &dummyi, w, z, &ln, work, &lldwork,
#  ifdef USE_COMPLEX
            rwork,
#  endif
            iwork, ifail, &linfo);
      CHKERR(Num_free_Sprimme(work, ctx));
   }

   /* Copy z to a */
   Num_copy_matrix_Sprimme(z, n, n, n, a, lda, ctx);

   CHKERR(Num_free_Sprimme(z, ctx)); 
#  ifdef USE_COMPLEX
   CHKERR(Num_free_Rprimme(rwork, ctx));
#  endif
   CHKERR(Num_free_iblasprimme(iwork, ctx));
   CHKERR(Num_free_iblasprimme(ifail, ctx));

   CHKERRM(linfo != 0, PRIMME_LAPACK_FAILURE, "Error in xheev with info %d",
          (int)linfo);
   return 0;

#else /* USE_XHEEV */

   PRIMME_BLASINT ln, llda;
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(lda, &llda));
   PRIMME_BLASINT lldwork = 0;
   PRIMME_BLASINT linfo = 0;
#  ifdef USE_COMPLEX
   REAL *rwork;
#  endif

   /* Zero dimension matrix may cause problems */
   if (n == 0) return 0;

#  ifdef USE_COMPLEX
   CHKERR(Num_malloc_Rprimme(3*n, &rwork, ctx));
#  endif

   /* Call to know the optimal workspace */

   lldwork = -1;
   SCALAR lwork0 = 0;
   XHEEV(jobz, uplo, &ln, a, &llda, w, &lwork0, &lldwork,
#     ifdef USE_COMPLEX
         rwork,
#     endif
         &linfo); 
   lldwork = REAL_PART(lwork0);

   if (linfo == 0) {
      SCALAR *work = NULL;
      CHKERR(Num_malloc_Sprimme(lldwork, &work, ctx));
      XHEEV(jobz, uplo, &ln, a, &llda, w, work, &lldwork,
#     ifdef USE_COMPLEX
            rwork,
#     endif
            &linfo); 
      CHKERR(Num_free_Sprimme(work, ctx));
   }

#  ifdef USE_COMPLEX
   CHKERR(Num_free_Rprimme(rwork, ctx));
#  endif
   
   CHKERRM(linfo != 0, PRIMME_LAPACK_FAILURE, "Error in xheev with info %d",
          (int)linfo);
   return 0;

#endif
}
#endif /* (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF) */


/*******************************************************************************
 * Subroutines for dense generalize eigenvalue decomposition
 * NOTE: xhegvx is used instead of xhegv because xhegv is not in ESSL
 ******************************************************************************/
 
#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
TEMPLATE_PLEASE
int Num_hegv_Sprimme(const char *jobz, const char *uplo, PRIMME_INT n,
      SCALAR *a, PRIMME_INT lda, SCALAR *b0, PRIMME_INT ldb0, REAL *w,
      primme_context ctx) {

   /* Call heev if b is null */
   if (b0 == NULL) {
      return Num_heev_Sprimme(jobz, uplo, n, a, lda, w, ctx);
   }

#ifndef USE_XHEGV

   PRIMME_BLASINT ln, llda;
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(lda, &llda));
   PRIMME_BLASINT lldwork = 0;
   PRIMME_BLASINT linfo = 0;
   PRIMME_BLASINT ONE = 1;
   SCALAR *z, *b;
   REAL abstol=0.0;
#  ifdef USE_COMPLEX
   REAL *rwork;
#  endif
   PRIMME_BLASINT *iwork, *ifail;
   REAL   dummyr=0;
   PRIMME_BLASINT dummyi=0;

   /* Zero dimension matrix may cause problems */
   if (n == 0) return 0;

   /* Allocate arrays */

   CHKERR(Num_malloc_Sprimme(n*n, &z, ctx)); 
   CHKERR(Num_malloc_Sprimme(n*n, &b, ctx)); 
#  ifdef USE_COMPLEX
   CHKERR(Num_malloc_Rprimme(7*n, &rwork, ctx));
#  endif
   CHKERR(Num_malloc_iblasprimme(5*n, &iwork, ctx));
   CHKERR(Num_malloc_iblasprimme(n, &ifail, ctx));

   CHKERR(Num_copy_trimatrix_Sprimme(b0, n, n, ldb0,
         *uplo == 'U' || *uplo == 'u' ? 0 : 1, 0, b, n,
         0 /*not to zero rest of the matrix */));

   /* Call to know the optimal workspace */

   lldwork = -1;
   SCALAR lwork0 = 0;
   XHEGVX(&ONE, jobz, "A", uplo, &ln, a, &llda, b, &ln, &dummyr, &dummyr,
         &dummyi, &dummyi, &abstol, &dummyi, w, z, &ln, &lwork0, &lldwork,
#  ifdef USE_COMPLEX
         rwork,
#  endif
         iwork, ifail, &linfo);
   lldwork = REAL_PART(lwork0);

   if (linfo == 0) {
      SCALAR *work = NULL;
      CHKERR(Num_malloc_Sprimme(lldwork, &work, ctx));
      XHEGVX(&ONE, jobz, "A", uplo, &ln, a, &llda, b, &ln, &dummyr, &dummyr,
            &dummyi, &dummyi, &abstol, &dummyi, w, z, &ln, work, &lldwork,
#  ifdef USE_COMPLEX
            rwork,
#  endif
            iwork, ifail, &linfo);
      CHKERR(Num_free_Sprimme(work, ctx));
   }

   /* Copy z to a */
   Num_copy_matrix_Sprimme(z, n, n, n, a, lda, ctx);

   CHKERR(Num_free_Sprimme(z, ctx)); 
   CHKERR(Num_free_Sprimme(b, ctx)); 
#  ifdef USE_COMPLEX
   CHKERR(Num_free_Rprimme(rwork, ctx));
#  endif
   CHKERR(Num_free_iblasprimme(iwork, ctx));
   CHKERR(Num_free_iblasprimme(ifail, ctx));

   CHKERRM(linfo != 0, PRIMME_LAPACK_FAILURE, "Error in xhegvx with info %d",
          (int)linfo);
   return 0;


#else /* USE_XHEGV */

   PRIMME_BLASINT ln, llda;
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(lda, &llda));
   PRIMME_BLASINT lldwork = 0;
   PRIMME_BLASINT linfo = 0;
   PRIMME_BLASINT ONE = 1;
   SCALAR *b;
#  ifdef USE_COMPLEX
   REAL *rwork;
#  endif

   /* Zero dimension matrix may cause problems */
   if (n == 0) return 0;

   CHKERR(Num_malloc_Sprimme(n*n, &b, ctx)); 
#  ifdef USE_COMPLEX
   CHKERR(Num_malloc_Rprimme(3*n, &rwork, ctx));
#  endif

   Num_copy_matrix_Sprimme(b0, n, n, ldb0, b, n, ctx);

   /* Call to know the optimal workspace */

   lldwork = -1;
   SCALAR lwork0 = 0;
   XHEGV(&ONE, jobz, uplo, &ln, a, &llda, b, &ln, w, &lwork0, &lldwork,
#  ifdef USE_COMPLEX
         rwork,
#  endif
         &linfo);
   lldwork = REAL_PART(lwork0);

   if (linfo == 0) {
      SCALAR *work = NULL;
      CHKERR(Num_malloc_Sprimme(lldwork, &work, ctx));
      XHEGV(&ONE, jobz, uplo, &ln, a, &llda, b, &ln, w, work, &lldwork,
#  ifdef USE_COMPLEX
            rwork,
#  endif
            &linfo);
      CHKERR(Num_free_Sprimme(work, ctx));
   }

   CHKERR(Num_free_Sprimme(b, ctx)); 
#  ifdef USE_COMPLEX
   CHKERR(Num_free_Rprimme(rwork, ctx));
#  endif
   return 0;
 
#endif
}
#endif /* (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF) */

/*******************************************************************************
 * Subroutines for dense Schur decomposition
 * NOTE: only for complex matrices
 ******************************************************************************/
 
#if defined(USE_COMPLEX) && ((!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF))

TEMPLATE_PLEASE
int Num_gees_Sprimme(const char *jobvs, PRIMME_INT n, SCALAR *a, PRIMME_INT lda,
      SCALAR *w, SCALAR *vs, PRIMME_INT ldvs, primme_context ctx) {

   PRIMME_BLASINT ln, llda, lldvs;
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(lda, &llda));
   CHKERR(to_blas_int(ldvs, &lldvs));
   PRIMME_BLASINT lldwork = 0;
   PRIMME_BLASINT linfo = 0;
   PRIMME_BLASINT lsdim = 0;
   REAL *rwork;
   PRIMME_BLASINT dummyi=0;

   /* Zero dimension matrix may cause problems */
   if (n == 0) return 0;

   /* Allocate arrays */

   CHKERR(Num_malloc_Rprimme(n, &rwork, ctx));

   /* Call to know the optimal workspace */

   lldwork = -1;
   SCALAR lwork0 = 0;
   XGEES(jobvs, "N", NULL, &ln, a, &llda, &lsdim, w, vs, &lldvs, &lwork0,
         &lldwork, rwork, &dummyi, &linfo);
   lldwork = REAL_PART(lwork0);

   if (linfo == 0) {
      SCALAR *work = NULL;
      CHKERR(Num_malloc_Sprimme(lldwork, &work, ctx));
      XGEES(jobvs, "N", NULL, &ln, a, &llda, &lsdim, w, vs, &lldvs, work,
            &lldwork, rwork, &dummyi, &linfo);
      CHKERR(Num_free_Sprimme(work, ctx));
   }

   CHKERR(Num_free_Rprimme(rwork, ctx));

   CHKERRM(linfo != 0, PRIMME_LAPACK_FAILURE, "Error in xgees with info %d",
          (int)linfo);
   return 0;
}
#endif /* defined(USE_COMPLEX) && ((!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)) */

/*******************************************************************************
 * Subroutines for dense singular value decomposition
 ******************************************************************************/
 
#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
TEMPLATE_PLEASE
int Num_gesvd_Sprimme(const char *jobu, const char *jobvt, PRIMME_INT m,
      PRIMME_INT n, SCALAR *a, PRIMME_INT lda, REAL *s, SCALAR *u,
      PRIMME_INT ldu, SCALAR *vt, PRIMME_INT ldvt, primme_context ctx) {

   PRIMME_BLASINT lm, ln, llda, lldu, lldvt;
   CHKERR(to_blas_int(m, &lm));
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(lda, &llda));
   CHKERR(to_blas_int(ldu, &lldu));
   CHKERR(to_blas_int(ldvt, &lldvt));
   PRIMME_BLASINT lldwork = 0;
   PRIMME_BLASINT linfo = 0;
#ifdef USE_COMPLEX
   REAL dummyr = 0.0;
#endif

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return 0;

   /* Call to know the optimal workspace */

   lldwork = -1;
   SCALAR lwork0 = 0;
   XGESVD(jobu, jobvt, &lm, &ln, a, &llda, s, u, &lldu, vt, &lldvt, &lwork0,
          &lldwork,
#ifdef USE_COMPLEX
          &dummyr,
#endif
         &linfo);
   lldwork = REAL_PART(lwork0);

   if (linfo == 0) {
      SCALAR *work = NULL;
      CHKERR(Num_malloc_Sprimme(lldwork, &work, ctx));
#  ifdef USE_COMPLEX
      REAL *rwork;
      CHKERR(Num_malloc_Rprimme(5*n, &rwork, ctx));
#  endif
      XGESVD(jobu, jobvt, &lm, &ln, a, &llda, s, u, &lldu, vt, &lldvt, work,
            &lldwork,
#ifdef USE_COMPLEX
            rwork,
#endif
            &linfo);
      CHKERR(Num_free_Sprimme(work, ctx));
#  ifdef USE_COMPLEX
      CHKERR(Num_free_Rprimme(rwork, ctx));
#  endif
   }

   CHKERRM(linfo != 0, PRIMME_LAPACK_FAILURE, "Error in xgesvd with info %d",
          (int)linfo);
   return 0;
}
#endif /* (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF) */

/*******************************************************************************
 * Subroutine Num_hetrf_Sprimme - LL^H factorization with pivoting
 ******************************************************************************/

#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
TEMPLATE_PLEASE
int Num_hetrf_Sprimme(const char *uplo, PRIMME_INT n, SCALAR *a, PRIMME_INT lda,
      int *ipivot, primme_context ctx) {
#ifndef USE_ZGESV

   PRIMME_BLASINT ln, llda;
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(lda, &llda));
   PRIMME_BLASINT *lipivot;
   PRIMME_BLASINT lldwork = 0;
   PRIMME_BLASINT linfo = 0; 

   /* Zero dimension matrix may cause problems */

   if (n == 0) return 0;

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      CHKERR(Num_malloc_iblasprimme(n, &lipivot, ctx));
   } else {
      lipivot = (PRIMME_BLASINT *)ipivot; /* cast avoid compiler warning */
   }

   /* Call to know the optimal workspace */
   
   lldwork = -1;
   SCALAR lwork0 = 0;
   XHETRF(uplo, &ln, a, &llda, lipivot, &lwork0, &lldwork, &linfo);
   lldwork = REAL_PART(lwork0);

   if (linfo == 0) {
      SCALAR *work;
      CHKERR(Num_malloc_Sprimme(lldwork, &work, ctx));
      XHETRF(uplo, &ln, a, &llda, lipivot, work, &lldwork, &linfo);
      CHKERR(Num_free_Sprimme(work, ctx));
   }

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      int i;
      for(i=0; i<n; i++)
         ipivot[i] = (int)lipivot[i];
      CHKERR(Num_free_iblasprimme(lipivot, ctx));
   }

   CHKERRM(linfo != 0, PRIMME_LAPACK_FAILURE, "Error in xhetrf with info %d",
          (int)linfo);
   return 0;

#else /* USE_ZGESV */

   (void)ipivot;
   (void)ctx;
   /* Lapack's R core library doesn't have zhetrf. The functionality is       */
   /* implemented by replacing the input matrix with a full general matrix.   */
   /* And Num_zhetrs_Sprimme will solve the general linear system.            */

   /* Copy the given upper/lower triangular part into the lower/upper part    */

   PRIMME_INT i, j;
   if (*uplo == 'L' || *uplo == 'l') {
      for (i=0; i<n; i++) {
         for (j=i+1; j<n; j++) {
            a[lda*j+i] = CONJ(a[lda*i+j]);
         }
      }
   }
   else {
      for (i=1; i<n; i++) {
         for (j=0; j<i; j++) {
            a[lda*j+i] = CONJ(a[lda*i+j]);
         }
      }
   }

   return 0;

#endif
}
#endif /* (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF) */

/*******************************************************************************
 * Subroutine Num_hetrs_Sprimme - b = A\b where A may store a LL^H factorization
 ******************************************************************************/
 
#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
TEMPLATE_PLEASE
int Num_hetrs_Sprimme(const char *uplo, PRIMME_INT n, PRIMME_INT nrhs,
      SCALAR *a, PRIMME_INT lda, int *ipivot, SCALAR *b, PRIMME_INT ldb,
      primme_context ctx) {

   PRIMME_BLASINT ln, llda, lldb, lnrhs;
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(nrhs, &lnrhs));
   CHKERR(to_blas_int(lda, &llda));
   CHKERR(to_blas_int(ldb, &lldb));
   PRIMME_BLASINT *lipivot;
   PRIMME_BLASINT linfo = 0; 

   /* Zero dimension matrix may cause problems */

   if (n == 0 || nrhs == 0) return 0;

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      CHKERR(Num_malloc_iblasprimme(n, &lipivot, ctx));
      PRIMME_INT i;
      for(i=0; i<n; i++) {
         lipivot[i] = (PRIMME_BLASINT)ipivot[i];
      }
   } else {
      lipivot = (PRIMME_BLASINT *)ipivot; /* cast avoid compiler warning */
   }

#ifndef USE_ZGESV
   XHETRS(uplo, &ln, &lnrhs, a, &llda, lipivot, b, &lldb, &linfo);
   CHKERRM(linfo != 0, PRIMME_LAPACK_FAILURE, "Error in xhetrs with info %d",
          (int)linfo);
#else /* USE_ZGESV */
   (void)uplo;
   XGESV(&ln, &lnrhs, a, &llda, lipivot, b, &lldb, &linfo);
   CHKERRM(linfo != 0, PRIMME_LAPACK_FAILURE, "Error in xgesv with info %d",
          (int)linfo);
#endif

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      CHKERR(Num_free_iblasprimme(lipivot, ctx));
   }

   return 0;
}
#endif /* (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF) */

/*******************************************************************************
 * Subroutine Num_potrf_Sprimme - Cholesky factorization
 ******************************************************************************/
 
#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
TEMPLATE_PLEASE
int Num_potrf_Sprimme(const char *uplo, PRIMME_INT n, SCALAR *a, PRIMME_INT lda,
      int *info, primme_context ctx) {

   (void)ctx;

   PRIMME_BLASINT ln, llda;
   CHKERR(to_blas_int(n, &ln));
   CHKERR(to_blas_int(lda, &llda));
   PRIMME_BLASINT linfo = 0; 

   /* Zero dimension matrix may cause problems */
   if (n == 0) {
      if (info) *info = 0;
      return 0;
   }

   XPOTRF(uplo, &ln, a, &llda, &linfo);
   CHKERRM(info == NULL && linfo != 0, PRIMME_LAPACK_FAILURE,
         "Error in xpotrf with info %d\n", (int)linfo);
   if (info) *info = (int)linfo;

   return 0;
}
#endif /* (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF) */


/*******************************************************************************
 * Subroutine Num_trsm_Sprimme - b = op(A)\b, where A is triangular
 ******************************************************************************/
 
#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
TEMPLATE_PLEASE
int Num_trsm_Sprimme(const char *side, const char *uplo, const char *transa,
      const char *diag, PRIMME_INT m, PRIMME_INT n, HSCALAR alpha, SCALAR *a,
      PRIMME_INT lda, SCALAR *b, PRIMME_INT ldb, primme_context ctx) {

   (void)ctx;
   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldb = ldb;

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return 0;

   XTRSM(side, uplo, transa, diag, &lm, &ln, (BLASSCALAR *)&alpha,
         (BLASSCALAR *)a, &llda, (BLASSCALAR *)b,
         &lldb STRING_LEN1 STRING_LEN1 STRING_LEN1 STRING_LEN1);

   return 0;
}
#endif /* (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF) */

TEMPLATE_PLEASE
int Num_trsm_hd_Sprimme(const char *side, const char *uplo, const char *transa,
      const char *diag, PRIMME_INT m, PRIMME_INT n, HSCALAR alpha, HSCALAR *a,
      PRIMME_INT lda, SCALAR *b, PRIMME_INT ldb, primme_context ctx) {

   /* Cast the matrix b to HSCALAR */

   HSCALAR *bf = NULL;
   PRIMME_INT ldbf;
   CHKERR(
         Num_matrix_astype_Sprimme(b, m, n, ldb, PRIMME_OP_SCALAR, (void **)&bf,
               &ldbf, PRIMME_OP_HSCALAR, 1 /* alloc */, 1 /* copy */, ctx));

   /* Call the kernel */

   CHKERR(Num_trsm_SHprimme(
         side, uplo, transa, diag, m, n, alpha, a, lda, bf, ldbf, ctx));

   /* Copy back and destroy the cast matrix */

   PRIMME_INT ldb_ = ldb;
   CHKERR(Num_matrix_astype_Sprimme(bf, m, n, ldbf, PRIMME_OP_HSCALAR,
         (void **)&b, &ldb_, PRIMME_OP_SCALAR, -1 /* destroy */, 1 /* copy */,
         ctx));

   return 0;
}

/*******************************************************************************
 * Subroutine Num_getrf_Sprimme - Factorize A=LU
 ******************************************************************************/
 
#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
TEMPLATE_PLEASE
int Num_getrf_Sprimme(PRIMME_INT m, PRIMME_INT n, SCALAR *a, PRIMME_INT lda,
      int *ipivot, primme_context ctx) {
   (void)ctx;

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT *lipivot = NULL;
   PRIMME_BLASINT linfo = 0;

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return 0;

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      CHKERR(MALLOC_PRIMME(n, &lipivot));
   } else {
      lipivot = (PRIMME_BLASINT *)ipivot; /* cast avoid compiler warning */
   }

   XGETRF(&lm, &ln, a, &llda, lipivot, &linfo);

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      PRIMME_INT i;
      if (ipivot) for(i=0; i<n; i++)
         ipivot[i] = (int)lipivot[i];
      free(lipivot);
   }

   CHKERRM(linfo != 0, PRIMME_LAPACK_FAILURE, "Error in xgesv with info %d",
         (int)linfo);

   return 0;
}
#endif /* (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF) */

/*******************************************************************************
 * Subroutine Num_getrs_Sprimme - Computes A\X where A=LU computed with getrf
 ******************************************************************************/
 
#if (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)
TEMPLATE_PLEASE
int Num_getrs_Sprimme(const char *trans, PRIMME_INT n, PRIMME_INT nrhs,
      SCALAR *a, PRIMME_INT lda, int *ipivot, SCALAR *b, PRIMME_INT ldb,
      primme_context ctx) {
   (void)ctx;

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lnrhs = nrhs;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldb = ldb;
   PRIMME_BLASINT *lipivot = NULL;
   PRIMME_BLASINT linfo = 0;

   /* Zero dimension matrix may cause problems */
   if (n == 0 || nrhs == 0) return 0;

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      CHKERR(MALLOC_PRIMME(n, &lipivot));
      PRIMME_INT i;
      for(i=0; i<n; i++) {
         lipivot[i] = (PRIMME_BLASINT)ipivot[i];
      }
   } else {
      lipivot = (PRIMME_BLASINT *)ipivot; /* cast avoid compiler warning */
   }

   XGETRS(trans, &ln, &lnrhs, a, &llda, lipivot, b, &lldb, &linfo);

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      free(lipivot);
   }

   CHKERRM(linfo != 0, PRIMME_LAPACK_FAILURE, "Error in xgetrs with info %d",
         (int)linfo);

   return 0;
}
#endif /* (!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF) */

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
int Num_compute_gramm_ddh_Sprimme(SCALAR *X, PRIMME_INT m, PRIMME_INT n,
      PRIMME_INT ldX, SCALAR *Y, PRIMME_INT ldY, HSCALAR alpha, HSCALAR *H,
      PRIMME_INT ldH, int isherm, primme_context ctx) {

   (void)isherm;

   CHKERR(Num_gemm_ddh_Sprimme(
         "C", "N", n, n, m, 1.0, X, ldX, Y, ldY, alpha, H, ldH, ctx));

   return 0;
}

#endif /* USE_HOST */

#endif /* SUPPORTED_TYPE */
