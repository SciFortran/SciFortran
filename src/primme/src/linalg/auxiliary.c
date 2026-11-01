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
 * File: auxiliary.c
 *
 * Purpose - Miscellanea functions to copy and permuting matrices.
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../linalg/auxiliary.c"
#endif

#include <string.h>   /* memmove */
#include <assert.h>
#include <math.h>
#include "numerical.h"

#ifdef SUPPORTED_TYPE


/******************************************************************************
 * Function Num_matrix_astype - Create a matrix, y (if asked) and copy the
 *    matrix x into y (if asked).
 *
 *    do_alloc do_copy xt==yt action
 *    --------------------------------------------
 *        0       0      *    do nothing
 *        0       1      *    copy x into y
 *        1       0      0    create y of type yt
 *        1       1      0    create y of type yt and copy x into y
 *        1       *      1    assign *y = x
 *       -1       *      1    do nothing
 *       -1       0      0    destroy x
 *       -1       1      0    copy x into y and destroy x
 *
 * PARAMETERS
 * ---------------------------
 * x           The source matrix
 * m           The number of rows of x
 * n           The number of columns of x
 * ldx         The leading dimension of x
 * xt          The datatype of x
 * y           On output y = x
 * ldy         The leading dimension of y
 * yt          The datatype of y
 *
 * NOTE: x and y *cannot* partially overlap
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_matrix_astype_Sprimme(void *x, PRIMME_INT m, PRIMME_INT n,
      PRIMME_INT ldx, primme_op_datatype xt, void **y, PRIMME_INT *ldy,
      primme_op_datatype yt, int do_alloc, int do_copy, primme_context ctx) {

   /* Replace primme_op_default */

   if (xt == primme_op_default) xt = PRIMME_OP_SCALAR;
   if (yt == primme_op_default) yt = PRIMME_OP_SCALAR;
   
   /* Call the function that y has the type of the SCALAR */

   if (yt != PRIMME_OP_SCALAR) {
      switch(yt) {
#ifdef SUPPORTED_HALF_TYPE
      case primme_op_half:   return Num_matrix_astype_Shprimme(x, m, n, ldx, xt, y, ldy, yt, do_alloc, do_copy, ctx);
#endif
#ifndef PRIMME_WITHOUT_FLOAT
      case primme_op_float:  return Num_matrix_astype_Ssprimme(x, m, n, ldx, xt, y, ldy, yt, do_alloc, do_copy, ctx);
#endif
      case primme_op_double: return Num_matrix_astype_Sdprimme(x, m, n, ldx, xt, y, ldy, yt, do_alloc, do_copy, ctx);
#ifdef PRIMME_WITH_NATIVE_COMPLEX_QUAD
      case primme_op_quad:   return Num_matrix_astype_Sqprimme(x, m, n, ldx, xt, y, ldy, yt, do_alloc, do_copy, ctx);
#endif
#ifdef USE_HOST
      case primme_op_int:    return Num_matrix_astype_iprimme(x, m, n, ldx, xt, y, ldy, yt, do_alloc, do_copy, ctx);
#endif
      default: CHKERR(PRIMME_FUNCTION_UNAVAILABLE);
      }
   }

   /* Quick exit */

   if (xt == PRIMME_OP_SCALAR && do_alloc) {
      *y = x;
      if (ldy) *ldy = ldx;
      return 0;
   }

   /* Create workspace for y and copy x on y */

   SCALAR *y0 = NULL;
   PRIMME_INT ldy0 = 0;
   if (do_alloc > 0) {
      Mem_keep_frame(
            ctx); /* The next allocation will not be freed in this function */
      CHKERR(Num_malloc_Sprimme(m * n, &y0, ctx));
      *y = (void*)y0;
      ldy0 = m;
      if (ldy) *ldy = m;
   } else {
      y0 = (SCALAR*)*y;
      ldy0 = (ldy ? *ldy : 1);
   }

   if (do_copy && x != NULL) {
      CHKERR(Num_copy_Tmatrix_Sprimme(x, xt, m, n, ldx, y0, ldy0, ctx));
   }

   /* Destroy x if asked */

   if (do_alloc < 0 && x != y0) CHKERR(Num_free_Sprimme((SCALAR *)x, ctx));

   return 0;
}

#ifdef USE_HOST
#ifdef USE_DOUBLE

TEMPLATE_PLEASE
int Num_matrix_astype_iprimme(void *x, PRIMME_INT m, PRIMME_INT n,
      PRIMME_INT ldx, primme_op_datatype xt, void **y, PRIMME_INT *ldy,
      primme_op_datatype yt, int do_alloc, int do_copy, primme_context ctx) {

   /* Replace primme_op_default */

   if (xt == primme_op_default) xt = primme_op_int;
   if (yt == primme_op_default) yt = primme_op_int;
   
   /* Call the function that y has the type of the SCALAR */

   if (yt != primme_op_int) {
      switch(yt) {
#ifdef SUPPORTED_HALF_TYPE
      case primme_op_half:   return Num_matrix_astype_Shprimme(x, m, n, ldx, xt, y, ldy, yt, do_alloc, do_copy, ctx);
#endif
#ifndef PRIMME_WITHOUT_FLOAT
      case primme_op_float:  return Num_matrix_astype_Ssprimme(x, m, n, ldx, xt, y, ldy, yt, do_alloc, do_copy, ctx);
#endif
      case primme_op_double: return Num_matrix_astype_Sdprimme(x, m, n, ldx, xt, y, ldy, yt, do_alloc, do_copy, ctx);
#ifdef PRIMME_WITH_NATIVE_COMPLEX_QUAD
      case primme_op_quad:   return Num_matrix_astype_Sqprimme(x, m, n, ldx, xt, y, ldy, yt, do_alloc, do_copy, ctx);
#endif
      default: CHKERR(PRIMME_FUNCTION_UNAVAILABLE);
      }
   }

   /* Quick exit */

   if (xt == primme_op_int && do_alloc) {
      *y = x;
      if (ldy) *ldy = ldx;
      return 0;
   }

   /* Create workspace for y and copy x on y */

   int *y0 = NULL;
   PRIMME_INT ldy0 = 0;
   if (do_alloc > 0) {
      Mem_keep_frame(
            ctx); /* The next allocation will not be freed in this function */
      CHKERR(Num_malloc_iprimme(m * n, &y0, ctx));
      *y = (void*)y0;
      ldy0 = m;
      if (ldy) *ldy = m;
   } else {
      y0 = (int*)*y;
      ldy0 = (ldy ? *ldy : 1);
   }

   if (do_copy && x != NULL) {
      CHKERR(Num_copy_Tmatrix_iprimme(x, xt, m, n, ldx, y0, ldy0, ctx));
   }

   /* Destroy x if asked */

   if (do_alloc < 0 && x != y0) CHKERR(Num_free_iprimme((int *)x, ctx));

   return 0;
}

#endif /* USE_DOUBLE */
#endif /* USE_HOST */

/******************************************************************************
 * Function Num_copy_matrix_astype - copy the matrix x into y.
 *
 * PARAMETERS
 * ---------------------------
 * x           The source matrix
 * xm0         The starting row
 * xn0         The starting column
 * m           The number of rows of x
 * n           The number of columns of x
 * ldx         The leading dimension of x
 * xt          The datatype of x
 * y           On output y = x
 * xm0         The starting row
 * xn0         The starting column
 * ldy         The leading dimension of y
 * yt          The datatype of y
 *
 * NOTE: x and y *cannot* partially overlap
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_copy_matrix_astype_Sprimme(void *x, PRIMME_INT xm0, PRIMME_INT xn0,
      PRIMME_INT m, PRIMME_INT n, PRIMME_INT ldx, primme_op_datatype xt,
      void *y, PRIMME_INT ym0, PRIMME_INT yn0, PRIMME_INT ldy,
      primme_op_datatype yt, primme_context ctx) {

   /* Replace primme_op_default */

   if (xt == primme_op_default) xt = PRIMME_OP_SCALAR;
   if (yt == primme_op_default) yt = PRIMME_OP_SCALAR;
   
   /* Call the function that y has the type of the SCALAR */

   if (yt != PRIMME_OP_SCALAR) {
      switch(yt) {
#ifdef SUPPORTED_HALF_TYPE
      case primme_op_half:   return Num_copy_matrix_astype_Shprimme(x, xm0, xn0, m, n, ldx, xt, y, ym0, yn0, ldy, yt, ctx);
#endif
#ifndef PRIMME_WITHOUT_FLOAT
      case primme_op_float:  return Num_copy_matrix_astype_Ssprimme(x, xm0, xn0, m, n, ldx, xt, y, ym0, yn0, ldy, yt, ctx);
#endif
      case primme_op_double: return Num_copy_matrix_astype_Sdprimme(x, xm0, xn0, m, n, ldx, xt, y, ym0, yn0, ldy, yt, ctx);
#ifdef PRIMME_WITH_NATIVE_COMPLEX_QUAD
      case primme_op_quad:   return Num_copy_matrix_astype_Sqprimme(x, xm0, xn0, m, n, ldx, xt, y, ym0, yn0, ldy, yt, ctx);
#endif
      default: CHKERR(PRIMME_FUNCTION_UNAVAILABLE);
      }
   }

   size_t xt_size;
   CHKERR(Num_sizeof_Sprimme(xt, &xt_size));
   return Num_copy_Tmatrix_Sprimme(&((char *)x)[xt_size * (xm0 + ldx * xn0)],
         xt, m, n, ldx, &((SCALAR *)y)[ym0 + ldy * yn0], ldy, ctx);
}


/******************************************************************************
 * Function Num_sizeof_Sprimme - Return the size of an element with the given
 *    type.
 *
 * INPUT/OUTPUT PARAMETERS
 * ---------------------------
 * t    Type
 * s    returned size of the type in bytes
 *
 * RETURN
 * ------
 * int  error code
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_sizeof_Sprimme(primme_op_datatype t, size_t *s) {

   if (t == primme_op_default) t = PRIMME_OP_SCALAR;

   *s = 0;

   switch(t) {
#  ifdef SUPPORTED_HALF_TYPE
   case primme_op_half:    *s = sizeof(PRIMME_HALF); break;
#  endif
   case primme_op_float:   *s = sizeof(float); break;
   case primme_op_double:  *s = sizeof(double); break;
#  ifdef PRIMME_WITH_NATIVE_QUAD
   case primme_op_quad:    *s = sizeof(PRIMME_QUAD); break;
#  endif
   case primme_op_int:     *s = sizeof(int); break;
   default:                return PRIMME_FUNCTION_UNAVAILABLE;
   }

#ifdef USE_COMPLEX
   *s *= 2;
#endif

   return 0;
}

/******************************************************************************
 * Function Num_machine_epsilon_Sprimme - Return the machine epsilon of the
 *    type.
 *
 * INPUT/OUTPUT PARAMETERS
 * ---------------------------
 * t    Type
 * eps  The returned machine epsilon
 *
 * RETURN
 * ------
 * int  error code
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_machine_epsilon_Sprimme(primme_op_datatype t, double *eps) {

   /* Replace primme_op_default */

   if (t == primme_op_default) t = PRIMME_OP_SCALAR;
   
   /* Call the function that t has the type of the SCALAR */

   if (t != PRIMME_OP_SCALAR) {
      switch(t) {
#ifdef SUPPORTED_HALF_TYPE
      case primme_op_half:   return Num_machine_epsilon_Shprimme(t, eps);
#endif
#ifndef PRIMME_WITHOUT_FLOAT
      case primme_op_float:  return Num_machine_epsilon_Ssprimme(t, eps);
#endif
      case primme_op_double: return Num_machine_epsilon_Sdprimme(t, eps);
#ifdef PRIMME_WITH_NATIVE_COMPLEX_QUAD
      case primme_op_quad:   return Num_machine_epsilon_Sqprimme(t, eps);
#endif
      default: return PRIMME_FUNCTION_UNAVAILABLE;
      }
   }

   if (eps) *eps = MACHINE_EPSILON;

   return 0;
}

#ifdef USE_HOST

/******************************************************************************
 * Function Num_copy_matrix_conj - Copy the matrix x' into y
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
 * NOTE: x and y *can't* overlap
 *
 ******************************************************************************/

#if !defined(USE_HALFCOMPLEX) || defined(PRIMME_WITH_NATIVE_COMPLEX_HALF)
TEMPLATE_PLEASE
int Num_copy_matrix_conj_Sprimme(SCALAR *x, PRIMME_INT m, PRIMME_INT n,
      PRIMME_INT ldx, SCALAR *y, PRIMME_INT ldy, primme_context ctx) {
   (void)ctx;

   PRIMME_INT i, j;

   assert(m == 0 || n == 0 || (ldx >= m && ldy >= n));

   /* TODO: assert x and y don't overlap */
   for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
         y[j * ldy + i] = CONJ(x[i * ldx + j]);

   return 0;
}
#endif

#ifdef USE_DOUBLE
TEMPLATE_PLEASE
int Num_zero_matrix_Tprimme(void *x, primme_op_datatype xt, PRIMME_INT m,
      PRIMME_INT n, PRIMME_INT ldx, primme_context ctx) {

   switch (xt) {
#  ifdef SUPPORTED_HALF_TYPE
      case primme_op_half:   return Num_zero_matrix_hprimme((PRIMME_HALF*)x, m, n, ldx, ctx);
#  endif
#  ifndef PRIMME_WITHOUT_FLOAT
      case primme_op_float:  return Num_zero_matrix_sprimme((float*)      x, m, n, ldx, ctx);
#  endif
      case primme_op_double: return Num_zero_matrix_dprimme((double*)     x, m, n, ldx, ctx);
#  ifdef PRIMME_WITH_NATIVE_COMPLEX_QUAD
      case primme_op_quad:   return Num_zero_matrix_qprimme((PRIME_QUAD*) x, m, n, ldx, ctx);
#  endif
      case primme_op_int: {
         PRIMME_INT i, j;
         int *xi = (int *)x;
         for (i = 0; i < n; i++)
            for (j = 0; j < m; j++) xi[i * ldx + j] = 0;
         break;
      }
      default: CHKERR(PRIMME_FUNCTION_UNAVAILABLE);
   }

   return 0;
}
#endif /* USE_DOUBLE */

/******************************************************************************
 * Function Num_copy_trimatrix - Copy the upper/lower triangular part of the
 *    matrix x into y
 *
 * PARAMETERS
 * ---------------------------
 * x           The source matrix
 * m           The number of rows of x
 * n           The number of columns of x
 * ldx         The leading dimension of x
 * ul          if 0, copy the upper part; otherwise copy the lower part 
 * i0          The row index that diagonal starts from
 * y           On output y = x
 * ldy         The leading dimension of y
 * zero        If nonzero, zero the triangular part not copied
 *
 * NOTE: the distance between x and y can be less than ldx, or
 *       x and y *cannot* overlap at all
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_copy_trimatrix_Sprimme(SCALAR *x, int m, int n, int ldx, int ul,
      int i0, SCALAR *y, int ldy, int zero) {

   int i, j, jm;

   assert(m == 0 || n == 0 || (ldx >= m && ldy >= m));
   if (x == y) return 0;
   if (ul == 0) {
      /* Copy upper part */

      if (ldx == ldy && (y > x ? y-x : x-y) < ldx) {
         /* x and y overlap */
         for (i=0; i<n; i++) {
            memmove(&y[i*ldy], &x[i*ldx], sizeof(SCALAR)*min(i0+i+1, m));
            /* zero lower part*/
            if (zero) for (j=min(i0+i+1, m); j<m; j++) SET_ZERO(y[i*ldy+j]);
         }
      }
      else {
         /* x and y don't overlap */
         for (i=0; i<n; i++) {
            for (j=0, jm=min(i0+i+1, m); j<jm; j++)
               y[i*ldy+j] = x[i*ldx+j];
            /* zero lower part*/
            if (zero) for (j=min(i0+i+1, m); j<m; j++) SET_ZERO(y[i*ldy+j]);
         }
      }
   }
   else {
      /* Copy lower part */

      if (ldx == ldy && (y > x ? y-x : x-y) < ldx) {
         /* x and y overlap */
         for (i=0; i<n; i++) {
            memmove(&y[i*ldy+i0+i], &x[i*ldx+i0+i], sizeof(SCALAR)*(m-min(i0+i, m)));
            /* zero upper part*/
            if (zero) for (j=0, jm=min(i0+i, m); j<jm; j++) SET_ZERO(y[i*ldy+j]);
         }
      }
      else {
         /* x and y don't overlap */
         for (i=0; i<n; i++) {
            for (j=i+i0; j<m; j++)
               y[i*ldy+j] = x[i*ldx+j];
            /* zero upper part*/
            if (zero) for (j=0, jm=min(i0+i, m); j<jm; j++) SET_ZERO(y[i*ldy+j]);
         }
      }
   }

   return 0;
}


/******************************************************************************
 * Function Num_copy_trimatrix_compact - Copy the upper triangular part of the matrix x
 *    into y contiguously, i.e., y has all columns of x row-stacked
 *
 * PARAMETERS
 * ---------------------------
 * x           The source upper triangular matrix
 * m           The number of rows of x
 * n           The number of columns of x
 * ldx         The leading dimension of x
 * i0          The row index that diagonal starts from
 * y           On output y = x and nonzero elements of y are contiguous
 * ly          Output the final length of y
 *
 * NOTE: x and y *cannot* overlap
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_copy_trimatrix_compact_Sprimme(SCALAR *x, PRIMME_INT m, int n,
      PRIMME_INT ldx, int i0, SCALAR *y, int *ly) {

   int i, j, k;

   assert(m == 0 || n == 0 || ldx >= m);
   if (m < n) return -1;

   for (i=0, k=0; i<n; i++)
      for (j=0; j<=i+i0; j++)
         y[k++] = x[i*ldx+j];
   if (ly) *ly = k;

   return 0;
}

/******************************************************************************
 * Function Num_copy_compact_trimatrix - Copy y into the upper triangular part of the
 *    matrix x
 *
 * PARAMETERS
 * ---------------------------
 * x           The source vector
 * m           The number of rows of y
 * n           The number of columns of y
 * i0          The row index that diagonal starts from
 * y           On output the upper triangular part of y has x
 * ldy         The leading dimension of y
 *
 * NOTE: x and y *cannot* overlap
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_copy_compact_trimatrix_Sprimme(SCALAR *x, PRIMME_INT m, int n, int i0,
      SCALAR *y, int ldy) {

   int i, j, k;

   assert(m == 0 || n == 0 || (ldy >= m && m >= n));
   if (m < n) return -1;

   for (i=n-1, k=(n+1)*n/2+i0*n-1; i>=0; i--)
      for (j=i+i0; j>=0; j--)
         y[i*ldy+j] = x[k--];

   return 0;
}

/*******************************************************************************
 * Subroutine compute_submatrix - This subroutine computes the nX x nX submatrix
 *    R = X'*H*X, where H stores the upper triangular part of a symmetric matrix,
 *    or H is a full non-Hermitian matrix.
 *    
 * Input parameters
 * ----------------
 * X        The coefficient vectors retained from the previous iteration
 *
 * nX       Number of columns of X
 *
 * H        Matrix
 *
 * nH       Dimension of H
 *
 * ldH      Leading dimension of H
 *
 * isherm   whether H is Hermitian
 *
 * ldR      Leading dimension of R
 *
 * Output parameters
 * -----------------
 * R - nX x nX matrix computed 
 *
 ******************************************************************************/

#if !defined(USE_HALF) && !defined(USE_HALFCOMPLEX)
TEMPLATE_PLEASE
int compute_submatrix_Sprimme(SCALAR *X, int nX, int ldX, SCALAR *H, int nH,
      int ldH, int isherm, SCALAR *R, int ldR, primme_context ctx) {

   if (nH == 0 || nX == 0) return 0;

   SCALAR *rwork;
   CHKERR(Num_malloc_Sprimme((size_t)nH * (size_t)nX, &rwork, ctx));

   /* rwork = H * X */

   Num_zero_matrix_Sprimme(rwork, nH, nX, nH, ctx);
   if (isherm) {
      CHKERR(Num_hemm_Sprimme(
            "L", "U", nH, nX, 1.0, H, ldH, X, ldX, 0.0, rwork, nH, ctx));
   } else {
      CHKERR(Num_gemm_Sprimme(
            "N", "N", nH, nX, nH, 1.0, H, ldH, X, ldX, 0.0, rwork, nH, ctx));
   }

   /* R = X' * rwork */

   Num_zero_matrix_Sprimme(R, nX, nX, ldR, ctx);
   CHKERR(Num_gemm_Sprimme(
         "C", "N", nX, nX, nH, 1.0, X, ldX, rwork, nH, 0.0, R, ldR, ctx));
   CHKERR(Num_free_Sprimme(rwork, ctx));

  return 0;
}
#endif /* !defined(USE_HALF) && !defined(USE_HALFCOMPLEX) */

#endif /* USE_HOST */

/******************************************************************************
 * Function Num_copy_matrix_columns - Copy the matrix x(:,xin) into y(:,yin)
 *
 * PARAMETERS
 * ---------------------------
 * x           The source matrix
 * m           The number of rows of x
 * xin         The column indices to copy
 * n           The number of columns of x
 * ldx         The leading dimension of x
 * y           On output y(:,yin) = x(:,xin)
 * yin         The column indices of y to be modified
 * ldy         The leading dimension of y
 *
 * NOTE: x(:,xin) and y(:,yin) *cannot* overlap
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_copy_matrix_columns_Sprimme(SCALAR *x, PRIMME_INT m, int *xin, int n,
                                     PRIMME_INT ldx, SCALAR *y, int *yin,
                                     PRIMME_INT ldy, primme_context ctx) {

   int i;

   /* TODO: assert x and y don't overlap */
   for (i = 0; i < n; i++) {
      Num_copy_Sprimme(m, &x[(xin ? xin[i] : i) * ldx], 1,
            &y[(yin ? yin[i] : i) * ldy], 1, ctx);
   }

   return 0;
}

/******************************************************************************
 * Function Num_copy_matrix_rows - Copy the matrix x(xin,:) into y(yin,:)
 *
 * PARAMETERS
 * ---------------------------
 * x           The source matrix
 * xim         The row indices to copy
 * m           The number of rows of x
 * n           The number of columns of x
 * ldx         The leading dimension of x
 * y           On output y(yin,:) = x(xin,:)
 * yim         The row indices of y to be modified
 * ldy         The leading dimension of y
 *
 * NOTE: x(xin,:) and y(yin,:) *cannot* overlap
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_copy_matrix_rows_Sprimme(SCALAR *x, int *xim, int m, int n,
                                     PRIMME_INT ldx, SCALAR *y, int *yim,
                                     PRIMME_INT ldy, primme_context ctx) {

   int i;

   /* TODO: assert x and y don't overlap */
   for (i = 0; i < m; i++) {
      Num_copy_Sprimme(
            n, &x[xim ? xim[i] : i], ldx, &y[yim ? yim[i] : i], ldy, ctx);
   }

   return 0;
}

/******************************************************************************
 * Subroutine permute_vecs - This routine permutes a set of vectors according
 *            to a permutation array perm.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * m, n, ld    The number of rows and columns and the leading dimension of
 *vecs perm        The permutation of the columns rwork       Temporary space
 *of size the number of rows iwork       Temporary space of size the number
 *of columns
 *
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * vecs        The matrix whose columns will be reordered
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int permute_vecs_Sprimme(SCALAR *vecs, PRIMME_INT m, int n, PRIMME_INT ld,
                         int *perm_, primme_context ctx) {

   int currentIndex;     /* Index of vector in sorted order                   */
   int sourceIndex;      /* Position of out-of-order vector in original order */
   int destinationIndex; /* Position of out-of-order vector in sorted order */
   int tempIndex;        /* Used to swap                                      */
   int *perm;            /* A copy of perm_                                   */
   SCALAR *rwork;        /* vector buffer */

   CHKERR(Num_malloc_iprimme(n, &perm, ctx));
   CHKERR(Num_malloc_Sprimme(m, &rwork, ctx));

   /* Check perm_ is a permutation */

#ifndef NDEBUG
   for (tempIndex = 0; tempIndex < n; tempIndex++)
      perm[tempIndex] = 0;
   for (tempIndex = 0; tempIndex < n; tempIndex++) {
      assert(0 <= perm_[tempIndex] && perm_[tempIndex] < n);
      perm[perm_[tempIndex]] = 1;
   }
   for (tempIndex = 0; tempIndex < n; tempIndex++)
      assert(perm[tempIndex] == 1);
#endif

   /* Copy of perm_ into perm, to avoid to modify the input permutation */

   for (tempIndex = 0; tempIndex < n; tempIndex++)
      perm[tempIndex] = perm_[tempIndex];

   /* Continue until all vectors are in the sorted order */

   currentIndex = 0;
   while (1) {

      /* Find a vector that does not belong in its original position */
      while ((currentIndex < n) && (perm[currentIndex] == currentIndex)) {
         currentIndex++;
      }

      /* Return if they are in the sorted order */
      if (currentIndex >= n) {
         break;
      }

      /* Copy the vector to a buffer for swapping */
      Num_copy_Sprimme(m, &vecs[currentIndex * ld], 1, rwork, 1, ctx);

      destinationIndex = currentIndex;
      /* Copy vector perm[destinationIndex] into position destinationIndex */

      while (perm[destinationIndex] != currentIndex) {

         sourceIndex = perm[destinationIndex];
         Num_copy_Sprimme(m, &vecs[sourceIndex * ld], 1,
                          &vecs[destinationIndex * ld], 1, ctx);
         tempIndex = perm[destinationIndex];
         perm[destinationIndex] = destinationIndex;
         destinationIndex = tempIndex;
      }

      /* Copy the vector from the buffer to where it belongs */
      Num_copy_Sprimme(m, rwork, 1, &vecs[destinationIndex * ld], 1, ctx);
      perm[destinationIndex] = destinationIndex;

      currentIndex++;
   }

   /* Check permutation */
   for (currentIndex = 0; currentIndex < n; currentIndex++)
      assert(perm[currentIndex] == currentIndex);

   CHKERR(Num_free_iprimme(perm, ctx));
   CHKERR(Num_free_Sprimme(rwork, ctx));

   return 0;
}

#ifdef USE_DOUBLE
TEMPLATE_PLEASE
int permute_vecs_iprimme(int *vecs, int n, int *perm_, primme_context ctx) {

   int currentIndex;     /* Index of vector in sorted order                   */
   int sourceIndex;      /* Position of out-of-order vector in original order */
   int destinationIndex; /* Position of out-of-order vector in sorted order   */
   int tempIndex;        /* Used to swap                                      */
   int *perm;            /* A copy of perm_                                   */
   int aux;

   /* Check that perm_ and iwork do not overlap */

   CHKERR(Num_malloc_iprimme(n, &perm, ctx));

   /* Check perm_ is a permutation */

#ifndef NDEBUG
   for (tempIndex=0; tempIndex<n; tempIndex++) perm[tempIndex] = 0;
   for (tempIndex=0; tempIndex<n; tempIndex++) {
      assert(0 <= perm_[tempIndex] && perm_[tempIndex] < n);
      perm[perm_[tempIndex]] = 1;
   }
   for (tempIndex=0; tempIndex<n; tempIndex++) assert(perm[tempIndex] == 1);
#endif

   /* Copy of perm_ into perm, to avoid to modify the input permutation */

   for (tempIndex=0; tempIndex<n; tempIndex++)
      perm[tempIndex] = perm_[tempIndex];

   /* Continue until all vectors are in the sorted order */

   currentIndex = 0;
   while (1) {

      /* Find a vector that does not belong in its original position */
      while ((currentIndex < n) && (perm[currentIndex] == currentIndex)) {
         currentIndex++;
      }

      /* Return if they are in the sorted order */
      if (currentIndex >= n) {
         break;
      }

      /* Copy the vector to a buffer for swapping */
      aux = vecs[currentIndex];

      destinationIndex = currentIndex;
      /* Copy vector perm[destinationIndex] into position destinationIndex */

      while (perm[destinationIndex] != currentIndex) {

         sourceIndex = perm[destinationIndex];
         vecs[destinationIndex] = vecs[sourceIndex];
         tempIndex = perm[destinationIndex];
         perm[destinationIndex] = destinationIndex;
         destinationIndex = tempIndex;
      }

      /* Copy the vector from the buffer to where it belongs */
      vecs[destinationIndex] = aux;
      perm[destinationIndex] = destinationIndex;

      currentIndex++;
   }

   /* Check permutation */
   for (currentIndex=0; currentIndex < n; currentIndex++)
      assert(perm[currentIndex] == currentIndex);

   CHKERR(Num_free_iprimme(perm, ctx));

   return 0;
}
#endif


/******************************************************************************
 * Subroutine Num_compact_vecs - copy certain columns of matrix into another
 *       matrix, i.e., work = vecs(perm). If avoidCopy and perm indices are
 *       consecutive the routine returns a reference in vecs and doesn't copy.
 *            
 *
 * PARAMETERS
 * ---------------------------
 * 
 * vecs        The matrix
 * m           The number of rows of vecs
 * n           The number of columns of to copy
 * ld          The leading dimension of vecs
 * perm        The indices of columns to copy
 * work        The columns are copied to this matrix
 * ldwork      The leading dimension of work
 * avoidCopy   If nonzero, the copy is avoid
 *
 * return      Reference of a matrix with the columns ordered as perm
 *
 ******************************************************************************/

TEMPLATE_PLEASE
SCALAR* Num_compact_vecs_Sprimme(SCALAR *vecs, PRIMME_INT m, int n, 
      PRIMME_INT ld, int *perm, SCALAR *work, PRIMME_INT ldwork,
      int avoidCopy, primme_context ctx) {

   int i;

   if (avoidCopy) {
      for (i=0; i<n-1 && perm[i]+1 == perm[i+1]; i++);
      if (i >= n-1) return &vecs[ld*perm[0]];
   }

   for (i=0; i < n; i++) {
     Num_copy_matrix_Sprimme(&vecs[perm[i] * ld], m, 1, ld, &work[i * ldwork],
                             ld, ctx);
   }
   return work;
}

/******************************************************************************
 * Function Num_scale_matrix - Scale the column of the matrix x into y
 *
 * PARAMETERS
 * ---------------------------
 * x           The source matrix
 * m           The number of rows of x
 * n           The number of columns of x
 * ldx         The leading dimension of x
 * s           The scales
 * y           On output y = x
 * ldy         The leading dimension of y
 *
 * NOTE: x and y *can* overlap
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_scale_matrix_Sprimme(SCALAR *x, PRIMME_INT m, PRIMME_INT n,
      PRIMME_INT ldx, HREAL *s, SCALAR *y, PRIMME_INT ldy, primme_context ctx) {

   PRIMME_INT i;

   assert(m == 0 || n == 0 || (ldx >= m && ldy >= m));

   /* Copy the matrix some columns backward, and other cases */
   /* TODO: assert x and y don't overlap */
   for (i=0; i<n; i++) {
      Num_copy_Sprimme(m, &x[ldx*i], 1, &y[ldy*i], 1, ctx);
      Num_scal_Sprimme(m, s[i], &y[ldy*i], 1, ctx);
   }

   return 0;
}

#endif /* SUPPORTED_TYPE */
