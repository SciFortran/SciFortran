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
 * File: factorize.c
 *
 * Purpose - Functions to factorize and back-solve a hermitian matrix M.
 *           and A'Kinv*B*X
 *  
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/factorize.c"
#endif

#include "numerical.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "update_projection.h"
#include "factorize.h"
#endif

#ifdef SUPPORTED_TYPE

#if defined(USE_HOST) && ((!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF))

/******************************************************************************
 * Function UDUDecompose - This function computes an UDU decomposition of the
 *   matrix M.  See LAPACK routine dsytrf for more information on how the
 *   decomposition is performed.
 *
 *
 * Input Parameters
 * ----------------
 * M  A (numOrthoConst+numEvals) x (numOrthoConst+numEvals) array that contains 
 *    the upper triangular portion of a dimM x dimM hermitian matrix.  
 *    The leading dimension of the array is numOrthoConst+numEvals.
 *
 * dimM  The dimension of the matrix M
 *
 * rwork Real work array of dimension at least dimM.  Optimal size is dimM*NB 
 *       where NB is the block size returned by LAPACK routine ilaenv.
 * 
 *
 * Output Parameters
 * -----------------
 * UDU  Array of dimension dimM x dimM containing the UDU decomposition of M.
 *
 * ipivot  Integer array of length dimM containing pivot mapping
 *
 *
 * Return Value
 * ------------
 * int error code
 ******************************************************************************/

TEMPLATE_PLEASE
int UDUDecompose_Sprimme(SCALAR *M, int ldM, SCALAR *UDU, int ldUDU,
      int *ipivot, int dimM, primme_context ctx) {

   /* Quick return for M with dimension 1 */

   if (dimM == 1) {
      *UDU = *M;
   }
   else {
      /* Copy the upper triangular portion of M into UDU */

      CHKERR(Num_copy_trimatrix_Sprimme(
            M, dimM, dimM, ldM, 0 /* up */, 0, UDU, ldUDU, 0));

      /* Perform the decomposition */

      CHKERR(Num_hetrf_Sprimme("U", dimM, UDU, ldUDU, ipivot, ctx));
   }

   return 0;
}

/******************************************************************************
 * Function UDUSolve - This function solves a dense hermitian linear system
 *   given a right hand side (rhs) and a UDU factorization.
 *
 *
 * Input Parameters
 * ----------------
 * UDU     Two-dimensional of dimension dim and leading dimension dim.
 *         Contains block diagonal and multipliers necessary to construct
 *         the upper triangular matrix U.  See LAPACK routine dsytrf for more
 *         details.
 *
 * ipivot  Permutation array that determines how rows and columns of the
 *         factorization were permuted for stability.
 *
 * dim     The dimension of the linear system
 *
 * rhs     The right hand side of the linear system
 *
 * primme  Structure containing various solver parameters
 *
 *
 * Output Parameters
 * -----------------
 * sol     The solution of the linear system 
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int UDUSolve_Sprimme(SCALAR *UDU, int *ipivot, int dim, SCALAR *rhs, int nrhs,
      int ldrhs, SCALAR *sol, int ldsol, primme_context ctx) {

   if (dim == 1) {
      int i;
      for (i=0; i<nrhs; i++) {
         sol[ldsol*i] = rhs[ldrhs*i]/(*UDU); 
      }
   }
   else {
      CHKERR(Num_copy_matrix_SHprimme(rhs, dim, nrhs, ldrhs, sol, ldsol, ctx));
      CHKERR(Num_hetrs_Sprimme("U", dim, nrhs, UDU, dim, ipivot, sol, ldsol,
                  ctx));
   }

   return 0;

}

#endif /* defined(USE_HOST) && ((!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)) */

/*******************************************************************************
 * Subroutine update_XKinvBX - Updates the matrix M=X'*Kinv*B*X and returns
 *    its factorization.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * X           Matrix with size nLocal x numCols+blockSize
 * ldX         The leading dimension of X
 * KinvBX      Matrix with size nLocal x numCols+blockSize
 * ldKinvBX    The leading dimension of KinvBX
 * numCols     The number of columns that haven't changed
 * blockSize   The number of columns that have changed
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * M           Matrix with size numCols+blockSize with value X'*Y, if it's symmetric
 *             only Z(1:numCols+blockSize,numCols:numCols+blockSize) will be updated
 * ldM         The leading dimension of M,
 * Mfact       Factorization of M
 * ldMfact     The leading dimension of Mfact; if it is zero, it is used numCols+blockSize
 * ipivot      Integer array of length numCols+blockSize containing pivot mapping
 *
 * Return Value
 * ------------
 * int error code
 ******************************************************************************/

TEMPLATE_PLEASE
int update_XKinvBX_Sprimme(SCALAR *X, PRIMME_INT ldX, SCALAR *KinvBX,
      PRIMME_INT ldKinvBX, HSCALAR *M, int ldM, int numCols,
      int blockSize, HSCALAR *Mfact, int ldMfact, int *ipivot,
      primme_context ctx) {

   primme_params *primme = ctx.primme;

   /* Update M. Note that if B=I, then M is Hermitian */

   CHKERR(update_projection_Sprimme(X, ldX, KinvBX, ldKinvBX, M, ldM,
         primme->nLocal, numCols, blockSize,
         primme->massMatrixMatvec ? 0 /* no Hermitian */ : 1 /* Hermitian */,
         ctx));

   int nM = numCols + blockSize; /* dimension of M after updating */

   /* Quick return for M with dimension 0 */

   if (nM == 0) return 0;

   /* Quick return for M with dimension 1 */

   if (nM == 1) {
      *Mfact = *M;
      return 0;
   }

   /* If B == I, M is Hermitian, so M is factorize as LDL^T. Otherwise M is   */
   /* factorize as LU                                                         */

   if (primme->massMatrixMatvec == NULL) {

      /* Copy the upper triangular portion of M into Mfact */

      CHKERR(Num_copy_trimatrix_SHprimme(
            M, nM, nM, ldM, 0 /* up */, 0, Mfact, ldMfact, 0));

      /* Perform the LDL^T decomposition */

      CHKERR(Num_hetrf_SHprimme("U", nM, Mfact, ldMfact, ipivot, ctx));
   }
   else {
      /* Copy M into Mfact */

      CHKERR(Num_copy_matrix_SHprimme(M, nM, nM, ldM, Mfact, ldMfact, ctx));

      /* Perform the LU decomposition */

      CHKERR(Num_getrf_SHprimme(nM, nM, Mfact, ldMfact, ipivot, ctx));
   }

   return 0;
}

#if defined(USE_HOST) && ((!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF))

/******************************************************************************
 * Function MSolve - This function solves a dense hermitian linear system
 *   given a right hand side (rhs) and a UDU factorization.
 *
 *
 * Input Parameters
 * ----------------
 * UDU     Two-dimensional of dimension dim and leading dimension dim.
 *         Contains block diagonal and multipliers necessary to construct
 *         the upper triangular matrix U.  See LAPACK routine dsytrf for more
 *         details.
 *
 * ipivot  Permutation array that determines how rows and columns of the
 *         factorization were permuted for stability.
 *
 * dim     The dimension of the linear system
 *
 * rhs     The right hand side of the linear system
 *
 * primme  Structure containing various solver parameters
 *
 *
 * Output Parameters
 * -----------------
 * sol     The solution of the linear system
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int MSolve_Sprimme(SCALAR *Mfact, int *ipivot, int dim, SCALAR *rhs, int nrhs,
      int ldrhs, SCALAR *sol, int ldsol, primme_context ctx) {

   primme_params *primme = ctx.primme;

   /* Quick return for M with dimension 0 and 1 */

   if (dim == 0) {
      return 0;
   }
   else if (dim == 1) {
      int i;
      for (i = 0; i < nrhs; i++) sol[i * ldsol] = rhs[i * ldrhs] / *Mfact;
      return 0;
   }

   /* Copy rhs into sol */

   CHKERR(Num_copy_matrix_Sprimme(rhs, dim, nrhs, ldrhs, sol, ldsol, ctx));

   /* Solve with the proper decomposition */

   if (primme->massMatrixMatvec == NULL) {
      CHKERR(Num_hetrs_Sprimme(
            "U", dim, nrhs, Mfact, dim, ipivot, sol, ldsol, ctx));
   }
   else {
      CHKERR(Num_getrs_Sprimme(
            "N", dim, nrhs, Mfact, dim, ipivot, sol, ldsol, ctx));
   }

   return 0;
}

#endif /* defined(USE_HOST) && ((!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)) */

#endif /* SUPPORTED_TYPE */
