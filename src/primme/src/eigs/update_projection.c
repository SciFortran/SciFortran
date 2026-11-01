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
 * File: update_projection.c
 *
 * Purpose - Adds blockSize new columns and rows to H.
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/update_projection.c"
#endif


#include <stdlib.h>
#include <assert.h>
#include "numerical.h"
#include "common_eigs.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "auxiliary_eigs.h"
#include "update_projection.h"
#endif

#ifdef SUPPORTED_TYPE

/*******************************************************************************
 * Subroutine update_projection - Z = X'*Y. If Z is a Hermitian matrix 
 *    whose columns will be updated with blockSize vectors, even though space 
 *    for the entire Z is allocated, only the upper triangular portion is 
 *    stored. 
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * X           Matrix with size nLocal x numCols+blockSize
 * ldX         The leading dimension of X
 * Y           Matrix with size nLocal x numCols+blockSize
 * ldY         The leading dimension of Y
 * Z           Matrix with size numCols+blockSize x numCols+blockSize
 * numCols     The number of columns that haven't changed
 * blockSize   The number of columns that have changed
 * rwork       Workspace
 * lrwork      Size of rwork
 * isSymmetric Nonzero if Z is symmetric/Hermitian
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * Z           Matrix with size numCols+blockSize with value X'*Y, if it's symmetric
 *             only Z(1:numCols+blockSize,numCols:numCols+blockSize) will be updated
 * ldZ         The leading dimension of Z
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int update_projection_Sprimme(SCALAR *X, PRIMME_INT ldX, SCALAR *Y,
      PRIMME_INT ldY, HSCALAR *Z, PRIMME_INT ldZ, PRIMME_INT nLocal,
      int numCols, int blockSize, int isSymmetric, primme_context ctx) {

   primme_params *primme = ctx.primme;

   assert(ldX >= nLocal && ldY >= nLocal && ldZ >= numCols+blockSize);

   /* ------------ */
   /* Quick return */
   /* ------------ */

   if (blockSize <= 0) return 0;

   /* --------------------------------------------------------------------- */
   /* Grow Z by blockSize number of rows and columns all at once            */
   /* --------------------------------------------------------------------- */

   int m = numCols+blockSize;
   Num_zero_matrix_SHprimme(&Z[ldZ * numCols], m, blockSize, ldZ, ctx);
   CHKERR(Num_gemm_ddh_Sprimme("C", "N", m, blockSize, nLocal, 1.0, X, ldX,
         &Y[ldY * numCols], ldY, 0.0, &Z[ldZ * numCols], ldZ, ctx));

   /* -------------------------------------------------------------- */
   /* Alternative to the previous call:                              */
   /*    Compute next the additional rows of each new column vector. */
   /*    Only the upper triangular portion is computed and stored.   */
   /* -------------------------------------------------------------- */

   /*
   for (j = numCols; j < numCols+blockSize; j++) {
      Num_gemv_Sprimme("C", primme->nLocal, j-numCols+1, 1.0,
         &X[primme->nLocal*numCols], primme->nLocal, &Y[primme->nLocal*j], 1, 
         0.0, &rwork[maxCols*(j-numCols)+numCols], 1);  
   }
   */

   if (!isSymmetric) {
      Num_zero_matrix_SHprimme(&Z[numCols], blockSize, numCols, ldZ, ctx);
      CHKERR(Num_gemm_ddh_Sprimme("C", "N", blockSize, numCols, nLocal, 1.0,
            &X[ldX * numCols], ldX, Y, ldY, 0.0, &Z[numCols], ldZ, ctx));
   }

   if (primme->numProcs > 1 && isSymmetric) {
      /* --------------------------------------------------------------------- */
      /* Reduce the upper triangular part of the new columns in Z.             */
      /* --------------------------------------------------------------------- */

      HSCALAR *rwork;
      CHKERR(Num_malloc_SHprimme((numCols+blockSize)*blockSize, &rwork, ctx));
      int count;
      CHKERR(Num_copy_trimatrix_compact_SHprimme(
            &Z[ldZ * numCols], m, blockSize, ldZ, numCols, rwork, &count));
      assert(count <= (numCols+blockSize)*blockSize);

      CHKERR(globalSum_SHprimme(rwork, count, ctx));

      CHKERR(Num_copy_compact_trimatrix_SHprimme(
            rwork, m, blockSize, numCols, &Z[ldZ * numCols], ldZ));
      CHKERR(Num_free_SHprimme(rwork, ctx));
   }
   else if (primme->numProcs > 1 && !isSymmetric) {
      /* --------------------------------------------------------------------- */
      /* Reduce Z(:,numCols:end) and Z(numCols:end,:).                         */
      /* --------------------------------------------------------------------- */

      HSCALAR *rwork;
      int count = m * blockSize + blockSize * numCols;
      CHKERR(Num_malloc_SHprimme(count, &rwork, ctx));
      Num_copy_matrix_SHprimme(&Z[ldZ*numCols], m, blockSize, ldZ,
            rwork, m, ctx);
      Num_copy_matrix_SHprimme(&Z[numCols], blockSize, numCols, ldZ,
            &rwork[m*blockSize], blockSize, ctx);

      CHKERR(globalSum_SHprimme(rwork, count, ctx));

      Num_copy_matrix_SHprimme(rwork, m, blockSize, m, &Z[ldZ*numCols],
            ldZ, ctx);
      Num_copy_matrix_SHprimme(&rwork[m*blockSize], blockSize, numCols,
            blockSize, &Z[numCols], ldZ, ctx);
      CHKERR(Num_free_SHprimme(rwork, ctx));
   }

   return 0;
}

#endif /* SUPPORTED_TYPE */
