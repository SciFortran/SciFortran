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
 * File: update_W.c
 *
 * Purpose - Computes A*V(:,0) through A*V(:,blockSize-1) after
 *           V has been expanded by blksze correction vectors.
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/update_W.c"
#endif


#include "numerical.h"
#include "template_normal.h"
#include "common_eigs.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "update_W.h"
#include "auxiliary_eigs.h"
#include "auxiliary_eigs_normal.h"
#include "ortho.h"
#endif

#ifdef SUPPORTED_TYPE

/*******************************************************************************
 * Subroutine update_QR - Computes the QR factorization (A-targetShift*B)*V
 *    updating only the columns nv:nv+blockSize-1 of Q and R.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * BV         B*V
 * nLocal     Number of rows of each vector stored on this node
 * ldBV       The leading dimension of BV
 * W          A*V
 * ldW        The leading dimension of W
 * basisSize  Number of vectors in V
 * blockSize  The current block size
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * Q          The Q factor
 * R          The R factor
 * QtQ        Q'Q
 * fQtQ       The Cholesky factor of QtQ
 ******************************************************************************/

TEMPLATE_PLEASE
int update_Q_Sprimme(SCALAR *BV, PRIMME_INT nLocal, PRIMME_INT ldBV, SCALAR *W,
      PRIMME_INT ldW, SCALAR *Q, PRIMME_INT ldQ, HSCALAR *R, int ldR,
      HSCALAR *QtQ, int ldQtQ, HSCALAR *fQtQ, int ldfQtQ, double targetShift,
      int basisSize, int blockSize, int *nQ, primme_context ctx) {

   int i;

   /* Quick exit */

   if (blockSize <= 0 || R == NULL) return 0;

   assert(ldBV >= nLocal && ldW >= nLocal && ldQ >= nLocal &&
          ldR >= basisSize + blockSize);

   /* Q(:,c) = W(:,c) - BV(:,c)*target for c = basisSize:basisSize+blockSize-1 */

   HEVAL *t;
   CHKERR(KIND(Num_malloc_RHprimme, Num_malloc_SHprimme)(blockSize, &t, ctx));
   for (i=0; i<blockSize; i++) t[i] = targetShift;
   CHKERR(Num_compute_residuals_Sprimme(nLocal, blockSize, t,
         &BV[ldBV * basisSize], ldBV, &W[ldW * basisSize], ldW,
         &Q[ldQ * basisSize], ldQ, ctx));
   CHKERR(KIND(Num_free_RHprimme, Num_free_SHprimme)(t, ctx));

   /* Ortho Q(:,c) for c = basisSize:basisSize+blockSize-1 */

   CHKERR(ortho_block_Sprimme(Q, ldQ, QtQ, ldQtQ, fQtQ, ldfQtQ, R, ldR, *nQ,
         *nQ + blockSize - 1, NULL, 0, 0, NULL, 0, nLocal,
         ctx.primme->maxBasisSize, nQ, ctx));

   /* Zero the lower-left part of R */

   Num_zero_matrix_SHprimme(&R[basisSize], blockSize, basisSize, ldR, ctx);

   return 0;
}

#endif /* SUPPORTED_TYPE */
