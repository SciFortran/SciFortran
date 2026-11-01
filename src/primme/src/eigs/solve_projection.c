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
 * File: solve_H.c
 * 
 * Purpose - Solves the eigenproblem for the matrix V'*A*V.
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/solve_projection.c"
#endif

#include <math.h>
#include <assert.h>
#include "common_eigs.h"
#include "numerical.h"
#include "template_normal.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "solve_projection.h"
#include "ortho.h"
#include "auxiliary_eigs.h"
#include "auxiliary_eigs_normal.h"
#endif

#ifdef SUPPORTED_TYPE

#if defined(USE_HOST) && ((!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF))

/*******************************************************************************
 * Subroutine solve_H - This procedure solves the project problem and return
 *       the projected vectors (hVecs) and values (hVals) in the order according
 *       to primme.target.
 *        
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * H              The matrix V'*A*V
 * basisSize      The dimension of H, R, QtV and hU
 * ldH            The leading dimension of H
 * VtBV           The matrix V'*B*V
 * ldVtBV         The leading dimension of VtBV
 * R              The factor R for the QR decomposition of (A - target*I)*V
 * ldR            The leading dimension of R
 * QtV            Q'*V
 * ldQtV          The leading dimension of QtV
 * numConverged   Number of eigenvalues converged to determine ordering shift
 * lrwork         Length of the work array rwork
 * primme         Structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * hU             The left singular vectors of R or the eigenvectors of QtV/R
 * ldhU           The leading dimension of hU
 * hVecs          The coefficient vectors such as V*hVecs will be the Ritz vectors
 * ldhVecs        The leading dimension of hVecs
 * hVals          The Ritz values
 * hSVals         The singular values of R
 * rwork          Workspace
 * iwork          Workspace in integers
 *
 * Return Value
 * ------------
 * int -  0 upon successful return
 *     - -1 Num_dsyev/zheev was unsuccessful
 ******************************************************************************/

TEMPLATE_PLEASE
int solve_H_Sprimme(SCALAR *H, int basisSize, int ldH, SCALAR *VtBV, int ldVtBV,
      SCALAR *R, int ldR, SCALAR *QtV, int ldQtV, SCALAR *QtQ, int ldQtQ,
      SCALAR *hU, int ldhU, SCALAR *hVecs, int ldhVecs, EVAL *hVals,
      REAL *hSVals, int numConverged, primme_context ctx) {

   int i;

   /* In parallel (especially with heterogeneous processors/libraries) ensure */
   /* that every process has the same hVecs and hU. Only processor 0 solves   */
   /* the projected problem and broadcasts the resulting matrices to the rest */

   if (ctx.primme->procID == 0) {
      switch (ctx.primme->projectionParams.projection) {
         case primme_proj_RR:
            CHKERR(solve_H_RR_Sprimme(H, ldH, VtBV, ldVtBV, hVecs, ldhVecs,
                     hVals, basisSize, numConverged, ctx))
            break;

         case primme_proj_harmonic:
            CHKERR(solve_H_Harm_Sprimme(H, ldH, QtV, ldQtV, R, ldR, QtQ, ldQtQ,
                  VtBV, ldVtBV, hVecs, ldhVecs, hU, ldhU, hVals, basisSize,
                  numConverged, ctx));
            break;

         case primme_proj_refined:
            CHKERR(solve_H_Ref_Sprimme(H, ldH, hVecs, ldhVecs, hU, ldhU, hSVals,
                  R, ldR, QtQ, ldQtQ, VtBV, ldVtBV, hVals, basisSize,
                  numConverged, ctx));
            break;

         default:
            assert(0);
      }
   }

   /* Broadcast hVecs, hU, hVals, hSVals */

   CHKERR(solve_H_brcast_Sprimme(basisSize, hU, ldhU, hVecs, ldhVecs, hVals,
            hSVals, ctx));
 
   /* Return memory requirements */

   if (H == NULL) {
      return 0;
   }

   /* -------------------------------------------------------- */
   /* Update the leftmost and rightmost Ritz values ever seen  */
   /* -------------------------------------------------------- */
   for (i=0; i<basisSize; i++) {
     ctx.primme->stats.estimateMinEVal =
         min(ctx.primme->stats.estimateMinEVal, EVAL_REAL_PART(hVals[i]));
     ctx.primme->stats.estimateMaxEVal =
         max(ctx.primme->stats.estimateMaxEVal, EVAL_REAL_PART(hVals[i])); 
     ctx.primme->stats.estimateLargestSVal =
         max(ctx.primme->stats.estimateLargestSVal, EVAL_ABS(hVals[i])); 
   }

   return 0;
}


/*******************************************************************************
 * Subroutine solve_H_RR - This procedure solves the eigenproblem for the
 *            matrix H.
 *        
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * H              The matrix V'*A*V
 * basisSize      The dimension of H, R, hU
 * ldH            The leading dimension of H
 * VtBV           The matrix V'*B*V
 * ldVtBV         The leading dimension of VtBV
 * numConverged   Number of eigenvalues converged to determine ordering shift
 * lrwork         Length of the work array rwork
 * primme         Structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * hVecs          The eigenvectors of H or the right singular vectors
 * ldhVecs        The leading dimension of hVecs
 * hVals          The Ritz values
 * hSVals         The singular values of R
 * rwork          Workspace
 * iwork          Workspace in integers
 *
 * Return Value
 * ------------
 * int -  0 upon successful return
 *     - -1 Num_dsyev/zheev was unsuccessful
 ******************************************************************************/

TEMPLATE_PLEASE
int solve_H_RR_Sprimme(SCALAR *H, int ldH, SCALAR *VtBV, int ldVtBV,
      SCALAR *hVecs, int ldhVecs, EVAL *hVals, int basisSize, int numConverged,
      primme_context ctx) {

   /* Some LAPACK implementations don't like zero-size matrices */
   if (basisSize == 0) return 0;

#ifdef USE_HERMITIAN
   primme_params *primme = ctx.primme;
   int i, j; /* Loop variables    */
   int index;

   /* ------------------------------------------------------------------- */
   /* Copy the upper triangular portion of H into hvecs.  We need to do   */
   /* this since DSYEV overwrites the input matrix with the eigenvectors. */  
   /* Note that H is maxBasisSize-by-maxBasisSize and the basisSize-by-   */
   /* basisSize submatrix of H is copied into hvecs.                      */
   /* ------------------------------------------------------------------- */

   if (primme->target != primme_largest) {
      CHKERR(Num_copy_trimatrix_Sprimme(H, basisSize, basisSize, ldH,
            0 /* upper part */, 0, hVecs, ldhVecs, 0));
   } else { /* (primme->target == primme_largest) */
      for (j = 0; j < basisSize; j++) {
         for (i = 0; i <= j; i++) hVecs[ldhVecs * j + i] = -H[ldH * j + i];
      }
   }

   CHKERR(Num_hegv_Sprimme("V", "U", basisSize, hVecs, ldhVecs, VtBV, ldVtBV,
                           hVals, ctx));

   /* ---------------------------------------------------------------------- */
   /* ORDER the eigenvalues and their eigenvectors according to the desired  */
   /* target:  smallest/Largest or interior closest abs/leq/geq to a shift   */
   /* ---------------------------------------------------------------------- */

   if (primme->target == primme_smallest) 
      return 0;

   if (primme->target == primme_largest) {
      for (i = 0; i < basisSize; i++) {
         hVals[i] = -hVals[i];
      }
   }
   else { 
      int *permu;  /* permutation of the pairs */
      CHKERR(Num_malloc_iprimme(basisSize, &permu, ctx));

      /* ---------------------------------------------------------------- */
      /* Select the interior shift. Use the first unlocked shift, and not */
      /* higher ones, even if some eigenpairs in the basis are converged. */
      /* Then order the ritz values based on the closeness to the shift   */
      /* from the left, from right, or in absolute value terms            */
      /* ---------------------------------------------------------------- */

      /* TODO: order properly when numTargetShifts > 1 */

      double targetShift = 
        primme->targetShifts[min(primme->numTargetShifts-1, numConverged)];

      if (primme->target == primme_closest_geq) {
   
         /* ---------------------------------------------------------------- */
         /* find hVal closest to the right of targetShift, i.e., closest_geq */
         /* ---------------------------------------------------------------- */
         for (j=0;j<basisSize;j++) 
              if (hVals[j]>=targetShift) break;
           
         /* figure out this ordering */
         index = 0;
   
         for (i=j; i<basisSize; i++) {
            permu[index++]=i;
         }
         for (i=0; i<j; i++) {
            permu[index++]=i;
         }
      }
      else if (primme->target == primme_closest_leq) {
         /* ---------------------------------------------------------------- */
         /* find hVal closest_leq to targetShift                             */
         /* ---------------------------------------------------------------- */
         for (j=basisSize-1; j>=0 ;j--) 
             if (hVals[j]<=targetShift) break;
           
         /* figure out this ordering */
         index = 0;
   
         for (i=j; i>=0; i--) {
            permu[index++]=i;
         }
         for (i=basisSize-1; i>j; i--) {
            permu[index++]=i;
         }
      }
      else if (primme->target == primme_closest_abs) {

         /* ---------------------------------------------------------------- */
         /* find hVal closest but geq than targetShift                       */
         /* ---------------------------------------------------------------- */
         for (j=0;j<basisSize;j++) 
             if (hVals[j]>=targetShift) break;

         i = j-1;
         index = 0;
         while (i>=0 && j<basisSize) {
            if (fabs(hVals[i]-targetShift) < fabs(hVals[j]-targetShift)) 
               permu[index++] = i--;
            else 
               permu[index++] = j++;
         }
         if (i<0) {
            for (i=j;i<basisSize;i++) 
                    permu[index++] = i;
         }
         else if (j>=basisSize) {
            for (j=i;j>=0;j--)
                    permu[index++] = j;
         }
      }
      else if (primme->target == primme_largest_abs) {

         j = 0;
         i = basisSize-1;
         index = 0;
         while (i>=j) {
            if (fabs(hVals[i]-targetShift) > fabs(hVals[j]-targetShift)) 
               permu[index++] = i--;
            else 
               permu[index++] = j++;
         }

      }

      /* ---------------------------------------------------------------- */
      /* Reorder hVals and hVecs according to the permutation             */
      /* ---------------------------------------------------------------- */
      permute_vecs_Rprimme(hVals, 1, basisSize, 1, permu, ctx);
      permute_vecs_Sprimme(hVecs, basisSize, basisSize, ldhVecs, permu, ctx);

      CHKERR(Num_free_iprimme(permu, ctx));
   }

   return 0;   
#else

   /* Copy H into T */

   SCALAR *T;
   CHKERR(Num_malloc_Sprimme(basisSize * basisSize, &T, ctx));
   CHKERR(Num_copy_matrix_Sprimme(
         H, basisSize, basisSize, ldH, T, basisSize, ctx));

   /* Cholesky factorization of VtBV = U'U, and apply T = U'\T/U */

   SCALAR *U=NULL; /* Cholesky factor of VtBV */
   if (VtBV) {
      CHKERR(Num_malloc_Sprimme(basisSize*basisSize, &U, ctx));
      Num_copy_matrix_Sprimme(
            VtBV, basisSize, basisSize, ldVtBV, U, basisSize, ctx);
      CHKERR(Num_potrf_Sprimme("U", basisSize, U, basisSize, NULL, ctx));
      CHKERR(Num_trsm_Sprimme("R", "U", "N", "N", basisSize, basisSize, 1.0,
            U, basisSize, T, basisSize, ctx));
      CHKERR(Num_trsm_Sprimme("L", "U", "C", "N", basisSize, basisSize, 1.0,
            U, basisSize, T, basisSize, ctx));
   }


   /* Compute the Schur factorization of T */

   CHKERR(Num_gees_Sprimme(
         "V", basisSize, T, basisSize, hVals, hVecs, ldhVecs, ctx));
   CHKERR(Num_free_Sprimme(T, ctx));

   /* Apply hVecs = U\hVecs */

   if (VtBV) {
      CHKERR(Num_trsm_Sprimme("L", "U", "N", "N", basisSize, basisSize, 1.0,
            U, basisSize, hVecs, ldhVecs, ctx));
      CHKERR(Num_free_Sprimme(U, ctx));
   }

   /* Order the eigenpairs following primme.target */

   int i;
   int *perm;
   CHKERR(Num_malloc_iprimme(basisSize, &perm, ctx));
   for (i = 0; i < basisSize; i++) perm[i] = i;

   for (i = 1; i < basisSize; i++) {
      CHKERR(insertionSort_Sprimme(hVals[i], hVals, 0.0, NULL, 0, NULL, perm, i,
            numConverged, ctx.primme));
   }
   permute_vecs_Sprimme(hVecs, basisSize, basisSize, ldhVecs, perm, ctx);

   CHKERR(Num_free_iprimme(perm, ctx));

   return 0;
#endif /* USE_HERMITIAN */
}

/*******************************************************************************
 * Subroutine solve_H_Harm - This procedure implements the harmonic extraction
 *    in a novelty way. In standard harmonic the next eigenproblem is solved:
 *       V'*(A-s*I)'*(A-s*I)*V*X = V'*(A-s*I)'*V*X*L,
 *    where (L_{i,i},X_i) are the harmonic-Ritz pairs. In practice, it is
 *    computed (A-s*I)*V = Q*R and it is solved instead:
 *       R*X = Q'*V*X*L,
 *    which is a generalized non-Hermitian problem. Instead of dealing with
 *    complex solutions, which are unnatural in context of Hermitian problems,
 *    we propose the following. Note that,
 *       (A-s*I)*V = Q*R -> Q'*V*inv(R) = Q'*inv(A-s*I)*Q.
 *    And note that Q'*V*inv(R) is Hermitian if A is, and also that
 *       Q'*V*inv(R)*Y = Y*inv(L) ->  Q'*V*X*L = R*X,
 *    with Y = R*X. So this routine computes X by solving the Hermitian problem
 *    Q'*V*inv(R).
 *        
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * H             The matrix V'*A*V
 * ldH           The leading dimension of H
 * R             The R factor for the QR decomposition of (A - target*I)*V
 * ldR           The leading dimension of R
 * basisSize     Current size of the orthonormal basis V
 * lrwork        Length of the work array rwork
 * primme        Structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * hVecs         The orthogonal basis of inv(R) * eigenvectors of QtV/R
 * ldhVecs       The leading dimension of hVecs
 * hU            The eigenvectors of QtV/R
 * ldhU          The leading dimension of hU
 * hVals         The Ritz values of the vectors in hVecs
 * rwork         Workspace
 *
 * Return Value
 * ------------
 * int -  0 upon successful return
 *     - -1 Num_dsyev/zheev was unsuccessful
 ******************************************************************************/

STATIC int solve_H_Harm_Sprimme(SCALAR *H, int ldH, SCALAR *QtV, int ldQtV,
      SCALAR *R, int ldR, SCALAR *QtQ, int ldQtQ, SCALAR *VtBV, int ldVtBV,
      SCALAR *hVecs, int ldhVecs, SCALAR *hU, int ldhU, EVAL *hVals,
      int basisSize, int numConverged, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, ret;
   double *oldTargetShifts, zero=0.0;
   primme_target oldTarget;

   (void)numConverged; /* unused parameter */

   /* Some LAPACK implementations don't like zero-size matrices */
   if (basisSize == 0) return 0;

   SCALAR *rwork;
   CHKERR(Num_malloc_Sprimme(basisSize*basisSize, &rwork, ctx));

   /* Factorize R */

   SCALAR *fR = rwork;
   int *pivots;
   CHKERR(Num_malloc_iprimme(basisSize, &pivots, ctx));
   Num_copy_matrix_Sprimme(R, basisSize, basisSize, ldR, fR, basisSize, ctx);
   CHKERR(Num_getrf_Sprimme(basisSize, basisSize, fR, basisSize, pivots, ctx));

   /* QAQ = QtV*inv(R) = R'\QtV' */

   Num_copy_matrix_conj_Sprimme(
         QtV, basisSize, basisSize, ldQtV, hVecs, ldhVecs, ctx);
   CHKERR(Num_getrs_Sprimme("C", basisSize, basisSize, fR, basisSize, pivots,
         hVecs, ldhVecs, ctx));

   /* Compute eigenpairs of (Q'AQ, Q'Q) */

   oldTargetShifts = primme->targetShifts;
   oldTarget = primme->target;
   primme->targetShifts = &zero;
   switch(primme->target) {
      case primme_closest_geq:
         primme->target = primme_largest;
         break;
      case primme_closest_leq:
         primme->target = primme_smallest;
         break;
      case primme_closest_abs:
         primme->target = primme_largest_abs;
         break;
      default:
         assert(0);
   }
   ret = solve_H_RR_Sprimme(
         hVecs, ldhVecs, QtQ, ldQtQ, hVecs, ldhVecs, hVals, basisSize, 0, ctx);
   primme->targetShifts = oldTargetShifts;
   primme->target = oldTarget;
   CHKERRM(ret, ret, "Error calling solve_H_RR_Sprimme");

   Num_copy_matrix_Sprimme(hVecs, basisSize, basisSize, ldhVecs, hU, ldhU, ctx);

   /* Transfer back the eigenvectors to V, hVecs = R\hVecs */

   CHKERR(Num_getrs_Sprimme("N", basisSize, basisSize, fR, basisSize, pivots,
         hVecs, ldhVecs, ctx));
   CHKERR(Num_free_iprimme(pivots, ctx));
   CHKERR(Bortho_local_SHprimme(hVecs, ldhVecs, NULL, 0, 0, basisSize - 1, NULL,
         0, 0, basisSize, VtBV, ldVtBV, primme->iseed, ctx));

   /* Compute Rayleigh quotient lambda_i = x_i'*H*x_i */

   CHKERR(Num_hemm_Sprimme("L", "U", basisSize, basisSize, 1.0, H,
      ldH, hVecs, ldhVecs, 0.0, rwork, basisSize, ctx));

   for (i=0; i<basisSize; i++) {
      hVals[i] = KIND(REAL_PART, )(Num_dot_Sprimme(
            basisSize, &hVecs[ldhVecs * i], 1, &rwork[basisSize * i], 1, ctx));
   }
   CHKERR(Num_free_Sprimme(rwork, ctx));

   return 0;
}

/*******************************************************************************
 * Subroutine solve_H_Ref - This procedure solves the singular value
 *            decomposition of matrix R
 *        
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * H             The matrix V'*A*V
 * ldH           The leading dimension of H
 * R             The R factor for the QR decomposition of (A - target*I)*V
 * ldR           The leading dimension of R
 * basisSize     Current size of the orthonormal basis V
 * lrwork        Length of the work array rwork
 * primme        Structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * hVecs         The right singular vectors of R
 * ldhVecs       The leading dimension of hVecs
 * hU            The left singular vectors of R
 * ldhU          The leading dimension of hU
 * hSVals        The singular values of R
 * hVals         The Ritz values of the vectors in hVecs
 * rwork         Workspace
 *
 * Return Value
 * ------------
 * int -  0 upon successful return
 *     - -1 was unsuccessful
 ******************************************************************************/

STATIC int solve_H_Ref_Sprimme(SCALAR *H, int ldH, SCALAR *hVecs, int ldhVecs,
      SCALAR *hU, int ldhU, REAL *hSVals, SCALAR *R, int ldR, SCALAR *QtQ,
      int ldQtQ, SCALAR *VtBV, int ldVtBV, EVAL *hVals, int basisSize,
      int targetShiftIndex, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, j; /* Loop variables    */

   (void)targetShiftIndex; /* unused parameter */

   /* Some LAPACK implementations don't like zero-size matrices */
   if (basisSize == 0) return 0;

   /* Copy R into hVecs */
   Num_copy_matrix_Sprimme(R, basisSize, basisSize, ldR, hVecs, ldhVecs, ctx);

   if (QtQ) {
      /* Factorize QtQ */
      Num_copy_matrix_Sprimme(
            QtQ, basisSize, basisSize, ldQtQ, hU, ldhU, ctx);
      CHKERR(Num_potrf_Sprimme("U", basisSize, hU, ldhU, NULL, ctx));
      CHKERR(Num_trmm_Sprimme("L", "U", "N", "N", basisSize, basisSize, 1.0, hU,
            ldhU, hVecs, ldhVecs, ctx));
   }

   SCALAR *U_VtBV=NULL; /* Cholesky factor of VtBV */
   if (VtBV) {
      CHKERR(Num_malloc_Sprimme(basisSize*basisSize, &U_VtBV, ctx));
      Num_copy_matrix_Sprimme(
            VtBV, basisSize, basisSize, ldVtBV, U_VtBV, basisSize, ctx);
      CHKERR(Num_potrf_Sprimme("U", basisSize, U_VtBV, basisSize, NULL, ctx));
      CHKERR(Num_trsm_Sprimme("R", "U", "N", "N", basisSize, basisSize, 1.0,
            U_VtBV, basisSize, hVecs, basisSize, ctx));
   }

   /* Note gesvd returns transpose(V) rather than V and sorted in descending */
   /* order of the singular values */

   CHKERR(Num_gesvd_Sprimme("S", "O", basisSize, basisSize, hVecs, ldhVecs,
         hSVals, hU, ldhU, hVecs, ldhVecs, ctx));

   /* Transpose back V */

   SCALAR *rwork;
   CHKERR(Num_malloc_Sprimme((size_t)basisSize * basisSize, &rwork, ctx));
   for (j = 0; j < basisSize; j++) {
      for (i = 0; i < basisSize; i++) {
         rwork[basisSize * j + i] = CONJ(hVecs[ldhVecs * i + j]);
      }
   }
   Num_copy_matrix_Sprimme(
         rwork, basisSize, basisSize, basisSize, hVecs, ldhVecs, ctx);

   if (VtBV) {
      CHKERR(Num_trsm_Sprimme("L", "U", "N", "N", basisSize, basisSize, 1.0,
            U_VtBV, basisSize, hVecs, ldhVecs, ctx));
      CHKERR(Num_free_Sprimme(U_VtBV, ctx));
   }

   /* Rearrange V, hSVals and hU in ascending order of singular value   */
   /* if target is not largest abs.                                     */

   if (primme->target == primme_closest_abs ||
         primme->target == primme_closest_leq ||
         primme->target == primme_closest_geq) {
      int *perm;
      CHKERR(Num_malloc_iprimme(basisSize, &perm, ctx));
      for (i = 0; i < basisSize; i++) perm[i] = basisSize - 1 - i;
      permute_vecs_Rprimme(hSVals, 1, basisSize, 1, perm, ctx);
      permute_vecs_Sprimme(hVecs, basisSize, basisSize, ldhVecs, perm, ctx);
      permute_vecs_Sprimme(hU, basisSize, basisSize, ldhU, perm, ctx);
      CHKERR(Num_free_iprimme(perm, ctx));
   }

   /* compute Rayleigh quotient lambda_i = x_i'*H*x_i */

   CHKERR(Num_hemm_Sprimme("L", "U", basisSize, basisSize, 1.0, H, ldH, hVecs,
         ldhVecs, 0.0, rwork, basisSize, ctx));

   for (i=0; i<basisSize; i++) {
      hVals[i] = KIND(REAL_PART, )(Num_dot_Sprimme(
            basisSize, &hVecs[ldhVecs * i], 1, &rwork[basisSize * i], 1, ctx));
   }
   CHKERR(Num_free_Sprimme(rwork, ctx));

   return 0;
}

/*******************************************************************************
 * Subroutine solve_H_brcast - This procedure broadcast the solution of the
 *       projected problem (hVals, hSVals, hVecs, hU) from process 0 to the rest.
 *
 * NOTE: the optimal implementation will use an user-defined broadcast function.
 *       To ensure backward compatibility, we used globalSum instead.
 * 
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * basisSize      The dimension of hVecs, hU, hVals and hSVals
 * lrwork         Length of the work array rwork
 * primme         Structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * hU             The left singular vectors of R or the eigenvectors of QtV/R
 * ldhU           The leading dimension of hU
 * hVecs          The coefficient vectors such as V*hVecs will be the Ritz vectors
 * ldhVecs        The leading dimension of hVecs
 * hVals          The Ritz values
 * hSVals         The singular values of R
 * rwork          Workspace
 *
 * Return Value
 * ------------
 * error code                          
 ******************************************************************************/

STATIC int solve_H_brcast_Sprimme(int basisSize, SCALAR *hU, int ldhU,
                                  SCALAR *hVecs, int ldhVecs, EVAL *hVals,
                                  REAL *hSVals, primme_context ctx) {


   SCALAR *rwork0;             /* next SCALAR free */
   SCALAR *rwork;
   const size_t c = sizeof(SCALAR)/sizeof(REAL);

   /* Quick exit */

   if (basisSize <= 0) return 0;

   /* Allocate memory */

   int n=0;  /* number of SCALAR packed */
   if (hVecs) n += basisSize*basisSize;
   if (hU) n += basisSize*basisSize;
   if (hVals) n += KIND((basisSize + c - 1) / c, basisSize);
   if (hSVals) n += (basisSize + c-1)/c;
   CHKERR(Num_malloc_Sprimme(n, &rwork, ctx));
   rwork0 = rwork;

   /* Pack hVecs */

   if (hVecs) {
      if (ctx.primme->procID == 0) {
         Num_copy_matrix_Sprimme(hVecs, basisSize, basisSize, ldhVecs, rwork0,
               basisSize, ctx);
      }
      rwork0 += basisSize*basisSize;
   }

   /* Pack hU */

   if (hU) {
      if (ctx.primme->procID == 0) {
         Num_copy_matrix_Sprimme(hU, basisSize, basisSize, ldhU, rwork0,
               basisSize, ctx);
      }
      rwork0 += basisSize*basisSize;
   }

   /* Pack hVals */

   if (hVals) {
      if (ctx.primme->procID == 0) {
#ifdef USE_HERMITIAN
         //  When complex, avoid to reduce with an uninitialized value
         rwork0[(basisSize + c-1)/c-1] = 0.0; 
         Num_copy_matrix_Rprimme(hVals, basisSize, 1, basisSize, (REAL*)rwork0,
               basisSize, ctx);
         rwork0 += (basisSize + c - 1) / c;
#else
         Num_copy_matrix_Sprimme(
               hVals, basisSize, 1, basisSize, rwork0, basisSize, ctx);
         rwork0 += basisSize;
#endif
      }
   }

   /* Pack hSVals */

   if (hSVals) {
      if (ctx.primme->procID == 0) {
         rwork0[(basisSize + c-1)/c-1] = 0.0; /* When complex, avoid to reduce with an   */
                                              /* uninitialized value                     */
         Num_copy_matrix_Rprimme(hSVals, basisSize, 1, basisSize, (REAL*)rwork0,
               basisSize, ctx);
      }
      rwork0 += (basisSize + c-1)/c;
   }

   /* Perform the broadcast */

   CHKERR(broadcast_Sprimme(rwork, n, ctx));
   rwork0 = rwork;

   /* Unpack hVecs */

   if (hVecs) {
      Num_copy_matrix_Sprimme(rwork0, basisSize, basisSize, basisSize, hVecs,
            ldhVecs, ctx);
      rwork0 += basisSize*basisSize;
   }

   /* Unpack hU */

   if (hU) {
      Num_copy_matrix_Sprimme(rwork0, basisSize, basisSize, basisSize, hU,
            ldhU, ctx);
      rwork0 += basisSize*basisSize;
   }

   /* Unpack hVals */

   if (hVals) {
#ifdef USE_HERMITIAN
      Num_copy_matrix_Rprimme((REAL*)rwork0, basisSize, 1, basisSize, hVals,
            basisSize, ctx);
      rwork0 += (basisSize + c-1)/c;
#else
      Num_copy_matrix_Sprimme(
            rwork0, basisSize, 1, basisSize, hVals, basisSize, ctx);
      rwork0 += basisSize;
#endif
   }

   /* Unpack hSVals */

   if (hSVals) {
      Num_copy_matrix_Rprimme((REAL*)rwork0, basisSize, 1, basisSize, hSVals,
               basisSize, ctx);
      rwork0 += (basisSize + c-1)/c;
   }

   CHKERR(Num_free_Sprimme(rwork, ctx));

   return 0;
}

/*******************************************************************************
 * Function prepare_vecs - This subroutine checks that the
 *    conditioning of the coefficient vectors are good enough to converge
 *    with the requested accuracy. For now refined extraction is the only one that
 *    may present problems: two similar singular values in the projected problem
 *    may correspond to distinct eigenvalues in the original problem. If that is
 *    the case, the singular vector may have components of both eigenvectors,
 *    which prevents the residual norm be lower than some degree. Don't dealing
 *    with this may lead into stagnation.
 *
 *    It is checked that the next upper bound about the angle of the right
 *    singular vector v of A and the right singular vector vtilde of A+E,
 *
 *      sin(v,vtilde) <= sqrt(2)*||E||/sval_gap <= sqrt(2)*||A||*machEps/sval_gap,
 *
 *    is less than the upper bound about the angle of exact eigenvector u and
 *    the approximate eigenvector utilde,
 *
 *      sin(u,utilde) <= ||r||/eval_gap <= ||A||*eps/eval_gap.
 *
 *    (see pp. 211 in Matrix Algorithms vol. 2 Eigensystems, G. W. Steward).
 *
 *    If the inequality doesn't hold, do Rayleigh-Ritz onto the subspace
 *    spanned by both vectors.
 *
 *    we have found cases where this is not enough or the performance improves
 *    if Rayleigh-Ritz is also done when the candidate vector has a small
 *    angle with the last vector in V and when the residual norm is larger than
 *    the singular value.
 *
 *    When only one side of the shift is targeted (primme_closest_leq/geq), we
 *    allow to take eigenvalue of the other side but close to the shift. In the
 *    current heuristic they shouldn't be farther than the smallest residual
 *    norm in the block. This heuristic obtained good results in solving the
 *    augmented problem with shifts from solving the normal equations.
 *
 * NOTE: this function assumes hSVals are arranged in increasing order.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * basisSize    Projected problem size
 * i0           Index of the first pair to check
 * blockSize    Number of candidates wanted
 * H            The matrix V'*A*V
 * ldH          The leading dimension of H
 * hVals        The Ritz values
 * hSVals       The singular values of R
 * hVecs        The coefficient vectors
 * ldhVecs      The leading dimension of hVecs
 * targetShiftIndex The target shift used in (A - targetShift*B) = Q*R
 * arbitraryVecs The number of vectors modified (input/output)
 * smallestResNorm The smallest residual norm in the block
 * flags        Array indicating the convergence of the Ritz vectors
 * RRForAll     If false compute Rayleigh-Ritz only in clusters with
 *              candidates. If true, compute it in every cluster.
 * machEps      Machine precision
 * rworkSize    The length of rwork
 * rwork        Workspace
 * iwork        Integer workspace
 * primme       Structure containing various solver parameters
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int prepare_vecs_Sprimme(int basisSize, int i0, int blockSize, SCALAR *H,
      int ldH, EVAL *hVals, REAL *hSVals, SCALAR *hVecs, int ldhVecs,
      int targetShiftIndex, int *arbitraryVecs, double smallestResNorm,
      int *flags, int RRForAll, SCALAR *hVecsRot, int ldhVecsRot,
      primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, j, k;         /* Loop indices */
   int candidates;      /* Number of eligible pairs */
   int someCandidate;   /* If there is an eligible pair in the cluster */
   double aNorm, eps;

   /* Quick exit */

   if (primme->projectionParams.projection != primme_proj_refined ||
         basisSize == 0) {
      /* Find next candidates, starting from iev(*blockSize)+1 */

      if (flags) {
         for (i = i0, j = 0; i < basisSize && j < blockSize; i++) {
            if (flags[i] == UNCONVERGED) j++;
         }
      }

      return 0;
   }

   /* Quick exit */

   if (blockSize == 0) {
      return 0;
   }

   aNorm = (primme->aNorm <= 0.0) ?
      primme->stats.estimateLargestSVal : primme->aNorm;

   /* Before the first eigenpair converges, there's no information about the  */
   /* requested tolerance. In that case eps is set as ten times less than the */
   /* the current smallest residual norm, if smallestResNorm provides that    */
   /* information.                                                            */
   eps = primme->stats.maxConvTol > 0.0 ? primme->stats.maxConvTol : (
         smallestResNorm < HUGE_VAL ? smallestResNorm/10.0 : 0.0);
   double eps_orth;
   CHKERR(machineEpsOrth_Sprimme(&eps_orth, ctx));
   /* NOTE: the constant 6.28 is needed to pass                               */
   /* testi-100-LOBPCG_OrthoBasis-2-primme_closest_abs-primme_proj_refined.F  */
   eps = max(6.28 * eps_orth, eps);

   for (candidates=0, i=min(*arbitraryVecs,basisSize), j=i0;
         j < basisSize && candidates < blockSize; ) {

      double ip;
      /* -------------------------------------------------------------------- */
      /* Count all eligible values (candidates) from j up to i.               */
      /* -------------------------------------------------------------------- */

      for ( ; j < i; j++) {
         if (!flags || flags[j] == UNCONVERGED)
            candidates++;
      }
     
      if (candidates >= blockSize) break;
 
      /* -------------------------------------------------------------------- */
      /* Find the first i-th vector i>j with enough good conditioning, ie.,   */
      /* the singular value is separated enough from the rest (see the header */
      /* comment in this function). Also check if there is an unconverged     */
      /* value in the block.                                                  */
      /* -------------------------------------------------------------------- */

      for (i=j+1, someCandidate=0, ip=0.0; i<basisSize; i++) {

         /* Check that this approximation:                                    */
         /* sin singular vector: max(hSVals)*machEps/(hSvals[i]-hSVals[i+1])  */
         /* is less than the next one:                                        */
         /* sin eigenvector    : aNorm*eps/(hVals[i]-hVals[i+1]).             */
         /* Also try to include enough coefficient vectors into the cluster   */
         /* so that the cluster is close the last included vectors into the   */
         /* basis.                                                            */
         /* TODO: check the angle with all vectors added to the basis in the  */
         /* previous iterations; for now only the last one is considered.     */
         /* NOTE: we don't want to check hVecs(end,i) just after restart, so  */
         /* we don't use the value when it is zero.                           */

         double minDiff = sqrt(2.0)*hSVals[basisSize-1]*MACHINE_EPSILON/
            (aNorm*eps/EVAL_ABS(hVals[i]-hVals[i-1]));
         double ip0 = ABS(hVecs[(i-1)*ldhVecs+basisSize-1]);
         double ip1 = ((ip += ip0*ip0) != 0.0) ? ip : HUGE_VAL;

         someCandidate = 1;

         if (fabs(hSVals[i]-hSVals[i-1]) >= minDiff
               && (smallestResNorm >= HUGE_VAL
                  || sqrt(ip1) >= smallestResNorm/aNorm/3.16)) 
            break;
      }
      i = min(i, basisSize);

      /* ----------------------------------------------------------------- */
      /* If the cluster j:i-1 is larger than one vector and there is some  */
      /* unconverged pair in there, compute the approximate eigenvectors   */
      /* with Rayleigh-Ritz. If RRForAll do also this when there is no     */
      /* candidate in the cluster.                                         */
      /* ----------------------------------------------------------------- */

      if (i-j > 1 && (someCandidate || RRForAll)) {
         SCALAR *aH, *ahVecs;
         int aBasisSize = i-j;
         CHKERR(Num_malloc_Sprimme((size_t)basisSize*aBasisSize, &aH, ctx));
         ahVecs = &hVecsRot[ldhVecsRot*j+j];

         /* Zero hVecsRot(:,arbitraryVecs:i-1) */
         Num_zero_matrix_Sprimme(&hVecsRot[ldhVecsRot*(*arbitraryVecs)],
               primme->maxBasisSize, i-*arbitraryVecs, ldhVecsRot, ctx);

         /* hVecsRot(:,arbitraryVecs:i-1) = I */
         for (k=*arbitraryVecs; k<i; k++)
            hVecsRot[ldhVecsRot*k+k] = 1.0;
 
         /* aH = hVecs(:,j:i-1)'*H*hVecs(:,j:i-1) */
         CHKERR(compute_submatrix_Sprimme(&hVecs[ldhVecs * j], aBasisSize,
               ldhVecs, H, basisSize, ldH,
               KIND(1 /* Hermitian */, 0 /* non-Hermitian */), aH, aBasisSize,
               ctx));

         /* Compute and sort eigendecomposition aH*ahVecs = ahVecs*diag(hVals(j:i-1)) */
         CHKERR(solve_H_RR_Sprimme(aH, aBasisSize, NULL, 0, ahVecs, ldhVecsRot,
               &hVals[j], aBasisSize, targetShiftIndex, ctx));

         /* hVecs(:,j:i-1) = hVecs(:,j:i-1)*ahVecs */
         Num_zero_matrix_Sprimme(aH, basisSize, aBasisSize, basisSize, ctx);
         CHKERR(Num_gemm_Sprimme("N", "N", basisSize, aBasisSize, aBasisSize,
               1.0, &hVecs[ldhVecs*j], ldhVecs, ahVecs, ldhVecsRot, 0.0,
               aH, basisSize, ctx));
         Num_copy_matrix_Sprimme(aH, basisSize, aBasisSize, basisSize,
               &hVecs[ldhVecs*j], ldhVecs, ctx);
         CHKERR(Num_free_Sprimme(aH, ctx));

         /* Indicate that before i may not be singular vectors */
         *arbitraryVecs = i;
      }
   }

   return 0;
}

/*******************************************************************************
 * Function map_vecs: given two basis V and W the function returns the
 *    permutation p such as p[i] is the column in V closest in angle to
 *    column ith on W.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V,W          The orthonormal bases
 *
 * m,n          Dimensions on V and W            
 *
 * n0           Update p[n0:n-1]
 *
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * perm         The returned permutation
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int map_vecs_Sprimme(HSCALAR *V, int m, int nV, int ldV, HSCALAR *W, int n0,
      int n, int ldW, int *p, primme_context ctx) {

   int i;         /* Loop variable                                     */

   /* Compute the norm of the columns V(n0:n-1) */

   HREAL *Vnorms = NULL;
   CHKERR(Num_malloc_RHprimme(nV, &Vnorms, ctx));
   for (i = 0; i < nV; i++) {
      Vnorms[i] = sqrt(REAL_PART(
            Num_dot_SHprimme(m, &V[ldV * i], 1, &V[ldV * i], 1, ctx)));
   }
      
   /* Compute V'*W[n0:n-1] */

   HSCALAR *ip = NULL;
   CHKERR(Num_malloc_SHprimme(nV * (n - n0), &ip, ctx));
   Num_zero_matrix_SHprimme(ip, nV, n - n0, nV, ctx);
   CHKERR(Num_gemm_SHprimme("C", "N", nV, n - n0, m, 1.0, V, ldV, &W[ldW * n0],
         ldW, 0.0, ip, nV, ctx));

   for (i = n0; i < n; i++) {
      /* Find the j that maximizes ABS(V[j]'*W[i]/Vnorms[j]) and is not */
      /* in p(0:i-1)                                                    */

      int j, jmax=-1;
      HREAL ipmax = -1;
      for (j = 0; j < nV; j++) {
         HREAL ipij = ABS(ip[nV * (i - n0) + j]);
         if (ipij > ipmax * Vnorms[j]) {
            /* Check that j is not in p(0:i-1) */
            int k;
            for (k = 0; k < i && p[k] != j; k++)
               ;
            if (k < i) continue;

            /* Update ipmax and jmax */
            ipmax = fabs(ipij / Vnorms[j]);
            jmax = j;
         }
      }
      if (jmax < 0) {
         jmax = i;
      }

      /* Assign the map */

      p[i] = jmax;
   }
   
   CHKERR(Num_free_RHprimme(Vnorms, ctx));
   CHKERR(Num_free_SHprimme(ip, ctx));
    
   return 0;
}

#endif /* defined(USE_HOST) && ((!defined(USE_HALF) && !defined(USE_HALFCOMPLEX)) || defined(BLASLAPACK_WITH_HALF)) */

#endif /* SUPPORTED_TYPE */
