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
 * File: ortho.c
 *
 * Purpose - Orthonormalizes a block of vectors, vector by vector, 
 *           against two bases and among themselves. Gram-Schmidt is used 
 *           with reorthogonalization based on Daniel's test. 
 *           For the purpose of the test, the norm of the resulting vector 
 *           is computed without synchronizations. Because of floating point,
 *           this norm is not accurate if reortho is needed, but it is more 
 *           than sufficient for the test. If at least one reortho has been
 *           performed, the actual norm must be computed. 
 *
 * Note on numerical properties
 *           The relative error in the implicit s1 is 
 *                    O( machineEps * (s0/s1)^2 ), 
 *           so when Daniel's test succeeds error(s1) = O(2*machineEps)).
 *
 *           For almost linearly dependent vectors a large error in s1
 *           might fail to detect a problem during the second reortho.
 *           Let s1 = s1hat*(1+err), with err = O(macheps*(s0/s1)^2).
 *              eg. s0/s1 = 1e8 => (s1-s1hat)/s1hat = O(1)
 *              eg. s0/s1 = 1e7 => (s1-s1hat)/s1hat = O(1e-2)
 *              eg. s0/s1 = 1e4 => (s1-s1hat)/s1hat = O(1e-8)
 *           When s0 is reset to s00 = s1, and perform again the Daniel's test
 *              s11/s00 = s11/s1hat(1+err) < 0.75 <=> s11/s1hat < .75*(1+err)
 *           if err is around 1e-2 it does not affect Daniels test a lot.
 *           Thus, if s0/s1 > 1.49e7, we must compute the actual s1.
 *
 * Note on linear dependency
 *           If the original vector norm is reduced more than machEps times,
 *           it is considered linear dependent. If R is not returned, the
 *           vector is replaced by a random vector. Otherwise the vector is
 *           zeroed.
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/ortho.c"
#endif

#include <assert.h>
#include "numerical.h"
#include "common_eigs.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "factorize.h"
#include "auxiliary_eigs.h"
#include "ortho.h"
#endif
 
#ifdef SUPPORTED_TYPE

/**********************************************************************
 * Function ortho - This routine orthonormalizes
 * a block of of vectors (from b1 to including b2 in basis)
 * against other vectors in the same array (from 0 to b1-1 in basis),
 * against a set of locked vectors (from 0 to numLocked-1 in locked),
 * and themselves.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * ldV        Leading dimension of the basis
 * b1, b2     Range of indices of vectors to be orthonormalized
 *            (b1 can be zero, but b1 must be <= b2)
 * ldR        Leading dimension in R
 * locked     Array that holds locked vectors if they are in-core
 * ldLocked   Leading dimension of locked
 * numLocked  Number of vectors in locked
 * RLocked    R for locked vectors
 * ldRLocked  The leading dimension of RLocked
 * nLocal     Number of rows of each vector stored on this node
 *
 * rworkSize  Length of rwork array
 * primme     Primme struct. Contains globalSumReal and Parallelism info
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * V       Basis vectors
 * R       Rotations done in the basis: V_input = V_output * R
 * iseed   Seeds used to generate random vectors
 * rwork   Contains buffers and other necessary work arrays
 *
 * Return Value
 * ------------
 * >0  - Insufficient workspace provided.  Positive number indicates
 *       appropriate amount.
 *  0  - success
 * -1  - some size or leading dimension was < 0
 * -2  - b1 > b2
 * -3  - A limit number of randomizations has been performed without
 *       yielding an orthogonal direction
 * 
 **********************************************************************/

TEMPLATE_PLEASE
int Bortho_gen_Sprimme(SCALAR *V, PRIMME_INT ldV, HSCALAR *R, int ldR,
      int b1, int b2, SCALAR *locked, PRIMME_INT ldLocked, int numLocked,
      HSCALAR *RLocked, int ldRLocked, PRIMME_INT nLocal,
      int (*B)(SCALAR *, PRIMME_INT, SCALAR *, PRIMME_INT, int, void *),
      void *Bctx, PRIMME_INT *iseed, int *b2_out, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i;                   /* Loop indices */
   /* TODO: replace by a dynamic criterion when to stop orthogonalizing local */
   /* vectors. Observed performance improvement when maxNumOrthos increases.  */
   int maxNumOrthos = primme?3:7; /* We let 2 reorthogonalizations before randomize */
                                  /* for nLocal length vectors, and 6 orthogonalisations */
                                  /* for the rest */
   int maxNumRandoms = 10;  /* We do not allow more than 10 randomizations */
   double tol = sqrt(2.0L)/2.0L; /* We set Daniel et al. test to .707 */
   double t0;

   /*----------------------------------*/
   /* input and workspace verification */
   /*----------------------------------*/
   assert(nLocal >= 0 && numLocked >= 0 &&
          ldV >= nLocal && (numLocked == 0 || ldLocked >= nLocal) &&
          (R == NULL || ldR > b2) && (!RLocked || ldRLocked >= numLocked));

   tol = sqrt(2.0L)/2.0L;

   /* Zero the columns from b1 to b2 of R */
   if (R) {
      CHKERR(Num_zero_matrix_SHprimme(
            &R[b1 * ldR], b2 + 1, b2 - b1 + 1, ldR, ctx));
   }
   if (RLocked) {
      CHKERR(Num_zero_matrix_SHprimme(
            RLocked, numLocked, b2 - b1 + 1, ldRLocked, ctx));
   }

   /*---------------------------------------------------*/
   /* main loop to orthogonalize new vectors one by one */
   /*---------------------------------------------------*/

   t0 = primme_wTimer();

   // Allocate overlaps and Bx

   HSCALAR *overlaps;
   CHKERR(Num_malloc_SHprimme(b2+1 + numLocked, &overlaps, ctx));
   Num_zero_matrix_SHprimme(overlaps, 1, b2 + 1 + numLocked, 1, ctx);
   SCALAR *Bx = NULL;
   CHKERR(Num_malloc_Sprimme(B?nLocal:0, &Bx, ctx));
   if (b2_out) *b2_out = b1;
   double eps_orth;
   CHKERR(machineEpsOrth_Sprimme(&eps_orth, ctx));

   for(i=b1; i <= b2; i++) {
    
      int nOrth;  // number of orthogonalizations of the current vector
      int randomizations;  // times the current vector has been replaced by a
                           // random vector
      int reorth = 1;   // flag to keep iterating
      int updateR =
            (R || RLocked ? 1 : 0); // flag to keep updating R, after the first
      int Bx_update = 0;   // flag indicating if Bx = B*V[i]
                                 // randomization it is set to zero
      HREAL s0=0.0;   // B-norm of the current vector before deflating V and locked
      HREAL s02=0.0;  // s0 squared
      HREAL s1=0.0;   // B-norm of the current vector after deflating V and locked
      HREAL s12=0.0;  // s1 squared

      for (nOrth=0, randomizations=0; reorth; ) {

         if (nOrth >= maxNumOrthos) {
            /* Stop updating R when replacing one of the columns of the V */
            /* with a random vector                                           */

            if (updateR) {
               if (R) R[ldR*i + i] = 0.0;
               updateR = 0;
            }

            if (randomizations >= maxNumRandoms) {
               goto clean;
            }
            PRINTF(5, "Randomizing in ortho: %d, vector size of %" PRIMME_INT_P,
                  i, nLocal);

            CHKERR(Num_larnv_Sprimme(2, iseed, nLocal, &V[ldV*i], ctx)); 
            randomizations++;
            nOrth = 0;
            Bx_update = 0;    // V[i] has changed, so Bx != B*V[i]
         }

         nOrth++;

         // Compute B*V[i]

         if (B && !Bx_update) {
            CHKERR(B(&V[ldV*i], ldV, Bx, nLocal, 1, Bctx));
         }
         else if (!B) {
            Bx = &V[ldV*i];
         }

         // Compute the B norm of the current vector, V[i], if it wasn't computed
         // in previous iteration

         if (nOrth == 1) {
            s02 = REAL_PART(Num_dot_Sprimme(nLocal, &V[ldV*i], 1, Bx, 1, ctx));
            if (primme) primme->stats.numOrthoInnerProds += 1;
         }

         // Compute overlaps = [V[0:i-1]'*B*V[i] 
 
         if (i > 0) {
            CHKERR(Num_gemv_ddh_Sprimme(
                  "C", nLocal, i, 1.0, V, ldV, Bx, 1, 0.0, overlaps, 1, ctx));
            if (primme) primme->stats.numOrthoInnerProds += i;
         }

         if (numLocked > 0) {
            CHKERR(Num_gemv_ddh_Sprimme("C", nLocal, numLocked, 1.0, locked,
                  ldLocked, Bx, 1, 0.0, &overlaps[i], 1, ctx));
            if (primme) primme->stats.numOrthoInnerProds += numLocked;
         }

         overlaps[i+numLocked] = s02;
         CHKERR(globalSum_SHprimme(overlaps, i + numLocked + 1, ctx));

         if (updateR) {
            if (R) {
               CHKERR(Num_axpy_SHprimme(
                     i, 1.0, overlaps, 1, &R[ldR * i], 1, ctx));
            }
            if (RLocked) {
               CHKERR(Num_axpy_SHprimme(numLocked, 1.0, overlaps + i, 1,
                     &RLocked[ldRLocked * (i - b1)], 1, ctx));
             }
         }

         if (numLocked > 0) { /* locked array most recently accessed */
            // Compute V[i] = V[i] - locked'*overlaps[i:i+numLocked-1]
            CHKERR(Num_gemv_dhd_Sprimme("N", nLocal, numLocked, -1.0, locked,
                  ldLocked, &overlaps[i], 1, 1.0, &V[ldV * i], 1, ctx));
            if (primme) primme->stats.numOrthoInnerProds += numLocked;
         }

         if (i > 0) {
            CHKERR(Num_gemv_dhd_Sprimme("N", nLocal, i, -1.0, V, ldV, overlaps,
                  1, 1.0, &V[ldV * i], 1, ctx));
            if (primme) primme->stats.numOrthoInnerProds += i;
         }

         Bx_update = 0;    // V[i] has changed, so Bx != B*V[i]
 
         if (nOrth == 1) {
            s0 = sqrt(s02 = REAL_PART(overlaps[i+numLocked]));
         }

         /* Compute the norm s1 explicitly */
         
         if (B) {
            CHKERR(B(&V[ldV*i], ldV, Bx, nLocal, 1, Bctx));
            Bx_update = 1;
         }
         s12 = REAL_PART(
               Num_dot_Sprimme(nLocal, &V[ldV * i], 1, Bx, 1, ctx));
         if (primme) primme->stats.numOrthoInnerProds += 1;
         CHKERR(globalSum_RHprimme(&s12, 1, ctx));
         s1 = sqrt(s12);

         if (!ISFINITE(s0) || !ISFINITE(s1) || s1 <= eps_orth * s0) {
            PRINTF(5, "Vector %d lost all significant digits in ortho",
                  i - b1);
            nOrth = maxNumOrthos;
         }
         else if (s1 <= tol*s0 || (!primme && nOrth < maxNumOrthos)) {
            /* No numerical benefit in normalizing the vector before reortho */
            s0 = s1;
            s02 = s12;
         }
         else {
            if (updateR && R) R[ldR * i + i] = s1;

            if (ISFINITE((HREAL)(1.0/s1))) {
               CHKERR(Num_scal_Sprimme(nLocal, 1.0/s1, &V[ldV*i], 1, ctx));
               break;
            }
            else {
               PRINTF(5, "Vector %d lost all significant digits in ortho",
                     i - b1);
               nOrth = maxNumOrthos;
            }
         } 
      }
      if (b2_out) *b2_out = i+1;
   }

clean:
   if (primme) primme->stats.timeOrtho += primme_wTimer() - t0;

   CHKERR(Num_free_SHprimme(overlaps, ctx));
   if (B) CHKERR(Num_free_Sprimme(Bx, ctx));

   /* Check orthogonality */
   /*
   if (numLocked) {
      SCALAR *H = (SCALAR*)malloc(sizeof(SCALAR)*numLocked*numLocked);
      Num_gemm_Sprimme("C", "N", numLocked, numLocked, nLocal, 1.0, locked,
            ldLocked, locked, ldLocked, 0.0, H, numLocked);
      for(i=0; i < numLocked; i++) {
         for(j=0; j < i; j++) assert(ABS(H[numLocked*i+j]) < 1e-13);
         assert(fabs(1 - ABS(H[numLocked*i+i])) < 1e-13);
      }
      free(H);
   }
   if (b2+1) {
      SCALAR *H = (SCALAR*)malloc(sizeof(SCALAR)*(b2+1)*(b2+1));
      Num_gemm_Sprimme("C", "N", b2+1, b2+1, nLocal, 1.0, V,
            ldV, V, ldV, 0.0, H, b2+1);
      for(i=0; i < b2+1; i++) {
         for(j=0; j < i; j++) assert(ABS(H[(b2+1)*i+j]) < 1e-13);
         assert(H[(b2+1)*i+i] == 0.0 || fabs(1 - ABS(H[(b2+1)*i+i])) < 1e-13);
      }
      free(H);
   }
   if (numLocked) {
      SCALAR *H = (SCALAR*)malloc(sizeof(SCALAR)*(b2+1)*numLocked);
      Num_gemm_Sprimme("C", "N", numLocked, b2+1, nLocal, 1.0, locked,
            ldLocked, V, ldV, 0.0, H, numLocked);
      for(i=0; i < b2+1; i++) {
         for(j=0; j < numLocked; j++) assert(ABS(H[numLocked*i+j]) < 1e-13);
      }
      free(H);
   }
   */

   return 0;
}

TEMPLATE_PLEASE
int ortho_Sprimme(SCALAR *V, PRIMME_INT ldV, HSCALAR *R, int ldR, int b1, int b2,
                  SCALAR *locked, PRIMME_INT ldLocked, int numLocked,
                  PRIMME_INT nLocal, PRIMME_INT *iseed, primme_context ctx) {

   int b2_out;
   CHKERR(Bortho_gen_Sprimme(V, ldV, R, ldR, b1, b2, locked, ldLocked,
         numLocked, NULL, 0, nLocal, NULL, NULL, iseed, &b2_out, ctx));
   return b2_out == b2+1 ? 0 : -3;
}

#ifdef USE_HOST

#ifndef ORTHO__PRIVATE_H
#define ORTHO__PRIVATE_H
struct local_matvec_ctx { /* HSCALAR */
   void *B;
   int n, ldB;
   primme_context ctx;
};
#endif /* ORTHO__PRIVATE_H */

STATIC int local_matvec(HSCALAR *x, PRIMME_INT ldx, HSCALAR *y, PRIMME_INT ldy,
      int bs, void *Bctx_) {
   struct local_matvec_ctx *Bctx = (struct local_matvec_ctx*)Bctx_;
   primme_context ctx = Bctx->ctx;
   CHKERR(Num_zero_matrix_SHprimme(y, Bctx->n, 1, Bctx->n, ctx));
   CHKERR(Num_hemm_SHprimme("L", "U", Bctx->n, bs, 1.0, (HSCALAR *)Bctx->B,
         Bctx->ldB, x, ldx, 0.0, y, ldy, ctx));
   return 0;
}

TEMPLATE_PLEASE
int Bortho_local_Sprimme(HSCALAR *V, int ldV, HSCALAR *R,
      int ldR, int b1, int b2, HSCALAR *locked, int ldLocked,
      int numLocked, PRIMME_INT nLocal, HSCALAR *B, int ldB, PRIMME_INT *iseed,
      primme_context ctx) {

   /* Remove MPI communications from the context */

   ctx.numProcs = 1;
   ctx.procID = 0;
   ctx.mpicomm = NULL;
   ctx.primme = NULL;

   /* Call orthogonalization */

   struct local_matvec_ctx Bctx = {(void *)B, (int)nLocal, ldB, ctx};
   int b2_out;
   CHKERR(Bortho_gen_SHprimme(V, ldV, R, ldR, b1, b2, locked, ldLocked,
         numLocked, NULL, 0, nLocal, B ? local_matvec : NULL, &Bctx, iseed,
         &b2_out, ctx));
   return b2_out == b2 + 1 ? 0 : -3;
}

#endif /* USE_HOST */

STATIC int B_matvec(SCALAR *x, PRIMME_INT ldx, SCALAR *y, PRIMME_INT ldy,
      int bs, void *ctx_) {
   primme_context ctx = *(primme_context*)ctx_;
   CHKERR(massMatrixMatvec_Sprimme(
         x, ctx.primme->nLocal, ldx, y, ldy, 0, bs, ctx));
   return 0;
}


TEMPLATE_PLEASE
int Bortho_block_Sprimme(SCALAR *V, PRIMME_INT ldV, HSCALAR *VLtBVL,
      int ldVLtBVL, HSCALAR *fVLtBVL, int ldfVLtBVL, HSCALAR *R, PRIMME_INT ldR,
      int b1, int b2, SCALAR *locked, PRIMME_INT ldLocked, int numLocked,
      SCALAR *BV, PRIMME_INT ldBV, HSCALAR *RLocked, int ldRLocked,
      PRIMME_INT nLocal, int maxRank, int *b2_out, primme_context ctx) {

   return Bortho_block_gen_Sprimme(V, ldV, VLtBVL, ldVLtBVL, fVLtBVL, ldfVLtBVL,
         R, ldR, b1, b2, locked, ldLocked, numLocked,
         ctx.primme->massMatrixMatvec ? B_matvec : NULL, &ctx, BV, ldBV,
         RLocked, ldRLocked, nLocal, maxRank, b2_out, ctx);
}

TEMPLATE_PLEASE
int ortho_block_Sprimme(SCALAR *V, PRIMME_INT ldV, HSCALAR *VLtBVL,
      int ldVLtBVL, HSCALAR *fVLtBVL, int ldfVLtBVL, HSCALAR *R, PRIMME_INT ldR,
      int b1, int b2, SCALAR *locked, PRIMME_INT ldLocked, int numLocked,
      HSCALAR *RLocked, int ldRLocked, PRIMME_INT nLocal, int maxRank,
      int *b2_out, primme_context ctx) {

   return Bortho_block_gen_Sprimme(V, ldV, VLtBVL, ldVLtBVL, fVLtBVL, ldfVLtBVL,
         R, ldR, b1, b2, locked, ldLocked, numLocked, NULL, NULL, NULL, 0,
         RLocked, ldRLocked, nLocal, maxRank, b2_out, ctx);
}


/**********************************************************************
 * Function Bortho_gen_block - This routine orthonormalizes
 * a block of of vectors (from b1 to including b2 in basis)
 * against other vectors in the same array (from 0 to b1-1 in basis),
 * against a set of locked vectors (from 0 to numLocked-1 in locked),
 * and themselves.
 *
 * If given, the returning R and RLocked satisfy
 *
 *    input_V = [locked output_V] * [RLocked; R]
 *
 * INPUT/OUTPUT PARAMETERS
 * -----------------------
 * V          Basis vectors
 * ldV        Leading dimension of the basis
 * VLtBVL     [locked V]' * B * [locked V], referring V(0:b1-1) on input,
 *            and V(0:b2-1) on output
 * ldVLtBVL   The leading dimension of VLtBVL
 * fVLtBVL    The Cholesky factor of VLtBVL
 * ldfVLtBVL  The leading dimension of fVLtBVL
 * R          Rotations done in the basis regarding V
 * ldR        The leading dimension of R
 * b1, b2     Range of V columns to be orthonormalized
 * locked     Array that holds locked vectors if they are in-core
 * ldLocked   Leading dimension of locked
 * numLocked  Number of vectors in locked
 * B          Inner product
 * Bctx       Context for function B
 * BV         B*V
 * ldBV       The leading dimension of BV
 * RLocked    Rotations done in the basis regarding locked
 * ldRLocked  The leading dimension of RLocked
 * nLocal     Number of rows of each vector stored on this node
 * maxRank    largest rank of the basis being orthogonalized
 * b2_out     The number of linear independent columns
 * ctx        primme context
 *
 * Return Value
 * ------------
 *  error code
 * 
 **********************************************************************/

STATIC int Bortho_block_gen_Sprimme(SCALAR *V, PRIMME_INT ldV, HSCALAR *VLtBVL,
      int ldVLtBVL, HSCALAR *fVLtBVL, int ldfVLtBVL, HSCALAR *R, PRIMME_INT ldR,
      int b1, int b2, SCALAR *locked, PRIMME_INT ldLocked, int numLocked,
      int (*B)(SCALAR *, PRIMME_INT, SCALAR *, PRIMME_INT, int, void *),
      void *Bctx, SCALAR *BV, PRIMME_INT ldBV, HSCALAR *RLocked, int ldRLocked,
      PRIMME_INT nLocal, int maxRank, int *b2_out, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, j;               /* loop indices */
   HSCALAR *A, *C, *Y;      /* auxiliary local matrices */
   HSCALAR *VLtBVLdA;        /* auxiliary local matrices */
   HREAL *D, *N;            /* singular values */
   int ldA;                /* leading dimension of A */
   b2++; /* NOTE: Let's use C range convention */

   /* Quick exit */

   if (b2 <= b1) {
      *b2_out = b2;
      return 0;
   }

   double eps_orth;
   CHKERR(machineEpsOrth_Sprimme(&eps_orth, ctx));

   if (VLtBVL == NULL) {
      CHKERR(Bortho_gen_Sprimme(V, ldV, R, ldR, b1, b2 - 1, locked, ldLocked,
            numLocked, RLocked, ldRLocked, nLocal, B, Bctx, primme->iseed,
            b2_out, ctx));
      if (B && BV) {
         CHKERR(B(&V[ldV * b1], ldV, &BV[ldBV * b1], ldBV, *b2_out - b1, Bctx));
      }
      return 0;
   }

   // TEMP!!!
   // CHKERR(Bortho_gen_Sprimme(V, ldV, R, ldR, b1, b2 - 1, locked, ldLocked,
   //          numLocked, RLocked, ldRLocked, nLocal, NULL, NULL, primme->iseed,
   //          ctx));
   // Num_zero_matrix_SHprimme(&VLtBVL[ldVLtBVL*(numLocked+b1)], numLocked+b2, b2-b1, ldVLtBVL, ctx);
   // for (i=numLocked+b1; i<numLocked+b2; i++) {
   //    VLtBVL[ldVLtBVL*i+i] = 1.0;
   // }
   // *b2_out = b2;
   // return 0;
   
 
   double t0 = primme_wTimer();

   /* input and workspace verification */

   assert(nLocal >= 0 && numLocked >= 0 && ldV >= nLocal &&
          (numLocked == 0 || ldLocked >= nLocal) && (R == NULL || ldR >= b2) &&
          (!R || numLocked == 0) && (!RLocked || ldRLocked >= numLocked));

   /* Zero the columns from b1 to b2 of R */

   HSCALAR *r = NULL; /* The diagonal block of R(b1:b2) */
   int ldr = 0;

   if (R) {
      CHKERR(Num_zero_matrix_SHprimme(&R[ldR * b1], b1, b2 - b1, ldR, ctx));
      r = &R[ldR * b1 + b1];
      ldr = ldR;
   }

   if (RLocked) {
      CHKERR(Num_zero_matrix_SHprimme(
            RLocked, numLocked, b2 - b1, ldRLocked, ctx));
      if (!r) {
         CHKERR(Num_malloc_SHprimme((b2 - b1) * (b2 - b1), &r, ctx));
         ldr = b2 - b1;
      }
   }
   if (r) {
      CHKERR(Num_zero_matrix_SHprimme(r, b2 - b1, b2 - b1, ldr, ctx));
      for (i=0; i<b2 - b1; i++) {
         r[ldr*i+i] = 1.0;
      }
   }

   /* Allocate workspace */

   CHKERR(Num_malloc_RHprimme(b2-b1, &D, ctx));
   CHKERR(Num_malloc_RHprimme(b2-b1, &N, ctx));
   A = &VLtBVL[ldVLtBVL*(b1+numLocked)];
   ldA = ldVLtBVL;
   CHKERR(Num_malloc_SHprimme((numLocked+b1)*(b2-b1), &VLtBVLdA, ctx));
   CHKERR(Num_malloc_SHprimme((b2-b1)*(b2-b1), &Y, ctx));
   int nVL = b1+numLocked;
   CHKERR(Num_malloc_SHprimme((b2-b1)*(b2-b1), &C, ctx));
   SCALAR *BX;
   PRIMME_INT ldBX;
   if (B) {
      if (BV) {
         ldBX = ldBV;
         BX = &BV[b1 * ldBV];
      }
      else {
         ldBX = primme && nLocal == primme->nLocal ? primme->ldOPs : nLocal;
         CHKERR(Num_malloc_Sprimme(ldBX * (b2 - b1), &BX, ctx));
      }
   }
   else {
      ldBX = ldV;
      BX = &V[b1 * ldV];
   }

   /* Main loop to orthogonalize new vectors one by one. Just kidding. Although
    * seeming complicated it is just doing iterative CG-SVQB */

   *b2_out = b2; 
   int its;
   int maxits = 5, plus1 = 5;
   int Yortho = 1;
   for (its=0; its<maxits; its++) {
      /* Notation:                                          */
      /* Vp = [locked V(0:b1-1)]; Xp = V(b1:b2-1) */
      /* Xc = V(b1:b2-1)                                    */
      /* Xn is Xp in the next iteration                     */

      /* Do Xp <= (Xp - Vp*VLtBVLdA)/N*Y/D if it isn't the first iteration     */
      /* Do A <= [locked V(0:b2)]'*B*Xc always                                  */

      if (B) {
         if (its > 0) {
            CHKERR(Num_ortho_kernel(locked, nLocal, numLocked, ldLocked, V, b1,
                  b2, ldV, VLtBVLdA, nVL, D, Y, b2 - b1, Yortho, NULL, 0, NULL,
                  0, ctx));
         }

         /* Update BX */
         CHKERR(B(&V[ldV * b1], ldV, BX, ldBX, b2 - b1, Bctx));

         CHKERR(Num_ortho_kernel(locked, nLocal, numLocked, ldLocked, V, b1, b2,
               ldV, NULL, 0, NULL, NULL, 0, 0, BX, ldBX, A, ldA, ctx));
      } else {
         CHKERR(Num_ortho_kernel(locked, nLocal, numLocked, ldLocked, V, b1, b2,
               ldV, VLtBVLdA, nVL, its == 0 ? NULL : D, Y, b2 - b1, Yortho,
               &V[ldV * b1], ldV, A, ldA, ctx));
      }
      if (primme) {
         primme->stats.numOrthoInnerProds += (numLocked + b1) * (b2 - b1) +
#ifdef USE_HOST
                                             (b2 - b1) * (b2 - b1);
#else
                                             (b2 - b1 + 1) / 2 * (b2 - b1);
#endif
      }

      /* Check orthogonality. We stop the SVQB iteration one iteration after the
       * condition number of the basis is less than 3. This is check by
       * rank_estimation. The first time the basis passes the condition, plus1
       * is set to the index of the next iteration. */

      if (rank_estimation(VLtBVL, numLocked + b1,
                   numLocked + b2, maxRank, ldVLtBVL) == numLocked + b2) {
         if (its >= plus1) {
            /* Pass the check when the norm of V(b1:b2-1) are close to one */
            for (i = b1;
                  i < b2 &&
                  ABS(VLtBVL[ldVLtBVL * (numLocked + i) + numLocked + i] -
                        (HSCALAR)1.0) < .8;
                  i++)
               ;
            if (i >= b2) break;
         }
         else plus1 = min(its + 1, plus1);
      }

      if (ctx.procID == 0) {
         // If the norm of the vector overflows, discard the inner product with
         // other vectors and set the maximum value as the norm
         for (i = 0; i < b2 - b1; i++) {
            if (REAL_PART(A[ldA * i + i]) < MACHINE_MAX) continue;
            CHKERR(Num_zero_matrix_SHprimme(&VLtBVL[ldVLtBVL * (numLocked + i)],
                  numLocked + i, 1, ldVLtBVL, ctx));
            A[ldA * i + i] = (HSCALAR)MACHINE_MAX;
         }

         /* C = Xn'*B*Xn = Xc'*(I-B*Vc/(Vc'*B*Vc)*Vc'*B)*Xc = Xc'*B*Xc -
            (Xc'*B*Vc)*((Vc'*B*Vc)\(Vc'*B*Xc))       */
         /*   = A(b1:b2-1,:) -
            (A(0:numLocked+b1-1,:)'*(VLtBVL\A(0:numLocked+b1-1,:))   */

         CHKERR(Num_copy_matrix_SHprimme(
               &A[numLocked + b1], b2 - b1, b2 - b1, ldA, C, b2 - b1, ctx));
         CHKERR(Num_copy_matrix_SHprimme(
               A, nVL, b2 - b1, ldA, VLtBVLdA, nVL, ctx));
         CHKERR(Num_trsm_SHprimme("L", "U", "C", "N", nVL, b2 - b1, 1.0,
               fVLtBVL, ldfVLtBVL, VLtBVLdA, nVL, ctx));
         CHKERR(Num_gemm_SHprimme("C", "N", b2 - b1, b2 - b1, numLocked + b1,
               -1.0, VLtBVLdA, nVL, VLtBVLdA, nVL, 1.0, C, b2 - b1, ctx));
         CHKERR(Num_trsm_SHprimme("L", "U", "N", "N", nVL, b2 - b1, 1.0,
               fVLtBVL, ldfVLtBVL, VLtBVLdA, nVL, ctx));

         /* N(i) = ||Xc(:,i)||; normXc = max ||N(i)|| */

         for (i=0; i<b2-b1; i++) {
            N[i] = sqrt(max(ABS(C[(b2-b1)*i+i]), eps_orth));
         }
         for (i = 0; i < b2 - b1; i++)
            for (j = 0; j <= i; j++) C[(b2 - b1) * i + j] /= N[i] * N[j];

         /* [D, Y] = eig(C) or Y = chol(C) */

         CHKERR(decomposition(C, b2 - b1, b2 - b1, Y, b2 - b1, D, &Yortho, ctx));

         /* D = sqrt(D) */

         for (i=0; i<b2-b1; i++) {
            D[i] = sqrt(max(D[i], eps_orth * (b2 - b1)));
         }
      }

      CHKERR(broadcast_SHprimme(VLtBVLdA, nVL*(b2-b1), ctx));
      CHKERR(broadcast_RHprimme(N, b2-b1, ctx));
      CHKERR(broadcast_RHprimme(D, b2-b1, ctx));
      CHKERR(broadcast_SHprimme(Y, (b2 - b1) * (b2 - b1), ctx));
      CHKERR(broadcast_iprimme(&Yortho, 1, ctx));

      if (RLocked) {
         /* R(0:b1-1,:) += VtBV\Vc'*Xc*gi = VtBV\A(0:b1-1,:)*R(b1:b2-1,:) */

         CHKERR(Num_gemm_SHprimme("N", "N", numLocked, b2 - b1, b2 - b1, 1.0,
               VLtBVLdA, nVL, r, ldr, 1.0, RLocked, ldRLocked, ctx));
      }
      if (R) {
         /* R(0:b1-1,:) += VtBV\Vc'*Xc*gi = VtBV\A(0:b1-1,:)*R(b1:b2-1,:) */

         CHKERR(Num_gemm_SHprimme("N", "N", b1, b2 - b1, b2 - b1, 1.0,
               &VLtBVLdA[numLocked], nVL, r, ldr, 1.0, &R[ldR * b1], ldR, ctx));
      }
      if (r) {
         /* R(b1:b2-1,b1:b2-1) = N * R(b1:b2-1,b1:b2-1) */

         for (i = 0; i < b2 - b1; i++) {
            for (j = 0; j < b2 - b1; j++) { r[ldr * i + j] *= N[j]; }
         }

         /* C = Y' * R(b1:b2-1,b1:b2-1) */

         if (Yortho) {
            CHKERR(Num_gemm_SHprimme("C", "N", b2 - b1, b2 - b1, b2 - b1, 1.0,
                  Y, b2 - b1, r, ldr, 0.0, C, b2 - b1, ctx));
         } else {
            CHKERR(Num_copy_matrix_SHprimme(
                  r, b2 - b1, b2 - b1, ldr, C, b2 - b1, ctx));
            CHKERR(Num_trmm_SHprimme("L", "U", "N", "N", b2 - b1, b2 - b1, 1.0,
                  Y, b2 - b1, C, b2 - b1, ctx));
         }

         /* R(b1:b2-1,b1:b2-1) = D * C */

         for (i = 0; i < b2 - b1; i++) {
            for (j = 0; j < b2 - b1; j++) {
               r[ldr * i + j] = D[j] * C[(b2 - b1) * i + j];
            }
         }
      }

      /* Y = N \ Y or Y = Y * N */

      if (Yortho) {
         for (i = 0; i < b2 - b1; i++) {
            for (j = 0; j < b2 - b1; j++) {
               Y[(b2 - b1) * i + j] /= N[j];
            }
         }
      } else {
         for (i = 0; i < b2 - b1; i++) {
            for (j = 0; j < b2 - b1; j++) {
               Y[(b2 - b1) * i + j] *= N[i];
            }
         }
      }
   } /* end while, I hope you enjoyed the loop */

   /* Return the number of linearly independent columns */

   if (ctx.procID == 0) {
      b2 = rank_estimation(
                 VLtBVL, numLocked + b1, numLocked + b2, maxRank, ldVLtBVL) -
           numLocked;
   }
   CHKERR(broadcast_iprimme(&b2, 1, ctx));
   *b2_out = b2;

   /* Free workspaces */
  
   if (!R && RLocked) CHKERR(Num_free_SHprimme(r, ctx));
   CHKERR(Num_free_RHprimme(D, ctx));
   CHKERR(Num_free_RHprimme(N, ctx));
   CHKERR(Num_free_SHprimme(VLtBVLdA, ctx));
   CHKERR(Num_free_SHprimme(Y, ctx));
   CHKERR(Num_free_SHprimme(C, ctx));
   if (B && !BV) CHKERR(Num_free_Sprimme(BX, ctx));

   /* Update fVLtBVL */

   CHKERR(update_cholesky_Sprimme(VLtBVL, ldVLtBVL, fVLtBVL, ldfVLtBVL,
         numLocked + b1, numLocked + b2, ctx));

   if (primme) primme->stats.timeOrtho += primme_wTimer() - t0;

   return 0;
}

/**********************************************************************
 * Function ortho_single_iteration -- This function orthogonalizes
 *    applies ones the projector (I-BQQ') on X. Also returns
 *    the norms of ||(I-BQQ')X(i)||.
 *   
 * ARRAYS AND PARAMETERS
 * ----------------
 * Q               The basis of the projector I-BQQ'.
 * mQ, nQ, ldQ     Rows, columns and leading dimension of Q.
 * X               The vectors to apply the projector I-QQ'.
 * inX             Column indices to apply the projector (optional).
 * nX              Number of columns to apply the projector.
 * ldX             The leading dimension of X.
 * overlaps        The norms of Q'X(i) (optional).
 * norms           The norms of (I-QQ')X(i) (optional).
 * rwork           Auxiliary space
 * lrwork          Available rwork
 *
 ****************************************************************************/

TEMPLATE_PLEASE
int ortho_single_iteration_Sprimme(SCALAR *Q, int nQ, PRIMME_INT ldQ,
      SCALAR *BQ, PRIMME_INT ldBQ, HSCALAR *QtBQ, int ldQtBQ, SCALAR *X,
      int *inX, int nX, PRIMME_INT ldX, HREAL *norms, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, j, M=PRIMME_BLOCK_SIZE;
   PRIMME_INT mQ = primme->nLocal;

   double t0 = primme_wTimer();

   /* Check if the indices of inX are contiguous */

   if (inX && nX > 0) {
      for (i=0, j=inX[0]; i<nX && j==inX[i]; i++, j++);
      if (i >= nX) {
         X = &X[inX[0]*ldX];
         inX = NULL;
      }
   }

   /* y = Q'*X */

   HSCALAR *y;
   CHKERR(Num_malloc_SHprimme(nQ * nX, &y, ctx));
   CHKERR(Num_zero_matrix_SHprimme(y, nQ, nX, nQ, ctx));
   if (!inX) {
      CHKERR(Num_gemm_ddh_Sprimme(
            "C", "N", nQ, nX, mQ, 1.0, Q, ldQ, X, ldX, 0.0, y, nQ, ctx));
   }
   else {
      int m=min(M, mQ);
      SCALAR *X0;
      CHKERR(Num_malloc_Sprimme(m*nX, &X0, ctx));
      CHKERR(Num_zero_matrix_SHprimme(y, nQ, nX, nQ, ctx));
      for (i=0, m=min(M,mQ); i < mQ; i+=m, m=min(m,mQ-i)) {
         CHKERR(Num_copy_matrix_columns_Sprimme(
               &X[i], m, inX, nX, ldX, X0, NULL, m, ctx));
         CHKERR(Num_gemm_ddh_Sprimme(
               "C", "N", nQ, nX, m, 1.0, &Q[i], ldQ, X0, m, 1.0, y, nQ, ctx));
      }
      CHKERR(Num_free_Sprimme(X0, ctx));
   }
   primme->stats.numOrthoInnerProds += nQ*nX;

   /* Reduction on y */

   CHKERR(globalSum_SHprimme(y, nQ*nX, ctx));
   
   /* z = QtBQ\y */

   HSCALAR *z = NULL;
   if (QtBQ) {
      HSCALAR *fQtBQ;
      int *pQtBQ;
      CHKERR(Num_malloc_SHprimme(nQ * nX, &z, ctx));
      CHKERR(Num_malloc_SHprimme(nQ * nQ, &fQtBQ, ctx));
      CHKERR(Num_malloc_iprimme(nQ, &pQtBQ, ctx));
      CHKERR(UDUDecompose_SHprimme(QtBQ, ldQtBQ, fQtBQ, nQ, pQtBQ, nQ, ctx));
      CHKERR(UDUSolve_SHprimme(fQtBQ, pQtBQ, nQ, y, nX, nQ, z, nQ, ctx));
      CHKERR(Num_free_SHprimme(fQtBQ, ctx));
      CHKERR(Num_free_iprimme(pQtBQ, ctx));
   } else {
      z = y;
   }

   if (QtBQ) CHKERR(Num_free_SHprimme(z, ctx));

   /* X = X - BQ*(QtBQ\y); norms(i) = norm(X(i)) */

   SCALAR *X0 = NULL;
   int m=min(M, mQ);
   if (inX) {
      CHKERR(Num_malloc_Sprimme(m*nX, &X0, ctx));
   }
   if (norms) for (i=0; i<nX; i++) norms[i] = 0.0;
   for (i=0; i < mQ; i+=m, m=min(m,mQ-i)) {
      if (inX) {
         CHKERR(Num_copy_matrix_columns_Sprimme(
               &X[i], m, inX, nX, ldX, X0, NULL, m, ctx));
      }
      CHKERR(Num_gemm_dhd_Sprimme("N", "N", m, nX, nQ, -1.0, &BQ[i], ldBQ, y, nQ,
            1.0, inX ? X0 : &X[i], inX ? m : ldX, ctx));
      if (inX) {
         CHKERR(Num_copy_matrix_columns_Sprimme(
               X0, m, NULL, nX, m, &X[i], inX, ldX, ctx));
      }
      if (norms) {
         for (j = 0; j < nX; j++) {
            SCALAR *x = inX ? &X0[j * m] : &X[j * ldX + i];
            norms[j] += REAL_PART(Num_dot_Sprimme(m, x, 1, x, 1, ctx));
         }
      }
   }

   if (norms) {
      /* Store the reduction of norms */
      CHKERR(globalSum_RHprimme(norms, nX, ctx));
 
      for (i=0; i<nX; i++) norms[i] = sqrt(norms[i]);
      primme->stats.numOrthoInnerProds += nX;
   }

   CHKERR(Num_free_SHprimme(y, ctx));
   CHKERR(Num_free_Sprimme(X0, ctx));

   primme->stats.timeOrtho += primme_wTimer() - t0;

   return 0;
}

/******************************************************************************
 * Function Num_ortho_kernel - This subroutine performs the next operations:
 *
 *    V(b1:b2) = (V(b1:b2) - [Q V]*A)*Y/D if A, D and Y provided
 *    B = [Q V]'*W
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * Q           input matrix
 * m           number of rows on Q, V
 * nQ          number of columns of Q
 * ldQ         leading dimension of Q
 * V           input/output matrix of size m x b2
 * ldV         leading dimension of V
 * A           input matrix of size nQ+b2 x (b2-b1)
 * ldA         leading dimension of A
 * D           diagonal of a diagonal matrix (optional)
 * Y           input matrix of size (b2-b1) x (b2-b1) (optional)
 * ldY         leading dimension of Y (optional)
 * Yortho      Whether Y is orthonormal matrix (1) or a Cholesky factor (0)
 * W           input matrix of size m x b2-b1
 * ldW         leading dimension of W
 * B           output matrix of size (nQ+b2) x (b2-b1)
 * ldB         leading dimension of B
 *
 ******************************************************************************/

STATIC int Num_ortho_kernel(SCALAR *Q, PRIMME_INT M, int nQ, PRIMME_INT ldQ,
      SCALAR *V, int b1, int b2, PRIMME_INT ldV, HSCALAR *A, int ldA, HREAL *D,
      HSCALAR *Y, int ldY, int Yortho, SCALAR *W, PRIMME_INT ldW, HSCALAR *B,
      int ldB, primme_context ctx) {

   PRIMME_INT i;     /* Loop variables */
   PRIMME_INT m;
   /* Set X = V(b1:b2) */
   int nX = b2 - b1;
   SCALAR *X = &V[ldV * b1];
   SCALAR *Xo = NULL;

   if (D && Y) {
      if (nX == 1) {
         m = min(PRIMME_BLOCK_SIZE, M);   /* Number of rows in the cache */
         if (!Yortho) {
            Y[0] = (HSCALAR)1.0/Y[0];
            Yortho = 1;
         }
      } else {
         m = min(PRIMME_BLOCK_SIZE, M);   /* Number of rows in the cache */
      }
      if (Yortho) {
         CHKERR(Num_malloc_Sprimme(m * nX, &Xo, ctx));
         CHKERR(Num_zero_matrix_Sprimme(Xo, m, nX, m, ctx));
      }
   }
   else {
      m = M;
   }

   HSCALAR *Bo;
   int ldBo;
   if (ctx.numProcs == 1) {
      Bo = B;
      ldBo = ldB;
   } else {
      CHKERR(Num_malloc_SHprimme((nQ + b2) * nX, &Bo, ctx));
      ldBo = nQ + b2;
   }

   /* Zero Bo */
   if (Bo) CHKERR(Num_zero_matrix_SHprimme(Bo, nQ + b2, nX, ldBo, ctx));

   /* Y(:,i) = Y(:,i)/D[i] */
   if (D && Y) {
      if (Yortho) {
         for (i=0; i < nX; i++) D[i] = 1.0/D[i];
         CHKERR(Num_scale_matrix_SHprimme(Y, nX, nX, ldY, D, Y, ldY, ctx));
      }
   }

   for (i=0; i < M; i+=m, m=min(m,M-i)) {
      if (D && Y) {
         /* X = X - Q*A(0:nQ,:) */
         if (nQ > 0) {
            CHKERR(Num_gemm_dhd_Sprimme("N", "N", m, nX, nQ, -1.0, &Q[i], ldQ,
                  A, ldA, 1.0, &X[i], ldV, ctx));
         }

         /* X = X - V*A(nQ:nQ+b1) */
         if (b1 > 0) {
            CHKERR(Num_gemm_dhd_Sprimme("N", "N", m, nX, b1, -1.0, &V[i], ldV,
                  A + nQ, ldA, 1.0, &X[i], ldV, ctx));
         }

         /* Xo = X*Y */
         if (Yortho) {
            CHKERR(Num_gemm_dhd_Sprimme("N", "N", m, nX, nX, 1.0, &X[i], ldV, Y,
                  ldY, 0.0, Xo, m, ctx));
            /* X = Xo*D */
            CHKERR(Num_copy_matrix_Sprimme(Xo, m, nX, m, &X[i], ldV, ctx));
         } else {
            CHKERR(Num_trsm_hd_Sprimme(
                  "R", "U", "N", "N", m, nX, 1.0, Y, ldY, &X[i], ldV, ctx));
         }
      }

      if (!W || !Bo) continue;

      /* Bo(0:nQ-1,:) += Q'*W */
      CHKERR(Num_gemm_ddh_Sprimme("C", "N", nQ, nX, m, 1.0, &Q[i], ldQ, &W[i],
            ldW, i == 0 ? 0.0 : 1.0, Bo, ldBo, ctx));
#ifdef USE_HOST
      /* Bo(nQ:nQ+b2-1,:) += V(:b2-1)'*W */
      CHKERR(Num_gemm_ddh_Sprimme("C", "N", b2, nX, m, 1.0, &V[i], ldV, &W[i],
            ldW, i == 0 ? 0.0 : 1.0, Bo + nQ, ldBo, ctx));
#else
      /* Bo(nQ:nQ+b1-1,:) += V(:b1-1)'*W */
      CHKERR(Num_gemm_ddh_Sprimme("C", "N", b1, nX, m, 1.0, &V[i], ldV, &W[i],
            ldW, i == 0 ? 0.0 : 1.0, Bo + nQ, ldBo, ctx));

      /* Bo(nQ+b1:nQ+b2-1,:) += V(b1:b2-1)'*W */
      CHKERR(Num_compute_gramm_ddh_Sprimme(&V[ldV * b1 + i], m, nX, ldV, &W[i],
            ldW, i == 0 ? 0.0 : 1.0, Bo + nQ + b1, ldBo, 1 /* symmetric */,
            ctx));
#endif
   }

   /* B = globalSum(Bo) */
   if (ctx.numProcs > 1) {
      CHKERR(globalSum_SHprimme(Bo, (nQ + b2) * nX, ctx));
      CHKERR(Num_copy_matrix_SHprimme(Bo, nQ+b2, nX, nQ+b2, B, ldB, ctx));
   }

   if (ctx.numProcs > 1) CHKERR(Num_free_SHprimme(Bo, ctx));
   if (D && Y && Yortho) CHKERR(Num_free_Sprimme(Xo, ctx));

   return 0;
}

/*******************************************************************************
 * Subroutine decomposition - This procedure compute the Cholesky factor of H.
 *    If it fails, it returns instead the eigendecomposition of H.
 *        
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * H             The input matrix of size n x n
 * n             Size of H
 * ldH           The leading dimension of H
 * primme        Structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * Y             The Cholesky factor of H or the right singular vectors of H
 * ldY           The leading dimension of Y
 * Yortho        Return whether Y is orthonormal matrix (1) or a Cholesky factor (0)
 * evals         All ones if Yortho==1; otherwise the eigenvalues of H, non-increasing order
 *
 * Return Value
 * ------------
 * error code
 ******************************************************************************/

STATIC int decomposition(HSCALAR *H, int n, int ldH, HSCALAR *Y, int ldY,
      HREAL *evals, int *Yortho, primme_context ctx) {

   int i, j; /* Loop variables    */

   /* Try Cholesky factorization */

   CHKERR(Num_copy_matrix_SHprimme(H, n, n, ldH, Y, ldY, ctx));
   int info;
   CHKERR(Num_potrf_SHprimme("U", n, Y, ldY, &info, ctx));
   if (info == 0) {
      *Yortho = 0;

      /* evals = 1 */

      for (i=0; i < n; i++) { 
         evals[i] = 1.0;
      }
   } else {
      /* Otherwise use eigenvalue decomposition */

      /* XHEEV returns the eigenvalues in ascending order. This function will */
      /* return them in descending order. An easy way is to compute the       */
      /* eigenpairs of -H instead. Do Y = -H.                                 */

      for (i=0; i < n; i++) {
         for (j = 0; j <= i; j++) { Y[ldY * i + j] = -H[ldH * i + j]; }
      }
 
      CHKERR(Num_heev_SHprimme("V", "U", n, Y, ldY, evals, ctx));

      /* evals = -evals */

      for (i=0; i < n; i++) { 
         evals[i] *= -1.0;
      }
      *Yortho = 1;
   }

   return 0;
}

/*******************************************************************************
 * Subroutine rank_estimation - Return the number of consecutive linearly
 *    independent columns in Q, given the inner-product matrix, V=Q'Q.
 *
 *    The routine checks that ||D*V*D - I||_F^2 < 0.8, where D is the diagonal
 *    matrix that makes the diagonal of D*V*D all ones. If that is the case,
 *    Q is not rank deficient.
 *
 *    The routine starts inspecting from column n0 of V, which corresponds to
 *    the inner products for the column n0 on the basis Q, and checks that the
 *    inner-product with itself is nonzero (that is, the column Q is not null),
 *    and that the cosines of the angles to the previous columns are less than
 *    .8/n. Under that condition, the condition number of Q should be 3 at most. 
 *
 * INPUT PARAMETERS
 * ----------------
 * V             The inner-product matrix, of size n1 x n1
 * ldV           The leading dimension of V
 * n0, n1        The range of columns to check
 * n             The maximum size that V will have
 *
 * Return Value
 * ------------
 * int       The number of consecutive columns that are linearly independent
 ******************************************************************************/

STATIC int rank_estimation(HSCALAR *V, int n0, int n1, int n, int ldV) {

   int i, j;

   for(i=n0; i<n1; i++) {
      HREAL Vii = REAL_PART(V[i * ldV + i]);
      if (!ISFINITE(Vii) || Vii <= (HREAL)0.0) break;
      for (j = 0; j < i; j++) {
         if (ABS(V[i * ldV + j]) >
               .8 / n * sqrt(Vii * REAL_PART(V[j * ldV + j])))
            break;
      }
      if (j < i) break;
   }

   return i;
}

/*******************************************************************************
 * Subroutine update_cholesky - update the Cholesky factor given new columns
 *    and rows on the original matrix.
 *
 * INPUT PARAMETERS
 * ----------------
 * V             The inner-product matrix, of size n1 x n1
 * ldV           The leading dimension of V
 * n0, n1        The range of columns to check
 * n             The maximum size that V will have
 *
 * Return Value
 * ------------
 * int       The number of consecutive columns that are linearly independent
 ******************************************************************************/

TEMPLATE_PLEASE
int update_cholesky_Sprimme(HSCALAR *VtV, int ldVtV, HSCALAR *fVtV, int ldfVtV,
      int n0, int n, primme_context ctx) {

   HSCALAR *A;
   CHKERR(Num_malloc_SHprimme(n * (n - n0), &A, ctx));
   if (ctx.procID == 0) {
      CHKERR(Num_copy_matrix_SHprimme(&VtV[ldVtV * n0], n,
            n - n0, ldVtV, A, n, ctx));
      CHKERR(Num_trsm_SHprimme("L", "U", "C", "N", n0, n - n0, 1.0, fVtV,
            ldfVtV, A, n, ctx));
      CHKERR(Num_gemm_SHprimme("C", "N", n - n0, n - n0, n0, -1.0, A, n, A, n,
            1.0, &A[n0], n, ctx));
      CHKERR(Num_potrf_SHprimme("U", n - n0, &A[n0], n, NULL, ctx));
   }
   CHKERR(broadcast_SHprimme(A, n * (n - n0), ctx));
   CHKERR(Num_copy_matrix_SHprimme(
         A, n, n - n0, n, &fVtV[ldfVtV * n0], ldfVtV, ctx));
   CHKERR(Num_free_SHprimme(A, ctx));

   return 0;
}

#endif /* SUPPORTED_TYPE */
