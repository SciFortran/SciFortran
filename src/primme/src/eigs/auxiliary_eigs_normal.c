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
 * File: auxiliary_eigs_normal.c
 *
 * Purpose - Miscellanea functions used by PRIMME EIGS dependent on the
 *           Hermiticity/normality of the operator
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/auxiliary_eigs_normal.c"
#endif

#include <string.h> /* memset */
#include "common_eigs.h"
#include "numerical.h"
#include "template_normal.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "auxiliary_eigs_normal.h"
#include "auxiliary_eigs.h"
#endif

#ifdef SUPPORTED_TYPE

/******************************************************************************
 * Function Num_compute_residuals - This subroutine performs the next operation
 *    in a cache-friendly way:
 *
 *    r = Ax - eval*Bx
 *
 * PARAMETERS
 * ---------------------------
 * n           The number of rows of x, Ax and r
 * eval        The value to compute the residual vector r
 * Bx          The vector Bx
 * Ax          The vector Ax
 * r           On output r = Ax - eval*Bx
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_compute_residuals_Sprimme(PRIMME_INT m, int n, HEVAL *eval,
      SCALAR *Bx, PRIMME_INT ldBx, SCALAR *Ax, PRIMME_INT ldAx, SCALAR *r,
      PRIMME_INT ldr, primme_context ctx) {

#ifdef USE_HOST
   int j;
   for (j = 0; j < n; j++) {
      int k, M = min(m, PRIMME_BLOCK_SIZE);
      for (k = 0; k < m; k += M, M = min(M, m - k)) {
         CHKERR(Num_copy_Sprimme(
               M, &Ax[ldAx * j + k], 1, &r[ldr * j + k], 1, ctx));
         CHKERR(Num_axpy_Sprimme(
               M, -eval[j], &Bx[ldBx * j + k], 1, &r[ldr * j + k], 1, ctx));
      }
   }

#else
   // Cache is not exploit for GPU; also Num_copy_Sprimme has a lot of overhead

   int j;

   CHKERR(Num_copy_matrix_Sprimme(Ax, m, n, ldAx, r, ldr, ctx));
   for (j = 0; j < n; j++) {
      CHKERR(Num_axpy_Sprimme(
            m, -eval[j], &Bx[ldBx * j], 1, &r[ldr * j], 1, ctx));
   }
#endif

   return 0;
}

/******************************************************************************
 * Function Num_update_VWXR - This subroutine performs the next operations:
 *
 *    X0 = V*h(nX0b+1:nX0e), X1 = V*h(nX1b+1:nX1e), X2 = V*h(nX2b+1:nX2e)
 *    Wo = W*h(nWob+1:nWoe),
 *    BX0 = BV*h(nBX0b+1:nBX0e), BX1 = BV*h(nBX1b+1:nBX1e)
 *    R = W*h(nRb+1:nRe) - BV*h(nRb+1:nRe)*diag(hVals(nRb+1:nRe)),
 *    Rnorms = norms(R),
 *    rnorms = norms(W*h(nrb+1:nre) - BV*h(nrb+1:nre)*diag(hVals(nrb+1:nre)))
 *    xnorms = norms(V*h(nxb+1:nxe))
 *    G = (V*h(1:nG))'*BV*h(1:nG)
 *    H = (V*h(1:nH))'*W*h(1:nH)
 *
 * NOTE: if Rnorms and rnorms are requested, nRb-nRe+nrb-nre < mV
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V, W, BV    input basis
 * mV,nV,ldV   number of rows and columns and leading dimension of V, W and BV
 * h           input rotation matrix
 * nh          Number of columns of h
 * ldh         The leading dimension of h
 * hVals       Array of values
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * X0             Output matrix V*h(nX0b:nX0e-1) (optional)
 * nX0b, nX0e     Range of columns of h
 * X1             Output matrix V*h(nX1b:nX1e-1) (optional)
 * nX1b, nX1e     Range of columns of h
 * X2             Output matrix V*h(nX2b:nX2e-1) (optional)
 * nX2b, nX2e     Range of columns of h
 * Wo             Output matrix W*h(nWob:nWoe-1) (optional)
 * nWob, nWoe     Range of columns of h
 * BX0            Output matrix BV*h(nBX0b:nBX0e-1) (optional)
 * nBX0b, nBX0e   Range of columns of h
 * BX1            Output matrix BV*h(nBX1b:nBX1e-1) (optional)
 * nBX1b, nBX1e   Range of columns of h
 * BX2            Output matrix BV*h(nBX2b:nBX2e-1) (optional)
 * nBX2b, nBX2e   Range of columns of h
 * R              Output matrix (optional)
 * nRb, nRe       Range of columns of h and hVals
 * Rnorms         Output array with the norms of R (optional)
 * rnorms         Output array with the extra residual vector norms (optional)
 * nrb, nre       Columns of residual vector to compute the norm
 * xnorms         Output array with V*h(nxb:nxe-1) vector norms (optional)
 * nxb, nxe       Columns of V*h to compute the norm
 * 
 * NOTE: n*e, n*b are zero-base indices of ranges where the first value is
 *       included and the last isn't.
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_update_VWXR_Sprimme(SCALAR *V, SCALAR *W, SCALAR *BV, PRIMME_INT mV,
      int nV, PRIMME_INT ldV, HSCALAR *h, int nh, int ldh, HEVAL *hVals,
      SCALAR *X0, int nX0b, int nX0e, PRIMME_INT ldX0,
      SCALAR *X1, int nX1b, int nX1e, PRIMME_INT ldX1,
      SCALAR *X2, int nX2b, int nX2e, PRIMME_INT ldX2,
      SCALAR *Wo, int nWob, int nWoe, PRIMME_INT ldWo,
      SCALAR *R, int nRb, int nRe, PRIMME_INT ldR, HREAL *Rnorms,
      SCALAR *BX0, int nBX0b, int nBX0e, PRIMME_INT ldBX0,
      SCALAR *BX1, int nBX1b, int nBX1e, PRIMME_INT ldBX1,
      SCALAR *BX2, int nBX2b, int nBX2e, PRIMME_INT ldBX2,
      HREAL *rnorms, int nrb, int nre,
      HSCALAR *G, int nG, int ldG,
      HSCALAR *H, int nH, int ldH,
      HREAL *xnorms, int nxb, int nxe,
      primme_context ctx) {

   (void)nh;
 
   PRIMME_INT i;     /* Loop variables */
   int j;            /* Loop variables */
   int m=min(PRIMME_BLOCK_SIZE, mV);   /* Number of rows in the cache */
   int nXb, nXe, nYb, nYe, nBXb, nBXe, ldG0=0, ldH0=0;
   PRIMME_INT ldX, ldY, ldBX;
   SCALAR *X, *Y, *BX;
   HSCALAR *G0=NULL, *H0=NULL, *workGH = NULL;

   assert(mV <= ldV && nh <= ldh && (!G || nG <= ldG) && (!H || nH <= ldH));
   assert(Rnorms == NULL || nRb >= nRe || nrb >= nre || Rnorms != rnorms);

   /* Figure out which columns of V*h, W*h and BV*h to compute */

   nXb = nYb = nBXb = INT_MAX;
   nXe = nYe = nBXe = 0;

   if (X0 && nX0b < nX0e) nXb = min(nXb, nX0b), nXe = max(nXe, nX0e);
   if (X1 && nX1b < nX1e) nXb = min(nXb, nX1b), nXe = max(nXe, nX1e);
   if (X2 && nX2b < nX2e) nXb = min(nXb, nX2b), nXe = max(nXe, nX2e);
   if (R && nRb < nRe && !BV) nXb = min(nXb, nRb), nXe = max(nXe, nRe);
   if (G && nG > 0) nXb = min(nXb, 0), nXe = max(nXe, nG);
   if (H && nH > 0) nXb = min(nXb, 0), nXe = max(nXe, nH);
   if (rnorms && nrb < nre && !BV) nXb = min(nXb, nrb), nXe = max(nXe, nre);
   if (xnorms && nxb < nxe) nXb = min(nXb, nxb), nXe = max(nXe, nxe);
   nXe = max(nXe, nXb);

   if (Wo && nWob < nWoe) nYb = min(nYb, nWob), nYe = max(nYe, nWoe);
   if (R && nRb < nRe) nYb = min(nYb, nRb), nYe = max(nYe, nRe);
   if (rnorms && nrb < nre) nYb = min(nYb, nrb), nYe = max(nYe, nre);
   if (H && nH > 0) nYb = min(nYb, 0), nYe = max(nYe, nH);
   nYe = max(nYe, nYb);

   if (BV) {
      if (BX0 && nBX0b < nBX0e) nBXb = min(nBXb, nBX0b), nBXe = max(nBXe, nBX0e);
      if (BX1 && nBX1b < nBX1e) nBXb = min(nBXb, nBX1b), nBXe = max(nBXe, nBX1e);
      if (BX2 && nBX2b < nBX2e) nBXb = min(nBXb, nBX2b), nBXe = max(nBXe, nBX2e);
      if (R && nRb < nRe) nBXb = min(nBXb, nRb), nBXe = max(nBXe, nRe);
      if (G && nG > 0) nBXb = min(nBXb, 0), nBXe = max(nBXe, nG);
      if (rnorms && nrb < nre) nBXb = min(nBXb, nrb), nBXe = max(nBXe, nre);
   }
   nBXe = max(nBXe, nBXb);

   assert(nXe <= nh || nXb >= nXe); /* Check dimension */
   assert(nYe <= nh || nYb >= nYe); /* Check dimension */
   assert(nBXe <= nh || nBXb >= nXe); /* Check dimension */

   double t0 = primme_wTimer();

   CHKERR(Num_malloc_Sprimme(m * (nXe - nXb), &X, ctx));
   CHKERR(Num_malloc_Sprimme(m * (nYe - nYb), &Y, ctx));
   CHKERR(Num_malloc_Sprimme(m * (nBXe - nBXb), &BX, ctx));
   ldX = ldY = ldBX = m;
   CHKERR(Num_zero_matrix_Sprimme(X, m, nXe - nXb, ldX, ctx));
   CHKERR(Num_zero_matrix_Sprimme(Y, m, nYe - nYb, ldY, ctx));
   CHKERR(Num_zero_matrix_Sprimme(BX, m, nBXe - nBXb, ldBX, ctx));

   int nGH = (G ? nG * nG : 0) + (H ? nH * nH : 0);
   if (ctx.numProcs > 1) {
      CHKERR(Num_malloc_SHprimme(nGH, &workGH, ctx));
      if (G) {
         G0 = workGH;
         ldG0 = nG;
      }
      if (H) {
         H0 = workGH + (G ? nG * nG : 0);
         ldH0 = nH;
      }
   }
   else {
      G0 = G; ldG0 = ldG;
      H0 = H; ldH0 = ldH;
   }

   if (Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = 0.0;
   if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = 0.0;
   if (G) CHKERR(Num_zero_matrix_SHprimme(G0, nG, nG, ldG0, ctx));
   if (H) CHKERR(Num_zero_matrix_SHprimme(H0, nH, nH, ldH0, ctx));
   if (xnorms) for (i=nxb; i<nxe; i++) xnorms[i-nxb] = 0.0;

   for (i=0; i < mV; i+=m, m=min(m,mV-i)) {
      /* X = V*h(nXb:nXe-1) */
      CHKERR(Num_gemm_dhd_Sprimme("N", "N", m, nXe-nXb, nV, 1.0,
         &V[i], ldV, &h[nXb*ldh], ldh, 0.0, X, ldX, ctx));

      /* X0 = X(nX0b-nXb:nX0e-nXb-1) */
      if (X0)
         CHKERR(Num_copy_matrix_Sprimme(
               &X[ldX * (nX0b - nXb)], m, nX0e - nX0b, ldX, &X0[i], ldX0, ctx));

      /* X1 = X(nX1b-nXb:nX1e-nXb-1) */
      if (X1) CHKERR(Num_copy_matrix_Sprimme(&X[ldX*(nX1b-nXb)], m, nX1e-nX1b,
            ldX, &X1[i], ldX1, ctx));

      /* X2 = X(nX2b-nXb:nX2e-nXb-1) */
      if (X2) CHKERR(Num_copy_matrix_Sprimme(&X[ldX*(nX2b-nXb)], m, nX2e-nX2b,
            ldX, &X2[i], ldX2, ctx));

      /* Y = W*h(nYb:nYe-1) */
      if (nYb < nYe) CHKERR(Num_gemm_dhd_Sprimme("N", "N", m, nYe-nYb, nV,
            1.0, &W[i], ldV, &h[nYb*ldh], ldh, 0.0, Y, ldY, ctx));

      /* Wo = Y(nWob-nYb:nWoe-nYb-1) */
      if (Wo) CHKERR(Num_copy_matrix_Sprimme(&Y[ldY*(nWob-nYb)], m, nWoe-nWob,
            ldY, &Wo[i], ldWo, ctx));

      /* BX = BV*h(nBXb:nBXe-1) */

      if (BV) CHKERR(Num_gemm_dhd_Sprimme("N", "N", m, nBXe-nBXb, nV, 1.0,
         &BV[i], ldV, &h[nBXb*ldh], ldh, 0.0, BX, ldBX, ctx));

      /* BX0 = BX(nX0b-nXb:nX0e-nXb-1) */
      if (BX0)
         CHKERR(Num_copy_matrix_Sprimme(&BX[ldBX * (nBX0b - nBXb)], m,
               nBX0e - nBX0b, ldBX, &BX0[i], ldBX0, ctx));

      /* BX1 = BX(nBX1b-nBXb:nBX1e-nBXb-1) */
      if (BX1) CHKERR(Num_copy_matrix_Sprimme(&BX[ldBX*(nBX1b-nBXb)], m, nBX1e-nBX1b,
            ldBX, &BX1[i], ldBX1, ctx));

      /* BX2 = BX(nBX2b-nBXb:nBX2e-nBXb-1) */
      if (BX2) CHKERR(Num_copy_matrix_Sprimme(&BX[ldBX*(nBX2b-nBXb)], m, nBX2e-nBX2b,
            ldBX, &BX2[i], ldBX2, ctx));

      /* G += X(:,0:nG-1)'*X(:,0:nG-1) */

      if (G) {
         CHKERR(Num_compute_gramm_ddh_Sprimme(X, m, nG, ldX, BV ? BX : X, ldX,
               i == 0 ? 0.0 : 1.0, G0, ldG0, 1 /* symmetric */, ctx));
      }

      /* H = X(:,0:nH-1)'*Y(:,0:nH-1) */

      if (H) {
         CHKERR(Num_compute_gramm_ddh_Sprimme(X, m, nH, ldX, Y, ldY,
               i == 0 ? 0.0 : 1.0, H0, ldH0,
               KIND(1 /*symmetric*/, 0 /*not symmetric*/), ctx));
      }

      /* xnorms = norm(X(nxb-nXb:nxe-nXb-1)) */
      if (xnorms) for (j=nxb; j<nxe; j++) {
            xnorms[j - nxb] += REAL_PART(Num_dot_Sprimme(
                  m, &X[ldX * (j - nXb)], 1, &X[ldX * (j - nXb)], 1, ctx));
      }

      /* R = Y(nRb-nYb:nRe-nYb-1) - BX(nRb-nYb:nRe-nYb-1)*diag(nRb:nRe-1) */
      if (R) {
         CHKERR(Num_compute_residuals_Sprimme(m, nRe - nRb, &hVals[nRb],
               BV ? &BX[ldBX * (nRb - nBXb)] : &X[ldX * (nRb - nXb)],
               BV ? ldBX : ldX, &Y[ldY * (nRb - nYb)], ldY, &R[i], ldR, ctx));
         if (Rnorms) {
            for (j = nRb; j < nRe; j++) {
               Rnorms[j - nRb] +=
                     REAL_PART(Num_dot_Sprimme(m, &R[i + ldR * (j - nRb)], 1,
                           &R[i + ldR * (j - nRb)], 1, ctx));
            }
         }
      }

      /* rnorms = Y(nrb-nYb:nre-nYb-1) - BX(nrb-nYb:nre-nYb-1)*diag(nrb:nre-1) */
      if (rnorms)  {
         CHKERR(Num_compute_residuals_Sprimme(m, nre - nrb, &hVals[nrb],
               BV ? &BX[ldBX * (nrb - nBXb)] : &X[ldX * (nrb - nXb)],
               BV ? ldBX : ldX, &Y[ldY * (nrb - nYb)], ldY,
               &Y[ldY * (nrb - nYb)], ldY, ctx));
         for (j = nrb; j < nre; j++) {
            rnorms[j - nrb] += REAL_PART(Num_dot_Sprimme(
                  m, &Y[ldY * (j - nYb)], 1, &Y[ldY * (j - nYb)], 1, ctx));
         }
      }
   }

   /* Reduce and copy back G0 and H0 */

   if (ctx.numProcs > 1) {
      CHKERR(globalSum_SHprimme(workGH, nGH, ctx));
   }
   if (G) CHKERR(Num_copy_matrix_SHprimme(G0, nG, nG, ldG0, G, ldG, ctx));
   if (H) CHKERR(Num_copy_matrix_SHprimme(H0, nH, nH, ldH0, H, ldH, ctx));

   /* Reduce Rnorms, rnorms and xnorms and sqrt the results */

   if (ctx.numProcs > 1) {
      HREAL *tmp;
      CHKERR(Num_malloc_RHprimme(nRe - nRb + nre - nrb + nxe - nxb, &tmp, ctx));
      j = 0;
      if (Rnorms) for (i=nRb; i<nRe; i++) tmp[j++] = Rnorms[i-nRb];
      if (rnorms) for (i=nrb; i<nre; i++) tmp[j++] = rnorms[i-nrb];
      if (xnorms) for (i=nxb; i<nxe; i++) tmp[j++] = xnorms[i-nxb];
      if (j) CHKERR(globalSum_RHprimme(tmp, j, ctx));
      j = 0;
      if (Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = sqrt(tmp[j++]);
      if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = sqrt(tmp[j++]);
      if (xnorms) for (i=nxb; i<nxe; i++) xnorms[i-nxb] = sqrt(tmp[j++]);
      CHKERR(Num_free_RHprimme(tmp, ctx));
   }
   else {
      if (Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = sqrt(Rnorms[i-nRb]);
      if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = sqrt(rnorms[i-nrb]);
      if (xnorms) for (i=nxb; i<nxe; i++) xnorms[i-nxb] = sqrt(xnorms[i-nxb]);
   }

   CHKERR(Num_free_Sprimme(X, ctx));
   CHKERR(Num_free_Sprimme(Y, ctx));
   CHKERR(Num_free_Sprimme(BX, ctx));
   CHKERR(Num_free_SHprimme(workGH, ctx));

   if (ctx.primme) ctx.primme->stats.timeDense += primme_wTimer() - t0;
   ctx.primme->stats.flopsDense +=
         mV * (nXe - nXb) * nV + mV * (nYe - nYb) * nV +
         (BV ? mV * (nBXe - nBXb) * nV : 0) + (G ? mV * nG * nG : 0) +
         (H ? mV * nH * nH : 0) + (xnorms ? (nxe - nxb) * mV : 0) +
         (R ? (nRe - nRb) * mV : 0) + ((R && Rnorms) ? (nRe - nRb) * mV : 0) +
         (R ? (nre - nrb) * mV : 0);

   return 0; 
}

/*******************************************************************************
 * Subroutine convTestFun - wrapper around primme.convTestFun; evaluate if the
 *    the approximate eigenpair eval, evec with given residual norm is
 *    considered as converged.
 *
 * INPUT PARAMETERS
 * ----------------
 * eval       the eigenvalue
 * evec       the eigenvector
 * givenEvec  whether eigenvector is provided
 * rNorm      the residual vector norm
 * 
 * OUTPUT
 * ------
 * isconv   if non-zero, the pair is considered converged.
 ******************************************************************************/

TEMPLATE_PLEASE
int convTestFun_Sprimme(HEVAL eval, SCALAR *evec, int givenEvec, HREAL rNorm,
      int *isconv, primme_context ctx) {

   primme_params *primme = ctx.primme;

   /* Cast eval and rNorm */

   KIND(double, PRIMME_COMPLEX_DOUBLE) evald = eval;
   double rNormd = rNorm;

   /* Cast evec if given */

   void *evec0 = NULL;
   if (evec && givenEvec)
      CHKERR(Num_matrix_astype_Sprimme(evec, primme->nLocal, 1, primme->nLocal,
            PRIMME_OP_SCALAR, &evec0, NULL, primme->convTestFun_type,
            1 /* alloc */, 1 /* copy */, ctx));

   /* If an evec is going to be passed to convTestFun, but nLocal is 0,       */
   /* then fake the evec with a nonzero pointer in order to not be mistaken   */
   /* by not passing a vector.                                                */
   SCALAR dummy;

   if (primme->nLocal == 0 && givenEvec) evec0 = &dummy;

   int ierr=0;
   CHKERRM((primme->convTestFun((double *)&evald, givenEvec ? evec : NULL,
                  &rNormd, isconv, primme, &ierr),
                 ierr),
         -1, "Error returned by 'convTestFun' %d", ierr);

   if (primme->nLocal > 0 && evec && givenEvec && evec != (SCALAR *)evec0)
      CHKERR(Num_free_Sprimme((SCALAR*)evec0, ctx));

   return 0;
}

TEMPLATE_PLEASE
int monitorFun_Sprimme(HEVAL *basisEvals, int basisSize, int *basisFlags,
      int *iblock, int blockSize, HREAL *basisNorms, int numConverged,
      HEVAL *lockedEvals, int numLocked, int *lockedFlags, HREAL *lockedNorms,
      int inner_its, HREAL LSRes, const char *msg, double time,
      primme_event event, double startTime, primme_context ctx) {

   /* Quick exit */

   primme_params *primme = ctx.primme;
   if (!primme->monitorFun) return 0;

   /* Cast basisEvals, basisNorms, lockedEvals, lockedNorms and LSRes */
   /* to monitorFun_type                                              */

   void *basisEvals0, *basisNorms0, *lockedEvals0, *lockedNorms0, *LSRes0;
   CHKERR(KIND(Num_matrix_astype_RHprimme, Num_matrix_astype_SHprimme)(
         basisEvals, 1, basisSize, 1, PRIMME_OP_HREAL, (void **)&basisEvals0,
         NULL, primme->monitorFun_type, 1 /* alloc */, 1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_RHprimme(basisNorms, 1, basisSize, 1,
         PRIMME_OP_HREAL, (void **)&basisNorms0, NULL, primme->monitorFun_type,
         1 /* alloc */, 1 /* copy */, ctx));
   CHKERR(KIND(Num_matrix_astype_RHprimme, Num_matrix_astype_SHprimme)(
         lockedEvals, 1, numLocked, 1, PRIMME_OP_HREAL, (void **)&lockedEvals0,
         NULL, primme->monitorFun_type, 1 /* alloc */, 1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_RHprimme(lockedNorms, 1, numLocked, 1,
         PRIMME_OP_HREAL, (void **)&lockedNorms0, NULL, primme->monitorFun_type,
         1 /* alloc */, 1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_RHprimme(&LSRes, 1, 1, 1,
         PRIMME_OP_HREAL, (void **)&LSRes0, NULL, primme->monitorFun_type,
         1 /* alloc */, 1 /* copy */, ctx));

   /* Call the user-defined functions */

   primme->stats.elapsedTime = primme_wTimer() - startTime;
   int err = 0;
   CHKERRM(
         (primme->monitorFun(basisEvals0, &basisSize, basisFlags, iblock,
                &blockSize, basisNorms0, &numConverged, lockedEvals0, &numLocked,
                lockedFlags, lockedNorms0, inner_its >= 0 ? &inner_its : NULL,
                LSRes >= 0 ? LSRes0 : NULL, msg, &time, &event, primme, &err),
               err),
         -1, "Error returned by monitorFun: %d", err);

   if (basisEvals != (HEVAL *)basisEvals0)
      CHKERR(KIND(Num_free_RHprimme, Num_free_SHprimme)(
            (HEVAL *)basisEvals0, ctx));
   if (basisNorms != (HREAL *)basisNorms0)
      CHKERR(Num_free_RHprimme((HREAL *)basisNorms0, ctx));
   if (lockedEvals != (HEVAL *)lockedEvals0)
      CHKERR(KIND(Num_free_RHprimme,Num_free_SHprimme)((HEVAL *)lockedEvals0, ctx));
   if (lockedNorms != (HREAL *)lockedNorms0)
      CHKERR(Num_free_RHprimme((HREAL *)lockedNorms0, ctx));
   if (&LSRes != (HREAL *)LSRes0)
      CHKERR(Num_free_RHprimme((HREAL *)LSRes0, ctx));

   return 0;
}

/******************************************************************************
 * Subroutine insertionSort -- This subroutine locks a converged Ritz value
 *   by insertion sorting it into the evals array.  A permutation array, perm,
 *   is maintained to keep track of the position the value would have been
 *   placed in had sorting not been performed.  This allows the locked Ritz
 *   vectors to be sorted at a later time upon return to the user.
 *   The order is ascending or descending for smallest/largest respectively.
 *   For interior, it is the order of convergence except for the same shifts.
 *   In that case, Ritz values that satisfy the criterion closer come first.
 *
 *
 * Input parameters
 * ----------------
 * newVal   The Ritz value to be locked
 *
 * newNorm  The residual norm of the Ritz value to be locked
 *
 * newFlag     The current flag
 *
 * n        The length of the input arrays: evals, resNorms, flags and perm
 *
 * initialShift  The index of the first targetShift to consider 
 * 
 * primme  Structure containing various solver parameters
 *
 *
 * Input/Output parameters
 * -----------------------
 * evals    The sorted list of locked Ritz values
 *
 * resNorms The residual norms corresponding to the locked Ritz values
 *
 * flags    The flags of the locked pairs
 *
 * perm     The permutation array indicating each Ritz values original
 *          unsorted position.
 *
 ******************************************************************************/

TEMPLATE_PLEASE int insertionSort_Sprimme(HEVAL newVal, HEVAL *evals,
      HREAL newNorm, HREAL *resNorms, int newFlag, int *flags, int *perm, int n,
      int initialShift, primme_params *primme) {

   int i, current; /* Indices used for sorting */
   HREAL ithShift, currentShift;

   /* ------------------------------------------------------------------ */
   /* Find smallest index to insert the Ritz value. The eigenvalue order */
   /* depends on how we target eigenvalues.                              */
   /* ------------------------------------------------------------------ */

#ifdef USE_HERMITIAN
   if (primme->target == primme_smallest) {

      for (i = n; i > 0; i--) {
         if (newVal >= evals[i - 1]) break;
      }
   } else if (primme->target == primme_largest) {

      for (i = n; i > 0; i--) {
         if (newVal <= evals[i - 1]) break;
      }
   } else
#endif /* USE_HERMITIAN */
   {
      /* For interior cases maintain convergence order except for the same shift
       * * Example: s1 s2 s2 s3 s4. Only eigenvalues converged for s2 may
       * switch.  *
       * Therefore, we only need to look back as long as the shift is the same.
       */

      currentShift = primme->targetShifts[min(
            primme->numTargetShifts - 1, initialShift + n)];

#ifdef USE_HERMITIAN
      if (primme->target == primme_closest_geq) {
         for (i = n; i > 0; i--) {
            ithShift = primme->targetShifts[min(
                  primme->numTargetShifts - 1, initialShift + i - 1)];
            if (ithShift != currentShift ||
                  newVal - currentShift >= evals[i - 1] - currentShift)
               break;
         }
      } else if (primme->target == primme_closest_leq) {
         for (i = n; i > 0; i--) {
            ithShift = primme->targetShifts[min(
                  primme->numTargetShifts - 1, initialShift + i - 1)];
            if (ithShift != currentShift ||
                  currentShift - newVal >= currentShift - evals[i - 1])
               break;
         }
      } else
#endif /* USE_HERMITIAN */
            if (primme->target == primme_closest_abs) {
         for (i = n; i > 0; i--) {
            ithShift = primme->targetShifts[min(
                  primme->numTargetShifts - 1, initialShift + i - 1)];
            if (ithShift != currentShift ||
                  EVAL_ABS(newVal - currentShift) >=
                        EVAL_ABS(evals[i - 1] - currentShift))
               break;
         }
      } else if (primme->target == primme_largest_abs) {
         for (i = n; i > 0; i--) {
            ithShift = primme->targetShifts[min(
                  primme->numTargetShifts - 1, initialShift + i - 1)];
            if (ithShift != currentShift ||
                  EVAL_ABS(newVal - currentShift) <=
                        EVAL_ABS(evals[i - 1] - currentShift))
               break;
         }
      } else {
         /* This should never happen */
         return PRIMME_FUNCTION_UNAVAILABLE;
      }
   }

   /* Shift the array to make room for the new Ritz value */

   for (current = n-1; current >= i; current--) {
      evals[current+1] = evals[current];
      if (resNorms) resNorms[current+1] = resNorms[current];
      if (perm)     perm[current+1] = perm[current];
      if (flags)    flags[current+1] = flags[current];
   }

   /* Insert the new value */

   evals[i] = newVal;
   if (resNorms) resNorms[i] = newNorm;
   if (perm)     perm[i] = n;
   if (flags)    flags[i] = newFlag;

   return 0;
}
 
#endif /* SUPPORTED_TYPE */
