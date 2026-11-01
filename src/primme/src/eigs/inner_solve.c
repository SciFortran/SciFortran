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
 * File: inner_solve.c
 *
 * Purpose - Solves the correction equation using hermitian simplified QMR.
 *  
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/inner_solve.c"
#endif

#include "numerical.h"
#include "template_normal.h"
#include "common_eigs.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "inner_solve.h"
#include "factorize.h"
#include "update_W.h"
#include "auxiliary_eigs.h"
#include "auxiliary_eigs_normal.h"
#endif

#ifdef SUPPORTED_TYPE

#ifdef USE_HERMITIAN

/*******************************************************************************
 * Function inner_solve - This subroutine solves the correction equation
 *    
 *           (I-BQQ)(I-Bxx)(A-shift*B)(I-xx'B)(I-QQ'B)sol = -r 
 *
 *    with Q = evecs, using hermitian simplified QMR.
 *    A preconditioner may be applied to this system to accelerate convergence.
 *    The preconditioner is assumed to approximate (A-shift*B)^{-1}.  The
 *    classical JD method as described in Templates for the Solution of 
 *    Eigenvalue Problems by Bai, et. al. requires that the preconditioner is 
 *    computed orthogonally to x and evecs. This code implements all 
 *    possible variations for projectors as defined by user parameters
 *    and setup_JD_projectors(). The QMR transparently calls the resulting
 *    projected matrix and preconditioner.
 *
 *
 * Input parameters
 * ----------------
 * x              The current Ritz vector for which the correction is being solved.
 *
 * r              The residual with respect to the Ritz vector.
 *
 * evecs          The converged Ritz vectors
 *
 * evecsHat       K^{-1}*B*evecs where K is a hermitian preconditioner.
 *
 * Mfact          The factors of the hermitian projection (evecs'*evecsHat). 
 *
 * ipivot         The pivoting for the Mfact factorization
 *
 * xKinvBx        The value x'*Kinv*B*x needed if skew-X projection
 *
 * LprojectorQ,   Points to an array that includes all the left projector.
 * LprojectorX    Can be any combination of [evecs x], [evecs], [x], NULL.
 *
 * RprojectorQ    Points to an array that includes the right skew projector for Q:
 *                It can be [evecsHat] or Null
 *
 * RprojectorX    Points to an array that includes the right skew projector for x:
 *                It can be [Kinvx] or Null
 *
 * sizeLprojectorQ   Number of colums of Lprojector
 *
 * sizeLprojectorX   Number of colums of Lprojector
 *
 * sizeRprojectorQ   Number of colums of RprojectorQ
 *
 * sizeRprojectorX   Number of colums of LprojectorX
 *
 * eval           The current Ritz value 
 *
 * shift          Correction eq. shift. The closer the shift is to the target 
 *                eigenvalue, the more accurate the correction will be.
 *
 * ctx         Structure containing various solver parameters
 *
 *
 * Input/Output parameters
 * -----------------------
 * r       The residual with respect to the Ritz vector.  May be altered upon
 *         return.
 * rnorm   On input, the 2 norm of r. No need to recompute it initially.
 *         On output, the estimated 2 norm of the updated eigenvalue residual
 * touch   Parameter used in inner solve stopping criteria
 * 
 * Output parameters
 * -----------------
 * sol   The solution (correction) of the correction equation
 *
 * Return Value
 * ------------
 * Error code
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int inner_solve_Sprimme(int blockSize, SCALAR *x, PRIMME_INT ldx, SCALAR *Bx,
      PRIMME_INT ldBx, SCALAR *r, PRIMME_INT ldr, HREAL *rnorm, SCALAR *evecs,
      PRIMME_INT ldevecs, HSCALAR *Mfact, int *ipivot, HSCALAR *xKinvBx,
      SCALAR *LprojectorQ, PRIMME_INT ldLprojectorQ, SCALAR *LprojectorX,
      PRIMME_INT ldLprojectorX, SCALAR *LprojectorBQ, PRIMME_INT ldLprojectorBQ,
      SCALAR *LprojectorBX, PRIMME_INT ldLprojectorBX, SCALAR *RprojectorQ,
      PRIMME_INT ldRprojectorQ, SCALAR *RprojectorX, PRIMME_INT ldRprojectorX,
      int sizeLprojectorQ, int sizeLprojectorX, int sizeRprojectorQ,
      int sizeRprojectorX, SCALAR *sol, PRIMME_INT ldsol, HEVAL *eval,
      KIND(double, PRIMME_COMPLEX_DOUBLE) * shift, int *touch, double startTime,
      primme_context ctx) {

   primme_params *primme = ctx.primme;
   int maxIterations; /* The maximum # iterations allowed. Depends on primme */

   /* QMR parameters */

   PRIMME_INT nLocal = primme->nLocal;
   SCALAR *g, *d, *delta, *w;
   CHKERR(Num_malloc_Sprimme(nLocal * blockSize, &g, ctx));
   CHKERR(Num_malloc_Sprimme(nLocal * blockSize, &d, ctx));
   CHKERR(Num_malloc_Sprimme(nLocal * blockSize, &delta, ctx));
   CHKERR(Num_malloc_Sprimme(nLocal * blockSize, &w, ctx));
   HREAL *sigma_prev, *rho_prev, *rho;
   CHKERR(Num_malloc_RHprimme(blockSize, &sigma_prev, ctx));
   CHKERR(Num_malloc_RHprimme(blockSize, &rho_prev, ctx));
   CHKERR(Num_malloc_RHprimme(blockSize, &rho, ctx));
   double *alpha_prev;
   CHKERR(Num_malloc_dprimme(blockSize, &alpha_prev, ctx));
   HREAL *Theta_prev, *Theta, *tau_init, *tau_prev, *tau; 
   CHKERR(Num_malloc_RHprimme(blockSize, &Theta_prev, ctx));
   CHKERR(Num_malloc_RHprimme(blockSize, &Theta, ctx));
   CHKERR(Num_malloc_RHprimme(blockSize, &tau_init, ctx));
   CHKERR(Num_malloc_RHprimme(blockSize, &tau_prev, ctx));
   CHKERR(Num_malloc_RHprimme(blockSize, &tau, ctx));

   /* Parameters used to dynamically update eigenpair */
   double *Beta_prev, *Delta_prev, *Psi_prev, *eta;
   double *eval_prev, *eres_updated;
   double *Gamma_prev, *Phi_prev;
   double *gamma;
   HREAL *normBx;
   CHKERR(Num_malloc_dprimme(blockSize, &Beta_prev, ctx));
   CHKERR(Num_malloc_dprimme(blockSize, &Delta_prev, ctx));
   CHKERR(Num_malloc_dprimme(blockSize, &Psi_prev, ctx));
   CHKERR(Num_malloc_dprimme(blockSize, &eta, ctx));
   CHKERR(Num_malloc_dprimme(blockSize, &eval_prev, ctx));
   CHKERR(Num_malloc_dprimme(blockSize, &eres_updated, ctx));
   CHKERR(Num_malloc_dprimme(blockSize, &Gamma_prev, ctx));
   CHKERR(Num_malloc_dprimme(blockSize, &Phi_prev, ctx));
   CHKERR(Num_malloc_dprimme(blockSize, &gamma, ctx));
   CHKERR(Num_malloc_RHprimme(blockSize, &normBx, ctx));

   /* Auxiliary arrays */
   HREAL *dot_sol = NULL, *Bnormsol = NULL;
   if (primme->correctionParams.convTest == primme_adaptive ||
         primme->correctionParams.convTest == primme_adaptive_ETolerance) {
      if (primme->massMatrixMatvec) {
         CHKERR(Num_malloc_RHprimme(blockSize, &Bnormsol, ctx));
      } else {
         CHKERR(Num_malloc_RHprimme(blockSize, &dot_sol, ctx));
      }
   }
   int *p, *p0; /* permutation of the right-hand-sides and auxiliary permutation */
   CHKERR(Num_malloc_iprimme(blockSize, &p, ctx));
   CHKERR(Num_malloc_iprimme(blockSize, &p0, ctx));
    
   double LTolerance, ETolerance, LTolerance_factor, ETolerance_factor;
   int isConv;
   int i;

   /* -----------------------------------------*/
   /* Set up convergence criteria by Tolerance */
   /* -----------------------------------------*/

   for (i=0; i<blockSize; i++) {
      tau_prev[i] = tau_init[i] = rnorm[i];       /* Assumes zero initial guess */
   }

   /* NOTE: In any case stop when linear system residual is less than         */
   /*       max(machEps,eps)*aNorm.                                           */
   double eps_matrix, eps_orth;
   CHKERR(machineEpsMatrix_Sprimme(&eps_matrix, ctx));
   CHKERR(machineEpsOrth_Sprimme(&eps_orth, ctx));
   LTolerance = max(eps_matrix, eps_orth) * problemNorm_Sprimme(1, primme);
   LTolerance_factor = 1.0;
   ETolerance = 0.0;
   ETolerance_factor = 0.0;

   switch(primme->correctionParams.convTest) {
   case primme_full_LTolerance:
      /* stop when linear system residual norm is less than problemNorm*eps.  */
      /* NOTE: the criterion is covered by the default values set before.     */
       break;
   case primme_decreasing_LTolerance:
      /* stop when linear system residual norm is less than relTolBase^-its   */
      LTolerance = max(LTolerance,
            pow(primme->correctionParams.relTolBase, 
               -(double)*touch));
      (*touch)++;
      break;
   case primme_adaptive:
      /* stop when estimate eigenvalue residual norm is less than             */  
      /* problemNorm*eps. Eigenresidual tol may not be achievable, because it */
      /* iterates on  P(A-s)P not on (A-s). But tau reflects the residual norm*/
      /* on P(A-s)P. So stop when linear system residual norm or the estimate */
      /* eigenvalue residual norm is less than problemNorm*eps/1.8.           */
      LTolerance_factor = pow(1.8, -(double)*touch);
      ETolerance_factor = pow(1.8, -(double)*touch);
      break; 
   case primme_adaptive_ETolerance:
      /* Besides the primme_adaptive criteria, stop when estimate eigenvalue  */
      /* residual norm is less than tau_init*0.1                              */
      LTolerance_factor = pow(1.8, -(double)*touch);
      ETolerance_factor = pow(1.8, -(double)*touch);
      ETolerance = 0.1;
     }
   
   /* --------------------------------------------------------*/
   /* Set up convergence criteria by max number of iterations */
   /* --------------------------------------------------------*/

   /* compute first total number of remaining matvecs */

   if (primme->maxMatvecs > 0) {
      maxIterations = primme->maxMatvecs - primme->stats.numMatvecs;
   }
   else {
      maxIterations = INT_MAX;
   }

   /* Perform primme.maxInnerIterations, but do not exceed total remaining */
   if (primme->correctionParams.maxInnerIterations > 0) {

      maxIterations = min(primme->correctionParams.maxInnerIterations, 
                          maxIterations);
   }

   /* --------------------------------------------------------*/
   /* Rest of initializations                                 */
   /* --------------------------------------------------------*/

   /* Assume zero initial guess */
   CHKERR(Num_copy_matrix_Sprimme(r, nLocal, blockSize, ldr, g, nLocal, ctx));

   CHKERR(apply_projected_preconditioner(g, nLocal, evecs, ldevecs, RprojectorQ,
         ldRprojectorQ, x, ldx, RprojectorX, ldRprojectorX, sizeRprojectorQ,
         sizeRprojectorX, xKinvBx, Mfact, ipivot, d, nLocal, blockSize, ctx));

   for (i=0; i<blockSize; i++) Theta_prev[i] = 0.0L;
   for (i=0; i<blockSize; i++) eval_prev[i] = eval[i];
   CHKERR(Num_dist_dots_real_Sprimme(
         g, nLocal, d, nLocal, nLocal, blockSize, rho_prev, ctx));

   /* Initialize recurrences used to dynamically update the eigenpair */

   for (i=0; i<blockSize; i++) Beta_prev[i] = Delta_prev[i] = Psi_prev[i] = 0.0;
   for (i=0; i<blockSize; i++) Gamma_prev[i] = Phi_prev[i] = eres_updated[i] = 0.0;

   /* other initializations */

   Num_zero_matrix_Sprimme(delta, nLocal, blockSize, nLocal, ctx);
   Num_zero_matrix_Sprimme(sol, nLocal, blockSize, ldsol, ctx);

   if ((primme->correctionParams.convTest == primme_adaptive ||
             primme->correctionParams.convTest == primme_adaptive_ETolerance) &&
         primme->massMatrixMatvec) {
      CHKERR(Num_dist_dots_real_Sprimme(
            Bx, ldBx, Bx, ldBx, nLocal, blockSize, normBx, ctx));
   } else {
      for (i = 0; i < blockSize; i++) normBx[i] = 1.0;
   }

   /*----------------------------------------------------------------------*/
   /*------------------------ Begin Inner Loop ----------------------------*/
   /*----------------------------------------------------------------------*/

   for (i=0; i<blockSize; i++) p[i] = i;
   int numIts;        /* Number of inner iterations                          */
   for (numIts = 0; numIts < maxIterations && blockSize > 0; numIts++) {

      CHKERR(apply_projected_matrix(d, nLocal, shift, LprojectorQ,
            ldLprojectorQ, sizeLprojectorQ, LprojectorBQ, ldLprojectorBQ,
            LprojectorX, ldLprojectorX, LprojectorBX,
            ldLprojectorBX, sizeLprojectorX, blockSize, w, nLocal, ctx));
      CHKERR(Num_dist_dots_real_Sprimme(
            d, nLocal, w, nLocal, nLocal, blockSize, sigma_prev, ctx));

      int conv;
      for (i=0; i<blockSize; i++) p0[i] = i;
      for (i = conv = 0; i < blockSize; i++) {
         if (!ISFINITE(sigma_prev[p[i]]) || sigma_prev[p[i]] == 0.0L) {
            PRINTF(5, "Exiting because SIGMA %e in block vector %d",
                  sigma_prev[p[i]], p[i]);

            /* sol = r if first iteration */
            if (numIts == 0) {
               CHKERR(Num_copy_matrix_Sprimme(
                     &r[ldr * i], nLocal, 1, ldr, &sol[ldsol * i], ldsol, ctx));
            }
            CHKERR(perm_set_value_on_pos(p0, i, blockSize - ++conv, blockSize));
            continue;
         }

         alpha_prev[p[i]] = rho_prev[p[i]]/sigma_prev[p[i]];
         if (!ISFINITE(alpha_prev[p[i]]) ||
               fabs(alpha_prev[p[i]]) < MACHINE_EPSILON ||
               fabs(alpha_prev[p[i]]) > 1.0L / MACHINE_EPSILON) {
            PRINTF(5, "Exiting because ALPHA %e in block vector %d",
                  alpha_prev[p[i]], p[i]);

            /* sol = r if first iteration */
            if (numIts == 0) {
               CHKERR(Num_copy_matrix_Sprimme(
                     &r[ldr * i], nLocal, 1, ldr, &sol[ldsol * i], ldsol, ctx));
            }
            CHKERR(perm_set_value_on_pos(p0, i, blockSize - ++conv, blockSize));
            continue;
         }

         Num_axpy_Sprimme(nLocal, -alpha_prev[p[i]], &w[nLocal * i], 1,
               &g[nLocal * i], 1, ctx);
      }

      /* Apply permutation p0 and shrink blockSize */
      CHKERR(permute_vecs_iprimme(p, blockSize, p0, ctx));
      CHKERR(permute_vecs_dprimme(shift, 1, blockSize, 1, p0, ctx));
      CHKERR(permute_vecs_SHprimme(xKinvBx, 1, blockSize, 1, p0, ctx));
      CHKERR(permute_vecs_Sprimme(g, nLocal, blockSize, nLocal, p0, ctx));
      CHKERR(permute_vecs_Sprimme(d, nLocal, blockSize, nLocal, p0, ctx));
      CHKERR(permute_vecs_Sprimme(delta, nLocal, blockSize, nLocal, p0, ctx));
      if (sizeLprojectorX) CHKERR(permute_vecs_Sprimme(LprojectorX, nLocal, blockSize, nLocal, p0, ctx));
      if (sizeRprojectorX) CHKERR(permute_vecs_Sprimme(RprojectorX, nLocal, blockSize, nLocal, p0, ctx));
      CHKERR(permute_vecs_Sprimme(r, nLocal, blockSize, ldr, p0, ctx));
      CHKERR(permute_vecs_Sprimme(x, nLocal, blockSize, ldx, p0, ctx));
      CHKERR(permute_vecs_Sprimme(sol, nLocal, blockSize, ldsol, p0, ctx));
      blockSize -= conv;
      if (sizeLprojectorX) sizeLprojectorX -= conv;
      if (sizeRprojectorX) sizeRprojectorX -= conv;
      if (blockSize <= 0) break;

      CHKERR(Num_dist_dots_real_Sprimme(
            g, nLocal, g, nLocal, nLocal, blockSize, Theta, ctx));

      for (i = 0; i < blockSize; i++) {
         Theta[p[i]] = sqrt(Theta[p[i]]) / tau_prev[p[i]];
         double c = 1.0/sqrt(1+Theta[p[i]]*Theta[p[i]]);
         tau[p[i]] = tau_prev[p[i]]*Theta[p[i]]*c;

         gamma[p[i]] = c*c*Theta_prev[p[i]]*Theta_prev[p[i]];
         eta[p[i]] = alpha_prev[p[i]] * c * c;

#ifdef USE_HOST
         int j;
         if (dot_sol) dot_sol[i] = 0.0;
         for (j = 0; j < nLocal; j++) {
            SET_COMPLEX(delta[i * nLocal + j],
                  TO_COMPLEX(delta[i * nLocal + j]) * (HSCALAR)gamma[p[i]] +
                        TO_COMPLEX(d[nLocal * i + j]) * (HSCALAR)eta[p[i]]);
            SET_COMPLEX(
                  sol[ldsol * i + j], TO_COMPLEX(delta[nLocal * i + j]) +
                                            TO_COMPLEX(sol[ldsol * i + j]));
            if (dot_sol)
               dot_sol[i] += REAL_PART(CONJ(TO_COMPLEX(sol[ldsol * i + j])) *
                                       TO_COMPLEX(sol[ldsol * i + j]));
         }
#else
         Num_scal_Sprimme(
               nLocal, (HSCALAR)gamma[p[i]], &delta[i * nLocal], 1, ctx);
         Num_axpy_Sprimme(nLocal, (HSCALAR)eta[p[i]], &d[i * nLocal], 1,
               &delta[i * nLocal], 1, ctx);
         Num_axpy_Sprimme(
               nLocal, 1.0, &delta[i * nLocal], 1, &sol[i * ldsol], 1, ctx);
         if (dot_sol)
            dot_sol[i] = REAL_PART(Num_dot_Sprimme(
                  nLocal, &sol[i * ldsol], 1, &sol[i * ldsol], 1, ctx));
#endif
      }

      if (dot_sol) CHKERR(globalSum_RHprimme(dot_sol, blockSize, ctx));

      /* Compute B-norm of sol if adapting stopping and a generalized problem is
       * being solved */

      if (Bnormsol) {
         CHKERR(massMatrixMatvec_Sprimme(
               sol, ldsol, nLocal, w, nLocal, 0, blockSize, ctx));

         CHKERR(Num_dist_dots_real_Sprimme(
               sol, ldsol, w, nLocal, nLocal, blockSize, Bnormsol, ctx));
      }

      for (i=0; i<blockSize; i++) p0[i] = i;
      for (i = conv = 0; i < blockSize; i++) {
         if (fabs(rho_prev[p[i]]) == 0.0L ) {
            PRINTF(5, "Exiting because abs(rho) %e in block vector %d",
                  fabs(rho_prev[p[i]]), p[i]);
            CHKERR(perm_set_value_on_pos(p0, i, blockSize - ++conv, blockSize));
            continue;
         }
      
         if (numIts > 0 && tau[p[i]] < LTolerance) {
            PRINTF(5, " tau < LTol %e %e in block vector %d", tau[p[i]],
                  LTolerance, p[i]);
            CHKERR(perm_set_value_on_pos(p0, i, blockSize - ++conv, blockSize));
            continue;
         }
         if (ETolerance > 0.0 || ETolerance_factor > 0.0) {
            /* --------------------------------------------------------*/
            /* Adaptive stopping based on dynamic monitoring of eResid */
            /* --------------------------------------------------------*/

            /* Update the Ritz value and eigenresidual using the */
            /* following recurrences.                            */

            double Delta = gamma[p[i]] * Delta_prev[p[i]] + eta[p[i]] * rho_prev[p[i]];
            double Beta = Beta_prev[p[i]] - Delta;
            double Phi = gamma[p[i]] * gamma[p[i]] * Phi_prev[p[i]] + eta[p[i]] * eta[p[i]] * sigma_prev[p[i]];
            double Psi = gamma[p[i]] * Psi_prev[p[i]] + gamma[p[i]] * Phi_prev[p[i]];
            double Gamma = Gamma_prev[p[i]] + 2.0L * Psi + Phi;

            /* Perform the update: update the eigenvalue and the square of the
             * residual norm */

            double Bnorm_x_plus_sol =
                  1.0 + (primme->massMatrixMatvec ? Bnormsol[i] : dot_sol[i]);
            double eval_updated =
                  shift[i] +
                  (eval[p[i]] - shift[i] + 2 * Beta + Gamma) / Bnorm_x_plus_sol;
            double eres2_updated =
                  ((double)tau[p[i]] * tau[p[i]]) / Bnorm_x_plus_sol +
                  (normBx[p[i]]*((double)eval[p[i]] - shift[i] + Beta) *
                        ((double)eval[p[i]] - shift[i] + Beta)) /
                        Bnorm_x_plus_sol -
                  (eval_updated - shift[i]) * (eval_updated - shift[i]);

            /* If numerical problems, let eres about the same as tau */
            double eres_prev = eres_updated[p[i]];
            if (eres2_updated < 0) {
               eres_updated[p[i]] =
                     sqrt(((double)tau[p[i]] * tau[p[i]]) / Bnorm_x_plus_sol);
            }
            else 
               eres_updated[p[i]] = sqrt(eres2_updated);

         
            Delta_prev[p[i]] = Delta;
            Beta_prev[p[i]] = Beta;
            Phi_prev[p[i]] = Phi;
            Psi_prev[p[i]] = Psi;
            Gamma_prev[p[i]] = Gamma;

            assert(ISFINITE(Delta) && ISFINITE(Beta) && ISFINITE(Phi)
                  && ISFINITE(Psi) && ISFINITE(Gamma) && ISFINITE(eval_updated)
                  && ISFINITE(eres2_updated) && ISFINITE(eres_updated[p[i]]));

            /* --------------------------------------------------------*/
            /* Stopping criteria                                       */
            /* --------------------------------------------------------*/

            if (numIts > 0 && (tau_prev[p[i]] <= eres_updated[p[i]] ||
                                    eres_prev <= tau[p[i]])) {
               PRINTF(5, " tau < R eres for block vector %d", p[i]);
               CHKERR(perm_set_value_on_pos(p0, i, blockSize - ++conv, blockSize));
               continue;
            }

            if (primme->target == primme_smallest &&
                  eval_updated > eval_prev[p[i]]) {
               PRINTF(5, "eval_updated > eval_prev in block vector %d", p[i]);
               CHKERR(perm_set_value_on_pos(p0, i, blockSize - ++conv, blockSize));
               continue;
            }
            else if (primme->target == primme_largest && eval_updated < eval_prev[p[i]]){
               PRINTF(5, "eval_updated < eval_prev in block vector %d", p[i]);
               CHKERR(perm_set_value_on_pos(p0, i, blockSize - ++conv, blockSize));
               continue;
            } else if (primme->target == primme_closest_abs &&
                       fabs(eval[p[i]] - eval_updated) >
                             tau_init[p[i]] + eres_updated[p[i]]) {
               PRINTF(5, "|eval-eval_updated| > tau0+eres in block vector %d",
                     p[i]);
               CHKERR(perm_set_value_on_pos(p0, i, blockSize - ++conv, blockSize));
               continue;
            }
          
            if (numIts > 0 && eres_updated[p[i]] < ETolerance*tau_init[p[i]]) {
               PRINTF(5, "eres < eresTol %e in block vector %d",
                     eres_updated[p[i]], p[i]);
               CHKERR(perm_set_value_on_pos(p0, i, blockSize - ++conv, blockSize));
               continue;
            }

            /* Check if some of the next conditions is satisfied:             */
            /* a) estimate eigenvalue residual norm (eres_updated) is less    */
            /*    than eps*aNorm*Etolerance_factor                            */
            /* b) linear system residual norm is less                         */
            /*    than eps*aNorm*LTolerance_factor                            */
            /* The result is to check if eps*aNorm is less than               */
            /* max(tau/LTolerance_factor, eres_updated/ETolerance_factor).    */

            double tol = min(tau[p[i]] / LTolerance_factor,
                  eres_updated[p[i]] / ETolerance_factor);
            CHKERR(convTestFun_Sprimme(eval_updated, NULL,
                  0 /* evec not given */, tol, &isConv, ctx));

            if (numIts > 0 && isConv) {
               PRINTF(5,
                     "eigenvalue and residual norm passed convergence "
                     "criterion in block vector %d",
                     p[i]);
               (*touch)++;
               CHKERR(perm_set_value_on_pos(p0, i, blockSize - ++conv, blockSize));
               continue;
            }

            eval_prev[p[i]] = eval_updated;

            /* Report inner iteration */
            if (primme->monitorFun) {
               int ZERO = 0, UNCO = UNCONVERGED;
               HREAL evalr = eval_updated, resr = eres_updated[p[i]];

               CHKERR(monitorFun_Sprimme(&evalr, 1, &UNCO, &ZERO, 1, &resr, -1,
                     NULL, -1, NULL, NULL, numIts, tau[p[i]], NULL, 0.0,
                     primme_event_inner_iteration, startTime, ctx));
            }

           /* --------------------------------------------------------*/
         } /* End of if adaptive JDQMR section                        */
           /* --------------------------------------------------------*/
         else
         {
            /* Check if the linear system residual norm (tau) is less         */
            /* than eps*aNorm*LTolerance_factor. Note that QMR residual can   */
            /* be sqrt(iterations) times away from the actual residual.       */

            CHKERR(convTestFun_Sprimme(eval[p[i]], NULL, 0 /* evec not given */,
                  tau[p[i]] / LTolerance_factor * sqrt((double)numIts), &isConv,
                  ctx));

            if (numIts > 0 && isConv) {
               PRINTF(5,
                     "eigenvalue and residual norm "
                     "passed convergence criterion in block vector %d",
                     p[i]);
               CHKERR(perm_set_value_on_pos(p0, i, blockSize - ++conv, blockSize));
               continue;
            }

            else if (primme->monitorFun) {
               /* Report for non adaptive inner iterations */
               int ZERO = 0, UNCO = UNCONVERGED;
               CHKERR(monitorFun_Sprimme(&eval[p[i]], 1, &UNCO, &ZERO, 1,
                     &rnorm[p[i]], 0, NULL, 0, NULL, NULL, numIts, tau[p[i]],
                     NULL, 0.0, primme_event_inner_iteration, startTime, ctx));
            }
         }
      }

      /* Apply permutation p0 and shrink blockSize */
      CHKERR(permute_vecs_iprimme(p, blockSize, p0, ctx));
      CHKERR(permute_vecs_dprimme(shift, 1, blockSize, 1, p0, ctx));
      CHKERR(permute_vecs_SHprimme(xKinvBx, 1, blockSize, 1, p0, ctx));
      CHKERR(permute_vecs_Sprimme(g, nLocal, blockSize, nLocal, p0, ctx));
      CHKERR(permute_vecs_Sprimme(d, nLocal, blockSize, nLocal, p0, ctx));
      CHKERR(permute_vecs_Sprimme(delta, nLocal, blockSize, nLocal, p0, ctx));
      if (sizeLprojectorX) CHKERR(permute_vecs_Sprimme(LprojectorX, nLocal, blockSize, nLocal, p0, ctx));
      if (sizeRprojectorX) CHKERR(permute_vecs_Sprimme(RprojectorX, nLocal, blockSize, nLocal, p0, ctx));
      CHKERR(permute_vecs_Sprimme(r, nLocal, blockSize, ldr, p0, ctx));
      CHKERR(permute_vecs_Sprimme(x, nLocal, blockSize, ldx, p0, ctx));
      CHKERR(permute_vecs_Sprimme(sol, nLocal, blockSize, ldsol, p0, ctx));
      blockSize -= conv;
      if (sizeLprojectorX) sizeLprojectorX -= conv;
      if (sizeRprojectorX) sizeRprojectorX -= conv;
      if (blockSize <= 0) break;

      if (numIts + 1 < maxIterations) {

         CHKERR(apply_projected_preconditioner(g, nLocal, evecs, ldevecs,
               RprojectorQ, ldRprojectorQ, x, ldx, RprojectorX, ldRprojectorX,
               sizeRprojectorQ, sizeRprojectorX, xKinvBx, Mfact, ipivot, w, nLocal,
               blockSize, ctx));

         CHKERR(Num_dist_dots_real_Sprimme(
               g, nLocal, w, nLocal, nLocal, blockSize, rho, ctx));

         for (i=0; i< blockSize; i++) {
            HREAL beta = rho[p[i]]/rho_prev[p[i]];
            Num_axpy_Sprimme(
                  nLocal, beta, &d[nLocal * i], 1, &w[nLocal * i], 1, ctx);

            rho_prev[p[i]] = rho[p[i]];
            tau_prev[p[i]] = tau[p[i]];
            Theta_prev[p[i]] = Theta[p[i]];
         }

         /* Alternate between w and d buffers in successive iterations
          * This saves a memory copy. */
         SCALAR *ptmp = d;
         d = w;
         w = ptmp;
      }

     /* --------------------------------------------------------*/
   } /* End of QMR main while loop                              */
     /* --------------------------------------------------------*/

   CHKERR(Num_free_Sprimme(g, ctx));
   CHKERR(Num_free_Sprimme(d, ctx));
   CHKERR(Num_free_Sprimme(delta, ctx));
   CHKERR(Num_free_Sprimme(w, ctx));
   CHKERR(Num_free_RHprimme(sigma_prev, ctx));
   CHKERR(Num_free_RHprimme(rho_prev, ctx));
   CHKERR(Num_free_RHprimme(rho, ctx));
   CHKERR(Num_free_dprimme(alpha_prev, ctx));
   CHKERR(Num_free_RHprimme(Theta_prev, ctx));
   CHKERR(Num_free_RHprimme(Theta, ctx));
   CHKERR(Num_free_RHprimme(tau_init, ctx));
   CHKERR(Num_free_RHprimme(tau_prev, ctx));
   CHKERR(Num_free_RHprimme(tau, ctx));
   CHKERR(Num_free_dprimme(Beta_prev, ctx));
   CHKERR(Num_free_dprimme(Delta_prev, ctx));
   CHKERR(Num_free_dprimme(Psi_prev, ctx));
   CHKERR(Num_free_dprimme(eta, ctx));
   CHKERR(Num_free_dprimme(eval_prev, ctx));
   CHKERR(Num_free_dprimme(eres_updated, ctx));
   CHKERR(Num_free_dprimme(Gamma_prev, ctx));
   CHKERR(Num_free_dprimme(Phi_prev, ctx));
   CHKERR(Num_free_dprimme(gamma, ctx));
   CHKERR(Num_free_RHprimme(normBx, ctx));
   CHKERR(Num_free_RHprimme(Bnormsol, ctx));
   CHKERR(Num_free_RHprimme(dot_sol, ctx));
   CHKERR(Num_free_iprimme(p, ctx));
   CHKERR(Num_free_iprimme(p0, ctx));

   return 0;
}
   

/*******************************************************************************
 * Function apply_projected_preconditioner - This routine applies the
 *    projected preconditioner to a vector v by computing:
 *
 *         result = (I-KinvBx/xKinvBx*x') (I - Qhat (Q'*Qhat)^{-1}Q') Kinv*v
 *
 *    First we apply the preconditioner Kinv*v, and then the two projectors 
 *    are computed one after the other.
 *    
 * Input Parameters
 * ----------------
 * v      The vector the projected preconditioner will be applied to.
 *
 * Q      The matrix evecs where evecs are the locked/converged eigenvectors
 *
 * RprojectorQ     The matrix K^{-1}BQ (often called Qhat), Q, or nothing,
 *                 as determined by setup_JD_projectors.
 *
 * x               The current Ritz vector.
 *
 * RprojectorX     The matrix K^{-1}Bx (if needed)
 *
 * sizeRprojectorQ The number of columns in RprojectorQ
 *
 * sizeRprojectorX The number of columns in RprojectorX
 *
 * xKinvBx The value x^T (Kinv*B*x). It is computed in the setup_JD_projectors
 *
 * Mfact  The factorization of (Q'*K^{-1}*B*Q).
 *
 * ipivot Permutation array indicating how the rows of the Mfact decomposition
 *        have been pivoted.
 *
 * primme   Structure containing various solver parameters.
 *
 *
 * Output parameters
 * -----------------
 * result The result of the application.
 *
 ******************************************************************************/

STATIC int apply_projected_preconditioner(SCALAR *v, PRIMME_INT ldv, SCALAR *Q,
      PRIMME_INT ldQ, SCALAR *RprojectorQ, PRIMME_INT ldRprojectorQ, SCALAR *x,
      PRIMME_INT ldx, SCALAR *RprojectorX, PRIMME_INT ldRprojectorX,
      int sizeRprojectorQ, int sizeRprojectorX, HSCALAR *xKinvBx,
      HSCALAR *Mfact, int *ipivot, SCALAR *result, PRIMME_INT ldresult,
      int blockSize, primme_context ctx) {

   assert(sizeRprojectorX == 0 || sizeRprojectorX == blockSize);

   /* Place K^{-1}v in result */
   primme_params *primme = ctx.primme;
   CHKERR(applyPreconditioner_Sprimme(v, primme->nLocal, ldv, result,
            ldresult, blockSize, ctx));

   CHKERR(apply_skew_projector(Q, ldQ, RprojectorQ, ldRprojectorQ, Mfact, ipivot,
            sizeRprojectorQ, result, ldresult, blockSize, ctx));

   if (sizeRprojectorX <= 0) return 0;

   int i;
   for (i=0; i<blockSize; i++) {
      CHKERR(apply_skew_projector(&x[ldx * i], ldx,
            &RprojectorX[ldRprojectorX * i], ldRprojectorX, &xKinvBx[i], NULL,
            1, &result[ldresult * i], ldresult, 1, ctx));
   }

   return 0;
}

/*******************************************************************************
 * Subroutine apply_skew_projector - Apply the skew projector to a vector v:
 *
 *     v = (I-Qhat*inv(Q'Qhat)*Q') v
 *
 *   The result is placed back in v.  Q is the matrix of converged Ritz 
 *   vectors or the current Ritz vector.
 *
 * Input Parameters
 * ----------------
 * Q       The matrix of converged Ritz vectors and the current Ritz vector
 *
 * Qhat    The matrix of K^{-1}BQ
 *
 * Mfact   The factorization of the (Q'*Qhat) matrix
 *
 * ipivot  The pivot array for the Mfact factorization
 *
 * numCols Number of columns of Q and Qhat
 *
 * Input/Output Parameters
 * -----------------------
 * v       The vector to be skewed orthogonalized 
 * 
 ******************************************************************************/

STATIC int apply_skew_projector(SCALAR *Q, PRIMME_INT ldQ, SCALAR *Qhat,
      PRIMME_INT ldQhat, HSCALAR *Mfact, int *ipivot, int numCols, SCALAR *v,
      PRIMME_INT ldv, int blockSize, primme_context ctx) {

   primme_params *primme = ctx.primme;

   if (numCols <= 0 || blockSize <= 0) return 0;

   double t0 = primme_wTimer();

   HSCALAR *overlaps; /* overlaps of v with columns of Q   */
   CHKERR(Num_malloc_SHprimme(numCols * blockSize, &overlaps, ctx));

   /* Compute workspace = Q'*v */
   CHKERR(Num_gemm_ddh_Sprimme("C", "N", numCols, blockSize, primme->nLocal,
         1.0, Q, ldQ, v, ldv, 0.0, overlaps, numCols, ctx));
   if (primme) primme->stats.numOrthoInnerProds += numCols * blockSize;

   /* Global sum: overlaps = Q'*v */
   CHKERR(globalSum_SHprimme(overlaps, numCols * blockSize, ctx));

   /* Backsolve only if there is a skew projector */
   if (Mfact != NULL) {
      /* Solve (Q'Qhat)^{-1}*overlaps = overlaps = Q'*v for alpha by */
      /* backsolving  with Mfact.                 */

      CHKERR(MSolve_SHprimme(Mfact, ipivot, numCols, overlaps, blockSize,
            numCols, overlaps, numCols, ctx));
   }

   /* Compute v=v-Qhat*overlaps */
   CHKERR(Num_gemm_dhd_Sprimme("N", "N", primme->nLocal, blockSize, numCols,
         -1.0, Qhat, ldQhat, overlaps, numCols, 1.0, v, ldv, ctx));

   CHKERR(Num_free_SHprimme(overlaps, ctx));

   if (primme) primme->stats.timeOrtho += primme_wTimer() - t0;

   return 0;
}


/*******************************************************************************
 * Subroutine apply_projected_matrix - This subroutine applies the 
 *    projected matrix (I-BX*X)*(I-BQ*Q)*(A-shift*B) to a vector v
 *
 * Input Parameters
 * ----------------
 * v      The vector the projected matrix will be applied to
 *
 * shift  The amount the matrix is shifted by.
 *
 * Q      The converged Ritz vectors and the current Ritz vector
 *
 * BQ     B*Q
 *
 * dimQ   The number of columns of Q
 * 
 * rwork  Workspace of size 2*dimQ
 *
 * primme   Structure containing various solver parameters
 *
 *
 * Output Parameters
 * -----------------
 * result The result of the application.
 *
 ******************************************************************************/

STATIC int apply_projected_matrix(SCALAR *v, PRIMME_INT ldv, double *shift,
      SCALAR *Q, PRIMME_INT ldQ, int nQ, SCALAR *BQ, PRIMME_INT ldBQ, SCALAR *X,
      PRIMME_INT ldX, SCALAR *BX, PRIMME_INT ldBX, int nX, int blockSize,
      SCALAR *result, PRIMME_INT ldresult, primme_context ctx) {

   assert(nX == 0 || nX == blockSize);
   primme_params *primme = ctx.primme;

   /* result = A * v */

   CHKERR(matrixMatvec_Sprimme(
         v, primme->nLocal, ldv, result, ldresult, 0, blockSize, ctx));

   /* Bv = B * v */

   SCALAR *Bv;
   PRIMME_INT ldBv;
   if (primme->massMatrixMatvec) {
      ldBv = primme->ldOPs;
      CHKERR(Num_malloc_Sprimme(ldBv * blockSize, &Bv, ctx));
      CHKERR(massMatrixMatvec_Sprimme(
            v, primme->nLocal, ldv, Bv, ldBv, 0, blockSize, ctx));
   } else {
      ldBv = ldv;
      Bv = v;
   }

   /* result -= shift * Bv */

   int i;
   for (i = 0; i < blockSize; i++) {
      Num_axpy_Sprimme(primme->nLocal, -shift[i], &Bv[ldBv * i], 1,
            &result[ldresult * i], 1, ctx);
   }

   if (primme->massMatrixMatvec) {
      CHKERR(Num_free_Sprimme(Bv, ctx));
   }

   /* result = (I-BQ*Q')*result */

   CHKERR(apply_skew_projector(
         Q, ldQ, BQ, ldBQ, NULL, NULL, nQ, result, ldresult, blockSize, ctx));

   /* result = (I-BX*X)*result for each vector */

   if (nX <= 0) return 0;
   for (i = 0; i < blockSize; i++) {
      CHKERR(apply_skew_projector(&X[ldX * i], ldX, &BX[ldBX * i], ldBX, NULL,
            NULL, 1, &result[ldresult * i], ldresult, 1, ctx));
   }

   return 0;
}

/*******************************************************************************
 * Subroutine perm_set_value_on_pos - find the position of the value on the
 *    vector and swap with the given position.
 *
 * Input/Output Parameters
 * -----------------------
 * p      The input vector (it is modified)
 *
 * val    The value to seek
 *
 * pos    The position where the value is going to be moved
 *
 * n      The length of the vector
 *
 * Output Parameters
 * -----------------
 * error code
 *
 ******************************************************************************/

STATIC int perm_set_value_on_pos(int *p, int val, int pos, int n) {

   int i;
   for (i=0; i<n; i++) {
      if (p[i] == val) {
         p[i] = p[pos];
         p[pos] = val;
         return 0;
      }
   }

   return -1;
}

#endif /* USE_HERMITIAN */

#endif /* SUPPORTED_TYPE */
