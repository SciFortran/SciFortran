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
 * File: init.c
 *
 * Purpose - Generate the basis that will be used during the first
 *           iteration of the method.
 *  
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/init.c"
#endif

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "numerical.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "init.h"
#include "update_W.h"
#include "ortho.h"
#include "factorize.h"
#include "auxiliary_eigs.h"
#endif

#ifdef SUPPORTED_TYPE

/*******************************************************************************
 * subroutine init_basis - This subroutine is used to 
 *    initialize the basis V.
 *
 * Here we initialize the basis V as well as the array of locked
 * vectors if locking has been enabled.  The cases handled are:
 * 
 *    I. No initial vectors are provided.  BlockSize orthonormal initial
 *       vectors are created and used to form an orthonormal block Krylov
 *       basis.  The size of the basis will be minRestartSize.
 * 
 *   II. Initial vectors are provided.
 *
 *       1.a if locking is disabled, minRestartSize or greater number of initial
 *          vectors are provided and they are orthonormalized by calling ortho.
 * 
 *       1.b if locking is enable, up to minRestartSize of initial vectors are
 *          copied to V and orthonormalized. 
 *       2. There are fewer than minRestartSize initial vectors provided.
 *          A Krylov subspace of dimension restartSize - initSize vectors
 *          is created so that restartSize initial vectors will be available.
 * 
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * evecs      Array containing the orthogonalization constraints and initial 
 *            guesses. Holds as many as primme->numOrthoConst + primme->numEvals
 *            vectors.
 *
 * machEps    double machine precision
 * 
 * ctx          Structure containing various solver parameters
 * 
 * 
 * 
 * OUTPUT PARAMETERS
 * -----------------
 * *basisSize   The size of the resulting basis V
 *
 * *nextGuess   The index of the next initial guess stored in the evecs array
 *
 * *numGuesses  When locking is enabled, the number of remaining initial guesses
 * 
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * V            The orthonormal basis
 *
 * W            A*V
 *
 * evecsHat     K^{-1}*B*evecs, given a preconditioner K
 *
 * M            evecs'*evecsHat.  Its dimension is as large as 
 *              (primme->numOrthoConst + primme->numEvals).
 *
 * Mfact        The factorization of M
 *
 * ipivot       The pivots of the Mfact factorization
 *
 * VtBV         V'*B*V (used by Bortho_block)
 *
 * fVtBV        The Cholesky factor of VtBV (used by Bortho_block)
 *
 * Return value
 * ------------
 *  error code
 ******************************************************************************/

TEMPLATE_PLEASE
int init_basis_Sprimme(SCALAR *V, PRIMME_INT nLocal, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, SCALAR *BV, PRIMME_INT ldBV, SCALAR *evecs,
      PRIMME_INT ldevecs, SCALAR *Bevecs, PRIMME_INT ldBevecs, SCALAR *evecsHat,
      PRIMME_INT ldevecsHat, HSCALAR *M, int ldM, HSCALAR *Mfact, int ldMfact,
      int *ipivot, HSCALAR *VtBV, int ldVtBV, HSCALAR *fVtBV, int ldfVtBV,
      int maxRank, int *basisSize, int *nextGuess, int *numGuesses,
      primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i;
   int initSize;
   int random = 0;

   /*-----------------------------------------------------------------------*/
   /* Orthogonalize the orthogonalization constraints provided by the user. */
   /* If a preconditioner is given and inner iterations are to be           */
   /* performed, then initialize M.                                         */
   /*-----------------------------------------------------------------------*/

   if (primme->numOrthoConst > 0) {
      int nV;
      CHKERR(Bortho_block_Sprimme(evecs, ldevecs, VtBV, ldVtBV, fVtBV, ldfVtBV,
            NULL, 0, 0, primme->numOrthoConst - 1, NULL, 0, 0, Bevecs, ldBevecs,
            NULL, 0, nLocal, maxRank, &nV, ctx));
      CHKERRM(nV != primme->numOrthoConst, PRIMME_ORTHO_CONST_FAILURE,
            "The given orthogonal constrains are not full rank");

      /* Initialize evecsHat, M, and its factorization Mfact,ipivot. This */
      /* allows the orthogonalization constraints to be included in the   */
      /* projector (I-BQQ'). Only needed if there is preconditioning, and */
      /* JDqmr inner iterations with a right, skew projector. Only in     */
      /* that case, M and Mfact are not NULL                              */

      if (M != NULL) {

         primme->ShiftsForPreconditioner = NULL;

         CHKERR(applyPreconditioner_Sprimme(Bevecs ? Bevecs : evecs,
               primme->nLocal, Bevecs ? ldBevecs : ldevecs, evecsHat,
               ldevecsHat, primme->numOrthoConst, ctx));

         CHKERR(update_XKinvBX_Sprimme(evecs, ldevecs, evecsHat, ldevecsHat, M,
               ldM, 0, primme->numOrthoConst, Mfact, ldMfact, ipivot, ctx));

      } /* if evecsHat and M=evecs'evecsHat, UDU are needed */

   } /* if numOrthoCont >0 */

   /* Handle case when some or all initial guesses are provided by */
   /* the user                                                     */
   if (!primme->locking) {
      initSize = min(primme->maxBasisSize, primme->initSize);
   } else {
      initSize = min(primme->minRestartSize, primme->initSize);
   }
   initSize = max(0, min(primme->n - primme->numOrthoConst, initSize));
   *numGuesses = primme->initSize - initSize;
   *nextGuess = primme->numOrthoConst + initSize;

   /* Copy over the initial guesses provided by the user */
   CHKERR(Num_copy_matrix_Sprimme(&evecs[primme->numOrthoConst * ldevecs],
         nLocal, initSize, ldevecs, V, ldV, ctx));

   switch (primme->initBasisMode) {
   case primme_init_krylov: random = 0; break;
   case primme_init_random:
      random = max(0, primme->minRestartSize - initSize);
      break;
   case primme_init_user:
      random = max(primme->maxBlockSize - initSize, 0);
      break;
   default: assert(0);
   }
   random = max(0, min(primme->n - primme->numOrthoConst - initSize, random));
   for (i = 0; i < random; i++) {
      Num_larnv_Sprimme(
            2, primme->iseed, nLocal, &V[ldV * (initSize + i)], ctx);
   }
   *basisSize = initSize + random;

   /* Orthonormalize the guesses provided by the user */
   CHKERR(Bortho_block_Sprimme(V, ldV, VtBV, ldVtBV, fVtBV, ldfVtBV, NULL, 0, 0,
         *basisSize - 1, evecs, ldevecs, primme->numOrthoConst, BV, ldBV, NULL,
         0, nLocal, maxRank, basisSize, ctx));

   CHKERR(matrixMatvec_Sprimme(V, nLocal, ldV, W, ldW, 0, *basisSize, ctx));

   if (primme->initBasisMode == primme_init_krylov) {
      int minRestartSize =
            min(primme->minRestartSize, primme->n - primme->numOrthoConst);
      CHKERR(init_block_krylov(V, nLocal, ldV, W, ldW, BV, ldBV, *basisSize,
            minRestartSize - 1, evecs, ldevecs, primme->numOrthoConst, VtBV,
            ldVtBV, fVtBV, ldfVtBV, maxRank, ctx));

      *basisSize = minRestartSize;
   }

   return 0;
}


/*******************************************************************************
 * Subroutine init_block_krylov - Initializes the basis as an orthonormal 
 *    block Krylov subspace.  
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * dv1, dv2    Range of indices over which the basis will be generated
 * 
 * locked      The array of locked Ritz vectors
 * 
 * numLocked   The number of vectors in the locked array
 *
 * machEps     machine precision needed in ortho()
 *
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * V  The orthonormal basis
 * 
 * W  A*V
 *
 * VtBV         V'*B*V (used by Bortho_block)
 *
 * fVtBV        The Cholesky factor of VtBV (used by Bortho_block)
 *
 * Return value
 * ------------
 * int -  0 upon success
 *       -1 if orthogonalization failed
 * 
 ******************************************************************************/

STATIC int init_block_krylov(SCALAR *V, PRIMME_INT nLocal, PRIMME_INT ldV,
      SCALAR *W, PRIMME_INT ldW, SCALAR *BV, PRIMME_INT ldBV, int dv1, int dv2,
      SCALAR *locked, PRIMME_INT ldlocked, int numLocked, HSCALAR *VtBV,
      int ldVtBV, HSCALAR *fVtBV, int ldfVtBV, int maxRank, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i;               /* Loop variables */
   int numNewVectors;   /* Number of vectors to be generated */
   int blockSize;       /* blockSize used in practice */
   
   numNewVectors = dv2 - dv1 + 1;

   /* Quick exit */

   if (numNewVectors <= 0) return 0;
 
   /*----------------------------------------------------------------------*/
   /* Generate a single Krylov space if there are only a few vectors to be */
   /* generated, else generate a block Krylov space with                   */
   /* primme->maxBlockSize as the block Size.                              */ 
   /*----------------------------------------------------------------------*/

   blockSize = numNewVectors <= primme->maxBlockSize ? 1 : primme->maxBlockSize;

   /*----------------------------------------------------------------------*/
   /* Generate the initial vectors.                                        */
   /*----------------------------------------------------------------------*/

   for (i=dv1; i<dv1+blockSize; i++) {
      CHKERR(Num_larnv_Sprimme(2, primme->iseed, nLocal, &V[ldV*i], ctx));
   }
   int nV=0;
   CHKERR(Bortho_block_Sprimme(V, ldV, VtBV, ldVtBV, fVtBV, ldfVtBV, NULL, 0,
         dv1, dv1 + blockSize - 1, locked, ldlocked, numLocked, BV, ldBV, NULL,
         0, nLocal, maxRank, &nV, ctx));
   CHKERRM(nV != dv1+blockSize, -1, "Random basis is not full rank");

   /* Generate the remaining vectors in the sequence */

   int m = blockSize;
   for (i = dv1 + blockSize, m = min(m, dv2 - i + 1); i <= dv2;
         i += m, m = min(m, dv2 - i + 1)) {
      CHKERR(matrixMatvec_Sprimme(&V[ldV*(i-blockSize)], nLocal, ldV,
               &V[ldV*i], ldV, 0, m, ctx));

      CHKERR(Num_copy_matrix_Sprimme(
            &V[ldV * i], nLocal, m, ldV, &W[ldW * (i - blockSize)], ldW, ctx));

      CHKERR(Bortho_block_Sprimme(V, ldV, VtBV, ldVtBV, fVtBV, ldfVtBV, NULL, 0,
            i, i + m - 1, locked, ldlocked, numLocked, BV, ldBV, NULL, 0,
            nLocal, maxRank, &nV, ctx));
      int j;
      for (j = nV; j < i + m; j++) {
         Num_larnv_Sprimme(2, primme->iseed, nLocal, &V[ldV * j], ctx);
      }
      CHKERR(Bortho_block_Sprimme(V, ldV, VtBV, ldVtBV, fVtBV, ldfVtBV, NULL, 0,
            nV, i + m - 1, locked, ldlocked, numLocked, BV, ldBV, NULL, 0,
            nLocal, maxRank, &nV, ctx));
      CHKERRM(nV != i+m, -1, "Random basis in not full rank");
   }

   CHKERR(matrixMatvec_Sprimme(V, nLocal, ldV, W, ldW, dv2-blockSize+1,
            blockSize, ctx));

   return 0;
}

#endif /* SUPPORTED_TYPE */
