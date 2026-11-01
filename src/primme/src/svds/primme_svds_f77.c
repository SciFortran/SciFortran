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
 * File: primme_svds_f77.c
 *
 * Purpose - Implementation of PRIMME_SVDS F77 interface functions.
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../svds/primme_svds_f77.c"
#endif


#include <stdlib.h>   /* free */
#include "numerical.h"
#include "primme_svds_interface.h"


EXTERN_C void AS_FORTRAN(Sprimme_svds)(XREAL *svals, XSCALAR *svecs,
      XREAL *resNorms, primme_svds_params **primme_svds, int *ierr) {
   *ierr = Sprimme_svds(svals, svecs, resNorms, *primme_svds);
}

/* Only define these functions ones */
#ifdef USE_DOUBLE

EXTERN_C void AS_FORTRAN(primme_svds_initialize)(primme_svds_params **primme_svds) {
   *primme_svds = NULL;
   if (MALLOC_PRIMME(1, primme_svds) == 0)
      primme_svds_initialize(*primme_svds);
}

EXTERN_C void AS_FORTRAN(primme_svds_set_method)(primme_svds_preset_method *method,
      primme_preset_method *methodStage1, primme_preset_method *methodStage2,
      primme_svds_params **primme_svds, int *ierr) {
   *ierr = primme_svds_set_method(*method, *methodStage1, *methodStage2,
         *primme_svds);
}

EXTERN_C void AS_FORTRAN(primme_svds_display_params)(primme_svds_params **primme_svds) {
   primme_svds_display_params(**primme_svds);
}

EXTERN_C void AS_FORTRAN(primme_svds_free)(primme_svds_params **primme_svds) {
   primme_svds_free(*primme_svds);
   free(*primme_svds);
   *primme_svds = NULL;
}

EXTERN_C void AS_FORTRAN(primme_svds_set_member)(primme_svds_params **primme_svds,
      int *label, void *v, int *ierr) {
   *ierr = primme_svds_set_member(*primme_svds,
         (primme_svds_params_label)*label, v);
}

EXTERN_C void AS_FORTRAN(primme_svdstop_get_member)(
      primme_svds_params **primme_svds, int *label, void *v, int *ierr) {
   *ierr = primme_svds_get_member(*primme_svds,
         (primme_svds_params_label)*label, v);
}

EXTERN_C void AS_FORTRAN(primme_svds_get_member)(primme_svds_params *primme_svds,
      int *label, void *v, int *ierr) {
   *ierr = primme_svds_get_member(primme_svds,
         (primme_svds_params_label)*label, v);
}

#endif /* USE_DOUBLE */
