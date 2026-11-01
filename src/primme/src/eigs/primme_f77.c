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
 * File: primme_f77.c
 *
 * Purpose - Implementation of PRIMME F77 interface functions.
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/primme_f77.c"
#endif

#include <stdlib.h>   /* free */
#include "numerical.h"
#include "template_normal.h"
#include "common_eigs.h"
#include "primme_interface.h" /* for Xprimme */


/******************************************************************************/
/* The following differ only in the name mangling from C                      */
/******************************************************************************/

/*****************************************************************************
 * Wrapper for calling xprimme from Fortran 77/90. 
 * The only difference from primme: the return value passed as parameter 
 *****************************************************************************/

EXTERN_C void AS_FORTRAN(Xprimme)(XEVAL *evals, XSCALAR *evecs, XREAL *rnorms,
      primme_params **primme, int *ierr) {

   *ierr = Xprimme(evals, evecs, rnorms, *primme);

} /* end of xprimme_f77 wrapper for calling from Fortran */


/* Only define these functions ones */
#ifdef USE_DOUBLE

/*****************************************************************************
 * Initialize handles also the allocation of primme structure 
 *****************************************************************************/
EXTERN_C void AS_FORTRAN(primme_initialize)(primme_params **primme) {

   *primme = NULL;
   if (MALLOC_PRIMME(1, primme) == 0)
      primme_initialize(*primme);
}

/*****************************************************************************
 *  * Free the internally allocated work arrays of the primme structure 
 *   *****************************************************************************/
EXTERN_C void AS_FORTRAN(primme_free)(primme_params **primme) {

   primme_free(*primme);
   free(*primme);
   *primme = NULL;
}

/*****************************************************************************
 * Wrapper for displaying the primme parameters 
 *****************************************************************************/
EXTERN_C void AS_FORTRAN(primme_display_params)(primme_params **primme) {

   primme_display_params(*(*primme));
}

/*************************************************************************
 * set method takes "int method" as input. On return:
 *    returnValue = 0 successful completion 
 *    returnValue < 0 no such method exists. If not defined by user, defaults 
 *    have been set for maxBasisSize, minRestartSize, and maxBlockSize   
 *************************************************************************/
EXTERN_C void AS_FORTRAN(primme_set_method)(primme_params **primme,
      primme_preset_method *method, int *returnValue) {
  *returnValue = primme_set_method(*method, *primme);
}

/*************************************************************************
 * primme_set_member_f77(primme_params *primme, int *label, void *ptr)
 *
 * Sets any of the members of the primme data structure from a Fortran call.
 *
 * primme is a pointer to the primme data structure.
 *
 * void *ptr generic type pointer is used for input, and is cast to the 
 * required pointer type before it is assigned to a data structure member.
 *
 * The choice of member is through an integer label. The integers describing
 * the member are prefixed with PRIMMEF77_ followed by the actual name of 
 * the member. If the member is a structure, underscores are used for dots,
 *
 * Eg the label int for:    primme->correctionParams.projectors.LeftX
 * is:                    PRIMMEF77_correctionParams_projectors_LeftX
 *
 *************************************************************************/

EXTERN_C void AS_FORTRAN(primme_set_member)(primme_params **primme, int *label,
      void *v, int *ierr) {
   *ierr = primme_set_member(*primme, (primme_params_label)*label, v);
}

/*************************************************************************
 * primme_get_member_f77(primme_params *primme, int *label, void *ptr)
 *
 * Gets the value of any of the members of the primme data structure from 
 * a Fortran call.
 *
 * To be called ONLY from a f77 subroutine that is called by C PRIMME, 
 * (eg.,MV, PRECOND, GLOBALSUM) ---NOT from F77 driver.
 *
 *   F77 driver --> C PRIMME --> F77 function (MV) --> primme_get_member_f77
 *
 * For obtaining information from the driver, see below. The reason is that 
 * such f77 functions (MV, etc) receive the primme structure from PRIMME as
 * (primme_params *primme). However, the f77 driver has a pointer to **primme. 
 *
 * void *ptr generic type pointer is used for output, and is cast to the 
 * required pointer type before it is assigned the data structure member.
 *
 * The choice of member is through an integer label. The integers describing
 * the member are prefixed with PRIMMEF77_ followed by the actual name of 
 * the member. If the member is a structure, underscores are used for dots,
 *
 * Eg the label int for:    primme->correctionParams.projectors.LeftX
 * is:                    PRIMMEF77_correctionParams_projectors_LeftX
 *
 *************************************************************************/

EXTERN_C void AS_FORTRAN(primme_get_member)(primme_params *primme, int *label,
      void *v, int *ierr) {
   *ierr = primme_get_member(primme, (primme_params_label)*label, v);
}

/*************************************************************************
 *   If the user requires the current shift suggested by PRIMME to 
 *   be used in inverting (M-shift_i I), for the i-th eigenvalue
 *   this is at location (i-1) in C   (i-th location in Fortran)
 *   of the array primme->ShiftsForPreconditioner
 *************************************************************************/
EXTERN_C void AS_FORTRAN(primme_get_prec_shift)(primme_params *primme, int *i,
      double *shift) {
   *shift = primme->ShiftsForPreconditioner[*i-1];
}

/*************************************************************************/
/* If the information from primme is needed from the calling driver in F77
 * the primme is passed as  (primme_params **), while any f77 function such
 * as MV or preconditioning that needs primme information it passes it as
 * (primme_params *). The reason is that PRIMME C code works with the pointer
 * to primme and that's what it calls MV/GLOBALSUM with.
 *
 *   F77 driver --> other F77 subs --> primmetop_get_member_f77
 * 
 * If the user needs to obtain structure info from the top function that 
 * initialized primme they should call the wrappers:
 *
 *      primmetop_get_member_f77
 *  or
 *      primmetop_get_prec_shift_f77   (probably not needed ever)
 *
 * From within MV, PRECONDITION, or GLOBALSUMDOUBLE call the above function
 *      primme_get_member_f77
 *
 *************************************************************************/
EXTERN_C void AS_FORTRAN(primmetop_get_member)(primme_params **primme, int *label,
      void *ptr, int *ierr) {
   *ierr = primme_get_member(*primme, (primme_params_label)*label, ptr);
}

EXTERN_C void AS_FORTRAN(primmetop_get_prec_shift)(primme_params **primme, int *i,
      double *shift) {
   AS_FORTRAN(primme_get_prec_shift)(*primme, i, shift);
}

#endif /* USE_DOUBLE */
