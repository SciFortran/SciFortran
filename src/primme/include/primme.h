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
 * File: primme.h
 * 
 * Purpose - Main header with the PRIMME C interface functions.
 * 
 ******************************************************************************/

#ifndef PRIMME_H
#define PRIMME_H

/* Define version */

#define PRIMME_VERSION_MAJOR      3
#define PRIMME_VERSION_MINOR      2

#if defined(__clang__) && defined(__FLT16_EPSILON__)
#  define PRIMME_HALF __fp16
#  define PRIMME_WITH_NATIVE_HALF
#else
   struct _primme_half {int short a;};
#  define PRIMME_HALF struct _primme_half
#endif

#define PRIMME_QUAD double long

/* A C99 code with complex type is not a valid C++ code. However C++          */
/* compilers usually can take it. Nevertheless in order to avoid the warnings */
/* while compiling in pedantic mode, we use the proper complex type for C99   */
/* (complex double and complex float) and C++ (std::complex<double> and       */
/* std::complex<float>). Of course both complex types are binary compatible.  */

#ifdef __cplusplus
#  include <complex>
   struct _primme_complex_half {PRIMME_HALF r; PRIMME_HALF i;};
#  define PRIMME_COMPLEX_HALF struct _primme_complex_half
#  define PRIMME_COMPLEX_FLOAT std::complex<float>
#  define PRIMME_COMPLEX_DOUBLE std::complex<double>
#  define PRIMME_COMPLEX_QUAD std::complex<PRIMME_QUAD>
#else
#  include <complex.h>
/* Complex half is not part of C99 */
   struct _primme_complex_half {PRIMME_HALF r; PRIMME_HALF i;};
#  define PRIMME_COMPLEX_HALF struct _primme_complex_half
#  define PRIMME_COMPLEX_FLOAT float complex
#  define PRIMME_COMPLEX_DOUBLE double complex
#  define PRIMME_COMPLEX_QUAD long double complex
#endif

/* Required by some C++ compilers when including inttypes.h */
#if defined(__cplusplus) && !defined(__STDC_FORMAT_MACROS)
#  define __STDC_FORMAT_MACROS
#endif

#if !defined(PRIMME_INT_SIZE) || PRIMME_INT_SIZE == 64
#  include <stdint.h>
#  include <inttypes.h>
#  define PRIMME_INT int64_t
#  define PRIMME_INT_P PRId64
#  define PRIMME_INT_MAX INT64_MAX
#  ifndef PRId64
#     define PRId64 "ld"
#  endif
#elif PRIMME_INT_SIZE == 0
#  include <limits.h>
#  define PRIMME_INT int
#  define PRIMME_INT_P "d"
#  define PRIMME_INT_MAX INT_MAX
#elif PRIMME_INT_SIZE == 32
#  include <stdint.h>
#  include <inttypes.h>
#  define PRIMME_INT int32_t
#  define PRIMME_INT_P PRId32
#  define PRIMME_INT_MAX INT32_MAX
#  ifndef PRId64
#     define PRId64 "d"
#  endif
#else
#  define PRIMME_INT PRIMME_INT_SIZE
#  define PRIMME_INT_P "d"
#  define PRIMME_INT_MAX INT_MAX
#endif

#include "primme_eigs.h"
#include "primme_svds.h"

/* Error messages */

#define PRIMME_UNEXPECTED_FAILURE         -1
#define PRIMME_MALLOC_FAILURE             -2
#define PRIMME_MAIN_ITER_FAILURE          -3
#define PRIMME_LAPACK_FAILURE             -40
#define PRIMME_USER_FAILURE               -41
#define PRIMME_ORTHO_CONST_FAILURE        -42
#define PRIMME_PARALLEL_FAILURE           -43
#define PRIMME_FUNCTION_UNAVAILABLE       -44
#endif /* PRIMME_H */
