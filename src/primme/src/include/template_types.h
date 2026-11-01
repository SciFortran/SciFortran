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
 * File: template_types.h
 *
 * Purpose - Force a compilation for every type: half, float, double, quad
 *           and the complex variants.
 *
 ******************************************************************************/

#ifndef TEMPLATE_TYPES
#define TEMPLATE_TYPES

#ifdef THIS_FILE

/* Host types */

// #define SHOW_TYPE

#include "template_undef.h"

#ifdef SHOW_TYPE
#warning compiling double
#endif
#define USE_DOUBLE
#include THIS_FILE
#include "template_undef.h"
#undef USE_DOUBLE

#ifdef SHOW_TYPE
#warning compiling half
#endif
#define USE_HALF
#include THIS_FILE
#include "template_undef.h"
#undef USE_HALF

#ifdef SHOW_TYPE
#warning compiling half complex
#endif
#define USE_HALFCOMPLEX
#include THIS_FILE
#include "template_undef.h"
#undef USE_HALFCOMPLEX

#ifdef SHOW_TYPE
#warning compiling float
#endif
#define USE_FLOAT
#include THIS_FILE
#include "template_undef.h"
#undef USE_FLOAT

#ifdef SHOW_TYPE
#warning compiling float complex
#endif
#define USE_FLOATCOMPLEX
#include THIS_FILE
#include "template_undef.h"
#undef USE_FLOATCOMPLEX

#ifdef SHOW_TYPE
#warning compiling double complex
#endif
#define USE_DOUBLECOMPLEX
#include THIS_FILE
#include "template_undef.h"
#undef USE_DOUBLECOMPLEX

// #define USE_QUAD
// #include THIS_FILE
// #include "template_undef.h"
// #undef USE_QUAD
// #define USE_QUADCOMPLEX
// #include THIS_FILE
// #include "template_undef.h"
// #undef USE_QUADCOMPLEX

/* MAGMA types */


#ifdef SHOW_TYPE
#warning compiling half magma
#endif
#define USE_HALF_MAGMA
#include THIS_FILE
#include "template_undef.h"
#undef USE_HALF_MAGMA

#ifdef SHOW_TYPE
#warning compiling half complex magma
#endif
#define USE_HALFCOMPLEX_MAGMA
#include THIS_FILE
#include "template_undef.h"
#undef USE_HALFCOMPLEX_MAGMA

#ifdef SHOW_TYPE
#warning compiling float magma
#endif
#define USE_FLOAT_MAGMA
#include THIS_FILE
#include "template_undef.h"
#undef USE_FLOAT_MAGMA

#ifdef SHOW_TYPE
#warning compiling float complex magma
#endif
#define USE_FLOATCOMPLEX_MAGMA
#include THIS_FILE
#include "template_undef.h"
#undef USE_FLOATCOMPLEX_MAGMA

#ifdef SHOW_TYPE
#warning compiling double magma
#endif
#define USE_DOUBLE_MAGMA
#include THIS_FILE
#include "template_undef.h"
#undef USE_DOUBLE_MAGMA

#ifdef SHOW_TYPE
#warning compiling double complex magma
#endif
#define USE_DOUBLECOMPLEX_MAGMA
#include THIS_FILE
#include "template_undef.h"
#undef USE_DOUBLECOMPLEX_MAGMA


/* CUBLAS types */


#ifdef SHOW_TYPE
#warning compiling half cublas
#endif
#define USE_HALF_CUBLAS
#include THIS_FILE
#include "template_undef.h"
#undef USE_HALF_CUBLAS

#ifdef SHOW_TYPE
#warning compiling half complex cublas
#endif
#define USE_HALFCOMPLEX_CUBLAS
#include THIS_FILE
#include "template_undef.h"
#undef USE_HALFCOMPLEX_CUBLAS

#ifdef SHOW_TYPE
#warning compiling float cublas
#endif
#define USE_FLOAT_CUBLAS
#include THIS_FILE
#include "template_undef.h"
#undef USE_FLOAT_CUBLAS

#ifdef SHOW_TYPE
#warning compiling float complex cublas
#endif
#define USE_FLOATCOMPLEX_CUBLAS
#include THIS_FILE
#include "template_undef.h"
#undef USE_FLOATCOMPLEX_CUBLAS

#ifdef SHOW_TYPE
#warning compiling double cublas
#endif
#define USE_DOUBLE_CUBLAS
#include THIS_FILE
#include "template_undef.h"
#undef USE_DOUBLE_CUBLAS

#ifdef SHOW_TYPE
#warning compiling double complex cublas
#endif
#define USE_DOUBLECOMPLEX_CUBLAS
#include "template.h"
// #include THIS_FILE
// #include "template_undef.h"
// #undef USE_DOUBLECOMPLEX_CUBLAS
// #define USE_QUAD_CUBLAS
// #include THIS_FILE
// #include "template_undef.h"
// #undef USE_QUAD_CUBLAS
// #define USE_QUADCOMPLEX_CUBLAS
// #include THIS_FILE
// #undef USE_QUADCOMPLEX_CUBLAS
// #include "template_undef.h"

#endif /* THIS_FILE */
#endif /* TEMPLATE_TYPES */
