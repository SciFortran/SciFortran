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
 * File: template.h
 *
 * Purpose - Contains definitions of macros used along PRIMME.
 *    In short source files, *.c, are compiled several times, every time for a
 *    different type, referred as SCALAR. Examples of types are float, double,
 *    complex float, complex double, and the corresponding GPU versions. All
 *    types are described in section "Arithmetic". Other macros are defined to
 *    refer derived types and functions. An example is REAL, which is
 *    defined as the real (non-complex) version of SCALAR. For instance,
 *    when SCALAR is complex double, REAL is double. Similarly, it is possible
 *    to call the real version of a function. For instance,
 *    permute_vecs_Sprimme permutes vectors with type SCALAR and
 *    permute_vecs_Rprimme permutes vectors with type REAL.
 *
 *    When SCALAR is a GPU type, the pointers SCALAR* are supposed to point out
 *    memory allocated on GPUs, also called devices. For instance,
 *    Num_malloc_Sprimme allocates GPU memory when SCALAR is a GPU type. Use
 *    HSCALAR and HREAL as the non-GPU, also called host, versions of SCALAR and
 *    REAL. Also to use the non-GPU version use the suffices _SHprimme and
 *    _RHprimme. For instance,
 *    Num_malloc_SHprimme allocates memory on the host, and
 *    permute_vecs_RHprimme permute REAL vectors on the host.
 *
 ******************************************************************************/

#ifndef TEMPLATE_H
#define TEMPLATE_H

#include "template_types.h"

/**********************************************************************
 * Macros USE_FLOAT, USE_FLOATCOMPLEX, USE_DOUBLE and USE_DOUBLECOMPLEX -
 *    only one of them is defined at the same time, and identifies the
 *    type of SCALAR, one of float, complex float, double or complex double.
 *
 * Macro USE_COMPLEX - only defined when SCALAR is a complex type.
 *
 * Macro USE_HOST - CPU version
 *
 * Macro USE_MAGMA - MAGMA version
 *
 * Macro SUPPORTED_TYPE - defined if functions with the current type
 *    are going to be built.
 *
 * Macro SUPPORTED_HALF_TYPE - defined if functions with the current type
 *    has a version in half.
 **********************************************************************/

/* Helper macros and types used to define SCALAR and REAL and their variants */

#define HOST_STEM

#if defined(USE_HALF) || defined(USE_HALFCOMPLEX) || defined(USE_FLOAT) ||     \
      defined(USE_FLOATCOMPLEX) || defined(USE_DOUBLE) ||                      \
      defined(USE_DOUBLECOMPLEX)
#  define USE_HOST
#  define STEM
#  define IMPL(BL, MA) BL
#elif defined(USE_FLOAT_MAGMA) || defined(USE_FLOATCOMPLEX_MAGMA) ||           \
      defined(USE_DOUBLE_MAGMA) || defined(USE_DOUBLECOMPLEX_MAGMA) ||         \
      defined(USE_HALF_MAGMA) || defined(USE_HALFCOMPLEX_MAGMA)
#  define USE_MAGMA
#  define STEM magma_
#  define IMPL(BL, MA) MA
#elif defined(USE_FLOAT_CUBLAS) || defined(USE_FLOATCOMPLEX_CUBLAS) ||         \
      defined(USE_DOUBLE_CUBLAS) || defined(USE_DOUBLECOMPLEX_CUBLAS) ||       \
      defined(USE_HALF_CUBLAS) || defined(USE_HALFCOMPLEX_CUBLAS)
#  define USE_CUBLAS
#  define STEM cublas_
#  define IMPL(BL, MA) MA
#else
#  error 
#endif

#if !defined(CHECK_TEMPLATE) && !defined(STEM_C)
#   define STEM_C STEM
#endif

#if defined(USE_HALF) || defined(USE_HALF_MAGMA) ||                            \
      defined(USE_HALF_CUBLAS) || defined(USE_FLOAT) ||                        \
      defined(USE_FLOAT_MAGMA) || defined(USE_FLOAT_CUBLAS) ||                 \
      defined(USE_DOUBLE) || defined(USE_DOUBLE_MAGMA) ||                      \
      defined(USE_DOUBLE_CUBLAS) || defined(USE_QUAD) ||                       \
      defined(USE_QUAD_MAGMA) || defined(USE_QUAD_CUBLAS)
#  define USE_REAL
#elif defined(USE_HALFCOMPLEX) || defined(USE_HALFCOMPLEX_MAGMA) ||            \
      defined(USE_HALFCOMPLEX_CUBLAS) || defined(USE_FLOATCOMPLEX) ||          \
      defined(USE_FLOATCOMPLEX_MAGMA) || defined(USE_FLOATCOMPLEX_CUBLAS) ||   \
      defined(USE_DOUBLECOMPLEX) || defined(USE_DOUBLECOMPLEX_MAGMA) ||        \
      defined(USE_DOUBLECOMPLEX_CUBLAS) || defined(USE_QUADCOMPLEX) ||         \
      defined(USE_QUADCOMPLEX_MAGMA) || defined(USE_QUADCOMPLEX_CUBLAS)
#  define USE_COMPLEX
#else
#  error 
#endif

#if   defined(USE_DOUBLE)        || defined(USE_DOUBLE_MAGMA)          || defined(USE_DOUBLE_CUBLAS)
#  define      ARITH(H,K,S,C,D,Z,Q,W) D
#  define REAL_ARITH(H,K,S,C,D,Z,Q,W) D
#elif defined(USE_DOUBLECOMPLEX) || defined(USE_DOUBLECOMPLEX_MAGMA)   || defined(USE_DOUBLECOMPLEX_CUBLAS)
#  define      ARITH(H,K,S,C,D,Z,Q,W) Z
#  define REAL_ARITH(H,K,S,C,D,Z,Q,W) D
#elif defined(USE_FLOAT)         || defined(USE_FLOAT_MAGMA)           || defined(USE_FLOAT_CUBLAS)
#  define      ARITH(H,K,S,C,D,Z,Q,W) S
#  define REAL_ARITH(H,K,S,C,D,Z,Q,W) S
#elif defined(USE_FLOATCOMPLEX)  || defined(USE_FLOATCOMPLEX_MAGMA)    || defined(USE_FLOATCOMPLEX_CUBLAS)      
#  define      ARITH(H,K,S,C,D,Z,Q,W) C
#  define REAL_ARITH(H,K,S,C,D,Z,Q,W) S
#elif defined(USE_HALF)          || defined(USE_HALF_MAGMA)            || defined(USE_HALF_CUBLAS) 
#  define      ARITH(H,K,S,C,D,Z,Q,W) H
#  define REAL_ARITH(H,K,S,C,D,Z,Q,W) H
#elif defined(USE_HALFCOMPLEX)   || defined(USE_HALFCOMPLEX_MAGMA)     || defined(USE_HALFCOMPLEX_CUBLAS)          
#  define      ARITH(H,K,S,C,D,Z,Q,W) K
#  define REAL_ARITH(H,K,S,C,D,Z,Q,W) H
#elif defined(USE_QUAD)          || defined(USE_QUAD_MAGMA)            || defined(USE_QUAD_CUBLAS) 
#  define      ARITH(H,K,S,C,D,Z,Q,W) Q
#  define REAL_ARITH(H,K,S,C,D,Z,Q,W) Q
#elif defined(USE_QUADCOMPLEX)   || defined(USE_QUADCOMPLEX_MAGMA)     || defined(USE_QUADCOMPLEX_CUBLAS)          
#  define      ARITH(H,K,S,C,D,Z,Q,W) W
#  define REAL_ARITH(H,K,S,C,D,Z,Q,W) Q
#else
#  error
#endif

/* For host types, define SUPPORTED_HALF_TYPE when the compiler supports half
 * precision. For MAGMA, define the macro if MAGMA also supports half precision.
 *
 * Define SUPPORTED_TYPE when inspecting the signature functions to generate
 * the signature for all possible functions. Also define the macro for any
 * setting without half precision, and for half precision if the compiler
 * supports half precision.
 */

#if defined(PRIMME_WITH_HALF) && defined(PRIMME_WITH_NATIVE_HALF) &&           \
      (defined(USE_HOST) ||                                                    \
            (defined(PRIMME_WITH_MAGMA) && defined(USE_MAGMA) &&               \
                  defined(MAGMA_WITH_HALF)) ||                                 \
            ((defined(PRIMME_WITH_CUBLAS) || defined(PRIMME_WITH_HIPBLAS)) &&  \
                  defined(USE_CUBLAS)))
#  define SUPPORTED_HALF_TYPE
#endif

// Undefine SUPPORTED_TYPE when the current type is not supported. That is if
// one the next applies:
// - USE_HALF/COMPLEX/_MAGMA|CUBLAS is defined but SUPPORTED_HALF_TYPE is not.
// - USE_FLOAT/COMPLEX/_MAGMA|CUBLAS is defined but PRIMME_WITHOUT_FLOAT is defined.
// - USE_MAGMA is defined but PRIMME_WITH_MAGMA is not.

#define SUPPORTED_TYPE
#if !defined(CHECK_TEMPLATE) &&                                                \
      (((defined(USE_HALF) || defined(USE_HALFCOMPLEX) ||                      \
              defined(USE_HALF_MAGMA) || defined(USE_HALFCOMPLEX_MAGMA) ||     \
              defined(USE_HALF_CUBLAS) || defined(USE_HALFCOMPLEX_CUBLAS)) &&  \
             !defined(SUPPORTED_HALF_TYPE)) ||                                 \
            ((defined(USE_HALF_MAGMA) || defined(USE_HALFCOMPLEX_MAGMA)) &&    \
                  !defined(MAGMA_WITH_HALF)) ||                                \
            (defined(USE_MAGMA) && !defined(PRIMME_WITH_MAGMA)) ||             \
            (defined(USE_CUBLAS) && !(defined(PRIMME_WITH_CUBLAS) ||           \
                                          defined(PRIMME_WITH_HIPBLAS))) ||     \
                  ((defined(USE_FLOAT) || defined(USE_FLOATCOMPLEX) ||         \
                         defined(USE_FLOAT_MAGMA) ||                           \
                         defined(USE_FLOATCOMPLEX_MAGMA) ||                    \
                         defined(USE_HALF_CUBLAS) ||                           \
                         defined(USE_HALFCOMPLEX_CUBLAS)) &&                   \
                        defined(PRIMME_WITHOUT_FLOAT)))
#  undef SUPPORTED_TYPE
#endif

/* A C99 code with complex type is not a valid C++ code. However C++          */
/* compilers usually can take it. Nevertheless in order to avoid the warnings */
/* while compiling in pedantic mode, we use the proper complex type for C99   */
/* (complex double and complex float) and C++ (std::complex<double> and       */
/* std::complex<float>). Of course both complex types are binary compatible.  */

#ifdef USE_COMPLEX
#  ifndef __cplusplus
#     define REAL_PART(x) (creal(x))
#     define IMAGINARY_PART(x) (cimag(x))
#     define ABS(x) (cabs(x))
#     define CONJ(x) (conj(x))
#  else
#     define REAL_PART(x) (std::real(x))
#     define IMAGINARY_PART(x) (std::imag(x))
#     define ABS(x) (std::abs(x))
#     define CONJ(x) (std::conj(x))
#  endif
#else
#  define REAL_PART(x) (x)
#  define IMAGINARY_PART(x) 0
#  define ABS(x) (fabs(x))
#  define CONJ(x) (x)
#endif

/* Helper macros to support complex arithmetic for types without complex   */
/* support in C99. For now, only half precision has this problem. The      */
/* approach is to cast the unsupported complex type to a supported type    */
/* with more precision. For instance, complex half precision is cast to    */
/* complex single precision.                                               */
/*                                                                         */
/* NOTE: 'A' is an unsupported complex type and 'B' is a supported type    */
/* SET_ZERO(A)       : set A = 0                                           */
/* SET_COMPLEX(A, B) : set A = B                                           */
/* TO_COMPLEX(A)     : cast A to a supported complex type                  */
/* PLUS_EQUAL(A, B)  : set A += B                                          */
/* MULT_EQUAL(A, B)  : set A *= B                                          */

#if (defined(USE_HALFCOMPLEX) || defined(USE_HALFCOMPLEX_MAGMA) || defined(USE_HALFCOMPLEX_CUBLAS)) && !defined(PRIMME_WITH_NATIVE_COMPLEX_HALF)
#  define SET_ZERO(A) {(A).r = 0; (A).i = 0;}
#  define SET_COMPLEX(A,B) {(A).r = REAL_PART(B); (A).i = IMAGINARY_PART(B);}
#  ifndef __cplusplus
#     define TO_COMPLEX(A) ((A).r + (A).i * _Complex_I)
#  else
#     define TO_COMPLEX(A) (std::complex<HREAL>((HREAL)((A).r), (HREAL)((A).i)))
#  endif
#  define PLUS_EQUAL(A,B) {(A).r += REAL_PART(B); (A).i += IMAGINARY_PART(B);}
#  define MULT_EQUAL(A, B)                                                     \
   {                                                                           \
      HSCALAR C = TO_COMPLEX(A) * (B);                                         \
      (A).r += REAL_PART(C);                                                   \
      (A).i += IMAGINARY_PART(C);                                              \
   }
#else
#  define SET_ZERO(A) {(A) = 0.0;}
#  define SET_COMPLEX(A,B) (A) = (B)
#  if defined(USE_HALFCOMPLEX) && defined(__cplusplus)
#     define TO_COMPLEX(A) (HSCALAR(REAL_PART(A), IMAGINARY_PART(A)))
#     define PLUS_EQUAL(A, B) (A) = TO_COMPLEX(A) + (B)
#     define MULT_EQUAL(A, B) (A) = TO_COMPLEX(A) * (B)
#  else
#     define TO_COMPLEX(A) (A)
#     define PLUS_EQUAL(A, B) (A) += (B)
#     define MULT_EQUAL(A, B) (A) *= (B)
#  endif
#endif


/* TEMPLATE_PLEASE tags the functions whose prototypes depends on macros and  */
/* are used in other files. The macro has value only when the tool ctemplate  */
/* is inspecting the source files, which is indicated by the macro            */
/* CHECK_TEMPLATE being defined. See Makefile and tools/ctemplate.            */
/*                                                                            */
/* When SCALAR is not a complex type (e.g., float and double) the function    */
/* will be referred as _Sprimme and _Rprimme. Otherwise it will be referred   */
/* only as _Sprimme. The term TEMPLATE_PLEASE should prefix every function    */
/* that will be instantiated with different values for SCALAR and REAL.       */

#ifndef TEMPLATE_H_PRIVATE
#define TEMPLATE_H_PRIVATE

#  define USE_ARITH(Re,Co) ARITH(Re,Co,Re,Co,Re,Co,Re,Co)

#  define USE_SR(Re,Co,T,XH,STEM,POST) \
      USE(CONCAT(CONCAT(CONCAT(S,XH),T),primme), STR0(CONCAT(CONCAT(CONCAT(STEM,USE_ARITH(Re,Co)),primme),POST))) \
      USE(CONCAT(CONCAT(CONCAT(R,XH),T),primme), STR0(CONCAT(CONCAT(CONCAT(STEM,Re),primme),POST)))

#  define USE_TYPE(H,K,S,C,D,Z,Q,W,XH,STEM,POST)  \
      USE_SR(H,K,h,XH,STEM,POST) \
      USE_SR(S,C,s,XH,STEM,POST) \
      USE_SR(D,Z,d,XH,STEM,POST) \
      USE_SR(Q,W,q,XH,STEM,POST)

#endif /* TEMPLATE_H_PRIVATE */

#ifdef CHECK_TEMPLATE
#  define TEMPLATE_PLEASE \
      APPEND_FUNC(Sprimme,SCALAR_SUF) USE(Sprimme,"SCALAR_SUF") \
      USE(Rprimme,"REAL_SUF") USE(SHprimme,"HOST_SCALAR_SUF") \
      USE(RHprimme,"HOST_REAL_SUF") USE(SXprimme,"XSCALAR_SUF") \
      USE(RXprimme,"XREAL_SUF") USE_TYPE(h,k,s,c,d,z,q,w, , STEM_C, ) \
      USE_TYPE(h,k,s,c,d,z,q,w, X, HOST_STEM, ) \
      USE_TYPE(s,c,s,c,d,z,q,w, H, HOST_STEM, )

#  define STATIC APPEND_FUNC(,SCALAR_SUF) USE(,"SCALAR_SUF")

#else
#  define TEMPLATE_PLEASE
#  define STATIC
#endif /* CHECK_TEMPLATE */

/* Avoid to use the final type for integers and complex in generated       */
/* headers file. Instead use PRIMME_COMPLEX_FLOAT, _HALF and _DOUBLE.      */

#ifdef CHECK_TEMPLATE
#  undef PRIMME_HALF
#  undef PRIMME_COMPLEX_HALF
#  undef PRIMME_COMPLEX_FLOAT
#  undef PRIMME_COMPLEX_DOUBLE
#  undef PRIMME_QUAD
#  undef PRIMME_COMPLEX_QUAD
#  undef PRIMME_INT
#endif /* CHECK_TEMPLATE */

#endif /* TEMPLATE_H */
