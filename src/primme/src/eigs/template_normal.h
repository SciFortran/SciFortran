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
 * File: template_normal.h
 *
 * Purpose - Force a compilation for Hermitian and normal operator
 *
 ******************************************************************************/

#ifndef TEMPLATE_NORMAL_H
#define TEMPLATE_NORMAL_H

#ifndef WITH_KIND
#define WITH_KIND(X) CONCAT(X,KIND(,_normal))
#endif

#ifdef THIS_FILE

#ifdef CHECK_TEMPLATE
#  undef TEMPLATE_PLEASE
#  undef STATIC
#  define TEMPLATE_PLEASE \
      APPEND_FUNC(Sprimme,WITH_KIND(SCALAR_SUF)) \
      USE(Sprimme, STR0(WITH_KIND(SCALAR_SUF))) \
      USE(Rprimme, STR0(WITH_KIND(REAL_SUF))) \
      USE(SHprimme,STR0(WITH_KIND(HOST_SCALAR_SUF))) \
      USE(RHprimme,STR0(WITH_KIND(HOST_REAL_SUF))) \
      USE(SXprimme,STR0(WITH_KIND(XSCALAR_SUF))) \
      USE(RXprimme,STR0(WITH_KIND(XREAL_SUF))) \
      USE_TYPE(h,k,s,c,d,z,q,w,  , STEM_C, KIND_C) \
      USE_TYPE(h,k,s,c,d,z,q,w, X, HOST_STEM, KIND_C) \
      USE_TYPE(s,c,s,c,d,z,q,w, H, HOST_STEM, KIND_C)

#  define STATIC APPEND_FUNC(,WITH_KIND(SCALAR_SUF)) USE(,STR0(WITH_KIND(SCALAR_SUF)))
#elif !defined(KIND_C)
#  define KIND_C WITH_KIND()
#endif

#undef KIND
#undef USE_HERMITIAN
#undef USE_NORMAL

// #define SHOW_TYPE

#ifdef USE_COMPLEX
#  ifdef SHOW_TYPE
#     warning compiling normal
#  endif
#  define USE_NORMAL
#  define KIND(H,N) N
#  include THIS_FILE
#  undef USE_NORMAL
#  undef KIND
#endif

#ifdef SHOW_TYPE
#warning compiling Hermitian
#endif
#define USE_HERMITIAN
#define KIND(H,N) H
// #include THIS_FILE
// #undef USE_HERMITIAN
// #undef KIND

#undef TEMPLATE_NORMAL_H

#endif /* THIS_FILE */
#endif /* TEMPLATE_NORMAL_H */
