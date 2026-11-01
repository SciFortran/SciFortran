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
 * File: wtime.c
 *
 * Purpose - Time functions.
 *
 ******************************************************************************/

#define THIS_FILE "../linalg/wtime.c"

#include <stdlib.h>
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#  include <sys/time.h>
#  include <sys/resource.h>
#endif

#ifndef CHECK_TEMPLATE
#include "wtime.h"
#endif

#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))

double primme_wTimer() {
   struct timeval tv;
   gettimeofday(&tv, NULL);
   return ((double)tv.tv_sec) + ((double)tv.tv_usec) / (double)1E6;
}

/* In the unlikely event that gettimeofday() is not available, but POSIX is, 
 * we can use the following alternative definition for primme_wTimer, 
 * after including time.h at the top.
 */
/*
#include <time.h>
double primme_wTimer(int zeroTimer) {
   static struct timespec ts;
   static double StartingTime;
   
   if (zeroTimer) {
      clock_gettime(CLOCK_REALTIME, &ts);
      StartingTime = ((double) ts.tv_sec) + ((double) ts.tv_nsec )/(double) 1E9;
      return StartingTime;
   }
   else {
      clock_gettime(CLOCK_REALTIME, &ts);
      return ((double) ts.tv_sec) + ((double) ts.tv_nsec ) / (double) 1E9;
           - StartingTime;
   }
}
*/

#else

#include <windows.h>
double primme_wTimer() {
   return GetTickCount();
}

#endif
