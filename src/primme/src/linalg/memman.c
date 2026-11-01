/*******************************************************************************
 * Copyright (c) 2017, College of William & Mary
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
 * File: memman.c
 *
 * Purpose - Define functions to track allocated memory and free them in
 *           case of error.
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../linalg/memman.c"
#endif

#include <stdlib.h>   /* free */
#include <assert.h>
#include <math.h>
#include "common.h"
#include "memman.h"


/*******************************************************************************
 * Subroutine Mem_push_frame - Push a new frame in the context.
 * 
 * INPUT PARAMETERS
 * ----------------------------------
 * ctx      context
 *
 ******************************************************************************/

static int free_dummy(void *p, primme_context ctx) {
   (void)ctx;
   free(p);
   return 0;
}

int Mem_push_frame(primme_context *ctx) {

   /* Quick exit */

   if (!ctx) return 0;

   primme_frame *f = NULL;
   primme_alloc *a = NULL;
   if (MALLOC_PRIMME(1, &f) == 0 && MALLOC_PRIMME(1, &a) == 0) {
      f->prev_alloc = a;
      f->keep_frame = 0;
      f->prev = ctx->mm;
      a->p = f;
      a->free_fn = free_dummy;
      a->prev = NULL;
#ifndef NDEBUG
      a->debug = NULL;
#endif
      ctx->mm = f;
   } else {
      if (f) free(f);
      if (a) free(a);
      return -1;
   }

   return 0;
}

/*******************************************************************************
 * Subroutine Mem_pop_frame - Remove the last frame pushed in the context.
 * 
 * INPUT PARAMETERS
 * ----------------------------------
 * ctx      context
 *
 ******************************************************************************/

int Mem_pop_frame(primme_context *ctx) {

   /* Quick exit */

   if (!ctx || !ctx->mm) return 0;

   /* Show message if there is no previous frame  and they want to keep the
    * frame. */

   if (ctx->mm->keep_frame && !ctx->mm->prev && !ctx->mm->prev_alloc) {
      PRINTFALLCTX(*ctx, 1, "Warning: no frame where to keep allocations");
      return -1;
   }

   /* Store the reference to the previous frame. The current frame may be freed
    * in the followings commands. */

   primme_frame *mm_prev = ctx->mm->prev;

   /* If it is asked to keep the frame, transfer all registers to the   */
   /* previous frame.                                                   */

   if (ctx->mm->keep_frame && ctx->mm->prev) {
      primme_alloc *a = ctx->mm->prev_alloc;
      while (a) {
         primme_alloc *a_prev = a->prev;
         a->prev = ctx->mm->prev->prev_alloc;
         ctx->mm->prev->prev_alloc = a;
         a = a_prev;
      }
   }

   /* If not, the function should have freed all allocations. */

   else {

      /* Warning about allocations not made by Mem_push_frame */

#ifndef NDEBUG
      primme_alloc *a = ctx->mm->prev_alloc;
      while (a) {
         if (a->free_fn != free_dummy) {
            PRINTFALLCTX(*ctx, 1,
                  "Warning: the allocation at %s has not been freed",
               a->debug ? a->debug : "unknown");
            assert(0);
         }
         a = a->prev;
      }
#endif 

      /* If there are allocations, just free them */
 
      Mem_pop_clean_frame(*ctx);
   }

   /* Set the current frame as the previous one */

   ctx->mm = mm_prev;

   return 0;
}

/*******************************************************************************
 * Subroutine Mem_pop_clean_frame - Free all allocations registered already.
 *
 * INPUT PARAMETERS
 * ----------------------------------
 * ctx      context
 *
 ******************************************************************************/

int Mem_pop_clean_frame(primme_context ctx) {

   primme_alloc *a = ctx.mm ? ctx.mm->prev_alloc : NULL;
   if (ctx.mm) ctx.mm->prev_alloc = NULL;
   while (a) {
      primme_alloc *a_prev = a->prev;
      if (a->p) a->free_fn(a->p, ctx);
      free(a);
      a = a_prev;
   }

   return 0;
}

/*******************************************************************************
 * Subroutine Mem_keep_frame - Ask to not remove the last frame pushed.
 * 
 * INPUT/OUTPUT PARAMETERS
 * ----------------------------------
 * ctx  context
 *
 ******************************************************************************/

int Mem_keep_frame(primme_context ctx) {

   assert(ctx.mm);
   ctx.mm->keep_frame = 1;

   return 0;
}

/*******************************************************************************
 * Subroutine Mem_register_alloc - Register a pointer been allocated and the
 *    function to free the pointer.
 * 
 * INPUT PARAMETERS
 * ----------------------------------
 * p        Pointer been allocated
 * free_fn  Function to free the pointer
 * ctx      context
 *
 ******************************************************************************/

int Mem_register_alloc(void *p, free_fn_type free_fn, primme_context ctx) {

   assert(ctx.mm);

   primme_alloc *prev_alloc = ctx.mm->prev_alloc, *a;
   CHKERR(MALLOC_PRIMME(1, &a)); 
   a->p = p;
   a->free_fn = free_fn;
   a->prev = prev_alloc;
#ifndef NDEBUG
   a->debug = NULL;
#endif
   ctx.mm->prev_alloc = a;

   return 0;
}

/*******************************************************************************
 * Subroutine Mem_deregister_alloc - Remove the pointer from the current frame
 * 
 * INPUT PARAMETERS
 * ----------------------------------
 * p        Pointer been removed
 * ctx      context
 *
 ******************************************************************************/

int Mem_deregister_alloc(void *p, primme_context ctx) {

   if (!p) return 0;

   assert(ctx.mm);

   /* Find the register with the pointer p and who points out that register */

   primme_frame *f = ctx.mm;
   primme_alloc *a = NULL, **prev = NULL;
   while (f) {
      a = f->prev_alloc;
      prev = &f->prev_alloc;
      while (a && a->p != p) {
         prev = &a->prev;
         a = a->prev;
      }
      if (a) break;
      f = f->prev;
   }

   /* Remove the register and link the list properly */

   assert(a);
   *prev = a->prev;
   free(a);

   return 0;
}

/*******************************************************************************
 * Subroutine Mem_register_alloc - Register a pointer been allocated and the
 *    function to free the pointer.
 * 
 * INPUT PARAMETERS
 * ----------------------------------
 * p        Pointer been allocated
 * free_fn  Function to free the pointer
 * ctx      context
 *
 ******************************************************************************/

int Mem_debug_frame(const char *debug, primme_context ctx) {

   /* Quick exit */

   if (!ctx.mm) return 0;

#ifndef NDEBUG
   primme_alloc *a = ctx.mm->prev_alloc;
   while(a) {
      if (!a->debug) a->debug = debug;
      a = a->prev;
   }
#else
   (void)debug;
#endif

   return 0;
}
