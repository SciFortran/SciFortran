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
 * File: memman.h
 *
 * Purpose - Header file for memman.c
 *
 ******************************************************************************/


#ifndef memman_H
#define memman_H

struct primme_context_str;

typedef struct primme_alloc_str {
   void *p;                         /* Allocated pointer */
   int (*free_fn)(void *, struct primme_context_str);
                                    /* Function to free pointer */
   struct primme_alloc_str *prev;   /* Previous allocation */
#ifndef NDEBUG
   const char* debug;               /* String identifying the code that */
                                    /* generated the allocation         */
#endif
} primme_alloc;

typedef struct primme_frame_str {
   primme_alloc *prev_alloc;        /* Pointer to last allocation */
   int keep_frame;                  /* Flag to merge allocation in this frame */
                                    /* into the previous frame */
   struct primme_frame_str *prev;   /* Pointer to previous frame */
                        
} primme_frame;

int Mem_push_frame(struct primme_context_str *ctx);
int Mem_pop_frame(struct primme_context_str *ctx);
int Mem_pop_clean_frame(struct primme_context_str ctx);
int Mem_keep_frame(struct primme_context_str ctx);
int Mem_debug_frame(const char *debug, struct primme_context_str ctx);
typedef int (*free_fn_type)(void *, struct primme_context_str);
int Mem_register_alloc(void *p, free_fn_type free_fn, struct primme_context_str ctx);
int Mem_deregister_alloc(void *p, struct primme_context_str ctx);

#endif
