C*******************************************************************************
C  Copyright (c) 2018, College of William & Mary                                   
C  All rights reserved.                                                            
C                                                                                  
C  Redistribution and use in source and binary forms, with or without
C  modification, are permitted provided that the following conditions are met:     
C      * Redistributions of source code must retain the above copyright
C        notice, this list of conditions and the following disclaimer.             
C      * Redistributions in binary form must reproduce the above copyright         
C        notice, this list of conditions and the following disclaimer in the       
C        documentation and/or other materials provided with the distribution.      
C      * Neither the name of the College of William & Mary nor the
C        names of its contributors may be used to endorse or promote products      
C        derived from this software without specific prior written permission.     
C                                                                                  
C  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
C  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
C  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          
C  DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY       
C  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      
C  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
C  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
C  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
C  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C  
C  PRIMME: https://github.com/primme/primme
C  Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
C*******************************************************************************
C  File: primme_eigs_f77.h
C  
C  Purpose - Main header with the PRIMME EIGS F77 interface functions.
C  
C******************************************************************************

C-------------------------------------------------------
C     Defining easy to remember labels for setting the 
C     method in primme_set_method from Fortran
C-------------------------------------------------------
      integer
     : PRIMME_DEFAULT_METHOD,
     : PRIMME_DYNAMIC,
     : PRIMME_DEFAULT_MIN_TIME,
     : PRIMME_DEFAULT_MIN_MATVECS,
     : PRIMME_Arnoldi,
     : PRIMME_GD,
     : PRIMME_GD_plusK,
     : PRIMME_GD_Olsen_plusK,
     : PRIMME_JD_Olsen_plusK,
     : PRIMME_RQI,
     : PRIMME_JDQR,
     : PRIMME_JDQMR,
     : PRIMME_JDQMR_ETol,
     : PRIMME_STEEPEST_DESCENT,
     : PRIMME_LOBPCG_OrthoBasis,
     : PRIMME_LOBPCG_OrthoBasis_Window

      parameter(
     : PRIMME_DEFAULT_METHOD = 0,
     : PRIMME_DYNAMIC = 1,
     : PRIMME_DEFAULT_MIN_TIME = 2,
     : PRIMME_DEFAULT_MIN_MATVECS = 3,
     : PRIMME_Arnoldi = 4,
     : PRIMME_GD = 5,
     : PRIMME_GD_plusK = 6,
     : PRIMME_GD_Olsen_plusK = 7,
     : PRIMME_JD_Olsen_plusK = 8,
     : PRIMME_RQI = 9,
     : PRIMME_JDQR = 10,
     : PRIMME_JDQMR = 11,
     : PRIMME_JDQMR_ETol = 12,
     : PRIMME_STEEPEST_DESCENT = 13,
     : PRIMME_LOBPCG_OrthoBasis = 14,
     : PRIMME_LOBPCG_OrthoBasis_Window = 15
     :)

C-------------------------------------------------------
C     Defining easy to remember labels for setting the 
C     members of the primme structure from Fortran
C-------------------------------------------------------
      INTEGER
     : PRIMME_n                                      , 
     : PRIMME_matrixMatvec                           , 
     : PRIMME_matrixMatvec_type                      , 
     : PRIMME_applyPreconditioner                    , 
     : PRIMME_applyPreconditioner_type               , 
     : PRIMME_massMatrixMatvec                       , 
     : PRIMME_massMatrixMatvec_type                  , 
     : PRIMME_numProcs                               , 
     : PRIMME_procID                                 , 
     : PRIMME_commInfo                               ,
     : PRIMME_nLocal                                 ,
     : PRIMME_globalSumReal                          ,
     : PRIMME_globalSumReal_type                     ,
     : PRIMME_broadcastReal                          ,
     : PRIMME_broadcastReal_type                     ,
     : PRIMME_numEvals                               ,
     : PRIMME_target                                 ,
     : PRIMME_numTargetShifts                        ,
     : PRIMME_targetShifts                           ,
     : PRIMME_locking                                ,
     : PRIMME_initSize                               ,
     : PRIMME_numOrthoConst                          ,
     : PRIMME_maxBasisSize                           ,
     : PRIMME_minRestartSize                         ,
     : PRIMME_maxBlockSize                           ,
     : PRIMME_maxMatvecs                             ,
     : PRIMME_maxOuterIterations                     ,
     : PRIMME_iseed                                  ,
     : PRIMME_aNorm                                  ,
     : PRIMME_BNorm                                  ,
     : PRIMME_invBNorm                               ,
     : PRIMME_eps                                    ,
     : PRIMME_orth                                   ,
     : PRIMME_internalPrecision                      ,
     : PRIMME_printLevel                             ,
     : PRIMME_outputFile                             ,
     : PRIMME_matrix                                 ,
     : PRIMME_massMatrix                             ,
     : PRIMME_preconditioner                         ,
     : PRIMME_ShiftForPreconditioner                 ,
     : PRIMME_initBasisMode                          ,
     : PRIMME_projectionParams_projection            ,
     : PRIMME_restartingParams_maxPrevRetain         ,
     : PRIMME_correctionParams_precondition          ,
     : PRIMME_correctionParams_robustShifts          ,
     : PRIMME_correctionParams_maxInnerIterations    ,
     : PRIMME_correctionParams_projectors_LeftQ      ,
     : PRIMME_correctionParams_projectors_LeftX      ,
     : PRIMME_correctionParams_projectors_RightQ     ,
     : PRIMME_correctionParams_projectors_RightX     ,
     : PRIMME_correctionParams_projectors_SkewQ      ,
     : PRIMME_correctionParams_projectors_SkewX      ,
     : PRIMME_correctionParams_convTest              ,
     : PRIMME_correctionParams_relTolBase            ,
     : PRIMME_stats_numOuterIterations               ,
     : PRIMME_stats_numRestarts                      ,
     : PRIMME_stats_numMatvecs                       ,
     : PRIMME_stats_numPreconds                      ,
     : PRIMME_stats_numGlobalSum                     ,
     : PRIMME_stats_volumeGlobalSum                  ,
     : PRIMME_stats_numBroadcast                     ,
     : PRIMME_stats_volumeBroadcast                  ,
     : PRIMME_stats_flopsDense                       ,
     : PRIMME_stats_numOrthoInnerProds               ,
     : PRIMME_stats_elapsedTime                      ,
     : PRIMME_stats_timeMatvec                       ,
     : PRIMME_stats_timePrecond                      ,
     : PRIMME_stats_timeOrtho                        ,
     : PRIMME_stats_timeGlobalSum                    ,
     : PRIMME_stats_timeBroadcast                    ,
     : PRIMME_stats_timeDense                        ,
     : PRIMME_stats_estimateMinEVal                  ,
     : PRIMME_stats_estimateMaxEVal                  ,
     : PRIMME_stats_estimateLargestSVal              ,
     : PRIMME_stats_estimateBNorm                    ,
     : PRIMME_stats_estimateInvBNorm                 ,
     : PRIMME_stats_maxConvTol                       ,
     : PRIMME_stats_lockingIssue                     ,
     : PRIMME_dynamicMethodSwitch                    ,
     : PRIMME_convTestFun                            ,
     : PRIMME_convTestFun_type                       ,
     : PRIMME_convtest                               ,
     : PRIMME_ldevecs                                ,
     : PRIMME_ldOPs                                  ,
     : PRIMME_monitorFun                             ,
     : PRIMME_monitorFun_type                        ,
     : PRIMME_monitor                                ,
     : PRIMME_queue                                  ,
     : PRIMME_profile                                

      parameter(
     : PRIMME_n                                      = 1  ,
     : PRIMME_matrixMatvec                           = 2  ,
     : PRIMME_matrixMatvec_type                      = 3  ,
     : PRIMME_applyPreconditioner                    = 4  ,
     : PRIMME_applyPreconditioner_type               = 5  ,
     : PRIMME_massMatrixMatvec                       = 6  ,
     : PRIMME_massMatrixMatvec_type                  = 7  ,
     : PRIMME_numProcs                               = 8  ,
     : PRIMME_procID                                 = 9  ,
     : PRIMME_commInfo                               = 10  ,
     : PRIMME_nLocal                                 = 11  ,
     : PRIMME_globalSumReal                          = 12  ,
     : PRIMME_globalSumReal_type                     = 13  ,
     : PRIMME_broadcastReal                          = 14  ,
     : PRIMME_broadcastReal_type                     = 15  ,
     : PRIMME_numEvals                               = 16  ,
     : PRIMME_target                                 = 17  ,
     : PRIMME_numTargetShifts                        = 18  ,
     : PRIMME_targetShifts                           = 19  ,
     : PRIMME_locking                                = 20  ,
     : PRIMME_initSize                               = 21  ,
     : PRIMME_numOrthoConst                          = 22  ,
     : PRIMME_maxBasisSize                           = 23  ,
     : PRIMME_minRestartSize                         = 24  ,
     : PRIMME_maxBlockSize                           = 25  ,
     : PRIMME_maxMatvecs                             = 26  ,
     : PRIMME_maxOuterIterations                     = 27  ,
     : PRIMME_iseed                                  = 28  ,
     : PRIMME_aNorm                                  = 29  ,
     : PRIMME_BNorm                                  = 30  ,
     : PRIMME_invBNorm                               = 31  ,
     : PRIMME_eps                                    = 32  ,
     : PRIMME_orth                                   = 33  ,
     : PRIMME_internalPrecision                      = 34  ,
     : PRIMME_printLevel                             = 35  ,
     : PRIMME_outputFile                             = 36  ,
     : PRIMME_matrix                                 = 37  ,
     : PRIMME_massMatrix                             = 38  ,
     : PRIMME_preconditioner                         = 39  ,
     : PRIMME_ShiftForPreconditioner                 = 40  ,
     : PRIMME_initBasisMode                          = 41  ,
     : PRIMME_projectionParams_projection            = 42  ,
     : PRIMME_restartingParams_maxPrevRetain         = 43  ,
     : PRIMME_correctionParams_precondition          = 44  ,
     : PRIMME_correctionParams_robustShifts          = 45  ,
     : PRIMME_correctionParams_maxInnerIterations    = 46  ,
     : PRIMME_correctionParams_projectors_LeftQ      = 47  ,
     : PRIMME_correctionParams_projectors_LeftX      = 48  ,
     : PRIMME_correctionParams_projectors_RightQ     = 49  ,
     : PRIMME_correctionParams_projectors_RightX     = 50  ,
     : PRIMME_correctionParams_projectors_SkewQ      = 51  ,
     : PRIMME_correctionParams_projectors_SkewX      = 52  ,
     : PRIMME_correctionParams_convTest              = 53  ,
     : PRIMME_correctionParams_relTolBase            = 54  ,
     : PRIMME_stats_numOuterIterations               = 55  ,
     : PRIMME_stats_numRestarts                      = 56  ,
     : PRIMME_stats_numMatvecs                       = 57  ,
     : PRIMME_stats_numPreconds                      = 58  ,
     : PRIMME_stats_numGlobalSum                     = 59  ,
     : PRIMME_stats_volumeGlobalSum                  = 60  ,
     : PRIMME_stats_numBroadcast                     = 61  ,
     : PRIMME_stats_volumeBroadcast                  = 62  ,
     : PRIMME_stats_flopsDense                       = 63  ,
     : PRIMME_stats_numOrthoInnerProds               = 64  ,
     : PRIMME_stats_elapsedTime                      = 65  ,
     : PRIMME_stats_timeMatvec                       = 66  ,
     : PRIMME_stats_timePrecond                      = 67  ,
     : PRIMME_stats_timeOrtho                        = 68  ,
     : PRIMME_stats_timeGlobalSum                    = 69  ,
     : PRIMME_stats_timeBroadcast                    = 70  ,
     : PRIMME_stats_timeDense                        = 71  ,
     : PRIMME_stats_estimateMinEVal                  = 72  ,
     : PRIMME_stats_estimateMaxEVal                  = 73  ,
     : PRIMME_stats_estimateLargestSVal              = 74  ,
     : PRIMME_stats_estimateBNorm                    = 75  ,
     : PRIMME_stats_estimateInvBNorm                 = 76  ,
     : PRIMME_stats_maxConvTol                       = 77  ,
     : PRIMME_stats_lockingIssue                     = 78  ,
     : PRIMME_dynamicMethodSwitch                    = 79  ,
     : PRIMME_convTestFun                            = 80  ,
     : PRIMME_convTestFun_type                       = 81  ,
     : PRIMME_convtest                               = 82  ,
     : PRIMME_ldevecs                                = 83  ,
     : PRIMME_ldOPs                                  = 84  ,
     : PRIMME_monitorFun                             = 85  ,
     : PRIMME_monitorFun_type                        = 86  ,
     : PRIMME_monitor                                = 87  ,
     : PRIMME_queue                                  = 88  ,
     : PRIMME_profile                                = 89  
     : )

C-------------------------------------------------------
C    Defining easy to remember labels for setting the 
C    enum members for targeting, restarting and innertest
C-------------------------------------------------------

      integer*8 
     : primme_smallest,
     : primme_largest,
     : primme_closest_geq,
     : primme_closest_leq,
     : primme_closest_abs,
     : primme_largest_abs,
     : primme_proj_default,
     : primme_proj_RR,
     : primme_proj_harmonic,
     : primme_proj_refined,
     : primme_init_default,
     : primme_init_krylov,
     : primme_init_random,
     : primme_init_user,
     : primme_full_LTolerance,
     : primme_decreasing_LTolerance,
     : primme_adaptive_ETolerance,
     : primme_adaptive,
     : primme_event_outer_iteration,
     : primme_event_inner_iteration,
     : primme_event_restart,
     : primme_event_reset,
     : primme_event_converged,
     : primme_event_locked,
     : primme_op_default, 
     : primme_op_half,
     : primme_op_float,
     : primme_op_double,
     : primme_op_quad,
     : primme_op_int

      parameter(
     : primme_smallest = 0,
     : primme_largest = 1,
     : primme_closest_geq = 2,
     : primme_closest_leq = 3,
     : primme_closest_abs = 4,
     : primme_largest_abs = 5,
     : primme_proj_default = 0,
     : primme_proj_RR = 1,
     : primme_proj_harmonic = 2,
     : primme_proj_refined = 3,
     : primme_init_default = 0,
     : primme_init_krylov = 1,
     : primme_init_random = 2,
     : primme_init_user = 3,
     : primme_full_LTolerance = 0,
     : primme_decreasing_LTolerance = 1,
     : primme_adaptive_ETolerance = 2,
     : primme_adaptive = 3,
     : primme_event_outer_iteration = 0,
     : primme_event_inner_iteration = 1,
     : primme_event_restart = 2,
     : primme_event_reset = 3,
     : primme_event_converged = 4,
     : primme_event_locked = 5,
     : primme_op_default = 0, 
     : primme_op_half = 1,
     : primme_op_float = 2,
     : primme_op_double = 3,
     : primme_op_quad = 4,
     : primme_op_int = 5
     : )
