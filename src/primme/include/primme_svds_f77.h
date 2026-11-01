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
C  File: primme_svds_f77.h
C  
C  Purpose - Main header with the PRIMME SVDS F77 interface functions.
C  
C******************************************************************************

C-------------------------------------------------------
C     Defining easy to remember labels for setting the 
C     method in primme_svds_set_method from Fortran
C-------------------------------------------------------
      integer
     : primme_svds_default,
     : primme_svds_hybrid,
     : primme_svds_normalequations,
     : primme_svds_augmented

      parameter(
     : primme_svds_default = 0,
     : primme_svds_hybrid = 1,
     : primme_svds_normalequations = 2,
     : primme_svds_augmented = 3
     :)

C-------------------------------------------------------
C     Defining easy to remember labels for setting the 
C     members of the primme_svds structure from Fortran
C-------------------------------------------------------
 
      INTEGER*8
     : PRIMME_SVDS_primme                       , 
     : PRIMME_SVDS_primmeStage2                 , 
     : PRIMME_SVDS_m                            , 
     : PRIMME_SVDS_n                            , 
     : PRIMME_SVDS_matrixMatvec                 , 
     : PRIMME_SVDS_matrixMatvec_type            , 
     : PRIMME_SVDS_applyPreconditioner          , 
     : PRIMME_SVDS_applyPreconditioner_type     , 
     : PRIMME_SVDS_numProcs                     , 
     : PRIMME_SVDS_procID                       ,
     : PRIMME_SVDS_mLocal                       ,
     : PRIMME_SVDS_nLocal                       ,
     : PRIMME_SVDS_commInfo                     ,
     : PRIMME_SVDS_globalSumReal                ,
     : PRIMME_SVDS_globalSumReal_type           ,
     : PRIMME_SVDS_broadcastReal                ,
     : PRIMME_SVDS_broadcastReal_type           ,
     : PRIMME_SVDS_numSvals                     ,
     : PRIMME_SVDS_target                       ,
     : PRIMME_SVDS_numTargetShifts              ,
     : PRIMME_SVDS_targetShifts                 ,
     : PRIMME_SVDS_method                       ,
     : PRIMME_SVDS_methodStage2                 ,
     : PRIMME_SVDS_matrix                       ,
     : PRIMME_SVDS_preconditioner               ,
     : PRIMME_SVDS_locking                      ,
     : PRIMME_SVDS_numOrthoConst                ,
     : PRIMME_SVDS_aNorm                        ,
     : PRIMME_SVDS_eps                          ,
     : PRIMME_SVDS_precondition                 ,
     : PRIMME_SVDS_initSize                     ,
     : PRIMME_SVDS_maxBasisSize                 ,
     : PRIMME_SVDS_maxBlockSize                 ,
     : PRIMME_SVDS_maxMatvecs                   ,
     : PRIMME_SVDS_iseed                        ,
     : PRIMME_SVDS_printLevel                   ,
     : PRIMME_SVDS_internalPrecision            ,
     : PRIMME_SVDS_outputFile                   ,
     : PRIMME_SVDS_stats_numOuterIterations     ,
     : PRIMME_SVDS_stats_numRestarts            ,
     : PRIMME_SVDS_stats_numMatvecs             ,
     : PRIMME_SVDS_stats_numPreconds            ,
     : PRIMME_SVDS_stats_numGlobalSum           ,
     : PRIMME_SVDS_stats_volumeGlobalSum        ,
     : PRIMME_SVDS_stats_numBroadcast           ,
     : PRIMME_SVDS_stats_volumeBroadcast        ,
     : PRIMME_SVDS_stats_numOrthoInnerProds     ,
     : PRIMME_SVDS_stats_elapsedTime            ,
     : PRIMME_SVDS_stats_timeMatvec             ,
     : PRIMME_SVDS_stats_timePrecond            ,
     : PRIMME_SVDS_stats_timeOrtho              ,
     : PRIMME_SVDS_stats_timeGlobalSum          ,
     : PRIMME_SVDS_stats_timeBroadcast          ,
     : PRIMME_SVDS_stats_lockingIssue           ,
     : PRIMME_SVDS_convTestFun                  ,
     : PRIMME_SVDS_convTestFun_type             ,
     : PRIMME_SVDS_convtest                     ,
     : PRIMME_SVDS_monitorFun                   ,
     : PRIMME_SVDS_monitorFun_type              ,
     : PRIMME_SVDS_monitor                      ,
     : PRIMME_SVDS_queue                        ,
     : PRIMME_SVDS_profile                       

      parameter(
     : PRIMME_SVDS_primme                       = 1,
     : PRIMME_SVDS_primmeStage2                 = 2,
     : PRIMME_SVDS_m                            = 3,
     : PRIMME_SVDS_n                            = 4,
     : PRIMME_SVDS_matrixMatvec                 = 5,
     : PRIMME_SVDS_matrixMatvec_type            = 6,
     : PRIMME_SVDS_applyPreconditioner          = 7,
     : PRIMME_SVDS_applyPreconditioner_type     = 8,
     : PRIMME_SVDS_numProcs                     = 9,
     : PRIMME_SVDS_procID                       = 10,
     : PRIMME_SVDS_mLocal                       = 11,
     : PRIMME_SVDS_nLocal                       = 12,
     : PRIMME_SVDS_commInfo                     = 13,
     : PRIMME_SVDS_globalSumReal                = 14,
     : PRIMME_SVDS_globalSumReal_type           = 15,
     : PRIMME_SVDS_broadcastReal                = 16,
     : PRIMME_SVDS_broadcastReal_type           = 17,
     : PRIMME_SVDS_numSvals                     = 18,
     : PRIMME_SVDS_target                       = 19,
     : PRIMME_SVDS_numTargetShifts              = 20,
     : PRIMME_SVDS_targetShifts                 = 21,
     : PRIMME_SVDS_method                       = 22,
     : PRIMME_SVDS_methodStage2                 = 23,
     : PRIMME_SVDS_matrix                       = 24,
     : PRIMME_SVDS_preconditioner               = 25,
     : PRIMME_SVDS_locking                      = 26,
     : PRIMME_SVDS_numOrthoConst                = 27,
     : PRIMME_SVDS_aNorm                        = 28,
     : PRIMME_SVDS_eps                          = 29,
     : PRIMME_SVDS_precondition                 = 30,
     : PRIMME_SVDS_initSize                     = 31,
     : PRIMME_SVDS_maxBasisSize                 = 32,
     : PRIMME_SVDS_maxBlockSize                 = 33,
     : PRIMME_SVDS_maxMatvecs                   = 34,
     : PRIMME_SVDS_iseed                        = 35,
     : PRIMME_SVDS_printLevel                   = 36,
     : PRIMME_SVDS_internalPrecision            = 37,
     : PRIMME_SVDS_outputFile                   = 38,
     : PRIMME_SVDS_stats_numOuterIterations     = 39,
     : PRIMME_SVDS_stats_numRestarts            = 40,
     : PRIMME_SVDS_stats_numMatvecs             = 41,
     : PRIMME_SVDS_stats_numPreconds            = 42,
     : PRIMME_SVDS_stats_numGlobalSum           = 43,
     : PRIMME_SVDS_stats_volumeGlobalSum        = 44,
     : PRIMME_SVDS_stats_numBroadcast           = 45,
     : PRIMME_SVDS_stats_volumeBroadcast        = 46,
     : PRIMME_SVDS_stats_numOrthoInnerProds     = 47,
     : PRIMME_SVDS_stats_elapsedTime            = 48,
     : PRIMME_SVDS_stats_timeMatvec             = 49,
     : PRIMME_SVDS_stats_timePrecond            = 50,
     : PRIMME_SVDS_stats_timeOrtho              = 51,
     : PRIMME_SVDS_stats_timeGlobalSum          = 52,
     : PRIMME_SVDS_stats_timeBroadcast          = 53,
     : PRIMME_SVDS_stats_lockingIssue           = 54,
     : PRIMME_SVDS_convTestFun                  = 55,
     : PRIMME_SVDS_convTestFun_type             = 56,
     : PRIMME_SVDS_convtest                     = 57,
     : PRIMME_SVDS_monitorFun                   = 58,
     : PRIMME_SVDS_monitorFun_type              = 59,
     : PRIMME_SVDS_monitor                      = 60,
     : PRIMME_SVDS_queue                        = 61,
     : PRIMME_SVDS_profile                      = 62 
     :)

C-------------------------------------------------------
C    Defining easy to remember labels for setting the 
C    enum members for targeting and operator
C-------------------------------------------------------

      integer 
     : primme_svds_largest,
     : primme_svds_smallest,
     : primme_svds_closest_abs,
     : primme_svds_op_none,
     : primme_svds_op_AtA,
     : primme_svds_op_AAt,
     : primme_svds_op_augmented

      parameter(
     : primme_svds_largest = 0,
     : primme_svds_smallest = 1,
     : primme_svds_closest_abs = 2,
     : primme_svds_op_none = 0,
     : primme_svds_op_AtA = 1,
     : primme_svds_op_AAt = 2,
     : primme_svds_op_augmented = 3
     :)
