##################################################
# SCIFOR PROJECT
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.10.0)
PROJECT(scifor Fortran)

MESSAGE(STATUS "OS: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}")
SET(BUILD_TYPE "RELEASE" CACHE STRING "Build type")
STRING(TOUPPER "${BUILD_TYPE}" BUILD_TYPE)
MESSAGE(STATUS "BUILD_TYPE=${BUILD_TYPE}")
SET(CMAKE_BUILD_TYPE ${BUILD_TYPE} CACHE INTERNAL "Build type" FORCE)

# Add our local CMAKE modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)

# Set some options the user may choose
OPTION(USE_MPI "Use the MPI library for parallelization" ON)
OPTION(VERBOSE "Build a verbose Makefile" OFF)
OPTION(WITH_BLAS_LAPACK "Build using internal Blas/Lapack" OFF)
OPTION(WITH_SCALAPACK "Build ScaLapack support" ON)
IF(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE ON)
ELSE(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE OFF)
ENDIF(VERBOSE)



OPTION(LONG_PREFIX "Set long (ON) or short (OFF) installation prefix DEFAULT_PREFIX(/<fc_id>/[<git_branch>{/<debug>}]/<version>):  " ON)
if(LONG_PREFIX)
  message(STATUS "LONG_PREFIX is ON")
else()
  message(STATUS "LONG_PREFIX is OFF")
endif()

############################################################
# Look for MPI if needed
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/MpiConfig.cmake)


############################################################
# Check for Blas/Lapack/ScaLapack from different sources
############################################################
# Try to find MKL: IF found use it ELSE try finding Blas/Lapack
# in the OS, ELSE using internal Blas/Lapack
IF(NOT WITH_BLAS_LAPACK)
  FIND_PACKAGE(MKL QUIET)
  IF(MKL_FOUND)
    MESSAGE(STATUS "MKL found at: ${MKL_ROOT_DIR}")
    SET(BLA_LIBRARIES ${MKL_LIBRARIES})    
  ELSE()
    FIND_PACKAGE(BLAS QUIET)
    FIND_PACKAGE(LAPACK QUIET)
    IF(BLAS_FOUND AND LAPACK_FOUND)
      MESSAGE(STATUS "BLAS found as: ${BLAS_LIBRARIES}")
      MESSAGE(STATUS "LAPACK found as: ${LAPACK_LIBRARIES}")
      SET(BLA_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    ELSE()
      SET(WITH_BLAS_LAPACK "ON")
    ENDIF()
  ENDIF()
ENDIF()

SET(SCALAPACK_FOUND FALSE)
IF(WITH_SCALAPACK)
  IF(MKL_SCALAPACK_FOUND)
    MESSAGE(STATUS "MKL supports Scalapack+Blacs")
    SET(SCALAPACK_FOUND TRUE)
  ELSE()
    FIND_PACKAGE(SCALAPACK QUIET)	#TRY TO FIND SCALAPACK IN THE SYSTEM
    IF(SCALAPACK_FOUND)
      MESSAGE(STATUS "Scalapack found: ${SCALAPACK_LIBRARIES}")
    ENDIF()
  ENDIF()  
ENDIF()





############################################################
# Configure Fortran Compiler Main options
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/FortranCompilerConfig.cmake)
IF((WITH_SCALAPACK) AND (SCALAPACK_FOUND))
  ADD_DEFINITIONS(-D_SCALAPACK)
ENDIF()




############################################################
# Check for .git and extract the current SHA1 revision.
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/GitConfig.cmake)




############################################################
# Define the actual files and folders that make up the build
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/MainConfig.cmake)



############################################################
############################################################
############################################################
# DEFINE HERE THE DETAILS OF THE LIBRARY CONFIGURATION 

# Define some directories
SET(LIB_SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB_ETC ${CMAKE_SOURCE_DIR}/etc)
SET(LIB_ENV ${CMAKE_SOURCE_DIR}/etc/environment_modules)
SET(LIB_BIN ${CMAKE_SOURCE_DIR}/bin)
SET(LIB_TMP_INC ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME})
SET(LIB_TMP_ETC ${CMAKE_BINARY_DIR}/etc)
SET(LIB_TMP_BIN ${CMAKE_BINARY_DIR}/bin)
SET(LIB_TMP_VER ${CMAKE_BINARY_DIR}/version)
SET(LIB_TARGET_LIB ${CMAKE_INSTALL_PREFIX}/lib)
SET(LIB_TARGET_INC ${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME})
SET(LIB_TARGET_ETC ${CMAKE_INSTALL_PREFIX}/etc)
SET(LIB_TARGET_BIN ${CMAKE_INSTALL_PREFIX}/bin)
SET(LIB_TARGET_DIR ${CMAKE_INSTALL_PREFIX})

# Have the .mod files placed in the INC folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB_TMP_INC})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_Fortran_MODULE_DIRECTORY})


IF(WITH_BLAS_LAPACK)
  SET(SRC_DIR_BLAS ${LIB_SRC}/blas)
  ADD_SUBDIRECTORY(${SRC_DIR_BLAS})

  SET(SRC_DIR_LAPACK ${LIB_SRC}/lapack)
  ADD_SUBDIRECTORY(${SRC_DIR_LAPACK})
ENDIF()


SET(SRC_DIR_ARPACK ${LIB_SRC}/arpack)
ADD_SUBDIRECTORY(${SRC_DIR_ARPACK})


SET(SRC_DIR_FFTPACK ${LIB_SRC}/fftpack)
ADD_SUBDIRECTORY(${SRC_DIR_FFTPACK})


SET(SRC_DIR_MINPACK ${LIB_SRC}/minpack)
ADD_SUBDIRECTORY(${SRC_DIR_MINPACK})


SET(SRC_DIR_QUADPACK ${LIB_SRC}/quadpack)
ADD_SUBDIRECTORY(${SRC_DIR_QUADPACK})


SET(SRC_DIR_SF_MPI ${LIB_SRC}/SF_MPI)
ADD_SUBDIRECTORY(${SRC_DIR_SF_MPI})

#########################################

SET(SRC_DIR_SF_ARRAYS ${LIB_SRC}/SF_ARRAYS)
ADD_SUBDIRECTORY(${SRC_DIR_SF_ARRAYS})


SET(SRC_DIR_SF_COLORS ${LIB_SRC}/SF_COLORS)
ADD_SUBDIRECTORY(${SRC_DIR_SF_COLORS})


SET(SRC_DIR_SF_CONSTANTS ${LIB_SRC}/SF_CONSTANTS)
ADD_SUBDIRECTORY(${SRC_DIR_SF_CONSTANTS})


SET(SRC_DIR_SF_FONTS ${LIB_SRC}/SF_FONTS)
ADD_SUBDIRECTORY(${SRC_DIR_SF_FONTS})


SET(SRC_DIR_SF_MISC ${LIB_SRC}/SF_MISC)
ADD_SUBDIRECTORY(${SRC_DIR_SF_MISC})

SET(SRC_DIR_SF_IOTOOLS ${LIB_SRC}/SF_IOTOOLS)
ADD_SUBDIRECTORY(${SRC_DIR_SF_IOTOOLS})

SET(SRC_DIR_SF_DERIVATE ${LIB_SRC}/SF_DERIVATE)
ADD_SUBDIRECTORY(${SRC_DIR_SF_DERIVATE})


SET(SRC_DIR_SF_INTEGRATE ${LIB_SRC}/SF_INTEGRATE)
ADD_SUBDIRECTORY(${SRC_DIR_SF_INTEGRATE})


SET(SRC_DIR_SF_INTERPOLATE ${LIB_SRC}/SF_INTERPOLATE)
ADD_SUBDIRECTORY(${SRC_DIR_SF_INTERPOLATE})


SET(SRC_DIR_SF_PARSE_INPUT ${LIB_SRC}/SF_PARSE_INPUT)
ADD_SUBDIRECTORY(${SRC_DIR_SF_PARSE_INPUT})


SET(SRC_DIR_SF_TIMER ${LIB_SRC}/SF_TIMER)   #<- MPI
ADD_SUBDIRECTORY(${SRC_DIR_SF_TIMER})


SET(SRC_DIR_SF_LINALG ${LIB_SRC}/SF_LINALG) #<- MPI(BLACS)
ADD_SUBDIRECTORY(${SRC_DIR_SF_LINALG})


SET(SRC_DIR_SF_OPTIMIZE ${LIB_SRC}/SF_OPTIMIZE) #<- CONSTANTS,LINALG
ADD_SUBDIRECTORY(${SRC_DIR_SF_OPTIMIZE})


SET(SRC_DIR_SF_SPECIAL ${LIB_SRC}/SF_SPECIAL) #<- INTEGRATE
ADD_SUBDIRECTORY(${SRC_DIR_SF_SPECIAL})


SET(SRC_DIR_SF_FFT ${LIB_SRC}/SF_FFT) #<- INTEGRATE,ARRAYS,CONSTANTS
ADD_SUBDIRECTORY(${SRC_DIR_SF_FFT})


SET(SRC_DIR_SF_RANDOM ${LIB_SRC}/SF_RANDOM) #<- INTEGRATE,ARRAYS,IOTOOLS,LINALG
ADD_SUBDIRECTORY(${SRC_DIR_SF_RANDOM})


SET(SRC_DIR_SF_SP_LINALG ${LIB_SRC}/SF_SP_LINALG) #<- MISC, RANDOM, LINALG, MPI
ADD_SUBDIRECTORY(${SRC_DIR_SF_SP_LINALG})



ADD_SUBDIRECTORY(${LIB_SRC})

  
ADD_DEPENDENCIES(SF_TIMERLIB
  SF_MPILIB
)

ADD_DEPENDENCIES(SF_LINALGLIB
  SF_MPILIB
)

ADD_DEPENDENCIES(SF_OPTIMIZELIB
  SF_CONSTANTSLIB
  SF_LINALGLIB
)

ADD_DEPENDENCIES(SF_SPECIALLIB
  SF_INTEGRATELIB
)

ADD_DEPENDENCIES(SF_FFTLIB
  SF_ARRAYSLIB
  SF_CONSTANTSLIB
  SF_INTEGRATELIB
)

ADD_DEPENDENCIES(SF_RANDOMLIB
  SF_INTEGRATELIB
  SF_IOTOOLSLIB
  SF_ARRAYSLIB
  SF_LINALGLIB
)


ADD_DEPENDENCIES(SF_SP_LINALGLIB
  SF_MISCLIB
  SF_RANDOMLIB
  SF_LINALGLIB
  SF_MPILIB
)


ADD_DEPENDENCIES(SCIFORLIB
  SF_ARRAYSLIB
  SF_COLORSLIB
  SF_CONSTANTSLIB
  SF_FONTSLIB
  SF_MISCLIB
  SF_TIMERLIB
  SF_MPILIB
  SF_IOTOOLSLIB
  SF_PARSE_INPUTLIB
  SF_DERIVATELIB
  SF_LINALGLIB
  SF_OPTIMIZELIB
  SF_SPECIALLIB 
  SF_INTEGRATELIB
  SF_INTERPOLATELIB
  SF_RANDOMLIB
  SF_FFTLIB
  SF_SP_LINALGLIB
)



  
# Build the  library and have it placed in the target lib folder
IF(WITH_BLAS_LAPACK)
  ADD_LIBRARY(scifor 
    $<TARGET_OBJECTS:BLASLIB>
    $<TARGET_OBJECTS:LAPACKLIB>
    $<TARGET_OBJECTS:ARPACKLIB>
    $<TARGET_OBJECTS:FFTPACKLIB>
    $<TARGET_OBJECTS:MINPACKLIB>
    $<TARGET_OBJECTS:QUADPACKLIB>
    $<TARGET_OBJECTS:SF_ARRAYSLIB>
    $<TARGET_OBJECTS:SF_COLORSLIB>
    $<TARGET_OBJECTS:SF_CONSTANTSLIB>
    $<TARGET_OBJECTS:SF_FONTSLIB>
    $<TARGET_OBJECTS:SF_MISCLIB>
    $<TARGET_OBJECTS:SF_TIMERLIB>    
    $<TARGET_OBJECTS:SF_MPILIB>
    $<TARGET_OBJECTS:SF_PARSE_INPUTLIB>
    $<TARGET_OBJECTS:SF_IOTOOLSLIB>
    $<TARGET_OBJECTS:SF_LINALGLIB>
    $<TARGET_OBJECTS:SF_DERIVATELIB>
    $<TARGET_OBJECTS:SF_OPTIMIZELIB>
    $<TARGET_OBJECTS:SF_INTEGRATELIB>
    $<TARGET_OBJECTS:SF_SPECIALLIB>
    $<TARGET_OBJECTS:SF_FFTLIB>
    $<TARGET_OBJECTS:SF_INTERPOLATELIB>
    $<TARGET_OBJECTS:SF_RANDOMLIB>
    $<TARGET_OBJECTS:SF_SP_LINALGLIB>
    $<TARGET_OBJECTS:SCIFORLIB>
    )
ELSE()
  ADD_LIBRARY(scifor 
    $<TARGET_OBJECTS:ARPACKLIB>
    $<TARGET_OBJECTS:FFTPACKLIB>
    $<TARGET_OBJECTS:MINPACKLIB>
    $<TARGET_OBJECTS:QUADPACKLIB>
    $<TARGET_OBJECTS:SF_ARRAYSLIB>
    $<TARGET_OBJECTS:SF_COLORSLIB>
    $<TARGET_OBJECTS:SF_CONSTANTSLIB>
    $<TARGET_OBJECTS:SF_FONTSLIB>
    $<TARGET_OBJECTS:SF_MISCLIB>
    $<TARGET_OBJECTS:SF_TIMERLIB>
    $<TARGET_OBJECTS:SF_MPILIB>
    $<TARGET_OBJECTS:SF_IOTOOLSLIB>
    $<TARGET_OBJECTS:SF_PARSE_INPUTLIB>
    $<TARGET_OBJECTS:SF_LINALGLIB>
    $<TARGET_OBJECTS:SF_DERIVATELIB>
    $<TARGET_OBJECTS:SF_OPTIMIZELIB>
    $<TARGET_OBJECTS:SF_INTEGRATELIB>
    $<TARGET_OBJECTS:SF_SPECIALLIB>
    $<TARGET_OBJECTS:SF_FFTLIB>
    $<TARGET_OBJECTS:SF_INTERPOLATELIB>
    $<TARGET_OBJECTS:SF_RANDOMLIB>
    $<TARGET_OBJECTS:SF_SP_LINALGLIB>
    $<TARGET_OBJECTS:SCIFORLIB>
    )
ENDIF()



#Build the PKG-CONFIG file
INCLUDE(${CMAKE_MODULE_PATH}/BuildPkgConfigFile.cmake)
SET(TMP_PKCONFIG_FILE ${LIB_TMP_ETC}/${PROJECT_NAME}.pc)
BUILD_PKCONFIG(${TMP_PKCONFIG_FILE})

############################################################
############################################################
############################################################

INCLUDE(${CMAKE_MODULE_PATH}/PostConfig.cmake)
