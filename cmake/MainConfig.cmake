SET(USER_HOME $ENV{HOME})
SET(USER $ENV{USER})

IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET(CMAKE_Fortran_ARCHIVE_CREATE   "<CMAKE_AR> -r -LS <TARGET> <LINK_FLAGS> <OBJECTS>")
  SET(CMAKE_Fortran_ARCHIVE_FINISH   "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
ENDIF()

STRING(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" FC_ID)
STRING(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)

#DEFAULT LOCATION for installation is $HOME/opt/<libname>
SET(DEFAULT_LOCATION "$ENV{HOME}/opt")


#IF CIP IS DEFAULT (i.e. /usr/local in Unix/Linux) then we use DEFAULT_LOCATION + possible addendum
#if LONG_PREFIX=on
IF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  IF (LONG_PREFIX)
    SET(PATH_ID "${FC_ID}" CACHE INTERNAL "Default Path ID" FORCE)
    #IF not master branch, include simplified branch name
    IF( (NOT GIT_BRANCH MATCHES "master") )
      SET(PATH_ID  "${PATH_ID}/${GIT_BRANCH}")
    ENDIF()
    #IF DEBUG, add /debug to
    IF("${BUILD_TYPE}" MATCHES "debug")
      SET(PATH_ID  "${PATH_ID}/${BUILD_TYPE}")
    ENDIF()
    #set to $HOME/opt/<libname>(/<fc_id>/[<git_branch>{/<debug>}]/<version>)
    #set module name accordingly
    SET(DEFAULT_PREFIX "${DEFAULT_LOCATION}/${PROJECT_NAME}/${PATH_ID}/${VERSION}")    
    SET(MODULE_NAME "${PROJECT_NAME}/${PATH_ID}/${VERSION}")
    SET(FULL_VER "${PATH_ID}/${VERSION}")
  ELSE()
    SET(PATH_ID "custom" CACHE INTERNAL "Custom Path ID" FORCE)
    #set to $HOME/opt/<libname>
    #set module name accordingly
    SET(DEFAULT_PREFIX "${DEFAULT_LOCATION}/${PROJECT_NAME}/${PATH_ID}")
    SET(MODULE_NAME "${PROJECT_NAME}/${PATH_ID}")
    SET(FULL_VER "${PATH_ID}")
  ENDIF(LONG_PREFIX)
ELSE()
  SET(PATH_ID "custom" CACHE INTERNAL "Custom Path ID" FORCE)
  #set to user defined CMAKE_INSTALL_PREFIX
  #set module name accordingly
  SET(DEFAULT_PREFIX "${CMAKE_INSTALL_PREFIX}")
  SET(MODULE_NAME "${PROJECT_NAME}/${PATH_ID}")
  SET(FULL_VER "${PATH_ID}")
ENDIF()


#SET the resulting DEFAULT_PREFIX to CMAKE_INSTALL_PREFIX
SET(CMAKE_INSTALL_PREFIX "${DEFAULT_PREFIX}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)


#Set the module .version file name to be used later in PostConfig.cmake
SET(VERSION_PATH "${PROJECT_NAME}/${PATH_ID}")


