#Set permission

#777 is potentially dangerous
SET(PERMISSION_777
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
  GROUP_WRITE GROUP_READ GROUP_EXECUTE
  WORLD_WRITE WORLD_READ WORLD_EXECUTE)

SET(PERMISSION_700
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
)


#Build the Library module for environment modules
FILE(REMOVE_RECURSE ${LIB_TMP_ETC}/modules)
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY ON)

#Create default .version file if following conditions are verified:
#1. master branch
#2. no debug
IF( (NOT "${BUILD_TYPE}" MATCHES "debug") AND (NOT "${PATH_ID}" MATCHES "custom") AND (GIT_BRANCH MATCHES "master"))
  SET(TMP_VER_MODULE_FILE ${LIB_TMP_ETC}/modules/${VERSION_PATH}/.version)
  CONFIGURE_FILE(${LIB_ENV}/version.in ${TMP_VER_MODULE_FILE}  @ONLY)
  MESSAGE(STATUS "${Red}Version file${ColourReset}: ${TMP_VER_MODULE_FILE}")
ELSE()
  if(DEFINED TMP_VER_MODULE_FILE)
    unset(TMP_VER_MODULE_FILE CACHE)
  endif() 
ENDIF()


#Build the library module for environment modules
SET(TMP_MODULE_NAME "${MODULE_NAME}" CACHE PATH "Prefix prepended to install directories")
SET(TMP_ENV_MODULE_FILE ${LIB_TMP_ETC}/modules/${TMP_MODULE_NAME})
CONFIGURE_FILE(${LIB_ENV}/module.in ${TMP_ENV_MODULE_FILE} @ONLY)
MESSAGE(STATUS "${Red}Module file${ColourReset}: ${MODULE_NAME}")


#Build the user CONFIG scripts (sourced in user shell config file, i.e. .bashrc)
SET(USER_CONFIG_SCRIPT ${PROJECT_NAME}_config_user.sh)
SET(TMP_CONFIGVARS_USER_FILE ${LIB_TMP_ETC}/${USER_CONFIG_SCRIPT})
CONFIGURE_FILE(${LIB_ETC}/${PROJECT_NAME}_config_user.sh.in ${TMP_CONFIGVARS_USER_FILE} @ONLY)

SET(GLOBAL_CONFIG_SCRIPT ${PROJECT_NAME}_config_global.sh)
SET(TMP_CONFIGVARS_GLOBAL_FILE ${LIB_TMP_ETC}/${GLOBAL_CONFIG_SCRIPT})
CONFIGURE_FILE(${LIB_ETC}/${PROJECT_NAME}_config_global.sh.in ${TMP_CONFIGVARS_GLOBAL_FILE} @ONLY)


FILE(WRITE  ${LIB_TMP_VER}  "${VERSION}\n")

#Install .mod files
INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${LIB_TARGET_INC})
#Install environment module file

INSTALL(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E rm -f ${LIB_VERSION_FILE})")

INSTALL(TARGETS ${PROJECT_NAME} DESTINATION ${LIB_TARGET_LIB})

INSTALL(DIRECTORY ${LIB_TMP_ETC}/ DESTINATION ${LIB_TARGET_ETC})

#Install config scripts:
INSTALL(FILES ${TMP_CONFIGVARS_USER_FILE} DESTINATION ${LIB_TARGET_BIN}/
  PERMISSIONS ${PERMISSION_700} SETUID)

INSTALL(FILES ${TMP_CONFIGVARS_GLOBAL_FILE} DESTINATION ${LIB_TARGET_BIN}/
  PERMISSIONS ${PERMISSION_700} SETUID)

#Install version file
INSTALL(FILES ${LIB_TMP_VER} DESTINATION ${LIB_TARGET_DIR} 
  PERMISSIONS ${PERMISSION_700} SETUID)

INSTALL(FILES ${LIB_TARGET_ETC}/${PROJECT_NAME}.pc DESTINATION $ENV{HOME}/.pkgconfig.d/
  PERMISSIONS ${PERMISSION_700} SETUID)

INSTALL(DIRECTORY ${LIB_TARGET_ETC}/modules/ DESTINATION $ENV{HOME}/.modules.d)

get_filename_component(BARE_MAKE_PROGRAM ${CMAKE_MAKE_PROGRAM} NAME)
MESSAGE( STATUS "${Red}Library version:${ColourReset} ${VERSION}")
MESSAGE( STATUS "${Red}Library will be installed in:${ColourReset} ${CMAKE_INSTALL_PREFIX}")
MESSAGE( STATUS "
>> ${Red}TO CONCLUDE INSTALLATION${ColourReset} <<
Compile with:
$ ${BARE_MAKE_PROGRAM}
Install with:
$ ${BARE_MAKE_PROGRAM} install

Uninstall with:
$ ${BARE_MAKE_PROGRAM} uninstall
")

INSTALL(CODE "MESSAGE(
\"
ADD LIBRARY TO YOUR SYSTEM: 
Pick ONE method below [or add it in your bash profile, e.g. ~/.bashrc]:
${Yellow}Method 1: use the provided ${PROJECT_NAME} environment module${ColourReset}:
   $ module use $HOME/.modules.d
   $ module load ${TMP_MODULE_NAME} 

${Yellow}Method 2: source the config script${ColourReset}:
   $ source ${LIB_TARGET_BIN}/scifor_config_user.sh

${Yellow}Method 3: use pkg-config with the provided ${PROJECT_NAME}.pc${ColourReset}:
   $ export PKG_CONFIG_PATH=${LIB_TARGET_ETC}/:$PKG_CONFIG_PATH
   $ pkg-config --cflags --libs ${PROJECT_NAME}

${Yellow}Method ADMIN: Add this line to the system shell configuration file, e.g. /etc/bash.bashrc${ColourReset}
   $ source ${LIB_TARGET_BIN}/scifor_config_global.sh
\")
")


SET(MODULE_FILE_NAME ${TMP_MODULE_NAME} CACHE PATH " ") 
SET(CONFIG_FILE_PATH ${LIB_TARGET_BIN} CACHE PATH " ")


# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_MODULE_PATH}/DistClean.cmake
)


# Uninstall target
if(NOT TARGET uninstall)
  CONFIGURE_FILE(
    "${CMAKE_MODULE_PATH}/ConfigUninstall.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
    "${PROJECT_NAME}/${TMP_MODULE_NAME}" "${PROJECT_NAME}.pc" )
  
ENDIF()
