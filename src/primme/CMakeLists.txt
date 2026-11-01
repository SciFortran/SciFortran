# PRIMME 3.2.3 C core, vendored from https://github.com/SciFortran/primme
# Source layout follows the upstream src/Makefile.

SET(PRIMME_SRC
  src/eigs/auxiliary_eigs.c
  src/eigs/auxiliary_eigs_normal.c
  src/eigs/convergence.c
  src/eigs/correction.c
  src/eigs/factorize.c
  src/eigs/init.c
  src/eigs/inner_solve.c
  src/eigs/main_iter.c
  src/eigs/ortho.c
  src/eigs/primme_c.c
  src/eigs/primme_f77.c
  src/eigs/primme_interface.c
  src/eigs/restart.c
  src/eigs/solve_projection.c
  src/eigs/update_projection.c
  src/eigs/update_W.c
  src/linalg/auxiliary.c
  src/linalg/blaslapack.c
  src/linalg/magma_wrapper.c
  src/linalg/cublas_wrapper.c
  src/linalg/memman.c
  src/linalg/wtime.c
  src/svds/primme_svds_c.c
  src/svds/primme_svds_f77.c
  src/svds/primme_svds_interface.c
)

SET(PRIMME_HEADERS
  include/primme.h
  include/primme_eigs.h
  include/primme_eigs_f77.h
  include/primme_eigs_f90.inc
  include/primme_f77.h
  include/primme_f90.inc
  include/primme_svds.h
  include/primme_svds_f77.h
  include/primme_svds_f90.inc
)

ADD_LIBRARY(PRIMMELIB OBJECT ${PRIMME_SRC})
SET_TARGET_PROPERTIES(PRIMMELIB PROPERTIES POSITION_INDEPENDENT_CODE ON)
TARGET_INCLUDE_DIRECTORIES(PRIMMELIB PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/src/include
  ${CMAKE_CURRENT_SOURCE_DIR}/src/eigs
  ${CMAKE_CURRENT_SOURCE_DIR}/src/linalg
  ${CMAKE_CURRENT_SOURCE_DIR}/src/svds
)
TARGET_COMPILE_DEFINITIONS(PRIMMELIB PRIVATE F77UNDERSCORE)

# Keep the public C and Fortran include files available to SciFortran users.
INSTALL(FILES ${PRIMME_HEADERS} DESTINATION ${LIB_TARGET_INC})
INSTALL(FILES COPYING.txt DESTINATION ${LIB_TARGET_ETC}/licenses/primme)
