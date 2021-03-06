include_directories(
  ${Boost_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  )

install(FILES 	
  solver.h
  linearOperatorSolver.h
  cgSolver.h
  nlcgSolver.h
  lbfgsSolver.h
  sbSolver.h
  sbcSolver.h
  cgCallback.h	
  cgPreconditioner.h
  lwSolver.h
  lbfgsSolver.h
  lsqrSolver.h
  gpSolver.h
  gpBbSolver.h
  eigenTester.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

IF(ARMADILLO_FOUND)
  add_subdirectory(cpu)
ENDIF(ARMADILLO_FOUND)

IF( CUDA_FOUND)
  add_subdirectory(gpu)
ENDIF (CUDA_FOUND)
