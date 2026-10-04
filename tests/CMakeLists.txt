# Regression tests: short fixed-step runs compared against stored per-step stats
# (tests/regression/reference/*.csv). Run with `ctest -L regression` from the build tree;
# needs gmsh on PATH (or GMSH=...) to generate the meshes on first use.
find_program(DOPAMINE_GMSH gmsh)
set(_regression_dir "${CMAKE_CURRENT_SOURCE_DIR}/regression")
set(_regression_work "${CMAKE_BINARY_DIR}/regression")

foreach(_case tgv channel rb duct dynsmag graded pdrop hydro hydrowall inflow inflow_jarrin)
  foreach(_np 1 2)
    add_test(NAME regression_${_case}_np${_np}
             COMMAND "${_regression_dir}/check_case.sh" ${_case} ${_np} "${_regression_work}")
    set(_env "DOPAMINE_BIN=$<TARGET_FILE:dopamine>" "MPIRUN=${MPIEXEC_EXECUTABLE}")
    if(DOPAMINE_GMSH)
      list(APPEND _env "GMSH=${DOPAMINE_GMSH}")
    endif()
    set_tests_properties(regression_${_case}_np${_np} PROPERTIES
      ENVIRONMENT "${_env}"
      LABELS regression
      PROCESSORS ${_np}
      TIMEOUT 900)
  endforeach()
endforeach()

# Unit test of the divergence-free synthetic eddy kernel (no mesh, serial).
add_executable(test_sem_dfsem unit/test_sem_dfsem.f90)
target_link_options(test_sem_dfsem PRIVATE -Wl,--no-as-needed)
target_link_libraries(test_sem_dfsem PRIVATE dopamine_lib MPI::MPI_Fortran MPI::MPI_C dl)
add_test(NAME unit_sem_dfsem COMMAND test_sem_dfsem)
set_tests_properties(unit_sem_dfsem PROPERTIES LABELS unit TIMEOUT 300)
