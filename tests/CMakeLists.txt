# Regression tests: short fixed-step runs compared against stored per-step stats
# (tests/regression/reference/*.csv). Run with `ctest -L regression` from the build tree;
# needs gmsh on PATH (or GMSH=...) to generate the meshes on first use.
find_program(DOPAMINE_GMSH gmsh)
set(_regression_dir "${CMAKE_CURRENT_SOURCE_DIR}/regression")
set(_regression_work "${CMAKE_BINARY_DIR}/regression")

foreach(_case tgv channel rb)
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
