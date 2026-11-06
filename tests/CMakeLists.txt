# Regression tests: short fixed-step runs compared against stored per-step stats
# (tests/regression/reference/*.csv). Run with `ctest -L regression` from the build tree;
# needs gmsh on PATH (or GMSH=...) to generate the meshes on first use.
find_program(DOPAMINE_GMSH gmsh)
set(_regression_dir "${CMAKE_CURRENT_SOURCE_DIR}/regression")
set(_regression_work "${CMAKE_BINARY_DIR}/regression")

foreach(_case tgv tgv_skew tgv_skew_chk tgv_skew_exact default_skew hydro_exact duct_exact duct_exact_gamg channel channel_gamg channel_most rb duct dynsmag graded pdrop hydro hydrowall inflow inflow_jarrin inflow_backflow)
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

# AMD + adaptive wall-model filter, serial only (np 2 cannot be run on every host yet).
add_test(NAME regression_amd_np1
         COMMAND "${_regression_dir}/check_case.sh" amd 1 "${_regression_work}")
set(_env "DOPAMINE_BIN=$<TARGET_FILE:dopamine>" "MPIRUN=${MPIEXEC_EXECUTABLE}")
if(DOPAMINE_GMSH)
  list(APPEND _env "GMSH=${DOPAMINE_GMSH}")
endif()
set_tests_properties(regression_amd_np1 PROPERTIES ENVIRONMENT "${_env}" LABELS regression PROCESSORS 1 TIMEOUT 900)

# EXACT projection: the written mean pressure must stay pinned by the Dirichlet outlet (no gauge drift).
foreach(_np 1 2)
  add_test(NAME regression_outlet_gauge_np${_np}
           COMMAND "${_regression_dir}/check_mean_p.sh" outlet_gauge ${_np} "${_regression_work}" 0.5)
  set(_env "DOPAMINE_BIN=$<TARGET_FILE:dopamine>" "MPIRUN=${MPIEXEC_EXECUTABLE}")
  set_tests_properties(regression_outlet_gauge_np${_np} PROPERTIES
    ENVIRONMENT "${_env}" LABELS regression PROCESSORS ${_np} TIMEOUT 900 SKIP_RETURN_CODE 77)
endforeach()

# Unit test of the divergence-free synthetic eddy kernel (no mesh, serial).
add_executable(test_sem_dfsem unit/test_sem_dfsem.f90)
target_link_options(test_sem_dfsem PRIVATE -Wl,--no-as-needed)
target_link_libraries(test_sem_dfsem PRIVATE dopamine_lib MPI::MPI_Fortran MPI::MPI_C dl)
add_test(NAME unit_sem_dfsem COMMAND test_sem_dfsem)
set_tests_properties(unit_sem_dfsem PROPERTIES LABELS unit TIMEOUT 300)

add_executable(test_amd unit/test_amd.f90)
target_link_options(test_amd PRIVATE -Wl,--no-as-needed)
target_link_libraries(test_amd PRIVATE dopamine_lib MPI::MPI_Fortran MPI::MPI_C dl)
add_test(NAME unit_amd COMMAND test_amd)
set_tests_properties(unit_amd PROPERTIES LABELS unit TIMEOUT 60)

add_executable(test_inflow_db unit/test_inflow_db.f90)
target_link_options(test_inflow_db PRIVATE -Wl,--no-as-needed)
target_link_libraries(test_inflow_db PRIVATE dopamine_lib MPI::MPI_Fortran MPI::MPI_C dl)
add_test(NAME unit_inflow_db COMMAND test_inflow_db)
set_tests_properties(unit_inflow_db PROPERTIES LABELS unit TIMEOUT 120)

add_executable(test_inflow_gen unit/test_inflow_gen.f90)
target_link_options(test_inflow_gen PRIVATE -Wl,--no-as-needed)
target_link_libraries(test_inflow_gen PRIVATE dopamine_lib MPI::MPI_Fortran MPI::MPI_C dl)
add_test(NAME unit_inflow_gen COMMAND test_inflow_gen)
set_tests_properties(unit_inflow_gen PROPERTIES LABELS unit TIMEOUT 300)

add_executable(test_inflow_planes unit/test_inflow_planes.f90)
target_link_options(test_inflow_planes PRIVATE -Wl,--no-as-needed)
target_link_libraries(test_inflow_planes PRIVATE dopamine_lib MPI::MPI_Fortran MPI::MPI_C dl)
add_test(NAME unit_inflow_planes COMMAND test_inflow_planes)
set_tests_properties(unit_inflow_planes PROPERTIES LABELS unit TIMEOUT 600)

# End-to-end inflow-plane chain (needs gmsh): export grid -> generate -> replay, stop checks,
# restart, precursor recorder. `ctest -L inflow`.
add_test(NAME inflow_chain COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/inflow/chain.sh" "${CMAKE_BINARY_DIR}/inflow_chain")
set_tests_properties(inflow_chain PROPERTIES
  ENVIRONMENT "DOPAMINE_BIN=$<TARGET_FILE:dopamine>;DOPAMINE_GEN=$<TARGET_FILE:dopamine-inflow-gen>;MPIRUN=${MPIEXEC_EXECUTABLE}"
  LABELS inflow PROCESSORS 2 TIMEOUT 1200)

# Polygon face geometry on concave faces (serial, no mesh file).
add_executable(test_mesh_geometry unit/test_mesh_geometry.f90)
target_link_options(test_mesh_geometry PRIVATE -Wl,--no-as-needed)
target_link_libraries(test_mesh_geometry PRIVATE dopamine_lib MPI::MPI_Fortran MPI::MPI_C dl)
add_test(NAME unit_mesh_geometry COMMAND test_mesh_geometry)
set_tests_properties(unit_mesh_geometry PROPERTIES LABELS unit TIMEOUT 60)

# gmsh reader must refuse a mesh whose boundary faces lack a physical surface group.
add_test(NAME unit_gmsh_open_cell
         COMMAND test_mesh_geometry "${CMAKE_CURRENT_SOURCE_DIR}/unit/mesh/gmsh_open_tet.msh")
set_tests_properties(unit_gmsh_open_cell PROPERTIES LABELS unit TIMEOUT 60
  PASS_REGULAR_EXPRESSION "unlabelled boundary faces")

# cos(theta) floor of the over-relaxed k-vector (serial, no mesh).
add_executable(test_nonortho_limiter unit/test_nonortho_limiter.f90)
target_link_options(test_nonortho_limiter PRIVATE -Wl,--no-as-needed)
target_link_libraries(test_nonortho_limiter PRIVATE dopamine_lib MPI::MPI_Fortran MPI::MPI_C dl)
add_test(NAME unit_nonortho_limiter COMMAND test_nonortho_limiter)
set_tests_properties(unit_nonortho_limiter PROPERTIES LABELS unit TIMEOUT 60)

# Convective outlet: order in space of the face derivative and in time of its RK3W integration (serial, no mesh).
add_executable(test_outlet_bc unit/test_outlet_bc.f90)
target_link_options(test_outlet_bc PRIVATE -Wl,--no-as-needed)
target_link_libraries(test_outlet_bc PRIVATE dopamine_lib MPI::MPI_Fortran MPI::MPI_C dl)
add_test(NAME unit_outlet_bc COMMAND test_outlet_bc)
set_tests_properties(unit_outlet_bc PROPERTIES LABELS unit TIMEOUT 60)

# Exact cell projection on a skewed box (mesh from mk_box_mesh.py, no gmsh needed).
add_executable(test_cell_projection unit/test_cell_projection.f90)
target_link_options(test_cell_projection PRIVATE -Wl,--no-as-needed)
target_link_libraries(test_cell_projection PRIVATE dopamine_lib MPI::MPI_Fortran MPI::MPI_C dl)
add_test(NAME unit_cell_projection
         COMMAND sh -c "python3 ${CMAKE_CURRENT_SOURCE_DIR}/regression/mk_box_mesh.py 6 0.03 cp_box.msh && $<TARGET_FILE:test_cell_projection> cp_box.msh"
         WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
set_tests_properties(unit_cell_projection PROPERTIES LABELS unit TIMEOUT 120)
