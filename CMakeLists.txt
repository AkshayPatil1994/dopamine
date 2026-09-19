cmake_minimum_required(VERSION 3.22)
cmake_policy(SET CMP0135 NEW)   # ExternalProject URL timestamps = extraction time
project(dopamine VERSION 0.2.0 LANGUAGES Fortran CXX C)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_CXX_STANDARD     17)
set(CMAKE_C_STANDARD       11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# CMake's own NVHPC/PGI compiler module defaults CMAKE_Fortran_FLAGS_RELEASE to "-fast -O3".
# -fast bundles -Mvect, which forces O3-equivalent vectorization even on files given an
# explicit lower -O via per-source COMPILE_OPTIONS (needed as an nvfortran ICE workaround --
# see mod_advection.f90 below). Drop the CMake default in favor of just -O3, matching what's
# already set explicitly per-compiler further down.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE STRING "" FORCE)
endif()

# ============================================================================
# User-facing options
# ============================================================================
option(DOPAMINE_ENABLE_OPENMP     "Enable OpenMP threading in Fortran loops"  ON)
option(DOPAMINE_ENABLE_PETSC      "Enable PETSc pressure backend when available" ON)
option(DOPAMINE_FETCH_PETSC       "Download and build PETSc automatically"    ON)
option(DOPAMINE_PETSC_WITH_HYPRE  "Build PETSc with hypre (enables BoomerAMG preconditioner, recommended for large meshes)" ON)
option(DOPAMINE_PETSC_WITH_CUDA   "Build PETSc with CUDA GPU support"         OFF)
option(DOPAMINE_PETSC_USE_GPU_AWARE_MPI "Request GPU-aware MPI in PETSc CUDA builds" ON)
option(DOPAMINE_ENABLE_OPENACC   "Enable OpenACC GPU offload (requires NVHPC/nvfortran)" OFF)
set(DOPAMINE_OPENACC_GPU_ARCH "" CACHE STRING "Optional -gpu= target for OpenACC (e.g. cc86); empty lets nvfortran pick defaults")
option(DOPAMINE_FETCH_NETCDF      "Auto-build NetCDF C/Fortran when missing"  ON)
option(DOPAMINE_FETCH_CGNS        "Auto-build CGNS (Fortran+HDF5) when missing" ON)
option(DOPAMINE_BUILD_TESTS       "Build unit tests"                          ON)
option(DOPAMINE_ENABLE_LTO        "Enable interprocedural/link-time optimization (GNU: -flto, Intel: -ipo)" OFF)

# External build directory
set(DOPAMINE_EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/external" CACHE PATH
  "Directory where external libraries are downloaded and built")

# Guard against empty cache values that would collapse to absolute root paths
# like /petsc-install and break imported target include directories.
if(NOT DOPAMINE_EXTERNAL_DIR)
  set(DOPAMINE_EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/external" CACHE PATH
      "Directory where external libraries are downloaded and built" FORCE)
endif()

set(DOPAMINE_PETSC_VERSION "3.25.0" CACHE STRING "PETSc version tag (v<version>)")
set(DOPAMINE_PETSC_ARCH "arch-dopamine-opt" CACHE STRING "PETSc ARCH build directory name")
set(DOPAMINE_PETSC_CUDA_ARCH "" CACHE STRING "Optional PETSc CUDA architecture (for example: 80 or sm_80)")
set(DOPAMINE_PETSC_CUDAFLAGS "" CACHE STRING "Optional CUDAFLAGS passed through to PETSc configure")
set(DOPAMINE_PETSC_CUDA_HOST_COMPILER "" CACHE FILEPATH "Optional host compiler passed to nvcc via --compiler-bindir")

set(DOPAMINE_PETSC_INSTALL_DIR "${DOPAMINE_EXTERNAL_DIR}/petsc-install" CACHE PATH
    "Path to PETSc install prefix (containing include/ and lib/)")
set(DOPAMINE_PETSC_ARCH_DIR "${DOPAMINE_EXTERNAL_DIR}/petsc-src/${DOPAMINE_PETSC_ARCH}" CACHE PATH
    "Path to PETSc ARCH directory (for additional headers/libs)")

# Repair stale or malformed cache entries from earlier configurations.
if(DOPAMINE_PETSC_INSTALL_DIR STREQUAL "/petsc-install" OR
   DOPAMINE_PETSC_INSTALL_DIR STREQUAL "petsc-install")
  set(DOPAMINE_PETSC_INSTALL_DIR "${DOPAMINE_EXTERNAL_DIR}/petsc-install" CACHE PATH
    "Path to PETSc install prefix (containing include/ and lib/)" FORCE)
endif()
if(DOPAMINE_PETSC_ARCH_DIR STREQUAL "/petsc-src/${DOPAMINE_PETSC_ARCH}" OR
   DOPAMINE_PETSC_ARCH_DIR STREQUAL "petsc-src/${DOPAMINE_PETSC_ARCH}")
  set(DOPAMINE_PETSC_ARCH_DIR "${DOPAMINE_EXTERNAL_DIR}/petsc-src/${DOPAMINE_PETSC_ARCH}" CACHE PATH
    "Path to PETSc ARCH directory (for additional headers/libs)" FORCE)
endif()

include(ExternalProject)
include(GNUInstallDirs)

function(dopamine_configure_petsc_imported_target target_name)
  set(_petsc_pc "${DOPAMINE_PETSC_INSTALL_DIR}/lib/pkgconfig/petsc.pc")
  if(EXISTS "${_petsc_pc}")
    file(STRINGS "${_petsc_pc}" _petsc_libs_private_line REGEX "^Libs.private:")
    if(_petsc_libs_private_line)
      string(REGEX REPLACE "^Libs.private:[ ]*" "" _petsc_libs_private "${_petsc_libs_private_line}")
      separate_arguments(_petsc_libs_private NATIVE_COMMAND "${_petsc_libs_private}")
      if(_petsc_libs_private)
        # Collect -L search paths and -l/-absolute-path library entries separately.
        # IMPORTANT: do NOT use INTERFACE_LINK_OPTIONS for -L flags — CMake emits
        # INTERFACE_LINK_OPTIONS after INTERFACE_LINK_LIBRARIES on the linker command
        # line, so -lFOO would precede the -L that lets ld find libFOO.a.
        # Use INTERFACE_LINK_DIRECTORIES instead (always placed before -l entries).
        # For -l libraries known to live in the PETSc install dir, prefer absolute
        # paths so the search path is irrelevant.
        set(_link_dirs "")
        set(_link_libs "")
        set(_link_opts "")
        foreach(_item IN LISTS _petsc_libs_private)
          if(_item MATCHES "^-L(.+)$")
            list(APPEND _link_dirs "${CMAKE_MATCH_1}")
          elseif(_item MATCHES "^-l(.+)$")
            set(_libname "${CMAKE_MATCH_1}")
            set(_abs "${DOPAMINE_PETSC_INSTALL_DIR}/lib/lib${_libname}.a")
            if(EXISTS "${_abs}")
              # Resolve to absolute path — no search path needed and ordering
              # is unambiguous for static archives.
              list(APPEND _link_libs "${_abs}")
            else()
              list(APPEND _link_libs "${_item}")
            endif()
          elseif(_item MATCHES "^-Wl,")
            list(APPEND _link_opts "${_item}")
          else()
            # Absolute archives / shared objects passed directly.
            list(APPEND _link_libs "${_item}")
          endif()
        endforeach()
        if(_link_dirs)
          set_property(TARGET ${target_name} APPEND PROPERTY
            INTERFACE_LINK_DIRECTORIES "${_link_dirs}")
        endif()
        if(_link_libs)
          set_property(TARGET ${target_name} APPEND PROPERTY
            INTERFACE_LINK_LIBRARIES "${_link_libs}")
        endif()
        if(_link_opts)
          set_property(TARGET ${target_name} APPEND PROPERTY
            INTERFACE_LINK_OPTIONS "${_link_opts}")
        endif()
      endif()
    endif()
  else()
    # First configure with ExternalProject fetch: petsc.pc does not exist yet.
    # PETSc is built static with --download-f2cblaslapack=1, so provide a
    # deterministic fallback to satisfy BLAS/LAPACK symbols at final link.
    set_property(TARGET ${target_name} APPEND PROPERTY
      INTERFACE_LINK_LIBRARIES
      "${DOPAMINE_PETSC_INSTALL_DIR}/lib/libf2clapack.a"
      "${DOPAMINE_PETSC_INSTALL_DIR}/lib/libf2cblas.a"
      m)
  endif()
endfunction()

# ============================================================================
# MPI — required; must be provided by the system (not auto-built)
# ============================================================================
find_package(MPI REQUIRED COMPONENTS Fortran CXX C)

if(DOPAMINE_ENABLE_PETSC)
  # dopamine never uses PETSc's native Fortran module interface (petscsys/petscksp) --
  # the whole integration goes through pressure_backend_petsc.c via ISO_C_BINDING -- so
  # generating them is unneeded work. Skip it under NVHPC specifically: nvfortran 23.3
  # fails PETSc's configure-time probe for F2018 pointer-initialization-in-declaration,
  # which only that binding-generation step requires.
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(DOPAMINE_PETSC_FORTRAN_BINDINGS 0)
  else()
    set(DOPAMINE_PETSC_FORTRAN_BINDINGS 1)
  endif()
  set(DOPAMINE_PETSC_CONFIGURE_ARGS
    --prefix=${DOPAMINE_PETSC_INSTALL_DIR}
    --PETSC_ARCH=${DOPAMINE_PETSC_ARCH}
    --with-debugging=0
    --with-shared-libraries=0
    --with-fortran-bindings=${DOPAMINE_PETSC_FORTRAN_BINDINGS}
    --with-x=0
    --download-f2cblaslapack=1
    --with-mpi=1
    --with-cc=${MPI_C_COMPILER}
    --with-fc=${MPI_Fortran_COMPILER})

  if(DOPAMINE_PETSC_WITH_HYPRE)
    list(APPEND DOPAMINE_PETSC_CONFIGURE_ARGS
      --with-cxx=${MPI_CXX_COMPILER}
      --download-hypre=1)
  else()
    list(APPEND DOPAMINE_PETSC_CONFIGURE_ARGS --with-cxx=0)
  endif()

  if(DOPAMINE_PETSC_WITH_CUDA)
    set(DOPAMINE_PETSC_CUDA_FLAGS_ARG "${DOPAMINE_PETSC_CUDAFLAGS}")
    find_program(DOPAMINE_PETSC_CUDAC nvcc)
    if(NOT DOPAMINE_PETSC_CUDAC)
      message(FATAL_ERROR "DOPAMINE_PETSC_WITH_CUDA=ON requires nvcc in PATH")
    endif()
    list(APPEND DOPAMINE_PETSC_CONFIGURE_ARGS --with-cuda=1)
    list(APPEND DOPAMINE_PETSC_CONFIGURE_ARGS --with-cudac=${DOPAMINE_PETSC_CUDAC})
    if(DOPAMINE_PETSC_CUDA_HOST_COMPILER)
      if(DOPAMINE_PETSC_CUDA_FLAGS_ARG)
        string(APPEND DOPAMINE_PETSC_CUDA_FLAGS_ARG " ")
      endif()
      string(APPEND DOPAMINE_PETSC_CUDA_FLAGS_ARG "--compiler-bindir=${DOPAMINE_PETSC_CUDA_HOST_COMPILER}")
    endif()
    if(DOPAMINE_PETSC_CUDA_FLAGS_ARG)
      list(APPEND DOPAMINE_PETSC_CONFIGURE_ARGS CUDAFLAGS=${DOPAMINE_PETSC_CUDA_FLAGS_ARG})
    endif()
    if(DOPAMINE_PETSC_CUDA_ARCH)
      list(APPEND DOPAMINE_PETSC_CONFIGURE_ARGS --with-cuda-arch=${DOPAMINE_PETSC_CUDA_ARCH})
    endif()
    if(DOPAMINE_PETSC_USE_GPU_AWARE_MPI)
      list(APPEND DOPAMINE_PETSC_CONFIGURE_ARGS --with-cuda-aware-mpi=1)
    endif()
  else()
    # Without this, PETSc's own configure auto-probes for nvcc on PATH and silently
    # enables CUDA even when DOPAMINE_PETSC_WITH_CUDA=OFF, on any machine that has it.
    list(APPEND DOPAMINE_PETSC_CONFIGURE_ARGS --with-cuda=0)
  endif()

  if(DOPAMINE_FETCH_PETSC)
    # Ensure include/lib paths referenced by the imported target exist already
    # during CMake generation on first configure.
    file(MAKE_DIRECTORY "${DOPAMINE_PETSC_INSTALL_DIR}/include" "${DOPAMINE_PETSC_INSTALL_DIR}/lib")

    # Deliberately no BUILD_BYPRODUCTS here: with the Makefiles generator, CMake's
    # `clean` target sweeps declared byproducts of custom commands project-wide,
    # so listing libpetsc.a/libHYPRE.a there means a plain `make clean` deletes an
    # externally-built library that took minutes to compile, while leaving the
    # ExternalProject step stamps in an inconsistent state (they aren't declared
    # byproducts, so they survive, and a subsequent build then sees "already done"
    # stamps pointing at missing files). Ordering is already guaranteed via
    # add_dependencies(PETSc::PETSc ext_petsc) below, which the Makefiles generator
    # doesn't need BUILD_BYPRODUCTS for (Ninja would, if this project ever adopts it).
    ExternalProject_Add(ext_petsc
      PREFIX               "${DOPAMINE_EXTERNAL_DIR}/petsc"
      GIT_REPOSITORY       "https://gitlab.com/petsc/petsc.git"
      GIT_TAG              "v${DOPAMINE_PETSC_VERSION}"
      GIT_SHALLOW          ON
      DOWNLOAD_NO_PROGRESS ON
      BUILD_IN_SOURCE      ON
      PATCH_COMMAND        ${CMAKE_COMMAND} -E remove -f "<SOURCE_DIR>/CLAUDE.md"
      CONFIGURE_COMMAND
        <SOURCE_DIR>/configure
          ${DOPAMINE_PETSC_CONFIGURE_ARGS}
      BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} PETSC_DIR=<SOURCE_DIR> PETSC_ARCH=${DOPAMINE_PETSC_ARCH} all
      INSTALL_COMMAND
        ${CMAKE_MAKE_PROGRAM} PETSC_DIR=<SOURCE_DIR> PETSC_ARCH=${DOPAMINE_PETSC_ARCH} install
      LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
      UPDATE_DISCONNECTED ON
    )

    add_library(PETSc::PETSc UNKNOWN IMPORTED GLOBAL)
    set_target_properties(PETSc::PETSc PROPERTIES
      IMPORTED_LOCATION "${DOPAMINE_PETSC_INSTALL_DIR}/lib/libpetsc.a"
      INTERFACE_INCLUDE_DIRECTORIES "${DOPAMINE_PETSC_INSTALL_DIR}/include")
    dopamine_configure_petsc_imported_target(PETSc::PETSc)
    # Declare HYPRE as part of PETSc::PETSc's interface using its installed path.
    # This is necessary on the first configure (before petsc.pc exists) and avoids
    # the unreliable find_library glob-path search that caused the link error.
    if(DOPAMINE_PETSC_WITH_HYPRE)
      set_property(TARGET PETSc::PETSc APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES "${DOPAMINE_PETSC_INSTALL_DIR}/lib/libHYPRE.a")
    endif()
    add_dependencies(PETSc::PETSc ext_petsc)
    set(DOPAMINE_HAVE_PETSC TRUE)
    message(STATUS "dopamine: PETSc will be fetched and built in ${DOPAMINE_EXTERNAL_DIR}/petsc")
  else()
    find_path(DOPAMINE_PETSC_INCLUDE_DIR petscksp.h
      PATHS
        "${DOPAMINE_PETSC_INSTALL_DIR}/include"
        "${DOPAMINE_PETSC_ARCH_DIR}/include"
      NO_DEFAULT_PATH)
    find_library(DOPAMINE_PETSC_LIBRARY NAMES petsc
      PATHS
        "${DOPAMINE_PETSC_INSTALL_DIR}/lib"
        "${DOPAMINE_PETSC_ARCH_DIR}/lib"
      NO_DEFAULT_PATH)

    if(DOPAMINE_PETSC_INCLUDE_DIR AND DOPAMINE_PETSC_LIBRARY)
      add_library(PETSc::PETSc UNKNOWN IMPORTED GLOBAL)
      set_target_properties(PETSc::PETSc PROPERTIES
        IMPORTED_LOCATION "${DOPAMINE_PETSC_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${DOPAMINE_PETSC_INCLUDE_DIR}")
      dopamine_configure_petsc_imported_target(PETSc::PETSc)
      set(DOPAMINE_HAVE_PETSC TRUE)
      message(STATUS "dopamine: PETSc enabled via ${DOPAMINE_PETSC_LIBRARY}")
    else()
      set(DOPAMINE_HAVE_PETSC FALSE)
      message(WARNING "dopamine: PETSc requested but not found; pressure solver backend unavailable")
    endif()
  endif()
else()
  set(DOPAMINE_HAVE_PETSC FALSE)
endif()

# ============================================================================
# External prefix: everything lands in external/install/
# ============================================================================
set(EXT_INSTALL "${DOPAMINE_EXTERNAL_DIR}/install")
file(MAKE_DIRECTORY "${EXT_INSTALL}/include" "${EXT_INSTALL}/lib")
list(PREPEND CMAKE_PREFIX_PATH "${EXT_INSTALL}")

# ============================================================================
# NetCDF (required for sampling/output): system package first, auto-build fallback
# ============================================================================
find_package(PkgConfig QUIET)
set(DOPAMINE_HAVE_SYSTEM_NETCDF FALSE)
set(DOPAMINE_NETCDF_TARGET "")
set(DOPAMINE_NETCDF_INCLUDE_DIRS "")

# Fortran .mod files are compiler-specific (format/version), so a system netcdf-fortran
# package -- built with GNU Fortran -- can't be read by a different Fortran compiler
# (observed: "Corrupt or Old Module file" under nvfortran). Only trust it under GNU;
# any other compiler falls through to the auto-build path below, which compiles
# netcdf-fortran with whatever Fortran compiler this configure is actually using.
if(PkgConfig_FOUND AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  pkg_check_modules(NETCDF_FORTRAN QUIET IMPORTED_TARGET netcdf-fortran)
  if(TARGET PkgConfig::NETCDF_FORTRAN)
    set(DOPAMINE_HAVE_SYSTEM_NETCDF TRUE)
    set(DOPAMINE_NETCDF_TARGET PkgConfig::NETCDF_FORTRAN)
    if(DEFINED NETCDF_FORTRAN_INCLUDE_DIRS)
      list(APPEND DOPAMINE_NETCDF_INCLUDE_DIRS ${NETCDF_FORTRAN_INCLUDE_DIRS})
    endif()
    if(DEFINED NETCDF_FORTRAN_INCLUDEDIRS)
      list(APPEND DOPAMINE_NETCDF_INCLUDE_DIRS ${NETCDF_FORTRAN_INCLUDEDIRS})
    endif()
    find_path(NETCDF_FORTRAN_MOD_DIR netcdf.mod
      HINTS
        ${DOPAMINE_NETCDF_INCLUDE_DIRS}
        /usr/include
        /usr/local/include)
    if(NETCDF_FORTRAN_MOD_DIR)
      list(APPEND DOPAMINE_NETCDF_INCLUDE_DIRS "${NETCDF_FORTRAN_MOD_DIR}")
    endif()
    message(STATUS "dopamine: using system netcdf-fortran (pkg-config)")
  endif()
endif()

if(NOT DOPAMINE_HAVE_SYSTEM_NETCDF)
  if(DOPAMINE_FETCH_NETCDF)
    message(STATUS "dopamine: netcdf-fortran not found; auto-building NetCDF C/Fortran in ${EXT_INSTALL}")

    # The output writer creates NetCDF-4 files and XDMF references HDF datasets.
    # Require HDF5 for fallback builds to preserve that behavior.
    find_package(HDF5 QUIET COMPONENTS C)
    if(NOT HDF5_FOUND)
      message(FATAL_ERROR
        "NetCDF auto-build requires HDF5 development libraries (HDF5 C). "
        "Install HDF5 dev packages or provide a system netcdf-fortran installation.")
    endif()

    ExternalProject_Add(ext_netcdf_c
      PREFIX                  "${DOPAMINE_EXTERNAL_DIR}/netcdf-c"
      URL                     "https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.tar.gz"
      DOWNLOAD_NO_PROGRESS    ON
      CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${EXT_INSTALL}
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_C_COMPILER=${MPI_C_COMPILER}
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DBUILD_SHARED_LIBS=ON
        -DBUILD_TESTING=OFF
        -DENABLE_NETCDF_4=ON
        -DENABLE_DAP=OFF
        -DENABLE_BYTERANGE=OFF
        -DENABLE_NCZARR=OFF
      BUILD_BYPRODUCTS
        "${EXT_INSTALL}/lib/libnetcdf${CMAKE_SHARED_LIBRARY_SUFFIX}"
      LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
      UPDATE_DISCONNECTED ON
    )

    ExternalProject_Add(ext_netcdf_fortran
      PREFIX                  "${DOPAMINE_EXTERNAL_DIR}/netcdf-fortran"
      URL                     "https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.1.tar.gz"
      DOWNLOAD_NO_PROGRESS    ON
      CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${EXT_INSTALL}
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_C_COMPILER=${MPI_C_COMPILER}
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DBUILD_SHARED_LIBS=ON
        -DBUILD_TESTING=OFF
        -DENABLE_TESTS=OFF
        -DENABLE_EXAMPLES=OFF
        -DNETCDF_C_LIBRARY=${EXT_INSTALL}/lib/libnetcdf${CMAKE_SHARED_LIBRARY_SUFFIX}
        -DNETCDF_C_INCLUDE_DIR=${EXT_INSTALL}/include
      BUILD_BYPRODUCTS
        "${EXT_INSTALL}/lib/libnetcdff${CMAKE_SHARED_LIBRARY_SUFFIX}"
        "${EXT_INSTALL}/include/netcdf.mod"
      DEPENDS ext_netcdf_c
      LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
      UPDATE_DISCONNECTED ON
    )

    add_library(NetCDF::Fortran UNKNOWN IMPORTED GLOBAL)
    set_target_properties(NetCDF::Fortran PROPERTIES
      IMPORTED_LOCATION "${EXT_INSTALL}/lib/libnetcdff${CMAKE_SHARED_LIBRARY_SUFFIX}"
      INTERFACE_INCLUDE_DIRECTORIES "${EXT_INSTALL}/include")
    add_dependencies(NetCDF::Fortran ext_netcdf_fortran)

    set(DOPAMINE_NETCDF_TARGET NetCDF::Fortran)
    list(APPEND DOPAMINE_NETCDF_INCLUDE_DIRS "${EXT_INSTALL}/include")
  else()
    message(FATAL_ERROR
      "NetCDF Fortran not found (pkg-config: netcdf-fortran). "
      "Install netcdf-fortran or enable DOPAMINE_FETCH_NETCDF=ON.")
  endif()
endif()

# ============================================================================
# CGNS (optional, for mesh_type='cgns'): system package first, auto-build fallback
# ============================================================================
set(DOPAMINE_HAVE_SYSTEM_CGNS FALSE)
set(DOPAMINE_CGNS_TARGET "")
set(DOPAMINE_CGNS_INCLUDE_DIRS "")

# Same Fortran .mod compatibility caveat as NetCDF above: only trust a system
# cgns-fortran package under GNU Fortran.
if(PkgConfig_FOUND AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  pkg_check_modules(CGNS_PKG QUIET IMPORTED_TARGET cgns)
  if(TARGET PkgConfig::CGNS_PKG)
    set(DOPAMINE_HAVE_SYSTEM_CGNS TRUE)
    set(DOPAMINE_CGNS_TARGET PkgConfig::CGNS_PKG)
    if(DEFINED CGNS_PKG_INCLUDE_DIRS)
      list(APPEND DOPAMINE_CGNS_INCLUDE_DIRS ${CGNS_PKG_INCLUDE_DIRS})
    endif()
    find_path(CGNS_FORTRAN_MOD_DIR cgns.mod
      HINTS ${DOPAMINE_CGNS_INCLUDE_DIRS} /usr/include /usr/local/include)
    if(CGNS_FORTRAN_MOD_DIR)
      list(APPEND DOPAMINE_CGNS_INCLUDE_DIRS "${CGNS_FORTRAN_MOD_DIR}")
    endif()
    message(STATUS "dopamine: using system CGNS (pkg-config)")
  endif()
endif()

if(NOT DOPAMINE_HAVE_SYSTEM_CGNS)
  if(DOPAMINE_FETCH_CGNS)
    message(STATUS "dopamine: CGNS not found; auto-building CGNS (Fortran+HDF5) in ${EXT_INSTALL}")

    # CGNS's HDF5 backend needs the same HDF5 dev libraries NetCDF-4 already requires above.
    find_package(HDF5 QUIET COMPONENTS C)
    if(NOT HDF5_FOUND)
      message(FATAL_ERROR
        "CGNS auto-build requires HDF5 development libraries (HDF5 C). "
        "Install HDF5 dev packages or provide a system CGNS installation.")
    endif()

    ExternalProject_Add(ext_cgns
      PREFIX                  "${DOPAMINE_EXTERNAL_DIR}/cgns"
      URL                     "https://github.com/CGNS/CGNS/archive/refs/tags/v4.4.0.tar.gz"
      DOWNLOAD_NO_PROGRESS    ON
      CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${EXT_INSTALL}
        -DCMAKE_PREFIX_PATH=${EXT_INSTALL}
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_C_COMPILER=${MPI_C_COMPILER}
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DCGNS_ENABLE_FORTRAN=ON
        -DCGNS_ENABLE_HDF5=ON
        -DCGNS_BUILD_SHARED=ON
        -DCGNS_ENABLE_TESTS=OFF
        -DCGNS_BUILD_CGNSTOOLS=OFF
      BUILD_BYPRODUCTS
        "${EXT_INSTALL}/lib/libcgns${CMAKE_SHARED_LIBRARY_SUFFIX}"
        "${EXT_INSTALL}/include/cgns.mod"
      LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
      UPDATE_DISCONNECTED ON
    )

    add_library(CGNS::Fortran UNKNOWN IMPORTED GLOBAL)
    set_target_properties(CGNS::Fortran PROPERTIES
      IMPORTED_LOCATION "${EXT_INSTALL}/lib/libcgns${CMAKE_SHARED_LIBRARY_SUFFIX}"
      INTERFACE_INCLUDE_DIRECTORIES "${EXT_INSTALL}/include")
    add_dependencies(CGNS::Fortran ext_cgns)

    set(DOPAMINE_CGNS_TARGET CGNS::Fortran)
    list(APPEND DOPAMINE_CGNS_INCLUDE_DIRS "${EXT_INSTALL}/include")
  else()
    message(FATAL_ERROR
      "CGNS not found (pkg-config: cgns). "
      "Install a CGNS Fortran package or enable DOPAMINE_FETCH_CGNS=ON.")
  endif()
endif()

if(FALSE)  # ParMETIS removed — Morton (Z-order) SFC partitioner is the sole partitioner
  ExternalProject_Add(ext_gklib
    PREFIX                  "${DOPAMINE_EXTERNAL_DIR}/gklib"
    GIT_REPOSITORY          "https://github.com/KarypisLab/GKlib.git"
    GIT_TAG                 "master"
    GIT_SHALLOW             ON
    DOWNLOAD_NO_PROGRESS    ON
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${EXT_INSTALL}
      -DCMAKE_INSTALL_LIBDIR=lib
      -DCMAKE_BUILD_TYPE=Release
      -DCMAKE_C_COMPILER=${MPI_C_COMPILER}
      -DSHARED=OFF
    BUILD_BYPRODUCTS
      "${EXT_INSTALL}/lib/libGKlib.a"
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
    UPDATE_DISCONNECTED ON
    PATCH_COMMAND ${CMAKE_COMMAND} -E remove -f "<SOURCE_DIR>/CLAUDE.md"
  )

  add_library(GKLIB::GKLIB STATIC IMPORTED GLOBAL)
  set_target_properties(GKLIB::GKLIB PROPERTIES
    IMPORTED_LOCATION "${EXT_INSTALL}/lib/libGKlib.a"
    INTERFACE_INCLUDE_DIRECTORIES "${EXT_INSTALL}/include")
  add_dependencies(GKLIB::GKLIB ext_gklib)

  ExternalProject_Add(ext_metis
    PREFIX                  "${DOPAMINE_EXTERNAL_DIR}/metis"
    GIT_REPOSITORY          "https://github.com/KarypisLab/METIS.git"
    GIT_TAG                 "master"
    GIT_SHALLOW             ON
    BUILD_IN_SOURCE         ON
    DOWNLOAD_NO_PROGRESS    ON
    CONFIGURE_COMMAND       ${CMAKE_MAKE_PROGRAM} config prefix=${EXT_INSTALL} gklib_path=${EXT_INSTALL} cc=${MPI_C_COMPILER}
    BUILD_COMMAND           ${CMAKE_MAKE_PROGRAM} -C build
    INSTALL_COMMAND         ${CMAKE_MAKE_PROGRAM} -C build install
    BUILD_BYPRODUCTS
      "${EXT_INSTALL}/lib/libmetis.a"
    DEPENDS ext_gklib
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
    UPDATE_DISCONNECTED ON
    PATCH_COMMAND ${CMAKE_COMMAND} -E remove -f "<SOURCE_DIR>/CLAUDE.md"
  )

  # ParMETIS ships its own CMakeLists.txt (>= 4.0.3).
  # We patch in shared=OFF and enable the Fortran interface.
  ExternalProject_Add(ext_parmetis
    PREFIX                  "${DOPAMINE_EXTERNAL_DIR}/parmetis"
    GIT_REPOSITORY          "https://github.com/KarypisLab/ParMETIS.git"
    GIT_TAG                 "${PARMETIS_GIT_TAG}"
    GIT_SHALLOW             ON
    DOWNLOAD_NO_PROGRESS    ON
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${EXT_INSTALL}
      -DCMAKE_BUILD_TYPE=Release
      -DCMAKE_C_COMPILER=${MPI_C_COMPILER}
      -DGKLIB_PATH=${EXT_INSTALL}
      -DMETIS_PATH=${EXT_INSTALL}
      -DSHARED=OFF
    BUILD_BYPRODUCTS
      "${EXT_INSTALL}/lib/libparmetis.a"
    DEPENDS ext_metis
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
    UPDATE_DISCONNECTED ON
    PATCH_COMMAND ${CMAKE_COMMAND} -E remove -f "<SOURCE_DIR>/CLAUDE.md"
  )

endif()

# ============================================================================
# 7. OpenMP (optional, system-provided)
#    AppleClang does not ship libomp; install via: brew install libomp
# ============================================================================
if(DOPAMINE_ENABLE_OPENMP)
  # On macOS with AppleClang, point CMake at the Homebrew libomp if present
  if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang")
    find_program(BREW_CMD brew)
    if(BREW_CMD)
      execute_process(
        COMMAND ${BREW_CMD} --prefix libomp
        OUTPUT_VARIABLE LIBOMP_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
      if(LIBOMP_PREFIX)
        set(OpenMP_CXX_FLAGS         "-Xpreprocessor -fopenmp")
        set(OpenMP_CXX_LIB_NAMES     "omp")
        set(OpenMP_omp_LIBRARY       "${LIBOMP_PREFIX}/lib/libomp.dylib")
        set(OpenMP_CXX_INCLUDE_DIR   "${LIBOMP_PREFIX}/include")
        message(STATUS "macOS: using Homebrew libomp at ${LIBOMP_PREFIX}")
      else()
        message(STATUS "macOS: libomp not found via Homebrew — disabling OpenMP for C++")
        message(STATUS "       Install with: brew install libomp")
      endif()
    endif()
  endif()

  find_package(OpenMP COMPONENTS Fortran)
  find_package(OpenMP OPTIONAL_COMPONENTS CXX)   # non-fatal if C++ OpenMP missing
endif()

if(DOPAMINE_ENABLE_OPENACC)
  if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    message(FATAL_ERROR "DOPAMINE_ENABLE_OPENACC=ON requires nvfortran (NVHPC); "
                         "got Fortran compiler ID '${CMAKE_Fortran_COMPILER_ID}'")
  endif()
  find_package(OpenACC REQUIRED COMPONENTS Fortran)
endif()

add_library(petsc_pressure_backend OBJECT
  src/solver/pressure_backend_petsc.c)

target_link_libraries(petsc_pressure_backend PRIVATE MPI::MPI_C)
if(DOPAMINE_HAVE_PETSC)
  target_link_libraries(petsc_pressure_backend PRIVATE PETSc::PETSc)
  target_compile_definitions(petsc_pressure_backend PRIVATE FVM_HAVE_PETSC)
endif()

# ============================================================================
# Fortran source list (dependency order for module compilation)
# ============================================================================
set(DOPAMINE_FORTRAN_SOURCES
  src/mod_kinds.f90
  src/mod_parameters.f90
  src/mod_scalar.f90
  src/io/mod_warn_log.f90
  src/stats/mod_sampling_types.f90
  src/mesh/mod_mesh.f90
  src/mesh/mod_gmsh_reader.f90
  src/mesh/mod_foam_reader.f90
  src/mesh/mod_icem_reader.f90
  src/mesh/mod_cgns_reader.F90
  src/mesh/mod_partitioner.f90
  src/mesh/mod_mesh_io.f90
  src/mesh/mod_mesh_validate.f90
  src/mesh/mod_nonconforming.f90
  src/mesh/mod_streaming_decompose.f90
  src/io/mod_checkpoint_io.f90
  src/numerics/mod_nonortho.f90
  src/numerics/mod_gradient.f90
  src/numerics/mod_advection.f90
  src/numerics/mod_viscous.f90
  src/sgs/mod_vreman.f90
  src/sgs/mod_dynamic_smagorinsky.f90
  src/sgs/mod_wall_model.f90
  src/abl/mod_abl_types.f90
  src/abl/mod_abl_profiles.f90
  src/abl/mod_sem_inflow.f90
  src/abl/mod_abl_weather.f90
  src/abl/mod_abl_wall.f90
  src/solver/mod_pressure_solver.f90
  src/bc/mod_bc.f90
  src/bc/mod_bc_setup.f90
  src/solver/mod_fv_solver.f90
  src/init/mod_initconditions.f90
  src/stats/mod_sampling.f90
  src/stats/mod_stats.f90
  src/io/mod_input_parser.f90
  src/io/mod_runtime_env.f90
  src/io/mod_output_writer.f90
)

# ============================================================================
# Library target
# ============================================================================
add_library(dopamine_lib STATIC
  ${DOPAMINE_FORTRAN_SOURCES}
  $<TARGET_OBJECTS:petsc_pressure_backend>)

set_target_properties(dopamine_lib PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod")

target_include_directories(dopamine_lib PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_BINARY_DIR}/mod"
  "${EXT_INSTALL}/include"
  ${DOPAMINE_NETCDF_INCLUDE_DIRS}
  ${DOPAMINE_CGNS_INCLUDE_DIRS})

target_link_libraries(dopamine_lib PUBLIC
  MPI::MPI_Fortran
  MPI::MPI_C)


if(DOPAMINE_HAVE_PETSC)
  target_link_libraries(dopamine_lib PUBLIC PETSc::PETSc)
  if(DOPAMINE_PETSC_WITH_CUDA)
    # PETSc's CUDA/cupm device layer is C++ (e.g. mark_dcontext.o); the Fortran link
    # driver for the executables below doesn't pull in libstdc++ on its own for that.
    target_link_libraries(dopamine_lib PUBLIC stdc++)

    # petsc.pc's Libs.private lists -lcublas/-lcusparse/-lcusolver/-lcurand but, under
    # NVHPC, no matching -L for the directory those actually live in -- NVHPC ships them
    # under math_libs/, not next to nvfortran/nvcc -- so the bare -l names leave the
    # final link with undefined cublasDtrsv_v2 etc. Derive that directory from the
    # Fortran compiler's own install root and add it explicitly; a no-op on non-NVHPC
    # CUDA installs, where these libs already sit in a default search path (see the
    # "GPU (CUDA) Build" section, unaffected by this).
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
      get_filename_component(_nvhpc_bin_dir "${CMAKE_Fortran_COMPILER}" DIRECTORY)   # .../compilers/bin
      get_filename_component(_nvhpc_ver_dir "${_nvhpc_bin_dir}" DIRECTORY)           # .../compilers
      get_filename_component(_nvhpc_ver_dir "${_nvhpc_ver_dir}" DIRECTORY)           # .../26.5
      set(_nvhpc_math_libs_dir "${_nvhpc_ver_dir}/math_libs/lib64")
      if(EXISTS "${_nvhpc_math_libs_dir}")
        target_link_directories(dopamine_lib PUBLIC "${_nvhpc_math_libs_dir}")
      else()
        message(WARNING "dopamine: expected NVHPC math_libs at ${_nvhpc_math_libs_dir} "
                         "but it doesn't exist; PETSc's CUDA link (cublas/cusparse/cusolver) "
                         "may fail with undefined references. Check your NVHPC install layout.")
      endif()
      # PETSc's own petsc.pc Libs.private also carries -lcudart/-lcuda/-lnvidia-ml (the CUDA
      # runtime + driver, needed by HYPRE's device kernels and CUB's RDC-launch stubs), but
      # dopamine_configure_petsc_imported_target() above only parses that file if it already
      # exists at CMake *configure* time -- true when PETSc happens to already be built, but
      # not on a from-scratch DOPAMINE_FETCH_PETSC=ON build, where petsc.pc is only written
      # by ext_petsc's own build step, long after this configure runs. That leaves the final
      # link silently missing -lcudart on a first build, failing with "undefined reference to
      # __cudaRegisterFatBinary" etc. from HYPRE's/CUB's .o files. Add it unconditionally here
      # instead of depending on Libs.private timing.
      set(_nvhpc_cuda_libs_dir "${_nvhpc_ver_dir}/cuda/lib64")
      if(EXISTS "${_nvhpc_cuda_libs_dir}")
        target_link_directories(dopamine_lib PUBLIC "${_nvhpc_cuda_libs_dir}")
        # libcuda.so/libnvidia-ml.so are driver libraries; NVHPC ships link-time stubs for
        # them under cuda/lib64/stubs (the real .so.1 comes from the installed GPU driver
        # at runtime), same as petsc.pc's own cudalib recipe.
        if(EXISTS "${_nvhpc_cuda_libs_dir}/stubs")
          target_link_directories(dopamine_lib PUBLIC "${_nvhpc_cuda_libs_dir}/stubs")
        endif()
        target_link_libraries(dopamine_lib PUBLIC cudart cuda nvidia-ml)
      else()
        message(WARNING "dopamine: expected NVHPC CUDA libs at ${_nvhpc_cuda_libs_dir} "
                         "but it doesn't exist; PETSc's CUDA link (cudart/cuda) may fail "
                         "with undefined references. Check your NVHPC install layout.")
      endif()
    endif()
  endif()
endif()

target_link_libraries(dopamine_lib PUBLIC ${DOPAMINE_NETCDF_TARGET})
target_link_libraries(dopamine_lib PUBLIC ${DOPAMINE_CGNS_TARGET})

if(OpenMP_Fortran_FOUND AND DOPAMINE_ENABLE_OPENMP)
  target_link_libraries(dopamine_lib PUBLIC OpenMP::OpenMP_Fortran)
endif()
if(DOPAMINE_ENABLE_OPENACC AND OpenACC_Fortran_FOUND)
  target_link_libraries(dopamine_lib PUBLIC OpenACC::OpenACC_Fortran)
endif()
target_compile_options(dopamine_lib PRIVATE
  $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:GNU>>:
    -O3 -march=native -Wall -Wextra
    -ffree-line-length-none -fbacktrace
    $<$<CONFIG:Debug>:-fcheck=all -g>>
  $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:Intel,IntelLLVM>>:
    -O3 -xHost -warn all -traceback
    $<$<CONFIG:Debug>:-check bounds -g>>
  $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:NVHPC>>:
    -O3 -Minfo=accel
    $<$<BOOL:${DOPAMINE_ENABLE_OPENMP}>:-mp>
    $<$<AND:$<BOOL:${DOPAMINE_ENABLE_OPENACC}>,$<BOOL:${DOPAMINE_OPENACC_GPU_ARCH}>>:-gpu=${DOPAMINE_OPENACC_GPU_ARCH}>>
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:
    -O3 -march=native -Wall>
  $<$<AND:$<COMPILE_LANGUAGE:C>,$<C_COMPILER_ID:GNU,Clang,AppleClang>>:
    -O3 -march=native -Wall>)

# nvfortran 23.3 hits an internal compiler error (ICE) in resolve_reconstruction_context
# at -O2/-O3 ("NVFORTRAN-S-0000-Internal compiler error. bad cc 11") on a trivial array-
# broadcast assignment -- a compiler bug, not a code issue (confirmed: -O1 compiles clean,
# gfortran/ifx have no issue at any level). Drop just this file to -O1 under NVHPC rather
# than losing optimization project-wide; append after the target-level -O3 so it wins.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  # -Mvect (implied by the target-level -O3 above) stays sticky even after a later -O1 on
  # the same command line -- nvfortran-Info-Switch confirms "-Mvect forces -O3" -- so the
  # ICE workaround needs -Mnovect explicitly, not just a lower -O.
  set_source_files_properties(src/numerics/mod_advection.f90 PROPERTIES
    COMPILE_OPTIONS "-O1;-Mnovect")

  # Separate, persistent nvfortran ICE (confirmed present in both 23.3 and 26.5 -- not
  # fixed by upgrading): "Lowering Error: bad ast optype in expression" whenever this file
  # is compiled with -mp, independent of optimization level. Triggered by ordinary code
  # (a plain do-loop copying derived-type fields into arrays), not anything unusual here.
  # This file is the offline dopamine-decompose preprocessing tool's writer-thread loop,
  # not part of the solver hot path, so disabling OpenMP for just this file under NVHPC
  # only costs threading in mesh preprocessing, not simulation performance.
  set_source_files_properties(src/mesh/mod_streaming_decompose.f90 PROPERTIES
    COMPILE_OPTIONS "-Mnoopenmp")
endif()

# GPU-offloaded files: dual !$omp (CPU, default) / !$acc (GPU, -DDOPAMINE_ENABLE_OPENACC=ON)
# code paths on the same loops, selected via #ifdef _OPENACC (defined automatically when compiling
# with -acc; never defined otherwise) so the two directive sets never stack on one loop. Needs
# preprocessing under every compiler, not just NVHPC -- gfortran also has to see and strip the
# #ifdef/#else/#endif lines correctly, even though _OPENACC is never defined for it.
set_source_files_properties(
  src/sgs/mod_vreman.f90
  src/sgs/mod_dynamic_smagorinsky.f90
  src/numerics/mod_gradient.f90
  src/numerics/mod_viscous.f90
  src/numerics/mod_advection.f90
  src/solver/mod_fv_solver.f90
  src/main.f90
  PROPERTIES Fortran_PREPROCESS ON)

# Ensure external libraries are built before the Fortran sources
if(TARGET ext_netcdf_fortran)
  add_dependencies(dopamine_lib ext_netcdf_fortran)
endif()

# ============================================================================
# Interprocedural/link-time optimization (opt-in; untried, benchmark before
# defaulting on -- see the CPU performance audit notes)
# ============================================================================
if(DOPAMINE_ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT DOPAMINE_LTO_SUPPORTED OUTPUT DOPAMINE_LTO_ERROR LANGUAGES C Fortran)
  if(DOPAMINE_LTO_SUPPORTED)
    set_target_properties(dopamine_lib PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    message(STATUS "dopamine: LTO/IPO enabled for dopamine_lib")
  else()
    message(WARNING "DOPAMINE_ENABLE_LTO=ON but not supported by this toolchain: ${DOPAMINE_LTO_ERROR}")
  endif()
endif()

# ============================================================================
# Executable
# ============================================================================
add_executable(dopamine src/main.f90)
target_link_options(dopamine PRIVATE -Wl,--no-as-needed)
target_link_libraries(dopamine PRIVATE
  dopamine_lib
  MPI::MPI_Fortran
  MPI::MPI_C)

# Ensure libdl is linked for PETSc dynamic loading
target_link_libraries(dopamine PRIVATE dl)

# ============================================================================
# dopamine-decompose (offline mesh partitioner)
# ============================================================================
add_executable(dopamine-decompose src/tools/decompose_main.f90)
target_link_options(dopamine-decompose PRIVATE -Wl,--no-as-needed)
target_link_libraries(dopamine-decompose PRIVATE
  dopamine_lib
  MPI::MPI_Fortran
  MPI::MPI_C
  dl)

# ============================================================================
# dopamine-checkpoint-partition (checkpoint repartitioner)
# ============================================================================
add_executable(dopamine-checkpoint-partition src/tools/checkpoint_partition_main.f90)
target_link_options(dopamine-checkpoint-partition PRIVATE -Wl,--no-as-needed)
target_link_libraries(dopamine-checkpoint-partition PRIVATE
  dopamine_lib
  MPI::MPI_Fortran
  MPI::MPI_C
  dl)

if(DOPAMINE_ENABLE_LTO AND DOPAMINE_LTO_SUPPORTED)
  # LTO must also apply at the final link step (not just dopamine_lib's compile),
  # or the linker just uses each object's fallback native code without the
  # cross-TU merge pass that's the actual point of enabling this.
  set_target_properties(dopamine dopamine-decompose dopamine-checkpoint-partition
    PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# ============================================================================
# Tests
# ============================================================================
if(DOPAMINE_BUILD_TESTS)
  enable_testing()
  if(EXISTS "${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt")
    add_subdirectory(tests EXCLUDE_FROM_ALL)
  endif()
endif()

# ============================================================================
# cleanall — removes fetched/built external libraries (PETSc, hypre, ...) too.
# Plain `make clean` only removes dopamine's own build artifacts; rebuilding
# PETSc+hypre from scratch takes minutes, so that stays opt-in and explicit.
# ============================================================================
add_custom_target(cleanall
  COMMAND ${CMAKE_COMMAND} --build "${CMAKE_BINARY_DIR}" --target clean
  COMMAND ${CMAKE_COMMAND} -E rm -rf "${DOPAMINE_EXTERNAL_DIR}"
  COMMENT "Removing solver build artifacts and fetched external libraries (${DOPAMINE_EXTERNAL_DIR})"
  VERBATIM)

# ============================================================================
# Install
# ============================================================================
install(TARGETS dopamine dopamine-decompose dopamine-checkpoint-partition RUNTIME DESTINATION bin)
install(TARGETS dopamine_lib ARCHIVE DESTINATION lib)
install(DIRECTORY src/ DESTINATION include/dopamine
        FILES_MATCHING PATTERN "*.f90" PATTERN "*.hpp")

# ============================================================================
# Summary
# ============================================================================
set(DOPAMINE_SUMMARY_BOX_INNER_WIDTH 50)
string(REPEAT "═" ${DOPAMINE_SUMMARY_BOX_INNER_WIDTH} DOPAMINE_SUMMARY_BOX_RULE)
set(DOPAMINE_SUMMARY_TITLE " dopamine v${PROJECT_VERSION} configuration ")
string(LENGTH "${DOPAMINE_SUMMARY_TITLE}" DOPAMINE_SUMMARY_TITLE_LEN)
math(EXPR DOPAMINE_SUMMARY_TITLE_PAD "${DOPAMINE_SUMMARY_BOX_INNER_WIDTH} - ${DOPAMINE_SUMMARY_TITLE_LEN}")
if(DOPAMINE_SUMMARY_TITLE_PAD LESS 0)
  set(DOPAMINE_SUMMARY_TITLE_PAD 0)
endif()
string(REPEAT " " ${DOPAMINE_SUMMARY_TITLE_PAD} DOPAMINE_SUMMARY_TITLE_FILL)

message(STATUS "")
message(STATUS "╔${DOPAMINE_SUMMARY_BOX_RULE}╗")
message(STATUS "║${DOPAMINE_SUMMARY_TITLE}${DOPAMINE_SUMMARY_TITLE_FILL}║")
message(STATUS "╠${DOPAMINE_SUMMARY_BOX_RULE}╣")
message(STATUS "  Build type     : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Fortran        : ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")
message(STATUS "  C++            : ${CMAKE_CXX_COMPILER}")
message(STATUS "  MPI            : ${MPI_Fortran_COMPILER}")
message(STATUS "  ParMETIS       : fetch=${DOPAMINE_FETCH_PARMETIS}")
message(STATUS "  PETSc          : enabled=${DOPAMINE_ENABLE_PETSC}, fetch=${DOPAMINE_FETCH_PETSC}, available=${DOPAMINE_HAVE_PETSC}")
message(STATUS "  PETSc install  : ${DOPAMINE_PETSC_INSTALL_DIR}")
message(STATUS "  PETSc arch     : ${DOPAMINE_PETSC_ARCH}")
message(STATUS "  hypre/BoomerAMG: ${DOPAMINE_PETSC_WITH_HYPRE} (pressure solver default)")
message(STATUS "  PETSc CUDA     : ${DOPAMINE_PETSC_WITH_CUDA}")
message(STATUS "  GPU-aware MPI  : ${DOPAMINE_PETSC_USE_GPU_AWARE_MPI}")
message(STATUS "  CUDA arch      : ${DOPAMINE_PETSC_CUDA_ARCH}")
message(STATUS "  CUDAFLAGS      : ${DOPAMINE_PETSC_CUDAFLAGS}")
message(STATUS "  CUDA host cxx  : ${DOPAMINE_PETSC_CUDA_HOST_COMPILER}")
message(STATUS "  OpenMP         : ${OpenMP_Fortran_FOUND}")
message(STATUS "  External dir   : ${DOPAMINE_EXTERNAL_DIR}")
message(STATUS "╚${DOPAMINE_SUMMARY_BOX_RULE}╝")
message(STATUS "")
message(STATUS "First build: external libraries will be downloaded and compiled.")
message(STATUS "Subsequent builds: only dopamine sources recompile.")
message(STATUS "To skip auto-build: set -DDOPAMINE_FETCH_<LIB>=OFF and provide")
message(STATUS "  -DCMAKE_PREFIX_PATH=<path> pointing to installed libraries.")
message(STATUS "")
