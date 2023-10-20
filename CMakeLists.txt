cmake_minimum_required(VERSION 3.8...3.12 FATAL_ERROR)

# Extra policies (not covering all introduced after 3.12)
foreach(policy
    CMP0083 # CMake 3.14
    CMP0127
    # Add more here
    )
  if(POLICY ${policy})
    cmake_policy(SET ${policy} NEW)
  endif()
endforeach()

include(CMakeDependentOption)

project(ICP
  LANGUAGES
  CXX)

set(ICP_MAJOR_VERSION 1)
set(ICP_MINOR_VERSION 0)
# Allow custom cmake overrides - not present right now
set(icp_cmake_dir ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
set(CMAKE_MODULE_PATH ${icp_cmake_dir} ${CMAKE_MODULE_PATH})

find_package(Git QUIET)

# Use ctest for dashboard submissions
include(CTest)
include(VTKIcpDetermineVersion) # TODO: Considering doing similarly

find_package(VTK 9.3 REQUIRED)

option(BUILD_SHARED_LIBS "Build shared libraries." ON)

if(VTK_WRAP_PYTHON)
  set(ICP_WRAP_PYTHON ON)
  option(BUILD_PYTHON_WRAPPERS "Build python wrappers" ON)

  get_property(vtk_soabi GLOBAL PROPERTY _vtk_python_soabi)
  cmake_dependent_option(VTK_WINDOWS_PYTHON_DEBUGGABLE "Append `_d` to Python module names" OFF
    "WIN32;VTK_WRAP_PYTHON" OFF)
  mark_as_advanced(VTK_WINDOWS_PYTHON_DEBUGGABLE)
else()
  unset(BUILD_PYTHON_WRAPPERS CACHE)
endif()

# Defines CMAKE_INSTALL_LIBDIR etc
include(GNUInstallDirs)

option(ICP_RELOCATABLE_INSTALL "Do not embed hard-coded paths into the install" ON)
mark_as_advanced(ICP_RELOCATABLE_INSTALL)

set(ICP_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
set(CMAKE_CXX_VISIBILITY_PRESET "hidden")
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)


# Set up our directory structure for output libraries and binaries
# (Note: these are the build locations, not the install locations)
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${ICP_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
endif()

if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  if(UNIX)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${ICP_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
  else()
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${ICP_BINARY_DIR}${CMAKE_INSTALL_BINDIR}")
  endif()
endif()

if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${ICP_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
endif()

set(icp_requested_modules
  "VTK::ICP"
)

option(ICP_BUILD_ALL_MODULES "Build all modules by default" OFF)
mark_as_advanced(ICP_BUILD_ALL_MODULES)

set(icp_source_directories
  "${CMAKE_CURRENT_SOURCE_DIR}/Source"
)

vtk_module_find_modules(icp_module_files ${icp_source_directories})

# Consider introducing KITS (right now none are found).
vtk_module_find_kits(icp_kit_files ${icp_source_directories})

vtk_module_scan(
  MODULE_FILES        "${CMAKE_CURRENT_SOURCE_DIR}/Source/vtk.module"
  KIT_FILES           ${icp_kit_files} # empty right now
  REQUEST_MODULES     ${icp_requested_modules}
  REJECT_MODULES      ${icp_rejected_modules}
  PROVIDES_MODULES    icp_modules
  PROVIDES_KITS       icp_kits
  WANT_BY_DEFAULT     "${ICP_BUILD_ALL_MODULES}"
  ENABLE_TESTS        "${BUILD_TESTING}")

if (NOT DEFINED icp_cmake_destination)
  set(icp_cmake_destination
    "${CMAKE_INSTALL_LIBDIR}/cmake/icp${icp_version_suffix}")
endif()

if (VTK_WRAP_PYTHON)
  if (NOT DEFINED ICP_PYTHON_SITE_PACKAGES_SUFFIX)
    vtk_module_python_default_destination(ICP_PYTHON_SITE_PACKAGES_SUFFIX)
  endif()
endif()

option(ICP_VERSIONED_INSTALL "Install artifacts with versioned names." ON)
mark_as_advanced(ICP_VERSIONED_INSTALL)
set(vtk_version_suffix "")
if (ICP_VERSIONED_INSTALL)
  set(vtk_version_suffix "-${ICP_MAJOR_VERSION}.${ICP_MINOR_VERSION}")
endif ()

set(ICP_CUSTOM_LIBRARY_SUFFIX "<DEFAULT>"
  CACHE STRING "Custom library file name suffix (defaults to the version number)")
mark_as_advanced(ICP_CUSTOM_LIBRARY_SUFFIX)
if (ICP_CUSTOM_LIBRARY_SUFFIX STREQUAL "<DEFAULT>")
  if (ICP_VERSIONED_INSTALL)
    set(ICP_CUSTOM_LIBRARY_SUFFIX "${ICP_MAJOR_VERSION}.${ICP_MINOR_VERSION}")
  else ()
    set(ICP_CUSTOM_LIBRARY_SUFFIX "")
  endif ()
endif ()

vtk_module_build(
  MODULES               ${icp_modules}
  KITS                  ${icp_kits}
  INSTALL_EXPORT        "${icp_install_export}"
  ARCHIVE_DESTINATION   "lib"
  HEADERS_DESTINATION   "${CMAKE_INSTALL_INCLUDEDIR}/icp${icp_version_suffix}"
  CMAKE_DESTINATION     "${icp_cmake_destination}"
  LICENSE_DESTINATION   "${CMAKE_INSTALL_LIBDIR}/icp"
  HIERARCHY_DESTINATION "${CMAKE_INSTALL_LIBDIR}/vtk/hierarchy/${CMAKE_PROJECT_NAME}"
  LIBRARY_NAME_SUFFIX   "${ICP_CUSTOM_LIBRARY_SUFFIX}"
  VERSION               "${ICP_VERSION}"
  SOVERSION             "1" # If the library is small and should work as a remote module to VTK
  TEST_DATA_TARGET      IcpData
  INSTALL_HEADERS       "${ICP_INSTALL_SDK}"
  USE_EXTERNAL          "${ICP_USE_EXTERNAL}"
  TEST_INPUT_DATA_DIRECTORY   "${CMAKE_CURRENT_SOURCE_DIR}/Testing"
  TEST_OUTPUT_DATA_DIRECTORY  "${CMAKE_CURRENT_BINARY_DIR}/ExternalData/Testing")

if(BUILD_PYTHON_WRAPPERS)
  find_package(PythonInterp ${VTK_PYTHON_VERSION} QUIET)

  vtk_module_wrap_python(
    MODULES         ${icp_modules}
    INSTALL_EXPORT  ICPPython
    PYTHON_PACKAGE  "icpmodules" # Directories under site-packages, vtk/modules
    LIBRARY_DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    MODULE_DESTINATION  "${ICP_PYTHON_SITE_PACKAGES_SUFFIX}"
    CMAKE_DESTINATION   "${icp_cmake_destination}"
    INSTALL_HEADERS     "${ICP_INSTALL_SDK}"
    SOABI               "${vtk_soabi}"
    #USE_DEBUG_SUFFIX    "${VTK_WINDOWS_PYTHON_DEBUGGABLE}"
    #UTILITY_TARGET      "VTK::vtkbuild"
    WRAPPED_MODULES icp_python_wrapped_modules
    BUILD_STATIC    OFF) # TODO: Consider supporting static builds

  message(${icp_python_wrapped_modules})
  add_subdirectory(Wrapping/Python)

  export(
    EXPORT    ICPPython
    NAMESPACE VTK::
    FILE      "${CMAKE_BINARY_DIR}/${icp_cmake_destination}/ICPPython-targets.cmake")
  install(
    EXPORT      ICPPython
    NAMESPACE   VTK::
    FILE        ICPPython-targets.cmake
    DESTINATION "${icp_cmake_destination}"
    COMPONENT   "development")

endif()
