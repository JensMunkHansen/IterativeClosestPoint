# NOTE: We do not support static Python wrappers like VTK.
#       Static wrappers means linking libraries directly into
#       the python executable.

# This file is copied and adapted from VTK/Wrapping/Python/CMakeLists.txt
# and there may be some rendundant material.

set(_icpmodules_all)
set(_icp_python_imports)

foreach (_icp_python_module IN LISTS icp_python_wrapped_modules)
  get_property(_icp_python_library_name
    TARGET    "${_icp_python_module}"
    PROPERTY  "INTERFACE_vtk_module_library_name")
  string(APPEND _icpmodules_all
    "    \'${_icp_python_library_name}\',\n")
  string(APPEND _icp_python_imports
    "from .${_icp_python_library_name} import *\n")
endforeach ()

set(_icp_python_files)
foreach (_icp_module IN LISTS icp_modules)
  _vtk_module_get_module_property("${_icp_module}"
    PROPERTY "python_modules"
    VARIABLE _icp_module_python_files)
  list(APPEND _icp_python_files
    ${_icp_module_python_files})
endforeach()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/icpmodules/all.py.in"
  "${CMAKE_BINARY_DIR}/${ICP_PYTHON_SITE_PACKAGES_SUFFIX}/icpmodules/all.py"
  @ONLY)
list(APPEND _icp_python_files
  "${CMAKE_BINARY_DIR}/${ICP_PYTHON_SITE_PACKAGES_SUFFIX}/icpmodules/all.py")

if (BUILD_SHARED_LIBS)
  install(
    FILES       "${CMAKE_BINARY_DIR}/${ICP_PYTHON_SITE_PACKAGES_SUFFIX}/icpmodules/all.py"
    DESTINATION "${ICP_PYTHON_SITE_PACKAGES_SUFFIX}/icpmodules"
    COMPONENT   "python")
endif ()

# Python files to add
set(python_files
  icp.py
  icpmodules/__init__.py.in
  icpmodules/util/__init__.py
  icpmodules/util/io.py
)

foreach (python_file IN LISTS python_files)
  # add submodules that don't begin with "_" to _icpmodules_all
  if (python_file MATCHES "^icpmodules\/([A-Za-z][A-Za-z0-9_]*)(\/__init__)?[.]py")
    string(APPEND _icpmodules_all "    \'${CMAKE_MATCH_1}\',\n")
  endif ()
endforeach ()

set(ICP_PYTHON_CONFIGURATION_DEBUG "Debug" CACHE STRING "")
set(ICP_PYTHON_CONFIGURATION_RELEASE "Release" CACHE STRING "")

set(python_copied_modules)
foreach (python_file IN LISTS python_files)
  set(output_python_file
    "${CMAKE_BINARY_DIR}/${ICP_PYTHON_SITE_PACKAGES_SUFFIX}/${python_file}")
  if (python_file MATCHES "\\.in$")
    string(REPLACE ".in" "" output_python_file "${output_python_file}")
    configure_file(
      "${python_file}"
      "${output_python_file}"
      @ONLY)
  else ()
    add_custom_command(
      OUTPUT  "${output_python_file}"
      DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${python_file}"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
              "${CMAKE_CURRENT_SOURCE_DIR}/${python_file}"
              "${output_python_file}"
      COMMENT "Copying ${python_file} to the binary directory")
  endif ()
  if (BUILD_SHARED_LIBS)
    get_filename_component(python_file_directory "${python_file}" DIRECTORY)
    install(
      FILES       "${output_python_file}"
      DESTINATION "${ICP_PYTHON_SITE_PACKAGES_SUFFIX}/${python_file_directory}"
      COMPONENT   "python")
  endif ()
  list(APPEND python_copied_modules
    "${output_python_file}")
endforeach ()

set(icp_build_paths_module
  "${CMAKE_BINARY_DIR}/${ICP_PYTHON_SITE_PACKAGES_SUFFIX}/icpmodules/_build_paths.py")

# Figure out the relevant DLL paths to add on Windows
if (WIN32 AND ICP_DLL_PATHS)
  string(REPLACE ";" "',\n    '" icp_dll_paths_python "${ICP_DLL_PATHS}")
  file(WRITE "${icp_build_paths_module}"
    "paths = [\n    '${icp_dll_paths_python}',\n]\n")
  # Add the file to the zip only if required in the install tree.
  if (NOT ICP_RELOCATABLE_INSTALL AND NOT ICP_UNIFIED_INSTALL_TREE)
    list(APPEND _icp_python_files
      "${icp_build_paths_module}")
  endif ()
elseif (EXISTS "${icp_build_paths_module}")
  file(REMOVE "${icp_build_paths_module}")
endif ()

list(APPEND _icp_python_files
  ${python_copied_modules})

set(_icp_python_zip)

# TODO: Populate if we decide to support STATIC wrappers (we don't)
add_custom_target(icp_python_copy ALL
  DEPENDS
    ${python_copied_modules}
    ${_icp_python_zip})

