# Download dependencies by using FetchContent_Declare
# Use FetchContent_MakeAvailable only in those code parts where the dependency is actually needed

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)
set(FETCHCONTENT_BASE_DIR ${CMAKE_BINARY_DIR}/Thirdparty)

FetchContent_Declare(LBFGSpp
  GIT_REPOSITORY https://github.com/yixuan/LBFGSpp.git
  GIT_TAG master
  GIT_SHALLOW ON
  GIT_PROGRESS ON)
FetchContent_MakeAvailable(LBFGSpp)
