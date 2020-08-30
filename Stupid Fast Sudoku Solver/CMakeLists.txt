# CMakeList.txt : CMake project for Stupid Fast Sudoku Solver, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

find_package(fmt)

include_directories(
	/
	"${fmt_DIR}/../../../include"
)

set(CMakeTarget_HEADERS
)
set(CMakeTarget_SOURCES
	"Stupid Fast Sudoku Solver.cpp"
)

# Add source to this project's executable.
add_executable (CMakeTarget ${CMakeTarget_HEADERS} ${CMakeTarget_SOURCES})
set_property(TARGET CMakeTarget PROPERTY CXX_STANDARD 20)

target_link_libraries(CMakeTarget
	PRIVATE fmt::fmt-header-only)

# TODO: Add tests and install targets if needed.
