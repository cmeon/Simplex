# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/cmeon/SimplexLP/eigen

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cmeon/SimplexLP/lib

# Include any dependencies generated for this target.
include test/CMakeFiles/array_for_matrix_1.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/array_for_matrix_1.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/array_for_matrix_1.dir/flags.make

test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o: test/CMakeFiles/array_for_matrix_1.dir/flags.make
test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o: /home/cmeon/SimplexLP/eigen/test/array_for_matrix.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cmeon/SimplexLP/lib/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o -c /home/cmeon/SimplexLP/eigen/test/array_for_matrix.cpp

test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.i"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cmeon/SimplexLP/eigen/test/array_for_matrix.cpp > CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.i

test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.s"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cmeon/SimplexLP/eigen/test/array_for_matrix.cpp -o CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.s

test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o.requires:
.PHONY : test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o.requires

test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o.provides: test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/array_for_matrix_1.dir/build.make test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o.provides.build
.PHONY : test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o.provides

test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o.provides.build: test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o

# Object files for target array_for_matrix_1
array_for_matrix_1_OBJECTS = \
"CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o"

# External object files for target array_for_matrix_1
array_for_matrix_1_EXTERNAL_OBJECTS =

test/array_for_matrix_1: test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o
test/array_for_matrix_1: test/CMakeFiles/array_for_matrix_1.dir/build.make
test/array_for_matrix_1: test/CMakeFiles/array_for_matrix_1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable array_for_matrix_1"
	cd /home/cmeon/SimplexLP/lib/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/array_for_matrix_1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/array_for_matrix_1.dir/build: test/array_for_matrix_1
.PHONY : test/CMakeFiles/array_for_matrix_1.dir/build

test/CMakeFiles/array_for_matrix_1.dir/requires: test/CMakeFiles/array_for_matrix_1.dir/array_for_matrix.cpp.o.requires
.PHONY : test/CMakeFiles/array_for_matrix_1.dir/requires

test/CMakeFiles/array_for_matrix_1.dir/clean:
	cd /home/cmeon/SimplexLP/lib/test && $(CMAKE_COMMAND) -P CMakeFiles/array_for_matrix_1.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/array_for_matrix_1.dir/clean

test/CMakeFiles/array_for_matrix_1.dir/depend:
	cd /home/cmeon/SimplexLP/lib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cmeon/SimplexLP/eigen /home/cmeon/SimplexLP/eigen/test /home/cmeon/SimplexLP/lib /home/cmeon/SimplexLP/lib/test /home/cmeon/SimplexLP/lib/test/CMakeFiles/array_for_matrix_1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/array_for_matrix_1.dir/depend

