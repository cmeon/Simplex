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
include bench/spbench/CMakeFiles/test_sparseLU.dir/depend.make

# Include the progress variables for this target.
include bench/spbench/CMakeFiles/test_sparseLU.dir/progress.make

# Include the compile flags for this target's objects.
include bench/spbench/CMakeFiles/test_sparseLU.dir/flags.make

bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o: bench/spbench/CMakeFiles/test_sparseLU.dir/flags.make
bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o: /home/cmeon/SimplexLP/eigen/bench/spbench/test_sparseLU.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cmeon/SimplexLP/lib/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o"
	cd /home/cmeon/SimplexLP/lib/bench/spbench && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o -c /home/cmeon/SimplexLP/eigen/bench/spbench/test_sparseLU.cpp

bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.i"
	cd /home/cmeon/SimplexLP/lib/bench/spbench && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cmeon/SimplexLP/eigen/bench/spbench/test_sparseLU.cpp > CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.i

bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.s"
	cd /home/cmeon/SimplexLP/lib/bench/spbench && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cmeon/SimplexLP/eigen/bench/spbench/test_sparseLU.cpp -o CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.s

bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o.requires:
.PHONY : bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o.requires

bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o.provides: bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o.requires
	$(MAKE) -f bench/spbench/CMakeFiles/test_sparseLU.dir/build.make bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o.provides.build
.PHONY : bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o.provides

bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o.provides.build: bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o

# Object files for target test_sparseLU
test_sparseLU_OBJECTS = \
"CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o"

# External object files for target test_sparseLU
test_sparseLU_EXTERNAL_OBJECTS =

bench/spbench/test_sparseLU: bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o
bench/spbench/test_sparseLU: bench/spbench/CMakeFiles/test_sparseLU.dir/build.make
bench/spbench/test_sparseLU: /usr/lib/x86_64-linux-gnu/librt.so
bench/spbench/test_sparseLU: bench/spbench/CMakeFiles/test_sparseLU.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable test_sparseLU"
	cd /home/cmeon/SimplexLP/lib/bench/spbench && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_sparseLU.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
bench/spbench/CMakeFiles/test_sparseLU.dir/build: bench/spbench/test_sparseLU
.PHONY : bench/spbench/CMakeFiles/test_sparseLU.dir/build

bench/spbench/CMakeFiles/test_sparseLU.dir/requires: bench/spbench/CMakeFiles/test_sparseLU.dir/test_sparseLU.cpp.o.requires
.PHONY : bench/spbench/CMakeFiles/test_sparseLU.dir/requires

bench/spbench/CMakeFiles/test_sparseLU.dir/clean:
	cd /home/cmeon/SimplexLP/lib/bench/spbench && $(CMAKE_COMMAND) -P CMakeFiles/test_sparseLU.dir/cmake_clean.cmake
.PHONY : bench/spbench/CMakeFiles/test_sparseLU.dir/clean

bench/spbench/CMakeFiles/test_sparseLU.dir/depend:
	cd /home/cmeon/SimplexLP/lib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cmeon/SimplexLP/eigen /home/cmeon/SimplexLP/eigen/bench/spbench /home/cmeon/SimplexLP/lib /home/cmeon/SimplexLP/lib/bench/spbench /home/cmeon/SimplexLP/lib/bench/spbench/CMakeFiles/test_sparseLU.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : bench/spbench/CMakeFiles/test_sparseLU.dir/depend

