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
include test/CMakeFiles/sparse_product_2.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/sparse_product_2.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/sparse_product_2.dir/flags.make

test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o: test/CMakeFiles/sparse_product_2.dir/flags.make
test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o: /home/cmeon/SimplexLP/eigen/test/sparse_product.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cmeon/SimplexLP/lib/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o -c /home/cmeon/SimplexLP/eigen/test/sparse_product.cpp

test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sparse_product_2.dir/sparse_product.cpp.i"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cmeon/SimplexLP/eigen/test/sparse_product.cpp > CMakeFiles/sparse_product_2.dir/sparse_product.cpp.i

test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sparse_product_2.dir/sparse_product.cpp.s"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cmeon/SimplexLP/eigen/test/sparse_product.cpp -o CMakeFiles/sparse_product_2.dir/sparse_product.cpp.s

test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o.requires:
.PHONY : test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o.requires

test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o.provides: test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/sparse_product_2.dir/build.make test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o.provides.build
.PHONY : test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o.provides

test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o.provides.build: test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o

# Object files for target sparse_product_2
sparse_product_2_OBJECTS = \
"CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o"

# External object files for target sparse_product_2
sparse_product_2_EXTERNAL_OBJECTS =

test/sparse_product_2: test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o
test/sparse_product_2: test/CMakeFiles/sparse_product_2.dir/build.make
test/sparse_product_2: test/CMakeFiles/sparse_product_2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable sparse_product_2"
	cd /home/cmeon/SimplexLP/lib/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sparse_product_2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/sparse_product_2.dir/build: test/sparse_product_2
.PHONY : test/CMakeFiles/sparse_product_2.dir/build

test/CMakeFiles/sparse_product_2.dir/requires: test/CMakeFiles/sparse_product_2.dir/sparse_product.cpp.o.requires
.PHONY : test/CMakeFiles/sparse_product_2.dir/requires

test/CMakeFiles/sparse_product_2.dir/clean:
	cd /home/cmeon/SimplexLP/lib/test && $(CMAKE_COMMAND) -P CMakeFiles/sparse_product_2.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/sparse_product_2.dir/clean

test/CMakeFiles/sparse_product_2.dir/depend:
	cd /home/cmeon/SimplexLP/lib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cmeon/SimplexLP/eigen /home/cmeon/SimplexLP/eigen/test /home/cmeon/SimplexLP/lib /home/cmeon/SimplexLP/lib/test /home/cmeon/SimplexLP/lib/test/CMakeFiles/sparse_product_2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/sparse_product_2.dir/depend

