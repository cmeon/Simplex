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
include test/CMakeFiles/product_selfadjoint_6.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/product_selfadjoint_6.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/product_selfadjoint_6.dir/flags.make

test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o: test/CMakeFiles/product_selfadjoint_6.dir/flags.make
test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o: /home/cmeon/SimplexLP/eigen/test/product_selfadjoint.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cmeon/SimplexLP/lib/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o -c /home/cmeon/SimplexLP/eigen/test/product_selfadjoint.cpp

test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.i"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cmeon/SimplexLP/eigen/test/product_selfadjoint.cpp > CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.i

test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.s"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cmeon/SimplexLP/eigen/test/product_selfadjoint.cpp -o CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.s

test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o.requires:
.PHONY : test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o.requires

test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o.provides: test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/product_selfadjoint_6.dir/build.make test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o.provides.build
.PHONY : test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o.provides

test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o.provides.build: test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o

# Object files for target product_selfadjoint_6
product_selfadjoint_6_OBJECTS = \
"CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o"

# External object files for target product_selfadjoint_6
product_selfadjoint_6_EXTERNAL_OBJECTS =

test/product_selfadjoint_6: test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o
test/product_selfadjoint_6: test/CMakeFiles/product_selfadjoint_6.dir/build.make
test/product_selfadjoint_6: test/CMakeFiles/product_selfadjoint_6.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable product_selfadjoint_6"
	cd /home/cmeon/SimplexLP/lib/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/product_selfadjoint_6.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/product_selfadjoint_6.dir/build: test/product_selfadjoint_6
.PHONY : test/CMakeFiles/product_selfadjoint_6.dir/build

test/CMakeFiles/product_selfadjoint_6.dir/requires: test/CMakeFiles/product_selfadjoint_6.dir/product_selfadjoint.cpp.o.requires
.PHONY : test/CMakeFiles/product_selfadjoint_6.dir/requires

test/CMakeFiles/product_selfadjoint_6.dir/clean:
	cd /home/cmeon/SimplexLP/lib/test && $(CMAKE_COMMAND) -P CMakeFiles/product_selfadjoint_6.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/product_selfadjoint_6.dir/clean

test/CMakeFiles/product_selfadjoint_6.dir/depend:
	cd /home/cmeon/SimplexLP/lib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cmeon/SimplexLP/eigen /home/cmeon/SimplexLP/eigen/test /home/cmeon/SimplexLP/lib /home/cmeon/SimplexLP/lib/test /home/cmeon/SimplexLP/lib/test/CMakeFiles/product_selfadjoint_6.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/product_selfadjoint_6.dir/depend

