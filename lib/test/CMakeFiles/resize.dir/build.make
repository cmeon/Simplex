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
include test/CMakeFiles/resize.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/resize.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/resize.dir/flags.make

test/CMakeFiles/resize.dir/resize.cpp.o: test/CMakeFiles/resize.dir/flags.make
test/CMakeFiles/resize.dir/resize.cpp.o: /home/cmeon/SimplexLP/eigen/test/resize.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cmeon/SimplexLP/lib/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object test/CMakeFiles/resize.dir/resize.cpp.o"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/resize.dir/resize.cpp.o -c /home/cmeon/SimplexLP/eigen/test/resize.cpp

test/CMakeFiles/resize.dir/resize.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/resize.dir/resize.cpp.i"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cmeon/SimplexLP/eigen/test/resize.cpp > CMakeFiles/resize.dir/resize.cpp.i

test/CMakeFiles/resize.dir/resize.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/resize.dir/resize.cpp.s"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cmeon/SimplexLP/eigen/test/resize.cpp -o CMakeFiles/resize.dir/resize.cpp.s

test/CMakeFiles/resize.dir/resize.cpp.o.requires:
.PHONY : test/CMakeFiles/resize.dir/resize.cpp.o.requires

test/CMakeFiles/resize.dir/resize.cpp.o.provides: test/CMakeFiles/resize.dir/resize.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/resize.dir/build.make test/CMakeFiles/resize.dir/resize.cpp.o.provides.build
.PHONY : test/CMakeFiles/resize.dir/resize.cpp.o.provides

test/CMakeFiles/resize.dir/resize.cpp.o.provides.build: test/CMakeFiles/resize.dir/resize.cpp.o

# Object files for target resize
resize_OBJECTS = \
"CMakeFiles/resize.dir/resize.cpp.o"

# External object files for target resize
resize_EXTERNAL_OBJECTS =

test/resize: test/CMakeFiles/resize.dir/resize.cpp.o
test/resize: test/CMakeFiles/resize.dir/build.make
test/resize: test/CMakeFiles/resize.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable resize"
	cd /home/cmeon/SimplexLP/lib/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/resize.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/resize.dir/build: test/resize
.PHONY : test/CMakeFiles/resize.dir/build

test/CMakeFiles/resize.dir/requires: test/CMakeFiles/resize.dir/resize.cpp.o.requires
.PHONY : test/CMakeFiles/resize.dir/requires

test/CMakeFiles/resize.dir/clean:
	cd /home/cmeon/SimplexLP/lib/test && $(CMAKE_COMMAND) -P CMakeFiles/resize.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/resize.dir/clean

test/CMakeFiles/resize.dir/depend:
	cd /home/cmeon/SimplexLP/lib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cmeon/SimplexLP/eigen /home/cmeon/SimplexLP/eigen/test /home/cmeon/SimplexLP/lib /home/cmeon/SimplexLP/lib/test /home/cmeon/SimplexLP/lib/test/CMakeFiles/resize.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/resize.dir/depend

