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
include test/CMakeFiles/mapstaticmethods_5.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/mapstaticmethods_5.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/mapstaticmethods_5.dir/flags.make

test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o: test/CMakeFiles/mapstaticmethods_5.dir/flags.make
test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o: /home/cmeon/SimplexLP/eigen/test/mapstaticmethods.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cmeon/SimplexLP/lib/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o -c /home/cmeon/SimplexLP/eigen/test/mapstaticmethods.cpp

test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.i"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cmeon/SimplexLP/eigen/test/mapstaticmethods.cpp > CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.i

test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.s"
	cd /home/cmeon/SimplexLP/lib/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cmeon/SimplexLP/eigen/test/mapstaticmethods.cpp -o CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.s

test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o.requires:
.PHONY : test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o.requires

test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o.provides: test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/mapstaticmethods_5.dir/build.make test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o.provides.build
.PHONY : test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o.provides

test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o.provides.build: test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o

# Object files for target mapstaticmethods_5
mapstaticmethods_5_OBJECTS = \
"CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o"

# External object files for target mapstaticmethods_5
mapstaticmethods_5_EXTERNAL_OBJECTS =

test/mapstaticmethods_5: test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o
test/mapstaticmethods_5: test/CMakeFiles/mapstaticmethods_5.dir/build.make
test/mapstaticmethods_5: test/CMakeFiles/mapstaticmethods_5.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable mapstaticmethods_5"
	cd /home/cmeon/SimplexLP/lib/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mapstaticmethods_5.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/mapstaticmethods_5.dir/build: test/mapstaticmethods_5
.PHONY : test/CMakeFiles/mapstaticmethods_5.dir/build

test/CMakeFiles/mapstaticmethods_5.dir/requires: test/CMakeFiles/mapstaticmethods_5.dir/mapstaticmethods.cpp.o.requires
.PHONY : test/CMakeFiles/mapstaticmethods_5.dir/requires

test/CMakeFiles/mapstaticmethods_5.dir/clean:
	cd /home/cmeon/SimplexLP/lib/test && $(CMAKE_COMMAND) -P CMakeFiles/mapstaticmethods_5.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/mapstaticmethods_5.dir/clean

test/CMakeFiles/mapstaticmethods_5.dir/depend:
	cd /home/cmeon/SimplexLP/lib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cmeon/SimplexLP/eigen /home/cmeon/SimplexLP/eigen/test /home/cmeon/SimplexLP/lib /home/cmeon/SimplexLP/lib/test /home/cmeon/SimplexLP/lib/test/CMakeFiles/mapstaticmethods_5.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/mapstaticmethods_5.dir/depend

