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
include unsupported/test/CMakeFiles/matrix_function_6.dir/depend.make

# Include the progress variables for this target.
include unsupported/test/CMakeFiles/matrix_function_6.dir/progress.make

# Include the compile flags for this target's objects.
include unsupported/test/CMakeFiles/matrix_function_6.dir/flags.make

unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o: unsupported/test/CMakeFiles/matrix_function_6.dir/flags.make
unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o: /home/cmeon/SimplexLP/eigen/unsupported/test/matrix_function.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cmeon/SimplexLP/lib/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o"
	cd /home/cmeon/SimplexLP/lib/unsupported/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o -c /home/cmeon/SimplexLP/eigen/unsupported/test/matrix_function.cpp

unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix_function_6.dir/matrix_function.cpp.i"
	cd /home/cmeon/SimplexLP/lib/unsupported/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cmeon/SimplexLP/eigen/unsupported/test/matrix_function.cpp > CMakeFiles/matrix_function_6.dir/matrix_function.cpp.i

unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix_function_6.dir/matrix_function.cpp.s"
	cd /home/cmeon/SimplexLP/lib/unsupported/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cmeon/SimplexLP/eigen/unsupported/test/matrix_function.cpp -o CMakeFiles/matrix_function_6.dir/matrix_function.cpp.s

unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o.requires:
.PHONY : unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o.requires

unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o.provides: unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o.requires
	$(MAKE) -f unsupported/test/CMakeFiles/matrix_function_6.dir/build.make unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o.provides.build
.PHONY : unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o.provides

unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o.provides.build: unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o

# Object files for target matrix_function_6
matrix_function_6_OBJECTS = \
"CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o"

# External object files for target matrix_function_6
matrix_function_6_EXTERNAL_OBJECTS =

unsupported/test/matrix_function_6: unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o
unsupported/test/matrix_function_6: unsupported/test/CMakeFiles/matrix_function_6.dir/build.make
unsupported/test/matrix_function_6: unsupported/test/CMakeFiles/matrix_function_6.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable matrix_function_6"
	cd /home/cmeon/SimplexLP/lib/unsupported/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/matrix_function_6.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
unsupported/test/CMakeFiles/matrix_function_6.dir/build: unsupported/test/matrix_function_6
.PHONY : unsupported/test/CMakeFiles/matrix_function_6.dir/build

unsupported/test/CMakeFiles/matrix_function_6.dir/requires: unsupported/test/CMakeFiles/matrix_function_6.dir/matrix_function.cpp.o.requires
.PHONY : unsupported/test/CMakeFiles/matrix_function_6.dir/requires

unsupported/test/CMakeFiles/matrix_function_6.dir/clean:
	cd /home/cmeon/SimplexLP/lib/unsupported/test && $(CMAKE_COMMAND) -P CMakeFiles/matrix_function_6.dir/cmake_clean.cmake
.PHONY : unsupported/test/CMakeFiles/matrix_function_6.dir/clean

unsupported/test/CMakeFiles/matrix_function_6.dir/depend:
	cd /home/cmeon/SimplexLP/lib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cmeon/SimplexLP/eigen /home/cmeon/SimplexLP/eigen/unsupported/test /home/cmeon/SimplexLP/lib /home/cmeon/SimplexLP/lib/unsupported/test /home/cmeon/SimplexLP/lib/unsupported/test/CMakeFiles/matrix_function_6.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : unsupported/test/CMakeFiles/matrix_function_6.dir/depend

