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
include doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/depend.make

# Include the progress variables for this target.
include doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/progress.make

# Include the compile flags for this target's objects.
include doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/flags.make

doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o: doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/flags.make
doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o: doc/snippets/compile_MatrixBase_ones.cpp
doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o: /home/cmeon/SimplexLP/eigen/doc/snippets/MatrixBase_ones.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cmeon/SimplexLP/lib/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o"
	cd /home/cmeon/SimplexLP/lib/doc/snippets && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o -c /home/cmeon/SimplexLP/lib/doc/snippets/compile_MatrixBase_ones.cpp

doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.i"
	cd /home/cmeon/SimplexLP/lib/doc/snippets && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cmeon/SimplexLP/lib/doc/snippets/compile_MatrixBase_ones.cpp > CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.i

doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.s"
	cd /home/cmeon/SimplexLP/lib/doc/snippets && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cmeon/SimplexLP/lib/doc/snippets/compile_MatrixBase_ones.cpp -o CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.s

doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o.requires:
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o.requires

doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o.provides: doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o.requires
	$(MAKE) -f doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/build.make doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o.provides.build
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o.provides

doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o.provides.build: doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o

# Object files for target compile_MatrixBase_ones
compile_MatrixBase_ones_OBJECTS = \
"CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o"

# External object files for target compile_MatrixBase_ones
compile_MatrixBase_ones_EXTERNAL_OBJECTS =

doc/snippets/compile_MatrixBase_ones: doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o
doc/snippets/compile_MatrixBase_ones: doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/build.make
doc/snippets/compile_MatrixBase_ones: doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable compile_MatrixBase_ones"
	cd /home/cmeon/SimplexLP/lib/doc/snippets && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/compile_MatrixBase_ones.dir/link.txt --verbose=$(VERBOSE)
	cd /home/cmeon/SimplexLP/lib/doc/snippets && ./compile_MatrixBase_ones >/home/cmeon/SimplexLP/lib/doc/snippets/MatrixBase_ones.out

# Rule to build all files generated by this target.
doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/build: doc/snippets/compile_MatrixBase_ones
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/build

doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/requires: doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/compile_MatrixBase_ones.cpp.o.requires
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/requires

doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/clean:
	cd /home/cmeon/SimplexLP/lib/doc/snippets && $(CMAKE_COMMAND) -P CMakeFiles/compile_MatrixBase_ones.dir/cmake_clean.cmake
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/clean

doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/depend:
	cd /home/cmeon/SimplexLP/lib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cmeon/SimplexLP/eigen /home/cmeon/SimplexLP/eigen/doc/snippets /home/cmeon/SimplexLP/lib /home/cmeon/SimplexLP/lib/doc/snippets /home/cmeon/SimplexLP/lib/doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/snippets/CMakeFiles/compile_MatrixBase_ones.dir/depend

