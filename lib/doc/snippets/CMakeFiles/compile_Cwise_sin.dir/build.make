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
include doc/snippets/CMakeFiles/compile_Cwise_sin.dir/depend.make

# Include the progress variables for this target.
include doc/snippets/CMakeFiles/compile_Cwise_sin.dir/progress.make

# Include the compile flags for this target's objects.
include doc/snippets/CMakeFiles/compile_Cwise_sin.dir/flags.make

doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o: doc/snippets/CMakeFiles/compile_Cwise_sin.dir/flags.make
doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o: doc/snippets/compile_Cwise_sin.cpp
doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o: /home/cmeon/SimplexLP/eigen/doc/snippets/Cwise_sin.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cmeon/SimplexLP/lib/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o"
	cd /home/cmeon/SimplexLP/lib/doc/snippets && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o -c /home/cmeon/SimplexLP/lib/doc/snippets/compile_Cwise_sin.cpp

doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.i"
	cd /home/cmeon/SimplexLP/lib/doc/snippets && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cmeon/SimplexLP/lib/doc/snippets/compile_Cwise_sin.cpp > CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.i

doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.s"
	cd /home/cmeon/SimplexLP/lib/doc/snippets && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cmeon/SimplexLP/lib/doc/snippets/compile_Cwise_sin.cpp -o CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.s

doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o.requires:
.PHONY : doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o.requires

doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o.provides: doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o.requires
	$(MAKE) -f doc/snippets/CMakeFiles/compile_Cwise_sin.dir/build.make doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o.provides.build
.PHONY : doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o.provides

doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o.provides.build: doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o

# Object files for target compile_Cwise_sin
compile_Cwise_sin_OBJECTS = \
"CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o"

# External object files for target compile_Cwise_sin
compile_Cwise_sin_EXTERNAL_OBJECTS =

doc/snippets/compile_Cwise_sin: doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o
doc/snippets/compile_Cwise_sin: doc/snippets/CMakeFiles/compile_Cwise_sin.dir/build.make
doc/snippets/compile_Cwise_sin: doc/snippets/CMakeFiles/compile_Cwise_sin.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable compile_Cwise_sin"
	cd /home/cmeon/SimplexLP/lib/doc/snippets && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/compile_Cwise_sin.dir/link.txt --verbose=$(VERBOSE)
	cd /home/cmeon/SimplexLP/lib/doc/snippets && ./compile_Cwise_sin >/home/cmeon/SimplexLP/lib/doc/snippets/Cwise_sin.out

# Rule to build all files generated by this target.
doc/snippets/CMakeFiles/compile_Cwise_sin.dir/build: doc/snippets/compile_Cwise_sin
.PHONY : doc/snippets/CMakeFiles/compile_Cwise_sin.dir/build

doc/snippets/CMakeFiles/compile_Cwise_sin.dir/requires: doc/snippets/CMakeFiles/compile_Cwise_sin.dir/compile_Cwise_sin.cpp.o.requires
.PHONY : doc/snippets/CMakeFiles/compile_Cwise_sin.dir/requires

doc/snippets/CMakeFiles/compile_Cwise_sin.dir/clean:
	cd /home/cmeon/SimplexLP/lib/doc/snippets && $(CMAKE_COMMAND) -P CMakeFiles/compile_Cwise_sin.dir/cmake_clean.cmake
.PHONY : doc/snippets/CMakeFiles/compile_Cwise_sin.dir/clean

doc/snippets/CMakeFiles/compile_Cwise_sin.dir/depend:
	cd /home/cmeon/SimplexLP/lib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cmeon/SimplexLP/eigen /home/cmeon/SimplexLP/eigen/doc/snippets /home/cmeon/SimplexLP/lib /home/cmeon/SimplexLP/lib/doc/snippets /home/cmeon/SimplexLP/lib/doc/snippets/CMakeFiles/compile_Cwise_sin.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/snippets/CMakeFiles/compile_Cwise_sin.dir/depend

