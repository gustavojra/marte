# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /home/aroeira/psi4conda/bin/cmake

# The command to remove a file.
RM = /home/aroeira/psi4conda/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/aroeira/marte

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/aroeira/marte

# Include any dependencies generated for this target.
include CMakeFiles/marte.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/marte.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/marte.dir/flags.make

CMakeFiles/marte.dir/ccsd.cc.o: CMakeFiles/marte.dir/flags.make
CMakeFiles/marte.dir/ccsd.cc.o: ccsd.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aroeira/marte/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/marte.dir/ccsd.cc.o"
	/home/aroeira/psi4conda/bin/x86_64-conda_cos6-linux-gnu-g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/marte.dir/ccsd.cc.o -c /home/aroeira/marte/ccsd.cc

CMakeFiles/marte.dir/ccsd.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/marte.dir/ccsd.cc.i"
	/home/aroeira/psi4conda/bin/x86_64-conda_cos6-linux-gnu-g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aroeira/marte/ccsd.cc > CMakeFiles/marte.dir/ccsd.cc.i

CMakeFiles/marte.dir/ccsd.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/marte.dir/ccsd.cc.s"
	/home/aroeira/psi4conda/bin/x86_64-conda_cos6-linux-gnu-g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aroeira/marte/ccsd.cc -o CMakeFiles/marte.dir/ccsd.cc.s

# Object files for target marte
marte_OBJECTS = \
"CMakeFiles/marte.dir/ccsd.cc.o"

# External object files for target marte
marte_EXTERNAL_OBJECTS =

marte.so: CMakeFiles/marte.dir/ccsd.cc.o
marte.so: CMakeFiles/marte.dir/build.make
marte.so: /home/aroeira/psi4conda/lib/python3.6/site-packages/psi4/core.cpython-36m-x86_64-linux-gnu.so
marte.so: /home/aroeira/psi4conda/lib/libiomp5.so
marte.so: /home/aroeira/psi4conda/x86_64-conda_cos6-linux-gnu/sysroot/lib/libgomp.so
marte.so: /home/aroeira/psi4conda/x86_64-conda_cos6-linux-gnu/sysroot/usr/lib/libpthread.so
marte.so: CMakeFiles/marte.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/aroeira/marte/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module marte.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/marte.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/marte.dir/build: marte.so

.PHONY : CMakeFiles/marte.dir/build

CMakeFiles/marte.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/marte.dir/cmake_clean.cmake
.PHONY : CMakeFiles/marte.dir/clean

CMakeFiles/marte.dir/depend:
	cd /home/aroeira/marte && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/aroeira/marte /home/aroeira/marte /home/aroeira/marte /home/aroeira/marte /home/aroeira/marte/CMakeFiles/marte.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/marte.dir/depend
