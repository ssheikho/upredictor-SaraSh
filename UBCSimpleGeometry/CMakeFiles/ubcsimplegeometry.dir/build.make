# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor"

# Include any dependencies generated for this target.
include UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/depend.make

# Include the progress variables for this target.
include UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/progress.make

# Include the compile flags for this target's objects.
include UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/flags.make

UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o: UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/flags.make
UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o: UBCSimpleGeometry/UBCSimpleGeometry.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o"
	cd "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/UBCSimpleGeometry" && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o   -c "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/UBCSimpleGeometry/UBCSimpleGeometry.c"

UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.i"
	cd "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/UBCSimpleGeometry" && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/UBCSimpleGeometry/UBCSimpleGeometry.c" > CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.i

UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.s"
	cd "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/UBCSimpleGeometry" && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/UBCSimpleGeometry/UBCSimpleGeometry.c" -o CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.s

UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o.requires:

.PHONY : UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o.requires

UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o.provides: UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o.requires
	$(MAKE) -f UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/build.make UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o.provides.build
.PHONY : UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o.provides

UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o.provides.build: UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o


# Object files for target ubcsimplegeometry
ubcsimplegeometry_OBJECTS = \
"CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o"

# External object files for target ubcsimplegeometry
ubcsimplegeometry_EXTERNAL_OBJECTS =

UBCSimpleGeometry/libubcsimplegeometry.a: UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o
UBCSimpleGeometry/libubcsimplegeometry.a: UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/build.make
UBCSimpleGeometry/libubcsimplegeometry.a: UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library libubcsimplegeometry.a"
	cd "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/UBCSimpleGeometry" && $(CMAKE_COMMAND) -P CMakeFiles/ubcsimplegeometry.dir/cmake_clean_target.cmake
	cd "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/UBCSimpleGeometry" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ubcsimplegeometry.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/build: UBCSimpleGeometry/libubcsimplegeometry.a

.PHONY : UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/build

UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/requires: UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/UBCSimpleGeometry.o.requires

.PHONY : UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/requires

UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/clean:
	cd "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/UBCSimpleGeometry" && $(CMAKE_COMMAND) -P CMakeFiles/ubcsimplegeometry.dir/cmake_clean.cmake
.PHONY : UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/clean

UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/depend:
	cd "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor" "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/UBCSimpleGeometry" "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor" "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/UBCSimpleGeometry" "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : UBCSimpleGeometry/CMakeFiles/ubcsimplegeometry.dir/depend

