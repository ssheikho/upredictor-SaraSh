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
CMAKE_SOURCE_DIR = "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator"

# Include any dependencies generated for this target.
include CMakeFiles/ViconAnnotator.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ViconAnnotator.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ViconAnnotator.dir/flags.make

CMakeFiles/ViconAnnotator.dir/WristDisplacement.o: CMakeFiles/ViconAnnotator.dir/flags.make
CMakeFiles/ViconAnnotator.dir/WristDisplacement.o: WristDisplacement.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ViconAnnotator.dir/WristDisplacement.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ViconAnnotator.dir/WristDisplacement.o -c "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/WristDisplacement.cpp"

CMakeFiles/ViconAnnotator.dir/WristDisplacement.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ViconAnnotator.dir/WristDisplacement.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/WristDisplacement.cpp" > CMakeFiles/ViconAnnotator.dir/WristDisplacement.i

CMakeFiles/ViconAnnotator.dir/WristDisplacement.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ViconAnnotator.dir/WristDisplacement.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/WristDisplacement.cpp" -o CMakeFiles/ViconAnnotator.dir/WristDisplacement.s

CMakeFiles/ViconAnnotator.dir/WristDisplacement.o.requires:

.PHONY : CMakeFiles/ViconAnnotator.dir/WristDisplacement.o.requires

CMakeFiles/ViconAnnotator.dir/WristDisplacement.o.provides: CMakeFiles/ViconAnnotator.dir/WristDisplacement.o.requires
	$(MAKE) -f CMakeFiles/ViconAnnotator.dir/build.make CMakeFiles/ViconAnnotator.dir/WristDisplacement.o.provides.build
.PHONY : CMakeFiles/ViconAnnotator.dir/WristDisplacement.o.provides

CMakeFiles/ViconAnnotator.dir/WristDisplacement.o.provides.build: CMakeFiles/ViconAnnotator.dir/WristDisplacement.o


CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o: CMakeFiles/ViconAnnotator.dir/flags.make
CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o: NameCrossIndex.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o -c "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/NameCrossIndex.cpp"

CMakeFiles/ViconAnnotator.dir/NameCrossIndex.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ViconAnnotator.dir/NameCrossIndex.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/NameCrossIndex.cpp" > CMakeFiles/ViconAnnotator.dir/NameCrossIndex.i

CMakeFiles/ViconAnnotator.dir/NameCrossIndex.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ViconAnnotator.dir/NameCrossIndex.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/NameCrossIndex.cpp" -o CMakeFiles/ViconAnnotator.dir/NameCrossIndex.s

CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o.requires:

.PHONY : CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o.requires

CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o.provides: CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o.requires
	$(MAKE) -f CMakeFiles/ViconAnnotator.dir/build.make CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o.provides.build
.PHONY : CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o.provides

CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o.provides.build: CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o


CMakeFiles/ViconAnnotator.dir/ViconData.o: CMakeFiles/ViconAnnotator.dir/flags.make
CMakeFiles/ViconAnnotator.dir/ViconData.o: ViconData.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/ViconAnnotator.dir/ViconData.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ViconAnnotator.dir/ViconData.o -c "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconData.cpp"

CMakeFiles/ViconAnnotator.dir/ViconData.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ViconAnnotator.dir/ViconData.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconData.cpp" > CMakeFiles/ViconAnnotator.dir/ViconData.i

CMakeFiles/ViconAnnotator.dir/ViconData.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ViconAnnotator.dir/ViconData.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconData.cpp" -o CMakeFiles/ViconAnnotator.dir/ViconData.s

CMakeFiles/ViconAnnotator.dir/ViconData.o.requires:

.PHONY : CMakeFiles/ViconAnnotator.dir/ViconData.o.requires

CMakeFiles/ViconAnnotator.dir/ViconData.o.provides: CMakeFiles/ViconAnnotator.dir/ViconData.o.requires
	$(MAKE) -f CMakeFiles/ViconAnnotator.dir/build.make CMakeFiles/ViconAnnotator.dir/ViconData.o.provides.build
.PHONY : CMakeFiles/ViconAnnotator.dir/ViconData.o.provides

CMakeFiles/ViconAnnotator.dir/ViconData.o.provides.build: CMakeFiles/ViconAnnotator.dir/ViconData.o


CMakeFiles/ViconAnnotator.dir/ViconFunctions.o: CMakeFiles/ViconAnnotator.dir/flags.make
CMakeFiles/ViconAnnotator.dir/ViconFunctions.o: ViconFunctions.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/ViconAnnotator.dir/ViconFunctions.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ViconAnnotator.dir/ViconFunctions.o -c "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconFunctions.cpp"

CMakeFiles/ViconAnnotator.dir/ViconFunctions.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ViconAnnotator.dir/ViconFunctions.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconFunctions.cpp" > CMakeFiles/ViconAnnotator.dir/ViconFunctions.i

CMakeFiles/ViconAnnotator.dir/ViconFunctions.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ViconAnnotator.dir/ViconFunctions.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconFunctions.cpp" -o CMakeFiles/ViconAnnotator.dir/ViconFunctions.s

CMakeFiles/ViconAnnotator.dir/ViconFunctions.o.requires:

.PHONY : CMakeFiles/ViconAnnotator.dir/ViconFunctions.o.requires

CMakeFiles/ViconAnnotator.dir/ViconFunctions.o.provides: CMakeFiles/ViconAnnotator.dir/ViconFunctions.o.requires
	$(MAKE) -f CMakeFiles/ViconAnnotator.dir/build.make CMakeFiles/ViconAnnotator.dir/ViconFunctions.o.provides.build
.PHONY : CMakeFiles/ViconAnnotator.dir/ViconFunctions.o.provides

CMakeFiles/ViconAnnotator.dir/ViconFunctions.o.provides.build: CMakeFiles/ViconAnnotator.dir/ViconFunctions.o


CMakeFiles/ViconAnnotator.dir/ViconJoints.o: CMakeFiles/ViconAnnotator.dir/flags.make
CMakeFiles/ViconAnnotator.dir/ViconJoints.o: ViconJoints.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/ViconAnnotator.dir/ViconJoints.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ViconAnnotator.dir/ViconJoints.o -c "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconJoints.cpp"

CMakeFiles/ViconAnnotator.dir/ViconJoints.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ViconAnnotator.dir/ViconJoints.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconJoints.cpp" > CMakeFiles/ViconAnnotator.dir/ViconJoints.i

CMakeFiles/ViconAnnotator.dir/ViconJoints.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ViconAnnotator.dir/ViconJoints.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconJoints.cpp" -o CMakeFiles/ViconAnnotator.dir/ViconJoints.s

CMakeFiles/ViconAnnotator.dir/ViconJoints.o.requires:

.PHONY : CMakeFiles/ViconAnnotator.dir/ViconJoints.o.requires

CMakeFiles/ViconAnnotator.dir/ViconJoints.o.provides: CMakeFiles/ViconAnnotator.dir/ViconJoints.o.requires
	$(MAKE) -f CMakeFiles/ViconAnnotator.dir/build.make CMakeFiles/ViconAnnotator.dir/ViconJoints.o.provides.build
.PHONY : CMakeFiles/ViconAnnotator.dir/ViconJoints.o.provides

CMakeFiles/ViconAnnotator.dir/ViconJoints.o.provides.build: CMakeFiles/ViconAnnotator.dir/ViconJoints.o


CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o: CMakeFiles/ViconAnnotator.dir/flags.make
CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o: ViconTrajectories.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o -c "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconTrajectories.cpp"

CMakeFiles/ViconAnnotator.dir/ViconTrajectories.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ViconAnnotator.dir/ViconTrajectories.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconTrajectories.cpp" > CMakeFiles/ViconAnnotator.dir/ViconTrajectories.i

CMakeFiles/ViconAnnotator.dir/ViconTrajectories.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ViconAnnotator.dir/ViconTrajectories.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconTrajectories.cpp" -o CMakeFiles/ViconAnnotator.dir/ViconTrajectories.s

CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o.requires:

.PHONY : CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o.requires

CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o.provides: CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o.requires
	$(MAKE) -f CMakeFiles/ViconAnnotator.dir/build.make CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o.provides.build
.PHONY : CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o.provides

CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o.provides.build: CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o


CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o: CMakeFiles/ViconAnnotator.dir/flags.make
CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o: ViconAnnotator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o -c "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconAnnotator.cpp"

CMakeFiles/ViconAnnotator.dir/ViconAnnotator.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ViconAnnotator.dir/ViconAnnotator.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconAnnotator.cpp" > CMakeFiles/ViconAnnotator.dir/ViconAnnotator.i

CMakeFiles/ViconAnnotator.dir/ViconAnnotator.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ViconAnnotator.dir/ViconAnnotator.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/ViconAnnotator.cpp" -o CMakeFiles/ViconAnnotator.dir/ViconAnnotator.s

CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o.requires:

.PHONY : CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o.requires

CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o.provides: CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o.requires
	$(MAKE) -f CMakeFiles/ViconAnnotator.dir/build.make CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o.provides.build
.PHONY : CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o.provides

CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o.provides.build: CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o


# Object files for target ViconAnnotator
ViconAnnotator_OBJECTS = \
"CMakeFiles/ViconAnnotator.dir/WristDisplacement.o" \
"CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o" \
"CMakeFiles/ViconAnnotator.dir/ViconData.o" \
"CMakeFiles/ViconAnnotator.dir/ViconFunctions.o" \
"CMakeFiles/ViconAnnotator.dir/ViconJoints.o" \
"CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o" \
"CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o"

# External object files for target ViconAnnotator
ViconAnnotator_EXTERNAL_OBJECTS =

ViconAnnotator: CMakeFiles/ViconAnnotator.dir/WristDisplacement.o
ViconAnnotator: CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o
ViconAnnotator: CMakeFiles/ViconAnnotator.dir/ViconData.o
ViconAnnotator: CMakeFiles/ViconAnnotator.dir/ViconFunctions.o
ViconAnnotator: CMakeFiles/ViconAnnotator.dir/ViconJoints.o
ViconAnnotator: CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o
ViconAnnotator: CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o
ViconAnnotator: CMakeFiles/ViconAnnotator.dir/build.make
ViconAnnotator: CMakeFiles/ViconAnnotator.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable ViconAnnotator"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ViconAnnotator.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ViconAnnotator.dir/build: ViconAnnotator

.PHONY : CMakeFiles/ViconAnnotator.dir/build

CMakeFiles/ViconAnnotator.dir/requires: CMakeFiles/ViconAnnotator.dir/WristDisplacement.o.requires
CMakeFiles/ViconAnnotator.dir/requires: CMakeFiles/ViconAnnotator.dir/NameCrossIndex.o.requires
CMakeFiles/ViconAnnotator.dir/requires: CMakeFiles/ViconAnnotator.dir/ViconData.o.requires
CMakeFiles/ViconAnnotator.dir/requires: CMakeFiles/ViconAnnotator.dir/ViconFunctions.o.requires
CMakeFiles/ViconAnnotator.dir/requires: CMakeFiles/ViconAnnotator.dir/ViconJoints.o.requires
CMakeFiles/ViconAnnotator.dir/requires: CMakeFiles/ViconAnnotator.dir/ViconTrajectories.o.requires
CMakeFiles/ViconAnnotator.dir/requires: CMakeFiles/ViconAnnotator.dir/ViconAnnotator.o.requires

.PHONY : CMakeFiles/ViconAnnotator.dir/requires

CMakeFiles/ViconAnnotator.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ViconAnnotator.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ViconAnnotator.dir/clean

CMakeFiles/ViconAnnotator.dir/depend:
	cd "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator" "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator" "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator" "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator" "/home/sara/Dropbox/Reaching Study/GitRepo/ReachPredictor/upredictor/uPredictor/ViconAnnotator/CMakeFiles/ViconAnnotator.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/ViconAnnotator.dir/depend
