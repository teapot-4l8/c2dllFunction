# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.16

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "D:\JetBrains\CLion 2020.1\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "D:\JetBrains\CLion 2020.1\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\unityProject\PhysicsPrj\Assets\compileFromCpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\unityProject\PhysicsPrj\Assets\compileFromCpp\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/compileFromCpp.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/compileFromCpp.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/compileFromCpp.dir/flags.make

CMakeFiles/compileFromCpp.dir/double_slit2_FFT_tri.c.obj: CMakeFiles/compileFromCpp.dir/flags.make
CMakeFiles/compileFromCpp.dir/double_slit2_FFT_tri.c.obj: ../double_slit2_FFT_tri.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\unityProject\PhysicsPrj\Assets\compileFromCpp\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/compileFromCpp.dir/double_slit2_FFT_tri.c.obj"
	D:\Programs\mingw64\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\compileFromCpp.dir\double_slit2_FFT_tri.c.obj   -c D:\unityProject\PhysicsPrj\Assets\compileFromCpp\double_slit2_FFT_tri.c

CMakeFiles/compileFromCpp.dir/double_slit2_FFT_tri.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/compileFromCpp.dir/double_slit2_FFT_tri.c.i"
	D:\Programs\mingw64\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\unityProject\PhysicsPrj\Assets\compileFromCpp\double_slit2_FFT_tri.c > CMakeFiles\compileFromCpp.dir\double_slit2_FFT_tri.c.i

CMakeFiles/compileFromCpp.dir/double_slit2_FFT_tri.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/compileFromCpp.dir/double_slit2_FFT_tri.c.s"
	D:\Programs\mingw64\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\unityProject\PhysicsPrj\Assets\compileFromCpp\double_slit2_FFT_tri.c -o CMakeFiles\compileFromCpp.dir\double_slit2_FFT_tri.c.s

# Object files for target compileFromCpp
compileFromCpp_OBJECTS = \
"CMakeFiles/compileFromCpp.dir/double_slit2_FFT_tri.c.obj"

# External object files for target compileFromCpp
compileFromCpp_EXTERNAL_OBJECTS =

compileFromCpp.exe: CMakeFiles/compileFromCpp.dir/double_slit2_FFT_tri.c.obj
compileFromCpp.exe: CMakeFiles/compileFromCpp.dir/build.make
compileFromCpp.exe: CMakeFiles/compileFromCpp.dir/linklibs.rsp
compileFromCpp.exe: CMakeFiles/compileFromCpp.dir/objects1.rsp
compileFromCpp.exe: CMakeFiles/compileFromCpp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\unityProject\PhysicsPrj\Assets\compileFromCpp\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable compileFromCpp.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\compileFromCpp.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/compileFromCpp.dir/build: compileFromCpp.exe

.PHONY : CMakeFiles/compileFromCpp.dir/build

CMakeFiles/compileFromCpp.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\compileFromCpp.dir\cmake_clean.cmake
.PHONY : CMakeFiles/compileFromCpp.dir/clean

CMakeFiles/compileFromCpp.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\unityProject\PhysicsPrj\Assets\compileFromCpp D:\unityProject\PhysicsPrj\Assets\compileFromCpp D:\unityProject\PhysicsPrj\Assets\compileFromCpp\cmake-build-debug D:\unityProject\PhysicsPrj\Assets\compileFromCpp\cmake-build-debug D:\unityProject\PhysicsPrj\Assets\compileFromCpp\cmake-build-debug\CMakeFiles\compileFromCpp.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/compileFromCpp.dir/depend

