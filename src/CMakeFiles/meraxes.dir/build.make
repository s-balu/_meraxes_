# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_COMMAND = /fred/oz113/yqin/3rd_party/lib/python3.6/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /fred/oz113/yqin/3rd_party/lib/python3.6/site-packages/cmake/data/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /fred/oz025/yqin/bitbucket/meraxes/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /fred/oz025/yqin/bitbucket/meraxes/src

# Include any dependencies generated for this target.
include CMakeFiles/meraxes.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/meraxes.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/meraxes.dir/flags.make

CMakeFiles/meraxes.dir/core/meraxes.c.o: CMakeFiles/meraxes.dir/flags.make
CMakeFiles/meraxes.dir/core/meraxes.c.o: core/meraxes.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/fred/oz025/yqin/bitbucket/meraxes/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/meraxes.dir/core/meraxes.c.o"
	/apps/skylake/software/core/gcccore/7.3.0/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/meraxes.dir/core/meraxes.c.o   -c /fred/oz025/yqin/bitbucket/meraxes/src/core/meraxes.c

CMakeFiles/meraxes.dir/core/meraxes.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/meraxes.dir/core/meraxes.c.i"
	/apps/skylake/software/core/gcccore/7.3.0/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /fred/oz025/yqin/bitbucket/meraxes/src/core/meraxes.c > CMakeFiles/meraxes.dir/core/meraxes.c.i

CMakeFiles/meraxes.dir/core/meraxes.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/meraxes.dir/core/meraxes.c.s"
	/apps/skylake/software/core/gcccore/7.3.0/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /fred/oz025/yqin/bitbucket/meraxes/src/core/meraxes.c -o CMakeFiles/meraxes.dir/core/meraxes.c.s

# Object files for target meraxes
meraxes_OBJECTS = \
"CMakeFiles/meraxes.dir/core/meraxes.c.o"

# External object files for target meraxes
meraxes_EXTERNAL_OBJECTS =

bin/meraxes: CMakeFiles/meraxes.dir/core/meraxes.c.o
bin/meraxes: CMakeFiles/meraxes.dir/build.make
bin/meraxes: lib/libmeraxes.a
bin/meraxes: lib/libmlog.a
bin/meraxes: /apps/skylake/software/OpenMPI/3.0.0-GCC-7.3.0/lib/libmpi.so
bin/meraxes: /apps/skylake/software/compiler/gcc/7.3.0/gsl/2.4/lib/libgsl.so
bin/meraxes: /apps/skylake/software/compiler/gcc/7.3.0/gsl/2.4/lib/libgslcblas.so
bin/meraxes: /apps/skylake/software/mpi/gcc/7.3.0/openmpi/3.0.0/hdf5/1.10.1/lib/libhdf5.so
bin/meraxes: /apps/skylake/software/core/szip/2.1.1/lib/libsz.so
bin/meraxes: /usr/lib64/libz.so
bin/meraxes: /usr/lib64/libdl.so
bin/meraxes: /usr/lib64/libm.so
bin/meraxes: /usr/lib64/libpthread.so
bin/meraxes: /apps/skylake/software/mpi/gcc/7.3.0/openmpi/3.0.0/hdf5/1.10.1/lib/libhdf5_hl.so
bin/meraxes: /apps/skylake/software/mpi/gcc/7.3.0/openmpi/3.0.0/fftw/3.3.7/lib/libfftw3f_mpi.so
bin/meraxes: /apps/skylake/software/mpi/gcc/7.3.0/openmpi/3.0.0/fftw/3.3.7/lib/libfftw3f.so
bin/meraxes: CMakeFiles/meraxes.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/fred/oz025/yqin/bitbucket/meraxes/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable bin/meraxes"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/meraxes.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/meraxes.dir/build: bin/meraxes

.PHONY : CMakeFiles/meraxes.dir/build

CMakeFiles/meraxes.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/meraxes.dir/cmake_clean.cmake
.PHONY : CMakeFiles/meraxes.dir/clean

CMakeFiles/meraxes.dir/depend:
	cd /fred/oz025/yqin/bitbucket/meraxes/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /fred/oz025/yqin/bitbucket/meraxes/src /fred/oz025/yqin/bitbucket/meraxes/src /fred/oz025/yqin/bitbucket/meraxes/src /fred/oz025/yqin/bitbucket/meraxes/src /fred/oz025/yqin/bitbucket/meraxes/src/CMakeFiles/meraxes.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/meraxes.dir/depend
