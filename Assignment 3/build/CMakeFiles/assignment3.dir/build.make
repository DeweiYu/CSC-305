# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.24.2/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.24.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3/build

# Include any dependencies generated for this target.
include CMakeFiles/assignment3.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/assignment3.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/assignment3.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/assignment3.dir/flags.make

CMakeFiles/assignment3.dir/src/main.cpp.o: CMakeFiles/assignment3.dir/flags.make
CMakeFiles/assignment3.dir/src/main.cpp.o: /Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3/src/main.cpp
CMakeFiles/assignment3.dir/src/main.cpp.o: CMakeFiles/assignment3.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/assignment3.dir/src/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/assignment3.dir/src/main.cpp.o -MF CMakeFiles/assignment3.dir/src/main.cpp.o.d -o CMakeFiles/assignment3.dir/src/main.cpp.o -c /Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3/src/main.cpp

CMakeFiles/assignment3.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/assignment3.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3/src/main.cpp > CMakeFiles/assignment3.dir/src/main.cpp.i

CMakeFiles/assignment3.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/assignment3.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3/src/main.cpp -o CMakeFiles/assignment3.dir/src/main.cpp.s

# Object files for target assignment3
assignment3_OBJECTS = \
"CMakeFiles/assignment3.dir/src/main.cpp.o"

# External object files for target assignment3
assignment3_EXTERNAL_OBJECTS =

assignment3: CMakeFiles/assignment3.dir/src/main.cpp.o
assignment3: CMakeFiles/assignment3.dir/build.make
assignment3: CMakeFiles/assignment3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable assignment3"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/assignment3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/assignment3.dir/build: assignment3
.PHONY : CMakeFiles/assignment3.dir/build

CMakeFiles/assignment3.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/assignment3.dir/cmake_clean.cmake
.PHONY : CMakeFiles/assignment3.dir/clean

CMakeFiles/assignment3.dir/depend:
	cd /Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3 /Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3 /Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3/build /Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3/build /Users/liugongheng/Desktop/Trent/CSC305/CG-Spring-2023-main/Assignment_3/build/CMakeFiles/assignment3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/assignment3.dir/depend

