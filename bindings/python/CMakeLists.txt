# Copyright (c) 2025 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

cmake_minimum_required(VERSION 3.7)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# define the Python project
project(unopy)

include_directories("${CMAKE_SOURCE_DIR}/../../uno")

set(UNO_DIRECTORY ${CMAKE_SOURCE_DIR}/../../uno)
set(SOURCE_FILES ${UNO_DIRECTORY}/options/Options.cpp)
set(HEADER_FILES
	${UNO_DIRECTORY}/linear_algebra/SparseVector.hpp
	${UNO_DIRECTORY}/linear_algebra/RectangularMatrix.hpp
	${UNO_DIRECTORY}/linear_algebra/SymmetricMatrix.hpp
	${UNO_DIRECTORY}/options/Options.hpp)
set(PYTHON_FILES modules/SparseVector.cpp modules/RectangularMatrix.cpp modules/SymmetricMatrix.cpp modules/Options.cpp unopy.cpp)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 REQUIRED)

pybind11_add_module(unopy 
	${SOURCE_FILES}
	${HEADER_FILES}
	${PYTHON_FILES}
)