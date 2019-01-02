cmake_minimum_required(VERSION 3.6.3 FATAL_ERROR)
project(demuxEMS CXX)

add_compile_options(-std=c++11)

find_package(Boost 1.67.0 REQUIRED COMPONENTS iostreams)
if (NOT Boost_FOUND OR NOT Boost_IOSTREAMS_FOUND)
	message(FATAL_ERROR "Boost is not found. The iostreams library is required.")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

find_package(HTSLIB REQUIRED)
if (NOT HTSLIB_FOUND)
	message(FATAL_ERROR "htslib is not found. Required for parsing BAM files.")
endif()

add_subdirectory(src)
