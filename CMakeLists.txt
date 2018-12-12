cmake_minimum_required(VERSION 3.6.3 FATAL_ERROR)
project(demuxEMS CXX)

add_compile_options(-std=c++11)

find_package(ZLIB REQUIRED)
if (NOT ZLIB_FOUND)
	message(FATA_ERROR "zlib not found. Required for the SAMtools.")
endif()

find_package(Boost 1.61.0 REQUIRED COMPONENTS iostreams)
if (NOT Boost_FOUND OR NOT Boost_IOSTREAMS_FOUND)
	message(FATA_ERROR "boost not found. The iostreams library is required.")
endif()

# add_subdirectory(ext)
add_subdirectory(src)
