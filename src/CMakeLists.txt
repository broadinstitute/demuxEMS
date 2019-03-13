include_directories(${Boost_INCLUDE_DIRS})
include_directories(${HTSLIB_INCLUDE_DIRS})

file(GLOB sources *.cpp)
file(GLOB headers *.hpp *.h)
remove(sources scan_snp.cpp)

add_compile_options(-Wall -O3 -ffast-math)
add_library(scan_snp_core ${sources} ${headers})

add_executable(scan_snp scan_snp.cpp)
target_link_libraries(scan_snp scan_snp_core ${Boost_IOSTREAMS_LIBRARY} ${HTSLIB_LIBRARIES})
