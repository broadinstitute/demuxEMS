include_directories(${Boost_INCLUDE_DIRS})
include_directories(${HTSLIB_INCLUDE_DIRS})

file(GLOB sources *.cpp)
file(GLOB headers *.hpp)
list(REMOVE_ITEM sources scan_snp.cpp demuxEMS.cpp debug.cpp)
message(STATUS, ${sources})

add_compile_options(-Wall -O3 -ffast-math)
add_library(demuxEMS_core ${sources} ${headers})


# add_executable(demuxEMS demuxEMS.cpp)
# target_link_libraries(demuxEMS demuxEMS_core ${Boost_IOSTREAMS_LIBRARY} ${HTSLIB_LIBRARIES})
# add_executable(scan_snp scan_snp.cpp)
# target_link_libraries(scan_snp demuxEMS_core ${Boost_IOSTREAMS_LIBRARY} ${HTSLIB_LIBRARIES})
add_executable(debug debug.cpp)
target_link_libraries(debug demuxEMS_core ${Boost_IOSTREAMS_LIBRARY} ${HTSLIB_LIBRARIES})



# configure_file(utils.py utils.py COPYONLY)
# configure_file(PROBer-single-batch-estimate PROBer-single-batch-estimate COPYONLY)

# configure_file(PROBer PROBer COPYONLY)
# configure_file(PROBer-generateVariationPlot PROBer-generateVariationPlot COPYONLY)
 
# install(TARGETS PROBer-build-reference PROBer-parse-alignments PROBer-run-em  PROBer-simulate-reads PROBer-analyze-iCLIP PROBer-sample-iCLIP PROBer-run-gibbs PROBer-bootstrap DESTINATION bin)
# install(PROGRAMS PROBer PROBer-generateVariationPlot DESTINATION bin)
