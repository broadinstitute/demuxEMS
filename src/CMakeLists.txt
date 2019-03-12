include_directories(${Boost_INCLUDE_DIRS})
include_directories(${HTSLIB_INCLUDE_DIRS})

message(STATUS, ${Boost_IOSTREAMS_LIBRARY})

file(GLOB sources BamAlignment.cpp DataCollector.cpp SEQstring.cpp SamParser.cpp BamWriter.cpp SamHeaderText.cpp VCFLoader.cpp)
file(GLOB headers BamAlignment.hpp BarcodeUtils.hpp DataCollector.hpp QUALstring.hpp SamHeaderText.hpp SufficientStatistics.hpp BamWriter.hpp CIGARstring.hpp MDstring.hpp SEQstring.hpp SamParser.hpp VCFLoader.hpp)

add_compile_options(-Wall -O3 -ffast-math)
add_library(demuxEMS_core ${sources} ${headers})

add_executable(scan_snp scan_snp.cpp)
target_link_libraries(scan_snp demuxEMS_core ${Boost_IOSTREAMS_LIBRARY} ${HTSLIB_LIBRARIES})

# add_executable(demuxEMS demuxEMS.cpp)
# target_link_libraries(demuxEMS demuxEMS_core ${Boost_IOSTREAMS_LIBRARY} ${HTSLIB_LIBRARIES})
# add_executable(debug debug.cpp)
# target_link_libraries(debug demuxEMS_core ${Boost_IOSTREAMS_LIBRARY} ${HTSLIB_LIBRARIES})
