### This cmake file is modified from CMAKE's FindZLIB.cmake file.

set(_HTSLIB_SEARCHES)

find_path(HTSLIB_INCLUDE_DIR
	NAMES htslib/sam.h
	PATHS ${_HTSLIB_SEARCHES}
	HINTS ENV HTSLIB_ROOT
)

find_library(HTSLIB_LIBRARY
	NAMES libhts.so libhts.dylib libhts.a
	PATHS ${_HTSLIB_SEARCHES}
	HINTS ENV HTSLIB_ROOT)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HTSLIB REQUIRED_VARS HTSLIB_LIBRARY HTSLIB_INCLUDE_DIR)

if(HTSLIB_FOUND)
	get_filename_component(HTSLIB_LIBRARY_DIR ${HTSLIB_LIBRARY} DIRECTORY)


	find_path(LIBDEFLATE_INCLUDE_DIR
		NAMES libdeflate.h
		HINTS ${HTSLIB_INCLUDE_DIR}
	)

	find_library(LIBDEFLATE_LIBRARY
		NAMES libdeflate.so libdeflate.dylib libdeflate.a
		HINTS ${HTSLIB_LIBRARY_DIR}
	)
	
	if (NOT LIBDEFLATE_INCLUDE_DIR OR NOT LIBDEFLATE_LIBRARY)
		message(FATAL_ERROR "libdeflate not found. Required for the htslib.")
	endif()


	find_path(LIBLZMA_INCLUDE_DIR
		NAMES lzma.h
		HINTS ${HTSLIB_INCLUDE_DIR}
	)

	find_library(LIBLZMA_LIBRARY
		NAMES liblzma.so liblzma.dylib liblzma.a
		HINTS ${HTSLIB_LIBRARY_DIR}
	)
	
	if (NOT LIBLZMA_INCLUDE_DIR OR NOT LIBLZMA_LIBRARY)
		message(FATAL_ERROR "liblzma not found. Required for the htslib.")
	endif()


	set(HTSLIB_INCLUDE_DIRS ${LIBDEFLATE_INCLUDE_DIR} ${LIBLZMA_INCLUDE_DIR} ${HTSLIB_INCLUDE_DIR})
	set(HTSLIB_LIBRARIES ${LIBDEFLATE_LIBRARY} ${LIBLZMA_LIBRARY} ${HTSLIB_LIBRARY})
endif()
