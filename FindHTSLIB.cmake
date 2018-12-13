### This cmake file is modified from CMAKE's FindZLIB.cmake file.

set(_HTSLIB_SEARCHES)

find_path(HTSLIB_INCLUDE_DIR
	NAMES htslib/sam.h
	PATHS ${_HTSLIB_SEARCHES}
	HINTS ENV HTSLIB_ROOT
)

find_library(HTSLIB_LIBRARY
	NAMES libhts.a
	PATHS ${_HTSLIB_SEARCHES}
	HINTS ENV HTSLIB_ROOT)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HTSLIB REQUIRED_VARS HTSLIB_LIBRARY HTSLIB_INCLUDE_DIR)

if(HTSLIB_FOUND)
	find_path(ZLIB_INCLUDE_DIR
		NAMES zlib.h
		HINTS ${HTSLIB_INCLUDE_DIR}
	)

	get_filename_component(HTSLIB_LIBRARY_DIR ${HTSLIB_LIBRARY} DIRECTORY)

	find_library(ZLIB_LIBRARY
		NAMES libz.a
		HINTS ${HTSLIB_LIBRARY_DIR}
	)
	
	if (NOT ZLIB_INCLUDE_DIR OR NOT ZLIB_LIBRARY)
		message(FATAL_ERROR "zlib not found. Required for the htslib.")
	endif()

	set(HTSLIB_INCLUDE_DIRS ${ZLIB_INCLUDE_DIR} ${HTSLIB_INCLUDE_DIR})
	set(HTSLIB_LIBRARIES ${ZLIB_LIBRARY} ${HTSLIB_LIBRARY})
endif()
