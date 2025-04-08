cmake_minimum_required(VERSION 3.16)

project(Landslide)

if ( CMAKE_SYSTEM_NAME MATCHES "Windows" )
    set( my_std_pre "/std:" )
else()
    set( my_std_pre "-std=" )
endif()
 
set( basic_cxx17 "c++17" )
set( str_cxx17 "${my_std_pre}${basic_cxx17}" )
 
include( CheckCXXCompilerFlag )
check_cxx_compiler_flag( "${str_cxx17}" _cpp_17_flag_supported )
if ( _cpp_17_flag_supported )
    set( CMAKE_CXX_STANDARD 17 )
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(UNIX)
	find_package(GDAL REQUIRED)
endif(UNIX)

if(MSVC)
	set(GDAL_INCLUDE_DIRS C:/Users/MPCBBishop/anaconda3/Lib/site-packages/osgeo/include/gdal)
	set(GDAL_LIBRARIES C:/Users/MPCBBishop/anaconda3/Lib/site-packages/osgeo/lib/gdal_i.lib)
endif(MSVC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

include_directories(${GDAL_INCLUDE_DIRS})
link_libraries(${GDAL_LIBRARIES})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${PROJECT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROJECT_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
add_subdirectory(src)