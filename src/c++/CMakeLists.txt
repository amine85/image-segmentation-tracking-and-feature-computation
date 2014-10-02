cmake_minimum_required(VERSION 2.4)

if(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0005 OLD)
endif(COMMAND CMAKE_POLICY)


###################################
# set project's name
###################################
project( nano )

####################################################################
# Make sure that binaries aren't being built in the source directory
####################################################################
# if("${nano-project_BINARY_DIR}" MATCHES "^${nano-project_SOURCE_DIR}$")
#   message(FATAL_ERROR "nano-project must be built with separate source and binary directories.")
# endif("${nano-project_BINARY_DIR}" MATCHES "^${nano-project_SOURCE_DIR}$")

#########################################
# ADD NEW MODULES HERE
#########################################
option( BUILD_PREPROCESSING "Build Preprocessing Tools" ON )
option( BUILD_SEGMENTATION "Build Segmentation Tools" OFF )
option( BUILD_WELLCROPPING "Build WellCropping Tools" OFF )
option( BUILD_TRACKING "Build Tracking Tools" OFF )
option( BUILD_SCATTERPLOT "Build Scatterplot Tools" OFF )

include_directories( ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )
#########################################
# Building Installers
#########################################

set(INSTALL_BIN_DIR bin)
set(INSTALL_LIB_DIR lib)

set( EXE_DIR ${nano_BINARY_DIR}/exe )
set( LIB_DIR ${nano_BINARY_DIR}/libs )

# just debugging here
# message(${nano_SOURCE_DIR})
# message(${nano_BINARY_DIR})
# message(${EXE_DIR})
# message(${LIB_DIR})

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXE_DIR} )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_DIR} )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR} )

#########################################
# DEPENDENCIES
#########################################

# find_package( OpenCV REQUIRED )
# include_directories( ${OpenCV_INCLUDE_DIRS} )

FIND_PACKAGE(ITK REQUIRED)
INCLUDE( ${ITK_USE_FILE} )

# # FIND_PACKAGE(VTK)
# # 
# # IF(VTK_FOUND)
# #       INCLUDE(${VTK_USE_FILE})
# # ELSE(VTK_FOUND)
# #       MESSAGE(FATAL_ERROR "Cannot build without VTK")
# # ENDIF(VTK_FOUND)

if(MSVC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /openmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif (MSVC)



####################################################################################
# NOW WE CAN INCLUDE THE MODULES ACCORDING TO THE OPTIONS SELECTED ABOVE:
####################################################################################

if(BUILD_PREPROCESSING)
  add_subdirectory(preprocessing)
endif(BUILD_PREPROCESSING)

if(BUILD_SEGMENTATION)
  add_subdirectory(segmentation)
endif(BUILD_SEGMENTATION)

if(BUILD_WELLCROPPING)
  add_subdirectory(wellcropping)
endif(BUILD_WELLCROPPING)

if(BUILD_TRACKING)
  add_subdirectory(tracking)
endif(BUILD_TRACKING)

if(BUILD_SCATTERPLOT)
  add_subdirectory(scatterPlot)
endif(BUILD_SCATTERPLOT)



































