project( WellCropping )

cmake_minimum_required(VERSION 2.8)

# find_package( OpenCV REQUIRED )
# include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable( wellDetection wellDetection.cpp )
target_link_libraries( wellDetection ${OpenCV_LIBS} ${ITK_LIBRARIES})

# add_executable( wellDetectionv2 test.cpp )
# target_link_libraries( wellDetectionv2 ${OpenCV_LIBS} ${ITK_LIBRARIES})

# add_executable( wellDetectionv3 itkFFTWF_FFTTest.cxx itkFFTTest.h  )
# target_link_libraries( wellDetectionv3 ${OpenCV_LIBS} ${ITK_LIBRARIES})


# add_executable( wellDetection wellDetection_zeiss_v1.cpp )
# target_link_libraries( wellDetection ${OpenCV_LIBS} ${ITK_LIBRARIES})


# add_executable( wellDetection_zeiss_v1 wellDetection_zeiss_v1.cpp )
# target_link_libraries( wellDetection_zeiss_v1 ${OpenCV_LIBS} ${ITK_LIBRARIES})

# add_executable( wellDetection_single wellDetection_single.cpp )
# target_link_libraries( wellDetection_single ${OpenCV_LIBS} ${ITK_LIBRARIES})
