PROJECT(Preprocessing)

cmake_minimum_required(VERSION 2.0)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# adding a library
ADD_LIBRARY(helper helpers.cpp helpers.h)
TARGET_LINK_LIBRARIES(helper ${ITK_USE_LIBRARIES} )

# image stacking: 
ADD_EXECUTABLE(image_to_stack image_to_stack.cpp )
TARGET_LINK_LIBRARIES(image_to_stack helper ${ITK_LIBRARIES} )

# linear unmixing:
ADD_EXECUTABLE(unmix16 unmix16.cpp )
TARGET_LINK_LIBRARIES(unmix16 helper ${ITK_LIBRARIES} )

# background subtraction:
ADD_EXECUTABLE(background_subtraction background_subtraction.cpp )
TARGET_LINK_LIBRARIES(background_subtraction helper ${ITK_LIBRARIES} )

