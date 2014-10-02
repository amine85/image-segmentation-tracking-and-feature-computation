PROJECT(MixtureSegmentation)

cmake_minimum_required(VERSION 2.0)

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

SET(MIXTURE_SEGMENTATION_SRCS
    seg_helpers.cpp
)
SET(MIXTURE_SEGMENTATION_HDRS
    seg_helpers.h
)

ADD_LIBRARY(MixtureSegmentation ${MIXTURE_SEGMENTATION_SRCS} ${MIXTURE_SEGMENTATION_HDRS})

# kmeans clustering segmenation 
ADD_EXECUTABLE(mixture_segment main.cpp )
TARGET_LINK_LIBRARIES(mixture_segment MixtureSegmentation  ${ITK_LIBRARIES} )

# removing small components
#ADD_EXECUTABLE(clean_mixture_segmentation remove_small_components.cpp )
#TARGET_LINK_LIBRARIES(clean_mixture_segmentation MixtureSegmentation  ${ITK_LIBRARIES} )

# find seeds
#ADD_EXECUTABLE(detect_seeds detect_seeds.cpp )
#TARGET_LINK_LIBRARIES(detect_seeds MixtureSegmentation  ${ITK_LIBRARIES} )


# find number of cells and segment using spectral clustering
#ADD_EXECUTABLE(spectral_segment spectral_segment.cpp )
#TARGET_LINK_LIBRARIES(spectral_segment  MixtureSegmentation  ${ITK_LIBRARIES} )
