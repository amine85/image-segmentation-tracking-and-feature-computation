cmake_minimum_required(VERSION 2.8)

PROJECT( ScatterPlot )

ADD_EXECUTABLE(ScatterPlot ScatterPlot.cpp )
TARGET_LINK_LIBRARIES( ScatterPlot ${ITK_LIBRARIES} ) 
