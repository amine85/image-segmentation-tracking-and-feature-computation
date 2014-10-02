#ifndef __HELPERS_H
#define __HELPERS_H

// c includes
#include<stdio.h>

// c++ inlcudes
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
#include <iterator> 
#include <vector>

// itk includes
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>
#include <itkSubtractImageFilter.h>



//typedefs
namespace helpers{
	
typedef unsigned short InputPixelType16;
typedef unsigned short OutputPixelType16;
typedef float FloatPixelType;


typedef itk::Image<InputPixelType16,3> InputImageType16;
typedef itk::Image<InputPixelType16,2> Input2DImageType16;
typedef itk::Image<InputPixelType16,3> OutputImageType16;
typedef itk::Image<InputPixelType16,2> Output2DImageType16;
typedef itk::Image<float,3> FloatImageType;
typedef itk::Image<float,2> Float2DImageType;


}

#endif