#ifndef __SEG_HELPERS_H
#define __SEG_HELPERS_H
#define MPICH_IGNORE_CXX_SEEK

// headers 
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
#include <iterator> 
#include <vector>
#include <algorithm>
#include <math.h>
#include <fstream> 
#include <cmath>
// itk 
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkLineIterator.h>
#include <itkMedianImageFilter.h>
#include <itkBinaryMedianImageFilter.h>

#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryMedianImageFilter.h>
#include <itkOtsuThresholdImageFilter.h> 
#include <itkConnectedComponentImageFilter.h>
#include <itkScalarConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include "itkLabelGeometryImageFilter.h"
#include <itkLaplacianRecursiveGaussianImageFilter.h>
#include "itkDiscreteGaussianImageFilter.h"
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkGrayscaleErodeImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>


#include <itkMaximumEntropyThresholdImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>


// vxl/vnl
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_sparse_matrix.h>
#include <vnl/vnl_diag_matrix.h>

#include <vnl/vnl_hungarian_algorithm.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_sparse_symmetric_eigensystem.h>

// eigen headers 
// #include <Eigen/Dense>
//#include "EigendSover.h"

// Macros
#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))
#define SIGN(a)  (((a)>=0)?1:-1)
#define EPS 1e-15


namespace seg_helpers{
    // Type Defs :
    // Pixel types
    typedef unsigned char BinPixelType;
    typedef unsigned short InputPixelType;
    typedef unsigned short LabelPixelType;
    typedef float FloatPixelType;
    // Image types
    typedef itk::Image<BinPixelType,3> BinImageType;
    typedef itk::Image<InputPixelType,3> InputImageType;
    typedef itk::Image<LabelPixelType,3> LabelImageType;
    typedef itk::Image<FloatPixelType,3> FloatImageType;

    typedef itk::Image<BinPixelType,2> BinImageType2D;
    typedef itk::Image<InputPixelType,2> InputImageType2D;
    typedef itk::Image<LabelPixelType,2> LabelImageType2D;
    typedef itk::Image<FloatPixelType,2> FloatImageType2D;

    // Iterator types
    typedef itk::ImageRegionIterator<BinImageType> BinImageIteratorType;
    typedef itk::ImageRegionIterator<InputImageType> InputImageIteratorType;
    typedef itk::ImageRegionIterator<LabelImageType> LabelImageIteratorType;
    typedef itk::ImageRegionIterator<FloatImageType> FloatImageIteratorType;
    
    // Filters
    typedef itk::DiscreteGaussianImageFilter<InputImageType,InputImageType> DiscreteGaussianType;
    
    // templated functions
    template <typename T>
        typename T::Pointer GetITKImageOfSize(typename T::SizeType size)
        {
            typename T::Pointer outputImage = T::New();
            typename T::IndexType start;
            start[0] = 0;  // size along X
            start[1] = 0;  // size along Y
            if(T::ImageDimension > 2)
                start[2] = 0;  // size along time 

            typename T::RegionType region;
            region.SetSize( size );
            region.SetIndex( start );

            outputImage->SetRegions( region );
            outputImage->Allocate();
            outputImage->FillBuffer(0);
            try
            {
                outputImage->Update();
            }
            catch(itk::ExceptionObject &err)
            {
                std::cerr << "ExceptionObject caught!" <<std::endl;
                std::cerr << err << std::endl;
                //return EXIT_FAILURE;
            }
            return outputImage;
        };

    template <typename T>
        typename T::Pointer readImage(const char *filename)
        {
            printf("Reading %s ... \n",filename);
            typedef typename itk::ImageFileReader<T> ReaderType;
            typename ReaderType::Pointer reader = ReaderType::New();

            ReaderType::GlobalWarningDisplayOff();
            reader->SetFileName(filename);
            try
            {
                reader->Update();
            }
            catch(itk::ExceptionObject &err)
            {
                std::cerr << "ExceptionObject caught!" <<std::endl;
                std::cerr << err << std::endl;
                //return EXIT_FAILURE;
            }
            printf("Done\n");
            return reader->GetOutput();
        };
    template <typename T>
        int writeImage(typename T::Pointer im, const char* filename)
        {
            printf("Writing %s ... \n",filename);
            typedef typename itk::ImageFileWriter<T> WriterType;

            typename WriterType::Pointer writer = WriterType::New();
            writer->SetFileName(filename);
            writer->SetInput(im);
            try
            {
                writer->Update();
            }
            catch(itk::ExceptionObject &err)
            {
                std::cerr << "ExceptionObject caught!" <<std::endl;
                std::cerr << err << std::endl;
                return EXIT_FAILURE;
            }
            return EXIT_SUCCESS;
        };
        
    // Functions
    bool KmeansClustering(const InputPixelType * Image, BinPixelType * LabelImage, unsigned int  N, int K, int iternum);

	void MinimumErrorTresholding( const InputPixelType * Image, BinPixelType * LabelImage, unsigned int N );

    void GetDistanceResponsev3(LabelImageType::Pointer inputLog2DImage,BinImageType::Pointer binaryImage, \
            FloatImageType::Pointer responseImage,unsigned int levels);

	void GetMaxLaplacianOfGaussianResponse( InputImageType::Pointer image, BinImageType::Pointer binaryImage, \
			FloatImageType::Pointer responseImage, unsigned int min_radius, unsigned int max_radius, unsigned int N );

    void GetSeedImagev3(FloatPixelType * ResponseImagePtr, BinPixelType * SeedImagePtr,\
            const unsigned int nr, const unsigned int nc, const unsigned int window_size);

    // Auxiliary functions
    void RemoveSmallComponents(BinImageType::Pointer BinImage, unsigned int min_volume); 
    void RemoveLargeComponents(BinImageType::Pointer BinImage, unsigned int max_volume); 
    
    // Helpers
    InputPixelType GetMax(const InputPixelType * X, unsigned int M);
    void BubbleSortAscend(vnl_vector<unsigned int> &X, vnl_vector<unsigned int> &indices);
}
#endif

































