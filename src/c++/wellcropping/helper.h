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
#include "itkDiscreteGaussianImageFilter.h"
#include"itkSmoothingRecursiveGaussianImageFilter.h"

#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

#include "itkConfigure.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"

#if defined(ITK_USE_FFTWF) || defined(ITK_USE_FFTWD)
#include "itkFFTWInverseFFTImageFilter.h"
#include "itkFFTWForwardFFTImageFilter.h"
#endif

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkMultiThreader.h"
#include <omp.h>

#include "itkMultiplyImageFilter.h"



typedef unsigned short InputPixelType16;
typedef unsigned short OutputPixelType16;
typedef double FloatPixelType;


typedef itk::Image<InputPixelType16,3> InputImageType16;
typedef itk::Image<InputPixelType16,2> Input2DImageType16;
typedef itk::Image<InputPixelType16,3> OutputImageType16;
typedef itk::Image<InputPixelType16,2> Output2DImageType16;
typedef itk::Image<double,3> FloatImageType;
typedef itk::Image<double,2> Float2DImageType;

typedef itk::Image< unsigned char, 2>               ImageCHAR2;
typedef itk::Image< unsigned char, 3>               ImageCHAR3;

typedef itk::Image< unsigned short, 3>               ImageUS3;


typedef itk::Image< double, 1>               ImageF1;
typedef itk::Image< std::complex<double>, 1> ImageCF1;
typedef itk::Image< double, 2>               ImageF2;
typedef itk::Image< std::complex<double>, 2> ImageCF2;
typedef itk::Image< double, 3>               ImageF3;
typedef itk::Image< std::complex<double>, 3> ImageCF3;

typedef itk::Image< std::complex< double >, 2 >  ComplexImageTypeCF2;

typedef itk::MultiplyImageFilter <ImageF2, ImageF2 > MultiplyImageFilterTypeF2F2;
typedef itk::MultiplyImageFilter <ImageCF2, ImageCF2 > MultiplyImageFilterTypeCF2CF2;

void NCCFun( ImageF2::Pointer sliceImage, ImageF2::Pointer pattern_1, ImageF2::Pointer pattern_2, ImageF2::Pointer *outPat1, ImageF2::Pointer *outPat2 );


template <typename T>
typename T::Pointer readImage(const char *filename)
{
        //printf("Reading %s ... \n",filename);
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
        //printf("Done\n");
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

template <typename TIN, typename TOUT>
typename TOUT::Pointer FFT2( typename TIN::Pointer img )
{
    typedef typename itk::FFTWForwardFFTImageFilter< TIN > FFT2FilterType;
  // Real to complex pointer. This computes the forward FFT.
  typename FFT2FilterType::Pointer R2C = FFT2FilterType::New();
  R2C->SetInput( img );

  try
    {
    R2C->Update();
    }
  catch ( itk::ExceptionObject & err )
    {
                std::cerr << "ExceptionObject caught!" <<std::endl;
                std::cerr << err << std::endl;
    }
   return R2C->GetOutput();
}


template <typename TIN, typename TOUT>
typename TOUT::Pointer IFFT2( typename TIN::Pointer img )
{
    typedef typename itk::FFTWInverseFFTImageFilter< TIN, TOUT > IFFT2FilterType;
  // Real to complex pointer. This computes the forward FFT.
  typename IFFT2FilterType::Pointer R2C = IFFT2FilterType::New();
  R2C->SetInput( img );

  try
    {
    R2C->Update();
    }
  catch ( itk::ExceptionObject & err )
    {
                std::cerr << "ExceptionObject caught!" <<std::endl;
                std::cerr << err << std::endl;
    }
   return R2C->GetOutput();
}

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

template <typename TIN >
typename TIN::Pointer RotateImage180( typename TIN::Pointer imgIn )
{
    typename TIN::Pointer imgRot = GetITKImageOfSize<TIN>( imgIn->GetLargestPossibleRegion().GetSize() );
    
    typename TIN::PixelType * imgInPtr = imgIn->GetBufferPointer();
    typename TIN::PixelType * imgRotPtr = imgRot->GetBufferPointer();
    
    typename TIN::PixelType nr = imgIn->GetLargestPossibleRegion().GetSize()[0];
    typename TIN::PixelType nc = imgIn->GetLargestPossibleRegion().GetSize()[1];
//     typename TIN::PixelType ns = imgIn->GetLargestPossibleRegion().GetSize()[2];
    
    unsigned long long sz = (unsigned long long)nr*(unsigned long long)nc;
    
    for(size_t i_iter = 0; i_iter< sz; ++i_iter )
    {
        imgRotPtr[sz-i_iter-1] = imgInPtr[ i_iter ];
    }
    
    return imgRot;
}

template <typename TIN >
typename TIN::Pointer PadImagePost( typename TIN::Pointer imgIn, int xpad, int ypad  )
{
    typename TIN::PixelType nr = imgIn->GetLargestPossibleRegion().GetSize()[0];
    typename TIN::PixelType nc = imgIn->GetLargestPossibleRegion().GetSize()[1];
    
    typename TIN::PixelType nr_new = imgIn->GetLargestPossibleRegion().GetSize()[0]+xpad;
    typename TIN::PixelType nc_new = imgIn->GetLargestPossibleRegion().GetSize()[1]+ypad;
    
    typename TIN::SizeType size;
    size[0] = nr_new;       size[1] = nc_new;
    
    typename TIN::Pointer imgPad = GetITKImageOfSize<TIN>( size );
    
    typename TIN::PixelType * imgInPtr = imgIn->GetBufferPointer();
    typename TIN::PixelType * imgPadPtr = imgPad->GetBufferPointer();
    
    unsigned long long indx, indx_new;
    for( unsigned long long ii = 0; ii < nr_new; ++ ii ){
        for( unsigned long long jj = 0; jj < nc_new; ++ jj ){
            if( ii >= nr || jj >= nc ){
                indx_new = jj*nr_new + ii;
                imgPadPtr[ indx_new ] = 0;
            }
            else{
                indx = jj*nr + ii;
                indx_new = jj*nr_new + ii;
                imgPadPtr[ indx_new ] = imgInPtr[ indx ];
            }
        }
    }
    return imgPad;
}

template <typename TIN >
typename TIN::Pointer PadImagePrePost( typename TIN::Pointer imgIn, int xpad, int ypad  )
{
    typename TIN::PixelType nr = imgIn->GetLargestPossibleRegion().GetSize()[0];
    typename TIN::PixelType nc = imgIn->GetLargestPossibleRegion().GetSize()[1];
    
    typename TIN::PixelType nr_new = nr+2*xpad;
    typename TIN::PixelType nc_new = nc+2*ypad;
    
    typename TIN::SizeType size;
    size[0] = nr_new;       size[1] = nc_new;
    
    typename TIN::Pointer imgPad = GetITKImageOfSize<TIN>( size );
    
    typename TIN::PixelType * imgInPtr = imgIn->GetBufferPointer();
    typename TIN::PixelType * imgPadPtr = imgPad->GetBufferPointer();
    
    unsigned long long indx, indx_new;
    for( unsigned long long ii = 0; ii < nr_new; ++ ii ){
        for( unsigned long long jj = 0; jj < nc_new; ++ jj ){
            if( ii < xpad || ii >= xpad+nr || jj < ypad || jj>= ypad+nc ){
                indx_new = jj*nr_new + ii;
                imgPadPtr[ indx_new ] = 0;
            }
            else{
                indx = (jj-ypad)*nr + (ii-xpad);
                indx_new = jj*nr_new + ii;
                imgPadPtr[ indx_new ] = imgInPtr[ indx ];
            }
        }
    }
    return imgPad;
}

template <typename TIN >
typename TIN::Pointer PadImagePrePostRepeat( typename TIN::Pointer imgIn, int xpad, int ypad  )
{
    // Pad the image repeating the boundary pixels
    typename TIN::PixelType nr = imgIn->GetLargestPossibleRegion().GetSize()[0];
    typename TIN::PixelType nc = imgIn->GetLargestPossibleRegion().GetSize()[1];
    
    typename TIN::PixelType nr_new = nr+2*xpad;
    typename TIN::PixelType nc_new = nc+2*ypad;
    
    typename TIN::SizeType size;
    size[0] = nr_new;       size[1] = nc_new;
    
    typename TIN::Pointer imgPad = GetITKImageOfSize<TIN>( size );
    
    typename TIN::PixelType * imgInPtr = imgIn->GetBufferPointer();
    typename TIN::PixelType * imgPadPtr = imgPad->GetBufferPointer();
    
    unsigned long long indx, indx_new;
    for( unsigned long long ii = 0; ii < nr_new; ++ ii ){
        for( unsigned long long jj = 0; jj < nc_new; ++ jj ){
            indx_new = jj*nr_new + ii;
            long long ii_in = std::max( (long long)0, (long long)(ii-xpad) );
            ii_in = std::min( (long long)(ii_in), (long long)(nr-1) );
            long long jj_in = std::max( (long long)0, (long long)(jj-ypad) );
            jj_in = std::min( (long long)(jj_in), (long long)(nc-1) );
            
            indx = (jj_in)*nr + (ii_in);
            imgPadPtr[ indx_new ] = imgInPtr[ indx ];
        }
    }
    return imgPad;
}

template <typename TIN >
typename TIN::Pointer PowerImage( typename TIN::Pointer imgIn, int poww  )
{
    typename TIN::SizeType size = imgIn->GetLargestPossibleRegion().GetSize();
    
    typename TIN::Pointer imgOut = GetITKImageOfSize< TIN >( size );
    
    typename TIN::PixelType * imgInPtr = imgIn->GetBufferPointer();
    typename TIN::PixelType * imgOutPtr = imgOut->GetBufferPointer();
    
    unsigned long long indx;
    for( unsigned long long ii = 0; ii < size[0]; ++ ii ){
        for( unsigned long long jj = 0; jj < size[1]; ++ jj ){
            indx = ii + size[0]*jj;
            if( poww == 2 ){
                imgOutPtr[ indx ] = imgInPtr[ indx ]*imgInPtr[ indx ];
            }
        }
    }
    return imgOut;
}


template <typename TIN, typename TOUT >
typename TOUT::Pointer ExtractSlide( typename TIN::Pointer imgIn, int slice  )
{
    long long n1 = imgIn->GetLargestPossibleRegion().GetSize()[0];
    long long n2 = imgIn->GetLargestPossibleRegion().GetSize()[1];
    long long n3 = imgIn->GetLargestPossibleRegion().GetSize()[2];
    
    typename TOUT::SizeType size;
    size[0] = n1;       size[1] = n2;
    
    typename TOUT::Pointer imgOut = GetITKImageOfSize< TOUT >( size );
    
    typename TIN::PixelType * imgInPtr = imgIn->GetBufferPointer();
    typename TOUT::PixelType * imgOutPtr = imgOut->GetBufferPointer();
    
    unsigned long long indx, indx_org;
    unsigned long long kk = slice;
    for( unsigned long long ii = 0; ii < n1; ++ ii ){
        for( unsigned long long jj = 0; jj < n2; ++ jj ){
            indx = ii + n1*jj+n1*n2*kk;
            indx_org = ii + n1*jj;
            imgOutPtr[ indx_org ] = imgInPtr[ indx ];
        }
    }
    return imgOut;
}

template <typename TIN >
typename TIN::Pointer LocalSum( typename TIN::Pointer imgIn, typename TIN::SizeType sizeTem, int poww  )
{
    sizeTem[0] = (sizeTem[0]-1)/2;
    sizeTem[1] = (sizeTem[1]-1)/2;
    
    typename TIN::PixelType n1 = imgIn->GetLargestPossibleRegion().GetSize()[0];
    typename TIN::PixelType n2 = imgIn->GetLargestPossibleRegion().GetSize()[1];
    
    typename TIN::SizeType size;
    size[0] = n1 + 2*sizeTem[0];       
    size[1] = n2 + 2*sizeTem[1];       
    
    typename TIN::Pointer imgOut = GetITKImageOfSize< TIN >( size );
    
    typename TIN::PixelType * imgInPtr = imgIn->GetBufferPointer();
    typename TIN::PixelType * imgOutPtr = imgOut->GetBufferPointer();
    
    unsigned long long indx_in, indx;
    for( unsigned long long ii = 0; ii < size[0]; ++ii ){
        unsigned long long iimin = std::max( (unsigned long long)0, ii-2*sizeTem[0] );
        unsigned long long iimax = std::min( (unsigned long long)n1-1,ii );
        for( unsigned long long jj = 0; jj < size[1]; ++jj ){
            unsigned long long jjmin = std::max( (unsigned long long)0, jj-2*sizeTem[1] );
            unsigned long long jjmax = std::min( (unsigned long long)n2-1,jj );
            double temp = 0.0;
            for( unsigned long long ii_in = iimin; ii_in <= iimax; ++ii_in ){
                for( unsigned long long jj_in = jjmin; jj_in <= jjmax; ++jj_in ){
                    indx_in = ii_in + n1*jj_in;
                    temp += imgInPtr[ indx_in ];
                }
            }
            indx = ii + size[0]*jj;
            if( poww == 1 )
                imgOutPtr[ indx ] = temp;
            else
                imgOutPtr[ indx ] = temp*temp;
        }
    }
    return imgOut;
}

template <typename TIN >
typename TIN::Pointer LocalSum2( typename TIN::Pointer imgIn, typename TIN::SizeType sizeTem  )
{
    typename TIN::Pointer imgInPad = PadImagePrePost< TIN >( imgIn, sizeTem[0], sizeTem[1]);
    
    typename TIN::SizeType size = imgInPad->GetLargestPossibleRegion().GetSize();
    
//     typename TIN::Pointer imgCumSum1 = GetITKImageOfSize< TIN >( size );
    
    typename TIN::PixelType * imgInPadPtr = imgInPad->GetBufferPointer();
//     typename TIN::PixelType * imgCumSum1Ptr = imgCumSum1->GetBufferPointer();
    
    // Cumsum 1
    unsigned long long indx_j_post, indxPosi;
    unsigned long long indx_j_pre, indxNeg;
    unsigned long long indx;
    for( unsigned long long ii = 0; ii < size[0]; ++ ii ){
        for( unsigned long long jj = 0; jj < size[1]; ++ jj ){
            if( jj == 0)
                continue;
            indx_j_post = ii + size[0]*(jj);
            indx_j_pre = ii + size[0]*(jj-1);
            imgInPadPtr[ indx_j_post ] += imgInPadPtr[ indx_j_pre ];
        }
    }
    
    typename TIN::SizeType sizeC;
    sizeC[0] = size[0];
    sizeC[1] = size[1] - sizeTem[1] - 1;
    
    typename TIN::Pointer imgC = GetITKImageOfSize< TIN >( sizeC );
    typename TIN::PixelType * imgCPtr = imgC->GetBufferPointer();
    
    for( unsigned long long ii = 0; ii < sizeC[0]; ++ ii ){
        for( unsigned long long jj = 0; jj < sizeC[1]; ++ jj ){
            
            indx = ii + jj*sizeC[0];
            
            indxPosi = ii + (jj+sizeTem[1])*size[0];
            indxNeg = ii + (jj)*size[0];
            
            imgCPtr[ indx ] = imgInPadPtr[ indxPosi ] - imgInPadPtr[ indxNeg ];
        }
    }
    
    // Cumsum 2
    for( unsigned long long jj = 0; jj < sizeC[1]; ++ jj ){
        for( unsigned long long ii = 0; ii < sizeC[0]; ++ ii ){
            if( ii == 0)
                continue;
            indx_j_post = ii + sizeC[0]*(jj);
            indx_j_pre = ii-1 + sizeC[0]*(jj);
            imgCPtr[ indx_j_post ] += imgCPtr[ indx_j_pre ];
        }
    }
    
    typename TIN::SizeType sizeCN;
    sizeCN[0] = sizeC[0] - sizeTem[0] - 1;
    sizeCN[1] = sizeC[1];
    
    typename TIN::Pointer imgCN = GetITKImageOfSize< TIN >( sizeCN );
    typename TIN::PixelType * imgCNPtr = imgCN->GetBufferPointer();
    
    
        for( unsigned long long jj = 0; jj < sizeCN[1]; ++ jj ){
            for( unsigned long long ii = 0; ii < sizeCN[0]; ++ ii ){
            
            indx = ii + jj*sizeCN[0];
            
            indxPosi = ii+sizeTem[0] + (jj)*sizeC[0];
            indxNeg = ii + (jj)*sizeC[0];
            
            imgCNPtr[ indx ] = imgCPtr[ indxPosi ] - imgCPtr[ indxNeg ];
        }
    }
    return imgCN;
}

template <typename TIN >
double variance( typename TIN::Pointer imgIn )
{
    typename TIN::SizeType sizePattern_1 = imgIn->GetLargestPossibleRegion().GetSize();
    typename TIN::PixelType * pattern_1Ptr = imgIn->GetBufferPointer();
    
    double mean = 0, m2 = 0, variance = 0, n=0;
    double sizeeTemp = sizePattern_1[0]*sizePattern_1[1];
    for(size_t i=0; i<sizeeTemp; ++i) {
        n = n + 1;
        double delta = (double)pattern_1Ptr[i] - mean;
        mean += delta/n;
        m2 += delta*((double)pattern_1Ptr[i] - mean);
        variance = m2/(n - 1);
    }
    
    return variance;
}


#endif
