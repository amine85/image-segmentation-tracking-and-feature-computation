
#include<iostream>

#include"itkImageFileReader.h"
#include"itkImageFileWriter.h"
#include"itkImage.h"
#include"omp.h"
#include "itkMultiThreader.h"

typedef itk::Image< unsigned short,  3> ImageType_16bit;

template <typename T>
typename T::Pointer readImage(const char* filename)
{
// 	printf("Reading %s ... \n",filename);
	std::cout << std::endl << "Reading ... " << filename;
	typedef typename itk::ImageFileReader<T> ReaderType;

	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);
	try
	{
		reader->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	itk::Size<3> inputImageSize = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
	std::cout<<" done: Image size: "<<inputImageSize[0]<<"x"<<inputImageSize[1]<<"x"<<inputImageSize[2] <<std::flush;
	return reader->GetOutput();
}

template <typename T>
int writeImage(typename T::Pointer im, const char* filename)
{
// 	printf("Writing %s ... \n",filename);
	std::cout << std::endl << "Writing ... " << filename;
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
	itk::Size<3> inputImageSize = im->GetLargestPossibleRegion().GetSize();
	std::cout<<" done: Image size: "<<inputImageSize[0]<<"x"<<inputImageSize[1]<<"x"<<inputImageSize[2] <<std::flush;
	return EXIT_SUCCESS;
}


template <typename T>
typename T::Pointer createImage( itk::Size<3> size )
{
  	typename T::Pointer imageLabelMontage = T::New();
  	typename T::PointType originz;
  	originz[0] = 0;
  	originz[1] = 0;
  	originz[2] = 0;
  	imageLabelMontage->SetOrigin( originz );
  	typename T::IndexType indexStich;
  	indexStich.Fill(0);
  	typename T::RegionType regionz;
  	regionz.SetSize ( size  );
  	regionz.SetIndex( indexStich );
  	imageLabelMontage->SetRegions( regionz );
  	imageLabelMontage->Allocate();
  	imageLabelMontage->FillBuffer(0);
  	try
  	{
  		imageLabelMontage->Allocate();
  	}
  	catch(itk::ExceptionObject &err)
  	{
  		std::cerr << "ExceptionObject caught!" <<std::endl;
  		std::cerr << err << std::endl;
  	}
  	return imageLabelMontage;
}


int main( int argc, char* argv [] ) {
  if( argc < 5 ){
    std::cerr << "Need to specify two input parameters. <Input image name ch1 (x-axis)>, <Input image name ch2 (y-axis)>, <output image name (scatter plot)>, <OmpNumThreads>." << std::endl;
    exit(1);
  }
  std::cout << std::endl << "Input filenames";
  for( int i=0;i<argc;++i )
    std::cout << i << " " << argv[i] << std::endl;
  
  // Read the image in
  ImageType_16bit::Pointer InputImageCH1 = readImage< ImageType_16bit >( argv[1] );
  ImageType_16bit::Pointer InputImageCH2 = readImage< ImageType_16bit >( argv[2] );
  
  int numThreads = atoi(argv[4]);
        
    // Threads
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);
    omp_set_num_threads(numThreads);
  
  // Compare sizes
  itk::Size<3> SizeInputImageCH1 = InputImageCH1->GetLargestPossibleRegion().GetSize();
  itk::Size<3> SizeInputImageCH2 = InputImageCH2->GetLargestPossibleRegion().GetSize();
  
  if( (SizeInputImageCH1[0]!=SizeInputImageCH2[0]) || (SizeInputImageCH1[1]!=SizeInputImageCH2[1]) || (SizeInputImageCH1[2]!=SizeInputImageCH2[2]) ){
    std::cerr << "The image sizes are not the same." << std::endl;
    exit(1);
  }    
  
  // Create output image
  itk::Size<3> SizeOut;
  SizeOut[0] = 5000;
  SizeOut[1] = 5000;
  SizeOut[2] = 1;
  ImageType_16bit::Pointer OutputImage = createImage< ImageType_16bit >( SizeOut );
  
  
  unsigned long long sizeXY = SizeInputImageCH1[1] * SizeInputImageCH1[0];
  unsigned long long sizeX = SizeInputImageCH1[0];
  
  ImageType_16bit::PixelType * ArrayInputImageCH1 = InputImageCH1->GetBufferPointer();
  ImageType_16bit::PixelType * ArrayInputImageCH2 = InputImageCH2->GetBufferPointer();
  ImageType_16bit::PixelType * ArrayOutputImage = OutputImage->GetBufferPointer();
  
#pragma omp parallel for
  for(int itx = 0; itx < SizeInputImageCH1[0]; ++itx )
  {
    for(int ity = 0; ity < SizeInputImageCH1[1]; ++ity )
    {
      for(int itz = 0; itz < SizeInputImageCH1[2]; ++itz )
      {
	unsigned long long offsetIn = (itz*sizeXY)+(ity*sizeX)+itx;
	
	ImageType_16bit::PixelType pixelCh1 = ArrayInputImageCH1[ offsetIn ];
	ImageType_16bit::PixelType pixelCh2 = ArrayInputImageCH2[ offsetIn ];
	
	if( (pixelCh1>=5000) || (pixelCh2>=5000) )
	  continue;
	
	unsigned long long offsetOut = (pixelCh2*SizeOut[0])+pixelCh1;
	ArrayOutputImage[ offsetOut ] = ArrayOutputImage[ offsetOut ] + 1;
      }
    }
  }
  
  writeImage< ImageType_16bit >( OutputImage, argv[3] );
  std::cout << std::endl;
  return 0;
}





