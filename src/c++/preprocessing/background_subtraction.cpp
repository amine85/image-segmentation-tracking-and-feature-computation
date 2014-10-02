#include "helpers.h"
#include "omp.h"
#include "itkMultiThreader.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"
using namespace helpers;
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

}
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
}


#define MAX_TIME 1000
bool file_exists(char *filename)
{
	FILE * fp = fopen(filename,"r");
	if(fp!=NULL)
	{
		fclose(fp);
		return true;
	}
	return false;
}


int main(int argc, char**argv)
{

	std::cout<<"number of arguments is: "<<argc<<std::endl;
	if(argc <5)
	{
		std::cout<<"Usage: background_subtraction <InputImageFileName> <OutputImageFileName> <GaussFilterSize>\n";
		return 0;
	}
	
        int numThreads = atoi(argv[4]);
        // Threads
        itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
        itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);
        omp_set_num_threads(numThreads);
	
	std::string input = argv[1];
	std::string outfname = argv[2];
	std::cout<<"Input filename:"<<input<<std::endl;
	

	InputImageType16::Pointer inputImage = readImage<InputImageType16>(input.c_str());

	InputImageType16::SizeType  inputImageSize = inputImage->GetLargestPossibleRegion().GetSize();
	size_t numRows = inputImageSize[0];  // size along X
	size_t numCols = inputImageSize[1];  // size along Y
	size_t numTimes =inputImageSize[2];  // size along time	

// 	InputImageType16::Pointer outputImage = InputImageType16::New();
// 
// 	InputImageType16::IndexType start1;
// 	start1[0] = 0;  // size along X
// 	start1[1] = 0;  // size along Y
// 	start1[2] = 0;  // size along time 
// 	
// 	
// 	InputImageType16::SizeType  size1;
// 	size1[0] = numRows;//inputImage->GetLargestPossibleRegion()->GetSize()[0];  // size along X
// 	size1[1] = numCols;//inputImage->GetLargestPossibleRegion()->GetSize()[1];  // size along Y
// 	size1[2] = numTimes;  // size along time
// 
// 	InputImageType16::RegionType region1;
// 	region1.SetSize( size1 );
// 	region1.SetIndex( start1 );
// 
// 	outputImage->SetRegions( region1 );
// 	outputImage->Allocate();
// 	outputImage->FillBuffer(0);
// 	outputImage->Update();

	InputPixelType16 * inputImagePtr = inputImage->GetBufferPointer();
// 	InputPixelType16 * outputImagePtr = outputImage->GetBufferPointer();

	#pragma omp  parallel for schedule(dynamic, 1)
	for(size_t t = 0; t<numTimes; ++t)
	{
		Input2DImageType16::Pointer input_image = Input2DImageType16::New();

		Input2DImageType16::IndexType start;
		start[0] = 0;  // size along X
		start[1] = 0;  // size along Y
		
		
		Input2DImageType16::SizeType  size;
		size[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];  // size along X
		size[1] = inputImage->GetLargestPossibleRegion().GetSize()[1];  // size along Y

		Input2DImageType16::RegionType region;
		region.SetSize( size );
		region.SetIndex( start );

		input_image->SetRegions( region );
		input_image->Allocate();
		input_image->FillBuffer(0);
		input_image->Update();

		InputPixelType16 * input2DImagePtr = input_image->GetBufferPointer();
		size_t offset = t*(numRows*numCols);
		for(size_t index = 0; index < (numRows*numCols); ++index)
		{
			input2DImagePtr[index] = inputImagePtr[index+offset];
		}
		
		Input2DImageType16::SpacingType spacing;
		spacing[0] = 1;
		spacing[1] = 1;
		input_image->SetSpacing(spacing);
		typedef itk::SmoothingRecursiveGaussianImageFilter<Input2DImageType16,Input2DImageType16> FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput(input_image);
		filter->SetSigma(atof(argv[3]));
		try
		{
			filter->Update();
		}
		catch(itk::ExceptionObject &err)
		{
			std::cerr << "ExceptionObject caught!" <<std::endl;
			std::cerr << err << std::endl;
		}
		Input2DImageType16::Pointer backg_image = filter->GetOutput();

		
		InputPixelType16 * backgImagePtr = backg_image->GetBufferPointer();
		
		
		float diffVal = 0.0;
		for(size_t index = 0; index < (numRows*numCols); ++index)
		{
			diffVal = (float)input2DImagePtr[index]-(float)backgImagePtr[index];
			if(diffVal>0)
				inputImagePtr[index+offset] = (InputPixelType16) diffVal;
			else
				inputImagePtr[index+offset] = 0;
		}		
		
// 		typedef itk::SubtractImageFilter <Input2DImageType16, Input2DImageType16, Float2DImageType> SubtractImageFilterType;
// 		SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New ();
// 		subtractFilter->SetInput1(input_image);
// 		subtractFilter->SetInput2(backg_image);
// 		try
// 		{
// 			subtractFilter->Update();
// 		}
// 		catch(itk::ExceptionObject &err)
// 		{
// 			std::cerr << "ExceptionObject caught!" <<std::endl;
// 			std::cerr << err << std::endl;
// 		}
// 		Float2DImageType::Pointer diff_image = subtractFilter->GetOutput();
// 		FloatPixelType * diffImagePtr = diff_image->GetBufferPointer();		
// 
// 
// 
// 		for(size_t index = 0; index < (numRows*numCols); ++index)
// 		{
// 			if(diffImagePtr[index]>0)
// 				outputImagePtr[index+offset] = (InputPixelType16) diffImagePtr[index];
// 			else
// 				outputImagePtr[index+offset] = 0;
// 		}
		

	}

	
	std::cout<<outfname<<std::endl;
// 	writeImage<InputImageType16>(outputImage, outfname.c_str());
	writeImage<InputImageType16>(inputImage, outfname.c_str());
	return 0;
}
