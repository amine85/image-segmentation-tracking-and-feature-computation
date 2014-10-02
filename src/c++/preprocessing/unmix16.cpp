#include "helpers.h"
#include "omp.h"
#include "itkMultiThreader.h"
#include "itkMaximumEntropyThresholdImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#define EPSILON 1e-15
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
	if(argc <7)
	{
		std::cout<<"Usage: unmix16 <InputImageFileName1> <InputImageFileName2> <OutputImageFileName1> <OutputImageFileName2> <MixingMatrix(txt)> <NewScatterImage> <OmpNumThreads>\n";
		return 0;
	}
	std::string fnameNewScatter = argv[6];
	
	int numThreads = atoi(argv[7]);
    // Threads
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);
    omp_set_num_threads(numThreads);
	
	std::string fname01 = argv[1];
	std::string fname02 = argv[2];
	std::cout<<"Input channel 1:"<<fname01<<std::endl;
 	std::cout<<"Input channel 2:"<<fname02<<std::endl;

	InputImageType16::Pointer inputImage01 = readImage<InputImageType16>(fname01.c_str());
	InputImageType16::Pointer inputImage02 = readImage<InputImageType16>(fname02.c_str());

	/***************************************************************************************************************/
	// Estimate the fingerprint(mixing) matrix 
	/***************************************************************************************************************/
	// initialize the mixing matrix //
	std::vector< std::vector< float > > mixingMatrix;
	mixingMatrix.resize(2);
	mixingMatrix[0].resize(2);
	mixingMatrix[1].resize(2);
	
	mixingMatrix[0][0] = 1;
	mixingMatrix[0][1] = 0;
	mixingMatrix[1][0] = 0;
	mixingMatrix[1][1] = 1;

	// initialize the updated matrix //
	std::vector< std::vector< float > > updatedMixingMatrix;
	updatedMixingMatrix.resize(2);
	updatedMixingMatrix[0].resize(2);
	updatedMixingMatrix[1].resize(2);
	
	FILE *fp = fopen(argv[5],"r");
	if(fp==NULL)
	{
	  printf("Could not read mixing matrix file\n");
	  return 1;
	}
	else
	{
	    for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
		    float num = -1;
		    fscanf(fp,"%f",&num);
		    mixingMatrix[i][j] = num;
		    std::cout<<mixingMatrix[i][j]<<"\t";
		    
		  }
		  std::cout<<"\n";
	    }
	}
	fclose(fp);
	
	
	size_t numRows=inputImage01->GetLargestPossibleRegion().GetSize()[0];
	size_t numCols=inputImage01->GetLargestPossibleRegion().GetSize()[1];
	size_t numTimes=inputImage01->GetLargestPossibleRegion().GetSize()[2];
	
	/***************************************************************************************************************/
	// unmix by matrix inversion	
	/***************************************************************************************************************/  
	
	float a11 =  mixingMatrix[0][0];
	float a21 =  mixingMatrix[1][0];
	float a12 =  mixingMatrix[0][1];
	float a22 =  mixingMatrix[1][1];
	
	float norm = sqrt(a11*a11 + a21*a21);
	a11 /= (norm+EPSILON);
	a21 /= (norm+EPSILON);
	norm = sqrt(a12*a12 + a22*a22);
	a12 /= (norm+EPSILON);
	a22 /= (norm+EPSILON);
	
	float detM = a11*a22 - a21*a12;
	updatedMixingMatrix[0][0] =  a22/(detM+EPSILON);
	updatedMixingMatrix[1][0] = -a21/(detM+EPSILON);
	updatedMixingMatrix[0][1] = -a12/(detM+EPSILON);
	updatedMixingMatrix[1][1] =  a11/(detM+EPSILON);
	
	std::cout<<"Normalized Matrix:\n";
	std::cout<<  a11 << "\t"<< a12 <<"\n";
	std::cout<<  a21 << "\t"<< a22 <<"\n";  
	std::cout<< std::endl;
	std::cout<<"Inverted Matrix:\n";
	std::cout<<  updatedMixingMatrix[0][0] << "\t"<< updatedMixingMatrix[0][1]<<"\n";
	std::cout<<  updatedMixingMatrix[1][0] << "\t"<< updatedMixingMatrix[1][1]<<"\n";  
	std::cout<< std::endl;
	

// 	// allocate output image space //
// 	InputImageType16::Pointer outputImage01 = InputImageType16::New();
// 
// 	InputImageType16::IndexType start;
// 	start[0] = 0;  // size along X
// 	start[1] = 0;  // size along Y
// 	start[2] = 0;  // size along time 
// 	
// 	InputImageType16::SizeType  size;
// 	size[0] = numRows;  // size along X
// 	size[1] = numCols;  // size along Y
// 	size[2] = numTimes;  // size along time
// 
// 	InputImageType16::RegionType region;
// 	region.SetSize( size );
// 	region.SetIndex( start );
// 
// 	outputImage01->SetRegions( region );
// 	outputImage01->Allocate();
// 	outputImage01->FillBuffer(0);
// 	outputImage01->Update();
// 
// 	InputImageType16::Pointer outputImage02 = InputImageType16::New();
// 	outputImage02->SetRegions( region );
// 	outputImage02->Allocate();
// 	outputImage02->FillBuffer(0);
// 	outputImage02->Update();
// 
	InputPixelType16 * inputImage01Ptr = inputImage01->GetBufferPointer();
	InputPixelType16 * inputImage02Ptr = inputImage02->GetBufferPointer();
// 	//copy the output image into the ITK image
// 	typedef itk::ImageRegionIteratorWithIndex< InputImageType16 > IteratorType2D16;
// 	InputPixelType16 * outputImage01Ptr = outputImage01->GetBufferPointer();
// 	InputPixelType16 * outputImage02Ptr = outputImage02->GetBufferPointer();
// 	
	std::cout<< "I am  starting inversion"<<std::endl;
	
	// declare loop variables
	size_t nPixels = (numRows*numCols*numTimes);
//	float outVal01 = 0.0;
//	float outVal02 = 0.0;
	
//	float inVal01 = 0.0;
//	float inVal02 = 0.0;

  float a_00 = updatedMixingMatrix[0][0];
  float a_10 = updatedMixingMatrix[1][0];
  float a_01 = updatedMixingMatrix[0][1];
  float a_11 = updatedMixingMatrix[1][1];

	
	//#pragma omp  parallel for schedule(dynamic, 1)
	for(size_t index=0; index<nPixels; ++index)
	{	
	  
// 		float outVal01 = (inputImage01Ptr[index]*updatedMixingMatrix[0][0] +  inputImage02Ptr[index]*updatedMixingMatrix[0][1]);
// 		float outVal02 = (inputImage01Ptr[index]*updatedMixingMatrix[1][0] +  inputImage02Ptr[index]*updatedMixingMatrix[1][1]);

		float inVal01 = (float) inputImage01Ptr[index];
		float inVal02 = (float) inputImage02Ptr[index];
		
    inputImage01Ptr[index] = (InputPixelType16)std::max( 0.0f, (inVal01*a_00 +  inVal02*a_01) );
    inputImage02Ptr[index] = (InputPixelType16)std::max( 0.0f, (inVal01*a_10 +  inVal02*a_11) );
	}
	
	std::cout<< "I am done ....................."<<std::endl;
	std::string outfName01 = argv[3];
	std::string outfName02 = argv[4];
	std::cout<< "I am writing: "<<outfName01<<std::endl<<std::flush;
	std::cout<< "I am writing: "<<outfName02<<std::endl<<std::flush;
	
//	writeImage<InputImageType16>(outputImage01,outfName01.c_str());
//	writeImage<InputImageType16>(outputImage02,outfName02.c_str());
	
	writeImage<InputImageType16>(inputImage01,outfName01.c_str());
	writeImage<InputImageType16>(inputImage02,outfName02.c_str());
    
    
    { // COMPUTE THE NEW MIXING MATRIX
    
        // Compare sizes
        itk::Size<3> SizeInputImageCH1 = inputImage01->GetLargestPossibleRegion().GetSize();
        itk::Size<3> SizeInputImageCH2 = inputImage02->GetLargestPossibleRegion().GetSize();
        
        // Create output image
        itk::Size<3> SizeOut;
        SizeOut[0] = 5000;
        SizeOut[1] = 5000;
        SizeOut[2] = 1;
        InputImageType16::Pointer OutputImage = createImage< InputImageType16 >( SizeOut );
        
        
        unsigned long long sizeXY = SizeInputImageCH1[1] * SizeInputImageCH1[0];
        unsigned long long sizeX = SizeInputImageCH1[0];
        
        InputImageType16::PixelType * ArrayInputImageCH1 = inputImage01->GetBufferPointer();
        InputImageType16::PixelType * ArrayInputImageCH2 = inputImage02->GetBufferPointer();
        InputImageType16::PixelType * ArrayOutputImage = OutputImage->GetBufferPointer();
        
//         #pragma omp parallel for
        nPixels = SizeInputImageCH1[0]*SizeInputImageCH1[1]*SizeInputImageCH1[2];
        for(size_t index=0; index<nPixels; ++index)
        {
            InputImageType16::PixelType pixelCh1 = ArrayInputImageCH1[ index ];
            InputImageType16::PixelType pixelCh2 = ArrayInputImageCH2[ index ];
            
            if( (pixelCh1>=5000) || (pixelCh2>=5000) )
            continue;
            
            unsigned long long offsetOut = (pixelCh2*SizeOut[0])+pixelCh1;
            ArrayOutputImage[ offsetOut ] = ArrayOutputImage[ offsetOut ] + 1;
        }
        
        writeImage< InputImageType16 >( OutputImage, fnameNewScatter.c_str() );
    }
    std::cout << std::endl;
	
}
