
#include "helpers.h"
#include "omp.h"
#include "itkMultiThreader.h"

using namespace helpers;
#if defined(_MSC_VER)
#pragma warning(disable: 4996)
#endif
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

double start_t,end_t,diff_t;

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
}




int main(int argc, char **argv)
{
  	std::cout<<"number of arguments is: "<<argc<<std::endl;
	if(argc <5)
	{
		std::cout<<"Usage: image_to_stack  <InputImageFileNames> <OutputDirectory> <OutputFileName> <OmpNumThreads>\n";
		return 0;
	}
	
  
	std::string inputImageFilename = argv[1];  
	std::string outputDirectory = argv[2];  
	std::string outputFileName = argv[3];
    int numThreads = atoi(argv[4]);
        
    // Threads
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);
    omp_set_num_threads(numThreads);
	
	std::string line;
	// Read Input Filenames:
	std::vector< std::string > infnames;
	
	std::ifstream myfile (inputImageFilename.c_str());
	if (myfile.is_open())
	  {
		while ( myfile.good() )
		{
		  std::getline (myfile,line);
		  if(!line.empty())
		  {
			infnames.push_back(line);
		  }
		  //std::cout << line << endl;
		}
		myfile.close();
	  }
	else std::cout << "Unable to open input files"; 
	

	// Read Images:
	std::vector< InputImageType16::Pointer > inputImages;

	int num_time_points = infnames.size();
	std::cout<< " number of time points :"<<infnames.size()<<std::endl;

	inputImages.resize( num_time_points );
// #pragma omp  parallel for schedule(dynamic, 1)
	for(int t =0; t<num_time_points; t++)
	{
		
		helpers::InputImageType16::Pointer currImage = readImage<helpers::InputImageType16>(infnames[t].c_str());	
// 		inputImages.push_back(currImage);
		inputImages[t] = currImage;
		
	}
	
	// Some images have 2049 pixels, need to remove one pixel
	int firsRow = 0;
    int lastRow = 0;
    int firsCol = 0;
    int lastCol = 0;
    
    helpers::InputPixelType16 * firstImagePtr = inputImages[0]->GetBufferPointer();
    InputImageType16::SizeType sizeInput = inputImages.at(0)->GetLargestPossibleRegion().GetSize();
//     for( long long jj = 0; jj<sizeInput[1]; ++jj ){
    long long ii = 0;
    long long jj = 0;
    long long indx = 0;
    
    jj = 0;
    for( ii = 0; ii<sizeInput[0]; ++ii ){
        indx = ii+jj*sizeInput[0];
        firsRow += firstImagePtr[ indx ];
    }
    jj = (int)sizeInput[1]-1;
    for( ii = 0; ii<sizeInput[0]; ++ii ){
        indx = ii+jj*sizeInput[0];
        lastRow += firstImagePtr[ indx ];
    }
    ii = 0;
    for( jj = 0; jj<sizeInput[1]; ++jj ){
        indx = ii+jj*sizeInput[0];
        firsCol += firstImagePtr[ indx ];
    }
    ii = (int)sizeInput[0]-1;
    for( jj = 0; jj<sizeInput[1]; ++jj ){
        indx = ii+jj*sizeInput[0];
        lastCol += firstImagePtr[ indx ];
    }
    firsRow = firsRow>0?0:1;
    lastRow = lastRow>0?0:1;
    firsCol = firsCol>0?0:1;
    lastCol = lastCol>0?0:1;
    
//     if( firsRow>0 && lastRow>0 ){
//         std::cout << std::endl << "Error, one of the sizes is more than 2049 row" << std::flush;
//         return 1;
//     }
// 
//     if( firsCol>0 && lastCol>0 ){
//         std::cout << std::endl << "Error, one of the sizes is more than 2049 col" << std::flush;
//         return 1;
//     }

	// Create stack image
	InputImageType16::Pointer stackImage = InputImageType16::New();
	InputImageType16::RegionType region;
	InputImageType16::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	InputImageType16::SizeType sizeOut;
    sizeOut[0] = sizeInput[0];
    sizeOut[1] = sizeInput[1];
// 	sizeOut[0] = (long long)sizeInput[0] - firsCol-lastCol;
// 	sizeOut[1] = (long long)sizeInput[1] - firsRow-lastRow;
	sizeOut[2] = num_time_points;
	region.SetSize(sizeOut);
	region.SetIndex(start);
	stackImage->SetRegions(region);
	stackImage->Allocate();
	stackImage->FillBuffer(0);
	try
	{
		stackImage->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}

	// use pointers //
	helpers::InputPixelType16 * stackImagePtr = stackImage->GetBufferPointer();
    
    double averageValue = 0;
    size_t rc = sizeInput[0]*sizeInput[1];
	for(int t=0; t<num_time_points; t++)
	{
	    size_t offset = t*rc;
	    
	    helpers::InputPixelType16 * currImagePtr = inputImages[t]->GetBufferPointer();
	    
	    for(unsigned int i = 0; i <rc; ++i){
		 stackImagePtr[i+offset] = currImagePtr[i];
         averageValue += currImagePtr[i];
        }
	}
	averageValue /= (num_time_points*sizeInput[0]*sizeInput[1]);
    
    // Some images have an empty value in the border, so this will replace tha value by the average.
	for(int kk=0; kk<num_time_points; kk++){
        if( firsRow == 1 ){
            jj = 0;
            for( ii = 0; ii<sizeInput[0]; ++ii ){
                indx = ii+jj*sizeInput[0]+kk*rc;
                stackImagePtr[ indx ] = averageValue;
            }
        }
        if( lastRow == 1 ){
            jj = (int)sizeInput[1]-1;
            for( ii = 0; ii<sizeInput[0]; ++ii ){
                indx = ii+jj*sizeInput[0]+kk*rc;
                stackImagePtr[ indx ] = averageValue;
            }
        }
        if( firsCol == 1 ){
            ii = 0;
            for( jj = 0; jj<sizeInput[1]; ++jj ){
                indx = ii+jj*sizeInput[0]+kk*rc;
                stackImagePtr[ indx ] = averageValue;
            }
        }
        if( lastCol == 1 ){
            ii = (int)sizeInput[0]-1;
            for( jj = 0; jj<sizeInput[1]; ++jj ){
                indx = ii+jj*sizeInput[0]+kk*rc;
                stackImagePtr[ indx ] = averageValue;
            }
        }
    }
	
	// save the output stack:
	std::string savefName = outputDirectory + "/" + outputFileName;
	writeImage<InputImageType16>(stackImage,savefName.c_str());



	return 0;
}
