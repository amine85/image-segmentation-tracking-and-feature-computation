#include "preprocessor.h"
#include "omp.h"
#include "itkMultiThreader.h"

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

int main(int argc, char**argv)
{

    std::cout<<"number of arguments is: "<<argc<<std::endl;
    if(argc <8)
    {
        std::cout<<"Usage: unmix16 <InputImageFileName1> <InputImageFileName2> "
                   "<OutputImageFileName1> <OutputImageFileName2> "
                   "<MixingMatrix(txt)> "
                   "<GaussianFitlerSize(bgsub)>"
                   "<OmpNumThreads>\n";
        return 0;
    }
    
    int sigma = atoi(argv[6]);
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
    
    size_t numRows=inputImage01->GetLargestPossibleRegion().GetSize()[0];
    size_t numCols=inputImage01->GetLargestPossibleRegion().GetSize()[1];
    size_t numTimes=inputImage01->GetLargestPossibleRegion().GetSize()[2];

    InputPixelType16 * inputImage01Ptr = inputImage01->GetBufferPointer();
    InputPixelType16 * inputImage02Ptr = inputImage02->GetBufferPointer();
    
    // initialize the mixing matrix //
    std::vector< std::vector< float > > mixingMatrix;
    mixingMatrix.resize(2);
    mixingMatrix[0].resize(2);
    mixingMatrix[1].resize(2);
    mixingMatrix[0][0] = 1;
    mixingMatrix[0][1] = 0;
    mixingMatrix[1][0] = 0;
    mixingMatrix[1][1] = 1;

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
    
    // preprocessing code //
    Preporcessor< InputPixelType16, Input2DImageType16 > * pr = new Preporcessor< InputPixelType16, Input2DImageType16 >;
    pr->setInputImage(inputImage01Ptr,inputImage02Ptr);
    pr->setGaussianSigma((InputPixelType16)sigma);
    pr->setMixingMatrix(mixingMatrix);
    pr->setLinearUnmixing();
    pr->setBackgroundSubtraction();
    pr->setDimensions(numRows,numCols,numTimes);
    pr->runPreprocessor();
    delete pr;
    
    
    std::cout<< "I am done ....................."<<std::endl;
    std::string outfName01 = argv[3];
    std::string outfName02 = argv[4];
    std::cout<< "I am writing: "<<outfName01<<std::endl<<std::flush;
    std::cout<< "I am writing: "<<outfName02<<std::endl<<std::flush;
    
    writeImage<InputImageType16>(inputImage01,outfName01.c_str());
    writeImage<InputImageType16>(inputImage02,outfName02.c_str());    
    
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



