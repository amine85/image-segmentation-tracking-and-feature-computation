#include "seg_helpers.h"
#include "omp.h"
#include "itkMultiThreader.h"
//#include "itkSmoothingRecursiveGaussianImageFilter.h"

using namespace seg_helpers;
#if defined(_MSC_VER)
#pragma warning(disable: 4996)
#endif

int optionsCreate(const char* optfile, std::map<std::string,std::string>& options)
{
    options.clear();
    std::ifstream fin(optfile);
    assert(fin.good());
    std::string name; 
    fin>>name;
    while(fin.good()) {
        char cont[100];	
        fin.getline(cont, 99);
        options[name] = std::string(cont);
        fin>>name;
    }
    fin.close();
    return 0;
}


int main(int argc, char **argv)
{

    if(argc <5)
    {
        std::cout<<"Usage: mixture_segment <InputImageFileName>  <LabelImageFileName> <parameters file> <OmpNumThreads>\n";
        return 0;
    }

    int numThreads = atoi(argv[4]);
    // Threads
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);
    omp_set_num_threads(numThreads);

    std::string ifName = argv[1];
    std::string ofName = argv[2];
    std::string opFileName = argv[3];
    
    // To save the seed image
    std::string seedImgfName;
    std::size_t found2;
    found2 = ofName.find_last_of(".");
    seedImgfName = ofName.substr(0,found2)+"seed.tif";

    // declare & define segmentation paramters ////////////////////////////////////
    //default values
    unsigned int K = 3;		// number of intensity mixtures
    unsigned int minVolume = 200;	// minimum object volume
    unsigned int maxVolume = 300;	// maximum object volume
    float gfSig = 2.0; 		// gauss filter parameters
    unsigned int levels = 10;	// number of thresholding levels
    unsigned int winSize = 5; 	// clustering window size
    unsigned int num_time_seeds = 7; // number of time points to estimate seeds
    unsigned int save_seed_image = 0; // flag to save the seed image (seed postfix) (0-no, 1-yes)
	unsigned int min_rad = 15; // minimum radius of a cell
	unsigned int max_rad = 25; // maximum radius of a cell
    
    std::map<std::string, std::string> opts; 
    optionsCreate( opFileName.c_str(), opts);
    std::map<std::string,std::string>::iterator mi;
    mi = opts.find("K"); 
    if(mi!=opts.end())
    {
        std::istringstream ss((*mi).second); 
        ss>>K; 
    }	
    mi = opts.find("minVolume"); 
    if(mi!=opts.end())
    {
        std::istringstream ss((*mi).second); 
        ss>>minVolume; 
    }	
    mi = opts.find("sigmaGfilter"); 
    if(mi!=opts.end())
    {
        std::istringstream ss((*mi).second); 
        ss>>gfSig; 
    }	
    mi = opts.find("numLevel"); 
    if(mi!=opts.end())
    {
        std::istringstream ss((*mi).second); 
        ss>>levels; 
    }	
    mi = opts.find("clusterXY"); 
    if(mi!=opts.end())
    {
        std::istringstream ss((*mi).second); 
        ss>>winSize; 
    }		
    mi = opts.find("num_time_seeds"); 
    if(mi!=opts.end())
    {
        std::istringstream ss((*mi).second); 
        ss>>num_time_seeds; 
    }
    mi = opts.find("save_seed_image"); 
    if(mi!=opts.end())
    {
        std::istringstream ss((*mi).second); 
        ss>>save_seed_image; 
    }

    mi = opts.find("min_rad"); 
    if(mi!=opts.end())
    {
        std::istringstream ss((*mi).second); 
        ss>>min_rad; 
    }

    mi = opts.find("max_rad"); 
    if(mi!=opts.end())
    {
        std::istringstream ss((*mi).second); 
        ss>>max_rad; 
    }
    mi = opts.find("maxVolume"); 
    if(mi!=opts.end())
    {
        std::istringstream ss((*mi).second); 
        ss>>maxVolume; 
    }
	else{
		maxVolume = 0;		// set this to zero if you do not want to use it
	}
	
    std::cout<<"Segmentation parameters:\n";
    std::cout<<"K: "<<K<<std::endl;
    std::cout<<"minVolume: "<< minVolume<<std::endl;
    std::cout<<"sigmaGfilter: "<< gfSig<<std::endl;
    std::cout<<"numLevel: "<< levels<<std::endl;
    std::cout<<"clusterXY: "<< winSize <<std::endl;
    std::cout<<"num_time_seeds: "<< num_time_seeds <<std::endl;
    std::cout<<"save_seed_image: "<< save_seed_image <<std::endl;
    std::cout<<"maxVolume: "<< maxVolume <<std::endl;
    std::cout<<"min_rad: "<< min_rad <<std::endl;
    std::cout<<"max_rad: "<< max_rad <<std::endl;
    
    unsigned found = ifName.find_last_of(".");
    // 	std::cout << " fullName: " << ifName.substr(0,found) << '\n';
    // 	std::cout << " extension: " << ifName.substr(found+1) << '\n';
    std::string seedfName = ifName.substr(0,found) + ".txt";
    std::cout << "Seed FileName: " << seedfName << '\n';
    FILE * fp = fopen(seedfName.c_str(),"w");

    /*******Read Input Image & Allocate Memory for Output Image************/

    // read input image:
    InputImageType::Pointer inputImage = readImage<InputImageType>(ifName.c_str());
    InputPixelType ncc = inputImage->GetLargestPossibleRegion().GetSize()[0];
    InputPixelType nrr = inputImage->GetLargestPossibleRegion().GetSize()[1];
    InputPixelType ns = inputImage->GetLargestPossibleRegion().GetSize()[2];
    InputPixelType * inputImagePtr = inputImage->GetBufferPointer();
    
    // generate output image
    LabelImageType::SizeType size;
    size[0] = ncc;	size[1] = nrr;	size[2] = ns;
    BinImageType::Pointer outputImage = GetITKImageOfSize<BinImageType>(size);
    BinPixelType * outputImagePtr = outputImage->GetBufferPointer();	

    // parameters for 2D image
    InputImageType::SizeType sz2;
    sz2[0] = ncc;	sz2[1] = nrr;	sz2[2] = 1;
    unsigned long long sz = (unsigned long long)ncc*(unsigned long long)nrr;

    // To save the seed image
    LabelImageType::SizeType sizeSeed;
    sizeSeed[0] = ncc;   sizeSeed[1] = nrr;   sizeSeed[2] = num_time_seeds;
    BinImageType::Pointer seedImageGlobal = GetITKImageOfSize<BinImageType>( sizeSeed );
    BinPixelType * seedImageGlobalPtr = seedImageGlobal->GetBufferPointer();    
    
    
    /**************Main Time Loop ****************************************/
#pragma omp  parallel for schedule(dynamic, 1)
    for(size_t it_time = 0; it_time<ns; ++it_time)
    {	
        // 	    std::cout<<"Binarizing Frame: "<<it_time<<std::endl;
        /************** Copy and Filter the Current Frame ****************************************/
        LabelImageType::Pointer iImage = GetITKImageOfSize<LabelImageType>(sz2);
        LabelPixelType * iImagePtr = iImage->GetBufferPointer();
        unsigned long long offset = (unsigned long long)it_time *sz; 

        for(size_t i = 0; i < sz; ++i)
        {
            iImagePtr[i] = inputImagePtr[i+offset];
        }

        if(gfSig>0)
        {
            DiscreteGaussianType::Pointer filter = DiscreteGaussianType::New();
            filter->SetInput(iImage);
            filter->SetVariance(gfSig);
            try
            {
                filter->Update();
            }
            catch(itk::ExceptionObject &err)
            {
                std::cerr << "ExceptionObject caught!" <<std::endl;
                std::cerr << err << std::endl;
            }
            iImage = filter->GetOutput();
            iImage->DisconnectPipeline();
            iImagePtr = iImage->GetBufferPointer();
        }

        /************** Run K means clustering on the intensity image ****************************************/
        BinImageType::Pointer oImage = GetITKImageOfSize<BinImageType>(sz2);
        BinPixelType * oImagePtr = oImage->GetBufferPointer();

		MinimumErrorTresholding(iImagePtr,oImagePtr,sz);
		K = 2;

 //       if(!KmeansClustering(iImagePtr,oImagePtr,sz,K,1))
 //       {
 //           printf("K-means could not converge\n");
 //           if(K>1){
 //               if(!KmeansClustering(iImagePtr,oImagePtr,sz,K-1,1))
 //                   printf("Image might be empty\n");
 //           }
 //       }
        // relabel by largest area: 
        std::map<unsigned int,unsigned int> count_per_label;
        for(unsigned long long i=0 ; i<K ; ++i)
            count_per_label[i] = 0;
        for(unsigned long long i=0 ; i<sz ; ++i)
            count_per_label[oImagePtr[i]]+=1;

        vnl_vector<unsigned int> count_temp(K);
        count_temp.fill(0);
        for(unsigned long long  i=0 ; i<K ; ++i)
            count_temp(i)=count_per_label[i];

        vnl_vector<unsigned int> indices(K);
        indices.fill(0);
        BubbleSortAscend(count_temp,indices);

        // Test for spurios cases, proportion of background
        float background_proportion = (float) count_temp(0) / (float)sz;
        bool isBad = false;
        printf("background_proportion: %f\n",background_proportion);
        if(background_proportion<0.9 ||background_proportion == 1.0 )
            isBad = true;

        std::map<unsigned int,unsigned int> label_correspondence;
        for(unsigned long long  i=0 ; i<K ; ++i)
            label_correspondence[i] = indices[i];
        //         	writeImage<LabelImageType>(oImage, "/data/amine/Data/test/debg-1.tif");

        if(!isBad)
        {
            std::cout << "Is not Bad, time: " << it_time <<std::endl << std::flush;
            /********** copy to binary image & remove the small components	*****************************************/
            for(unsigned long long  i=0 ; i<sz ; ++i)
            {
                if(label_correspondence[oImagePtr[i]]>0)
                    oImagePtr[i] = 1;
            }
            //         	writeImage<LabelImageType>(oImage, "/data/amine/Data/test/debg0.tif");

            RemoveSmallComponents(oImage, minVolume);	

	    if( maxVolume )
            	RemoveLargeComponents(oImage, maxVolume);	
			
            // 		writeImage<LabelImageType>(oImage, "/data/amine/Data/test/debg1.tif");
            std::cout<<"Finished Removing Small Components for Frame: "<<it_time<<std::endl << std::flush;

            /**************** detect the seeds *****************************************************************/
            // only detect seeds for few first time points
            if( it_time < num_time_seeds )
            {
                // seed Image: 
                BinImageType::Pointer seedImage = GetITKImageOfSize<BinImageType>(sz2);
                BinPixelType * seedImagePtr = seedImage->GetBufferPointer();
                // distance Image:
                FloatImageType::Pointer responseImage = GetITKImageOfSize<FloatImageType>(sz2);		
                FloatPixelType * responseImagePtr = responseImage->GetBufferPointer();

                std::cout<<"Computing Distance Map for Frame: "<<it_time<<std::endl << std::flush;
                //GetDistanceResponsev3(iImage, oImage, responseImage, levels); // run this at multiple thresholds
				GetMaxLaplacianOfGaussianResponse( iImage, oImage, responseImage, min_rad, max_rad, sz); 
                std::cout<<"Clustering: "<<it_time<<std::endl << std::flush;
                GetSeedImagev3(responseImagePtr, seedImagePtr,(unsigned int)ncc,(unsigned int)nrr,winSize);
                
                // DEBUG, to save seeds and response
//                     std::ostringstream oss3;
//                     oss3 << "/home/melquiades/aLab/00-nano-team/bin/aaaa" << it_time << ".nrrd";
//                     std::string name3 = oss3.str();
//                     writeImage< FloatImageType >( responseImage, name3.c_str() );
//                     
//                     std::ostringstream oss4;
//                     oss4 << "/home/melquiades/aLab/00-nano-team/bin/bbbb" << it_time << ".nrrd";
//                     std::string name4 = oss4.str();
//                     writeImage< BinImageType >( seedImage, name4.c_str() );

                // write the seed image file //
                std::cout<<"Finished Detecting Seeds for Frame: "<<it_time<<std::endl << std::flush;
#pragma omp critical 
                {
                    for(unsigned long long  i=0 ; i<sz ; ++i)
                    {
                        if(seedImagePtr[i] > 0)
                        {
                            unsigned int xcol = i % ncc;
                            unsigned int yrow = (unsigned int)((double)(i - xcol)/(double)(ncc)) ;
                            fprintf(fp,"%d\t%d\t%d\n",(unsigned int)it_time,xcol,yrow);
                        }
                    }
                    // DEBUG Save the seed image, and erode it by 1 pixel each side
                    if( save_seed_image == 1 ){
                        for(unsigned long long  i=0 ; i<sz ; ++i){
                            int xcol = i % ncc;
                            int yrow = (unsigned int)((double)(i - xcol)/(double)(ncc)) ;
                            if( seedImagePtr[i]>0 ){
                                for( int yy = std::max(0,yrow-1); yy <= std::min(yrow+1,(int)nrr); ++yy ){
                                    for( int xx = std::max(0,xcol-1); xx <= std::min(xcol+1,(int)ncc); ++xx ){
                                        int indx = xx+yy*nrr;
                                        seedImageGlobalPtr[ it_time*sz+indx ] = seedImagePtr[i]>0?255:0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            for(unsigned long long  i=0 ; i<sz ; ++i)
                outputImagePtr[i+offset] = oImagePtr[i]>0?255:0;
        }
        else
        {
            std::cout << "Is Bad, time: " << it_time << "\n";
            for(unsigned long long  i=0 ; i<sz ; ++i)
                outputImagePtr[i+offset] = 0;		// set label image to zeros if too much background has been segmented
        }// end of is bad segmentation if statement
    }// end of time loop
    
    fclose(fp);
    
    // Write output image
    writeImage<BinImageType>(outputImage, ofName.c_str());
    
    // Write the seed image
    if( save_seed_image == 1 )
        writeImage<BinImageType>(seedImageGlobal, seedImgfName.c_str());
    
    
    return 0;
}
