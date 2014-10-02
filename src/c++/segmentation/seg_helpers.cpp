#include "seg_helpers.h"

#if defined(_MSC_VER)
#pragma warning(disable: 4018)
#pragma warning(disable: 4996)
#pragma warning(disable: 4101)
#endif

namespace seg_helpers{
    //***********************************************************************************//
    bool KmeansClustering(const InputPixelType * Image, BinPixelType * BinImage, unsigned int  N, int K, int iternum)
    {
        unsigned int i,j,r;
        unsigned int max_reinitializtion = 100;
        unsigned int num_reinitializtion = 0;
        unsigned char * Label = new unsigned char[N]; // old labels
        InputPixelType max_pixel_val =  GetMax(Image,N);

        float GlobalCost = 0.0;
        std::vector<unsigned char> labels;
        labels.resize(N);

        std::cout<<"Max Pixel Value: "<< max_pixel_val<< std::endl;
        // Run iternum Number of iterations of k-means
        for( i = 0; i<iternum; ++i)
        {
            float Cost = 0.0;

            vnl_matrix<float> OldCentroid(K,1);
            OldCentroid.fill(0);
            vnl_matrix<float> NewCentroid(K,1);

            num_reinitializtion = 0;		// reset the number of reinitializations to zero at each iteration
            printf("centroid: ");

            // Initialize to random values
            for( j = 0; j<K; ++j)
            {
                OldCentroid.put(j,0,(float) (rand()% (max_pixel_val/K) + (j*(max_pixel_val/K))));
                std::cout<<OldCentroid.get(j,0)<<", ";
            }
            std::cout << std::endl;

            for( r=0; r<N; ++r)
                Label[r] = 0;

            // Loop until convergence
            bool converged = false;
            while(!converged)
            {
                for( r=0; r<N; ++r)
                    Label[r] = 0;
                // assignment step of kmeans 
                vnl_vector<double> distance(K);
                for( r=0; r<N; ++r)
                {			
                    distance.fill(0);
                    for(j=0 ; j<K; ++j) // calculate eucledian distance
                    {
                        float diff = (float)Image[r] -  OldCentroid.get(j,0);
                        distance[j] = (double) diff*diff;	
                    }
                    Label[r] = (unsigned int)distance.arg_min();	// store the labels
                    //clabel(r) = distance.arg_min();	// store the labels
                }

                // update step of kmeans  

                // Count the number of elements per each class
                vnl_vector<unsigned int> count(K);
                count.fill(0);

                for( r=0; r<N; ++r)
                    count(Label[r]) += 1;

                // 			for( r=0; r<N; ++r)
                // 			{
                // 				count(clabel(r))+=1;
                // 			}			
                /*			printf("count:\n");

                            for( r=0; r<K; ++r)
                            std::cout<<count[r]<<std::endl;*/

                if(count.min_value()==0)
                {
                    printf("found an empty cluster, reinitializing..\n");
                    //				printf("centroid:\n");
                    for( j = 0; j<K; ++j)
                    {
                        // 					unsigned int idx =  (unsigned int) rand() % N ;
                        // 					OldCentroid[j] = (float) Image[idx];
                        //OldCentroid[j] = (float) (rand()% (max_pixel_val/K) + (j*(max_pixel_val/K)));
                        OldCentroid.put(j,0,(float) (rand()% (max_pixel_val/K) + (j*(max_pixel_val/K))));

                        //std::cout<<OldCentroid[j]<<std::endl;
                    }
                    num_reinitializtion += 1;
                    if(num_reinitializtion>max_reinitializtion)
                        return false;
                    continue;
                }

                NewCentroid.fill(0);
                // update the means
                for( r=0; r<N; ++r)
                {
                    NewCentroid.put(Label[r],0,NewCentroid.get(Label[r],0)+ ((float)Image[r]/(float)count[Label[r]]));
                }

                // compute change and check for convergence
                //float change = (NewCentroid-OldCentroid).magnitude();
                float change = (NewCentroid-OldCentroid).frobenius_norm();

                if(change < 0.001)
                {
                    printf("converged with change: %f\n",change);
                    converged = true;
                }
                else
                {
                    OldCentroid = NewCentroid;
                }
            }// end of while loop

            if(i==0)
            {
                for( r=0; r<N; ++r)
                    BinImage[r] = Label[r];
                // calculate current cost
                Cost = 0.0;
                for( r=0; r<N; ++r)
                    Cost += (float) std::pow(2,(float)Image[r]- NewCentroid.get(Label[r],0));
                //			Cost += (float) std::pow(2,(float)Image[r]- NewCentroid.get(clabel(r),0));
                GlobalCost = Cost;
            }
            else
            {
                // calculate current cost
                Cost = 0.0;
                for( r=0; r<N; ++r)
                    Cost += (float) std::pow(2,(float)Image[r]- NewCentroid.get(Label[r],0));
                //			Cost += (float) std::pow(2,(float)Image[r]- NewCentroid.get(clabel(r),0));

                if(Cost == MIN(Cost,GlobalCost))
                {
                    for( r=0; r<N; ++r)
                        BinImage[r] = Label[r];
                    GlobalCost = Cost;
                }
            }
        }// end of iteration loop
        delete [] Label;
        return true;
    }


    //***********************************************************************************//
	void MinimumErrorTresholding( const InputPixelType * Image, BinPixelType * LabelImage, unsigned int N )
	{
		std::cout<< "In min error thresholding..."<<std::endl;
		/*
		 *	Compute max value to rescale the image data first
		 */

		 InputPixelType max = GetMax( Image, N );
		 std::cout<<"max: "<<max<<std::endl;
		 
		 /*
		  *	Compute probabilites
		  */
		 double* prob = new double[max+1];
		 for(unsigned int i = 0; i<=max; ++i)
		 	prob[i] = 0;

		 // histogram
		 for(unsigned int i = 0; i<N; ++i)
		 	prob[ Image[i] ] += 1;

		 // normalize
		 for(unsigned int i = 0; i<=max; ++i)
		 	prob[i]/= (double)N;
	
		/*
		 * Minimize the cost function
		 */
		

		double mu = 0.0;
		for(unsigned int i = 0; i<=max; ++i)
		 	mu += prob[i] * (double)i;

		std::cout<<"global mean: "<<mu<<std::endl;
		
		double mu_0 = 0.0;
		double p_0 = 0.0;
		double mu_1 = mu;
		double p_1 = 1.0;

		unsigned int th = max;
		double cost = std::numeric_limits<double>::max();
		
		for(unsigned int i = 0; i<=max; ++i){
			
			p_0 += prob[i];
			p_1 -= prob[i];

			mu_0 += prob[i]*(double)i;
			mu_1 -= prob[i]*(double)i;

			double m0 = mu_0/(p_0+EPS);
			double m1 = mu_1/(p_1+EPS);
			double c = mu - (p_0*( std::log(p_0) + m0 * std::log(m0))) - (p_1*(std::log(p_1) + m1* std::log(m1)));
			if( cost != MIN(c,cost) ){
				cost = c;
				th = i;
			}
		}
		std::cout<<" minimum error threshold: "<<th<<std::endl;
		/*
		 *	Treshold the image data
		 */
		
		 for(unsigned int i = 0; i<N; ++i)
		 	LabelImage[i] = ( Image[i]>=th ? 1:	0 );

		/*
		 *	Clean up
		 */
		delete[] prob;

	}
	
    //***********************************************************************************//
	void GetMaxLaplacianOfGaussianResponse( InputImageType::Pointer image, BinImageType::Pointer binaryImage, \
			FloatImageType::Pointer responseImage, unsigned int min_radius, unsigned int max_radius, unsigned int N)
	{	
		std::cout<<"starting log computation ...\n";
		
		/*
		 *	Copy the image to a 2d image because of itk log filter
		 */
		InputImageType2D::SizeType sz;
		sz[0] = image->GetLargestPossibleRegion().GetSize()[0];
    	sz[1] = image->GetLargestPossibleRegion().GetSize()[1];
		InputImageType2D::Pointer image2d = GetITKImageOfSize<InputImageType2D>(sz);
		InputPixelType* image2dPtr = image2d->GetBufferPointer();
		InputPixelType* imagePtr = image->GetBufferPointer();
		for( unsigned int i=0; i<N; ++i )
			image2dPtr[i] = imagePtr[i];		
		
		typedef itk::LaplacianRecursiveGaussianImageFilter< InputImageType2D, FloatImageType2D > logType;

		FloatPixelType* response = responseImage->GetBufferPointer();
		BinPixelType* binary = binaryImage->GetBufferPointer();

		float sqrt2 = std::sqrt(2);

		// let's store index of foreground pixels ( gain some speed )
		std::vector<unsigned int> index;
		for( unsigned int i=0; i<N; ++i){
			if( binary[i] )
				index.push_back(i);
		}
		
		for( unsigned int r = min_radius; r<=max_radius; ++r ){
			
			float sigma = (float)r/sqrt2;
			float sigma2 = sigma*sigma;
			
			/*
			 *	Apply log filter
			 */
			logType::Pointer logFilter = logType::New();
			logFilter->SetInput( image2d );
			logFilter->SetSigma( sigma );
			try{
				logFilter->Update();
			}
			catch( itk::ExceptionObject & err ) {
				std::cerr << "Error calculating log response: " << err << std::endl ;
			}	
			/*
			 * Get max response	
			 */
			FloatPixelType* logResponse = logFilter->GetOutput()->GetBufferPointer();
			
			for( unsigned int i=0; i<index.size(); ++i )
				response[index[i]] =  MAX( response[index[i]], -sigma2*logResponse[index[i]] );

		}

		// try to multiply by distmap
		typedef itk::SignedMaurerDistanceMapImageFilter<BinImageType,FloatImageType>  DistMapFType;
		// compute distance map:
		DistMapFType::Pointer distFilter = DistMapFType::New();
		distFilter->SetInput(binaryImage) ;
		distFilter->SetSquaredDistance( false );      
		distFilter->SetInsideIsPositive( true );
		try {
			distFilter->Update() ;
		}
		catch( itk::ExceptionObject & err ) {
			std::cerr << "Error calculating distance transform: " << err << std::endl ;
		}	

		FloatPixelType* distPtr = distFilter->GetOutput()->GetBufferPointer();
		for( unsigned int i=0; i<index.size(); ++i)
			response[index[i]]*= distPtr[index[i]];

		std::cout<<"finished log computation....\n";

	}
    //***********************************************************************************//
    InputPixelType GetMax(const InputPixelType * X, unsigned int M)
    {
        unsigned int i;
        InputPixelType max = X[0];
        for(i=1 ; i<M ; ++i)
        {
            max = MAX(X[i],max);
        }
        return max;
    }
    //***********************************************************************************//
    void BubbleSortAscend(vnl_vector<unsigned int> &X, vnl_vector<unsigned int> &indices)
    {
        bool swaped = false;
        unsigned int N = X.size();
        unsigned int i,j,pass;
        unsigned int temp,tempidx;

        for(i=0 ; i<N; ++i)
            indices[i] = i;

        for(pass=0; pass<N; ++pass)
        {
            for(i=0 ; i<N-1; ++i)
            {
                j = i+1;
                if(X[j]>X[i])
                {
                    temp = X[j];
                    tempidx = indices[j];

                    X[j] = X[i];
                    indices[j] = indices[i];				

                    X[i] = temp;
                    indices[i] = tempidx;
                    swaped = false;
                }
            }
        }
    }
    //***********************************************************************************//
    void RemoveLargeComponents(BinImageType::Pointer BinImage, unsigned int max_volume)
    {
        //printf("Removing small connected components ...\n");

        typedef itk::RelabelComponentImageFilter<LabelImageType,LabelImageType> RelabelFilterType;
        typedef itk::ScalarConnectedComponentImageFilter<BinImageType,LabelImageType> ConnCompFilterType;

        ConnCompFilterType::Pointer ccfilter = ConnCompFilterType::New();
        ccfilter->SetInput( BinImage );
        ccfilter->SetFullyConnected(1);
        ccfilter->SetDistanceThreshold(0);
        try {
            ccfilter->Update();
        }
        catch( itk::ExceptionObject & err ) {
            std::cerr << "Error ccmp filter: " << err << std::endl ;
            // return -1;
        }

        // make sure the background is zero //
        LabelImageIteratorType it(ccfilter->GetOutput(),ccfilter->GetOutput()->GetLargestPossibleRegion());
        BinImageIteratorType it1(BinImage,BinImage->GetLargestPossibleRegion());
        for(it.GoToBegin(),it1.GoToBegin();!it.IsAtEnd(); ++it,++it1)
        {
            // 		if(it.Get()==1)
            // 			it.Set(0);
            if(it1.Get()==0)
                it.Set(0);
        }

        RelabelFilterType::Pointer rfilter = RelabelFilterType::New();
        rfilter->SetInput(ccfilter->GetOutput());
        rfilter->InPlaceOn();
        try {
            rfilter->Update();
        }
        catch( itk::ExceptionObject & err ) {
            std::cerr << "Error relabel filter: " << err << std::endl ;
            // return -1;
        }
        
        RelabelFilterType::ObjectSizeInPixelsContainerType Volumes = rfilter->GetSizeOfObjectsInPixels();
        LabelPixelType labelValue = Volumes.size()+1;// just to be safe
        
        //printf("minimum volume: %d\n",min_volume);
        //printf("Volume.size: %d\n",Volumes.size());
        for(LabelPixelType i=0; i<Volumes.size(); ++i)
        {
            //printf("volumes[%d]:%d\n",i,Volumes[i]);
            if(Volumes[i] <  max_volume)
            {
                labelValue = i;
                //printf("found labelValue\n");
                break;
            }
        }

        LabelImageIteratorType filterIter(rfilter->GetOutput(),rfilter->GetOutput()->GetLargestPossibleRegion());
        BinImageIteratorType imageIter(BinImage,BinImage->GetLargestPossibleRegion());

        for(filterIter.GoToBegin(),imageIter.GoToBegin();!filterIter.IsAtEnd(); ++imageIter,++filterIter)
        {
            if(filterIter.Get() < labelValue)
                imageIter.Set(0);
        }
    }


    //***********************************************************************************//
    void RemoveSmallComponents(BinImageType::Pointer BinImage, unsigned int min_volume)
    {
        //printf("Removing small connected components ...\n");

        typedef itk::RelabelComponentImageFilter<LabelImageType,LabelImageType> RelabelFilterType;
        typedef itk::ScalarConnectedComponentImageFilter<BinImageType,LabelImageType> ConnCompFilterType;

        ConnCompFilterType::Pointer ccfilter = ConnCompFilterType::New();
        ccfilter->SetInput( BinImage );
        ccfilter->SetFullyConnected(1);
        ccfilter->SetDistanceThreshold(0);
        try {
            ccfilter->Update();
        }
        catch( itk::ExceptionObject & err ) {
            std::cerr << "Error ccmp filter: " << err << std::endl ;
            // return -1;
        }

        // make sure the background is zero //
        LabelImageIteratorType it(ccfilter->GetOutput(),ccfilter->GetOutput()->GetLargestPossibleRegion());
        BinImageIteratorType it1(BinImage,BinImage->GetLargestPossibleRegion());
        for(it.GoToBegin(),it1.GoToBegin();!it.IsAtEnd(); ++it,++it1)
        {
            if(it1.Get()==0)
                it.Set(0);
        }

        RelabelFilterType::Pointer rfilter = RelabelFilterType::New();
        rfilter->SetInput(ccfilter->GetOutput());
        rfilter->InPlaceOn();
        try {
            rfilter->Update();
        }
        catch( itk::ExceptionObject & err ) {
            std::cerr << "Error relabel filter: " << err << std::endl ;
            // return -1;
        }
        
        RelabelFilterType::ObjectSizeInPixelsContainerType Volumes = rfilter->GetSizeOfObjectsInPixels();
        LabelPixelType labelValue = Volumes.size()+1;// just to be safe
        
        for(LabelPixelType i=0; i<Volumes.size(); ++i)
        {
            if(Volumes[i] < min_volume)
            {
                labelValue = i;
                break;
            }
        }

        LabelImageIteratorType filterIter(rfilter->GetOutput(),rfilter->GetOutput()->GetLargestPossibleRegion());
        BinImageIteratorType imageIter(BinImage,BinImage->GetLargestPossibleRegion());

        for(filterIter.GoToBegin(),imageIter.GoToBegin();!filterIter.IsAtEnd(); ++imageIter,++filterIter)
        {
            if(filterIter.Get()>labelValue)
                imageIter.Set(0);
        }
    }
    //***********************************************************************************//
    void GetSeedImagev3(FloatPixelType * ResponseImagePtr, BinPixelType * SeedImagePtr,\
            const unsigned int nr, const unsigned int nc, const unsigned int window_size)
    {

        // iterate through the image
        printf("nr:%d\n",nr);
        printf("nc:%d\n",nc);
        printf("window_size:%d\n",window_size);

        int x,y,idx,i,j,k,xx,yy,dumy;
        FloatPixelType val,nval;
        int hws = (int)window_size;
        int minX = hws +1 ;
        int maxX = (int)nr - hws- 1;
        int minY = hws + 1 ;
        int maxY = (int)nc - hws -1;  
        int ws = 2*hws+1;

        printf("ws:%d\n",ws);

        for( y = minY; y < maxY; ++y)
        {    
            dumy = y*nr;
            for( x = minX; x < maxX; ++x)
            {
                idx = dumy + x;
                //printf("before indexing center\n");
                //printf("(x,y) : (%d,%d)\n",x,y);
                val = ResponseImagePtr[idx];           // pixel value       
                SeedImagePtr[idx] = 0;      // set to local maximum 


                if(val>0.0)
                {
                    //printf("val: %f\n",val);

                    // look in the neighborhood of the window 
                    int counter = 0;
                    for( i = -hws; i<=hws; ++i)
                    {
                        xx = x + i;
                        for( j = -hws; j<=hws; ++j)
                        {
                            // neighbor pixel locations
                            if( i==0 && j==0)                   // center pixel
                                continue;

                            yy = y + j;
                            k = yy*nr + xx;
                            nval = ResponseImagePtr[k];        // neighbor value
                            //printf("val: %f\n",val);
                            //printf("nval: %f\n",nval);
                            if(val>nval)
                            {
                                counter +=1;
                            }

                        }//end of window col loop  
                    }// end of window row loop   
                    //printf("counter:%d\n",counter);

                    if(counter==(ws*ws-1))
                    {
                        SeedImagePtr[idx] = std::numeric_limits<LabelPixelType>::max();
                        //printf("added maximum\n");
                    }
                    // if it is still a local maximum, supress the neighbors
                    if(SeedImagePtr[idx]>0)
                    {
                        for( i = -hws; i<=hws; ++i)
                        {
                            xx = x + i;
                            for( j = -hws; j<=hws; ++j)
                            {
                                // neighbor pixel locations
                                if( i==0 && j==0)               // center pixel
                                    continue;                                   
                                yy = y + j;
                                k = yy*nr + xx;
                                SeedImagePtr[k] = 0;
                                //printf("suppressing neigbors\n");
                            }//end of window col loop  
                        }// end of window row loop                
                    }// end of suppression if statement
                }// end of if(val>0) 

            }// end of row image loop
        }// end of col image loop

    }
    //**********************************************************************************************************//
    void GetDistanceResponsev3(LabelImageType::Pointer InputImage,BinImageType::Pointer BinaryImage, \
            FloatImageType::Pointer ResponseImage,unsigned int levels)
    {
        
        LabelPixelType * InputImagePtr = InputImage->GetBufferPointer();
        BinPixelType * BinaryImagePtr = BinaryImage->GetBufferPointer();
        // 1- threshold the images
        // 2- compute distance map 
        // 3- normalize 

        // set up filter types
        typedef itk::BinaryThresholdImageFilter<LabelImageType,LabelImageType>  ThresholdFType;
//         typedef itk::MedianImageFilter<LabelImageType,LabelImageType>  MedianFType;
        typedef itk::ScalarConnectedComponentImageFilter<BinImageType,LabelImageType> ConCompFType;
        typedef itk::SignedMaurerDistanceMapImageFilter<LabelImageType,FloatImageType>  DistMapFType;
        typedef itk::LabelGeometryImageFilter< LabelImageType > LabGeomFType;  
        typedef itk::DiscreteGaussianImageFilter<FloatImageType,FloatImageType> GaussianFType;

        // get image parameters //
        size_t nr = InputImage->GetLargestPossibleRegion().GetSize()[0];
        size_t nc = InputImage->GetLargestPossibleRegion().GetSize()[1];    
        size_t nt = InputImage->GetLargestPossibleRegion().GetSize()[2]; 
        size_t offset = nc*nt;

        // get connected components //
        ConCompFType::Pointer connFilter = ConCompFType::New();
        connFilter->SetInput(BinaryImage);
        connFilter->SetFullyConnected(1);
        connFilter->SetDistanceThreshold(0);
        try{
            connFilter->Update();
        }
        catch( itk::ExceptionObject & err ) {
            std::cerr << "Error calculating connected components: " << err << std::endl ;
            // return -1;
        }     

        // make sure the background is zero 
        LabelImageIteratorType it(connFilter->GetOutput(),connFilter->GetOutput()->GetLargestPossibleRegion());
        BinImageIteratorType it1(BinaryImage,BinaryImage->GetLargestPossibleRegion());
        for(it.GoToBegin(),it1.GoToBegin();!it.IsAtEnd(); ++it,++it1)
        {
            if(it.Get()==1)
                it.Set(0);
            if(it1.Get()==0)
                it.Set(0);
        }
        // get labels 
        LabGeomFType::Pointer labGeomFilter = LabGeomFType::New();
        labGeomFilter->SetInput(connFilter->GetOutput());
        labGeomFilter->CalculatePixelIndicesOn();
        try{
            labGeomFilter->Update();
        }
        catch( itk::ExceptionObject & err ) {
            std::cerr << "Error calculating pixel indices: " << err << std::endl ;
        }  

        float min_intensity = std::numeric_limits<float>::max();
        float max_intensity = -std::numeric_limits<float>::max();
        // first mixture 
        unsigned long long sz = nr*nc;
        for(unsigned long long  i = 0; i < sz; ++i)
        {
            if(BinaryImagePtr[i]>0)
            {
                min_intensity = MIN(InputImagePtr[i],min_intensity);
                max_intensity = MAX(InputImagePtr[i],max_intensity);
            }
        }
        
        // main threshold loop
        float step = (max_intensity-min_intensity)/(float)levels;
        FloatPixelType * reponsePtr = ResponseImage->GetBufferPointer();

        for(unsigned int i=0; i<levels; ++i)
            // for(unsigned int i=0; i<1; ++i)
        {
            float lower_threshold = min_intensity + i*step;

            // threshold 
            ThresholdFType::Pointer thresholdFilter = ThresholdFType::New();
            thresholdFilter->SetInput(InputImage);
            thresholdFilter->SetLowerThreshold(lower_threshold);
            thresholdFilter->SetUpperThreshold(max_intensity);
            thresholdFilter->SetInsideValue(1);
            thresholdFilter->SetOutsideValue(0);    
            thresholdFilter->Update();
            //     writeImage<LabelImageType>(thresholdFilter->GetOutput(), "/data/amine/Data/test/thresholded1.tif");

            // compute distance map:
            DistMapFType::Pointer distFilter = DistMapFType::New();
            distFilter->SetInput(thresholdFilter->GetOutput()) ;
            distFilter->SetSquaredDistance( false );      
            distFilter->SetInsideIsPositive( true );
            try {
                distFilter->Update() ;
            }
            catch( itk::ExceptionObject & err ) {
                std::cerr << "Error calculating distance transform: " << err << std::endl ;
            }
            // set the background pixels to zero
            FloatImageIteratorType distIterator(distFilter->GetOutput(),distFilter->GetOutput()->GetLargestPossibleRegion());
            distIterator.GoToBegin();

            while(!distIterator.IsAtEnd())
            {
                if(distIterator.Get()<=0)
                    distIterator.Set(0.0);
                ++distIterator;
            }   
            //writeImage<FloatImageType>(distFilter->GetOutput(), "/data/amine/Data/test/dist_image.nrrd");

            GaussianFType::Pointer gaussianFilter = GaussianFType::New();
            gaussianFilter->SetInput(distFilter->GetOutput());
            gaussianFilter->SetVariance(1.0);
            try
            {
                gaussianFilter->Update();
            }
            catch(itk::ExceptionObject &err)
            {
                std::cerr << "ExceptionObject caught!" <<std::endl;
                std::cerr << err << std::endl;
            }
            FloatImageType::Pointer DistImage = gaussianFilter->GetOutput();   

            //writeImage<FloatImageType>(DistImage, "/data/amine/Data/test/dist_image_gf.nrrd");

            FloatPixelType * distPtr = DistImage->GetBufferPointer();   
            LabGeomFType::LabelsType allLabels = labGeomFilter->GetLabels();
            LabGeomFType::LabelsType::iterator labelIterator;   

            for( labelIterator = allLabels.begin(); labelIterator != allLabels.end(); labelIterator++ )
            {   
                LabGeomFType::LabelPixelType labelValue = *labelIterator;
                //std::cout << "at label: " << labelValue << std::endl;
                if(labelValue==0)       // do not want to include the background
                    continue;

                std::vector<LabGeomFType::LabelIndexType> indices = labGeomFilter->GetPixelIndices( labelValue );
                //std::cout << "number of indices: " << indices.size() << std::endl;

                float max_dist = -std::numeric_limits<float>::max();
                for(unsigned int i = 0; i<indices.size(); ++i)
                {
                    size_t index =indices[i][0]+indices[i][1]*nr + indices[i][2]*offset;
                    max_dist  = MAX(max_dist,distPtr[index]);
                }// end of label index loop

                for(unsigned int i = 0; i<indices.size(); ++i)
                {
                    size_t index =indices[i][0]+indices[i][1]*nr + indices[i][2]*offset;
                    reponsePtr[index] += (distPtr[index]/(max_dist+EPS));
                }// end of label index loop           
            }// end of label loop
        }// end of threshold level loop
    }

}// end of namespace
