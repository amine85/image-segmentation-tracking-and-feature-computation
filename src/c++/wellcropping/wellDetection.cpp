
#include"helper.h"

ImageF2::Pointer SubsBack( ImageF2::Pointer imgIn, double sigma );
void FindInitialWells( ImageF2::Pointer INCCPat1, ImageF2::Pointer INCCPat2, int wellSize, int minBorder, std::vector< std::vector< double > > *locations, ImageCHAR2::Pointer* WellCenterImage );
void FindShift( ImageCF2::Pointer*sliceImageFFT, ImageCF2::Pointer*sliceImagePreviousFFT, std::vector< int >*shift, bool flagRotate );
std::vector< std::vector< std::vector< int > > > LabelWells( std::vector< std::vector< double > >*locations, int wellPerRow, ImageF2::SizeType size );

template <typename T>
void CropImages( typename T::Pointer imageStackCrop, std::vector< std::vector< std::vector< int > > > allWells, std::vector< std::vector< int > > shift, std::string outputPath, int wellSize, std::string channelName, int wellPerRow );

int main(int argc, char *argv[])
{
    if( argc<6 ){
        std::cout<<std::endl<<"At least one channel required"<<std::endl<<std::flush;
        std::cout<<"Usage: wellDetection <pattimgname_1>  <pattimgname_2> <outputPath> <initial_slice> <wellPerRow > <OmpNumThreads> <image...>\n";
        exit(1);
    }

    // Read the parameters
    const char * pattimgname_1 = argv[1];
    const char * pattimgname_2 = argv[2];
    std::string outputPath = argv[3];
    int initial_slice = atoi(argv[4]);
    int wellPerRow = atoi(argv[5]);

    int numThreads = atoi(argv[6]);
    // Threads
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(numThreads);
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(numThreads);
    omp_set_num_threads(1);
    
    // Read the images to be cropped, the first one is always the brightfield
    std::vector< const char* > imageNames;
    std::vector< std::string > imageCH;
    
    for( int i=7; i<argc; ++i ){
        imageNames.push_back(argv[i]);          
        std::string temp( argv[i] );
        std::string temp2( argv[i] );

        std::size_t found, found2;
        found = temp.find_last_of(".");
        temp = temp.substr(found-3,3);
        found = temp2.find("/bin_");
        found2 = temp2.find("/bg_");
        if ( (found==std::string::npos) && (found2==std::string::npos) )
            imageCH.push_back( temp );
        else if( found!=std::string::npos )
            imageCH.push_back( temp+"bin" );
        else if( found2!=std::string::npos )
            imageCH.push_back( temp+"bg" );
        else
            imageCH.push_back( temp );
        std::cout << " file: " << imageCH.at(i-7) <<std::endl<<std::flush;
    }
    // Read images
    std::vector<ImageCHAR3::Pointer> imageStackCropBin;
    imageStackCropBin.resize( imageNames.size() );
    std::vector<ImageUS3::Pointer> imageStackCrop;
    imageStackCrop.resize( imageNames.size() );

    for( int ii=0; ii<imageNames.size(); ++ii ){
        std::size_t foundBin = imageCH.at(ii).find("bin");
        if (foundBin!=std::string::npos){ // Bin is present in the chanel name
            imageStackCropBin.at(ii) = readImage< ImageCHAR3 >( imageNames[ii] );
        }
        else{
            imageStackCrop.at(ii) = readImage< ImageUS3 >( imageNames[ii] );
        }
    }

    // The brightfield image
    const char* imgname3d = imageNames[0];

    // Area cluster according to the expected well size
    int wellSize;
    if( wellPerRow == 8 ){
        wellSize = 253-20; // (2048/8+31) Each side will be 115, this is for the 8x8 //FIXME
    }
    else if( wellPerRow == 7 ){
        wellSize = 289-20; // (2048/7+31) Each side will be 115, this is for the 7x7 //FIXME
    }
    else if( wellPerRow == 6 ){
	wellSize = 281;
    }
    else{
        std::cout << std::endl << "FIXME: The number of wells per row is not 7 or 8.";
        exit(1);
    }

    int minBorder = 70; // Minimum distance to find a well with respect to the border

    double sigma = 30.0; // For the background substraction

    // -------------------------------------------------------------------------------------------------------------------------
    // FFT stuff
    //   // exercise the name-value conversion methods
    //   itk::FFTWGlobalConfiguration::GetPlanRigorValue("FFTW_EXHAUSTIVE");
    //   itk::FFTWGlobalConfiguration::GetPlanRigorName(FFTW_EXHAUSTIVE);
    //   
    //   itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_EXHAUSTIVE);
    //   itk::FFTWGlobalConfiguration::SetReadWisdomCache(true);
    //   itk::FFTWGlobalConfiguration::SetWriteWisdomCache(true);

    itk::FFTWGlobalConfiguration::GetPlanRigorValue("FFTW_ESTIMATE");
    itk::FFTWGlobalConfiguration::GetPlanRigorName(FFTW_ESTIMATE);

    itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_ESTIMATE);
    itk::FFTWGlobalConfiguration::SetReadWisdomCache(false);
    itk::FFTWGlobalConfiguration::SetWriteWisdomCache(false);

    //   if(argc>1)
    //     {
    //     itk::FFTWGlobalConfiguration::SetWisdomCacheBase(argv[1]);
    //     }
    std::cout << "WriteWisdomCache  " << itk::FFTWGlobalConfiguration::GetWriteWisdomCache() << std::endl;
    std::cout << "ReadWisdomCache  " << itk::FFTWGlobalConfiguration::GetReadWisdomCache() << std::endl;
    std::cout << "PlanRigor  " << itk::FFTWGlobalConfiguration::GetPlanRigor() << std::endl;
    std::cout << "WisdomCacheBase " << itk::FFTWGlobalConfiguration::GetWisdomCacheBase()  << std::endl;
    std::cout << "WisdomeFile     " << itk::FFTWGlobalConfiguration::GetWisdomFileDefaultBaseName() << std::endl;
    // -------------------------------------------------------------------------------------------------------------------------

    // READ IMAGES
    ImageUS3::Pointer imageStack = imageStackCrop.at(0);//readImage< ImageF3 >( imgname3d );
    ImageF2::Pointer pattern_1 = readImage< ImageF2 >( pattimgname_1 );
    ImageF2::Pointer pattern_2 = readImage< ImageF2 >( pattimgname_2 );

    // Get sizes
    ImageUS3::SizeType sizeImageStack = imageStack->GetLargestPossibleRegion().GetSize();
    ImageF2::SizeType sizePattern_1 = pattern_1->GetLargestPossibleRegion().GetSize();
    ImageF2::SizeType sizePattern_2 = pattern_2->GetLargestPossibleRegion().GetSize();

    // Check that sizes matches, and that they are odd.    
    if( ( sizePattern_1[0] != sizePattern_2[0] ) || ( sizePattern_1[1] != sizePattern_2[1] ) ){
        std::cout << std::endl << "The pattern sizes do not match";
        exit(1);
    }
    if( ((sizePattern_1[0] & 1) == 0) || ((sizePattern_1[1]&1)==0) ){
        std::cout << std::endl << "The sizes of the pattern are not odd";
        exit(1);
    }

    // Substract the backgoud, to avoid response due to cells
    ImageF2::Pointer pattern_1Subs = SubsBack( pattern_1, sigma );
    ImageF2::Pointer pattern_2Subs = SubsBack( pattern_2, sigma );

    // Variables to store the output
    std::vector< std::vector< double > > locations;
    std::vector< std::vector< int > > shift;
    shift.resize( sizeImageStack[2] );
    std::vector< std::vector< std::vector< int > > > allWells; // Center of wells, plus the well type
    ImageCF2::Pointer sliceImagePreviousFFT, sliceImageFFT;

    for( long long ii=0; ii<sizeImageStack[2]; ++ii ){
        clock_t t1 = clock();
        std::cout << std::endl << "At time: " << ii << std::flush;

        ImageF2::Pointer sliceImage = ExtractSlide< ImageUS3, ImageF2 >( imageStack, ii );
        if( ii == 0 ){
            ImageF2::Pointer sliceImageBs = SubsBack( sliceImage, sigma );
            ImageF2::Pointer INCCPat1, INCCPat2;
            NCCFun( sliceImageBs, pattern_1Subs, pattern_2Subs, &INCCPat1, &INCCPat2 );

            // Find maximas
            ImageCHAR2::Pointer WellCenterImage;
            FindInitialWells( INCCPat1, INCCPat2, wellSize, minBorder, &locations, &WellCenterImage );

            std::ostringstream oss3;
            oss3 << "/../a_wellDet.nrrd";
            std::string name3 = outputPath+oss3.str();
            writeImage< ImageCHAR2 >( WellCenterImage, name3.c_str() );

            // Constrain the maxima, acording to the expected location of the FindInitialWells

            // Label the different wells according to row or coll. The variable locations has two empty values by default
            ImageF2::SizeType size = sliceImage->GetLargestPossibleRegion().GetSize();
            allWells = LabelWells( &locations, wellPerRow, size );

            // Initialize
            sliceImagePreviousFFT = FFT2< ImageF2, ImageCF2 >( sliceImage );
            continue;
        }

        // Find the shift
        bool flagRotate = ii%2;
        if( flagRotate == 1 )
            sliceImageFFT = FFT2< ImageF2, ImageCF2 >( RotateImage180< ImageF2 > ( sliceImage ) );
        else
            sliceImageFFT = FFT2< ImageF2, ImageCF2 >( sliceImage );

        FindShift( &sliceImageFFT, &sliceImagePreviousFFT, &shift.at(ii), flagRotate );
        // Correct shift according to rotation

        sliceImagePreviousFFT = sliceImageFFT;
    }

    // Crop
    for( int ii=0; ii<imageNames.size(); ++ii ){

        std::size_t foundBin = imageCH.at(ii).find("bin");
        if (foundBin!=std::string::npos){ // Bin is present in the chanel name
            CropImages<ImageCHAR3>( imageStackCropBin.at(ii), allWells, shift, outputPath, wellSize, imageCH.at(ii), wellPerRow );
        }
        else{
            CropImages<ImageUS3>( imageStackCrop.at(ii), allWells, shift, outputPath, wellSize, imageCH.at(ii), wellPerRow );
        }
    }

    // Save information of the crops
    std::string centerFile = outputPath+"/a_centers.txt";
    FILE * centerFilefp = fopen( centerFile.c_str(), "w" );
    for( int row=0; row<wellPerRow; ++row ){
        for( int col=0; col<wellPerRow; ++col ){
            int xx_cen = allWells.at(row).at(col).at(0);
            int yy_cen = allWells.at(row).at(col).at(1);
            fprintf(centerFilefp,"%d\t%d\t%d\t%d\n",row,col,xx_cen,yy_cen);
        }
    }
    fclose(centerFilefp);

    // FIXME SHIFT, is provided in relative value with respect to the previous timepoint 
    // The shift is applied center-shift
    std::string shiftFile = outputPath+"/a_shifts.txt";
    FILE * shiftFilefp = fopen( shiftFile.c_str(), "w" );

    int xx_shift, yy_shift;
    for( int time=0; time<sizeImageStack[2]; ++time ){
        if( time == 0 ){
            xx_shift = 0;
            yy_shift = 0;
        }
        else{
            xx_shift = shift.at(time).at(0);
            yy_shift = shift.at(time).at(1);
        }
        fprintf( shiftFilefp,"%d\t%d\t%d\n",time,xx_shift,yy_shift);
    }
    fclose( shiftFilefp );

    return 0;

}

template <typename T>
void CropImages( typename T::Pointer imageStackCrop, std::vector< std::vector< std::vector< int > > > allWells, std::vector< std::vector< int > > shift, std::string outputPath, int wellSize, std::string channelName, int wellPerRow ){

    typename T::SizeType size = imageStackCrop->GetLargestPossibleRegion().GetSize();
    typename T::PixelType * imageStackCropPtr = imageStackCrop->GetBufferPointer();

    // Crop image
    typename T::SizeType sizeWell;
    sizeWell[0] = wellSize;
    sizeWell[1] = wellSize;
    sizeWell[2] = size[2];
    typename T::Pointer WellImage = GetITKImageOfSize<T>( sizeWell );
    typename T::PixelType * WellImagePtr = WellImage->GetBufferPointer();

    for( int row=0; row<wellPerRow; ++row ){
        for( int col=0; col<wellPerRow; ++col ){
            // Clear output
            for( int k=0; k<size[2]; ++k ){
                for( int j=0; j<wellSize; ++j ){
                    for( int i=0; i<wellSize; ++i ){
                        WellImagePtr[ i + j*wellSize + k*wellSize*wellSize ] = 0;
                    }
                }
            }
            // Not a valid well. FIXME, this can be done easily, just look for the max on that are, but I think is better to keep it like this
            if( allWells.at(row).at(col).at(0) == -1000 || allWells.at(row).at(col).at(1) == -1000 ){
                // In case the well is not detected, it is initialized with -1000
                std::size_t found = channelName.find("CH0");
                if (found!=std::string::npos)
                    std::cout << std::endl << "FIXME: well missing at R:" << row+1 << ", C:" << col+1 << ". IN: " << outputPath << std::flush;
                continue;
            }
            bool flagToSave = 0;
            int xx_cen = allWells.at(row).at(col).at(0)-(wellSize-1)/2;
            int yy_cen = allWells.at(row).at(col).at(1)-(wellSize-1)/2;
            // Copy data
            for( int k=0; k<size[2]; ++k ){
                for( int j=0; j<wellSize; ++j ){
                    if( yy_cen+j < 0 || yy_cen+j >= size[1] )
                        continue;
                    for( int i=0; i<wellSize; ++i ){
                        if( xx_cen+i < 0 || xx_cen+i >= size[0] )
                            continue;
                        typename T::PixelType value = imageStackCropPtr[ xx_cen+i + (yy_cen+j)*size[0] + k*size[0]*size[1] ];
                        WellImagePtr[ i + j*wellSize + k*wellSize*wellSize ] = value;
                        if( value > 0 ){
                            flagToSave = 1;
                        }
                    }
                }
                // Apply the shift
                if( k < size[2]-1 ){
                    //                     std::cout << std::endl << "SHIFT: " << shift.at(k+1).at(0) << " " << shift.at(k+1).at(1);
                    xx_cen = xx_cen - shift.at(k+1).at(0); ///////////////////////
                    yy_cen = yy_cen - shift.at(k+1).at(1); ///////////////////////
                }
            }

            // Only save the wells that are not empty
            if( flagToSave == 0 )
                continue;
            // Save well
            std::ostringstream oss1;
            oss1 << "/imgR" << row+1 << "C" << col+1 << channelName << ".tif";
            std::string name1 = outputPath+oss1.str();
            writeImage< T >( WellImage, name1.c_str() );
        }
    }
}



std::vector< std::vector< std::vector< int > > > LabelWells( std::vector< std::vector< double > >*locations, int wellPerRow, ImageF2::SizeType size ){
    std::vector< std::vector< std::vector< int > > > allWells;
    allWells.resize( wellPerRow );
    for( int row=0; row<wellPerRow; ++row ){
        allWells.at(row).resize( wellPerRow );
        for( int col=0; col<wellPerRow; ++col ){
            allWells.at(row).at(col).resize( 4 );
            allWells.at(row).at(col).at(0) = -1000; // Not valid by default
            allWells.at(row).at(col).at(1) = -1000; // Not valid by default
            allWells.at(row).at(col).at(2) = -1000; // Not valid by default
            allWells.at(row).at(col).at(3) = -1000; // Not valid by default
        }
    }
    //     if( (*locations).size() > wellPerRow*wellPerRow )
    //     {
    //         std::cout << std::endl << "FIXME: we have more wells thatn expected";
    // //         int tt;
    // //         std::cin >> tt;
    //         exit(1);
    //     }
    // Find the locations but only keep the maxima of each region
    for( int i=0; i<(*locations).size(); ++i ){
        int xx = (*locations).at(i).at(0);
        int yy = (*locations).at(i).at(1);
        int type = (*locations).at(i).at(2);
        int value = (*locations).at(i).at(3); // to implement in case of multiple wells
        int col = (double)xx/ ((double)size[0]/(double)wellPerRow);
        int row = (double)yy/ ((double)size[1]/(double)wellPerRow);

        // I will beter clean the max imag, bef finding maxima
        //         if( xx < minBorder || yy < minBorder || xx > (int)size[0]-(int)minBorder || yy > (int)size[1]-(int)minBorder )
        //             continue;

        if( value > allWells.at(row).at(col).at(3) || allWells.at(row).at(col).at(3) == -1000 ){
            allWells.at(row).at(col).at(0) = xx;
            allWells.at(row).at(col).at(1) = yy;
            allWells.at(row).at(col).at(2) = type;
            allWells.at(row).at(col).at(3) = value;
            std::cout << std::endl << "WELL, row: " << row << ", col: " << col << ", x: " << xx << ", yy: " << yy;
        }
    }
    return allWells;
}


void FindShift( ImageCF2::Pointer*sliceImageFFT, ImageCF2::Pointer*sliceImagePreviousFFT, std::vector< int >*shift, bool flagRotate ){

    MultiplyImageFilterTypeCF2CF2::Pointer multiplyFilterCF2CF2Pat1 = MultiplyImageFilterTypeCF2CF2::New ();
    multiplyFilterCF2CF2Pat1->SetInput1( (*sliceImageFFT) );
    multiplyFilterCF2CF2Pat1->SetInput2( (*sliceImagePreviousFFT) );
    multiplyFilterCF2CF2Pat1->Update();

    // ICorr
    ImageF2::Pointer imgCorr = IFFT2< ImageCF2, ImageF2 >( multiplyFilterCF2CF2Pat1->GetOutput() );

    ImageF2::SizeType size = imgCorr->GetLargestPossibleRegion().GetSize();
    ImageF2::PixelType * imgCorrPtr = imgCorr->GetBufferPointer();

    // Find max of the two
    double max= -10000;
    unsigned long long indx;
    int ii_max, jj_max;
    for( int j=0; j<size[1]; ++j )
    {
        for( int i=0; i<size[0]; ++i )
        {
            indx = i+j*size[0];
            if( imgCorrPtr[ indx ] > max ){
                max = imgCorrPtr[ indx ];
                ii_max = i;
                jj_max = j;
            }
        }
    }
    int minShift_ii, minShift_jj;
    if( ii_max < size[0]/2 )
    {
        if( flagRotate == 1 )
            minShift_ii = ii_max+1;
        else
            minShift_ii = -(ii_max+1);
    }
    else
    {
        if( flagRotate == 1 )
            minShift_ii = -((int)size[0]-1-ii_max);
        else
            minShift_ii = ((int)size[0]-1-ii_max);
    }
    if( jj_max < size[1]/2 )
    {
        if( flagRotate == 1 )
            minShift_jj = jj_max+1;
        else
            minShift_jj = -(jj_max+1);
    }
    else
    {
        if( flagRotate == 1 )
            minShift_jj = -((int)size[1]-1-jj_max);
        else
            minShift_jj = ((int)size[1]-1-jj_max);
    }


    (*shift).resize(2);
    (*shift).at(0) = minShift_ii;
    (*shift).at(1) = minShift_jj;
}


void FindInitialWells( ImageF2::Pointer INCCPat1, ImageF2::Pointer INCCPat2, int wellSize, int minBorder, std::vector< std::vector< double > > *locations, ImageCHAR2::Pointer* WellCenterImage ){
    // This is the main function, it locates the wells based on the NCC coefficient

    ImageF2::SizeType size = INCCPat1->GetLargestPossibleRegion().GetSize();
    int wellSizeRadius = (wellSize-1)/2;
    ImageF2::PixelType * INCCPat1Ptr = INCCPat1->GetBufferPointer();
    ImageF2::PixelType * INCCPat2Ptr = INCCPat2->GetBufferPointer();

    // Find max of the two
    ImageF2::Pointer INCCMax = GetITKImageOfSize<ImageF2>( size );
    ImageF2::PixelType * INCCMaxPtr = INCCMax->GetBufferPointer();
    unsigned long long indx;
    for( int j=0; j<size[1]; ++j )
    {
        for( int i=0; i<size[0]; ++i )
        {
            indx = i+j*size[0];
            if( i < minBorder || j < minBorder || i > (int)size[0]-(int)minBorder || j > (int)size[1]-(int)minBorder ){
                INCCPat2Ptr[ indx ] = -1000;
                INCCPat1Ptr[ indx ] = -1000;
            }

            if( INCCPat1Ptr[ indx ] > INCCPat2Ptr[ indx ] )
                INCCMaxPtr[ indx ] = INCCPat1Ptr[ indx ];
            else
                INCCMaxPtr[ indx ] = INCCPat2Ptr[ indx ];
        }
    }

    double* max_row;
    max_row = new double[size[0]*size[1]];

    //     #pragma omp  parallel for collapse(2) schedule(dynamic, 1)
    for( int j=0; j<size[1]; ++j )
    {
        for( int i=0; i<size[0]; ++i )
        {
            double max = -1000;
            for( int i_in=std::max( i-wellSizeRadius, 0); i_in<= std::min( (int)(size[0]-1), i+wellSizeRadius ); ++i_in )
            {
                if( INCCMaxPtr[ j*size[0] + i_in ] > max )
                {
                    max = INCCMaxPtr[ j*size[0]+i_in ];
                }
            }
            max_row[ j*size[0]+i ] = max;
        }
    }

    double* max_col;
    max_col = new double[size[0]*size[1]];

    //     #pragma omp  parallel for collapse(2) schedule(dynamic, 1)
    for( int i=0; i<size[0]; ++i )
    {
        for( int j=0; j<size[1]; ++j )
        {
            double max = -1000;
            for( int j_in=std::max( j-wellSizeRadius, 0); j_in<= std::min( (int)(size[1]-1), j+wellSizeRadius ); ++j_in )
            {
                if( max_row[ j_in*size[0] + i ] > max )
                {
                    max = max_row[ j_in*size[0]+i ];
                }
            }
            max_col[ j*size[0]+i ] = max;
            //             INCCPat1Ptr[ j*size[0]+i ] = max;
        }
    }

    // Find the max of the two
    (*WellCenterImage) = GetITKImageOfSize<ImageCHAR2>( size );
    ImageCHAR2::PixelType * WellCenterPtr = (*WellCenterImage)->GetBufferPointer();
    for( int j=0+minBorder; j<size[1]-minBorder; ++j )
    {
        for( int i=0+minBorder; i<size[0]-minBorder; ++i )
        {
            if( max_col[j*size[0]+i] <= INCCMaxPtr[j*size[0]+i] )
            {
                // Locations
                std::vector< double > CoordWell;
                CoordWell.resize(4);
                CoordWell.at(0) = i; // x
                CoordWell.at(1) = j; // y
                if( INCCPat1Ptr[j*size[0]+i] > INCCPat2Ptr[j*size[0]+i] )
                {
                    CoordWell.at(2) = 1; // Type
                    CoordWell.at(3) = INCCPat1Ptr[j*size[0]+i]; // value
                    (*locations).push_back( CoordWell );
                    // This is just for debug
                    for( int i_in=i-40; i_in<i+40; ++i_in )
                    {
                        for( int j_in=j-40; j_in<j+40; ++j_in )
                        {
                            WellCenterPtr[ j_in*size[0]+i_in ] = 100;
                        }
                    }
                    //                     std::cout << std::endl << i << " " << j;
                }
                else
                {
                    CoordWell.at(2) = 2; // Type
                    CoordWell.at(3) = INCCPat2Ptr[j*size[0]+i]; // Value
                    (*locations).push_back( CoordWell );
                    // This is just for debug
                    for( int i_in=i-40; i_in<i+40; ++i_in )
                    {
                        for( int j_in=j-40; j_in<j+40; ++j_in )
                        {
                            WellCenterPtr[ j_in*size[0]+i_in ] = 200;
                        }
                    }
                    //                     std::cout << std::endl << i << " " << j;
                }
            }
        }
    }
    //     int yy;
    //     std::cin >> yy;
    delete [] max_row;
    delete [] max_col;

}




ImageF2::Pointer SubsBack( ImageF2::Pointer imgIn, double sigma ){

    //     FilterType::Pointer smoothingRecursiveGaussianImageFilter = FilterType::New();
    //     smoothingRecursiveGaussianImageFilter->SetInput(reader->GetOutput());
    //     smoothingRecursiveGaussianImageFilter->SetSigma(sigma);
    //     smoothingRecursiveGaussianImageFilter->Update();

    ImageF2::SizeType size = imgIn->GetLargestPossibleRegion().GetSize();
    typedef itk::SmoothingRecursiveGaussianImageFilter< ImageF2, ImageF2 >  filterType;
    // Create and setup a Gaussian filter
    filterType::Pointer gaussianFilter = filterType::New();
    gaussianFilter->SetInput( imgIn );
    gaussianFilter->SetSigma( sigma );
    try
    {
        gaussianFilter->Update();
    }
    catch(itk::ExceptionObject &err)
    {
        std::cerr << "ExceptionObject caught!" <<std::endl;
        std::cerr << err << std::endl;
    }

    ImageF2::Pointer imgOut = GetITKImageOfSize<ImageF2>( size );

    ImageF2::PixelType * imgInPtr = imgIn->GetBufferPointer();
    ImageF2::PixelType * imgOutPtr = imgOut->GetBufferPointer();
    ImageF2::PixelType * imgGaussPtr = gaussianFilter->GetOutput()->GetBufferPointer();


    unsigned long long indx;
    for( unsigned long long ii = 0; ii < size[0]; ++ ii ){
        for( unsigned long long jj = 0; jj < size[1]; ++ jj ){
            indx = ii + size[0]*jj;
            imgOutPtr[ indx ] = std::max( (double)0, imgInPtr[ indx ] - imgGaussPtr[ indx ] );
        }
    }
    return imgOut;

}


void NCCFun( ImageF2::Pointer sliceImage, ImageF2::Pointer pattern_1, ImageF2::Pointer pattern_2, ImageF2::Pointer *outPat1, ImageF2::Pointer *outPat2 ){

    clock_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13;

    t1 = clock();

    ImageF2::SizeType sizeSliceImage = sliceImage->GetLargestPossibleRegion().GetSize();

    ImageF2::SizeType sizePattern_1 = pattern_1->GetLargestPossibleRegion().GetSize();
    ImageF2::SizeType sizePattern_2 = pattern_2->GetLargestPossibleRegion().GetSize();

    // Multipy
    ImageF2::Pointer imgCorrPat1, imgCorrPat2;
    {
        ImageCF2::Pointer imageStackPadPat1FFT = FFT2< ImageF2, ImageCF2 >( PadImagePost< ImageF2 >( sliceImage, sizePattern_1[0]-1, sizePattern_1[1]-1 ) );
        MultiplyImageFilterTypeCF2CF2::Pointer multiplyFilterCF2CF2Pat1 = MultiplyImageFilterTypeCF2CF2::New ();
        multiplyFilterCF2CF2Pat1->SetInput1( imageStackPadPat1FFT );
        multiplyFilterCF2CF2Pat1->SetInput2( FFT2< ImageF2, ImageCF2 >( PadImagePost< ImageF2 >( RotateImage180< ImageF2 > ( pattern_1 ), sizeSliceImage[0]-1, sizeSliceImage[1]-1 ) ) );      
        multiplyFilterCF2CF2Pat1->Update();

        MultiplyImageFilterTypeCF2CF2::Pointer multiplyFilterCF2CF2Pat2 = MultiplyImageFilterTypeCF2CF2::New ();
        multiplyFilterCF2CF2Pat2->SetInput1( imageStackPadPat1FFT );
        multiplyFilterCF2CF2Pat2->SetInput2( FFT2< ImageF2, ImageCF2 >( PadImagePost< ImageF2 >( RotateImage180< ImageF2 > ( pattern_2 ), sizeSliceImage[0]-1, sizeSliceImage[1]-1 ) ) );      
        multiplyFilterCF2CF2Pat2->Update();

        // ICorr
        imgCorrPat1 = IFFT2< ImageCF2, ImageF2 >( multiplyFilterCF2CF2Pat1->GetOutput() );
        imgCorrPat1->DisconnectPipeline();
        imgCorrPat2 = IFFT2< ImageCF2, ImageF2 >( multiplyFilterCF2CF2Pat2->GetOutput() );
        imgCorrPat2->DisconnectPipeline();
    }

    // Local local_sum (LocalQSumI2, will contain the result)
    ImageF2::Pointer LocalQSumI2 = LocalSum2< ImageF2 > ( PowerImage< ImageF2 >( sliceImage, 2 ), sizePattern_1 );
    ImageF2::Pointer LocalSumI2 = LocalSum2< ImageF2 > ( sliceImage, sizePattern_1 );

    ImageF2::Pointer stdImg = LocalQSumI2;

    // Std (LocalQSumI2, will contain the result)
    ImageF2::PixelType * LocalQSumI2Ptr = LocalQSumI2->GetBufferPointer();
    ImageF2::PixelType * stdImgPtr = stdImg->GetBufferPointer();
    ImageF2::PixelType * LocalSumI2Ptr = LocalSumI2->GetBufferPointer();
    ImageF2::SizeType sizeLQS2 = stdImg->GetLargestPossibleRegion().GetSize();

    // stdI=sqrt(max(LocalQSumI2-(LocalSumI2.^2)/numel(T),0) );
    unsigned long long indx;
    for( unsigned long long ii = 0; ii < sizeLQS2[0]; ++ ii ){
        for( unsigned long long jj = 0; jj < sizeLQS2[1]; ++ jj ){
            indx = ii + sizeLQS2[0]*jj;
            stdImgPtr[ indx ] = std::sqrt( std::max( (double)LocalQSumI2Ptr[ indx ] - (double)LocalSumI2Ptr[ indx ]*(double)LocalSumI2Ptr[ indx ]/(double)(sizeLQS2[0]*sizeLQS2[1]), (double)0 ) );
            //             stdImgPtr[ indx ] = LocalQSumI2Ptr[ indx ];
        }
    }

    // stdT
    double stdTPat1 = std::sqrt( (double)sizePattern_1[0]*sizePattern_1[1]-1)*std::sqrt( variance< ImageF2 >( pattern_1 ) );
    double stdTPat2 = std::sqrt( (double)sizePattern_1[0]*sizePattern_1[1]-1)*std::sqrt( variance< ImageF2 >( pattern_2 ) );

    // meanIT (meanIT=LocalSumI2*sum(T(:))/numel(T);)
    // LocalSumI2, will be used to store the results
    ImageF2::PixelType * pattern_1Ptr = pattern_1->GetBufferPointer();
    ImageF2::PixelType * pattern_2Ptr = pattern_2->GetBufferPointer();
    double avrPattern1 = 0;
    double avrPattern2 = 0;
    for( unsigned long long ii = 0; ii < sizePattern_1[0]; ++ ii ){
        for( unsigned long long jj = 0; jj < sizePattern_1[1]; ++ jj ){
            indx = ii + sizePattern_1[0]*jj;
            avrPattern1 += pattern_1Ptr[indx];
            avrPattern2 += pattern_2Ptr[indx];
        }
    }
    avrPattern1 /= (double)sizePattern_1[0]*(double)sizePattern_1[1];
    avrPattern2 /= (double)sizePattern_1[0]*(double)sizePattern_1[1];

    ImageF2::Pointer INCCPat1 = GetITKImageOfSize<ImageF2>( sizeSliceImage );
    ImageF2::Pointer INCCPat2 = GetITKImageOfSize<ImageF2>( sizeSliceImage );

    ImageF2::PixelType * INCCPat1Ptr = INCCPat1->GetBufferPointer();
    ImageF2::PixelType * imgCorrPat1Ptr = imgCorrPat1->GetBufferPointer();

    ImageF2::PixelType * INCCPat2Ptr = INCCPat2->GetBufferPointer();
    ImageF2::PixelType * imgCorrPat2Ptr = imgCorrPat2->GetBufferPointer();

    unsigned long long indx_in;
    for( unsigned long long ii = 0; ii < sizeSliceImage[0]; ++ ii ){
        for( unsigned long long jj = 0; jj < sizeSliceImage[1]; ++ jj ){
            indx = ii + sizeSliceImage[0]*jj;
            indx_in = (ii+( sizePattern_1[0]-1)/2 ) + sizeLQS2[0]*(jj+( sizePattern_1[1]-1)/2 );
            INCCPat1Ptr[ indx ] = 0.5+ ( (double)imgCorrPat1Ptr[ indx_in ] - (double)LocalSumI2Ptr[ indx_in ]*avrPattern1) / ( 2*stdTPat1*std::max( (double)stdImgPtr[ indx_in ], (double)stdTPat1/(double)1.0e5 ) );
            INCCPat2Ptr[ indx ] = 0.5+ ( (double)imgCorrPat2Ptr[ indx_in ] - (double)LocalSumI2Ptr[ indx_in ]*avrPattern2) / ( 2*stdTPat2*std::max( (double)stdImgPtr[ indx_in ], (double)stdTPat2/(double)1.0e5 ) );
        }
    }

    INCCPat1->DisconnectPipeline();
    INCCPat2->DisconnectPipeline();

    (*outPat1) = INCCPat1;
    (*outPat2) = INCCPat2;

}

// #endif
