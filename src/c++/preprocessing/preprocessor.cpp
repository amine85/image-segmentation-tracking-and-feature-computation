#include "preprocessor.h"

Preporcessor::Preporcessor()
{
    this->backgroundSubtract = true;
    this->linearUnmix = false;    
}
Preporcessor::~Preporcessor()
{
}

void Preporcessor::setInputImage(T * image1, T * image2)
{
    this->image1 = image1;
    this->image2 = image2;
}
void Preporcessor::setDimensions(size_t nr,size_t nc, size_t nt)
{
    this->nr = nr;
    this->nc = nc;
    this->nt = nt;
}
void Preporcessor::runPreprocessor()
{
    if(this->linearUnmix)
        this->runLinearUnmixing();
    if(this->backgroundSubtract)
    {
        this->runBackgroundSubtraction(this->image1);
        this->runBackgroundSubtraction(this->image2);
    }
}
    
void Preporcessor::runLinearUnmixing()    
{
    if(!image1 || !image2)
    {
        std::cout<<" you did not provide one or two input image, the output of linear unmixing will be the same as the input "<<std::endl;
        return;        
    }
    if(mixing_matrix.empty())
    {
        std::cout<<" mixing matrix is empty, the output of linear unmixing will be the same as the input "<<std::endl;
        return;
    }
    
    /***************************************************************************************************************/
    // unmix by matrix inversion    
    /***************************************************************************************************************/  
    
    float a11 =  mixing_matrix[0][0];
    float a21 =  mixing_matrix[1][0];
    float a12 =  mixing_matrix[0][1];
    float a22 =  mixing_matrix[1][1];
    
    float norm = sqrt(a11*a11 + a21*a21);
    a11 /= (norm+EPSILON);
    a21 /= (norm+EPSILON);
    norm = sqrt(a12*a12 + a22*a22);
    a12 /= (norm+EPSILON);
    a22 /= (norm+EPSILON);
    
    float detM = a11*a22 - a21*a12;
    mixing_matrix[0][0] =  a22/(detM+EPSILON);
    mixing_matrix[1][0] = -a21/(detM+EPSILON);
    mixing_matrix[0][1] = -a12/(detM+EPSILON);
    mixing_matrix[1][1] =  a11/(detM+EPSILON);
    
    std::cout<<"Normalized Matrix:\n";
    std::cout<<  a11 << "\t"<< a12 <<"\n";
    std::cout<<  a21 << "\t"<< a22 <<"\n";  
    std::cout<< std::endl;
    std::cout<<"Inverted Matrix:\n";
    std::cout<<  mixing_matrix[0][0] << "\t"<< mixing_matrix[0][1]<<"\n";
    std::cout<<  mixing_matrix[1][0] << "\t"<< mixing_matrix[1][1]<<"\n";  
    std::cout<< std::endl;
        
    
    // declare loop variables
    size_t nPixels = (nr*nc*nt);
    float outVal01 = 0.0;
    float outVal02 = 0.0;
    
    float inVal01 = 0.0;
    float inVal02 = 0.0;

    
    #pragma omp  parallel for schedule(dynamic, 1)
    for(size_t index=0; index<nPixels; ++index)
    {   
      
        inVal01 = (float) image1[index];
        inVal02 = (float) image2[index];
        
        outVal01 = (inVal01*mixing_matrix[0][0] +  inVal02*mixing_matrix[0][1]);
        outVal02 = (inVal01*mixing_matrix[1][0] +  inVal02*mixing_matrix[1][1]);
        
        if(outVal01>0)
        {
          image1[index] = (InputPixelType16) outVal01;       
        }
        else
        {
          image1[index] = 0;
        }
        
        if(outVal02>0)
        {
          image2[index] = (InputPixelType16) outVal02;
        }
        else
        {
          image2[index] = 0;
        }
    }
        
    return;
}
  
void Preporcessor::runBackgroundSubtraction(T * image)
{
    
    std::cout<< "subtracting background ..." << std::endl;
    #pragma omp  parallel for schedule(dynamic, 1)
    for(size_t t = 0; t<nt; ++t)
    {
        
        typename U::Pointer input2Dimage = U::New();
        typename U::IndexType start;
        start[0] = 0;  // size along X
        start[1] = 0;  // size along Y
        
        
        typename U::SizeType  size;
        size[0] = this->nr;  // size along X
        size[1] = this->nc;  // size along Y

        typename U::RegionType region;
        region.SetSize( size );
        region.SetIndex( start );

        input2Dimage->SetRegions( region );
        input2Dimage->Allocate();
        input2Dimage->FillBuffer(0);
        input2Dimage->Update();

        T * input2DImagePtr = input2Dimage->GetBufferPointer();
        size_t offset = t*(numRows*numCols);
        for(size_t index = 0; index < (numRows*numCols); ++index)
        {
            input2DImagePtr[index] = image[index+offset];
        }
        
        typename U::SpacingType spacing;
        spacing[0] = 1;
        spacing[1] = 1;
        input2Dimage->SetSpacing(spacing);
        typedef itk::SmoothingRecursiveGaussianImageFilter<U,U> FilterType;
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput(input2Dimage);
        filter->SetSigma(this->sigmab);
        try
        {
            filter->Update();
        }
        catch(itk::ExceptionObject &err)
        {
            std::cerr << "ExceptionObject caught!" <<std::endl;
            std::cerr << err << std::endl;
        }
        typename U::Pointer backg_image = filter->GetOutput();

        
        T * backgImagePtr = backg_image->GetBufferPointer();
        
        
        float diffVal = 0.0;
        for(size_t index = 0; index < (numRows*numCols); ++index)
        {
            diffVal = (float)input2DImagePtr[index]-(float)backgImagePtr[index];
            if(diffVal>0)
                image[index+offset] = (InputPixelType16) diffVal;
            else
                image[index+offset] = 0;
        }  
    }
    std::cout<< "finished subtracting background ..." << std::endl;
}




