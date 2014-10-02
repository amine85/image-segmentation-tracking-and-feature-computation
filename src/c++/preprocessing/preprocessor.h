#ifndef __PREPROCESSOR_H
#define __PREPROCESSOR_H

#include "helpers.h"


template <typename T,typename U>
class Preporcessor
{
     
public:
    Preporcessor();
    ~Preporcessor();    
    
    // set processing steps and parameters
    void setBackgroundSubtraction(bool value = true){this->backgroundSubtract = value;};
    void setLinearUnmixing(bool value = true){this->linearUnmix = value;};
    void setMixingMatrix(std::vector< std::vector<float> > mixing_matrix){this->mixing_matrix = mixing_matrix;};
    void setGaussianSigma(float s = 40.0){this->sigmab = s;};    
    void setDimensions(size_t nr,size_t nc, size_t nt);
    
    // set input and run, output will overwrite the input image
    void setInputImage(T * image1, T * image2);    // this input means the processing will be ran in place
    void runPreprocessor(void);
    
private:
    // processing flags
    bool backgroundSubtract;
    bool linearUnmix;
    
    // parameters and inputs
    T * image1;
    T * image2;
    std::vector< std::vector<float> > mixing_matrix;
    float sigmab;
    size_t nr;
    size_t nc;
    size_t nt;
    
    // processing functions
    void runLinearUnmixing(void);
    void runBackgroundSubtraction(void);
    
}

#endif