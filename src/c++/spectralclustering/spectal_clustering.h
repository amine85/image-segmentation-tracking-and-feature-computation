#ifndef __SPECTRAL_CLUTERING_H
#define __SPECTRAL_CLUTERING_H

#include <iostream>

// arpackpp includes
#include "arlsmat.h"
#include "areig.h"



/* function that prints matrix arrays */
template< typename T >
void printMatrix(T a[], int n)
{
    for(int i=0; i<n; ++i){
        for(int j=0; j<n; ++j){
            std::cout<< a[i+n*j] << "\t";
        }
        std::cout<<std::endl;
    }   
}

/* function that prints m eigenvectors of dimension n */
template< typename T >
void printEigVector(T eigVec[], int n, int m)
{
    
    for(int i=0; i<n; ++i){
        for(int j=0; j<m; ++j){
            std::cout<< eigVec[i+n*j] << "\t";
        }
        std::cout<<std::endl;
    }   
}
/* function that prints 1-D arrays*/
template< typename T >
void printArray(T a[], int n)
{
    for(int i=0; i<n; ++i)
        std::cout<< a[i] << "\t";
    std::cout<<std::endl;
}


template< typename T >
class SpectralClustering
{
public:
    SpectralClustering();
    ~SpectralClustering();
    
    // input/output functions 
    void setSigma(float s = 4.0){this->sigma = s;};
    void setNumClusters(unsigned int K);
    void compute(void);
    
    // templated functions 
    void setInputImage(T * image, unsigned int N);
    void getOutputImage(T * image);
    
    
private:

    void normalMatrixToCCSFormat(double matrix[], int n, int irow[], int pcol[], double A[]);        // n = dimension of the matrix
    void degreeMatrixToCCSFormat(double degree_matrix[], int n, int irow[], int pcol[], double A[]);
    void getColRowIndex(unsigned int index, unsigned int * i, unsigned int * j, unsigned int dim);
    void computeDegree(double matrix[], double degree_matrix[], int n);
    void computeLaplacianInPlace(double matrix[], double degree_matrix[], int n)
        
    // templated functions 
    void computeDimension(T * image);
    void computeDistanceMatrix(T * image, unsigned int N, double distanceMatrix[]);    
    
    
    // private members 
    unsigned int mat_dim;     // dimension of the problem (distance matrix)
    unsigned int img_dim;     // input image dimension
    
    unsigned int num_clus;     // number of clusters
    float sigma; 
    T * image;    
};

















#endif