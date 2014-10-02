
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

/* function that converts a matrix in 1-D format to a compressed column storage (CCS) */
void matrixToCCSFormat(double matrix[], int n, int irow[], int pcol[], double A[])
{
	int nj;
	int idx1 = 0;
	int idx2 = 0;
	bool firstCol;

	for(int j=0; j<n; ++j){
		nj = n*j;
		firstCol = true;
		for(int i=j; i<n; ++i){
			if(matrix[i+nj]!=0){
				A[idx1] = matrix[i+nj];
				irow[idx1] = i;
				if(firstCol){
					firstCol = false;
					pcol[idx2] = idx1;
					++idx2;
				}
				++idx1;
			}
		}
	}
	
}


/* function that converts a degree matrix in (just diagonal) to a compressed column storage (CCS) */
void degreeToCCSFormat(double degree_matrix[], int n, int irow[], int pcol[], double A[])
{
	for(int i=0; i<n; ++i){
		pcol[i]=i;
		irow[i]=i;
		A[i]=degree_matrix[i];
	}
}

/* function that computes the degree matrix*/
template < typename T >
void computeDegree(T matrix[], T degree_matrix[], int n)
{
	for(int i=0; i<n; ++i){
		T sum = 0.0;
		for(int j=0; j<n; ++j)
			sum+=matrix[i+n*j];
		degree_matrix[i] = sum;
	}
}

/* function that computes the laplacian matrix*/
template < typename T >
void computeLaplacian(T matrix[], T laplacian_matrix[], T degree_matrix[], int n)
{
	for(int i=0; i<n*n; ++i)
		laplacian_matrix[i] = -matrix[i];
	for(int i=0; i<n; ++i)
		laplacian_matrix[i+i*n] += degree_matrix[i];	
}

/* function that computes the laplacian matrix in place*/
template < typename T >
void computeLaplacianInPlace(T matrix[], T degree_matrix[], int n)
{
	for(int i=0; i<n*n; ++i)
		matrix[i] = -matrix[i];
	for(int i=0; i<n; ++i)
		matrix[i+i*n] += degree_matrix[i];	
}
int main()
{
	// define our input array 4x4 matrix
	int N = 16;
	int n = 4;
	double* matrix = new double[N];
	for(int i=0; i<N; ++i) 
		matrix[i] = 0.0;
	
	// diagonal to 1
	for(int i=0; i<n; ++i)
		matrix[i+n*i] = 1.0;
	matrix[n]=0.8;
	matrix[1]=0.8;
	matrix[3+2*n]=0.8;
	matrix[2+3*n]=0.8;
	
	matrix[2+n]=0.01;
	matrix[1+2*n]=0.01;
	
	printMatrix<double>(matrix,n);
	
		
	// let's compute the degree and the laplacian
	double* degree_matrix = new double[n]; // just store the diagonal
	computeDegree<double>(matrix,degree_matrix,n);
	std::cout<<"degree matrix: \n";
	printArray<double>(degree_matrix,n);
	
// 	double* laplacian_matrix = new double[N];
// 	computeLaplacian<double>(matrix,laplacian_matrix,degree_matrix,n);
// 	std::cout<<"laplacian matrix: \n";
// 	printMatrix<double>(laplacian_matrix,n);
	
	computeLaplacianInPlace<double>(matrix,degree_matrix,n);
	std::cout<<"laplacian matrix in place: \n";
	printMatrix<double>(matrix,n);
		
	int nnz = 0; // number of non-zero elements
	for(int i=0; i<N; ++i) 
		if(matrix[i]!=0)
			++nnz;
	std::cout<<"number of non-zeros: "<<nnz<<std::endl;
    
   
	
	// let's initialize arrays for CCS format
	int nnzA = (nnz+n)/2;

	double* A = new double[nnzA];
	int* irowA = new int[nnzA];
	int* pcolA = new int[n];
	matrixToCCSFormat(matrix, n, irowA, pcolA, A);
	printArray<double>(A,nnzA);
	printArray<int>(irowA,nnzA);
	printArray<int>(pcolA,n);
	
	// degree matrix
	int nnzB = n;
	double* B = new double[nnzB];
	int* irowB = new int[nnzB];
	int* pcolB = new int[nnzB];
	degreeToCCSFormat(degree_matrix, n, irowB, pcolB, B);
	std::cout<<"degree matrix: \n";
	printArray<double>(B,nnzB);
	printArray<int>(irowB,nnzB);
	printArray<int>(pcolB,nnzB);
	
	// define and solve the eigen problem
	double* EigVal = new double[n];
	double* EigVec = new double[n*n];
	
	int temp = AREig<double>(EigVal,EigVec,n,nnzA,A,irowA,pcolA,\
							nnzB,B,irowB,pcolB,'L',2);
    
    
  
	printArray<double>(EigVal,2);
    printEigVector<double>(EigVec,n,2);
	
	
	
	
	
	
	
	
	
	
	

}