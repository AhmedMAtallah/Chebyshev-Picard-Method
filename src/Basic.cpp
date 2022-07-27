/*
 
     Copyright c SPACE Lab 2019
    https://space.sdsu.edu/ 
San Diego State University (SDSU), Aerospace department 


	This  software is provided on an "as is" basis without warranty of any kind, express or implied.
	Under no circumstances and under no legal theory, whether in tort, contract, or otherwise, shall
	SPACE lab, or SDSU be liable to you or to any other person for any indirect,
	special, incidental, or consequential damages of any character including, without limitation, damages
	for software errors, work stoppage, computer failure or malfunction, loss of goodwill or for any and
	all other damages or losses.

Code Deveopers : Ahmed Atallah and Ahmad Bani Younes

   File name     : Basic.ccp
   Subject       : Example is used in ....................
   Description   :  Source file for the Basic functions

   
   Date Modified : 09/23/2019
     
*/




#include "../include/Basic.h"
#include "../include/const.h"


int iDivUp(int a, int b)
{
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

double  deg2rad(double angleInDegrees)
{
/*% DEG2RAD Convert angles from degrees to radians.
%   DEG2RAD(X) converts angle units from degrees to radians for each
%   element of X.
%
%   See also RAD2DEG.

% Copyright 2015 The MathWorks, Inc.
*/
	double angleInRadians;
//if isfloat(angleInDegrees)
    angleInRadians = (C_PI/180) * angleInDegrees;

	return angleInRadians;
}
 void trans( double *A,double *AT,const int m, const int n, const int ld)
 {
	 for (int i =0 ; i<m; i++)
		for (int j=0;j<n; j++)
			AT[IDX2F(j+1,i+1,n)]=A[IDX2F(i+1,j+1,m)];
	 
 }

void matmul(double* A, double* B, double* OUT, int m, int n, int q)
{
	for(int i=0;i<m; i++){
		for(int j=0;j<q;j++){
			double sum = {0.0};
			for(int j1=0;j1<n;j1++)
				sum += A[IDX2F(i+1,j1+1,m)]*B[IDX2F(j1+1,j+1,n)];
			OUT[IDX2F(i+1,j+1,m)] = sum;
		 }
   	}
}

double pos(double a)
{

	if (a<0)
		return -a;
	else
		return a;

}

/*
  Print out the contents of a section of array for debugging

  "A"  is a pointer to the array
  "m"  is the number of rows from the top to print
  "n"  is the number of cols from the left to print
  "ld" is the number of rows of the data container

*/
void printArray( const double* A, const int m, const int n, const int ld )
{
    int ii, jj;
    for ( ii=1; ii<=m; ii++ ) {
        for ( jj=1; jj<=n; jj++ ) {
            printf( "%1.8e\t", A[IDX2F(ii,jj,ld)] );
        }
        printf( "\n" );
    }
}


/*
  Write the contents of an array to a file

  "A"  is a pointer to the array
  "m"  is the number of rows from the top to print
  "n"  is the number of cols from the left to print
  "ld" is the number of rows of the data container
  "fname" is what to call the output file

*/
void writeArray( const double* A, const int m, const int n, const int ld,
                 const char* fname )
{

    FILE* f = fopen( fname, "w" );
    int ii, jj;
    for ( ii=1; ii<=m; ii++ ) {
        for ( jj=1; jj<=n; jj++ ) {
            fprintf( f, "%1.10f\t", A[IDX2F(ii,jj,ld)] );
        }
        fprintf( f, "\n" );
    }
    fclose( f );
}


// Function to round - off the number 
double Round_off(double NN, double n) 
{ 
    int h; 
    double l, a, b, c, d, e, i, j, m, f, g; 
    b = NN; 
    c = floor(NN); 
  
    // Counting the no. of digits to the left of decimal point 
    // in the given no. 
    for (i = 0; b >= 1; ++i) 
        b = b / 10; 
  
    d = n - i; 
    b = NN; 
    b = b * pow(10, d); 
    e = b + 0.5; 
    if ((float)e == (float)ceil(b)) { 
        f = (ceil(b)); 
        h = f - 2; 
        if (h % 2 != 0) { 
            e = e - 1; 
        } 
    } 
    j = floor(e); 
    m = pow(10, d); 
    j = j / m; 
	return j;
}
