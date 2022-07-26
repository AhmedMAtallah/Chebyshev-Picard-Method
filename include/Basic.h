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

   File name     : Basic.h
   Subject       : Example is used in ....................
   Description   :  Header for the Basic functions

   
   Date Modified : 09/23/2019
     

     
*/

#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <time.h>
#include <stdlib.h>

// cute macro for looking up indices using "1" based addressing ( a la fortran and matlab )
// and column-major array format [Note ld = Num of rows]
#define IDX2F(i,j,ld) ((((j)-1)*(ld))+((i)-1))
#define IDX3F(i,j,k,ld1,ld2) ((((k)-1)*(ld1*ld2))+(((j)-1)*(ld1))+((i)-1))

//#define N  20			// number of nodes 
//const int M=21;			// number of nodes + 1
// define constants 
#define Pi  3.1415926535897932
#define MU  3.986004418000000be+014
#define REQ 6378137
#define omega  7292115.0e-011 // radians/sec (Earth's angular velocity)
//#define DEG 20
//#define Max_Degree 300 //110
// mcpi parameters

//#define Nstates  6
//#define MaxIter 100 // Maximum Iteration - when the loop does not converege
//#define mcpi_tol 1.0e-13 //1.0e-10;

using namespace std;

//Round a / b to nearest higher integer value
int iDivUp(int a, int b);
double  deg2rad(double angleInDegrees);
double pos(double a);
extern "C" void trans( double *A,double *AT,const int m, const int n, const int ld);
extern "C" void matmul(double* A, double* B, double* OUT, int m, int n, int q);

double Round_off(double NN, double n);

/*
  Write the contents of an array to a file

  "A"  is a pointer to the array
  "m"  is the number of rows from the top to print
  "n"  is the number of cols from the left to print
  "ld" is the number of rows of the data container
  "fname" is what to call the output file

*/

void printArray( const double* A, const int m, const int n, const int ld );


void writeArray( const double* A, const int m, const int n, const int ld,
                 const char* fname );
