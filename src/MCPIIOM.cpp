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



   File name     : MCPIOM.cpp
   Subject       : Example is used in ....................
   Description   :  Source Functions of MCPI

   
   Date Modified : 09/23/2019
     
*/


#include "../include/MCPIIOM.h"
#include "../include/const.h"
//#include "FEM.h"
#include "../include/Basic.h"

//#define N  40			// number of nodes 
//#define M  (N+1)				// number of nodes + 1
void MCPI_CoeffsI(int N, int M, double*Im)
{
	int i, j;
	double *W=new double [M*M];
	double *T=new double [(N+1)*M];
	double * TT=new double [M*(N+1)];
	double * T2=new double [(N+2)*(M)];
	double * T2Z= new double [(N+2)*M];
 	double * Z =new double [(M+2)*M];
	double * V =new double [(N+1)*(N+2)];
	double *I_N=new double [(N+1)*(N+2)];

/*
	double W[M*M] = {0.0};
	double T[(N+1)*M] = {0.0};
	double TT[M*(N+1)] = {0.0};
	double T2[(N+2)*(M)] = {0.0};
	double T2Z[(N+2)*M] = {0.0};
	double Z[(M+2)*M] = {0.0};
	double V[(N+1)*(N+2)] = {0.0};
	double I_N[(N+1)*(N+2)]={0.0};
*/	
	double * tau=new double [M];
	
	for (int i =0; i<M; i++)
		tau[i]=cos(i*C_PI/N+pi);
	
	// BUILD W MATRIX (symmetric)
	for (int i=0;i<M*M;i++)
		W[i]=0.0;

	W[IDX2F(1,1,M)] = 0.5;
	for ( i=2; i<=M; i++ ) {
		W[IDX2F(i,i,M)] = 1.0;
	}
	W[IDX2F(M,M,M)] = 0.5;
	
	//BUILD T MATRIX (symmetric)
	for ( j=0; j<M; j++ ) {
		for ( i=0; i<=N; i++ ) {
			T[IDX2F(i+1,j+1,N+1)] = cos((double)i*acos(tau[j]));
		}
	}
//	printf("   T  \n");
//	printArray(T,M,M,M);
	
	// BUILD T MATRIX (symmetric)
	for ( j=0; j<N+1; j++ ) 
	{
		for ( i=0; i<M; i++ ) 
		{
			TT[IDX2F(i+1,j+1,M)] = cos((double)j*acos(tau[i]));
		}
	}
//	printf("   TT  \n");
//	printArray(TT,M,M,M);
	
				// BUILD T MATRIX (symmetric)
	for ( j=0; j<M; j++ ) {
		for ( i=0; i<N+2; i++ ) {
			T2[IDX2F(i+1,j+1,N+2)] = cos((double)i*acos(tau[j]));
		}
	}
	
//	printf("   T2 \n");
//	printArray(T2,M+1,M,M+1);
	
			// BUILD T MATRIX (symmetric)
	for ( j=0; j<M; j++ ) {
		for ( i=0; i<N+2; i++ ) {
			T2Z[IDX2F(i+1,j+1,N+2)] = cos((double)i*acos(tau[j]))-pow(-1,i+2);
		}
	}
	
//	printf("   T2Z  \n");
//	printArray(T2Z,N+2,M,N+2);
				// BUILD T MATRIX (symmetric)
	for ( j=0; j<M; j++ ) {
		for ( i=0; i<N+2; i++ )  {
			Z[IDX2F(i+1,j+1,N+2)] = pow(-1,i+2);
		}
	}
	
	//	printf("  Z  \n");
//	printArray(Z,N+2,N+1,N+2);
	//void trans( double *A,const int m, const int n, const int ld);
	// BUILD V MATRIX (symmetric)
	double vElem = 1.0/(double)N;
	for (int  i =0; i<(N+1)*(N+1); i++)
		V[i]=0.0;

	V[IDX2F(1,1,M)] = vElem;
	V[IDX2F(M,M,M)] = vElem;
	for ( i=2; i<=N; i++ ) {
		V[IDX2F(i,i,M)]	= 2.0*vElem;
	}
	for (int i=0; i<(N+1)*(N+2); i++)
		I_N[i]=0.0;
	
	I_N[IDX2F(1,2,M)]=1.0;
	I_N[IDX2F(2,1,M)]=0.25;
	I_N[IDX2F(2,3,M)]=0.25;
	
	for (int ii=3;ii<=M; ii++)
	{
		I_N[IDX2F(ii,ii-1,N+1)]=-0.5/((double)ii-2);
		I_N[IDX2F(ii,ii+1,N+1)]=0.5/((double)ii);		
	}
	
	//printf( " I n \n");
	//printArray(I_N,N+1,N+2,N+1);
	double WT[M*(N+1)];
	// Building Cx & Ca matrices
	matmul(W,TT,WT,M,N+1,M);  // Cx = T*W
	double WTV[M*(N+1)];
	matmul(WT,V,WTV,M,N+1,M);    // TV  = T*V
	// printf(" WTV  \n");
   // printArray(WTV,M,N+1,M);
//	double * IT;
 //   matmul(I_N,,STV,M,M,M);  // STV = S*T*V
    double ITZ[(N+1)*M];
    matmul(I_N,T2Z,ITZ,N+1,N+2,M);  // Ca  = R*S*T*V
   // printf(" ITZ \n");
   // printArray(ITZ,N+1,M,N+1);
	matmul(WTV,ITZ,Im,M,M,N+1);  // CxCa  = Cx*Ca
//	printf(" Im \n");
//    printArray(Im,M,M,M);

	//double ImT[M*M];

//TZI*V*T*W
//	printf(" Im \n");
  //  printArray(Im,M,M,M);
	
}


void errorAndUpdate(int MM,double timeSub,int Nstates2,double * x0,double *Xo,double *Xn,double* xAdd,double temp)
{
			double Err=0;
		    for(int node=0; node<MM;node++){
				for(int state=0;state<Nstates2;state++){
					// initialize the arrays
				    int indx  = IDX2F(node+1,state+1,MM);
				    Xn[indx]  = 0.0;
				    Err = 0.0;
				    Xn[indx]  = x0[state] + timeSub*xAdd[indx];
				    Err = pos(Xn[indx] - Xo[indx])/max(1.0,pos(Xo[indx])); 
				
					//if (loopCount==14)
				//		printf(" Error  %e \n",Err[indx]);
					if (indx==0)
						temp=Err;
				    if(Err>temp)     // comment this line when "Q-For" loop is used
					{
						temp = Err;  // comment this line when "Q-For" loop is used
					//	printf(" temp %e \n",temp);
					}
				//	Xo[indx]  = 0.0;
				    Xo[indx]  = Xn[indx] ;
				}
			}
}
