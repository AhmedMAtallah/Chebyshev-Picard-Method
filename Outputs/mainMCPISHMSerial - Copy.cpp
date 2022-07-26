
#include "SHM.h"
#include "MCPIIOM.h"
#include "Basic.h"
// intial conditions 

/*//--LEO
double init[6] = { -1.921314641994234E6,  -6.382345673714583E6 ,  0.414385389068071E6, 7.363157552624654e3,  -2.277518101688332e3  ,-0.938656130766960e3};
const double tSpan =  5.513674008269939e+03; //5400/10;			// time span of segment integration
const int numSeg =5 ; //15;				// number of segments to integrate
std::string outputPath = "./resultMCPISHMSerialLEO.txt";

//-- SSO
double init[6] = { -6.016660754228067E6,-3.868710263502779E6,0, -0.593051230763435E3,   0.922319797128030E3,   7.383875272629300e3 };
const double tSpan =  6.020800633500418e+03; //5400/10;			// time span of segment integration
const int numSeg =5 ; //15;				// number of segments to integrate
std::string outputPath = "./resultMCPISHMSerialSSO.txt";
*//*
//---GEO
double init[6] = {42121836 ,0 ,0, 0  , 3.077742489283442e3,                  0 };
const double tSpan =  8.616357055057827e+04; //5400/10;			// time span of segment integration
const int numSeg =5 ; //15;				// number of segments to integrate
std::string outputPath = "./resultMCPISHMSerialGEO.txt";
*/

//---HEO

 // -------------Inputs---------------------------//
double init[6] = {6.903780000000001e6,0 ,0, 0,   4.487913207680045e3  , 8.962155583672804e3 }; //---Initial state vector
const double tSpan =  4.306072857605435e+04; // time span of segment integration
const int numSeg =13 ; //15;				// number of segments to integrate
std::string outputPath = "./resultMCPISHMSerialHEO.txt"; // Output File

int  N=169; //---Chebyshev polynomial degree
int M=N+1;  //---Number of nodes

const int Nstates= 6; //----Number of states
const int MaxIter= 100; // Maximum Iteration - when the loop does not converege
double mcpi_tol= 1.0e-17; //Tolerance


//-----------Integration---------------------------//
// main function
int main() {

	cout << "************************************************************************" << endl;
    cout << "*        Copyright (c) Ahmed Atallah &Ahmad Bani Younes 2019           *" << endl;
    cout << "*         San Diego State University - AEROSPACE department            *" << endl;
    cout << "************************************************************************" << endl;
    cout << "    Subject       : MCPI - trajectory propagator for EGM2008            " << endl;
    cout << "    Description   : Integrate perturbed satellite motion with EGM 2008  " << endl;
	cout << "                    Gravity Acceleration in spherical Harmonic.         " << endl;
	cout << "    Date Modified : 08/13/2019                                          " << endl;
	cout << "************************************************************************" << endl << endl << endl;



	// allocating host memory for constant arrays
	int i;
	double x0[Nstates  ] = {0.0};  // starting Cond. Vector
	double Xo[M*Nstates] = {0.0};  // starting guess Matrix
	for ( int kk=0; kk<Nstates; kk++ )
		x0[kk] = init[kk];
		
	for (int m=0; m<M;m++)
		for ( int kk=0; kk<Nstates; kk++ )
			Xo[IDX2F(m+1,kk+1,M)] = init[kk];
	
	// BUILD TAU VECTOR {MCPI printf( "- time domain}
	double TAU[M] = {0.0};
	for ( i=0; i<M; i++ ) 
		TAU[i] = cos((double)i*Pi/(double)(M-1) +Pi);
	// allocating new memory space for loop arrays

	double xAdd[M*Nstates]={0.0};
	double G[M*Nstates]={0.0};
	double Xn[M*Nstates]={0.0};
	
	// timer for loop operations
	double timeSubArr[numSeg] = {0.0};
	double timeAddArr[numSeg] = {0.0};
	double timearray[numSeg]={0.0};
	double Tsum = {0.0};



	// Loading the MCPI coeff's
	// allocating MCPI-loop memory
	printf( "Allocating MCPI-loop memory for matrix calculations...\n\n" );
	double Im[M*M];
    MCPI_CoeffsI(N,M,Im);
	double ImT[M*M];
	trans( Im,ImT,M, M,M);


		
	// OPEN output file to write the results
	FILE *fToOpen;
	if((fToOpen=fopen( outputPath.c_str( ),"w"))==NULL)
	{
		printf("Cannot open file to write result!");
		return EXIT_FAILURE;
	}
	printf( "MCPI- trajectory propagation....\n" );
	printf( "Convergence Report \n" );
	printf( "==========================================================================\n\n" );
	double timeArray[M*numSeg] = {0};
	double timeArraySeg[M]={0};
	for(int IT=0; IT<numSeg; IT++){	
		clock_t mcpistart = clock();       //     get initial time in sec's 
		// set up time variables for current segment
		double b = ( (double)IT + 1.0 )*tSpan/numSeg;
		double a = (double)IT *tSpan/numSeg;
		double timeSub = (b-a)/2.0;
	    double timeAdd = (b+a)/2.0;
		timeSubArr[IT] = timeSub;
		timeAddArr[IT] = timeAdd;

		// iterate until convergence on this segment
	    int loopCount = 0;
	    double temp  = 1.000000000000000;    // this is to find the maximum Error
	   
		for ( int k=0; k<M; k++) {
			timeArraySeg[k] = timeSubArr[IT] * TAU[k] + timeAddArr[IT];
		}
		
		
//	    while (temp>mcpi_tol)
		while (loopCount<MaxIter)
		{ 
			accGravity(M, 200,timeArraySeg,Xo,G);
			matmul(ImT,G,xAdd,M,M,Nstates);			
			errorAndUpdate(M,timeSub,Nstates,x0,Xo,Xn,xAdd,temp);
		
			if(loopCount >= MaxIter)
				break;
			loopCount++;
		} 

		// Update the Initial Conditions
		for ( int kk=0; kk<Nstates; kk++ )
			x0[kk] = Xn[IDX2F(M,kk+1,M)];
		for (int m=0; m<M;m++)
		for ( int kk=0; kk<Nstates; kk++ )
			Xo[IDX2F(m+1,kk+1,M)] = x0[kk];
		timearray[IT] = ((double)clock() - mcpistart)*1000.0 / CLOCKS_PER_SEC;
		printf( " Orbit # %i takes (ms): %lf, converged in %i iterations\n", IT+1, timearray[IT], loopCount );
		Tsum += timearray[IT];	


		for ( int k=0; k<M; k++ ) {
			for ( int kk=0; kk<Nstates; kk++ )
				fprintf( fToOpen, "%6.16lf\t", Xn[IDX2F(k+1,kk+1,M)] );

			fprintf( fToOpen, "%6.16lf\t", timeArraySeg[k] );
			fprintf( fToOpen, "\n");
		}
	} // end for IT loop
	
	//fclose(fToOpen);	

	cout<<"\n\nDONE..."<<"\n\n";
		printf( "==========================================================================\n\n" );
	printf( "Total computation time for integrating the whole trajectory (ms): %lf\n\n", Tsum );	

    
    return 0;
}
