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

   File name     : const.h
   Subject       : Example is used in ....................
   Description   :  Constants

   
   Date Modified : 09/23/2019
     


   File name     : SHM.c
   Subject       : Example is used in ....................
   Description   :  Source Functions to evaluate the acceleratio due to EGM2008 gravity perturbations.

   
   Date Modified : 09/23/2019
     
*/

#include "../include/SHM.h"
#include "EGM2008.cc"
#include "../include/Basic.h"
#include <omp.h>
//#include "const.h"

// Declare Needed Variables
/*!
 * This is the allocation for the maximum size associated
 * Legendre polynomical matrix
 */
double P [(Max_Degree+3)*(Max_Degree+3)];
/*!
 * This is the allocation for the maximum size associated
 * Legendre polynomical matrix scale factor
 */
double scaleFactor [(Max_Degree+3)*(Max_Degree+3)];

/*!
 * \brief Matrix Multiplication
 * This is a simple matrix multiplication function
 *
 * \param[in] A Vector representation of matrix A
 * \param[in] B Vector representation of matrix B
 * \param[in] m Column dimension of A
 * \param[in] n Shared dimension of A and B
 * \param[in] q Row dimension of B
 * \param[out] OUT Matrix Output
 */
void matmulEGM(double* A, double* B, double* OUT, int m, int n, int q)
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

/*!
 * \brief Gravity Evaluation
 * This is the function that evaluates the spherical harmonic
 * series to provide acceleration
 *
 * \param[in] p 3 element position in ECEF
 * \param[in] Deg Degree and order of the serries to be used
 * \param[out] Gxyz Gravitational Acceloration Output
 */
void EGM2008( double* p, double* Gxyz, int DEG)
{

    double r             = {0.0};
    double phic          = {0.0};
    double lambda        = {0.0};
    double slambda       = {0.0};
    double clambda       = {0.0};
    double x = {0.0};
    double y = {0.0};
    double z = {0.0};
    double smlambda[Max_Degree+1] = {0.0};
    double cmlambda[Max_Degree+1] = {0.0};

    // determine radius of this thread's node
    x = p[0];
    y = p[1];
    z = p[2];

    // Compute geocentric radius
    r = pow( x*x + y*y + z*z , 0.5 );
    // Compute geocentric latitude
    phic  = asin( z / r );
    // Compute lambda
    lambda  = atan2( y, x );
    while (lambda<0)
        lambda = lambda+2*C_PI;
    while (lambda>=2*C_PI)
        lambda = lambda-2*C_PI;

//	printf("  Lambda  %20.19f\n", lambda);
    slambda = sin(lambda);
    clambda = cos(lambda);
    smlambda[0] = 0.0;
    cmlambda[0] = 1.0;
//printf( "0  Cmlambda Smlambda %d %20.19f %20.19f\n\n",0,cmlambda[0],smlambda[0]);
    smlambda[1] = slambda;
    cmlambda[1] = clambda;
//	printf( "1  Cmlambda Smlambda %d %20.19f %20.19f\n\n",1,cmlambda[1],smlambda[1]);
    for(int m=2;m<DEG+1;m++){
        smlambda[m] = 2.0*clambda*smlambda[m-1] - smlambda[m-2];
        cmlambda[m] = 2.0*clambda*cmlambda[m-1] - cmlambda[m-2];
//	 printf( "m  Cmlambda Smlambda %d %20.19f %20.19f\n\n",m,cmlambda[m],smlambda[m]);
    }


    loc_gravLegendre( phic, scaleFactor, P , DEG);

    loc_gravityPCPF( p, P, DEG, smlambda, cmlambda, r, scaleFactor, Gxyz );


}

/*!
 * \brief Gravity Potential Evaluation
 * This is the function that evaluates the spherical harmonic
 * serries to provide gravitational potential
 *
 * \param[in] p 3 element position in ECEF
 * \param[in] Deg Degree and order of the serries to be used
 * \param[out] Pot Gravitational Potential Output
 */
void EGM2008Pot( double* p, double* Pot, int DEG)
{
	// determine radius of this thread's node
	double r             = {0.0};
	double phic          = {0.0};
	double lambda        = {0.0};
	double slambda       = {0.0};
	double clambda       = {0.0};
		double smlambda[Max_Degree+1]  ;
		double cmlambda[Max_Degree+1] ;

		double x = p[0];
		double y = p[1];
		double z = p[2];
		int m;
		double P [(Max_Degree+3)*(Max_Degree+3)] ;
		double scaleFactor [(Max_Degree+3)*(Max_Degree+3)] ;

		// Compute geocentric radius
		r = pow( x*x + y*y + z*z , 0.5 );
		// Compute geocentric latitude
		phic  = asin( z / r );
		// Compute lambda
		lambda  = atan2( y, x );
		while (lambda<0){
				lambda = lambda+2*C_PI;
		}
		while (lambda>=2*C_PI){
				lambda = lambda-2*C_PI;
		}

		slambda = sin(lambda);
		clambda = cos(lambda);
		smlambda[0] = 0.0;
		cmlambda[0] = 1.0;
		smlambda[1] = slambda;
		cmlambda[1] = clambda;


		for(m=2;m<DEG+1;m++){
				smlambda[m] = 2.0*clambda*smlambda[m-1] - smlambda[m-2];
				cmlambda[m] = 2.0*clambda*cmlambda[m-1] - cmlambda[m-2];
//	 printf( "m  Cmlambda Smlambda %d %20.19f %20.19f\n\n",m,cmlambda[m],smlambda[m]);			
	}
		// Compute normalized associated legendre polynomials

		loc_gravLegendre( phic, scaleFactor, P , DEG);

		loc_gravityPot( p, P, DEG, smlambda, cmlambda, r, scaleFactor, Pot );


}

/*!
 * \brief Legendre Polyniomial Evaluation
 * This is the function that computes the normalized associated
 * legendre polynomials based on geocentric latitude
 *
 * \param[in] phi Geocentric latitude
 * \param[in] Deg Degree and order of the serries to be used
 * \param[out] P associated Legendre polynomial matrix
 * \param[out] scaleFactor Legendre scale factor
 */
void loc_gravLegendre( double phi, double* scaleFactor, double* P, int DEG )
{

	int k, p;
	double cphi = {0.0};
	double sphi = {0.0};


    cphi = cos(0.5*C_PI - phi);
    sphi = sin(0.5*C_PI - phi);
    // Seeds for recursion formula
    P[IDX2F(1,1,Max_Degree+3)] = 1.0;            // n = 0, m = 0;
    scaleFactor[IDX2F(1,1, Max_Degree+3)] = 0.0;
    P[IDX2F(2,1, Max_Degree+3)] = sqrt(3.0)*cphi ; // n = 1, m = 0;
    scaleFactor[IDX2F(2,1,Max_Degree+3)]  = 1.0;
    P[IDX2F(2,2,Max_Degree+3)] = sqrt(3.0)*sphi; // n = 1, m = 1;
    scaleFactor[IDX2F(2,2,Max_Degree+3)] = 0.0;

// // New Method
//     int nn = 2;
//     int mm = 0;
//     double m, n;
//     int limit = (DEG+3)*(DEG+4)/2;
//     for (int counter = 3; counter <= limit;counter++){
//         k = nn + 1;
//         p = mm + 1;
//         n = coefMatrix[counter][0];
//         m = coefMatrix[counter][1];
//
//         switch ( (int)(coefMatrix[counter][4]) )
//        {
//           case 1:
//             P[IDX2F(k,k,Max_Degree+3)] = sqrt(2*n+1.0)/sqrt(2.0*n)*sphi*P[IDX2F(k-1,k-1, Max_Degree+3)];
//             scaleFactor[IDX2F(k,k,Max_Degree+3)] = 0.0;
//              nn++;
//              mm = 0;
//              break;
//           case 2:
//             P[IDX2F(k,p,Max_Degree+3)] = (sqrt(2*n+1)/n)*(sqrt(2*n-1.0)*cphi*P[IDX2F(k-1,p,Max_Degree+3)] - (n-1)/sqrt(2*n-3)* P[IDX2F(k-2,p, Max_Degree+3)] );
//             scaleFactor[IDX2F(k,p,Max_Degree+3)] = sqrt( (n+1)*(n)/2);
//             mm++;
//             break;
//           case 3:
//             P[IDX2F(k,p,Max_Degree+3)] = sqrt(2*n+1)/(sqrt(n+m)*sqrt(n-m))*(sqrt(2*n-1.0)*cphi*P[IDX2F(k-1,p,Max_Degree+3)] - sqrt(n+m-1.0)*sqrt(n-m-1)/sqrt(2*n-3)*P[IDX2F(k-2,p,Max_Degree+3)] );
//             scaleFactor[IDX2F(k,p,Max_Degree+3)] = sqrt( (n+m+1)*(n-m));
//             mm++;
//             break;
//
//        }
//
//     }

// Old Method
   for (int nn = 2; nn <= DEG+2;nn++){
       double n = (double)nn;
       k = nn + 1;
       for(int mm=0; mm<=n;mm++) {
           double m = (double)mm;
           p = mm + 1;
           // Compute normalized associated legendre polynomials, P, via recursion relations
           // Scale Factor needed for normalization of dUdphi partial derivative
           if (n == m){
               P[IDX2F(k,k,Max_Degree+3)] = sqrt(2*n+1.0)/sqrt(2.0*n)*sphi*P[IDX2F(k-1,k-1, Max_Degree+3)];
               scaleFactor[IDX2F(k,k,Max_Degree+3)] = 0.0;
           }
           else if (m == 0){
               P[IDX2F(k,p,Max_Degree+3)] = (sqrt(2*n+1)/n)*(sqrt(2*n-1.0)*cphi*P[IDX2F(k-1,p,Max_Degree+3)] - (n-1)/sqrt(2*n-3)* P[IDX2F(k-2,p, Max_Degree+3)] );
               scaleFactor[IDX2F(k,p,Max_Degree+3)] = sqrt( (n+1)*(n)/2);
           }
           else {
               P[IDX2F(k,p,Max_Degree+3)] = sqrt(2*n+1)/(sqrt(n+m)*sqrt(n-m))*(sqrt(2*n-1.0)*cphi*P[IDX2F(k-1,p,Max_Degree+3)] - sqrt(n+m-1.0)*sqrt(n-m-1)/sqrt(2*n-3)*P[IDX2F(k-2,p,Max_Degree+3)] );
               scaleFactor[IDX2F(k,p,Max_Degree+3)] = sqrt( (n+m+1)*(n-m));
           }
       }
   }
}

/*!
 * \brief Internal Gravitational Acceloration Evaluation
 * This is the function that computes the gravitational acceloration based on
 * the associated Legendre polynomials and the state
 *
 * \param[in] p Position vector in ECEF
 * \param[in] P associated Legendre polynomial matrix
 * \param[in] scaleFactor Legendre scale factor
 * \param[in] Deg Degree and order of the serries to be used
 * \param[in] r Position vector in ECEF
 * \param[in] smlambda Trigonometric function of longitude
 * \param[in] smlambda Trigonometric function of longitude
 * \param[out] Gxyz Gravitational Acceloration Output
 */
void loc_gravityPCPF( double* p, double* P, int DEG, double* smlambda, double* cmlambda, double r, double* scaleFactor, double* Gxyz )
{

	int k, j;
	double radu;
	double rRatio  = {0.0};
	double rRatio_n  = {0.0};

	// initialize summation of gravity in radial coordinates
	double dUdrSumN   = {1.0};
	double dUdphiSumN    = {0.0};
	double dUdlambdaSumN = {0.0};

	double dUdrSumM   = {0.0};
	double dUdphiSumM = {0.0};
	double dUdlambdaSumM = {0.0};

	double dUdr    = {0.0};
	double dUdphi  = {0.0};
	double dUdlambda = {0.0};



    double x = p[0];
    double y = p[1];
    double z = p[2];
    radu = r;
    rRatio = C_Req/radu;
    rRatio_n = rRatio;
    // summation of gravity in radial coordinates




//    // New Method
//        int nn = 2;
//        int mm = 0;
//        double m, n;
//        dUdrSumM      = 0.0;
//        dUdphiSumM    = 0.0;
//        dUdlambdaSumM = 0.0;
//        int limit = (DEG+1)*(DEG+2)/2;
//        for (int counter = 3; counter <= limit;counter++){
//            k = nn + 1;
//            j = mm + 1;
//            n = coefMatrix[counter][0];
//            m = coefMatrix[counter][1];

//            switch ( (int)(coefMatrix[counter][4]) )
//           {
//              case 1:
//                dUdrSumM      = dUdrSumM + P[IDX2F(k,j,Max_Degree+3)] *(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
//                dUdphiSumM    = dUdphiSumM + ((P[IDX2F(k,j+1,Max_Degree+3)]*scaleFactor[IDX2F(k,j,Max_Degree+3)]) - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree+3)])*(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
//                dUdlambdaSumM = dUdlambdaSumM + m*P[IDX2F(k,j,Max_Degree+3)]*(S[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] - C[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
//                rRatio_n = rRatio_n*rRatio;
//                dUdrSumN      = dUdrSumN      + dUdrSumM*rRatio_n*k;
//                dUdphiSumN    = dUdphiSumN    + dUdphiSumM*rRatio_n;
//                dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM*rRatio_n;
//                 nn++;
//                 mm = 0;
//                 dUdrSumM      = 0.0;
//                 dUdphiSumM    = 0.0;
//                 dUdlambdaSumM = 0.0;
//                 break;
//              default:
//                dUdrSumM      = dUdrSumM + P[IDX2F(k,j,Max_Degree+3)] *(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
//                dUdphiSumM    = dUdphiSumM + ((P[IDX2F(k,j+1,Max_Degree+3)]*scaleFactor[IDX2F(k,j,Max_Degree+3)]) - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree+3)])*(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
//                dUdlambdaSumM = dUdlambdaSumM + m*P[IDX2F(k,j,Max_Degree+3)]*(S[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] - C[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
//                mm++;

//           }

//        }

    // // Old Method
    for (int n = 2; n <= DEG; n++) {
        k = n+1;
        rRatio_n = rRatio_n*rRatio;
        dUdrSumM      = 0.0;
        dUdphiSumM    = 0.0;
        dUdlambdaSumM = 0.0;
        for (int m = 0; m <= n; m++){
            j = m+1;
//		printf(  " dUdrSumM  %f\n",dUdrSumM);
            dUdrSumM      = dUdrSumM + P[IDX2F(k,j,Max_Degree+3)] *(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);

            dUdphiSumM    = dUdphiSumM + ((P[IDX2F(k,j+1,Max_Degree+3)]*scaleFactor[IDX2F(k,j,Max_Degree+3)]) -
 z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree+3)])*(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);

            dUdlambdaSumM = dUdlambdaSumM + m*P[IDX2F(k,j,Max_Degree+3)]*(S[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] - C[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);

        }
        dUdrSumN      = dUdrSumN      + dUdrSumM*rRatio_n*k;

        dUdphiSumN    = dUdphiSumN    + dUdphiSumM*rRatio_n;
 dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM*rRatio_n;

    }

    // gravity in spherical coordinates
    dUdr      = -C_MU/(radu*radu)*dUdrSumN ;
    dUdphi    =  C_MU/radu*dUdphiSumN ;
    dUdlambda =  C_MU/radu*dUdlambdaSumN ;

    //gravity in ECEF coordinates
    Gxyz[0] = ((1.0/radu)*dUdr - (z/(radu*radu*sqrt(x*x + y*y)))*dUdphi)*x - (dUdlambda/(x*x + y*y))*y;
    Gxyz[1] = ((1.0/radu)*dUdr - (z/(radu*radu*sqrt(x*x + y*y)))*dUdphi)*y + (dUdlambda/(x*x + y*y))*x;
    Gxyz[2] = (1.0/radu)*dUdr*z + ((sqrt(x*x + y*y))/(radu*radu))*dUdphi;

    // special case for poles
    /*
     atPole = abs(atan2(p[IDX2F(i+1,3,M)],sqrt(p[IDX2F(i+1,1,M)]^2 + p[IDX2F(i+1,2,M)]^2)))==C_PI/2;
     if any(atPole){
     gx(atPole) = 0;
     gy(atPole) = 0;
     gz(atPole) = (1./r(atPole)).*dUdr(atPole).*p((atPole),3);
     }
     */

}

double round1(double var) 
{ 
    // 37.66666 * 100 =3766.66 
    // 3766.66 + .5 =3767.16    for rounding off value 
    // then type cast to int so value is 3767 
    // then divided by 100 so the value converted into 37.67 
    double value =  (int)(10000000000*var); 
    return (double)value / 10000000000; 
} 
  
void accGravity(int M, int Deg,double *t,double* Xo,double* hG)
{
	double x0[3];
	double Gxyz[3];

//	#pragma omp parallel for			
//	#pragma omp parallel for
	for (int i=0; i<M; i++)
	{
	//	double th_interp =7292115.0e-011*t[i];
		
		double th_interp =0;
			
//		x0[0]=cos(th_interp)*Xo[i+0*M]+sin(th_interp)*Xo[i+1*M];
//		x0[1]=-sin(th_interp)*Xo[i+0*M]+cos(th_interp)*Xo[i+1*M];
//		x0[2]=Xo[i+2*M];
		x0[0]=Xo[i+0*M];
	        x0[1]=Xo[i+1*M];
                x0[2]=Xo[i+2*M];

	//	printf(" x0  %f \n",x0[0]);
	//	printf(" x1  %f \n",x0[1]);
		EGM2008( x0,  Gxyz, Deg);

		hG[i+3*M]=Gxyz[0];
		hG[i+4*M]=Gxyz[1];
		hG[i+5*M]=Gxyz[2];
//		Gxyz[0]=Round_off(Gxyz[0],10);
//		Gxyz[1]=Round_off(Gxyz[1],10);
//		Gxyz[2]=Round_off(Gxyz[2],10);

	//	hG[i+3*M]=cos(th_interp)*(Gxyz[0])+sin(th_interp)*(Gxyz[1]);
//		hG[i+4*M]=-sin(th_interp)*(Gxyz[0])+cos(th_interp)*(Gxyz[1]);
//		hG[i+5*M]=(Gxyz[2]);

		hG[i+0*M]=Xo[i+3*M];//+0*Xo[i+1*M];
		hG[i+1*M]=Xo[i+4*M];//-0*Xo[i+0*M];
		hG[i+2*M]=Xo[i+5*M];
			
	}
}


/*!
 * \brief Internal Gravity Potential Evaluation
 * This is the function that computes the gravitational acceloration based on
 * the associated Legendre polynomials and the state
 *
 * \param[in] p Position vector in ECEF
 * \param[in] P associated Legendre polynomial matrix
 * \param[in] scaleFactor Legendre scale factor
 * \param[in] Deg Degree and order of the serries to be used
 * \param[in] r Position vector in ECEF
 * \param[in] smlambda Trigonometric function of longitude
 * \param[in] smlambda Trigonometric function of longitude
 * \param[out] Pot Gravitational Potential Output
 */
void loc_gravityPot( double* p, double* P, int DEG, double* smlambda, double* cmlambda, double r, double* scaleFactor, double* Pot )
{

	int k, j;
	double radu;
	double rRatio  = {0.0};
	double rRatio_n  = {0.0};

	// initialize summation of gravity in radial coordinates
	double dUdrSumN   = {1.0};
	double dUdphiSumN    = {0.0};
	double dUdlambdaSumN = {0.0};

	double dUdrSumM   = {0.0};
	double dUdphiSumM = {0.0};
	double dUdlambdaSumM = {0.0};

	double dUdr  = {0.0};
	double dUdphi  = {0.0};
	double dUdlambda = {0.0};
	double USum = {0.0};


		double x = p[0];
		double y = p[1];
		double z = p[2];
		int n;
		int m;

		radu = r;
		rRatio = C_Req/radu;
		rRatio_n = rRatio;
		// summation of gravity in radial coordinates

		for (n = 2; n <= DEG; n++) {
				k = n+1;
				rRatio_n = rRatio_n*rRatio;
				dUdrSumM      = 0.0;
				dUdphiSumM    = 0.0;
				dUdlambdaSumM = 0.0;
				for (m = 0; m <= n; m++){
						j = m+1;

						dUdrSumM      = dUdrSumM + P[IDX2F(k,j,Max_Degree+3)] *(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
					
						dUdphiSumM    = dUdphiSumM + ((P[IDX2F(k,j+1,Max_Degree+3)]*
scaleFactor[IDX2F(k,j,Max_Degree+3)]) - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree+3)])*
(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
						dUdlambdaSumM = dUdlambdaSumM + m*P[IDX2F(k,j,Max_Degree+3)]*(S[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] - C[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
						// printf("UsumM: %e\n ",dUdrSumM);
				}
				dUdrSumN      = dUdrSumN      + dUdrSumM*rRatio_n;
				dUdphiSumN    = dUdphiSumN    + dUdphiSumM*rRatio_n;
				dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM*rRatio_n;
		}
		// gravity in spherical coordinates
		dUdr      =  C_MU/(radu)*dUdrSumN ;
		dUdphi    =  C_MU/radu*dUdphiSumN ;
		dUdlambda =  C_MU/radu*dUdlambdaSumN ;

		// gravity in potential
		// *Pot = MU * sqrt(dUdr*dUdr + dUdphi*dUdphi + dUdlambda*dUdlambda)/radu;
		*Pot =  sqrt(dUdr*dUdr);


}

/*!
 * \brief Jacobi Integral
 * This is the function computes the Jacobi Integral based on position and
 * state vector
 *
 * This is a usefule tool for checking the accuracy of conservitive
 * orbit propigation
 *
 * \param[in] solN State (position and velocity) vector in ECEF
 * \param[in] Deg Degree and order of the serries to be used
 * \param[out] H Jacobi Integral Output
 */
void jacobiIntegral(double* solN, double* H, int Deg){

	double KE,PE,RotTerm;

		KE = 0.5*(solN[3]*solN[3] + solN[4]*solN[4] + solN[5]*solN[5]);
		EGM2008Pot(solN, &PE, Deg);
		PE = -PE;
		RotTerm = 0.5*C_omega*C_omega*(solN[0]*solN[0] + solN[1]*solN[1]);
		*H  = PE + KE - RotTerm; // Hamiltonian

		// printf("KE: %e\tPE: %e\tRT: %e\tSum: %e\n ",KE,PE,RotTerm,*H);
		// getchar();

}
