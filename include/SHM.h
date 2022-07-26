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
     


   File name     : SHM.h
   Subject       : Example is used in ....................
   Description   :  Header for Source Functions to evaluate the acceleratio due to EGM2008 gravity perturbations.

   
   Date Modified : 09/23/2019
     
*/
#ifndef __SHM__
#define __SHM__





#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "const.h"

/* Macro for looking up indices using "1" based addressing ( a la fortran and matlab )
 and column-major array format [Note ld = Num of rows]
 */


#define Max_Degree 210
//#define N  200			// number of nodes 
//const int M=201;			// number of nodes + 1

/**
  Round a / b to nearest higher integer value
 
inline int iDivUp(int a, int b)
{
    return (a % b != 0) ? (a / b + 1) : (a / b);
}
*/

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
void matmulEGM(double* A, double* B, double* OUT, int m, int n, int q);


/*!
 * \brief Gravity Evaluation
 * This is the function that evaluates the spherical harmonic
 * series to provide acceloration
 *
 * \param[in] p 3-element position in ECEF
 * \param[in] Deg Degree and order of the series to be used
 * \param[in] delV_del_r_phi_lambda 9-element vector of gravity partials used to compute Jacobian fo$
 * \param[out] Gxyz Gravitational Acceloration Output
 */
void EGM2008( double* p, double* Gxyz, int DEG);


/*!
 * \brief Gravity Potential Evaluation
 * This is the function that evaluates the spherical harmonic
 * series to provide gravitational potential
 *
 * \param[in] p 3 element position in ECEF
 * \param[in] Deg Degree and order of the series to be used
 * \param[out] Pot Gravitational Potential Output
 */
void EGM2008Pot( double* p, double* Pot, int DEG);


/*!
 * \brief Legendre Polyniomial Evaluation
 * This is the function that computes the normalized associated
 * legendre polynomials based on geocentric latitude
 *
 * \param[in] phi Geocentric latitude
 * \param[in] Deg Degree and order of the series to be used
 * \param[out] P associated Legendre polynomial matrix
 * \param[out] scaleFactor Legendre scale factor
 */
void loc_gravLegendre( double phi, double* scaleFactor, double* P, int DEG );


/*!
 * \brief Internal Gravitational Acceloration Evaluation
 * This is the function that computes the gravitational acceloration based on
 * the associated Legendre polynomials and the state
 *
 * \param[in] p Position vector in ECEF
 * \param[in] P associated Legendre polynomial matrix
 * \param[in] scaleFactor Legendre scale factor
 * \param[in] Deg Degree and order of the series to be used
 * \param[in] r Position vector in ECEF
 * \param[in] smlambda Trigonometric function of longitude
 * \param[in] smlambda Trigonometric function of longitude
 * \param[out] Gxyz Gravitational Acceloration Output
 */
void loc_gravityPCPF( double* p, double* P, int DEG , double* smlambda, double* cmlambda, double r, double *scaleFactor, double* Gxyz );


/*!
 * \brief Internal Gravity Potential Evaluation
 * This is the function that computes the gravitational acceloration based on
 * the associated Legendre polynomials and the state
 *
 * \param[in] p Position vector in ECEF
 * \param[in] P associated Legendre polynomial matrix
 * \param[in] scaleFactor Legendre scale factor
 * \param[in] Deg Degree and order of the series to be used
 * \param[in] r Position vector in ECEF
 * \param[in] smlambda Trigonometric function of longitude
 * \param[in] smlambda Trigonometric function of longitude
 * \param[out] Pot Gravitational Potential Output
 */
void loc_gravityPot( double* p, double* P, int DEG , double* smlambda, double* cmlambda, double r, double *scaleFactor, double* Pot );


/*!
 * \brief Jacobi Integral
 * This is the function computes the Jacobi Integral based on position and
 * state vector
 *
 * This is a useful tool for checking the accuracy of conservitive
 * orbit propagation
 *
 * \param[in] solN State (position and velocity) vector in ECEF
 * \param[in] Deg Degree and order of the series to be used
 * \param[out] H Jacobi Integral Output
 */
void jacobiIntegral(double* solN, double* H, int Deg);

void accGravity(int M, int Deg,double *t,double* Xo,double* hG);


#endif
