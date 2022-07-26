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
     
*/

#ifndef _CONSTANTS_
#define _CONSTANTS_

#define C_PI 3.1415926535897932      // Pi
#define C_MU 3.986004418e14           // Gravitational Constant [km^3/s^2]
#define C_MUCan 1                    // Gravitational Constant Canonical Units
#define C_omega 7292115.0e-011       // Angular Speed of Earth [rad/s]
#define C_Req 6378137               // Equatorial Radius of Earth [km]
#define DU C_Req
#define TU sqrt(pow(DU,3)/C_MU)
#define pi 3.1415926535897932

#define Nx 16                     // Number of CGL nodes
#define Nxp1 17
#define Ny   16                     // Number of CGL nodes
#define Nyp1 17                 // Number of CGL node   s

#define Nx2   13                     // Number of CGL nodes
#define Nx2p1 14
#define Ny2   13                     // Number of CGL nodes
#define Ny2p1 14                 // Number of CGL node   s

#define Re  6378137.0
#define dr  0.02*Re
#define rmin Re
#define Nr 9 
#define Nrp1  10 

#define dr2  5.98*Re
#define rmin2 1.02*Re
#define Nr2 37 
#define Nr2p1  38 

#define Pi  3.1415926535897932
#define MU_EGM  3.986004418000000e+014
#define REQ 6378137

#define GM  3.986004418000000e+014
#define J2 1082.63E-6   // % Perturbation Parameters J2
#define rmax  7*Re
//	double rmin    = 1.02*Re;
#define rmid   1.02*Re
#endif
