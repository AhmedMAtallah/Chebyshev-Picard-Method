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



   File name     : MCPIOM.h
   Subject       : Example is used in ....................
   Description   :  Header file for Source Functions of MCPI   

  Date Modified : 09/23/2019  
     
*/

//#include "MCPIIOM.h"
#include "const.h"
#include "Basic.h"

void MCPI_CoeffsI(int N, int M,double*Im);

void errorAndUpdate(int MM,double timeSub,int Nstates2,double * x0,double *Xo,double *Xn,double* xAdd,double temp);
