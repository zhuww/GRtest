/*******************************************************************************
*                                                                              
*  calcGalPdotRange.cc                                                              
*                                                                              
*  Calculates Pdot/P caused by differentia Galactic acceleration               
*  for a range of distances in direction (l,b)  
*                               
*  Uses potential PJM16_best.Tpot as default
*                                                                              
*  Based on C++ code written by Walter Dehnen, 1995-96,                        
*                      Paul McMillan, 2007-,                                   
*  Lund Observatory, Lund University.                                          
*  address: Box 43, SE-221 00 Lund, Sweden                                     
*  e-mail:  paul@astro.lu.se                                                   
*                                                                              
*******************************************************************************/

#include <fstream>
#include <iostream>
#include "GalPot.h"

using std::cout;
using std::cin;


int main(int argc,char *argv[])  
{
    ifstream file;
    
    double DEG2RAD = 0.017453292519943295;
    
    double Rsun=8.2,zsun=0.014;
    double xsun=-Rsun,ysun=0.0;
    
    double l,b,d,d1,d2,step;
    double nx,ny,nz;
    double R,x,y,z;
    double P,dPdR,dPdz;
    double gRsun,gxsun,gysun,gzsun,gR,gx,gy,gz;
    double gdiff,dDoD,dPoP;
    
    if (argc != 5 && argc != 6) {
        cerr << "\n USAGE: calcGalPdot.exe <l[deg]> <b[deg]> <d1[kpc]> <d2[kpc]> (<pot_file)\n\n";
        exit(1);
    }
    
    if (argc == 5) {
            //file.open("/Users/wex/Science/psrsoft/GalPot-McMillan2016/pot/PJM16_best.Tpot");
        //file.open("/home/zhuww/work/timing/1713all/GRtest/GalPotMcMillan2016/pot/PJM16_best.Tpot");
        file.open("/homes/zhuww/data/timing/GRtest/Lazaridis/GalPotMcMillan2016/pot/PJM16_best.Tpot");
    } else {
    	file.open(argv[5]);
    }
    GalaxyPotential Phi(file);
    file.close();

	//l  = std::stod(argv[1]);
	//b  = std::stod(argv[2]);
	//d1 = std::stod(argv[3]);
	//d2 = std::stod(argv[4]);
	l  = atof(argv[1]);
	b  = atof(argv[2]);
	d1 = atof(argv[3]);
	d2 = atof(argv[4]);
	
	//unit vector (towards pulsar) 
	nx = cos(b*DEG2RAD) * cos(l*DEG2RAD);
	ny = cos(b*DEG2RAD) * sin(l*DEG2RAD);
	nz = sin(b*DEG2RAD);
	
    /*
	 * Galactic acceleration of the Sun
     */
     
    P = Phi(Rsun,zsun,dPdR,dPdz); //[kpc, kpc, kpc/Myr^2, kpc/Myr^2]
    
    gRsun = -dPdR * 3.09843657844953E-6; //[cm/s^2]
    gzsun = -dPdz * 3.09843657844953E-6; //[cm/s^2]
    
    gxsun = -gRsun;
    gysun =  0.0;
    
    /*
     * position of the pulsar in the Galaxy 
     */
    
    step = (d2 - d1)/100.0;
    if(step > 0.01) step = 0.01; //maximum step size of 10 pc
    
    for(d = d1; d < d2 + 1.0E-12; d += step) 
    {
		x = xsun + d*nx;
		y = ysun + d*ny;
		z = zsun + d*nz;

		R = sqrt(x*x + y*y);
	 
		//*** Galactic acceleration of the Pulsar

		P = Phi(R,z,dPdR,dPdz); //[kpc, kpc, kpc/Myr^2, kpc/Myr^2]

		gR = -dPdR * 3.09843657844953E-6; //[cm/s^2]
		gz = -dPdz * 3.09843657844953E-6; //[cm/s^2]

		gx = gR*x/R;
		gy = gR*y/R;

		 //*** Differential acceleration
 
		gdiff = (gx - gxsun)*nx + (gy - gysun)*ny + (gz - gzsun)*nz;
		dDoD  = -gdiff / 2.99792458E10;
		dPoP  = -dDoD;
 
		cout << d << "\t" << dPoP << "\n";
    }
}
