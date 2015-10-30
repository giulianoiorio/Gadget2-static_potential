#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


double KeplerPotential(double x, double y, double z, double mass)
{
	double R= sqrt(x*x+y*y+z*z);
	return -mass/R;
}

double *KeplerForce(double x, double y, double z, double mass)
{
	static double F[3];
	double R2= x*x+y*y+z*z;
	double den=(pow(R2,1.5));
	
	F[0]= -(mass*x)/den;
	F[1]= -(mass*y)/den;
	F[2]= -(mass*z)/den;
	
	//printf("mass= %f, x= %f %f %f, den= %f, force=%f", mass,x,y,z,den,F[0]);
	
	return F;
}

double MyamotoNagaiPotential(double x, double y, double z, double a, double b, double mass)
{
	double R2= x*x+y*y;
	double z_cyl= (a*a + sqrt(b*b+z*z));
	
	return -mass/sqrt( R2 + z_cyl*z_cyl );
}

double *MyamotoNagaiForce(double x, double y, double z, double a, double b, double mass)
{
	static double F[3];
	double R2= x*x+y*y;
	double z_cyl= (a*a + sqrt(b*b+z*z));
	double common_den= pow( (R2 + z_cyl*z_cyl), 1.5 );
	
	
	F[0]=-(x*mass)/(common_den);
	F[1]=-(y*mass)/(common_den);
	F[2]=-(z*mass*(z_cyl))/(common_den * sqrt(b*b+z*z) );
	
	return F;
}

double LogartmicHaloPotential(double x, double y, double z, double v_h2, double R_c, double q)
{
	double R2= x*x+y*y;
	double z_cyl=z/q; 
	
	return  ( 0.5*v_h2*log(R2+ z_cyl*z_cyl + R_c*R_c) )/(All.G);
}

double *LogartmicHaloForce(double x, double y, double z, double v_h2, double R_c, double q)
{
	static double F[3];
	double R2= x*x+y*y;
	double z_cyl=z/q; 
	double den=R2 + z_cyl*z_cyl + R_c*R_c;
	
	F[0]=-(x*v_h2)/(den*All.G);
	F[1]=-(y*v_h2)/(den*All.G);
	F[2]=-(z*v_h2)/(den*q*q*All.G);
	
	return F;
}

double HernquistPotential(double x, double y, double z, double mass, double rc)
{	
	double r= sqrt(x*x+y*y+z*z);
	
	return - mass/(r + rc);
}

double *HernquistForce(double x, double y, double z, double mass, double rc)
{	
	double r= sqrt(x*x+y*y+z*z);
	static double F[3];
	double cost= mass/(r*(rc+r)*(rc+r));
	//double cost=0;
	
	if (r>0.0001) //To avoid divergence
	{
	F[0]=-x*cost;
	F[1]=-y*cost;
	F[2]=-z*cost;
	}
	else 
	{
	F[0]=0;
	F[1]=0;
	F[2]=0;	
	}
	
	
	return F;
}


#ifdef EXTPOTENTIAL
	double *Static_Force(double x, double y, double z)
	{	
	
		
		switch (EXTPOTENTIAL) 
		{
			case 0:
			{
				const double pointmass=1e11;
				double *p;
				p=KeplerForce(x, y, z, pointmass);
				
				return p;	
				break;
			}
			
			//P07
			case 1:
			{	
				//disc
				const double Mdisc=1e11;
				const double a=6.5;
				const double b=0.26;
				//bulge
				const double Mbulge=3.4e10;
				const double rc=0.7;
				//halo
				const double vh=128*128*2;
				const double rh=12;
				const double q=1;
				
				double *p;
				static double F[3];
				
                
				p=LogartmicHaloForce(x, y, z, vh, rh, q);
				F[0]= *(p);
				F[1]= *(p+1);
				F[2]= *(p+2);
				
               
                
				p=HernquistForce(x, y, z, Mbulge, rc);
				F[0]+= *(p);
				F[1]+= *(p+1);
				F[2]+= *(p+2);
                
                // printf('FOOOOrza - 0: %f %f',x,F[0]);
				
                
				p=MyamotoNagaiForce(x, y, z, a, b, Mdisc);
				F[0]+= *(p);
				F[1]+= *(p+1);
				F[2]+= *(p+2);
                
               // printf('FOOOOrza - 0: %f %f',x,F[0]);
				
                 
				return F;
				break;
				
			}
			
			case 2:
			{
				const double discmass=10.0;
				const double a=1.0;
				const double b=1.0;
				double *p;
				p=MyamotoNagaiForce(x, y, z, a, b, discmass);
				return p;
				break;
			}
			
			case 3:
			{
				const double vh=10.0;
				const double rc=2.0;
				const double q=1.0;
				double *p;
				p=LogartmicHaloForce(x, y, z, vh, rc, q);
				return p;
				break;
			}
			
		}
	}
	
	double Static_Potential(double x, double y, double z)
	{
			switch (EXTPOTENTIAL) 
		{
			case 0:
			{
				const double pointmass=1e11;
				return KeplerPotential(x, y, z, pointmass);
				break;
			}
			
			//P07
			case 1:
			{
				//disc
				const double Mdisc=1e11;
				const double a=6.5;
				const double b=0.26;
				//bulge
				const double Mbulge=3.4e10;
				const double rc=0.7;
				//halo
				const double vh=128*2;
				const double rh=12;
				const double q=1;
				

				
				
				return LogartmicHaloPotential(x, y, z, vh, rh, q) +  HernquistPotential(x, y, z, Mbulge, rc) + MyamotoNagaiPotential(x, y, z, a, b, Mdisc);
				break;
			}
		
			
		}
	}



#else
double *Static_Force(double x, double y, double z)
{	
	static double p[3]= {0.,0.,0.};
	return p;
}

double Static_Potential(double x, double y, double z)
{	
	return 0;
}
#endif







