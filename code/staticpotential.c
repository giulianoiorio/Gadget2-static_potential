#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


double KeplerPotential(double x, double y, double z, double mass)
{
	double R= sqrt(x*x+y*y+z*z);
	return -mass/R;
}

double *KeplerForce(double x, double y, double z, double mass)
{
	static double F[3];
	double R= x*x+y*y+z*z;
	double den=(pow(R,1.5));
	
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

double LogartmicHaloPotential(double x, double y, double z, double v_h, double R_c, double q)
{
	double R2= x*x+y*y;
	double z_cyl=z/q; 
	
	return -0.5*v_h*log(R2+ z_cyl*z_cyl + R_c*R_c);
}

double *LogartmicHaloForce(double x, double y, double z, double v_h, double R_c, double q)
{
	static double F[3];
	double R2= x*x+y*y;
	double z_cyl=z/q; 
	double den=R2 + z_cyl*z_cyl + R_c*R_c;
	
	F[0]=-(x*v_h)/den;
	F[1]=-(y*v_h)/den;
	F[2]=-(z*v_h)/(den*q*q);
	
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
		
			case 1:
			{
				const double discmass=10.0;
				const double a=1.0;
				const double b=1.0;
				double *p;
				p=MyamotoNagaiForce(x, y, z, a, b, discmass);
				return p;
				break;
			}
			
			case 2:
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







