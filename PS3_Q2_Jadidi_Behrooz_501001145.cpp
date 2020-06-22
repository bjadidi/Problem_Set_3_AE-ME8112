#include<iostream>
#include <math.h>
#include<cmath>
#include <fstream>
using namespace std;

/*////////////////Begin_Header///////////////////////////////////////////////////
================================================================================
// This C++ code was written by Behrooz Jadidi_Student ID:501001145
// The first version was generated on MM_DD_YYYY
// This code was written as a solution to AE/ME8112 - Problem_Set3_Question2
// This code solves dT/dt=alpha*(d2T/dx2 + d2T/dy2)+S
// This code Uses a timestep of 0.01 s, a fully implicit scheme, and 80 x 80 equispaced control volumes 
// Also it can solve with unequal dx and dy
================================================================================
/*///////////////End_Header//////////////////////////////////////////////////////

//Global variables
//Geometry
#define nx 20
#define ny 20
#define lx 0.05
#define ly 0.05
#define N nx*ny
float A[N][N],b[N],T[N],T0[N],TF[N];
float deltax,deltay,deltat,rx,ry,alpha,Landa,cp,p,h,h_side,t,Tinf,T00,TS,C1x,C2x,C1y,C2y,Tol_ss,Res,sum,time;
int time_step;


//Decleration of Functions
void Zero_A_b ();
void Initial_C ();
void Cout_Ax_b ();
void Cout_T ();
void Cout_T0 ();
void Midle_Points ();
void West_3_Points ();
void East_3_Points ();
void North_3_Points ();
void South_3_Points ();
void South_West_Point ();
void South_East_Point ();
void North_West_Point ();
void North_East_Point ();
void Linear_Solver ();
void Calculate_Res ();
void Update_T ();
void matmul(float rr[N], float AA[N][N], float TT[N],float bb[N]);
float dot_product(float rr_hat[N],float rr[N]);
void matmul_2(float vv[N], float AA[N][N], float yy[N]);
//

////////////////Begin_Main code///////////////////////////////////////////////////
int main()
{	
//Variables_Decleration
deltat=0.01;
time=0.0;
time_step=0;
t=0.0035;//thickness
Landa=14.6; //Landa
p=1716.0; //density
cp=4817.0; //Cp
h=472.0; //h
h_side=36.4/(p*cp);//h_side/pcp
T00=298.15;// initial temp
Tinf=298.15;// infinite temp
TS=423.15;// B.C in south boundary
Tol_ss=0.05;
Res=10.0;

alpha=Landa/(p*cp);
cout<<endl<<"alpha="<<alpha;
deltax=lx/nx;
deltay=ly/ny;
cout<<endl<<"deltax="<<deltax;
cout<<endl<<"deltay="<<deltay;

rx=alpha*(deltat/(deltax*deltax));
ry=alpha*(deltat/(deltay*deltay));

cout<<endl<<"rx="<<rx;
cout<<endl<<"ry="<<ry;

C1x=(h*Tinf/Landa)/((1.0/deltax)+h/(2.0*Landa));
C2x=((1.0/deltax)-h/(2.0*Landa))/((1.0/deltax)+h/(2.0*Landa));
C1y=(h*Tinf/Landa)/((1.0/deltay)+h/(2.0*Landa));
C2y=((1.0/deltay)-h/(2.0*Landa))/((1.0/deltay)+h/(2.0*Landa));
//

//Initialize
Zero_A_b ();
Initial_C ();
//
//Cout_T0();
while (Res>Tol_ss)
{
time_step=time_step+1;
time=time+deltat;
//Create Ax=b
Midle_Points ();
West_3_Points ();
East_3_Points ();
North_3_Points ();
South_3_Points ();
South_West_Point ();
South_East_Point ();
North_West_Point ();
North_East_Point ();
//

//Solve Ax=b
Linear_Solver ();
//

Calculate_Res ();

//cout<<endl<<"Res="<<Res;
//Cout_T0 ();
Update_T ();

}


//Output 
Cout_Ax_b ();
Cout_T ();
cout<<endl<<"time="<<time<<"	"<<"time_step="<<time_step;
//

return 0;
}
///////////////End_Main code//////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////
//===========================FUNCTIONS===========================================
/////////////////////////////////////////////////////////////////////////////////

/*////////////////Begin_#///////////////////////////////////////////////////

/*///////////////End_#//////////////////////////////////////////////////////

////////////////Begin_Zero N*N///////////////////////////////////////////////////
void Zero_A_b ()
{
	for (int i=0;i<N;i++)
	{
		for (int j=0;j<N;j++)
		{
			A[i][j]=0.0;
		}
		b[i]=0.0;
	}
	
}
///////////////End_Zero N*N//////////////////////////////////////////////////////

////////////////Begin_Initial_C///////////////////////////////////////////////////
void Initial_C ()
{
	for (int i=0;i<N;i++)
	{
		T0[i]=T00;
		T[i]=0;
	}
}
///////////////End_Initial_C//////////////////////////////////////////////////////

////////////////Begin_Cout Ax=b///////////////////////////////////////////////////
void Cout_Ax_b ()
{
	cout<<endl<<"Ax=b"<<endl;
	for (int i=0;i<N;i++)
	{
		for (int j=0;j<N;j++)
		{
			cout<<A[i][j]<<" ";
		}
		cout<<" "<<"x"<<"="<<" "<<b[i];
		cout<<endl;
	}
	
}
///////////////End_Cout Ax=b//////////////////////////////////////////////////////

////////////////Begin_Cout_T///////////////////////////////////////////////////
void Cout_T ()
{
	cout<<endl<<"T="<<endl;
	for (int i=0;i<N;i++)
	{
		cout<<T[i];
		cout<<endl;
	}
	
}
///////////////End_Cout_T//////////////////////////////////////////////////////

////////////////Begin_Cout_T0///////////////////////////////////////////////////
void Cout_T0 ()
{
	cout<<endl<<"T0="<<endl;
	for (int i=0;i<N;i++)
	{
		cout<<T0[i];
		cout<<endl;
	}
	
}
///////////////End_Cout_T0//////////////////////////////////////////////////////


////////////////Begin_Midle_Points///////////////////////////////////////////////////
void Midle_Points ()
{
	for (int i=1;i<nx-1;i++)
	{
		for (int j=1;j<ny-1;j++)
		{
			A[j+i*nx][j+i*nx]=1.0+2.0*rx+2.0*ry+deltat*(2.0*h_side/t);
			A[j+i*nx][j+i*nx+1]=-rx;
			A[j+i*nx][j+i*nx-1]=-rx;
			A[j+i*nx][j+i*nx+nx]=-ry;
			A[j+i*nx][j+i*nx-nx]=-ry;
			
			b[j+i*nx]=T0[j+i*nx]+deltat*(2.0*h_side/t)*Tinf;
		}
	}	
}
///////////////End_Midle_Points//////////////////////////////////////////////////////

////////////////Begin_West_3_Points///////////////////////////////////////////////////
void West_3_Points ()
{
	for (int i=1;i<nx-1;i++)
	{
		A[i*nx][i*nx]=1.0+(2.0-C2x)*rx+2.0*ry+deltat*(2.0*h_side/t);
		A[i*nx][i*nx+1]=-rx;
		A[i*nx][i*nx+nx]=-ry;
		A[i*nx][i*nx-nx]=-ry;
			
		b[i*nx]=T0[i*nx]+deltat*(2.0*h_side/t)*Tinf+C1x*rx;
	}
}
///////////////End_West_3_Points//////////////////////////////////////////////////////

////////////////Begin_East_3_Points///////////////////////////////////////////////////
void East_3_Points ()
{
	for (int i=1;i<nx-1;i++)
	{
		A[(i+1)*nx-1][(i+1)*nx-1]=1.0+(2.0-C2x)*rx+2.0*ry+deltat*(2.0*h_side/t);
		A[(i+1)*nx-1][(i+1)*nx-1-1]=-rx;
		A[(i+1)*nx-1][(i+1)*nx-1+nx]=-ry;
		A[(i+1)*nx-1][(i+1)*nx-1-nx]=-ry;
			
		b[(i+1)*nx-1]=T0[(i+1)*nx-1]+deltat*(2.0*h_side/t)*Tinf+C1x*rx;
	}
}
///////////////End_East_3_Points//////////////////////////////////////////////////////

////////////////Begin_North_3_Points///////////////////////////////////////////////////
void North_3_Points ()
{
	for (int i=1;i<nx-1;i++)
	{
		A[(nx-1)*nx+i][(nx-1)*nx+i]=1.0+2*rx+(2.0-C2y)*ry+deltat*(2.0*h_side/t);
		A[(nx-1)*nx+i][(nx-1)*nx+i-1]=-rx;
		A[(nx-1)*nx+i][(nx-1)*nx+i+1]=-ry;
		A[(nx-1)*nx+i][(nx-1)*nx+i-nx]=-ry;
			
		b[(nx-1)*nx+i]=T0[(nx-1)*nx+i]+deltat*(2.0*h_side/t)*Tinf+C1y*rx;
	}
}
///////////////End_North_3_Points//////////////////////////////////////////////////////

////////////////Begin_South_3_Points///////////////////////////////////////////////////
void South_3_Points ()
{
	for (int i=1;i<nx-1;i++)
	{
		A[i][i]=1.0+2.0*rx+3*ry+deltat*(2.0*h_side/t);
		A[i][i+1]=-rx;
		A[i][i-1]=-rx;
		A[i][i+nx]=-ry;
					
		b[i]=T0[i]+deltat*(2.0*h_side/t)*Tinf+2.0*TS*ry;
	}		
}
///////////////End_South_3_Points//////////////////////////////////////////////////////

////////////////Begin_South_West_Point///////////////////////////////////////////////////
void South_West_Point ()
{
		int i=0;
		
		A[i][i]=1.0+(2.0-C2x)*rx+3.0*ry+deltat*(2.0*h_side/t);
		A[i][i+1]=-rx;
		A[i][i+nx]=-ry;
					
		b[i]=T0[i]+deltat*(2.0*h_side/t)*Tinf+2.0*TS*ry+C1x*rx;		
}
///////////////End_South_West_Point//////////////////////////////////////////////////////

////////////////Begin_South_East_Point///////////////////////////////////////////////////
void South_East_Point ()
{
		int i=nx-1;
		
		A[i][i]=1.0+(2.0-C2x)*rx+3.0*ry+deltat*(2.0*h_side/t);
		A[i][i-1]=-rx;
		A[i][i+nx]=-ry;
					
		b[i]=T0[i]+deltat*(2.0*h_side/t)*Tinf+2.0*TS*ry+C1x*rx;		
}
///////////////End_South_East_Point//////////////////////////////////////////////////////

////////////////Begin_North_West_Point///////////////////////////////////////////////////
void North_West_Point ()
{
		int i=(nx-1)*nx;
		
		A[i][i]=1.0+(2.0-C2x)*rx+(2.0-C2y)*ry+deltat*(2.0*h_side/t);
		A[i][i+1]=-rx;
		A[i][i-nx]=-ry;
					
		b[i]=T0[i]+deltat*(2.0*h_side/t)*Tinf+C1x*rx+C1y*ry;		
}
///////////////End_North_West_Point//////////////////////////////////////////////////////

////////////////Begin_North_East_Point///////////////////////////////////////////////////
void North_East_Point ()
{
		int i=nx*nx-1;
		
		A[i][i]=1.0+(2.0-C2x)*rx+(2.0-C2y)*ry+deltat*(2.0*h_side/t);
		A[i][i-1]=-rx;
		A[i][i-nx]=-ry;
					
		b[i]=T0[i]+deltat*(2.0*h_side/t)*Tinf+C1x*rx+C1y*ry;		
}
///////////////End_North_East_Point//////////////////////////////////////////////////////

////////////////Begin_Linear_Solver///////////////////////////////////////////////////

void Linear_Solver ()
{
//declerimg parameters
	float r[N],r_hat[N],y[N],tt[N],s[N],k[N],p[N],v[N],z[N];
	float rho,rho_old,alpha,omega,beta,resid;
	int i,j;
	bool converged ;
	
	// initial guess according to lecture content

	for (int i=0;i<N;i++)
	{
		T[i]=b[i]/A[i][i];
	}
	
	// variable initialize
	matmul (r,A,T,b);
	
	for (int i=0;i<N;i++)
	{
		r_hat[i]=r[i];
	}
	rho=1.0;
	alpha=1.0;
	omega=1.0;
	for (int i=0;i<N;i++)
	{
		v[i]=0.0;
		p[i]=0.0;
	}
	
	// precondition vector K=diag(A)
	for (int i=0;i<N;i++)
	{
		k[i]=A[i][i];
	} 
	
	// Preconditioned BICGSTAB algorithm main body
	
	converged = false;
	
	while (converged == false)// check if the norm is satisfied 
	{
	
	rho_old=rho;
	
	rho=dot_product(r_hat,r);
	beta= (rho/rho_old) * (alpha/omega);
	
	for (int i=0;i<N;i++)
	{
		p[i]=r[i] + beta * (p[i] - omega*v[i]);
		y[i]=p[i]/k[i];
	}  
  
  	matmul_2(v,A,y);
  	alpha= rho/dot_product( r_hat,v);
  	
  	for (int i=0;i<N;i++)
	{
		s[i]=r[i]-alpha*v[i];
		z[i]=s[i]/k[i];
	}  
	
  	matmul_2(tt,A,z);
  	omega= dot_product(tt,s)/dot_product(tt,tt);
  	
  	for (int i=0;i<N;i++)
	{
		T[i]=T[i] + alpha *y[i] + omega * z[i];
		r[i]=s[i] - omega * tt[i];
	} 
 	
 	resid = 0.0;
  
 	for (int i=0;i<N;i++)
 	{
 		resid=resid + r[i]*r[i];
	}
  	resid = sqrt(resid)/N;
	if( resid < 1e-9)
     converged= true;

  
	}


	
}
////////////////End_Linear_Solver///////////////////////////////////////////////////

////////////////Begin_Calculate_Res///////////////////////////////////////////////////

void Calculate_Res ()
{
	sum=0.0;
	for (int i=0;i<N;i++)
	{
		sum=sum+(T0[i]-T[i])*(T0[i]-T[i]);	
	}	
	Res=sqrt(sum)/N;
}
////////////////Begin_Calculate_Res///////////////////////////////////////////////////

////////////////Begin_Update_T///////////////////////////////////////////////////

void Update_T ()
{
	for (int i=0;i<N;i++)
	{
		T0[i]=T[i];
	}
	
}
////////////////Begin_Update_T///////////////////////////////////////////////////

////////////////Begin_matmul///////////////////////////////////////////////////

void matmul(float rr[N], float AA[N][N], float TT[N],float bb[N])
{
	float sum[N];
	for (int i=0;i<N;i++)
	{
		sum[i]=0.0;
		for (int j=0;j<N;j++)
		{
		
		sum[i]=sum[i]+A[i][j]*T[j];
		
		}
		rr[i]=b[i]-sum[i];
	}
}

///////////////End_matmul//////////////////////////////////////////////////////

////////////////Begin_matmul_2///////////////////////////////////////////////////

void matmul_2(float vv[N], float AA[N][N], float yy[N])
{
	float sum[N];
	for (int i=0;i<N;i++)
	{
		sum[i]=0.0;
		for (int j=0;j<N;j++)
		{
		
		sum[i]=sum[i]+AA[i][j]*yy[j];
		
		}
		vv[i]=sum[i];
	}
}

///////////////End_matmul_2//////////////////////////////////////////////////////

////////////////Begin_dot_product///////////////////////////////////////////////////

float dot_product(float rr_hat[N],float rr[N])
{
	float sum=0.0;
	for (int i=0;i<N;i++)
		{
		
		sum=sum+rr_hat[i]*rr[i];
		
		}
	
	return sum;
}

///////////////End_dot_product//////////////////////////////////////////////////////
