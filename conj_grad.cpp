#include <iostream>
#include <fstream>
#include <cmath>
#include "Timer.h"
#include "grid.cpp"
#define PI 3.14159265358979323846

using namespace std;

double vcycle();


int main(){
//int argc,char*argv[]

	int nx=2048;
	int ny=2048;
	int c=100;
	int eps=-1;
	double omega=0.8;
	double hx=2.0/(nx-1);
	double hy=1.0/(ny-1);
	double const1=4.0*PI*PI;
	double const2=1.0/(hx*hx);
	double const3=1.0/(hy*hy);
	double const4=(2.0*const2)+(2.0*const3)+const1;
	double const5=1.0/const4;
	double c1= sinh(2*PI);


	double *u= new double[nx*ny];
	double *f=new double[nx*ny];
	double *r=new double[nx*ny];
	double *ep=new double[nx*ny];


setuftozero();

setf();

setboundary();

n1= norm_u();


for (int it=1; it< maxit; it++){

	vcycle(u, f);

	n2= norm_u();

	if(abs(n2 - n1) < 1e-8){
		cout << "convergence!" << endl;
		goto Finish1;
	}
	
	n1 = n2;
}

Finish1:

std::ofstream myfile;
//, std::ios::app)
myfile.open("solution.txt");
for (int i=0; i< (ny-1); i++){
			for (int j=0; j< (nx-1); j++){
            if(i==0 & j==0){
            myfile<< "# x y u(x,y)" << "\n";
            }
			myfile<< (hx*j) << "\t" << (hy*i) << "\t" << u[i*nx+j] << "\n";
			}
}
myfile.close();

double time = 100;
siwir::Timer timer;
time = std::min(time, timer.elapsed());

delete [] u;
delete [] f;
delete [] r;

return 0;
}


double vcycle(double u, double f){

	hx=2.0/(nx-1);
	hy=1.0/(ny-1);
	const1=4.0*PI*PI;
	const2=1.0/(hx*hx);
	const3=1.0/(hy*hy);
	const4=(2.0*const2)+(2.0*const3)+const1;
	const5=1.0/const4;
	c1= sinh(2*PI);

	
	if (nx == 64 | nx ==65){
	    cout << "direct solving" << endl;

		cg();
		deletedz();		
		goto Vfinish;
	}

	// pre-smoothing
	for (int it=1; it< nu1; it++){

	jacobi();
	
	}

	setuftozero();
	
	setr();

	// restriction
	nx = ((nx-1)/2) + 1;
	ny = ((ny-1)/2) + 1;
	
	setftozero();
	setf2();
	
	setutozero();

	vcycle(u, f);
	
	coarsenx=nx;
	coarseny=ny;
	nx=2*(coarsenx-1)+1;
	ny=2*(coarseny-1)+1;	

	setutozero();
	seteptozero();
	copyvtoep();
	interpolateep()
	v = v + ep;

	u=v;
	for (int it=1; it< nu1; it++){

		jacobi();
	}

u=v;

Vfinish:
}





