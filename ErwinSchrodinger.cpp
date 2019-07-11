
#include <iostream>
#include <cmath>
#include <fstream>
#include "ScalarField.hpp"
#include "ScalarField.cpp"
#include <string>
#include <omp.h>

#define DIM 2
//#define NX1 10
//#define NX2 10
//#define NGHOST 1
//#define L 1.0
//#define ISIGMASQ 10.


void rk4_step(const int NGHOST, const int NX1, const int NX2, const double dt, const double dx, const double dy,
				const double idx, const double idy, const double alpha, const double beta, ScalarField &SF,
				ScalarField &dummy1, ScalarField &dummy2, ScalarField &dummy3, ScalarField &dummy4);

double V(const int i, const int j, const int NX1, const int NX2, const double dx, const double dy, const double beta);

int main() {
	std::ifstream infile("./param.par");
	std::string line;
	std::getline(infile, line);
	const int NX1 = std::stoi(line);
	std::getline(infile, line);
	const int NX2 = std::stoi(line);
	std::getline(infile, line);
	const int NGHOST = std::stoi(line);
	std::getline(infile, line);
	const double L = std::stod(line);
	std::getline(infile, line);
	const double ISIGMASQ = std::stod(line);
	std::getline(infile, line);
	const double TIME = std::stod(line);

	const double dx = L/NX1;
	const double dy = L/NX2;
	const double idx = 1./dx;
	const double idy = 1./dy;
	ScalarField SF;
	ScalarField dummy1, dummy2, dummy3, dummy4;
	SF.NewScalarField(DIM, NX1+2*NGHOST, NX2+2*NGHOST, NGHOST);
	dummy1.NewScalarField(DIM, NX1+2*NGHOST, NX2+2*NGHOST, NGHOST);
	dummy2.NewScalarField(DIM, NX1+2*NGHOST, NX2+2*NGHOST, NGHOST);
	dummy3.NewScalarField(DIM, NX1+2*NGHOST, NX2+2*NGHOST, NGHOST);
	dummy4.NewScalarField(DIM, NX1+2*NGHOST, NX2+2*NGHOST, NGHOST);
	dummy1.Zero(); dummy2.Zero(); dummy3.Zero(); dummy4.Zero();

	for (int i=0; i<NX1+2*NGHOST; ++i) {
		for (int j=0; j<NX2+2*NGHOST; ++j) {
			double x = dx*(i+(NX1/4)+0.5-(NX1+2*NGHOST)/2.);
			double y = dy*(j+0.5-(NX2+2*NGHOST)/2.);
			double rr = x*x+y*y;
			double val = std::exp(-0.5*rr*ISIGMASQ);

			SF.SetRe(i, j, val*std::cos(10.*x));
			SF.SetIm(i, j, val*std::sin(10.*x));
		}
	}
	SF.PeriodicBoundary();

    int plotnum = 0;
	double dt = 0.03*dx*dy;
	double nt = TIME/dt;
	std::cout << "Number of steps: " << nt << std::endl;
	for (int n=0; n<nt; ++n) {
		if (n%300==0)
		{
			std::cout << n << std::endl;

			if (plotnum < 10)
			{
				SF.WriteData("modsqr_0"+std::to_string(plotnum)+".dat");
			}
			else
			{
				SF.WriteData("modsqr_"+std::to_string(plotnum)+".dat");
			}

			plotnum += 1;
		}
		rk4_step(NGHOST,NX1,NX2,dt,dx,dy,idx,idy,1.,ISIGMASQ,SF,dummy1,dummy2,dummy3,dummy4); // dt = 7.*dx*dx for courant i think
	    
	}

	


	return 0;
}




void rk4_step(const int NGHOST, const int NX1, const int NX2, const double dt, const double dx, const double dy,
				const double idx, const double idy, const double alpha, const double beta, ScalarField &SF,
				ScalarField &dummy1, ScalarField &dummy2, ScalarField &dummy3, ScalarField &dummy4) {

	double RHS_RE, RHS_IM, DELTA_RE, DELTA_IM;


    #pragma omp parallel for
	for (int i=NGHOST; i<NX1+NGHOST; ++i) {
		for (int j=NGHOST; j<NX2+NGHOST; ++j) {
			RHS_RE = -dt*alpha*SF.LaplacianIm(idx,idy,i,j) + dt*V(i,j,NX1,NX2,dx,dy,beta)*SF.GetIm(i,j);
			RHS_IM = dt*alpha*SF.LaplacianRe(idx,idy,i,j) - dt*V(i,j,NX1,NX2,dx,dy,beta)*SF.GetRe(i,j);
			dummy1.SetRe(i,j,SF.GetRe(i,j)+0.5*RHS_RE);
			dummy1.SetIm(i,j,SF.GetIm(i,j)+0.5*RHS_IM);
		}
	}
	dummy1.PeriodicBoundary();

    #pragma omp parallel for
	for (int i=NGHOST; i<NX1+NGHOST; ++i) {
		for (int j=NGHOST; j<NX2+NGHOST; ++j) {
			RHS_RE = -dt*alpha*dummy1.LaplacianIm(idx,idy,i,j) + dt*V(i,j,NX1,NX2,dx,dy,beta)*dummy1.GetIm(i,j);
			RHS_IM = dt*alpha*dummy1.LaplacianRe(idx,idy,i,j) - dt*V(i,j,NX1,NX2,dx,dy,beta)*dummy1.GetRe(i,j);
			dummy2.SetRe(i,j,SF.GetRe(i,j)+0.5*RHS_RE);
			dummy2.SetIm(i,j,SF.GetIm(i,j)+0.5*RHS_IM);
		}
	}
	dummy2.PeriodicBoundary();

    #pragma omp parallel for
	for (int i=NGHOST; i<NX1+NGHOST; ++i) {
		for (int j=NGHOST; j<NX2+NGHOST; ++j) {
			RHS_RE = -dt*alpha*dummy2.LaplacianIm(idx,idy,i,j) + dt*V(i,j,NX1,NX2,dx,dy,beta)*dummy2.GetIm(i,j);
			RHS_IM = dt*alpha*dummy2.LaplacianRe(idx,idy,i,j) - dt*V(i,j,NX1,NX2,dx,dy,beta)*dummy2.GetRe(i,j);
			dummy3.SetRe(i,j,SF.GetRe(i,j)+RHS_RE);
			dummy3.SetIm(i,j,SF.GetIm(i,j)+RHS_IM);
		}
	}
	dummy3.PeriodicBoundary();
        
    #pragma omp parallel for
	for (int i=NGHOST; i<NX1+NGHOST; ++i) {
		for (int j=NGHOST; j<NX2+NGHOST; ++j) {
			RHS_RE = -dt*alpha*dummy3.LaplacianIm(idx,idy,i,j) + dt*V(i,j,NX1,NX2,dx,dy,beta)*dummy3.GetIm(i,j);
			RHS_IM = dt*alpha*dummy3.LaplacianRe(idx,idy,i,j) - dt*V(i,j,NX1,NX2,dx,dy,beta)*dummy3.GetRe(i,j);
			dummy4.SetRe(i,j,RHS_RE);
			dummy4.SetIm(i,j,RHS_IM);
		}
	}
	dummy4.PeriodicBoundary();


    #pragma omp parallel for
	for (int i=NGHOST; i<NX1+NGHOST; ++i) {
		for (int j=NGHOST; j<NX2+NGHOST; ++j) {
			DELTA_RE = (2.*dummy1.GetRe(i,j)+4.*dummy2.GetRe(i,j)+2.*dummy3.GetRe(i,j)+dummy4.GetRe(i,j)-8.*SF.GetRe(i,j))/6.;
			DELTA_IM = (2.*dummy1.GetIm(i,j)+4.*dummy2.GetIm(i,j)+2.*dummy3.GetIm(i,j)+dummy4.GetIm(i,j)-8.*SF.GetIm(i,j))/6.;
			SF.AddToRe(i,j,DELTA_RE);
			SF.AddToIm(i,j,DELTA_IM);
		}
	}
	SF.PeriodicBoundary();
}

double V(const int i, const int j, const int NX1, const int NX2, const double dx, const double dy, const double beta) {
	double x = dx*(i-0.5-NX1/2.);
	double y = dy*(j-0.5-NX2/2.);
	double rr = x*x+y*y;
	if (x>-0.05 and x<0.05)
		if (y>0.05 or y < -0.05 or (0==1))
			return 500000. + rr*beta*beta;
	return rr*(beta*beta);
	//return rr*(beta*beta);
}
