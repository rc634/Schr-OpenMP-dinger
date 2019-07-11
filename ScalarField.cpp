#include <fstream>
#include <string>
#include <iostream>
#include "ScalarField.hpp"

int ScalarField::Index(const int i, const int j) {
	return i*nx1_+j;
}

double ScalarField::GetRe(const int i) {
	return re_[i];
}

double ScalarField::GetRe(const int i, const int j) {
	const int ind = Index(i, j);
	return re_[ind];
}

double ScalarField::GetIm(const int i) {
	return im_[i];
}

double ScalarField::GetIm(const int i, const int j) {
	const int ind = Index(i, j);
	return im_[ind];
}

void ScalarField::SetRe(const int i, const double val) {
	re_[i] = val;
}

void ScalarField::SetRe(const int i, const int j, const double val) {
	const int ind = Index(i, j);
	re_[ind] = val;
}

void ScalarField::SetIm(const int i, const double val) {
	im_[i] = val;
}

void ScalarField::SetIm(const int i, const int j, const double val) {
	const int ind = Index(i, j);
	im_[ind] = val;
}

void ScalarField::AddToRe(const int i, const double val) {
	re_[i] += val;
}

void ScalarField::AddToRe(const int i, const int j, const double val) {
	const int ind = Index(i, j);
	re_[ind] += val;
}

void ScalarField::AddToIm(const int i, const double val) {
	im_[i] += val;
}

void ScalarField::AddToIm(const int i, const int j, const double val) {
	const int ind = Index(i, j);
	im_[ind] += val;
}


double ScalarField::LaplacianRe(const double idx, const double idy, const int i, const int j) {
	double lap = 0.0;
	const int ind = Index(i,j);
	if (ng_==1) {
		lap += (re_[Index(i-1,j)] + re_[Index(i+1,j)] - 2.*re_[ind])*idx*idx;
		lap += (re_[Index(i,j-1)] + re_[Index(i,j+1)] - 2.*re_[ind])*idy*idy;
	}
	if (ng_==2) {
		lap += ( -ootw*re_[Index(i+2,j)] + fth*re_[Index(i+1,j)] - 2.5*re_[ind] + fth*re_[Index(i-1,j)] - ootw*re_[Index(i-2,j)])*idx*idx;
		lap += ( -ootw*re_[Index(i,j+2)] + fth*re_[Index(i,j+1)] - 2.5*re_[ind] + fth*re_[Index(i,j-1)] - ootw*re_[Index(i,j-2)])*idy*idy;
	}
	if (ng_==3) {
		lap += ( oonnt*re_[Index(i+3,j)] - tht*re_[Index(i+2,j)] + 1.5*re_[Index(i+1,j)] - fnoet*re_[ind])*idx*idx;
		lap += ( oonnt*re_[Index(i-3,j)] - tht*re_[Index(i-2,j)] + 1.5*re_[Index(i-1,j)])*idx*idx;
		lap += ( oonnt*re_[Index(i,j+3)] - tht*re_[Index(i,j+2)] + 1.5*re_[Index(i,j+1)] - fnoet*re_[ind])*idy*idy;
		lap += ( oonnt*re_[Index(i,j-3)] - tht*re_[Index(i,j-2)] + 1.5*re_[Index(i,j-1)])*idy*idy;
	}
	return lap;
}

double ScalarField::LaplacianIm(const double idx, const double idy, const int i, const int j) {
	double lap = 0.0;
	const int ind = Index(i,j);
	if (ng_==1) {
		lap += (im_[Index(i-1,j)] + im_[Index(i+1,j)] - 2.*im_[ind])*idx*idx;
		lap += (im_[Index(i,j-1)] + im_[Index(i,j+1)] - 2.*im_[ind])*idy*idy;
	}
	if (ng_==2) {
		lap += ( -ootw*im_[Index(i+2,j)] + fth*im_[Index(i+1,j)] - 2.5*im_[ind] + fth*im_[Index(i-1,j)] - ootw*im_[Index(i-2,j)])*idx*idx;
		lap += ( -ootw*im_[Index(i,j+2)] + fth*im_[Index(i,j+1)] - 2.5*im_[ind] + fth*im_[Index(i,j-1)] - ootw*im_[Index(i,j-2)])*idy*idy;
	}
	if (ng_==3) {
		lap += ( oonnt*im_[Index(i+3,j)] - tht*im_[Index(i+2,j)] + 1.5*im_[Index(i+1,j)] - fnoet*im_[ind])*idx*idx;
		lap += ( oonnt*im_[Index(i-3,j)] - tht*im_[Index(i-2,j)] + 1.5*im_[Index(i-1,j)])*idx*idx;
		lap += ( oonnt*im_[Index(i,j+3)] - tht*im_[Index(i,j+2)] + 1.5*im_[Index(i,j+1)] - fnoet*im_[ind])*idy*idy;
		lap += ( oonnt*im_[Index(i,j-3)] - tht*im_[Index(i,j-2)] + 1.5*im_[Index(i,j-1)])*idy*idy;
	}
	return lap;
}

void ScalarField::WriteData( const std::string namemodsqr) {
	// std::cout << nx1_ << ", " << nx2_ << std::endl;
	std::ofstream outfile(namemodsqr,std::ios::binary);
	for (int i=0; i<nx1_; ++i) {
		for (int j=0; j<nx2_; ++j) {
			double val = im_[Index(i,j)]*im_[Index(i,j)] + re_[Index(i,j)]*re_[Index(i,j)];
			outfile.write(reinterpret_cast<const char*>(&val), sizeof(val));
		}
	}
	outfile.close();
}

void ScalarField::PeriodicBoundary() {
	for (int i = 0; i < nx1_; ++i) {
		for (int g = 0; g < ng_; ++g) {
			re_[Index(i,g)] = re_[Index(i,nx2_-2*ng_+g)];
			im_[Index(i,g)] = im_[Index(i,nx2_-2*ng_+g)];

			re_[Index(i,nx2_-1-g)] = re_[Index(i,2*ng_-1-g)];
			im_[Index(i,nx2_-1-g)] = im_[Index(i,2*ng_-1-g)];
		}
	}


	for (int j = 0; j < nx2_; ++j) {
		for (int g = 0; g < ng_; ++g) {
			re_[Index(g,j)] = re_[Index(nx1_-2*ng_+g,j)];
			im_[Index(g,j)] = im_[Index(nx1_-2*ng_+g,j)];

			re_[Index(nx1_-1-g,j)] = re_[Index(2*ng_-1-g,j)];
			im_[Index(nx1_-1-g,j)] = im_[Index(2*ng_-1-g,j)];
		}
	}
}


void ScalarField::Zero() {
	for (int i=0; i<nx1_; ++i) {
		for (int j=0; j<nx2_; ++j) {
			re_[Index(i,j)] = 0.;
			im_[Index(i,j)] = 0.;
		}
	}
}
