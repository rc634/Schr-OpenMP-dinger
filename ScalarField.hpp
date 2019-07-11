#ifndef SCALARFIELD_HPP_
#define SCALARFIELD_HPP_
#include <string>

class ScalarField {
private:
	int dim_ = 0;
	int nx1_ = 0;
	int nx2_ = 0;
	int size_ = 0;
	int ng_ = 0;
	double* re_;
	double* im_;
	const double ootw = 1./12.;
	const double fth = 4./3.;
	const double oonnt = 1./90.;
	const double fnoet = 49./18.;
	const double tht = 3./20.;
public:
	ScalarField() {};
	~ScalarField() {
		delete[] re_;
		delete[] im_;
	};
	void NewScalarField(const int dim, const int nx1, const int ng) {
		dim_ = dim;
		nx1_ = nx1;
		size_ = nx1;
		ng_ = ng;
		re_ = new double[size_];
		im_ = new double[size_];
	};

	void NewScalarField(const int dim, const int nx1, const int nx2, const int ng) {
		dim_ = dim;
		nx1_ = nx1;
		nx2_ = nx2;
		size_ = nx1*nx2;
		ng_ = ng;
		re_ = new double[size_];
		im_ = new double[size_];
	};

	int Index(const int i, const int j);
	double GetRe(const int i);
	double GetRe(const int i, const int j);
	double GetIm(const int i);
	double GetIm(const int i, const int j);
	void SetRe(const int i, const double val);
	void SetRe(const int i, const int j, const double val);
	void SetIm(const int i, const double val);
	void SetIm(const int i, const int j, const double val);
	void AddToRe(const int i, const double val);
	void AddToRe(const int i, const int j, const double val);
	void AddToIm(const int i, const double val);
	void AddToIm(const int i, const int j, const double val);
	double LaplacianRe(const double idx, const double idy, const int i, const int j);
	double LaplacianIm(const double idx, const double idy, const int i, const int j);
	void WriteData(const std::string namemodsqr);
	void PeriodicBoundary();
	void Zero();
};

#endif // SCALARFIELD_HPP_
