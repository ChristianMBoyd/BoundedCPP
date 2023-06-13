#pragma once
#include <cmath>
#include <complex>
#include <vector> // CHECK IF DEPRECATED
#include "mkl.h" // MKL backend for Eigen, check if this is needed here
#define EIGEN_USE_MKL_ALL // test to use intel oneMKL subroutines in Eigen
#include "Eigen/Dense" // need linear algebra tools

class Loss
{
public:
	Loss();
	std::complex<double> posRoot(std::complex<double>);
	std::complex<double> Pi0(double q, double w, double delta, double Qn, double Qnp);
	std::complex<double> Pi0Qn(double q, double w, double delta, double Qn, double Qnp);
	std::complex<double> Pi0Qnp(double q, double w, double delta, double Qn, double Qnp);
	double xiScale(double qx, double qy, double ex, double ey, double ez);
	std::complex<double> sumPi0Interval(double q, double w, double delta, double Qn, double L, const int min, const int max);
	std::complex<double> sumPi0QnInterval(double q, double w, double delta, double Qn, double L, const int min, const int max);
	std::complex<double> sumPi0QnpInterval(double q, double w, double delta, double Qn, double L, const int min, const int max);
	std::complex<double> sumPi0(double q, double w, double delta, double Qn, double L);
	std::complex<double> sumPi0New(double q, double w, double delta, double Qn, double L);
	int nMax(double cutoff, double L);
	bool evenQ(int val);
	double qScale(double qx, double qy, double mx);
	double LScale(double L, double mz);
	Eigen::VectorXd Qlist(const double L, const int nMax, const bool evenPar, const bool evenMax);
	Eigen::MatrixXcd mChi0Diag(double q, double w, double delta, Eigen::VectorXd& Qlist, double L, const bool evenPar);
	Eigen::MatrixXcd mChi0DiagNew(double q, double w, double delta, Eigen::VectorXd& Qlist, double L, const bool evenPar);
	Eigen::MatrixXcd mChi0OffDiag(double q, double w, double delta, Eigen::VectorXd& Qlist);
	Eigen::MatrixXcd mChi0(double qs, double w, double delta, Eigen::VectorXd& Qlist, double Ls, const bool evenPar);
	Eigen::VectorXd vCoulomb(double xi, Eigen::VectorXd& Qlist, const bool evenPar);
	Eigen::MatrixXcd mCoulomb(double xi, double alpha, double parTerm, double L, Eigen::VectorXd& vCoulomb);
	Eigen::MatrixXd ImChi(double dimRPA, Eigen::MatrixXcd& mChi0, Eigen::MatrixXcd& mCoulomb);
	double parityLoss(double q, double qs, double xi, double eps, double alpha, double expTerm, double gTerm, double dimRPA,
		double w, double delta, double L, double Ls, const int nMax, const int parity);
	double loss(double qx, double qy, double mx, double mz, double ex, double ey, double ez, double w, double delta, double L,
		double cutoff, double Delta);

	// OLD/DEPRECATED --- kept if need index-gymnastics later
	Eigen::VectorXi posList(const int parity, const int nMax);
	Eigen::VectorXi intList(const int parity, const int nMax);
	Eigen::VectorXi posToIntList(const int parity, const int nMax);


	// consider putting internal functions into private
private:
	const std::complex<double> i; // sqrt[-1]
	const double pi;
};