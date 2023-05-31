#pragma once
#include <cmath>
#include <complex>
#include <vector> // need variable-size arrays for below
#include "Eigen/Dense" // need linear algebra tools

class Loss
{
public:
	Loss();
	std::complex<double> posRoot(std::complex<double>);
	std::complex<double> Pi0(double q, double w, double delta, double Qn, double Qnp);
	double xiScale(double qx, double qy, double ex, double ey, double ez);
	std::complex<double> sumPi0Interval(double q, double w, double delta, double Qn, double L, const int min, const int max);
	std::complex<double> sumPi0(double q, double w, double delta, double Qn, double L);
	int nMax(double cutoff, double L);
	bool evenQ(int val);
	Eigen::VectorXi posList(const int parity, const int nMax);
	Eigen::VectorXi intList(const int parity, const int nMax); // may be unnecessary
	Eigen::VectorXi posToIntList(const int parity, const int nMax);
	double qScale(double qx, double qy, double mx);
	double LScale(double L, double mz);
	Eigen::VectorXcd mChi0DiagList(double q, double w, double delta, Eigen::VectorXd& Qlist, double L);
	Eigen::MatrixXcd mChi0Diag(double q, double w, double delta, Eigen::VectorXd& Qlist, double L, const int parity); // NEW
	Eigen::MatrixXcd mChi0OffDiagList(double q, double w, double delta, Eigen::VectorXd& Qlist);
	Eigen::MatrixXcd mChi0OffDiag(double q, double w, double delta, Eigen::VectorXd& Qlist); // NEW
	bool zeroQ(double Q, double L);
	Eigen::MatrixXcd mChi0(double qs, double w, double delta, double Ls, const int nMax, Eigen::VectorXd& Qlist, 
		Eigen::VectorXi& p2iList);
	Eigen::MatrixXcd mChi0New(double qs, double w, double delta, Eigen::VectorXd& Qlist, double Ls, const int parity); // NEW
	Eigen::MatrixXcd mCoulomb(double xi, double alpha, double parTerm, double L, Eigen::VectorXd& Qlist, Eigen::VectorXi& p2iList);
	Eigen::MatrixXd ImChi(double dimRPA, Eigen::MatrixXcd& mChi0, Eigen::MatrixXcd& mCoulomb);
	Eigen::VectorXd vCoulomb(double xi, Eigen::VectorXd& Qlist, Eigen::VectorXi& p2iList);
	double parityLoss(double q, double qs, double xi, double eps, double alpha, double expTerm, double gTerm, double dimRPA,
		double w, double delta, double L, double Ls, const int nMax, const int parity);
	double loss(double qx, double qy, double mx, double mz, double ex, double ey, double ez, double w, double delta, double L,
		double cutoff, double Delta);

	// consider putting internal functions into private
private:
	const std::complex<double> i; // sqrt[-1]
	const double pi;
};