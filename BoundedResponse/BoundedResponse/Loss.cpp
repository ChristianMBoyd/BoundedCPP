#include "Loss.h"

Loss::Loss() : i(0.0, 1.0) , pi (3.141592653589793238463)
{
	// thus far, only defining const "i" and "pi" values for easy use
}

// moves the branch cut from along the negative real axis to the positive real axis
// i.e., all square roots return a positive imaginary part
std::complex<double> Loss::posRoot(std::complex<double> z)
{
	return i * std::sqrt(-z);
}

// builds the full 2D planar wavevector integral of Pi0 out by calling Pi0Part() on each term
std::complex<double> Loss::Pi0Int(double q, double w, double delta, double Qn, double Qnp)
{
	// this term is universal
	const double prefactor = 1 / (8 * pi * pow(q, 2));
	
	// placeholders for each term in the integral
	std::complex<double> nTerm;
	std::complex<double> npTerm;
	
	// placeholders for variables (re)used to construct nTerm and npTerm
	double QQ, kMax;
	std::complex<double> wVal, rootArg;

	// calculate nTerm
	QQ = pow(Qn, 2);
	if (QQ > 1) // check if integral can even be nonzero
	{
		nTerm = 0;
	}
	else
	{
		kMax = std::sqrt(1 - QQ); // avoid entire branch cut issues if sqrt() acts on real values
		wVal = w + i * delta - (pow(q, 2) + pow(Qnp, 2) - pow(Qn, 2));
		rootArg = pow(wVal, 2) - 4 * pow(kMax * q, 2);
		// combine to produce nTerm integral result
		nTerm = prefactor * (wVal - posRoot(rootArg));
	}

	// calculate npTerm
	QQ = pow(Qnp, 2);
	if (QQ > 1) // check if integral can even be nonzero
	{
		npTerm = 0;
	}
	else
	{
		kMax = std::sqrt(1 - QQ); // avoid entire branch cut issues if sqrt() acts on real values
		// This term is different! There is now a (-) sign on pow(q,2).
		wVal = w + i * delta - (-pow(q, 2) + pow(Qnp, 2) - pow(Qn, 2));
		rootArg = pow(wVal, 2) - 4 * pow(kMax * q, 2);
		// combine to produce npTerm integral result
		npTerm = prefactor * (wVal - posRoot(rootArg));
	}

	return nTerm - npTerm;
}

// the effective planar wavevector as seen within the dielectric system
double Loss::Xi(double qx, double qy, double ex, double ey, double ez)
{
	double val = ex * pow(qx, 2) + ey * pow(qy, 2);
	val = val / ez;
	val = std::sqrt(val); // real-valued

	return val;
}

// a sum over a single Qn in Pi0Int across the integer indices enumerated by min and max
std::complex<double> Loss::sumPi0Interval(double q, double w, double delta, double Qn, double L, const int min, const int max)
{
	std::complex<double> sum; // placeholder for sum return
	const double waveFactor = pi / L; // convert integer indices into wavevectors
	const int size = max - min + 1; // total number of integer indices to sum over, including start and end points.
	double Qnp; // wavevector to be summed over

	Eigen::VectorXcd vec(size); // vector of entries

	int counter = 0;
	while (counter < size) // loop to enumerate entries
	{
		Qnp = waveFactor * (counter + min); // starts at -inner, ends at +inner
		vec[counter] = (1 / L) * Pi0Int(q, w, delta, Qnp, Qnp + Qn); // Qn is the offset, Qnp is summed over
		counter++;
	}

	sum = vec.sum();
	return sum;
}

// the internal sum that contributes to the diagonal part of Pi0, built using sumPi0Interval
std::complex<double> Loss::sumPi0(double q, double w, double delta, double Qn, double L)
{
	// integer indices within (-inner, inner) enumerate wavevectors that always have non-zero entries
	const int inner = int(std::floor(L / pi)); // this will never exceed O(10^5)

	std::complex<double> sum = sumPi0Interval(q, w, delta, Qn, L, -inner, inner); // sum from -inner to inner integer indices

	// now need to include the non-zero entries where (Qn + Qnp) is within (-1,1)
	const int min = int(std::ceil(-(Qn + 1) * L / pi)); // for values smaller, now Qn+Qnp is below -1 and does not contribute
	const int max = std::min(-inner - 1, int(std::floor((1 - Qn) * L / pi))); // consider overlap with innerVec range

	sum = sum + sumPi0Interval(q, w, delta, Qn, L, min, max); // sum over this secondary range

	return sum;
}

// consistent nMax function for use based off input parameters cutoff and L
int Loss::nMax(double cutoff, double L)
{
	return int(std::ceil(cutoff * L / pi)); // this is the maximum (positive) integral wavevector index for cutoff and L
}

// boolean check (eventually, bool->int) of whether or not an input is even
bool Loss::evenQ(int val)
{
	if (val % 2 == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

// enumeration of the parity-restricted, positive wavevectors for a given parity and nMax
Eigen::VectorXi Loss::posList(const int parity, const int nMax)
{
	const bool evenPar = evenQ(parity);
	const bool evenMax = evenQ(nMax);
	const int tot = int(std::floor((nMax + 1) / 2)) + int(std::floor((evenMax + evenPar) / 2));

	Eigen::VectorXi posList(tot); // placeholder list
	int counter = 0;
	while (counter < tot)
	{
		// finish!
		counter++;
	}
}
