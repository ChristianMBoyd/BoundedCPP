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

// builds the full 2D planar wavevector integral of Pi0
std::complex<double> Loss::Pi0(double q, double w, double delta, double Qn, double Qnp)
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
double Loss::xiScale(double qx, double qy, double ex, double ey, double ez)
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
		vec[counter] = (1 / L) * Pi0(q, w, delta, Qnp, Qnp + Qn); // Qn is the offset, Qnp is summed over
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

	// size of parity-restricted, positive entries
	const int tot = int(std::floor((nMax + 1) / 2)) + int(std::floor((evenMax + evenPar) / 2));

	Eigen::VectorXi posList(tot); // placeholder list
	int counter = 0;
	while (counter < tot)
	{
		posList[counter] = 1 - evenPar + 2 * counter; // enumerate entries
		counter++;
	}

	return posList;
}

// enumerate the full list of parity-restricted integers
Eigen::VectorXi Loss::intList(const int parity, const int nMax)
{
	bool evenPar = evenQ(parity); // if even, avoid double-counting zero element
	Eigen::VectorXi posList = Loss::posList(parity, nMax); // build off posList

	// negate, toss 0 term if even (avoid double-counting), then reverse order 
	Eigen::VectorXi negList = -posList(Eigen::lastN(posList.size() - evenPar).reverse(), Eigen::all);

	Eigen::VectorXi intList(posList.size() + negList.size());
	intList << negList, posList; // increments from most negative to positive parity-restricted integers

	return intList;
}

// using inversion symmetry, fill in negative-integer indices with posList results
Eigen::VectorXi Loss::posToIntList(const int parity, const int nMax)
{
	bool evenMax = evenQ(nMax);
	bool evenPar = evenQ(parity);

	const int tot = int(std::floor((nMax + 1) / 2)) + int(std::floor((evenMax + evenPar) / 2)); // total size of list

	Eigen::VectorXi fList(tot); // list to hold entries counting from 0 to tot

	int counter = 0;
	while (counter < tot) // fill in index positions of posList() vector
	{
		fList[counter] = counter;
		counter++;
	}

	// store reverse-order of fList, removing the double-counting of the Q=0 term in the case that the parity is even
	Eigen::VectorXi rList = fList(Eigen::lastN(fList.size() - evenPar).reverse(), Eigen::all);

	Eigen::VectorXi posToIntList(fList.size() + rList.size()); // vector to hold the result of joining rList and fList
	posToIntList << rList, fList; // Eigen means of joining lists

	return posToIntList;
}

// planar scaling of q with mass anisotropy
double Loss::qScale(double qx, double qy, double mx)
{
	return std::sqrt(pow(qx, 2) / mx + mx * pow(qy, 2)); // real-valued root argument 
}

// anisotropic re-scaling of effective L on internal sums
double Loss::LScale(double L, double mz)
{
	return L * std::sqrt(mz); // real-valued root argument
}

// the unique "diagonal" values of mChi0
Eigen::VectorXcd Loss::mChi0DiagList(double q, double w, double delta, Eigen::VectorXd& Qlist, double L)
{
	// Qlist passed by reference as per Eigen suggestions
	const int size = Qlist.size();

	Eigen::VectorXcd diagList(size);

	int counter = 0;
	while (counter < size)
	{
		diagList[counter] = L * sumPi0(q, w, delta, Qlist[counter], L);
		counter++;
	}

	return diagList;
}

// a lower-triangular matrix holding the unique entries of the (symmetric) off-diagonal elements of mChi0
Eigen::MatrixXcd Loss::mChi0OffDiagList(double q, double w, double delta, Eigen::VectorXd& Qlist)
{
	const int size = Qlist.size();

	Eigen::MatrixXcd offDiagList(size, size);

	int n = 0; // inner loop
	int m = 0; // outer loop

	while (m < size)
	{
		n = 0; // reset for additional loops
		while (n <= m)
		{
			// symmetric combination of out-of-plane wavevector arguments of Pi0
			offDiagList(m, n) = Pi0(q, w, delta, (Qlist[m] + Qlist[n]) / 2, (Qlist[m] - Qlist[n]) / 2)
				+ Pi0(q, w, delta, (Qlist[m] - Qlist[n]) / 2, (Qlist[m] + Qlist[n]) / 2);
			n++;
		}
		m++;
	}

	return offDiagList;
}

// check if the Qlist entry is zero, requires same L (possibly scaled) fed into Qlist 
bool Loss::zeroQ(double Q, double L)
{

	return (Q < pi / L); // concern over "Q == 0" check, only value Q<pi/L is Q=0.
}

// the non-interacting density response matrix, spans negative and positive integral indices
Eigen::MatrixXcd Loss::mChi0(double qx, double qy, double mx, double mz, double w, double delta, double L, double cutoff, const int parity)
{
	const double qs = qScale(qx, qy, mx);
	const double Ls = LScale(L, mz);
	const int nMax = Loss::nMax(cutoff, L);
	
	// positive out-of-plane wavevector list, cast to double explicitly for VectorXd initialization
	Eigen::VectorXd Qlist = (pi / Ls) * posList(parity, nMax).cast<double>();
	const Eigen::VectorXi p2iList = posToIntList(parity, nMax); // map from strictly positive integral indices in Qlist to all indices
	const int size = p2iList.size(); // the size of mChi0 when including negative integral indices

	// enumerate unique entries from Q inversion symmetry + mChi0 being symmetric in Qlist indices
	Eigen::VectorXcd diagList = mChi0DiagList(qs, w, delta, Qlist, Ls); // pass scaled qs and Ls to sumPi0
	Eigen::MatrixXcd offDiagList = mChi0OffDiagList(qs, w, delta, Qlist); // pass scaled qs to Pi0

	// enumerate entries across all integral indices of mChi0
	Eigen::MatrixXcd mChi0(size, size); // mChi0 matrix to be filled in
	std::complex<double> offDiagVal; // holder for offDiagList(m,n) depending upon m<=n conditions
	std::complex<double> diagVal; // holder for the "diagonal" contribution
	int m = 0; // outer loop
	int n = 0; // inner loop
	while (m < size)
	{
		n = 0; // reset inner loop
		while (n < size)
		{
			// "diagonal" part, sign not important
			if (p2iList[n] == p2iList[m]) // Qlist and diagList are index by p2iList, not m or m directly
			{
				// if Qlist[n] = 0, need to be double-counted -- done by casting int -> complex<double> first
				diagVal = (std::complex<double>(1 + zeroQ(Qlist[p2iList[n]], Ls), 0.0)) * diagList[p2iList[n]];
			}
			else
			{
				diagVal = 0;
			}

			// off-diagonal part, entries stored in lower-triangular offDiagList
			if (p2iList[n] > p2iList[m]) // offDiagList is indexed through p2iList, not n or m directly
			{
				offDiagVal = offDiagList(p2iList[n], p2iList[m]); // in this case, interchange order of (m,n)
			}
			else
			{
				offDiagVal = offDiagList(p2iList[m], p2iList[n]); // follows lower-triangular restriction, normal order
			}

			// overall factor of (1/2) -- CAUTION: requires 1.0 to avoid integer division (i.e., 1/2==0)
			mChi0(m, n) = (std::complex<double>(1.0 / 2, 0.0)) * (diagVal - offDiagVal);
			n++;
		}
		m++;
	}

	return mChi0;
}
