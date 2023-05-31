#include "Loss.h"

Loss::Loss() : i(0.0, 1.0) , pi (3.141592653589793238463)
{
	// thus far, only defining const "i" and "pi" values for easy use
}

// moves the branch cut from along the negative real axis to the positive real axis
// i.e., all square roots return a positive imaginary part
std::complex<double> Loss::posRoot(std::complex<double> z)
{
	return i * std::sqrt(-z); // i.e., std::sqrt() handles std::complex<double>, just wrong branch cut
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
	std::complex<double> sum(0, 0); // placeholder for sum return, initialized to zero
	const double waveFactor = pi / L; // convert integer indices into wavevectors
	const int size = max - min + 1; // total number of integer indices to sum over, including start and end points.
	double Qnp; // wavevector to be summed over

	int counter = 0;
	while (counter < size) // loop to enumerate entries
	{
		Qnp = waveFactor * (counter + min); // starts at -inner, ends at +inner
		sum += (1 / L) * Pi0(q, w, delta, Qnp, Qnp + Qn); // Qn is the offset, Qnp is summed over
		counter++;
	}

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

	sum += sumPi0Interval(q, w, delta, Qn, L, min, max); // sum over this secondary range

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

// enumerate the full list of parity-restricted integers -- may be unused!
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

// the unique entries diagonal entries within mChi0 --- includes the "off-diagonal" contribution AND factor of (1/2)!
Eigen::MatrixXcd Loss::mChi0Diag(double q, double w, double delta, Eigen::VectorXd& Qlist, double L, const int parity)
{
	const int size = Qlist.size(); // only explicitly calculating the smaller, parity-restricted positive entries
	Eigen::VectorXcd diagVec(size); // holder for diagonal entries

	std::complex<double> diag; // sumPi0 innput
	std::complex<double> offDiag; // bare Pi0 input

	int counter = 0;

	if (evenQ(parity)) // edge case if Qlist[0] == 0
	{
		diag = 2 * L * sumPi0(q, w, delta, Qlist[counter], L); // "diagonal" part double-counted for Qlist[0] == 0
		offDiag = Pi0(q, w, delta, Qlist[counter], 0) + Pi0(q, w, delta, 0, Qlist[counter]); // off-diag NOT double-counted
		diagVec[counter] = 0.5 * (diag - offDiag); // overall factor of (1/2), unrelated to edge case
		counter++; // skip this in following loop if encountered
	}

	while (counter < size)
	{
		diag = L * sumPi0(q, w, delta, Qlist[counter], L); // sumPi0 has a 1/L scaling internally (may remove)
		offDiag = Pi0(q, w, delta, Qlist[counter], 0) + Pi0(q, w, delta, 0, Qlist[counter]); // still needs the symmetrized result
		diagVec[counter] = 0.5 * (diag - offDiag); // includes overall factor of (1/2)!
		counter++;
	}

	return diagVec.asDiagonal();
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

// the unique off-diagonal entries, already symmetrized and with 0s on the diagonal
Eigen::MatrixXcd Loss::mChi0OffDiag(double q, double w, double delta, Eigen::VectorXd& Qlist)
{
	const int size = Qlist.size();
	Eigen::MatrixXcd offDiag = Eigen::MatrixXd::Constant(size, size, 0.0); // initialize zero matrix

	int m = 0; // outer
	int n = 0; // inner
	while (m < size)
	{
		n = 0; // reset before new loop
		while (n < m) // only filling in unique entries manually up to symmetry, excluding diagonal
		{
			offDiag(m, n) = 0.5 * (Pi0(q, w, delta, (Qlist[m] + Qlist[n]) / 2, (Qlist[m] - Qlist[n]) / 2)
				+ Pi0(q, w, delta, (Qlist[m] - Qlist[n]) / 2, (Qlist[m] + Qlist[n]) / 2)); // symmetrized result, includes (1/2)!
			n++;
		}
		m++;
	}

	// fill in zero entries with symmetrized values, leaving 0s along diagonal
	offDiag += offDiag.transpose().eval(); // force in-place evaluation before addition, attempt to avoid aliasing issues

	return offDiag;
}

// mChi0 built from smaller pre-existing matrices holding only the unique entries
Eigen::MatrixXcd Loss::mChi0New(double qs, double w, double delta, Eigen::VectorXd& Qlist, double Ls, const int parity)
{
	// enumerate minimal entries
	Eigen::MatrixXcd diag = mChi0Diag(qs, w, delta, Qlist, Ls, parity);
	Eigen::MatrixXcd offDiag = mChi0OffDiag(qs, w, delta, Qlist);

	Eigen::MatrixXcd miniChi0 = diag - offDiag; // diag and offDiag construction includes overall factor of (1/2)

	// determine total size -- note: need parity call anyway for mChi0Diag()
	bool evenPar = evenQ(parity);
	const int size = 2 * Qlist.size() - evenPar; // spans all integral indices: avoid double-counting zero if even

	Eigen::MatrixXcd mChi0(size, size); // full matrix

	if (evenQ(parity)) // case of even matrix dimensions -- i.e., a middle row/column
	{
		const int large = Qlist.size(); // dimensions of largest unique bottom-right block
		const int small = large - 1; // dimensions of smaller sides in the rest

		mChi0.bottomRightCorner(large, large) = miniChi0; // largest block
		mChi0.bottomLeftCorner(large, small) =
			miniChi0.rightCols(small).rowwise().reverse(); // counts from bottom middle leftward, ignoring center line
		mChi0.topRightCorner(small, large) =
			miniChi0.bottomRows(small).colwise().reverse(); // counts from right middle upward, ignoring center line
		mChi0.topLeftCorner(small, small) =
			mChi0.topRightCorner(small, large).rightCols(small).rowwise().reverse(); // counts from upper middle leftward
	}
	else // case of odd matrix dimensions -- i.e., decomposable into 4 equal size block matrices
	{
		const int inner = Qlist.size(); // length of minimal set

		mChi0.bottomRightCorner(inner, inner) = miniChi0; // all positive, unchanged
		mChi0.bottomLeftCorner(inner, inner) = miniChi0.rowwise().reverse(); // counts from bottom middle leftward
		mChi0.topRightCorner(inner, inner) = miniChi0.colwise().reverse(); // counts from the right-middle upward
		mChi0.topLeftCorner(inner, inner) = 
			mChi0.topRightCorner(inner, inner).rowwise().reverse(); // counts from midpoint leftward and upward
	}
	
	return mChi0;
}

// check if the Qlist entry is zero, requires same L (possibly scaled) fed into Qlist -- may be deprecated
bool Loss::zeroQ(double Q, double L)
{

	return (Q < pi / L); // concern over "Q == 0" check, instead the only value satisfying Q<pi/L is Q=0.
}

// the parity-dep (p2iList and Qlist) non-interacting density response matrix, all integral indices (enumerated through positives by symmetry)
Eigen::MatrixXcd Loss::mChi0(double qs, double w, double delta, double Ls, const int nMax, Eigen::VectorXd& Qlist, Eigen::VectorXi& p2iList)
{
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
				diagVal = (std::complex<double>(1.0 + 1.0 * zeroQ(Qlist[p2iList[n]], Ls))) * diagList[p2iList[n]];
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

// the Coulomb matrix that enters the RPA equation at fixed parity (as set by Qlist and p2iList)
Eigen::MatrixXcd Loss::mCoulomb(double xi, double alpha, double parTerm, double L, Eigen::VectorXd& Qlist, Eigen::VectorXi& p2iList)
{
	const int size = p2iList.size(); // size reflects both 
	bool diag; // holder for if n=m
	double prefactor;
	double infactor;
	double term;
	Eigen::MatrixXcd mCoulomb(size, size);
	int m = 0; // inner loop
	int n = 0; // outer loop

	while (m < size)
	{
		n = 0; // reset inner loop
		prefactor = 1.0 / (pow(xi, 2) + pow(Qlist[p2iList[m]], 2)); // cast to double before division, only m-dependent
		while (n < size)
		{
			diag = (n == m); // identity offset
			infactor = 1.0 / (pow(xi, 2) + pow(Qlist[p2iList[n]], 2)); // cast to double before division, only n-dependent
			term = prefactor * (diag - (xi / L) * (1 - alpha) * parTerm * infactor); // actual factor
			mCoulomb(m, n) = std::complex<double>(term); // cast to complex
			n++;
		}
		m++;
	}
	
	return mCoulomb;
}

// returns the (parity-dep) imaginary part of mChi (no signs yet) as a double-valued matrix, dimRPA INCLUDES 1/L!
Eigen::MatrixXd Loss::ImChi(double dimRPA, Eigen::MatrixXcd& mChi0, Eigen::MatrixXcd& mCoulomb)
{
	const int size = mChi0.rows(); // either matrix can be used here for reference
	std::complex<double> cdimRPA = std::complex<double>(dimRPA); // cast to complex

	Eigen::MatrixXcd mRPA = Eigen::MatrixXcd::Identity(size, size) - cdimRPA * (mChi0 * mCoulomb); // RPA matrix

	Eigen::MatrixXcd mChi = (mRPA.inverse()) * mChi0; // the parity-restricted RPA result of the mChi matrix

	return mChi.imag(); // return only the imaginary part, which is a double-valued matrix.
}

// the parity-dep (through Qlist and p2iList) Coulomb vector used to evaluate the sum over ImChi
Eigen::VectorXd Loss::vCoulomb(double xi, Eigen::VectorXd& Qlist, Eigen::VectorXi& p2iList)
{
	const int size = p2iList.size(); // includes positive/negative indices

	Eigen::VectorXd vCoulomb(size);

	int n = 0;
	while (n < size)
	{
		vCoulomb[n] = 1.0 / (pow(xi, 2) + pow(Qlist[p2iList[n]], 2));
		n++;
	}

	return vCoulomb;
}

// the parity-dep sum over ImChi, gTerm includes 1/L^2 but not parity-depedence!
double Loss::parityLoss(double q, double qs, double xi, double eps, double alpha, double expTerm, double gTerm, double dimRPA,
	double w, double delta,	double L, double Ls, const int nMax, const int parity)
{
	// enumerate parity-dependent wavevectors and maps from positive values to all integral values
	Eigen::VectorXd Qlist = (pi / L) * posList(parity, nMax).cast<double>(); // unscaled positive wavevectors
	Eigen::VectorXi p2iList = posToIntList(parity, nMax); // conversion from all integral indices to the positive values in Qlist
	const int size = p2iList.size(); // size of matrices operating over all integral indices

	// build mChi0
	Eigen::VectorXd mQlist = (L / Ls) * Qlist; // mChi0 takes scaled wavevectors, but same p2iList map
	Eigen::MatrixXcd mChi0 = Loss::mChi0New(qs, w, delta,Qlist, Ls, parity);

	// build mCoulomb
	double parTerm = (1.0 - pow(-1.0, parity) * expTerm) / (1.0 - pow(1.0, parity) * alpha * expTerm); // parity-dependent Coulomb term
	Eigen::MatrixXcd mCoulomb = Loss::mCoulomb(xi, alpha, parTerm, L, Qlist, p2iList);

	// build ImChi
	Eigen::MatrixXd ImChi = Loss::ImChi(dimRPA, mChi0, mCoulomb);

	// build vCoulomb
	Eigen::VectorXd vCoulomb = Loss::vCoulomb(xi, Qlist, p2iList);


	double product = vCoulomb.transpose() * ImChi * vCoulomb; // sum over imChi converted to matrix product
	double gParTerm = pow(1.0 - pow(-1.0, parity) * expTerm, 2) * (1.0 + pow(-1.0, parity) * alpha * expTerm) 
		/ (1.0 - pow(-1.0, parity) * alpha * expTerm); // parity-dependent prefactor, not included in gTerm

	return gTerm * gParTerm * product; // includes overall and parity-dependent prefactors
}

// the total surface loss function, summed over both parity sectors
double Loss::loss(double qx, double qy, double mx, double mz, double ex, double ey, double ez, double w, double delta,
	double L, double cutoff, double Delta)
{
	// integral index cutoff
	const int nMax = Loss::nMax(cutoff, L);
	
	// momenta
	double q = std::sqrt(pow(qx, 2) + pow(qy, 2)); // unscaled planar magnitude
	double qs = qScale(qx, qy, mx); // scaled planar magnitude
	double xi = xiScale(qx, qy, ex, ey, ez); // dielectric-scaled planar magnitude
	double Ls = LScale(L, mz); // mass-scaled L

	// dielectric parameters
	double eps = ez * xi / q; // effectrive planar dielectric constant
	double alpha = (eps - 1.0) / (eps + 1.0); // dielectric contribution to finite-size effects
	double expTerm = std::exp(-xi * L); // scale of finite-size effects

	// prefactors
	double dim = 3.0 * pow(pi, 2) * pow(Delta, 2) / (2.0 * std::sqrt(mz)); // overall prefactor
	double dimRPA = dim / (ez * L); // prefactor in RPA call, includes L!
	double gTerm = pow(xi, 2) / (pow(L, 2) * q * (eps + 1.0) * (1.0 + alpha * pow(expTerm, 2))); // loss prefactor, includes 1/L^2!

	// even parity
	int parity = 0;
	double evenLoss = parityLoss(q, qs, xi, eps, alpha, expTerm, gTerm, dimRPA, w, delta, L, Ls, nMax, parity);

	// odd parity
	parity = 1;
	double oddLoss = parityLoss(q, qs, xi, eps, alpha, expTerm, gTerm, dimRPA, w, delta, L, Ls, nMax, parity);

	double loss = -dim * (evenLoss + oddLoss); // sum over both parity sectors and include dimensions through dim
	return loss;
}
