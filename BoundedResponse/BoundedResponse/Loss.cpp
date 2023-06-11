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

// Consider: how to pull the logic out of this --- i.e., decide from Qlist which entries to compute, then push that work to GPU
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

	for (int counter = 0; counter < size; counter++) // loop to enumerate entries
	{
		Qnp = waveFactor * (counter + min); // starts at -inner, ends at +inner
		sum += (1 / L) * Pi0(q, w, delta, Qnp, Qnp + Qn); // Qn is the offset, Qnp is summed over
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

// enumeration of the parity-restricted, positive wavevectors for a given parity and nMax -- DEPRECATED
Eigen::VectorXi Loss::posList(const int parity, const int nMax)
{
	const bool evenPar = evenQ(parity);
	const bool evenMax = evenQ(nMax);

	// size of parity-restricted, positive entries
	const int tot = int(std::floor((nMax + 1) / 2)) + int(std::floor((evenMax + evenPar) / 2));

	Eigen::VectorXi posList(tot); // placeholder list
	for(int counter = 0; counter < tot; counter++)
	{
		posList[counter] = 1 - evenPar + 2 * counter; // enumerate entries
	}

	return posList;
}

// enumerate the full list of parity-restricted integers -- DEPRECATED
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

// using inversion symmetry, fill in negative-integer indices with posList results -- DEPRECATED
Eigen::VectorXi Loss::posToIntList(const int parity, const int nMax)
{
	bool evenMax = evenQ(nMax);
	bool evenPar = evenQ(parity);

	const int tot = int(std::floor((nMax + 1) / 2)) + int(std::floor((evenMax + evenPar) / 2)); // total size of list

	Eigen::VectorXi fList(tot); // list to hold entries counting from 0 to tot

	for (int counter = 0; counter < tot; counter++) // fill in index positions of posList() vector
	{
		fList[counter] = counter;
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

// non-negative wavevector list for even OR odd parity --- parity logic is pushed out into function calls via evenPar and evenMax
Eigen::VectorXd Loss::Qlist(const double L, const int nMax, const bool evenPar, const bool evenMax)
{
	// size of parity-restricted, positive entries
	const int tot = int(std::floor((nMax + 1) / 2)) + int(std::floor((evenMax + evenPar) / 2));

	Eigen::VectorXd Qlist(tot); // placeholder list
	const double dQ = pi / L; // wavevector spacing
	for (int counter = 0; counter < tot; counter++)
	{
		Qlist[counter] = dQ * (1 - evenPar + 2 * counter); // non-negative wavevector entries
	}

	return Qlist;
}

// the unique entries diagonal entries within mChi0
//	additionall handles the "off-diagonal" contribution, factor of (1/2), *and* parity dependence via (external) evenPar
Eigen::MatrixXcd Loss::mChi0Diag(double q, double w, double delta, Eigen::VectorXd& Qlist, double L, const bool evenPar)
{
	const int size = Qlist.size(); // only explicitly calculating the smaller, parity-restricted positive entries
	Eigen::VectorXcd diagVec(size); // holder for diagonal entries

	std::complex<double> diag; // sumPi0 innput
	std::complex<double> offDiag; // bare Pi0 input

	// mChi0Diag[0] has an edge case for even parity, the evenPar term handles this below
	diag = (1+evenPar) * L * sumPi0(q, w, delta, Qlist[0], L); // "diagonal" part double-counted for Qlist[0] == 0
	offDiag = Pi0(q, w, delta, Qlist[0], 0) + Pi0(q, w, delta, 0, Qlist[0]); // off-diag NOT double-counted
	diagVec[0] = 0.5 * (diag - offDiag); // overall factor of (1/2), unrelated to edge case

	// fill in the rest of the entries, no more edge cases to worry about
	for (int counter = 1; counter < size; counter++)
	{
		diag = L * sumPi0(q, w, delta, Qlist[counter], L); // sumPi0 has a 1/L scaling internally (may remove)
		offDiag = Pi0(q, w, delta, Qlist[counter], 0) + Pi0(q, w, delta, 0, Qlist[counter]); // still needs the symmetrized result
		diagVec[counter] = 0.5 * (diag - offDiag); // includes overall factor of (1/2)!
	}

	return diagVec.asDiagonal();
}

// the unique off-diagonal entries, already symmetrized and with 0s on the diagonal
Eigen::MatrixXcd Loss::mChi0OffDiag(double q, double w, double delta, Eigen::VectorXd& Qlist)
{
	const int size = Qlist.size();
	Eigen::MatrixXcd offDiag = Eigen::MatrixXd::Constant(size, size, 0.0); // initialize zero matrix

	for (int m = 0; m < size; m++)
	{
		for (int n = 0; n < m; n++) // only filling in unique entries manually up to symmetry, excluding diagonal
		{
			offDiag(m, n) = 0.5 * (Pi0(q, w, delta, (Qlist[m] + Qlist[n]) / 2, (Qlist[m] - Qlist[n]) / 2)
				+ Pi0(q, w, delta, (Qlist[m] - Qlist[n]) / 2, (Qlist[m] + Qlist[n]) / 2)); // symmetrized result, includes (1/2)!
		}
	}

	// fill in zero entries with symmetrized values, leaving 0s along diagonal
	auto init = offDiag; // forced evaluation anyway, manual temporary
	offDiag.transposeInPlace(); // transpose

	return init + offDiag; // combine offDiag with its tranpose, maintaining 0s on the diagonal
}

// Consider: explicit calls to even/odd parity cases since will be done anyway --- remove logic
// mChi0 built from smaller pre-existing matrices holding only the unique entries
Eigen::MatrixXcd Loss::mChi0(double qs, double w, double delta, Eigen::VectorXd& Qlist, double Ls, const bool evenPar)
	{
		// enumerate minimal entries
	auto diag = mChi0Diag(qs, w, delta, Qlist, Ls, evenPar);
	auto offDiag = mChi0OffDiag(qs, w, delta, Qlist);

	auto miniChi0 = diag - offDiag; // unique (by symmetry) entries, makes up bottom-right corner of full mChi0

	// determine total size, parity logic off-loaded through evenPar
	const int size = 2 * Qlist.size() - evenPar; // spans all integral indices: avoid double-counting zero if even

	Eigen::MatrixXcd mChi0(size, size); // full matrix dimensions, parity-dependent via evenPar in size def
	const int large = Qlist.size(); // dimensions of bottom-right block
	const int small = large - evenPar; // dimensions of smaller side in the case of even parity, avoids double-counting diagonals

	// filling in mChi0 from miniChi0
	mChi0.bottomRightCorner(large, large) = miniChi0; // largest block
	mChi0.bottomLeftCorner(large, small) =
		miniChi0.rightCols(small).rowwise().reverse(); // counts from bottom middle leftward, ignoring center line
	mChi0.topRightCorner(small, large) =
		miniChi0.bottomRows(small).colwise().reverse(); // counts from right middle upward, ignoring center line
	mChi0.topLeftCorner(small, small) =
		mChi0.topRightCorner(small, large).rightCols(small).rowwise().reverse(); // counts from upper middle leftward

	return mChi0;
}

// the Coulomb entries vector used to evaluate the sum over ImChi -- parity-dependence called externally through evenPar
Eigen::VectorXd Loss::vCoulomb(double xi, Eigen::VectorXd& Qlist, const bool evenPar)
{
	const int size = Qlist.size(); // enumerating only the unique entries by symmetry

	Eigen::VectorXd fList(size); // forward counting list

	const double xi2 = xi * xi; // just saving space below
	for (int counter = 0;  counter < size; counter++)
	{
		double Qval = Qlist[counter];
		fList[counter] = 1.0 / (xi2 + Qval * Qval);
	}

	// store reverse-order of fList, avoiding double-counting Q=0 term if even parity
	Eigen::VectorXd rList = fList(Eigen::lastN(fList.size() - evenPar).reverse(), Eigen::all);

	Eigen::VectorXd vCoulomb(fList.size() + rList.size()); // vector to hold the result of joining rList and fList
	vCoulomb << rList, fList; // Eigen means of joining lists

	return vCoulomb;
}

// the Coulomb matrix that enters the RPA equation at fixed parity (as set by Qlist and p2iList)
Eigen::MatrixXcd Loss::mCoulomb(double xi, double alpha, double parTerm, double L, Eigen::VectorXd& vCoulomb)
{
	const double internalFactor = (xi / L) * (1 - alpha) * parTerm; // multiplies the non-diagonal part (i.e., result is *not* proportional)

	Eigen::MatrixXd mCoulomb = vCoulomb.asDiagonal();
	mCoulomb.noalias() -= internalFactor * vCoulomb * vCoulomb.transpose(); // attempted .noalias() optimization

	return mCoulomb; // slow, but Eigen->MKL call may require this double->std::complex<double> cast
}

// returns the (parity-dep) imaginary part of mChi (no signs yet) as a double-valued matrix, dimRPA INCLUDES 1/L!
Eigen::MatrixXd Loss::ImChi(double dimRPA, Eigen::MatrixXcd& mChi0, Eigen::MatrixXcd& mCoulomb)
{
	const int size = mChi0.rows(); // either matrix can be used here for reference
	// RPA matrix acting on mChi, double->std::complex<double> cast appears to be requied for Eigen->MKL optimization
	Eigen::MatrixXcd mRPA = Eigen::MatrixXcd::Identity(size, size) - std::complex<double>(dimRPA) * (mChi0 * mCoulomb);

	// this tells Eigen to solve the linear algebra problem mRPA*mChi = mChi0 for mChi
	Eigen::MatrixXcd mChi = mRPA.partialPivLu().solve(mChi0); // needs an explicit MatrixXcd cast for .imag() call on return

	// currently crashes above, mkl seems to not be linking -- good luck

	return mChi.imag(); // return only the imaginary part, which is a real-valued <double>
}

// the parity-dep sum over ImChi, gTerm includes 1/L^2 but not parity-depedence!
double Loss::parityLoss(double q, double qs, double xi, double eps, double alpha, double expTerm, double gTerm, double dimRPA,
	double w, double delta,	double L, double Ls, const int nMax, const int parity)
{
	// parity logic performed here, the remaining Loss:: functions are purely mathematical and call these values as integers
	bool evenPar = evenQ(parity);
	bool evenMax = evenQ(nMax);

	// enumerate non-negative wavevector entries -- negative values calculated by symmetry
	auto Qlist = Loss::Qlist(L, nMax, evenPar, evenMax);

	// build mChi0
	Eigen::VectorXd mChi0Qlist = (L / Ls) * Qlist; // mChi0 takes scaled wavevectors
	auto mChi0 = Loss::mChi0(qs, w, delta, mChi0Qlist, Ls, evenPar);

	// build vCoulomb
	auto vCoulomb = Loss::vCoulomb(xi, Qlist, evenPar);

	// build mCoulomb
	double parTerm = (1.0 - pow(-1.0, parity) * expTerm) / (1.0 - pow(1.0, parity) * alpha * expTerm); // passed to mCoulomb
	auto mCoulomb = Loss::mCoulomb(xi, alpha, parTerm, L, vCoulomb);

	// build ImChi
	auto ImChi = Loss::ImChi(dimRPA, mChi0, mCoulomb);


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
	
	// length scales
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
