#include <iostream>
#include <string>

using namespace std;

typedef long long ll;
typedef unsigned long long ull;

struct MyExc
{
	static string text;
	MyExc() {};
	static void terminate() {
		cout << text;
		exit(-1);
	}
	static void set (string str) {
		text = str;
		set_terminate(terminate);
	}
	static void trw (string str) {
		MyExc::set(str);
		throw;
	}
};
string MyExc::text;

/**
 * Convert value to [0, mod) half-interval.
 */
void _toModCircle (ll& value, ull mod) {
	bool sign = value < 0;
	value = abs(value);
	if (value >= mod)
		value %= mod;
	if (sign && value)
		value = mod - value;
}

/**
 * Fast modulo power
 */
ll powMod (ll value, ull exp, ull mod) {
	if (exp == 0)
		return 1;

	_toModCircle(value, mod);

	ll res = 1;
	while (true) {
		if (exp & 1) {
			res *= value;
			res %= mod;
			if (exp == 1)
				return res;
		}
		value *= value;
		value %= mod;
		exp >>= 1;
	}
}

const size_t MAX_M = 100;
const size_t MAX_ROWS = MAX_M * MAX_M / 2;

struct MagicComplex
{
	static const ull MOD = 1000000007;

	/// Maintain 0 <= a < MOD
	ll a;
	bool im;

	MagicComplex (ll value = 0, bool isImaginary = false): im(isImaginary){
		_toModCircle(value, MOD);
		this->a = value;
	};

	string toString() const {
		if (!this->a)
			return "0";
		else if (!this->im)
			return to_string(this->a);
		else if (this->a == 1)
			return "i";
		else
			return string("i*") + to_string(this->a);
	}


	MagicComplex& operator-= (MagicComplex const& ro);
	friend MagicComplex operator* (MagicComplex const& lo, MagicComplex const& ro);
	MagicComplex& operator*= (MagicComplex const& ro);
	MagicComplex inverse();
};

MagicComplex& MagicComplex::operator-= (MagicComplex const& ro) {
	if (this->a && ro.a && this->im ^ ro.im)
		MyExc::trw(this->toString() + " -= " + ro.toString());

	if (ro.a)
		this->im = ro.im;

	this->a -= ro.a;
	if (this->a < 0)
		this->a += MagicComplex::MOD;

	return *this;
}

MagicComplex operator* (MagicComplex const& lo, MagicComplex const& ro) {
	ll a = lo.a * ro.a;
	if (lo.im && ro.im)
		return MagicComplex(-a, lo.im ^ ro.im);
	else
		return MagicComplex(a, lo.im ^ ro.im);
}

MagicComplex& MagicComplex::operator*= (MagicComplex const& ro) {
	*this = *this * ro;

	return *this;
}

MagicComplex MagicComplex::inverse() {
	ll a = powMod(this->a, MOD - 2, MOD);
	if (this->im && a)
		a = MOD - a;

	return MagicComplex(a, this->im);
}


MagicComplex matrix[MAX_ROWS][MAX_M+1];

/**
 * Mark all cells and fill sparse neighbours matrix.
 * `m` must be even.
 * If both `m` and `n` are even, it is better to ensure m <= n.
 *
 * For m, n = 4, 7 we will mark field like that:
 *  0  0  1  1
 *  2  2  3  3
 *  4  4  5  5
 *  6  6  7  7
 *  8  8  9  9
 * 10 10 11 11
 * 12 12 13 13
 *
 * we will get full matrix like that:
 * 1 0 i 0 0 0 0 0 0 0 0 0 0 0
 * 1 1 0 i 0 0 0 0 0 0 0 0 0 0
 * i 0 1 1 i 0 0 0 0 0 0 0 0 0
 * 0 i 0 1 0 i 0 0 0 0 0 0 0 0
 * 0 0 i 0 1 0 i 0 0 0 0 0 0 0
 * 0 0 0 i 1 1 0 i 0 0 0 0 0 0
 * 0 0 0 0 i 0 1 1 i 0 0 0 0 0
 * 0 0 0 0 0 i 0 1 0 i 0 0 0 0
 * 0 0 0 0 0 0 i 0 1 0 i 0 0 0
 * 0 0 0 0 0 0 0 i 1 1 0 i 0 0
 * 0 0 0 0 0 0 0 0 i 0 1 1 i 0
 * 0 0 0 0 0 0 0 0 0 i 0 1 0 i
 * 0 0 0 0 0 0 0 0 0 0 i 0 1 0
 * 0 0 0 0 0 0 0 0 0 0 0 i 1 1
 *
 * and sparse matrix like that:
 * 0 0 1 0 i ...
 * 0 1 1 0 i ...
 * i 0 1 1 i ...
 * i 0 1 0 i ...
 * i 0 1 0 i ...
 * i 1 1 0 i ...
 * i 0 1 1 i ...
 * i 0 1 0 i ...
 * i 0 1 0 i ...
 * i 1 1 0 i ...
 * i 0 1 1 i ...
 * i 0 1 0 i ...
 * i 0 1 0 0 ...
 * i 1 1 0 0 ...
 * .............
 */
void fillMatrix (size_t m, size_t n) {
	/// m MUST be even
	size_t colCapacity = m / 2;
	size_t totalCols = m + 1;
	size_t totalRows = colCapacity * n;

	for (size_t rowI = 0; rowI < totalRows; ++rowI) {
		MagicComplex* row = matrix[rowI];

		if (rowI >= colCapacity) {
			row[0].a = 1;
			row[0].im = true;
		}

		row[colCapacity].a = 1;

		size_t rest = rowI % m;
		if (rest && rest != m-1) {
			if (rest < colCapacity)
				row[colCapacity - 1].a = 1;
			else
				row[colCapacity + 1].a = 1;
		}

		if (rowI + colCapacity < totalRows) {
			row[totalCols - 1].a = 1;
			row[totalCols - 1].im = true;
		}
	}
}

void showMatrix (size_t m, size_t n) {
	size_t colCapacity = m / 2;
	size_t totalCols = m + 1;
	size_t totalRows = colCapacity * n;

	for (size_t rowI = 0; rowI < totalRows; ++rowI) {
		for (size_t colI = 0; colI < totalCols; ++colI)
			cout << matrix[rowI][colI].toString() << " ";
		cout << endl;
	}
}

ll det (size_t m, size_t n, bool show) {
	size_t colCapacity = m / 2;
	size_t totalCols = m + 1;
	size_t totalRows = colCapacity * n;

	MagicComplex res = MagicComplex(1);
	/// For each column nullify all elements under diagonal
	for (size_t i = 0; i < totalRows; ++i) {
		MagicComplex* baseRow = matrix[i];
		MagicComplex& diagEl = baseRow[colCapacity];
		MagicComplex diagElInv = diagEl.inverse();

		res *= diagEl;

		/// Process needed rows under `diagEl`
		for (size_t offset = 1; offset < min(colCapacity + 1, totalRows - i); ++offset) {
			MagicComplex* row = matrix[i + offset];
			if (row[colCapacity - offset].a) {
				MagicComplex mult = row[colCapacity - offset] * diagElInv;
				for (
					size_t colI = colCapacity;
					colI < min(totalCols, totalRows - i + colCapacity);
					++colI
				)
					row[colI - offset] -= baseRow[colI] * mult;
			}
		}

		if (show)
			showMatrix(m, n);
	}

	return res.a;
}

ll solve (size_t m, size_t n, bool show = false) {
	/// Ensure that m is even and if both even: m <= n
	if (m & 1) {
		if (n & 1)
			return 0;

		swap(m, n);
	}
	else if ((n & 1) == 0 && m > n)
		swap(m, n);

	/// Prepare matrix
	fillMatrix(m, n);


	if (show)
		showMatrix(m, n);

	return det(m, n, show);
}

int main(int argc, char **argv) {
	if (argc != 3) {
		cout << "There must be exactly two arguments" << endl;
		return 0;
	}

	size_t m = atoi(argv[1]), n = atoi(argv[2]);
	cout << solve(m, n) << endl;

	return 0;
}

