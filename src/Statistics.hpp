/////////////////////////////////////////////////////////////////////
// Object to help handle monovariate statistical reductions.       //
// Note that, when counting individual values, they are truncated  //
// to integers.                                                    //
//                                                                 //
// Code modeled by John Cook                                       //
// https://www.johndcook.com/blog/skewness_kurtosis/               //
// (accessed 5/5/2022)                                             //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 05/03-05/2022 - Brennan Young                                   //
// - created.                                                      //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GIS_STATISTICS_20220503
#define YOUNG_GIS_STATISTICS_20220503

#include <cmath>
#include <map>
#include <vector>

namespace bygis { // Brennan Young GIS namespace

class Statistics {
	long long n; // population size
	double xmin;
	double xmax;
	double xsum;
	double xsum2;
	double M1;
	double M2;
	double M3;
	double M4;
	
	bool countValues; // true to count values for median or mode
	std::map<int, int> values;
	std::map<double, int> fvalues;
	int valMaxCount;
	int xmode;
public:
	// constructors / destructor
	Statistics (bool cv=false);
	Statistics (const Statistics&);
	~Statistics ();
	
	// operators
	Statistics& operator= (const Statistics&);
	
	// operations
	void push (double); // add a sample
	template<typename T> void push (const std::vector<T>&);
	void trackValues (bool); // set whether values are tracked
	bool trackValues () const; // get whether values are being tracked
	int count () const;
	double min () const;
	double max () const;
	double range () const;
	double sum () const;
	double sum2 () const;
	double mean () const;
	double var () const;
	double stdev () const;
	double skew () const;
	double kurtosis () const;
	double median () const;
	int mode () const;
}; // Statistics


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

Statistics::Statistics ( bool cv )
: n(0), xmin(0), xmax(0), xsum(0), xsum2(0),
	M1(0), M2(0), M3(0), M4(0),
	countValues(cv), valMaxCount(0), xmode(0)
{}

Statistics::Statistics ( const Statistics& s )
: n(s.n), xmin(s.xmin), xmax(s.xmax), xsum(s.xsum), xsum2(s.xsum2),
	M1(s.M1), M2(s.M2), M3(s.M3), M4(s.M4),
	countValues(s.countValues), valMaxCount(s.valMaxCount),
	xmode(s.xmode)
{
	values = s.values;
	fvalues = s.fvalues;
}

Statistics::~Statistics ()
{}


// OPERATORS ////////////////////////////////////////////////////////

Statistics& Statistics::operator= ( const Statistics& s )
{
	if ( &s == this ) return *this;
	
	n = s.n;
	xmin = s.xmin;
	xmax = s.xmax;
	xsum = s.xsum;
	xsum2 = s.xsum2;
	M1 = s.M1;
	M2 = s.M2;
	M3 = s.M3;
	M4 = s.M4;
	
	countValues = s.countValues;
	values = s.values;
	fvalues = s.fvalues;
	valMaxCount = s.valMaxCount;
	xmode = s.xmode;
	
	return *this;
}


// OPERATIONS ///////////////////////////////////////////////////////

// Add a sample to the population being evaluated.
void Statistics::push ( double x )
{
	++n;
	if ( n == 1 || x < xmin ) xmin = x;
	if ( n == 1 || x > xmax ) xmax = x;
	xsum += x;
	xsum2 += x*x;
	
	long long n1 = n-1;
	double delta = x - M1;
	double delta_n = delta / n;
	double delta_n2 = delta_n * delta_n;
	double term1 = delta * delta_n * n1;
	M1 += delta_n;
	M4 += term1 * delta_n2 * (n*n - 3*n + 3)
		+ 6 * delta_n2 * M2 - 4 * delta_n * M3;
	M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
	M2 += term1;
	
	if ( countValues ) {
		int a = (int) x;
		if ( values.find(a) == values.end() ) values[a] = 0;
		if ( fvalues.find(x) == fvalues.end() ) fvalues[x] = 0;
		++values[a];
		++fvalues[x];
		if ( values[a] > valMaxCount ) {
			valMaxCount = values[a];
			xmode = a;
		}
	}
}

template <typename T>
void Statistics::push ( const std::vector<T>& x )
{
	for ( size_t i = 0; i < x.size(); ++i ) push(x[i]);
}

// Set whether values are tracked.
void Statistics::trackValues ( bool cv )
{
	countValues = cv;
}

// Report whether values are being tracked.
bool Statistics::trackValues () const
{
	return countValues;
}

int Statistics::count() const
{
	return n;
}

double Statistics::min () const
{
	return xmin;
}

double Statistics::max () const
{
	return xmax;
}

double Statistics::range () const
{
	return xmax - xmin;
}

double Statistics::sum () const
{
	return xsum;
}

double Statistics::sum2 () const
{
	return xsum2;
}

double Statistics::mean () const
{
	return M1;
}

double Statistics::var () const
{
	//double u = mean();
	//return xsum2 / n - u * u;
	return M2 / (n-1);
}

double Statistics::stdev () const
{
	return sqrt(var());
}

// Adjusted Fisher-Pearson standardized moment coefficient
double Statistics::skew () const
{
	// sum[ (x - mean)^3 ] / [(n-1)*stdev)^3]
	return sqrt(n) * M3 / pow(M2, 1.5);
}

double Statistics::kurtosis () const
{
	// n * sum[ (x[i]-mean)^4] / sum[ (x[i]-mean^2)^2 ]
	return n*M4 / (M2*M2) - 3.0;
}

double Statistics::median () const
{
	if ( n < 2 ) return xmode;
	long long i = 0;
	double i1 = ((double) n) / 2;
	std::map<double, int>::const_iterator it = fvalues.begin();
	for ( ; it != fvalues.end(); ++it )
	{
		if ( i + it->second > i1 ) return it->first;
		if ( fabs(i + it->second - i1) < 0.0000001 ) {
			double x = it->first;
			++it;
			return 0.5 * (x + it->first);
		}
		i += it->second;
	}
	return 0;
}

int Statistics::mode () const
{
	return xmode;
}

} // namespace bygis

#endif // YOUNG_GIS_STATISTICS_20220503