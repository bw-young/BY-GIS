/////////////////////////////////////////////////////////////////////
// Computational window object for GIS analysis.                   //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/13/2022 - Brennan Young                                      //
// - created.                                                      //
// 04/21/2022 - Brennan Young                                      //
// - simplified the names of all member functions.                 //
// - added size().                                                 //
// 04/28/2022 - Brennan Young                                      //
// - added operator+ (CompWindow), operator+ (int), and operator-  //
//   (int). Great for processes with nested neighborhood           //
//   operations.                                                   //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GIS_COMPWINDOW_20220413
#define YOUNG_GIS_COMPWINDOW_20220413

#include <cmath>

namespace bygis { // Brennan Young GIS namespace

class CompWindow {
private:
	int di; // row offset
	int dj; // column offset
	int ri; // radius in rows dimension
	int rj; // radius in columns dimension
	int r2; // square of ri
	bool circ; // true for circular window, rectangular otherwise
public:
	// constructors / destructor
	CompWindow (int r=1, bool c=false);
	CompWindow (int, int, int, int);
	CompWindow (const CompWindow&);
	~CompWindow ();
	
	// operators
	CompWindow& operator= (const CompWindow&);
	bool operator() (int, int, int, int) const;
	CompWindow operator+ (const CompWindow&);
	CompWindow operator+ (int);
	CompWindow operator- (int);
	
	// operations
	int size() const;
	int r () const;
	int r_row () const;
	int r_col () const;
	int offset_row () const;
	int offset_col () const;
	bool circular () const;
	
	void offset (int, int);
	void r (int, int);
	void r (int);
	void circular (bool);
	bool contains (int, int, int, int) const;
}; // CompWindow


// CONSTRUCTORS / DESTRUCTORS ///////////////////////////////////////

CompWindow::CompWindow ( int r, bool c )
: di(0), dj(0), ri(r), rj(r), r2(r*r), circ(c)
{}

// Let rj be <= 0 for circular window
CompWindow::CompWindow ( int i, int j, int ri, int rj )
: di(i), dj(j), ri(ri), rj(rj), r2(ri*ri), circ(rj <= 0)
{}

CompWindow::CompWindow ( const CompWindow& win )
: di(win.di), dj(win.dj), ri(win.ri), rj(win.rj),
	r2(win.ri*win.ri), circ(win.circ)
{}

CompWindow::~CompWindow ()
{}


// OPERATORS ////////////////////////////////////////////////////////

CompWindow& CompWindow::operator= ( const CompWindow& win )
{
	di = win.di;
	dj = win.dj;
	ri = win.ri;
	rj = win.rj;
	circ = win.circ;
	
	return *this;
}

bool CompWindow::operator() ( int fi, int fj, int i, int j ) const
{
	return contains(fi, fj, i, j);
}

// Return a rectangular window large enough to fit this window plus
// a traversal using win, for identifying the maximum window
// encapsulating nested neighborhood operations.
CompWindow CompWindow::operator+ ( const CompWindow& win )
{
	int mini = di - ri + win.di - win.ri;
	//int maxi = di + ri + win.di + win.ri;
	int minj = dj - rj + win.dj - win.rj;
	//int maxj = dj + rj + win.dj + win.rj;
	int new_di = di + win.di;
	int new_dj = dj + win.dj;
	int new_ri = new_di - mini;
	int new_rj = new_dj - minj;
	return CompWindow (new_di, new_dj, new_ri, new_rj);
}

// Return a window with each radius increased by the given amount.
CompWindow CompWindow::operator+ ( int x )
{
	CompWindow win (*this);
	win.ri += x;
	win.rj += x;
	win.r2 = win.ri * win.ri;
	return win;
}

// Return a window with each radius decreased by the given amount.
CompWindow CompWindow::operator- ( int x )
{
	return this->operator+(-x);
}


// OPERATIONS ///////////////////////////////////////////////////////

int CompWindow::size () const
{
	if ( circ ) {
		int sum = 2*ri+1;
		for ( int d = 1; d <= ri; ++d )
			sum += 2 * (2 * floor(sqrt(ri*ri - d*d)) + 1);
		return sum;
	}
	else return (2*ri+1) * (2*rj+1);
}

int CompWindow::r () const
{
	return ri;
}

int CompWindow::r_row () const
{
	return ri;
}

int CompWindow::r_col () const
{
	return rj;
}

int CompWindow::offset_row () const
{
	return di;
}

int CompWindow::offset_col () const
{
	return dj;
}

bool CompWindow::circular () const
{
	return circ;
}

// Set window offset from the focus.
void CompWindow::offset ( int i, int j )
{
	di = i;
	dj = j;
}

// Set window radius.
void CompWindow::r ( int i, int j )
{
	ri = i;
	rj = j;
}

// Set window radius.
void CompWindow::r ( int d )
{
	ri = rj = d;
}

// Identify the window as circular.
void CompWindow::circular ( bool c )
{
	circ = c;
}

// Return true if (i,j) is within the window with (fi,fj) as the
// focus.
bool CompWindow::contains ( int fi, int fj, int i, int j ) const
{
	if ( circ ) {
		double dx = i - fi;
		double dy = j - fj;
		return (dx*dx + dy*dy) <= r2;
	}
	
	return i >= (fi + di - ri)
		&& i <= (fi + di + ri)
		&& j >= (fj + dj - rj)
		&& j <= (fj + dj + rj);
}

} // namespace gis

#endif // YOUNG_GIS_COMPWINDOW_20220413