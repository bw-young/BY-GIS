/////////////////////////////////////////////////////////////////////
// GridCell object to facilitate handling raster grid cells.       //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/28/2022 - Brennan Young                                      //
// - created.                                                      //
// 05/17/2022 - Brennan Young                                      //
// - added operator==() and operator!=().                          //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GIS_GRIDCELL_20220428
#define YOUNG_GIS_GRIDCELL_20220428

#include <cmath>

namespace bygis { // Brennan Young GIS namespace

class GridCell {
public:
	// neighborhood constants
	static const unsigned char QUEEN;
	static const unsigned char ROOK;
	static const unsigned char BISHOP;
	
	int i;
	int j;
	int k;
	
	// constructors / destructor
	GridCell (int i=-1, int j=-1, int k=0);
	GridCell (const GridCell&);
	~GridCell ();
	
	// operators
	GridCell& operator= (const GridCell&);
	bool operator== (const GridCell&) const;
	bool operator!= (const GridCell&) const;
	bool operator< (const GridCell&) const;
	
	// operations
	bool isNbr (const GridCell&, const unsigned char&) const;
	bool isNbr (const GridCell&) const;
}; // GridCell

const unsigned char GridCell::QUEEN = 0;
const unsigned char GridCell::ROOK = 1;
const unsigned char GridCell::BISHOP = 2;


// CONSTRUCTORS / DESTRUCTORS ///////////////////////////////////////

GridCell::GridCell ( int row, int col, int band )
: i(row), j(col), k(band)
{}

GridCell::GridCell ( const GridCell& cell )
: i(cell.i), j(cell.j), k(cell.k)
{}

GridCell::~GridCell ()
{}


// OPERATORS ////////////////////////////////////////////////////////

GridCell& GridCell::operator= ( const GridCell& cell )
{
	i = cell.i;
	j = cell.j;
	k = cell.k;
	return *this;
}

bool GridCell::operator== ( const GridCell& cell ) const
{
	return i == cell.i && j == cell.j && k == cell.k;
}

bool GridCell::operator!= ( const GridCell& cell ) const
{
	return i != cell.i || j != cell.j || k != cell.k;
}

bool GridCell::operator< ( const GridCell& cell ) const
{
	return k < cell.k
		|| (k == cell.k && i < cell.i)
		|| (k == cell.k && i == cell.i && j < cell.j);
}


// OPERATIONS ///////////////////////////////////////////////////////

bool GridCell::isNbr ( const GridCell& cell, const unsigned char& nbrhd ) const
{
	int di = abs(i - cell.i);
	int dj = abs(j - cell.j);
	return !(di+dj == 0 || di+dj > 2)
		|| nbrhd == QUEEN
		|| (nbrhd == ROOK && di+dj == 1)
		|| (nbrhd == BISHOP && di+dj == 2);
}

bool GridCell::isNbr ( const GridCell& cell ) const
{
	return isNbr(cell, QUEEN);
}


// NON-MEMBER FUNCTIONS /////////////////////////////////////////////

bool isNbr ( const GridCell& A, const GridCell& B,
	const unsigned char& nbrhd )
{
	return A.isNbr(B, nbrhd);
}

bool isNbr ( const GridCell& A, const GridCell& B )
{
	return A.isNbr(B);
}

} // namespace bygis

#endif //YOUNG_GIS_GRIDCELL_20220428