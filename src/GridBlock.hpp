/////////////////////////////////////////////////////////////////////
// GridBlock object for loaded a portion of a raster grid into     //
// memory.                                                         //
// Designed to be used as part of the Raster object.               //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/02/2022 - Brennan Young                                      //
// - created.                                                      //
// 04/08/2022 - Brennan Young                                      //
// - added data type codes, dtype member, and dataTypeSize static  //
//   method.                                                       //
// 04/26/2022 - Brennan Young                                      //
// - corrected issue where copy and assignment wouldn't copy the   //
//   entire data block.                                            //
// 04/27/2022 - Brennan Young                                      //
// - corrected issue where the BSQ index at band > 0 scaled with n //
//   instead of ni*nj.                                             //
// 05/23/2022 - Brennan Young                                      //
// - transitioned to using new DataType object.                    //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GIS_GRIDBLOCK_20220407
#define YOUNG_GIS_GRIDBLOCK_20220407

#include <iostream>

#include "DataType.hpp"

namespace bygis { // Brennan Young GIS namespace

class GridBlock {
private:
	int indexBSQ (int, int, int) const;
	int indexBIL (int, int, int) const;
	int indexBIP (int, int, int) const;
public:
	// format codes
	static const unsigned char BSQ; // band-sequential interleave
	static const unsigned char BIL; // band-interleave by line
	static const unsigned char BIP; // band-interleave by pixel
	
	int i0; // first loaded row
	int j0; // first loaded column
	int k0; // first loaded band
	int ni; // number of loaded rows
	int nj; // number of loaded columns
	int nk; // number of loaded bands
	int n;  // number of loaded elements
	DataType datatype; // data type
	unsigned char format; // interleave format
	char* data; // loaded data (binary)
	bool update; // true if changes may have been made since load
	
	// constructors / destructor
	GridBlock (int row0=0, int col0=0, int band0=0,
		int nrows=0, int ncols=0, int nbands=0,
		unsigned char datatype=4, unsigned char interleave=BSQ);
	GridBlock (const GridBlock&);
	~GridBlock ();
	
	// operators
	GridBlock& operator= (const GridBlock&);
	char* operator() (int, int, int);
	
	// operations
	bool contains (int, int, int) const;
	int index (int, int, int) const;
}; // GridBlock

const unsigned char GridBlock::BSQ = 0;
const unsigned char GridBlock::BIL = 1;
const unsigned char GridBlock::BIP = 2;


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

GridBlock::GridBlock ( int row0, int col0, int band0,
	int nrows, int ncols, int nbands,
	unsigned char dtype, unsigned char interleave )
: i0(row0), j0(col0), k0(band0), ni(nrows), nj(ncols), nk(nbands),
	format(interleave), update(false)
{
	datatype = DataType(dtype);
	if ( format != BIL && format != BIP ) format = BSQ;
	
	if ( ni < 1 || nj < 1 || nk < 1 ) ni = nj = nk = 0;
	n = ni * nj * nk;
	
	try {
		data = new char [n < 1 ? 1 : n * datatype.bytes];
	}
	catch (...) {
		std::cout << "ERR: memory allocation: constructing GridBlock\n";
		exit(1);
	}
}

GridBlock::GridBlock ( const GridBlock& block )
: i0(block.i0), j0(block.j0), k0(block.k0),
	ni(block.ni), nj(block.nj), nk(block.nk), n(block.n),
	datatype(block.datatype), format(block.format),
	update(block.update)
{
	// copy data over
	try {
		data = new char [n < 1 ? 1 : n * datatype.bytes];
	}
	catch (...) {
		std::cout << "ERR: memory allocation: copying GridBlock\n";
		exit(1);
	}
	
	for ( int i = 0; i < n * datatype.bytes; ++i )
		data[i] = block.data[i];
}

GridBlock::~GridBlock ()
{
	delete [] data;
}


// OPERATORS ////////////////////////////////////////////////////////

// Assignment.
GridBlock& GridBlock::operator= ( const GridBlock& block )
{
	if ( &block == this ) return *this;
	
	// members
	i0 = block.i0;
	j0 = block.j0;
	k0 = block.k0;
	ni = block.ni;
	nj = block.nj;
	nk = block.nk;
	n  = block.n;
	datatype = block.datatype;
	format = block.format;
	update = block.update;
	
	// copy data over
	delete [] data;
	try {
		data = new char [n < 1 ? 1 : n * datatype.bytes];
	}
	catch (...) {
		std::cout << "ERR: memory allocation: copying GridBlock\n";
		exit(1);
	}
	
	for ( int i = 0; i < n * datatype.bytes; ++i )
		data[i] = block.data[i];
	
	return *this;
}

// Element access
char* GridBlock::operator() ( int i, int j, int k )
{
	return &data[datatype.bytes * index(i,j,k)];
}


// OPERATIONS ///////////////////////////////////////////////////////

// Returns the index to access the data array.
int GridBlock::indexBSQ ( int i, int j, int k ) const
{
	return (k-k0)*ni*nj + (i-i0)*nj + (j-j0);
}

int GridBlock::indexBIL ( int i, int j, int k ) const
{
	return (i-i0)*nj*nk + (k-k0)*nj + (j-j0);
}

int GridBlock::indexBIP ( int i, int j, int k ) const
{
	return nk * ((i-i0)*nj + (j-j0)) + (k-k0);
}

int GridBlock::index ( int i, int j, int k ) const
{
	if ( format == BSQ ) return indexBSQ(i,j,k);
	if ( format == BIL ) return indexBIL(i,j,k);
	return indexBIP(i,j,k);
}

// Returns true if (i,j,k) is within the block, false otherwise
bool GridBlock::contains ( int i, int j, int k ) const
{
	return n > 0
		&& i >= i0 && i < i0 + ni
		&& j >= j0 && j < j0 + nj
		&& k >= k0 && k < k0 + nk;
}

} // namespace bygis

#endif //YOUNG_GIS_GRIDBLOCK_20220407