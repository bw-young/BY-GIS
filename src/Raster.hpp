/////////////////////////////////////////////////////////////////////
// Raster object for GIS analysis.                                 //
// Designed to enable analysis of very large grids by loading      //
//   blocks of memory at a time rather than the entire grid.       //
//                                                                 //
// Notes:                                                          //
// + Not thread-safe.                                              //
// + Susceptible to an unlikely race condition associated with the //
//   use of tmpnam.                                                //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/02-13/2022 - Brennan Young                                   //
// - created.                                                      //
// 04/19/2022 - Brennan Young                                      //
// - added convert().                                              //
// 04/21/2022 - Brennan Young                                      //
// - added getBlockSize() and getBlockOffset().                    //
// - added X_to_j(), Y_to_i(), j_to_X(), i_to_Y(), and isNullXY(). //
// - added sub().                                                  //
// 04/26/2022 - Brennan Young                                      //
// - moved tempfilename from FileInfo object to Raster object.     //
// 04/28/2022 - Brennan Young                                      //
// - incorporated GridCell object.                                 //
// - ij_to_XY now takes x and y in the correct order.              //
// - get is now a member of Raster instead of Raster::Raster.      //
// 04/04/2022 - Brennan Young                                      //
// - added convert (str).                                          //
// 05/17/2022 - Brennan Young                                      //
// - changed BlockPattern to be offset to corner and number of     //
//   elements, instead of offset to center and radius. This should //
//   make setting the pattern much more comprehensible. Adjusted   //
//   block extent in loadBlock() accordingly.                      //
// - replaced getBlockOffset(), getBlockSize(), setBlockOffset(),  //
//   and setBlockSize() with blockOffset_row(), blockOffset_col(), //
//   blockOffset_band(), blockSize_band(), blockSize_col(),        //
//   blockSize_band(), and setBlockPattern().                      //
// 05/19/2022 - Brennan Young                                      //
// - the location where temporary files and created can now be     //
//   controlled with the static tempFileLoc() method, which sets   //
//   or reports the temppath member.                               //
// 05/23/2022 - Brennan Young                                      //
// - transitioned to using new DataType object.                    //
// - loadBlock() now better optimizes the size of the loaded block //
//   to more closely match the target data volume suggested in     //
//   the loading pattern.                                          //
// 05/25/2022 - Brennan Young                                      //
// - added getDataType().                                          //
// - corrected an issue in saveHeaderENVI() where saving a file    //
//   writes the wrong data type to the header file.                //
// 06/13/2022 - Brennan Young                                      //
// - sub(i0,i1,j0,j1) method now properly wraps, rather than       //
//   calling itself forever.                                       //
// 06/22/2022 - Brennan Young                                      //
// - renamed GeogInfo member from geog to geography.               //
// - added geog() to access and set raster geographic info.        //
// 06/23/2022 - Brennan Young                                      //
// - changed default constructor to be able to define the size of  //
//   an empty raster without geographic information.               //
// 08/18/2022 - Brennan Young                                      //
// - adjusted default loaded-block pattern from 10 rows and 100k   //
//   cols (abt. 4 MB per raster) to 3 rows and 1M cols (abt. 12 MB //
//   per raster). This should help the basic raster process 3x     //
//   more quickly with a 3x3 or smaller neighborhood.              //
// 11/14/2022 - Brennan Young                                      //
// - sub() now correctly checks for no-data in the appropriate     //
//   band.                                                         //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GIS_RASTER_20220407
#define YOUNG_GIS_RASTER_20220407

#include <cstdio>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "../BYstdlib/BinManip.hpp"
#include "FileInfo.hpp"
#include "GeogInfo.hpp"
#include "GridBlock.hpp"
#include "GridCell.hpp"

namespace bygis { // Brennan Young GIS namespace

class Raster {
private:
	struct HeaderInfo {
		std::map<std::string, std::string> fields;
	}; // HeaderInfo
	
	struct BlockPattern {
		// pattern to load a block of data relative to some (i,j,k)
		// the loaded block MUST contain (i,j,k).
		int di; // offset from a given row for the upper-left corner
		int dj; // offset from a given column for the upper-left corner
		int dk; // offset from a given band for the first band to include
		int ni; // number of rows
		int nj; // number of columns
		int nk; // number of bands
	}; // BlockPattern
	
	FileInfo file;
	HeaderInfo header;
	GeogInfo geography;
	BlockPattern pattern;
	
	static std::string temppath;
	std::string tempfilename;
	std::fstream disk;
	GridBlock block;
	
	static void fit_to_extent (int*, int*, int, int);
	int createTempFile ();
	int copyTempFile (const Raster&);
	void saveBlockBSQ ();
	void saveBlockBIL ();
	void saveBlockBIP ();
	void loadBlockBSQ ();
	void loadBlockBIL ();
	void loadBlockBIP ();
public:
	// constructors / destructor
	Raster (int i=0, int j=0, int k=0);
	Raster (const std::string&);
	Raster (const Raster&);
	Raster (const Raster&, int);
	~Raster ();
	
	// operators
	Raster& operator= (const Raster&);
	
	// operations (element access)
	double i_to_Y (int) const;
	double j_to_X (int) const;
	void ij_to_XY (int, int, double*, double*) const;
	void cell_to_XY (const GridCell&, double*, double*) const;
	int X_to_j (double) const;
	int Y_to_i (double) const;
	void XY_to_ij (double, double, int*, int*) const;
	GridCell XY_to_cell (double, double) const;
	
	int nrows () const;
	int ncols () const;
	int nbands () const;
	int size () const;
	int volume () const;
	
	double x0 () const;
	double y0 () const;
	double dx () const;
	double dy () const;
	
	const GeogInfo& geog() const;
	void geog(const GeogInfo&);
	
	char* get (int, int, int);
	char* get (int, int);
	char* get (const GridCell&);
	char* getXY (double, double, int);
	char* getXY (double, double);
	
	long long geti (int, int, int);
	long long geti (int, int);
	long long geti (const GridCell&);
	long long getiXY (double, double, int);
	long long getiXY (double, double);
	
	double getf (int, int, int);
	double getf (int, int);
	double getf (const GridCell&);
	double getfXY (double, double, int);
	double getfXY (double, double);
	
	template<typename T> void set (int, int, int, const T&);
	template<typename T> void set (int, int, const T&);
	template<typename T> void set (const GridCell&, const T&);
	template<typename T> void setXY (double, double, int, const T&);
	template<typename T> void setXY (double, double, const T&);
	void setNull (int, int, int);
	void setNull (int, int);
	void setNull (const GridCell&);
	void setNullXY (double, double, int);
	void setNullXY (double, double);
	
	template<typename T> void setNull (const T&);
	template<typename T> bool isNull (const T&) const;
	bool isNull (int, int, int);
	bool isNull (int, int);
	bool isNull (const GridCell&);
	bool isNullXY (double, double, int);
	bool isNullXY (double, double);
	long long getNulli () const;
	double getNullf () const;
	
	// operations
	bool isValid () const;
	bool contains (int, int) const;
	bool contains (int, int, int) const;
	bool contains (const GridCell&) const;
	bool containsXY (double, double) const;
	bool isLoaded (int, int) const;
	bool isLoaded (int, int, int) const;
	bool isLoaded (const GridCell&) const;
	bool isLoadedXY (double, double) const;
	
	int blockSize_row() const;
	int blockSize_col() const;
	int blockSize_band() const;
	void blockSize_row(int);
	void blockSize_col(int);
	void blockSize_band(int);
	int blockOffset_row() const;
	int blockOffset_col() const;
	int blockOffset_band() const;
	void blockOffset_row(int);
	void blockOffset_col(int);
	void blockOffset_band(int);
	void setBlockPattern(int, int, int, int, int, int);
	
	Raster sub (int, int, int, int, int, int);
	Raster sub (int, int, int, int);
	Raster sub (int, int);
	Raster sub (int, int, int, int, const std::vector<int>&);
	Raster sub (const std::vector<int>&);
	
	const DataType& getDataType () const;
	Raster convert (const DataType&, int);
	Raster convert (const DataType&);
	Raster convert (unsigned char, int);
	Raster convert (unsigned char);
	Raster convert (const std::string&, int);
	Raster convert (const std::string&);
	
	// I/O
	static std::string tempFileLoc ();
	static void tempFileLoc (const std::string&);
	
	int saveHeaderENVI (const std::string&);
	int saveFileENVI (const std::string&);
	int saveFileENVI ();
	int commit ();
	
	int loadHeaderENVI (const std::string&);
	int loadFileENVI (const std::string&);
	void unloadBlock ();
	void loadBlock (int, int, int, const BlockPattern&);
	char* request (int, int, int);
	
	// debug
	std::string getTempFilename () const;
}; // Raster

std::string Raster::temppath = "./";


// CONSTRUCTORS / DESTRUCTORS ///////////////////////////////////////

Raster::Raster ( int ni, int nj, int nk )
{
	pattern.di = -1;
	pattern.dj = -1;
	pattern.dk = 0;
	pattern.ni = 3;
	pattern.nj = 1000000;
	pattern.nk = 1;
	
	file.nr = ni;
	file.nc = nj;
	file.nb = nk;
	file.datatype = DataType(4);
	setNull(-9999);
	
	if ( createTempFile() != 0 ) {
		std::cout << "ERR: creating temporary file:"
			<< " Raster constructor\n";
		exit(1);
	}
	
	for ( int i = 0; i < ni; ++i ) {
		for ( int j = 0; j < nj; ++j ) {
			for ( int k = 0; k < nk; ++k ) setNull(i,j,k);
		}
	}
}

Raster::Raster ( const std::string& fn )
{
	pattern.di = -1;
	pattern.dj = -1;
	pattern.dk = 0;
	pattern.ni = 3;
	pattern.nj = 1000000;
	pattern.nk = 1;
	
	if ( createTempFile() != 0 ) {
		std::cout << "ERR: creating temporary file:"
			<< " Raster constructor\n";
		exit(1);
	}
	
	if ( fn.size() > 0 ) loadFileENVI(fn);
}

Raster::Raster ( const Raster& grid )
{
	file = grid.file;
	header = grid.header;
	geography = grid.geography;
	pattern = grid.pattern;
	
	if ( createTempFile() != 0 ) {
		std::cout << "ERR: creating temporary file:"
			<< " Raster copy constructor\n";
		exit(1);
	}
	
	if ( copyTempFile(grid) != 0 ) {
		std::cout << "ERR: copying temporary file:"
			<< " Raster copy constructor\n";
	}
}

// Copy raster dimensions and geography, with k bands.
Raster::Raster ( const Raster& grid, int k )
{
	file = grid.file;
	file.nb = k < 1 ? 1 : k; // always at least one band
	header = grid.header;
	geography = grid.geography;
	pattern = grid.pattern;
	
	if ( createTempFile() != 0 ) {
		std::cout << "ERR: creating temporary file:"
			<< " Raster constructor\n";
		exit(1);
	}
	
	// build out temporary file to the grid's byte size.
	for ( int i = 0; i < file.nr; ++i ) {
		for ( int j = 0; j < file.nc; ++j ) {
			for ( int k = 0; k < file.nb; ++k )
				disk.write(file.nulval, file.datatype.bytes);
		}
	}
}

// Destructor
Raster::~Raster ()
{
	// close files and delete temporary file
	disk.close();
	remove(tempfilename.c_str());
}


// OPERATORS ////////////////////////////////////////////////////////

Raster& Raster::operator= ( const Raster& grid )
{
	if ( &grid == this ) return *this;
	
	file = grid.file;
	header = grid.header;
	geography = grid.geography;
	pattern = grid.pattern;
	
	if ( !disk || !disk.is_open() ) {
		if ( createTempFile() != 0 ) {
			std::cout << "ERR: creating temporary file:"
				<< " Raster::operator=\n";
			exit(1);
		}
	}
	
	if ( copyTempFile(grid) != 0 ) {
		std::cout << "ERR: copying temporary file:"
			<< " Raster::operator=\n";
	}
	
	return *this;
}


// OPERATIONS ///////////////////////////////////////////////////////

// Fits some window or block extent (i0, i1) to the given constraints
// (a, b).
void Raster::fit_to_extent ( int* i0, int* i1, int a, int b )
{
	if ( *i0 < a ) {
		int d = a - *i0;
		*i0 += d;
		*i1 += d;
		if ( *i1 >= b ) *i1 = b;
	}
	else if ( *i1 > b ) {
		int d = b - *i1;
		*i0 += d;
		*i1 += d;
		if ( *i0 < a ) *i0 = a;
	}
}

// Convert (row,column) to (x,y) at the center of the grid cell.
double Raster::i_to_Y ( int i ) const
{
	return geography.y0 - (i+0.5) * geography.dy;
}

double Raster::j_to_X ( int j ) const
{
	return geography.x0 + (j+0.5) * geography.dx;
}

void Raster::ij_to_XY ( int i, int j, double* x, double* y ) const
{
	*x = j_to_X(j);
	*y = i_to_Y(i);
}

void Raster::cell_to_XY ( const GridCell& cell, double* x, double* y ) const
{
	ij_to_XY(cell.i, cell.j, x, y);
}

// Convert (x,y) to (row,column) for the containing grid cell.
int Raster::X_to_j ( double x ) const
{
	return floor((x - geography.x0) / geography.dx);
}

int Raster::Y_to_i ( double y ) const
{
	return floor((geography.y0 - y) / geography.dy);
}

void Raster::XY_to_ij ( double x, double y, int* i, int* j ) const
{
	*i = Y_to_i(y);
	*j = X_to_j(x);
}

GridCell Raster::XY_to_cell ( double x, double y ) const
{
	GridCell cell;
	XY_to_ij(x, y, &cell.i, &cell.j);
	return cell;
}

// Get number of lines or rows in file.
int Raster::nrows () const
{
	return file.nr;
}

// Get number of columns in file.
int Raster::ncols () const
{
	return file.nc;
}

// Get number of bands in file.
int Raster::nbands () const
{
	return file.nb;
}

// Get size of each band in file.
int Raster::size () const
{
	return file.nr * file.nc;
}

// Get the total number of elements in file.
int Raster::volume () const
{
	return file.nr * file.nc * file.nb;
}

// Get the x-coordinate of the upper-left corner of (0,0).
double Raster::x0 () const
{
	return geography.x0;
}

// Get the y-coordinate of the upper-left corner of (0,0).
double Raster::y0 () const
{
	return geography.y0;
}

// Get the spatial resolution in the x-dimension.
double Raster::dx () const
{
	return geography.dx;
}

// Get the spatial resolution in the y-dimension.
double Raster::dy () const
{
	return geography.dy;
}

// Get the grid's spatial information.
const GeogInfo& Raster::geog () const
{
	return geography;
}

// Set the grid's spatial information.
void Raster::geog ( const GeogInfo& ginfo )
{
	geography = ginfo;
}

// Get the value of cell (i,j,k).
char* Raster::get ( int i, int j, int k )
{
	return request(i, j, k);
}

char* Raster::get ( int i, int j )
{
	return get(i, j, 0);
}

char* Raster::get ( const GridCell& cell )
{
	return get(cell.i, cell.j, cell.k);
}

char* Raster::getXY ( double x, double y, int k )
{
	int i, j;
	XY_to_ij(x, y, &i, &j);
	return get(i, j, k);
}

char* Raster::getXY ( double x, double y )
{
	return getXY(x, y, 0);
}

long long Raster::geti ( int i, int j, int k )
{
	char* v = request(i, j, k);
	if ( !v ) return getNulli();
	if ( file.datatype == DataType::INT1 )
		return *((char*) v);
	if ( file.datatype == DataType::INT2 )
		return *((short*) v);
	if ( file.datatype == DataType::INT4 )
		return *((long*) v);
	if ( file.datatype == DataType::INT8 )
		return *((long long*) v);
	if ( file.datatype == DataType::UINT2 )
		return *((unsigned short*) v);
	if ( file.datatype == DataType::UINT4 )
		return *((unsigned long*) v);
	if ( file.datatype == DataType::UINT8 )
		return *((unsigned long long*) v);
	if ( file.datatype == DataType::FLOAT4 )
		return *((float*) v);
	if ( file.datatype == DataType::FLOAT8 )
		return *((double*) v);
	return 0;
}

long long Raster::geti ( int i, int j )
{
	return geti(i, j, 0);
}

long long Raster::geti ( const GridCell& cell )
{
	return geti(cell.i, cell.j, cell.k);
}

long long Raster::getiXY ( double x, double y, int k )
{
	int i, j;
	XY_to_ij(x, y, &i, &j);
	return geti(i, j, k);
}

long long Raster::getiXY ( double x, double y )
{
	return getiXY(x, y, 0);
}

double Raster::getf ( int i, int j, int k )
{
	char* v = request(i, j, k);
	if ( !v ) return getNullf();
	if ( file.datatype == DataType::INT1 )
		return *((char*) v);
	if ( file.datatype == DataType::INT2 )
		return *((short*) v);
	if ( file.datatype == DataType::INT4 )
		return *((long*) v);
	if ( file.datatype == DataType::INT8 )
		return *((long long*) v);
	if ( file.datatype == DataType::UINT2 )
		return *((unsigned short*) v);
	if ( file.datatype == DataType::UINT4 )
		return *((unsigned long*) v);
	if ( file.datatype == DataType::UINT8 )
		return *((unsigned long long*) v);
	if ( file.datatype == DataType::FLOAT4 )
		return *((float*) v);
	if ( file.datatype == DataType::FLOAT8 )
		return *((double*) v);
	return 0;
}

double Raster::getf ( int i, int j )
{
	return getf(i, j, 0);
}

double Raster::getf ( const GridCell& cell )
{
	return getf(cell.i, cell.j, cell.k);
}

double Raster::getfXY ( double x, double y, int k )
{
	int i, j;
	XY_to_ij(x, y, &i, &j);
	return getf(i, j, k);
}

double Raster::getfXY ( double x, double y )
{
	return getfXY(x, y, 0);
}

// Set value at (i,j,k). If (i,j,k) is not valid, does nothing.
template<typename T>
void Raster::set ( int i, int j, int k, const T& x )
{
	block.update = true;
	char* v = request(i, j, k);
	if ( !v ) return;
	if ( file.datatype == DataType::INT1 )
		*((char*) v) = (char) x;
	if ( file.datatype == DataType::INT2 )
		*((short*) v) = (short) x;
	if ( file.datatype == DataType::INT4 )
		*((long*) v) = (long) x;
	if ( file.datatype == DataType::INT8 )
		*((long long*) v) = (long long) x;
	if ( file.datatype == DataType::UINT2 )
		*((unsigned short*) v) = (unsigned short) x;
	if ( file.datatype == DataType::UINT4 )
		*((unsigned long*) v) = (unsigned long) x;
	if ( file.datatype == DataType::UINT8 )
		*((unsigned long long*) v) = (unsigned long long) x;
	if ( file.datatype == DataType::FLOAT4 )
		*((float*) v) = (float) x;
	if ( file.datatype == DataType::FLOAT8 )
		*((double*) v) = (double) x;
}

template<typename T>
void Raster::set ( int i, int j, const T& v )
{
	set(i, j, 0, v);
}

template <typename T>
void Raster::set ( const GridCell& cell, const T& v )
{
	set(cell.i, cell.j, cell.k, v);
}

template<typename T>
void Raster::setXY ( double x, double y, int k, const T& v )
{
	set(Y_to_i(y), X_to_j(x), k, v);
}

template<typename T>
void Raster::setXY ( double x, double y, const T& v )
{
	setXY(x, y, 0, v);
}

// Set the given cell to the no-data value. If (i,j,k) is not valid,
// does nothing.
void Raster::setNull ( int i, int j, int k )
{
	block.update = true;
	char* v = request(i, j, k);
	if ( !v ) return;
	if ( file.datatype == DataType::INT1 )
		*((char*) v) = *((char*) file.nulval);
	if ( file.datatype == DataType::INT2 )
		*((short*) v) = *((short*) file.nulval);
	if ( file.datatype == DataType::INT4 )
		*((long*) v) = *((long*) file.nulval);
	if ( file.datatype == DataType::INT8 )
		*((long long*) v) = *((long long*) file.nulval);
	if ( file.datatype == DataType::UINT2 )
		*((unsigned short*) v) = *((unsigned short*) file.nulval);
	if ( file.datatype == DataType::UINT4 )
		*((unsigned long*) v) = *((unsigned long*) file.nulval);
	if ( file.datatype == DataType::UINT8 ) {
		*((unsigned long long*) v) =
			*((unsigned long long*) file.nulval);
	}
	if ( file.datatype == DataType::FLOAT4 )
		*((float*) v) = *((float*) file.nulval);
	if ( file.datatype == DataType::FLOAT8 )
		*((double*) v) = *((double*) file.nulval);
}

void Raster::setNull ( int i, int j )
{
	setNull(i, j, 0);
}

void Raster::setNull ( const GridCell& cell )
{
	setNull(cell.i, cell.j, cell.k);
}

void Raster::setNullXY ( double x, double y, int k )
{
	setNull(Y_to_i(y), X_to_j(x), k);
}

void Raster::setNullXY ( double x, double y )
{
	setNullXY(x, y, 0);
}

// Set the no-data value.
template<typename T>
void Raster::setNull ( const T& x )
{
	delete [] file.nulval;
	try {
		file.nulval = new char [file.datatype.bytes];
	}
	catch (...)
	{
		std::cout << "ERR: memory allocation: Raster::setNull\n";
		exit(0);
	}
	
	if ( file.datatype == DataType::INT1 )
		*((char*) file.nulval) = x;
	else if ( file.datatype == DataType::INT2 )
		*((short*) file.nulval) = x;
	else if ( file.datatype == DataType::INT4 )
		*((long*) file.nulval) = x;
	else if ( file.datatype == DataType::INT8 )
		*((long long*) file.nulval) = x;
	else if ( file.datatype == DataType::UINT2 )
		*((unsigned short*) file.nulval) = x;
	else if ( file.datatype == DataType::UINT4 )
		*((unsigned long*) file.nulval) = x;
	else if ( file.datatype == DataType::UINT8 )
		*((unsigned long long*) file.nulval) = x;
	else if ( file.datatype == DataType::FLOAT4 )
		*((float*) file.nulval) = x;
	else if ( file.datatype == DataType::FLOAT8 )
		*((double*) file.nulval) = x;
}

// Return true if the given value is the no-data value or invalid.
// Always returns true with unknown data types.
template<typename T>
bool Raster::isNull ( const T& x ) const
{
	if ( file.datatype == DataType::INT1 )
		return *((char*) file.nulval) == (char) x;
	if ( file.datatype == DataType::INT2 )
		return *((short*) file.nulval) == (short) x;
	if ( file.datatype == DataType::INT4 )
		return *((long*) file.nulval) == (long) x;
	if ( file.datatype == DataType::INT8 )
		return *((long long*) file.nulval) == (long long) x;
	if ( file.datatype == DataType::UINT2 )
		return *((unsigned short*) file.nulval) == (unsigned short) x;
	if ( file.datatype == DataType::UINT4 )
		return *((unsigned long*) file.nulval) == (unsigned long) x;
	if ( file.datatype == DataType::UINT8 ) {
		return *((unsigned long long*) file.nulval)
			== (unsigned long long) x;
	}
	if ( file.datatype == DataType::FLOAT4 )
		return *((float*) file.nulval) == (float) x;
	if ( file.datatype == DataType::FLOAT8 )
		return *((double*) file.nulval) == (double) x;
	return true;
}

bool Raster::isNull ( int i, int j, int k )
{
	char* x = request(i, j, k);
	if ( !x ) return true;
	if ( file.datatype == DataType::INT1 )
		return *((char*) file.nulval) == (char) *((char*) x);
	if ( file.datatype == DataType::INT2 )
		return *((short*) file.nulval) == (short) *((short*) x);
	if ( file.datatype == DataType::INT4 )
		return *((long*) file.nulval) == (long) *((long*) x);
	if ( file.datatype == DataType::INT8 ) {
		return *((long long*) file.nulval)
			== (long long) *((long long*) x);
	}
	if ( file.datatype == DataType::UINT2 ) {
		return *((unsigned short*) file.nulval)
			== (unsigned short) *((unsigned short*) x);
	}
	if ( file.datatype == DataType::UINT4 ) {
		return *((unsigned long*) file.nulval)
			== (unsigned long) *((unsigned long*) x);
	}
	if ( file.datatype == DataType::UINT8 ) {
		return *((unsigned long long*) file.nulval)
			== (unsigned long long) *((unsigned long long*) x);
	}
	if ( file.datatype == DataType::FLOAT4 )
		return *((float*) file.nulval) == (float) *((float*) x);
	if ( file.datatype == DataType::FLOAT8 )
		return *((double*) file.nulval) == (double) *((double*) x);
	return true;
}

bool Raster::isNull ( int i, int j )
{
	return isNull(i, j, 0);
}

bool Raster::isNull ( const GridCell& cell )
{
	return isNull(cell.i, cell.j, cell.k);
}

bool Raster::isNullXY ( double x, double y, int k )
{
	return isNull(Y_to_i(y), X_to_j(x), k);
}

bool Raster::isNullXY ( double x, double y )
{
	return isNullXY(x, y, 0);
}

long long Raster::getNulli () const
{
	if ( file.datatype == DataType::INT1 )
		return *((char*) file.nulval);
	if ( file.datatype == DataType::INT2 )
		return *((short*) file.nulval);
	if ( file.datatype == DataType::INT4 )
		return *((long*) file.nulval);
	if ( file.datatype == DataType::INT8 )
		return *((long long*) file.nulval);
	if ( file.datatype == DataType::UINT2 )
		return *((unsigned short*) file.nulval);
	if ( file.datatype == DataType::UINT4 )
		return *((unsigned long*) file.nulval);
	if ( file.datatype == DataType::UINT8 )
		return *((unsigned long long*) file.nulval);
	if ( file.datatype == DataType::FLOAT4 )
		return *((float*) file.nulval);
	if ( file.datatype == DataType::FLOAT8 )
		return *((double*) file.nulval);
	return 0;
}

double Raster::getNullf () const
{
	if ( file.datatype == DataType::INT1 )
		return *((char*) file.nulval);
	if ( file.datatype == DataType::INT2 )
		return *((short*) file.nulval);
	if ( file.datatype == DataType::INT4 )
		return *((long*) file.nulval);
	if ( file.datatype == DataType::INT8 )
		return *((long long*) file.nulval);
	if ( file.datatype == DataType::UINT2 )
		return *((unsigned short*) file.nulval);
	if ( file.datatype == DataType::UINT4 )
		return *((unsigned long*) file.nulval);
	if ( file.datatype == DataType::UINT8 )
		return *((unsigned long long*) file.nulval);
	if ( file.datatype == DataType::FLOAT4 )
		return *((float*) file.nulval);
	if ( file.datatype == DataType::FLOAT8 )
		return *((double*) file.nulval);
	return 0.0;
}

// Return true if the raster is in a valid (read/writable) state.
bool Raster::isValid () const
{
	return disk.is_open();
}

// Return true if the raster contains the given grid cell.
bool Raster::contains ( int i, int j, int k ) const
{
	return i >= 0 && i < file.nr
		&& j >= 0 && j < file.nc
		&& k >= 0 && k < file.nb;
}

bool Raster::contains ( int i, int j ) const
{
	return contains(i,j,0);
}

bool Raster::contains ( const GridCell& cell ) const
{
	return contains(cell.i, cell.j, cell.k);
}

// Return true if the raster contains the given (x,y) coordinates.
bool Raster::containsXY ( double x, double y ) const
{
	return contains(Y_to_i(y), X_to_j(x));
}

// Return true if the given grid cell is loaded into memory.
bool Raster::isLoaded ( int i, int j, int k ) const
{
	return block.contains(i,j,k);
}

bool Raster::isLoaded ( int i, int j ) const
{
	return isLoaded(i, j, block.k0);
}

bool Raster::isLoaded ( const GridCell& cell ) const
{
	return isLoaded(cell.i, cell.j, cell.k);
}

// Return true if the given (x,y) coordinates are loaded into memory;
// false otherwise.
bool Raster::isLoadedXY ( double x, double y ) const
{
	return isLoaded(Y_to_i(y), X_to_j(x));
}

// Get/set the size of the block pattern for loading.
int Raster::blockSize_row () const
{
	return pattern.ni;
}

int Raster::blockSize_col () const
{
	return pattern.nj;
}

int Raster::blockSize_band () const
{
	return pattern.nk;
}

void Raster::blockSize_row ( int x )
{
	pattern.ni = x < 1 ? 1 : x;
}

void Raster::blockSize_col ( int x )
{
	pattern.nj = x < 1 ? 1 : x;
}

void Raster::blockSize_band ( int x )
{
	pattern.nk = x < 1 ? 1 : x;
}

// Get/set the offet of the block pattern for loading.
int Raster::blockOffset_row () const
{
	return pattern.di;
}

int Raster::blockOffset_col () const
{
	return pattern.dj;
}

int Raster::blockOffset_band () const
{
	return pattern.dk;
}

void Raster::blockOffset_row ( int x )
{
	if ( x > 0 ) x = 0;
	if ( x + pattern.ni <= 0 ) x = 1 - pattern.ni;
	pattern.di = x;
}

void Raster::blockOffset_col ( int x )
{
	if ( x > 0 ) x = 0;
	if ( x + pattern.nj <= 0 ) x = 1 - pattern.nj;
	pattern.dj = x;
}

void Raster::blockOffset_band ( int x )
{
	if ( x > 0 ) x = 0;
	if ( x + pattern.nk <= 0 ) x = 1 - pattern.nk;
	pattern.dk = x;
}

// Set the size of the blocks that are loaded from memory.
// di,dj,dk must be <= 0, ni,nj,nk must be >= 1. If the box does
// not contain the target point (0,0), adjusts offset to fit.
void Raster::setBlockPattern (
	int di, int dj, int dk, int ni, int nj, int nk )
{
	blockSize_row(ni);
	blockSize_col(nj);
	blockSize_band(nk);
	blockOffset_row(di);
	blockOffset_col(dj);
	blockOffset_band(dk);
}

// Create a raster that is a subset of the this raster.
Raster Raster::sub ( int i0, int i1, int j0, int j1,
	const std::vector<int>& bands )
{
	std::vector<int> bnd;
	if ( bands.size() > 0 ) bnd.reserve(bands.size());
	for ( size_t k = 0; k < bands.size(); ++k ) {
		if ( bands[k] >= 0 && bands[k] < file.nb )
			bnd.push_back(bands[k]);
	}
	
	if ( bnd.size() == 0 ) {
		bnd.push_back(0);
		return sub(i0, i1, j0, j1, bnd);
	}
	
	if ( i0 < 0 ) i0 = 0;
	if ( i1 > file.nr ) i1 = file.nr;
	if ( i1 == i0 ) ++i1;
	if ( j0 < 0 ) j0 = 0;
	if ( j1 > file.nc ) j1 = file.nc;
	if ( j1 == j0 ) ++j1;
	
	Raster out;
	
	out.file = file;
	out.file.filename = "";
	out.file.nr = i1-i0;
	out.file.nc = j1-j0;
	out.file.nb = bands.size();
	out.file.bandnames = std::vector<std::string> ();
	out.file.bandnames.reserve(bnd.size());
	for ( size_t k = 0; k < bnd.size(); ++k )
		out.file.bandnames.push_back(file.bandnames[bnd[k]]);
	
	out.geography = geography;
	out.geography.x0 = j_to_X(j0) - 0.5 * geography.dx;
	out.geography.y0 = i_to_Y(i0) + 0.5 * geography.dy;
	out.pattern = pattern;
	
	out.disk.clear();
	out.disk.seekg(0);
	
	int j;
	size_t k;
	for ( int i = i0; i < i1; ++i ) {
		for ( j = j0; j < j1; ++j ) {
			for ( k = 0; k < bnd.size(); ++k )
				out.disk.write(out.file.nulval, out.file.datatype.bytes);
		}
	}
	
	char* x;
	for ( int i = i0; i < i1; ++i ) {
		for ( j = j0; j < j1; ++j ) {
			for ( k = 0; k < bnd.size(); ++k ) {
				if ( isNull(i,j,bnd[k]) ) continue;
				x = get(i,j,bnd[k]);
				if ( file.datatype == DataType::INT1 )
					out.set(i-i0,j-j0,k, *((char*) x));
				if ( file.datatype == DataType::INT2 )
					out.set(i-i0,j-j0,k, *((short*) x));
				if ( file.datatype == DataType::INT4 )
					out.set(i-i0,j-j0,k, *((long*) x));
				if ( file.datatype == DataType::INT8 )
					out.set(i-i0,j-j0,k, *((long long*) x));
				if ( file.datatype == DataType::UINT2 )
					out.set(i-i0,j-j0,k, *((unsigned short*) x));
				if ( file.datatype == DataType::UINT4 )
					out.set(i-i0,j-j0,k, *((unsigned long*) x));
				if ( file.datatype == DataType::UINT8 )
					out.set(i-i0,j-j0,k, *((unsigned long long*) x));
				if ( file.datatype == DataType::FLOAT4 )
					out.set(i-i0,j-j0,k, *((float*) x));
				if ( file.datatype == DataType::FLOAT8 )
					out.set(i-i0,j-j0,k, *((double*) x));
			}
		}
	}
	
	return out;
}

Raster Raster::sub ( const std::vector<int>& bands )
{
	return sub(0, file.nr, 0, file.nc, bands);
}

Raster Raster::sub ( int i0, int i1, int j0, int j1, int k0, int k1 )
{
	if ( k0 < 0 ) k0 = 0;
	if ( k1 > file.nb ) k1 = file.nb;
	if ( k1 == k0 ) ++k1;
	std::vector<int> bands;
	bands.reserve(k1-k0);
	for ( int k = k0; k < k1; ++k ) bands.push_back(k);
	return sub(i0, i1, j0, j1, bands);
}

Raster Raster::sub ( int i0, int i1, int j0, int j1 )
{
	return sub(i0, i1, j0, j1, 0, file.nb);
}

Raster Raster::sub ( int k0, int k1 )
{
	return sub(0, file.nr, 0, file.nc, k0, k1);
}

// Return the raster's data type.
const DataType& Raster::getDataType () const
{
	return file.datatype;
}

// Convert the grid to a new data type. If numBands is not specified
// or is negative, translates this raster's data to the new one.
// Otherwise, creates an empy raster of the same number of rows and
// columns but a different number of bands (minimum 1).
Raster Raster::convert ( const DataType& datatype, int numBands )
{
	Raster out; // (creates empty temporary file)
	
	// copy header and file organization information
	out.file = file;
	out.header = header;
	out.geography = geography;
	out.pattern = pattern;
	
	// set up for conversion
	out.file.filename = "";
	if ( numBands < 0 ) out.file.nb = file.nb;
	else if ( numBands < 1 ) out.file.nb = 1;
	else out.file.nb = numBands;
	out.file.datatype = datatype;
	out.setNull(file.datatype.isInt() ? getNulli() : getNullf());
	
	// initialize temporary file
	out.disk.clear();
	out.disk.seekg(0);
	
	int j, k;
	for ( int i = 0; i < out.file.nr; ++i ) {
		for ( j = 0; j < out.file.nc; ++j ) {
			for ( k = 0; k < out.file.nb; ++k )
				out.disk.write(out.file.nulval, out.file.datatype.bytes);
		}
	}
	
	// copy data and translate to new data type
	if ( numBands < 0 ) {
		out.disk.clear();
		out.disk.seekg(0);
		
		std::fstream tmp (tempfilename,
			std::fstream::in | std::fstream::binary);
		if ( !tmp || tmp.eof() ) {
			std::cout << "ERR: copying temporary file:"
				<< " Raster::convert\n";
			exit(1);
		}
		
		char* x;
		for ( int i = 0; i < file.nr; ++i ) {
			for ( j = 0; j < file.nc; ++j ) {
				for ( k = 0; k < file.nb; ++k ) {
					if ( isNull(i,j) ) continue;
					x = get(i,j,k);
					if ( file.datatype == DataType::INT1 )
						out.set(i,j,k, *((char*) x));
					if ( file.datatype == DataType::INT2 )
						out.set(i,j,k, *((short*) x));
					if ( file.datatype == DataType::INT4 )
						out.set(i,j,k, *((long*) x));
					if ( file.datatype == DataType::INT8 )
						out.set(i,j,k, *((long long*) x));
					if ( file.datatype == DataType::UINT2 )
						out.set(i,j,k, *((unsigned short*) x));
					if ( file.datatype == DataType::UINT4 )
						out.set(i,j,k, *((unsigned long*) x));
					if ( file.datatype == DataType::UINT8 )
						out.set(i,j,k, *((unsigned long long*) x));
					if ( file.datatype == DataType::FLOAT4 )
						out.set(i,j,k, *((float*) x));
					if ( file.datatype == DataType::FLOAT8 )
						out.set(i,j,k, *((double*) x));
				}
			}
		}
	}
	
	return out;
}

Raster Raster::convert ( const DataType& datatype )
{
	return convert(datatype, -1);
}

Raster Raster::convert ( unsigned char datatype, int numBands )
{
	return convert(DataType(datatype), numBands);
}

Raster Raster::convert ( unsigned char datatype )
{
	return convert(DataType(datatype), -1);
}

Raster Raster::convert ( const std::string& datatype, int numBands )
{
	return convert(DataType(datatype), numBands);
}

Raster Raster::convert ( const std::string& datatype )
{
	return convert(DataType(datatype), -1);
}


// I/O //////////////////////////////////////////////////////////////

// Set the static temporary file location for the program's future
// temporary files.
std::string Raster::tempFileLoc ()
{
	return temppath;
}

void Raster::tempFileLoc ( const std::string& s )
{
	if ( s.size() == 0 ) temppath = "./";
	else if ( s.back() != '/' ) temppath = s + "/";
	else temppath = s;
}

// Create a temporary file for this raster.
int Raster::createTempFile ()
{
	// name for temporary file
	tempfilename = std::tmpnam(nullptr);
	tempfilename = temppath + "_tmp.rast."
		+ tempfilename.substr(1, tempfilename.size()-1);
	
	// create file
	disk = std::fstream (tempfilename, std::fstream::out);
	if ( !disk || disk.eof() ) return 1;
	disk.close();
	
	// open file for processing
	disk = std::fstream (tempfilename,
		std::fstream::in | std::fstream::out | std::fstream::binary);
	if ( !disk || disk.eof() ) return 1;
	
	return 0;
}

// Copy the raster's temporary file to this raster's temporary file.
int Raster::copyTempFile ( const Raster& grid )
{
	file = grid.file;
	block = GridBlock();
	
	std::fstream tmp (grid.tempfilename,
		std::fstream::in | std::fstream::binary);
	if ( !tmp || tmp.eof() ) return 1;
	
	disk.clear();
	disk.seekg(0);
	disk << tmp.rdbuf();
	
	if ( grid.block.update ) {
		block = grid.block;
		commit();
	}
	
	return 0;
}

// Reads a block of data from the file.
void Raster::saveBlockBSQ ()
{
	int i, n = 0;
	
	// return to the beginning of the file
	disk.clear();
	disk.seekg(0);
	
	// number of bytes to ignore before/after block's extent
	int nbi_before = block.i0 * file.nc * file.datatype.bytes;
	int nbi_after = (file.nr - (block.i0 + block.ni))
		* file.nc * file.datatype.bytes;
	int nbj_before = block.j0 * file.datatype.bytes;
	int nbj_after = (file.nc - (block.j0 + block.nj))
		* file.datatype.bytes;
	
	// number of bytes to write at a time
	int nb_write = block.nj * file.datatype.bytes;
	
	// skip bands before k0
	disk.ignore(block.k0 * file.nr * file.nc * file.datatype.bytes);
	
	for ( int k = 0; k < block.nk; ++k ) {
		// skip rows before i0
		disk.ignore(nbi_before);
		
		for ( i = 0; i < block.ni; ++i ) {
			// skip columns before j0
			disk.ignore(nbj_before);
			
			// write row
			disk.write(&block.data[n * file.datatype.bytes], nb_write);
			n += block.nj;
			
			// skip columns after j1
			disk.ignore(nbj_after);
		}
		
		// skip rows after i1
		disk.ignore(nbi_after);
	}
}

// Reads a block of data from the file.
void Raster::saveBlockBIL ()
{
	int k, n = 0;
	
	// return to the beginning of the file
	disk.clear();
	disk.seekg(0);
	
	// number of bytes to ignore before/after block's extent
	int nbj_before = block.j0 * file.datatype.bytes;
	int nbj_after = (file.nc - (block.j0 + block.nj))
		* file.datatype.bytes;
	int nbk_before = block.k0 * file.nc * file.datatype.bytes;
	int nbk_after = (file.nb - (block.k0 + block.nk))
		* file.nc * file.datatype.bytes;
	
	// number of bytes to write at a time
	int nb_write = block.nj * file.datatype.bytes;
	
	// skip rows before i0
	disk.ignore(block.i0 * file.nc * file.nb * file.datatype.bytes);
	
	for ( int i = 0; i < block.ni; ++i ) {
		// skip bands before k0
		disk.ignore(nbk_before);
		
		for ( k = 0; k < block.nk; ++k ) {
			// skip columns before j0
			disk.ignore(nbj_before);
			
			// write row
			disk.write(&block.data[n * file.datatype.bytes], nb_write);
			n += block.nj;
			
			// skip columns after j1
			disk.ignore(nbj_after);
		}
		
		// skip bands after k1
		disk.ignore(nbk_after);
	}
}

// Reads a block of data from the file.
void Raster::saveBlockBIP ()
{
	int j, n = 0;
	
	// return to the beginning of the file
	disk.clear();
	disk.seekg(0);
	
	// number of bytes to ignore before/after block's extent
	int nbj_before = block.j0 * file.nb * file.datatype.bytes;
	int nbj_after = (file.nc - (block.j0 + block.nj))
		* file.nb * file.datatype.bytes;
	int nbk_before = block.k0 * file.datatype.bytes;
	int nbk_after = (file.nb - (block.k0 + block.nk))
		* file.datatype.bytes;
	
	// number of bytes to write at a time
	int nb_write = block.nk * file.datatype.bytes;
	
	// skip rows before i0
	disk.ignore(block.i0 * file.nc * file.nb * file.datatype.bytes);
	
	for ( int i = 0; i < block.ni; ++i ) {
		// skip columns before j0
		disk.ignore(nbj_before);
		
		for ( j = 0; j < block.nj; ++j ) {
			// skip bands before k0
			disk.ignore(nbk_before);
			
			// write bands
			disk.write(&block.data[n * file.datatype.bytes], nb_write);
			n += block.nk;
			
			// skip bands after k1
			disk.ignore(nbk_after);
		}
		
		// skip columns after j1
		disk.ignore(nbj_after);
	}
}

// Save header information.
int Raster::saveHeaderENVI ( const std::string& fn )
{
	std::fstream tmp (fn, std::fstream::out);
	if ( !tmp || tmp.eof() ) return 1;
	
	tmp << "ENVI"
		<< "\nsamples = " << file.nc
		<< "\nlines   = " << file.nr
		<< "\nbands   = " << file.nb
		<< "\nheader offset = 0"
		<< "\nfile type = ENVI Standard"
		<< "\ndata type = " << ((int) file.datatype.code)
		<< "\ninterleave = ";
	if ( file.interleave == GridBlock::BIP ) tmp << "bip";
	else if ( file.interleave == GridBlock::BIL ) tmp << "bil";
	else tmp << "bsq";
	tmp << "\nbyte order = " << (bystd::osLittleEndian() ? "0" : "1");
	tmp << "\ndata ignore value = ";
	if ( file.datatype == DataType::INT1
			|| file.datatype == DataType::INT2
			|| file.datatype == DataType::INT4
			|| file.datatype == DataType::INT8
			|| file.datatype == DataType::UINT2
			|| file.datatype == DataType::UINT4
			|| file.datatype == DataType::UINT8 )
		tmp << getNulli();
	else tmp << bystd::num_to_str(getNullf());
	tmp << "\nmap info = {" << geography.projectionName
		<< ", " << geography.tiePoint_x
		<< ", " << geography.tiePoint_y
		<< ", " << bystd::num_to_str(geography.x0)
		<< ", " << bystd::num_to_str(geography.y0)
		<< ", " << bystd::num_to_str(geography.dx)
		<< ", " << bystd::num_to_str(geography.dy)
		<< ", " << geography.projectionZone
		<< ", " << geography.zoneCode
		<< ", " << geography.datum
		<< ", " << geography.units
		<< "}";
	tmp << "\ncoordinate system string = {"
		<< geography.coordinateSys << "}";
	tmp << "\nband names = {";
	for ( size_t i = 0; i < file.bandnames.size(); ++i ) {
		if ( i > 0 ) tmp << ", ";
		if ( file.bandnames[i].size() == 0 ) tmp << "Band " << (i+1);
		else tmp << file.bandnames[i];
	}
	tmp << "}";
	
	tmp.close();
	
	return 0;
}

// Save header information, save to temporary file, remove file with
// filename if it exists, and rename temporary file to given
// filename.
int Raster::saveFileENVI ( const std::string& fn )
{
	if ( file.filename.size() == 0 ) file.filename = fn;
	
	if ( block.update ) commit();
	
	disk.clear();
	disk.seekg(0);
	
	if ( !disk || disk.eof() ) return 1;
	
	// save binary file
	std::fstream tmp (fn,
		std::fstream::out | std::fstream::binary);
	if ( !tmp || tmp.eof() ) return 1;
	
	// copy from temporary file
	// (rename in stdio does not ensure overwrite, would require
	// disk to close, and would remove the temporary working file)
	tmp << disk.rdbuf(); // copy to temporary file
	tmp.close();
	
	// save header file
	std::string dir = bystd::fn_dir(fn);
	std::string basefn = bystd::fn_base(fn);
	saveHeaderENVI(dir + basefn + ".hdr");
	
	return 0;
}

int Raster::saveFileENVI ()
{
	return saveFileENVI(file.filename);
}

// Commit the current block to the temporary file.
int Raster::commit ()
{
	disk.clear();
	disk.seekg(0);
	
	if ( !disk || disk.eof() ) return 1;
	
	// save block of data
	if ( file.interleave == GridBlock::BIL ) saveBlockBIL();
	else if ( file.interleave == GridBlock::BIP ) saveBlockBIP();
	else saveBlockBSQ();
	
	block.update = false;
	
	return 0;
}

// Load header information. Returns 0 if successful, non-0 otherwise.
int Raster::loadHeaderENVI ( const std::string& fn )
{
	// read header to memory
	std::fstream tmp (fn, std::fstream::in);
	if ( !tmp || tmp.eof() ) return 1;
	
	std::string line, field, value;
	std::vector<std::string> elements;
	while ( !tmp.eof() ) {
		std::getline(tmp, line);
		if ( line.find("=") < line.size() ) {
			// commit previous field
			if ( field.size() > 0 )
				header.fields[field] = bystd::str_trim(value);
			
			// get new field
			elements = bystd::str_split(line, "=");
			field = bystd::str_trim(elements[0]);
			value = elements[1];
		}
		else value = value + line;
	}
	
	// commit previous field
	if ( field.size() > 0 )
		header.fields[field] = bystd::str_trim(value);
	
	tmp.close();
	
	// interpret header fields
	std::string nodatavalue;
	std::map<std::string, std::string>::iterator it =
		header.fields.begin();
	for ( ; it != header.fields.end(); ++it ) {
		field = it->first;
		value = it->second;
		
		if ( value.size() > 0 && value[0] == '{' ) {
			size_t i = value.find_first_not_of("{");
			size_t j = value.find_last_not_of("}");
			value = bystd::str_trim(value.substr(i, j-i+1));
		}
		
		if ( field == "samples" ) file.nc = atoi(value.c_str());
		else if ( field == "lines" ) file.nr = atoi(value.c_str());
		else if ( field == "bands" ) file.nb = atoi(value.c_str());
		else if ( field == "data type" ) {
			file.datatype = DataType(atoi(value.c_str()));
			// initialize no-data value
		}
		else if ( field == "byte order" )
			file.byteorder = atoi(value.c_str());
		else if ( field == "interleave" ) {
			value = bystd::str_upper(value);
			if ( value == "BIL" ) file.interleave = GridBlock::BIL;
			else if ( value == "BIP" ) file.interleave = GridBlock::BIP;
			else file.interleave = GridBlock::BSQ;
		}
		else if ( field == "data ignore value" )
			nodatavalue = value;
		else if ( field == "map info" )
			geography.parseENVI(value);
		else if ( field == "coordinate system string" )
			geography.coordinateSys = value;
		else if ( field == "band names" ) {
			file.bandnames = bystd::str_split(value,",");
		}
	}
	
	if ( nodatavalue.size() > 0 ) setNull(atof(nodatavalue.c_str()));
	else setNull(-9999);
	
	file.bandoffset.resize(file.nb);
	for ( int i = 0; i < file.nb; ++i ) file.bandoffset[i] = 0;
	
	return 0;
}

// Read header information and open file. Returns 0 if successful,
// non-0 otherwise.
int Raster::loadFileENVI ( const std::string& fn )
{
	file = FileInfo();
	header = HeaderInfo();
	geography = GeogInfo();
	
	file.filename = fn;
	
	std::string dir = bystd::fn_dir(fn);
	std::string basefn = bystd::fn_base(fn);
	std::string ext = bystd::fn_ext (fn);
	
	// read header to memory
	loadHeaderENVI(dir + basefn + ".hdr");
	
	std::fstream tmp (file.filename,
		std::fstream::in | std::fstream::binary);
	if ( !tmp || tmp.eof() ) return 1;
	
	disk.clear();
	disk.seekg(0);
	disk << tmp.rdbuf(); // copy to temporary file
	tmp.close();
	
	return 0;
}

// Release the current block from memory.
void Raster::unloadBlock ()
{
	if ( block.update ) commit();
	block = GridBlock ();
}

// Load a block of data from the file.
void Raster::loadBlockBSQ ()
{
	int i, n = 0;
	
	// return to the beginning of the file
	disk.clear();
	disk.seekg(0);
	
	// number of bytes to ignore before/after block's extent
	int nbi_before = block.i0 * file.nc * file.datatype.bytes;
	int nbi_after = (file.nr - (block.i0 + block.ni))
		* file.nc * file.datatype.bytes;
	int nbj_before = block.j0 * file.datatype.bytes;
	int nbj_after = (file.nc - (block.j0 + block.nj))
		* file.datatype.bytes;
	
	// number of bytes to read at a time
	int nb_read = block.nj * file.datatype.bytes;
	
	// skip bands before k0
	disk.ignore(block.k0 * file.nr * file.nc * file.datatype.bytes);
	
	for ( int k = 0; k < block.nk; ++k ) {
		// skip rows before i0
		disk.ignore(nbi_before);
		
		for ( i = 0; i < block.ni; ++i ) {
			// skip columns before j0
			disk.ignore(nbj_before);
			
			// read row
			disk.read(&block.data[n * file.datatype.bytes], nb_read);
			n += block.nj;
			
			// skip columns after j1
			disk.ignore(nbj_after);
		}
		
		// skip rows after i1
		disk.ignore(nbi_after);
	}
	
	// byteswap if necessary
	bool swap = (bystd::osLittleEndian() && (file.byteorder == 1))
		|| (!bystd::osLittleEndian() && (file.byteorder == 0));
	if ( swap ) {
		for ( int i = 0; i < block.n; i++ ) {
			bystd::byteswap(&block.data[i * file.datatype.bytes],
				file.datatype.bytes);
		}
	}
}

void Raster::loadBlockBIL ()
{
	int k, n = 0;
	
	// return to the beginning of the file
	disk.clear();
	disk.seekg(0);
	
	// number of bytes to ignore before/after block's extent
	int nbj_before = block.j0 * file.datatype.bytes;
	int nbj_after = (file.nc - (block.j0 + block.nj))
		* file.datatype.bytes;
	int nbk_before = block.k0 * file.nc * file.datatype.bytes;
	int nbk_after = (file.nb - (block.k0 + block.nk))
		* file.nc * file.datatype.bytes;
	
	// number of bytes to read at a time
	int nb_read = block.nj * file.datatype.bytes;
	
	// skip rows before i0
	disk.ignore(block.i0 * file.nc * file.nb * file.datatype.bytes);
	
	for ( int i = 0; i < block.ni; ++i ) {
		// skip bands before k0
		disk.ignore(nbk_before);
		
		for ( k = 0; k < block.nk; ++k ) {
			// skip columns before j0
			disk.ignore(nbj_before);
			
			// read row
			disk.read(&block.data[n * file.datatype.bytes], nb_read);
			n += block.nj;
			
			// skip columns after j1
			disk.ignore(nbj_after);
		}
		
		// skip bands after k1
		disk.ignore(nbk_after);
	}
	
	// byteswap if necessary
	bool swap = (bystd::osLittleEndian() && (file.byteorder == 1))
		|| (!bystd::osLittleEndian() && (file.byteorder == 0));
	if ( swap ) {
		for ( int i = 0; i < block.n; i += file.datatype.bytes )
			bystd::byteswap(&block.data[i], file.datatype.bytes);
	}
}

// Loads a block of data from the file.
void Raster::loadBlockBIP ()
{
	int j, n = 0;
	
	// return to the beginning of the file
	disk.clear();
	disk.seekg(0);
	
	// number of bytes to ignore before/after block's extent
	int nbj_before = block.j0 * file.nb * file.datatype.bytes;
	int nbj_after = (file.nc - (block.j0 + block.nj))
		* file.nb * file.datatype.bytes;
	int nbk_before = block.k0 * file.datatype.bytes;
	int nbk_after = (file.nb - (block.k0 + block.nk))
		* file.datatype.bytes;
	
	// number of bytes to read at a time
	int nb_read = block.nk * file.datatype.bytes;
	
	// skip rows before i0
	disk.ignore(block.i0 * file.nc * file.nb * file.datatype.bytes);
	
	for ( int i = 0; i < block.ni; ++i ) {
		// skip columns before j0
		disk.ignore(nbj_before);
		
		for ( j = 0; j < block.nj; ++j ) {
			// skip bands before k0
			disk.ignore(nbk_before);
			
			// read bands
			disk.read(&block.data[n * file.datatype.bytes], nb_read);
			n += block.nk;
			
			// skip bands after k1
			disk.ignore(nbk_after);
		}
		
		// skip columns after j1
		disk.ignore(nbj_after);
	}
	
	// byteswap if necessary
	bool swap = (bystd::osLittleEndian() && (file.byteorder == 1))
		|| (!bystd::osLittleEndian() && (file.byteorder == 0));
	if ( swap ) {
		for ( int i = 0; i < block.n; i += file.datatype.bytes )
			bystd::byteswap(&block.data[i], file.datatype.bytes);
	}
}

// Sets up a block to read data from the file, and then reads that
// data.
void Raster::loadBlock (
	int i, int j, int k, const BlockPattern& pattern )
{
	// get extent of block
	int i0 = i + pattern.di;
	int j0 = j + pattern.dj;
	int k0 = k + pattern.dk;
	int i1 = i0 + pattern.ni;
	int j1 = j0 + pattern.nj;
	int k1 = k0 + pattern.nk;
	
	// fit block to extent of data
	fit_to_extent(&i0, &i1, 0, file.nr);
	fit_to_extent(&j0, &j1, 0, file.nc);
	fit_to_extent(&k0, &k1, 0, file.nb);
	
	// see if more columns can be loaded
	int Vp = pattern.ni * pattern.nj * pattern.nk;
	int Vn = (i1-i0) * (j1-j0) * (k1-k0);
	int L = Vp - Vn; // unused elements
	int Ln = L / ((i1-i0) * (k1-k0)); // more columns to load
	if ( L > 0 && Ln > 0 ) {
		j1 = j1 + Ln;
		if ( j1 > file.nc ) j1 = file.nc;
	}
	
	// see if more rows can be loaded
	Vn = (i1-i0) * (j1-j0) * (k1-k0);
	L = Vp - Vn;
	Ln = L / ((j1-j0) * (k1-k0)); // more rows that to load
	if ( L > 0 && Ln > 0 ) {
		i1 = i1 + Ln;
		if ( i1 > file.nr ) i1 = file.nr;
	}
	
	// prepare block to receive data
	block = GridBlock (i0, j0, k0, i1-i0, j1-j0, k1-k0,
		file.datatype.code, file.interleave);
	
	// load block of data
	if ( file.interleave == GridBlock::BIL ) loadBlockBIL();
	else if ( file.interleave == GridBlock::BIP ) loadBlockBIP();
	else loadBlockBSQ();
}

// Request data at (i,j,k). If it is within the file but not loaded
// into memory
char* Raster::request ( int i, int j, int k )
{
	if ( !contains(i, j, k) ) return NULL;
	if ( !isLoaded(i, j, k) ) {
		unloadBlock();
		loadBlock(i, j, k, pattern);
	}
	return block(i,j,k);
}


// DEBUG ////////////////////////////////////////////////////////////

std::string Raster::getTempFilename () const
{
	return tempfilename;
}

} // namespace bygis

#endif //YOUNG_GIS_RASTER_20220407