/////////////////////////////////////////////////////////////////////
// FileInfo object for tracking raster file information.           //
// Designed to be used as part of the Raster object.               //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/12/2022 - Brennan Young                                      //
// - created.                                                      //
// 04/26/2022 - Brennan Young                                      //
// - removed tempfilename member.                                  //
// 05/23/2022 - Brennan Young                                      //
// - transitioned to using new DataType object.                    //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GIS_FILEINFO_20220412
#define YOUNG_GIS_FILEINFO_20220412

#include <iostream>

#include "DataType.hpp"

namespace bygis { // Brennan Young GIS namespace

class FileInfo {
public:
	std::string filename;
	char* nulval; // no-data value (binary)
	int nr;       // number of lines or rows
	int nc;       // number of columns
	int nb;       // number of bands
	DataType datatype;
	unsigned char byteorder; // 0 = LSF
	unsigned char interleave; // interleave format
	std::vector<std::string> bandnames;
	std::vector<int> bandoffset; // bytes to ignore to get to band
	
	// constructors / destructor
	FileInfo (const std::string& fn="");
	FileInfo (const FileInfo&);
	~FileInfo ();
	
	// operators
	FileInfo& operator= (const FileInfo&);
}; // GridBlock


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

FileInfo::FileInfo ( const std::string& fn )
: filename(fn), nr(0), nc(0), nb(0),
	byteorder(0), interleave(0)
{
	datatype = DataType(0);
	try {
		nulval = new char[datatype.bytes < 1 ? 1 : datatype.bytes];
	}
	catch (...) {
		std::cout << "ERR: memory allocation: constructing FileInfo\n";
		exit(1);
	}
}

FileInfo::FileInfo ( const FileInfo& info )
: filename(info.filename), nr(info.nr), nc(info.nc), nb(info.nb),
	datatype(info.datatype),
	byteorder(info.byteorder), interleave(info.interleave)
{
	// copy data over
	try {
		nulval = new char [datatype.bytes < 1 ? 1 : datatype.bytes];
	}
	catch (...) {
		std::cout << "ERR: memory allocation: copying FileInfo\n";
		exit(1);
	}
	
	for ( int i = 0; i < datatype.bytes; ++i )
		nulval[i] = info.nulval[i];
	
	bandnames = info.bandnames;
	bandoffset = info.bandoffset;
}

FileInfo::~FileInfo ()
{
	delete [] nulval;
}


// OPERATORS ////////////////////////////////////////////////////////

// Assignment.
FileInfo& FileInfo::operator= ( const FileInfo& info )
{
	if ( &info == this ) return *this;
	
	// copy data over
	delete [] nulval;
	try {
		nulval = new char [datatype.bytes < 1 ? 1 : datatype.bytes];
	}
	catch (...) {
		std::cout << "ERR: memory allocation: copying FileInfo\n";
		exit(1);
	}
	
	for ( int i = 0; i < datatype.bytes; ++i )
		nulval[i] = info.nulval[i];
	
	filename = info.filename;
	nr = info.nr;
	nc = info.nc;
	nb = info.nb;
	datatype = info.datatype;
	byteorder = info.byteorder;
	interleave = info.interleave;
	bandnames = info.bandnames;
	bandoffset = info.bandoffset;
	
	return *this;
}


} // namespace bygis

#endif //YOUNG_GIS_FILEINFO_20220411