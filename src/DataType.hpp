/////////////////////////////////////////////////////////////////////
// Object to help with interpreting and communicating data types.  //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 05/23/2022 - Brennan Young                                      //
// - created.                                                      //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GIS_DATATYPE_20220523
#define YOUNG_GIS_DATATYPE_20220523

#include "../BYstdlib/StrManip.hpp"

namespace bygis { // Brennan Young GIS namespace

class DataType {
public:
	// data type codes
	static const unsigned char UNKNOWN;
	static const unsigned char INT1; // 8-bit int
	static const unsigned char INT2; // 16-bit int
	static const unsigned char INT4; // 32-bit int
	static const unsigned char FLOAT4; // 32-bit float
	static const unsigned char FLOAT8; // 64-bit float
	static const unsigned char COMPLEX8; // float complex
	static const unsigned char COMPLEX16; // double complex
	static const unsigned char UINT2; // 16-bit unsigned int
	static const unsigned char UINT4; // 32-bit unsigned int
	static const unsigned char INT8; // 64-bit int
	static const unsigned char UINT8; // 64-bit unsigned int
	static unsigned char dataTypeSize (unsigned char);
	
	// data type information
	std::string name;
	unsigned char code;
	unsigned char bytes;
	
	// constructors, destructor
	DataType (const unsigned char& x=0);
	DataType (const std::string&);
	DataType (const DataType&);
	~DataType ();
	
	// operators
	DataType& operator= (const DataType&);
	bool operator< (const DataType&) const;
	bool operator> (const DataType&) const;
	bool operator<= (const DataType&) const;
	bool operator>= (const DataType&) const;
	bool operator== (const DataType&) const;
	bool operator!= (const DataType&) const;
	bool operator== (const unsigned char&) const;
	
	// operations
	bool isInt () const;
	bool isFloat () const;
	bool isSigned () const;
	std::string typeName () const;
	const unsigned char& typeCode () const;
	int size () const;
}; // DataType

const unsigned char DataType::UNKNOWN = 0;
const unsigned char DataType::INT1 = 1;
const unsigned char DataType::INT2 = 2;
const unsigned char DataType::INT4 = 3;
const unsigned char DataType::FLOAT4 = 4;
const unsigned char DataType::FLOAT8 = 5;
const unsigned char DataType::COMPLEX8 = 6;
const unsigned char DataType::COMPLEX16 = 9;
const unsigned char DataType::UINT2 = 12;
const unsigned char DataType::UINT4 = 13;
const unsigned char DataType::INT8 = 14;
const unsigned char DataType::UINT8 = 15;


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

DataType::DataType ( const unsigned char& x )
: code(x)
{
	name = typeName();
	code = typeCode();
	bytes = size();
}

DataType::DataType ( const std::string& s )
{
	name = bystd::str_lower(s);
	code = typeCode();
	name = typeName();
	bytes = size();
}

DataType::DataType ( const DataType& dt )
: name(dt.name), code(dt.code), bytes(dt.bytes)
{}

DataType::~DataType ()
{}


// OPERATORS ////////////////////////////////////////////////////////

DataType& DataType::operator= ( const DataType& dt )
{
	name = dt.name;
	code = dt.code;
	bytes = dt.bytes;
	return *this;
}

// logical operators compare bytesize, then float > int,
// then signed > unsigned.
bool DataType::operator< ( const DataType& dt ) const
{
	return bytes < dt.bytes
		|| (bytes == dt.bytes && isInt() && dt.isFloat())
		|| (bytes == dt.bytes
			&& ((isInt() && dt.isInt()) | (isFloat() && dt.isFloat()))
			&& !isSigned() && dt.isSigned());
}

bool DataType::operator> ( const DataType& dt ) const
{
	return bytes > dt.bytes
		|| (bytes == dt.bytes && isFloat() && dt.isInt())
		|| (bytes == dt.bytes
			&& ((isInt() && dt.isInt()) || (isFloat() && dt.isFloat()))
			&& isSigned() && !dt.isSigned());
}

bool DataType::operator<= ( const DataType& dt ) const
{
	return code == dt.code || bytes < dt.bytes
		|| (bytes == dt.bytes && isInt() && dt.isFloat())
		|| (bytes == dt.bytes
			&& ((isInt() && dt.isInt()) | (isFloat() && dt.isFloat()))
			&& !isSigned() && dt.isSigned());
}

bool DataType::operator>= ( const DataType& dt ) const
{
	return code == dt.code || bytes > dt.bytes
		|| (bytes == dt.bytes && isFloat() && dt.isInt())
		|| (bytes == dt.bytes
			&& ((isInt() && dt.isInt()) || (isFloat() && dt.isFloat()))
			&& isSigned() && !dt.isSigned());
}

bool DataType::operator== ( const DataType& dt ) const
{
	return code == dt.code;
}

bool DataType::operator!= ( const DataType& dt ) const
{
	return code != dt.code;
}

bool DataType::operator== ( const unsigned char& x ) const
{
	return code == x;
}


// OPERATIONS ///////////////////////////////////////////////////////

bool DataType::isInt () const
{
	return code == INT1 || code == INT2 || code == INT4
		|| code == UINT2 || code == UINT4 || code == UINT8;
}

bool DataType::isFloat () const
{
	return code == FLOAT4 || code == FLOAT8
		|| code == COMPLEX8 || code == COMPLEX16;
}

bool DataType::isSigned () const
{
	return !(code == UNKNOWN
		|| code == UINT2 || code == UINT4 || code == UINT8);
}

std::string DataType::typeName () const
{
	if ( code == INT1 ) return "byte";
	if ( code == INT2 ) return "short";
	if ( code == INT4 ) return "int";
	if ( code == INT8 ) return "long int";
	if ( code == UINT2 ) return "unsigned short";
	if ( code == UINT4 ) return "unsigned int";
	if ( code == UINT8 ) return "unsigned long int";
	if ( code == FLOAT4 ) return "float";
	if ( code == FLOAT8 ) return "double";
	if ( code == COMPLEX8 ) return "complex";
	if ( code == COMPLEX16 ) return "long complex";
	return "unknown";
}

const unsigned char& DataType::typeCode () const
{
	if ( name == "byte" || name == "char" || name == "int1" )
		return INT1;
	if ( name == "short" || name == "short int" || name == "int2" )
		return INT2;
	if ( name == "int" || name == "long" || name == "int4" )
		return INT4;
	if ( name == "long int" || name == "long long" || name == "int8" )
		return INT8;
	if ( name == "unsigned short" || name == "unsigned short int"
			|| name == "uint2" )
		return UINT2;
	if ( name == "unsigned int" || name == "unsigned long"
			|| name == "uint4" )
		return UINT4;
	if ( name == "unsigned long int" || name == "unsigned long long"
			|| name == "uint8" )
		return UINT8;
	if ( name == "float" || name == "float4" )
		return FLOAT4;
	if ( name == "double" || name == "long float" || name == "float8" )
		return FLOAT8;
	if ( name == "complex" || name == "float complex"
			|| name == "complex8" )
		return COMPLEX8;
	if ( name == "long complex" || name == "double complex"
			|| name == "complex16" )
		return COMPLEX16;
	return UNKNOWN;
}

int DataType::size () const
{
	if ( code == INT1 ) return 1;
	if ( code == INT2 || code == UINT2 ) return 2;
	if ( code == INT4 || code == UINT4 || code == FLOAT4 ) return 4;
	if ( code == INT8 || code == UINT8 || code == FLOAT8
			|| code == COMPLEX8 )
		return 8;
	if ( code == COMPLEX16 ) return 16;
	return 0;
}

} // namespace bygis

#endif // YOUNG_GIS_DATATYPE_20220523