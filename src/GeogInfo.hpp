/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 07/24/2017 - Brennan Young                                      //
// - created.                                                      //
// 04/07/2020 - Brennan Young                                      //
// - corrected compiler warnings.                                  //
// 01/19/2022 - Brennan Young                                      //
// - removed unecessary #include libraries.                        //
// - optimized constructors by parameterizing.                     //
// - removed clear() member function (just invoke constructor).    //
// 04/07/2022 - Brennan Young                                      //
// - copied old Grid MapInfo to be part of Raster object, as a     //
//   place holder while the more robust GeogInfo object is         //
//   developed into something more robust.                         //
// - added coordinateSys string member.                            //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GIS_GEOGINFO_20220407
#define YOUNG_GIS_GEOGINFO_20220407

#include "../BYstdlib/StrManip.hpp"

#include <cstdlib>
#include <sstream>
#include <string>

namespace bygis { // Brennan Young GIS namespace

class GeogInfo {
public:
	// coordinate system
	std::string coordinateSys;
	std::string projectionName; // name of projected coordinate system
	int         projectionZone; // UTM zone
	std::string zoneCode;       // UTM zone code (e.g., N or S)
	std::string datum;          // underlying datum
	std::string units;          // units of measurement
	
	// reference location (in another grid)
	int         tiePoint_x;     // reference x in file coordinates
	int         tiePoint_y;     // reference y in file coordinates
	
	// geographic location
	double      x0;             // UTM easting (top left corner)
	double      y0;             // UTM northing (top left corner)
	double      dx;             // cell size in x direction
	double      dy;             // cell size in y direction
	
	// constructors / destructor
	GeogInfo ();
	GeogInfo (const std::string&);
	GeogInfo (const GeogInfo&);
	~GeogInfo ();
	
	// operators
	GeogInfo& operator= (const GeogInfo&);
	
	// operations
	void parseENVI (const std::string&);
	
	// I/O
	std::string asString () const;
}; // GeogInfo


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

// Default constructor.
GeogInfo::GeogInfo ()
: coordinateSys("undefined"),
	projectionName("undefined"), projectionZone(0), zoneCode("N"),
	datum("undefined"), units("undefined"), tiePoint_x(0),
	tiePoint_y(0), x0(0.0), y0(0.0), dx(1.0), dy(1.0)
{}

// Constructor.
GeogInfo::GeogInfo ( const std::string& s )
: coordinateSys("undefined"),
	projectionName("undefined"), projectionZone(0), zoneCode("N"),
	datum("undefined"), units("undefined"),
	tiePoint_x(0), tiePoint_y(0),
	x0(0.0), y0(0.0), dx(1.0), dy(1.0)
{
	parseENVI(s);
}

GeogInfo::GeogInfo ( const GeogInfo& info )
: coordinateSys(info.coordinateSys),
	projectionName(info.projectionName),
	projectionZone(info.projectionZone),
	zoneCode(info.zoneCode), datum(info.datum), units(info.units),
	tiePoint_x(info.tiePoint_x), tiePoint_y(info.tiePoint_y),
	x0(info.x0), y0(info.y0), dx(info.dx), dy(info.dy)
{}

GeogInfo::~GeogInfo ()
{}


// OPERATORS ////////////////////////////////////////////////////////

GeogInfo& GeogInfo::operator= ( const GeogInfo& info )
{
	if ( &info == this ) return *this;
	
	coordinateSys   = info.coordinateSys;
	projectionName  = info.projectionName;
	projectionZone  = info.projectionZone;
	zoneCode        = info.zoneCode;
	datum           = info.datum;
	units           = info.units;
	tiePoint_x      = info.tiePoint_x;
	tiePoint_y      = info.tiePoint_y;
	x0              = info.x0;
	y0              = info.y0;
	dx              = info.dx;
	dy              = info.dy;
	
	return *this;
}


// OPERATIONS ///////////////////////////////////////////////////////

// Parses a string for map info, based on the format of the 'map
// info' field in ENVI header files.
void GeogInfo::parseENVI ( const std::string& s )
{
	if ( s.length() == 0 ) {
		*this = GeogInfo();
		return;
	}
	
	// split each element
	size_t a = s.find_first_not_of('{');
	size_t b = s.find_last_not_of('}');
	if ( b < s.size() ) ++b;
	std::vector<std::string> ss =
		bystd::str_split(s.substr(a, b), ',');
	
	// trim white space
	size_t i = 0;
	for ( ; i < ss.size(); ++i ) ss[i] = bystd::str_trim(ss[i]);
	
	// assign values from input string
	i = 0;
	if ( i < ss.size() ) projectionName = ss[i++];
	if ( i < ss.size() ) tiePoint_x     = atoi(ss[i++].c_str());
	if ( i < ss.size() ) tiePoint_y     = atoi(ss[i++].c_str());
	if ( i < ss.size() ) x0             = atof(ss[i++].c_str());
	if ( i < ss.size() ) y0             = atof(ss[i++].c_str());
	if ( i < ss.size() ) dx             = atof(ss[i++].c_str());
	if ( i < ss.size() ) dy             = atof(ss[i++].c_str());
	if ( i < ss.size() ) projectionZone = atoi(ss[i++].c_str());
	if ( i < ss.size() ) zoneCode       = ss[i++];
	if ( i < ss.size() ) datum          = ss[i++];
	if ( i < ss.size() ) units          = ss[i++];
}

// Compiles a string representation of the map info object using the
// format of the 'map info' component of ENVI header files.
std::string GeogInfo::asString () const
{
	std::stringstream ss;
	ss << std::fixed << "{"
		<< projectionName << ", "
		<< tiePoint_x     << ", "
		<< tiePoint_y     << ", "
		<< x0             << ", "
		<< y0             << ", "
		<< dx             << ", "
		<< dy             << ", "
		<< projectionZone << ", "
		<< zoneCode       << ", "
		<< datum          << ", "
		<< units          << "}";
	return ss.str();
}

} // bygis namespace

#endif // YOUNG_GIS_GEOGINFO_20220407