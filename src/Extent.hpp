/////////////////////////////////////////////////////////////////////
// Extent object for tracking geographic objects' spatial extent.  //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/18/2022 - Brennan Young                                      //
// - created.                                                      //
// 10/27/2022 - Brennan Young                                      //
// - added buffer().                                               //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GIS_EXTENT_20220418
#define YOUNG_GIS_EXTENT_20220418

namespace bygis { // Brennan Young GIS namespace

class Extent {
	public:
		double xmin;
		double xmax;
		double ymin;
		double ymax;
		double zmin;
		double zmax;
	
	// constructors / destructor
	Extent (double x0=0, double x1=0, double y0=0, double y1=0,
		double z0=0, double z1=0);
	Extent (const Extent&);
	~Extent ();
	
	// operators
	Extent& operator= (const Extent&);
	Extent operator+ (const Extent&) const;
	Extent operator* (const Extent&) const;
	
	// operations
	double xlen () const;
	double ylen () const;
	double zlen () const;
	double maxlen () const;
	double areaXY () const;
	double areaXZ () const;
	double areaYZ () const;
	double area () const;
	double volume () const;
	
	Extent buffer (double, double, double) const;
	Extent buffer (double, double) const;
	Extent buffer (double) const;
	Extent intersect (const Extent&) const;
	Extent merge (const Extent&) const;
}; // Extent


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

Extent::Extent ( double x0, double x1, double y0, double y1,
	double z0, double z1 )
: xmin(x0), xmax(x1), ymin(y0), ymax(y1), zmin(z0), zmax(z1)
{}

Extent::Extent ( const Extent& ext )
: xmin(ext.xmin), xmax(ext.xmax), ymin(ext.ymin), ymax(ext.ymax),
	zmin(ext.zmin), zmax(ext.zmax)
{}

Extent::~Extent () {}


// OPERATORS ////////////////////////////////////////////////////////

Extent& Extent::operator= ( const Extent& ext )
{
	xmin = ext.xmin;
	xmax = ext.xmax;
	ymin = ext.ymin;
	ymax = ext.ymax;
	zmin = ext.zmin;
	zmax = ext.zmax;
	return *this;
}

// Extent union / merge.
Extent Extent::operator+ ( const Extent& ext ) const
{
	return merge(ext);
}

// Extent intersection.
Extent Extent::operator* ( const Extent& ext ) const
{
	return intersect(ext);
}


// OPERATIONS ///////////////////////////////////////////////////////

// Get the length of the extent.
double Extent::xlen () const
{
	return xmax - xmin;
}

double Extent::ylen () const
{
	return ymax - ymin;
}

double Extent::zlen () const
{
	return zmax - zmin;
}

double Extent::maxlen () const
{
	double x = xlen();
	double y = ylen();
	double z = zlen();
	if ( x > y ) return x > z ? x : z;
	return y > z ? y : z;
}

double Extent::areaXY () const
{
	return xlen() * ylen();
}

double Extent::areaXZ () const
{
	return xlen() * zlen();
}

double Extent::areaYZ () const
{
	return ylen() * zlen();
}

double Extent::area () const
{
	return areaXY();
}

double Extent::volume () const
{
	return xlen() * ylen() * zlen();
}

// Widen the extent in both directions.
Extent Extent::buffer ( double x, double y, double z ) const
{
	return Extent(xmin-x, xmax+x, ymin-y, ymax+y, zmin-z, zmax+z);
}
Extent Extent::buffer ( double x, double y ) const
{
	return buffer(x, y, 0);
}
Extent Extent::buffer ( double x ) const
{
	return buffer(x, x, x);
}

// Get the intersection of the two extents.
// If they do not overlap, max = min for all dimensions.
Extent Extent::intersect ( const Extent& ext ) const
{
	if ( xmin > ext.xmax || xmax < ext.xmin
			|| ymin > ext.ymax || ymax < ext.ymin
			|| zmin > ext.zmax || zmax < ext.zmin ) {
		return Extent();
	}
	else return Extent (
		xmin < ext.xmin ? ext.xmin : xmin,
		xmax > ext.xmax ? ext.xmax : xmax,
		ymin < ext.ymin ? ext.ymin : ymin,
		ymax > ext.ymax ? ext.ymax : ymax,
		zmin < ext.zmin ? ext.zmin : zmin,
		zmax > ext.zmax ? ext.zmax : zmax);
}

// Get the union of the two extents.
Extent Extent::merge ( const Extent& ext ) const
{
	return Extent (
		xmin < ext.xmin ? xmin : ext.xmin,
		xmax > ext.xmax ? xmax : ext.xmax,
		ymin < ext.ymin ? ymin : ext.ymin,
		ymax > ext.ymax ? ymax : ext.ymax,
		zmin < ext.zmin ? zmin : ext.zmin,
		zmax > ext.zmax ? zmax : ext.zmax);
}

} // namespace bygis

#endif // YOUNG_GIS_EXTENT_20220418