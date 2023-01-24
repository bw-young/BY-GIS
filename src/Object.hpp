/////////////////////////////////////////////////////////////////////
// Generic spatial object class, container of basic extent         //
// information and object properties.                              //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/18/2022 - Brennan Young                                      //
// - created.                                                      //
// 05/09/2022 - Brennan Young                                      //
// - get() can now return a member variable (id, x, y, or z).      //
// 05/17/2022 - Brennan Young                                      //
// - contains() now returns true for member variables (id,x,y,z).  //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GIS_OBJECT_20220418
#define YOUNG_GIS_OBJECT_20220418

#include <map>
#include <string>

#include "Extent.hpp"

namespace bygis { // Brennan Young GIS namespace

class Object {
private:
	std::map<std::string, double> prop;
public:
	int id;
	
	// object/centroid position
	double x;
	double y;
	double z;
	Extent ext;
	
	// constructors, destructor
	Object (int i=0, double X=0, double Y=0, double Z=0,
		const Extent& ex=Extent());
	Object (int, const Extent&);
	Object (const Object &);
	~Object ();
	
	// operators
	Object& operator=(const Object&);
	
	// operations
	std::vector<std::string> properties () const;
	bool contains (const std::string&) const;
	double get (const std::string&) const;
	void set (const std::string&, double);
	void clear (const std::string&);
	void clear ();
}; // Object


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

Object::Object (
	int i, double X, double Y, double Z, const Extent& ex )
: id(i), x(X), y(Y), z(Z), ext(ex)
{}

Object::Object ( int i, const Extent& ex )
: id(i), x(0.5 * (ex.xmin + ex.xmax)),
	y(0.5 * (ex.ymin + ex.ymax)),
	z(0.5 * (ex.zmin + ex.zmax)), ext(ex)
{}

Object::Object ( const Object & obj )
: id(obj.id), x(obj.x), y(obj.y), z(obj.z), ext(obj.ext)
{
	prop = obj.prop;
}

Object::~Object () {}


// OPERATORS ////////////////////////////////////////////////////////

Object& Object::operator= ( const Object& obj )
{
	id = obj.id;
	x = obj.x;
	y = obj.y;
	z = obj.z;
	ext = obj.ext;
	prop = obj.prop;
	return *this;
}


// OPERATIONS ///////////////////////////////////////////////////////

// Returns a vector of the object's properties.
std::vector<std::string> Object::properties () const
{
	std::vector<std::string> p;
	p.reserve(prop.size());
	std::map<std::string, double>::const_iterator it = prop.begin();
	for ( ; it != prop.end(); ++it ) p.push_back(it->first);
	return p;
}

// Returns true if the object has a value for the given property.
bool Object::contains ( const std::string& s ) const
{
	return prop.find(s) != prop.end()
		|| s == "id" || s == "x" || s == "y" || s == "z";
}

// Returns the value of the given property. If the property does
// not exist, returns -9999.
double Object::get ( const std::string& s ) const
{
	std::map<std::string, double>::const_iterator it = prop.find(s);
	if ( it == prop.end() ) {
		if ( s == "id" ) return id;
		if ( s == "x" ) return x;
		if ( s == "y" ) return y;
		if ( s == "z" ) return z;
		return -9999;
	}
	return it->second;
}

// Set the value for the given property. If it does not exist,
// creates it. If it already exists, overwrites it.
void Object::set ( const std::string& s, double x )
{
	prop[s] = x;
}

// Remove the given property.
void Object::clear ( const std::string& s )
{
	prop.erase(s);
}

// Remove all of the object's properties.
void Object::clear ()
{
	prop.clear();
}

} // namespace bygis

#endif // YOUNG_GIS_OBJECT_20220418