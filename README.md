# BY-GIS
Brennan Young's GIS library (the parts of it I feel comfortable sharing, anyway).

Yes, I put everything in the object's .hpp file. I get it, I'm doing in wrong, etc. etc., but linking .cpp and .hpp files gives me a headache, and I prefer to contain each object in one file.

These may depend on functions in my more general BYStdLib library. You may need to muck with the pathing in each header file's #include statements.

## Contents

Objects that the user is likely to directly interface with. Each of these objects is contained within a header file (.hpp) of the same name.

| Object | Description |
| --- | --- |
| [CompWindow](CompWindow) | Data structure for facilitating raster neighborhood operations. |
| [Graph](Graph) | Data structure for handling relationships between objects within one set of IDs. May be directed or undirected, and may be a multigraph. |
| [GridCell](GridCell) | Simple data structure containing a row, column, and band number, to facilitate communication of grid coordinates. |
| [Object](Object) | Geographic object with a location, extent, and set of properties. |
| [Raster](Raster) | Geographic object containing a rectangle-tesselated grid with a location, extent, and set of bands of diverse data types. |
| [Statistics](Statistics) | Data structure for efficiently taking samples and reporting the distribution along one dimension. |

Objects that the user is ***NOT*** likely to directly interface with. Each of these objects is contained within a header file (.hpp) of the same name.

| Object | Description |
| --- | --- |
| [DataType](DataType) | Data structure for communicating data type: short, int, long, unsigned short, unsigned int, unsigned long, float, double, complex, long complex. |
| [Extent](Extent) | Data structure for specifying the minimum and maximum x- and y-coordinates relevant to a geographic object. |
| [FileInfo](FileInfo) | Data structure for holding and managing file information for raster data. |
| [GeogInfo](GeogInfo) | Data structure for holding and managing geographic information for raster data. |
| [GridBlock](GridBlock) | Data structure for holding and managing blocks of data from raster data files. |

I will make operation libraries available as I become comfortable sharing them.

## Objects

### CompWindow

Computation window object for facilitating raster neighborhood operations, specifically culminating in the contains method.

```C++
// constructors
CompWindow (int r=1, bool c=false);          // construct a new window with radius r and circular if c=true, rectangular otherwise.
CompWindow (int di, int dj, int ri, int rj); // construct a new window where the focus is offset from the center of the window by di rows and dj columns, and is rectangular with a vertical (row) height 2*ri+1 and horizontal (column) width 2*rj+1. If rj<0, ri is instead taken as a circular radius.
CompWindow (const CompWindow& win);          // copy an existing window object.

// destructor
~CompWindow ();

// operators
CompWindow& operator= (const CompWindow& win);        // copy an existing window object.
bool operator() (int fi, int fj, int i, int j) const; // shorthand to call the contains method with focus (fi,fj) and test (i,j).
CompWindow operator+ (const CompWindow& win);         // get a window that's large enough to contain both windows being summed.
CompWindow operator+ (int r);                         // get window with r added to this window's radius.
CompWindow operator- (int r);                         // get a window with r subtracted from this window's radius.

// operations
int size() const;        // get the number of cells in the window.
int r () const;          // get the radius of the window (vertical, if different from horizontal).
int r_row () const;      // get the radius of the window in the vertical (rows) dimension.
int r_col () const;      // get the radius of the window in the horizontal (columns) dimension.
int offset_row () const; // get the row-offset of the focus from the center of the window.
int offset_col () const; // get the column-offset of the focus from the center of the window.
bool circular () const;  // get true if the window is circular, false otherwise.

void offset (int i, int j); // set the offset of the focus by i rows and j columns.
void r (int ri, int rj);    // set the vertical (row) radius to ri and horizontal (column) radius to rj.
void r (int ri);            // set the radius (rows and columns) of the window to ri.
void circular (bool);       // set whether the window is circular or not.
bool contains (int fi, int fj, int i, int j) const; // return true if (i,j) is in the window if the focus is located at (fi,fj).
```

### DataType

Object for containing information about a piece or collection of data's type. The data types recognized by this object are:

| Code | Constant Member | String | Data Type |
| --- | --- | --- | --- |
| 0 | UNKNOWN | | Unknown data type. |
| 1 | INT1 | byte, char, or int1 | Signed 1-byte integer. |
| 2 | INT2 | short, short int, or int2 | Signed 2-byte integer. |
| 3 | INT4 | int, long, or int4 | Signed 4-byte integer. |
| 14 | INT8 | long int, long long, or int8 | Signed 8-byte integer. |
| 12 | UINT2 | unsigned short or unsigned short int or uint2 | Unsigned 2-byte integer. |
| 13 | UINT4 | unsigned int, unsigned long, or uint4 | Unsigned 4-byte integer. |
| 15 | UINT8 | unsigned long int, unsigned long long, or uint8 | Unsigned 8-byte integer. |
| 4 | FLOAT4 | float or float4 | 4-byte floating-point number. |
| 5 | FLOAT8 | double, long float, or float8 | 8-byte floating-point number. |
| 6 | COMPLEX8 | complex or float complex or complex8 | complex number composed of two 8-byte floating-point numbers. |
| 9 | COMPLEX16 | long complex or double complex or complex16 | complex number composed of two 16-byte floating-point numbers. |

```C++
// constructors
DataType (const unsigned char& x=0); // construct a new data type object of type x.
DataType (const std::string& s);     // construct a new data type called s.
DataType (const DataType& d);        // copy an existing data type object.

// destructor
~DataType ();

// operators
DataType& operator= (const DataType&);        // copy an existing data type object.
bool operator< (const DataType&) const;       // compare data types, large byte size > small byte size, then float > int, then signed > unsigned.
bool operator> (const DataType&) const;
bool operator<= (const DataType&) const;
bool operator>= (const DataType&) const;
bool operator== (const DataType&) const;
bool operator!= (const DataType&) const;
bool operator== (const unsigned char&) const;

// operations
bool isInt () const;                    // return true if the type is for integer numbers.
bool isFloat () const;                  // return true if the type is for floating-point numbers.
bool isSigned () const;                 // return true if the type is signed.
std::string typeName () const;          // get the data type's name.
const unsigned char& typeCode () const; // get the data type's code.
int size () const;                      // get the data type's byte size.
```

### Extent

Extent object to track and manage the extent of geographic objects.

```C++
// constructors
Extent (double x0=0, double x1=0, double y0=0, double y1=0, double z0=0, double z1=0); // construct an extent that encompasses (x0-x1, y0-y1, z0-z1).
Extent (const Extent&); // copy an existing extent object.

// destructor
~Extent ();

// operators
Extent& operator= (const Extent&);      // copy existing extent.
Extent operator+ (const Extent&) const; // merge two extents into an extent large enough to contain them both.
Extent operator* (const Extent&) const; // intersect two extents into an extent large enough to contain only the overlapping area.

// operations
double xlen () const;   // length of the extent in the x-dimension.
double ylen () const;   // length of the extent in the y-dimension.
double zlen () const;   // length of the extent in the z-dimension.
double maxlen () const; // length of the extent's longest dimension.
double areaXY () const; // area of the extent box in the x-y plane.
double areaXZ () const; // area of the extent box in the x-z plane.
double areaYZ () const; // area of the extent box in the y-z plane.
double area () const;   // area of the extent box in the x-y plane.
double volume () const; // volume of the extent box.

Extent buffer (double x, double y, double z) const; // widen the extent in both directions by x in the x-dimension, by y in the y-dimension, and by z in the z-dimension.
Extent buffer (double x, double y) const;           // widen the extent in both directions by x in the x-dimension and by y in the y-dimension.
Extent buffer (double x) const;                     // widen the extent in both directions by x in all dimensions.
Extent intersect (const Extent&) const;             // merge two extents into an extent large enough to contain them both.
Extent merge (const Extent&) const;                 // intersect two extents into an extent large enough to contain only the overlapping area.
```

### FileInfo

Container to hold basic raster file information.

```C++
std::string filename;
char* nulval;                       // no-data value (binary), which needs to be externally deleted and then allocated to the appropriate size.
int nr;                             // number of lines or rows.
int nc;                             // number of columns.
int nb;                             // number of bands.
DataType datatype;                  // data type of the raster file.
unsigned char byteorder;            // 0 = LSF.
unsigned char interleave;           // interleave format.
std::vector<std::string> bandnames; // names of each raster band.
std::vector<int> bandoffset;        // bytes to ignore to get to each raster band.

// constructors
FileInfo ();                     // Construct a blank file information object.
FileInfo (const FileInfo& info); // copy an existing file information object.

// destructor
~FileInfo ();

// operators
FileInfo& operator= (const FileInfo& info); // copy an existing file information object.
```

### GeogInfo

Container to hold basic geographic information for a raster.

```C++
// coordinate system
std::string coordinateSys;  // name of coordinate system / datum.
std::string projectionName; // name of projected coordinate system.
int         projectionZone; // UTM zone.
std::string zoneCode;       // UTM zone code (e.g., N or S).
std::string datum;          // underlying datum.
std::string units;          // units of measurement.

// reference location (in another grid)
int         tiePoint_x;     // reference x in file coordinates.
int         tiePoint_y;     // reference y in file coordinates.

// geographic location
double      x0;             // UTM easting (top left corner).
double      y0;             // UTM northing (top left corner).
double      dx;             // cell size in x direction.
double      dy;             // cell size in y direction.

// constructors / destructor
GeogInfo ();                     // construct an empty geographic information object.
GeogInfo (const std::string& s); // construct a geographic information object from an ENVI-format string.
GeogInfo (const GeogInfo& info); // copy an existing geographic information object.
~GeogInfo ();

// operators
GeogInfo& operator= (const GeogInfo& info); // copy an existing geographic information object.

// operations
void parseENVI (const std::string& s); // infer geographic information from an ENVI-format string.

// I/O
std::string asString () const; // concatenates geographic information into an ENVI-format string.
```

### Graph

Graph object for representing networks of a set of objects which collectively have one set of IDs. Graphs may be directed or undirected and may be multigraphs. The graph is internally structured as a nested map of {object ID: {related object IDs: {relationship keys: relationship value}}}. That is, the "A" relationship from object 1 to 32 is found at [1][32]["A"].

```C++
// constructors
Graph (bool dir=true, float x=0); // construct a graph that is directed if dir=true and where 0 indicates a non-relationship value.
Graph (const Graph& G);           // copy an existing graph.

// destructor
~Graph ();

// operators
Graph& operator=(const Graph& G);     // copy an existing graph.
bool operator<(const Graph& G) const; // comparison for organization only.

// operations
size_t size () const;                                           // get the number of unique IDs represented in the graph.

std::set<int> nbrs (int i, const std::string& key) const;       // get the incoming and out-going neighbors of i for the "key" relationship, as if they graph were undirected.
std::set<int> nbrs (int i) const;                               // get the incoming and outgoing neighbors i, as if the graph were undirected.
std::set<int> nbrs_to (int i, const std::string& key) const;    // get the incoming neighbors of i for the "key" relationship.
std::set<int> nbrs_to (int i) const;                            // get the incoming neighbors of i.
std::set<int> nbrs_from (int i, const std::string& key) const;  // get the outgoing neighbors of i for the "key" relationship.
std::set<int> nbrs_from (int i) const;                          // get the outgoing neighbors of i.

std::set<int> vertices () const;                                // get the set of all IDs represented in the graph.
std::set<std::string> keys () const;                            // get the set of all relationship keys represented in the graph.
std::set<std::string> keys (int i) const;                       // get the set of all relationship keys associated with i.
std::set<std::string> keys (int i, int j) const;                // get the set of all relationship keys between i and j.

bool contains (int i, int j, const std::string& key, bool undir) const; // return true if the graph contains the "key" relationship between i and j, but only if it exists from i to j if undir=false.
bool contains_dir (int i, int j, const std::string& key) const;   // returns true if the graph contains the "key" relationship between i and j, as if the graph were directed.
bool contains_undir (int i, int j, const std::string& key) const; // returns true if the graph contains the "key" relationship between i and j, as if the graph were undirected.
bool contains (int i, int j, const std::string& key) const;       // returns true if the graph contains the "key" relationship between i and j, depending on the graph's directedness.
bool contains_dir (int i, int j) const;                           // returns true if the graph contains a relationship between i and j, as if the graph were directed.
bool contains_undir (int i, int j) const;                         // returns true if the graph contains a relationship between i and j, as if the graph were undirected.
bool contains (int i, int j) const;                               // returns true if the graph contains a relationship between i and j, depending on the graph's directedness.
bool contains_dir (int i) const;                                  // returns true if the graph contains any outgoing relationship from i.
bool contains_undir (int i) const;                                // returns true if the graph contains any relationships with i.
bool contains (int i) const;                                      // returns true if the graph contains any relationships with i, depending on the graph's directedness.

float get (int i, int j, const std::string& key) const; // get the value of the "key" relationship from i to j.
float get (int i, int j) const;                         // get the value of the "" relationship from i to j.

void set (int i, int j, const std::string& key, bool undir, float x); // set the value of the "key" relationship from i to j to x, and from j to i if undir=true.
void set (int i, int j, const std::string& key, float x);             // set the value of the "key" relationship from i to j to x, depending on the graph's directedness.
void set_dir (int i, int j, const std::string& key, float x);         // set the value of the "key" relationship from i to j to x.
void set_undir (int i, int j, const std::string& key, float x);       // set the value of the "key" relationship from i to j, and from j to i, to x.
void set (int i, int j, float x);                                     // set the value of the "" relationship from i to j to x, depending on the graph's directedness.
void set_dir (int i, int j, float x);                                 // set the value of the "" relationship from i to j to x.
void set_undir (int i, int j, float x);                               // set the value of the "" relationship from i to j, and from j to i, to x.

void clear_dir (int i, int j);                                 // remove all relationships from i to j.
void clear_undir (int i, int j);                               // remove all relationships between i and j.
void clear (int i, int j);                                     // remove all relationships between i and j, depending on the graph's directedness.
void clear (int i, int j, const std::string& key, bool undir); // remove all "key" relationships from i to j, and from j to i if undir=true.
void clear_dir (int i, int j, const std::string& key);         // remove all "key" relationships from i to j.
void clear_undir (int i, int j, const std::string& key);       // remove all "key" relationships between i and j.
void clear (int i, int j, const std::string& key);             // remove all "key" relationships between i and j, depending on the graph's directedness.
void clear_dir (int i, const std::string& key);                // remove all out-going "key" relationships from i.
void clear_undir (int i, const std::string& key);              // remove all "key" relationships associated with i.
void clear (int i, const std::string& key);                    // remove all "key" relationships associated with i, depending on the graph's directedness.
void clear (int i);                                            // remove all relationships associated with i.
void clear (const std::string& key);                           // remove all "key" relationships.
void clear ();                                                 // remove all relationships.
```

### GridBlock

Data structure for managing a dynamic block of data intended to be read from a data file and writte to a raster data file.

The grid block formalizes codes for different styles of raster data interleave:

| Code | Constant Member | Description | Equation for 1D Index |
| --- | --- | --- | --- |
| 0 | BSQ | Band-sequential. | (k-k0)*ni*nj + (i-i0)*nj + (j-j0) |
| 1 | BIL | Band-interleave by line. | (i-i0)*nj*nk + (k-k0)*nj + (j-j0) |
| 2 | BIP | Band-interleave by pixel. | nk * ((i-i0)*nj + (j-j0)) + (k-k0) |

```C++
// public members
int i0;               // first loaded row
int j0;               // first loaded column
int k0;               // first loaded band
int ni;               // number of loaded rows
int nj;               // number of loaded columns
int nk;               // number of loaded bands
int n;                // number of loaded elements
DataType datatype;    // data type
unsigned char format; // interleave format
char* data;           // loaded data (binary)
bool update;          // set to true if changes may have been made since load

// constructors
GridBlock (int row0=0, int col0=0, int band0=0, int nrows=0, int ncols=0, int nbands=0, unsigned char datatype=4, unsigned char interleave=BSQ); // constructor to assign enough memory for a block with nrows, ncols, and nbands of the given data type. Data are organized as indicated with the interleave. The index in the block is offset by row0, col0, band 0.
GridBlock (const GridBlock& block); // copy the block of data.

// destructor
~GridBlock ();

// operators
GridBlock& operator= (const GridBlock& block); // copy the block of data.
char* operator() (int i, int j, int k);        // get the address to the data located at row i, column j, band k.

// operations
bool contains (int i, int j, int k) const;     // return true if the block contains row i, column j, band k.
int index (int i, int j, int k) const;         // return the index of row i, column j, band k.
}; // GridBlock
```

### GridCell

Data structure for containing basic grid cell information (row, column, and band) within a raster, and basic adjacency relationships with other cells.

Adjacency relationships are determined with one of three rules:

| Code | Cosntant Member | Description |
| --- | --- | --- |
| 0 | QUEEN | Queen's rule (*): all 4 edge-adjacent and 4 corner-adjacent cells count as neighbors. |
| 1 | ROOK | Rook's rule (+): all 4 edge-adjacent cells count as neighbors. |
| 2 | BISHOP | Bishop's rule (X): all 4 corner-adjacent cells count as neighbors. |

```C++
// public members
int i;
int j;
int k;

// constructors
GridCell (int i=-1, int j=-1, int k=0); // construct new grid cell for row i, column j, band k.
GridCell (const GridCell& cell);        // copy grid cell.

// destructor
~GridCell ();

// operators
GridCell& operator= (const GridCell& cell);   // copy grid cell.
bool operator== (const GridCell& cell) const; // return true if this and the given cell have the same (i,j,k).
bool operator!= (const GridCell& cell) const; // return true if this and the given cell do not have the same (i,j,k).
bool operator< (const GridCell& cell) const;  // comparison, where band is checked first, then row, then column.

// operations
bool isNbr (const GridCell& cell, const unsigned char& rule) const; // return true if the two cells are neighbors by the given rule.
bool isNbr (const GridCell& cell) const;                            // return true if the two cells are neighbors by Queen's rule.
```

### Object

Geographic object, represented with a position, an extent, and a set of properties. The exact extents of objects are represented with an external raster or polygon object.

```C++
// public members
int id;     // object identifier.
double x;   // centroid x-position.
double y;   // centroid y-position.
double z;   // centroid z-position.
Extent ext; // object's Extent object.

// constructors
Object (int i=0, double X=0, double Y=0, double Z=0, const Extent& ex=Extent()); // construct a new object with ID=i, centroid coordinates=(x,y,z), and the given extent.
Object (int, const Extent&); // construct a new object with ID=i, the given extent, and centroid coordinates inferred from the extent.
Object (const Object &);     // copy an existing object.

// destructor
~Object ();

// operators
Object& operator=(const Object&); // copy an existing object.

// operations
std::vector<std::string> properties () const; // get a vector of names or keys for all of the object's properties.
bool contains (const std::string& key) const; // return true if the object has a property called "key".
double get (const std::string& key) const;    // get the value of the object's "key" property.
void set (const std::string& key, double x);  // set the value of the object's "key" property to x.
void clear (const std::string& key);          // remove the object's "key" property.
void clear ();                                // remove all of the object's properties.
```

### Raster

Raster object for processing large files. Using these objects isn't fast and takes up space on your hard drive. Creates a copy of the file you plan to manipulate (or creates a new file for new rasters) and reads/writes blocks of data from/to it. If you run out of disk space, your program will have undefined behavior unless you have checks for file and grid data validity. Aborting the program early for any reason will prevent the destructor from being called, which destroys the temporary files -- these will need to be destroyed manually. Temporary raster files have names that begin with "_tmp_rast" and can be found at the directory indicated with the static temppath member, which can be retrieved or set (for ALL raster objects in the program) with the tempFileLoc method. I strongly recommend setting the temppath only once at the beginning of the program. The raster object currently only works for ENVI-format images with binary .dat files with the .hdr header files (https://www.l3harrisgeospatial.com/docs/enviheaderfiles.html, accessed 24 Jan 2023).

```C++
// constructors
Raster (int i, int j, int k);    // create a new reaster with i rows, j columns, and k bands.
Raster (const std::string& fn);  // load a raster file.
Raster (const Raster& r);        // copy an existing raster object.
Raster (const Raster& r, int k); // create an empty copy of the raster object (same rows, columns, data type, no-data value) with k bands (minimum 1).

// destructor
~Raster (); // destroys the raster object and deletes its temporary file.

// operators
Raster& operator= (const Raster& r); // copy an existing raster object.

// operations (element access)
double i_to_Y (int i) const;                    // convert row to geographic y-coordinate.
double j_to_X (int j) const;                    // convert column to geographic x-coordinate.
void ij_to_XY (int i, int j, double* x, double* y) const; // convert row and column to geographic coordinates.
void cell_to_XY (const GridCell& cell, double* x, double* y) const; // convert a GridCell object to geographic coordinates.
int X_to_j (double x) const;                    // convert x-coordinate to column.
int Y_to_i (double y) const;                    // convert y-coordinate to row.
void XY_to_ij (double x, double y, int* i, int* j) const; // convert geographic coordinates to row and column.
GridCell XY_to_cell (double x, double y) const; // convert geographic coordinates to a GridCell object.

int nrows () const;  // get the number of rows in the grid.
int ncols () const;  // get the number of columns in the grid.
int nbands () const; // get the number of bands in the grid.
int size () const;   // get the number of elements in each band (nrows x ncols).
int volume () const; // get the number of elements in the grid (nrows x ncols x nbands).

double x0 () const; // get the geographic x-coordinate for the upper-left corner of cell (0,0).
double y0 () const; // get the geographic y-coordinate for the upper-left corner of cell (0,0).
double dx () const; // get the geographic length in the x-dimension of each cell.
double dy () const; // get the geographic length in the y-dimension of each cell.

const GeogInfo& geog() const; // get the raster's geographic information in a GeogInfo object.
void geog(const GeogInfo&);   // set the raster's geographic information to a given GeogInfo object.

char* get (int i, int j, int k);         // get the address for the value at row i, column j, band k.
char* get (int i, int j);                // get the address for the value at row i, column j, band 0.
char* get (const GridCell& cell);        // get the address for the value in the given cell.
char* getXY (double x, double y, int k); // get the address for the value at the given geographic coordinates in band k.
char* getXY (double x, double y);        // get the address for the value at the given goegraphic coordinates in band 0.

long long geti (int i, int j, int k);          // get the integer value at row i, column j, band k.
long long geti (int i, int j);                 // get the integer value at row i, column j, band 0.
long long geti (const GridCell& cell);         // get the integer value in the given cell.
long long getiXY (double x, double y, int k);  // get the integer value at the given geographic coordinates in band k.
long long getiXY (double x, double y);         // get the integer value at the given geographic coordinates in band k.

double getf (int i, int j, int k);         // get the floating-point value at row i, column j, band k.
double getf (int i, int j);                // get the floating-point value at row i, column j, band 0.
double getf (const GridCell& cell);        // get the floating-point value in the given cell.
double getfXY (double x, double y, int k); // get the floating-point value at the given geographic coordinates in band k.
double getfXY (double x, double y);        // get the floating-point value at the given geographic coordinates in band k.

void set (int i, int j, int k, const T& v);         // set row i, column j, band k to the value v.
void set (int i, int j, const T& v);                // set row i, column j, band 0 to the value v.
void set (const GridCell& cell, const T& v);        // set the cell to the value v.
void setXY (double x, double y, int k, const T& v); // set the cell at the given geographic coordinates in band k to the value v.
void setXY (double x, double y, const T& v);        // set the cell at the given geographic coordinates in band 0 to the value v.

void setNull (int i, int j, int k);         // set row i, column j, band k to no-data.
void setNull (int i, int j);                // set row i, column j, band 0 to no-data.
void setNull (const GridCell& cell);        // set the cell to no-data.
void setNullXY (double x, double y, int k); // set the cell at the given geographic coordinates in band k to no-data.
void setNullXY (double x, double y);        // set the cell at the given geographic coordinates in band 0 to no-data.

template<typename T> void setNull (const T& v);      // set the value recognized as no-data to v.
template<typename T> bool isNull (const T& v) const; // returns true if the value is the no-data value.
bool isNull (int i, int j, int k);                   // returns true if row i, column j, band k is no-data.
bool isNull (int i, int j);                          // returns true if row i, column j, band 0 is no-data.
bool isNull (const GridCell& cell);                  // returns true if the cell is no-data.
bool isNullXY (double x, double y, int k);           // returns true if the cell at the given geographic coordinates in band k is no-data.
bool isNullXY (double x, double y);                  // returns true if the cell at the given geographic coordinates in band 0 is no-data.
long long getNulli () const;                         // get the no-data value as an integer.
double getNullf () const;                            // get the no-data value as a floating-point number.

// operations
bool isValid () const;                      // returns true if the raster has an open temporary file. This does not guarantee that the temporary file was created large enough.
bool contains (int i, int j) const;         // returns true if row i, column j can be found within the raster.
bool contains (int i, int j, int k) const;  // returns true if row i, column j, band k can be found within the raster.
bool contains (const GridCell& cell) const; // returns true if the cell can be found within the raster.
bool containsXY (double x, double y) const; // returns true if the cell at the given geographic coordinates can be found within the raster.
bool isLoaded (int i, int j) const;         // returns true if row i, column j is in the currently loaded GridBlock.
bool isLoaded (int i, int j, int k) const;  // returns true if row i, column j, band k is in the currently loaded GridBlock.
bool isLoaded (const GridCell& cell) const; // returns true if the cell is in the currently loaded GridBlock.
bool isLoadedXY (double x, double y) const; // returns true if the cell at the given geographic coordinates is in the currently loaded GridBlock.

int blockSize_row() const;    // get the number of rows targeted to load for a GridBlock.
int blockSize_col() const;    // get the number of columns targeted to load for a GridBlock.
int blockSize_band() const;   // get the number of bands targeted to load for a GridBlock.
void blockSize_row(int i);    // set the number of rows to target for loading Gridblocks.
void blockSize_col(int j);    // set the number of columns to target for loading GridBlocks.
void blockSize_band(int k);   // set the number of bands to target for loading GridBlocks.
int blockOffset_row() const;  // get the row-offset for loading Gridblocks.
int blockOffset_col() const;  // get the column-offset for loading Gridblocks.
int blockOffset_band() const; // get the band-offset for loading Gridblocks.
void blockOffset_row(int i);  // set the row-offset for loading GridBlocks.
void blockOffset_col(int j);  // set the column-offset for loading Gridblocks.
void blockOffset_band(int k); // set the band-offset for loading Gridblocks.
void setBlockPattern(int di, int dj, int dk, int ni, int nj, int nk); // set the GridBlock loading target to row-offset di, column-offset dj, band-offset dk with numbers of rows ni, columns nj, and bands nk.

Raster sub (int i0, int i1, int j0, int j1, int k0, int k1); // clip the raster to rows i0-i1, columns j0-j1, and bands k0-k1 (i1, j1, k1 not inclusive).
Raster sub (int i0, int i1, int j0, int j1);                 // clip the raster to rows i0-i1, columns j0-j1 with all bands (i1, j1 not inclusive).
Raster sub (int k0, int k1);                                 // clip the raster to bands k0-k1 (k1 not inclusive).
Raster sub (int i0, int i1, int j0, int j1, const std::vector<int>& k); // clip the raster to rows i0-i1, j0-j1, and bands contained in k (i1, j1 not inclusive).
Raster sub (const std::vector<int>& k);                      // clip the raster to bands contained in k.

// data type conversion -- see data type codes and strings listed with the DataType object.
const DataType& getDataType () const;          // get the raster's data type indicated in a DataType object.
Raster convert (const DataType& dt, int k);    // create an emtpy copy of the raster (same dimensions) with k bands of the specified data type.
Raster convert (const DataType& dt);           // copy the raster and convert its contents to the specified data type.
Raster convert (unsigned char dt, int k);      // create an empty copy of the raster (same dimensions) with the data type indicated with the given code.
Raster convert (unsigned char dt);             // copy the raster and convert its contents to the data type indicated with the given code.
Raster convert (const std::string& dt, int k); // create an empty copy of the raster (same dimensions) with the data type indicated with the given string.
Raster convert (const std::string& dt);        // copy the raster and convert its contents to the data type indicated with the given string.

// I/O
static std::string tempFileLoc ();            // get the universal Raster object temporary file location.
static void tempFileLoc (const std::string&); // set the universal Raster object temporary file location.

int saveHeaderENVI (const std::string& fn); // save just the raster's header file.
int saveFileENVI (const std::string& fn);   // save the raster's data and header files.
int saveFileENVI ();                        // update the raster's existing file.
int commit ();                              // commit the currently loaded GridBlock to the raster's temporary file.

int loadHeaderENVI (const std::string& fn); // load the raster's header file.
int loadFileENVI (const std::string& fn);   // load the raster's header information and data.

void unloadBlock ();                        // unload the raster's currently loaded GridBlock. This happens invisibly with get/set methods.
void loadBlock (int i, int j, int k, const BlockPattern&); // load a GridBlock targeting row i, column j, band k. This happens invisibly with get/set methods.
char* request (int i, int j, int k);        // request the cell located at row i, column j, band k. Unloads the currently loaded block if it does not contain (i,j,k), and loads the relevant block containing (i,j,k). This happens invisibly with get/set methods.

// debug
std::string getTempFilename () const; // get the name of the raster's temporary file.
```

### Statistics

Object to help handle monovariate statistical reductions. If determining the median and/or mode, values are truncated to integers and values are tracked in a map structure.

```C++
// constructors
Statistics (bool cv=false);           // construct a new statistics object, which will count individual values (i.e., for median and mode determination) if cv=true.
Statistics (const Statistics& stats); // copy an existing statistics object.

// destructor
~Statistics ();

// operators
Statistics& operator= (const Statistics& stats); // copy an existing statistics object.

// operations
void push (double x);                                     // sample the value x into the object's internal aggregate.
template<typename T> void push (const std::vector<T>& x); // sammple all values x into the object's internal aggregate.
void trackValues (bool);                                  // set to true to track individual values.
bool trackValues () const;                                // return true if individual values are being tracked.
int count () const;                                       // return the number of samples taken.
double min () const;                                      // get the minimum value sampled.
double max () const;                                      // get the maximum number sampled.
double range () const;                                    // get the maximum-minimum difference.
double sum () const;                                      // get the sum of sampled values.
double sum2 () const;                                     // get the square of the sum of sampled values.
double mean () const;                                     // get the mean of sampled values.
double var () const;                                      // get the variance of sampled values.
double stdev () const;                                    // get the standard devation of sampled values.
double skew () const;                                     // get the skew of sampled values.
double kurtosis () const;                                 // get the kurtosis of sampled values.
double median () const;                                   // get the median of sampled values, for all values that were tracked.
int mode () const;                                        // get the mode of sampled values, for all values that were tracked.
```

## Operations

I'm not yet ready to make these public. I'm using these for my work, and I've been taught to be both paranoid and open with my code and data. This repository is my compromise between those two end-members.
