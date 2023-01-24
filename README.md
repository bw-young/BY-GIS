# BY-GIS
Brennan Young's GIS library (the parts of it I feel comfortable sharing, anyway).

Yes, I put everything in the object's .hpp file. I get it, I'm doing in wrong, etc. etc., but linking .cpp and .hpp files gives me a headache, and I prefer to contain each object in one file.

These may depend on functions in my more general BYStdLib library. You may need to muck with the pathing in each header file's #include statements.

## Contents

Objects that the user is likely to directly interface with:
| Object | Description |
| --- | --- |
| [CompWindow](CompWindow) | Data structure for facilitating raster neighborhood operations. |
| [Graph](Graph) | Data structure for handling relationships between objects within one set of IDs. May be directed or undirected, and may be a multigraph. |
| [GridCell](GridCell) | Simple data structure containing a row, column, and band number, to facilitate communication of grid coordinates. |
| [Object](Object) | Geographic object with a location, extent, and set of properties. |
| [Raster](Raster) | Geographic object containing a rectangle-tesselated grid with a location, extent, and set of bands of diverse data types. |
| [Statistics](Statistics) | Data structure for efficiently taking samples and reporting the distribution along one dimension. |

Objects that the user is ***NOT*** likely to directly interface with:
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

### DataType

### Environment

### Extent

### FileInfo

### GeogInfo

### Graph

### GridBlock

### GridCell

### Object

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
Raster& operator= (const Raster&); // copy an existing raster object.

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

## Operations

I'm not yet ready to make these public. I'm using these for my work, and I've been taught to be both paranoid and open with my code and data. This repository is my compromise between those two end-members.
