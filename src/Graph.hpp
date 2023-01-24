/////////////////////////////////////////////////////////////////////
// Graph network class, container of basic information about       //
// relationships between objects. Can represent directed and       //
// undirected graphs, and multigraphs.                             //
//                                                                 //
// Does not assume undirected behavior -- need to assign (i,j)     //
// as well as (j,i), or use the set method for undirected          //
// relationships.                                                  //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 11/11-21/2022 - Brennan Young                                   //
// - created.                                                      //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GIS_GRAPH_20221111
#define YOUNG_GIS_GRAPH_20221111

#include <map>
#include <cmath>
#include <set>
#include <string>

namespace bygis { // Brennan Young GIS namespace

class Graph {
private:
	static const unsigned char UNDIRECTED;
	static const unsigned char FROM;
	static const unsigned char TO;
	
	// bool in pair is false if the relationship only exists *to*
	// vertex i from vertex j
	std::map<int, std::map<int, std::map<std::string,
		std::pair<bool, float> > > > data;
	
	std::set<int> nbrs (
		int, unsigned char, bool, const std::string&) const;
	void update (int, int, const std::string&, bool, float);
public:
	bool directed;
	float no_relationship;
	
	// constructors, destructor
	Graph (bool dir=true, float x=0);
	Graph (const Graph&);
	~Graph ();
	
	// operators
	Graph& operator=(const Graph&);
	bool operator<(const Graph&) const;
	
	// operations
	size_t size () const;
	std::set<int> nbrs (int, const std::string&) const;
	std::set<int> nbrs (int) const;
	std::set<int> nbrs_to (int, const std::string&) const;
	std::set<int> nbrs_to (int) const;
	std::set<int> nbrs_from (int, const std::string&) const;
	std::set<int> nbrs_from (int) const;
	std::set<int> vertices () const;
	std::set<std::string> keys () const;
	std::set<std::string> keys (int) const;
	std::set<std::string> keys (int, int) const;
	bool contains (int, int, const std::string&, bool) const;
	bool contains_dir (int, int, const std::string&) const;
	bool contains_undir (int, int, const std::string&) const;
	bool contains (int, int, const std::string&) const;
	bool contains_dir (int, int) const;
	bool contains_undir (int, int) const;
	bool contains (int, int) const;
	bool contains_dir (int) const;
	bool contains_undir (int) const;
	bool contains (int) const;
	float get (int, int, const std::string&) const;
	float get (int, int) const;
	void set (int, int, const std::string&, bool, float);
	void set (int, int, const std::string&, float);
	void set_dir (int, int, const std::string&, float);
	void set_undir (int, int, const std::string&, float);
	void set (int, int, float);
	void set_dir (int, int, float);
	void set_undir (int, int, float);
	void clear_dir (int, int);
	void clear_undir (int, int);
	void clear (int, int);
	void clear (int, int, const std::string&, bool);
	void clear_dir (int, int, const std::string&);
	void clear_undir (int, int, const std::string&);
	void clear (int, int, const std::string&);
	void clear_dir (int, const std::string&);
	void clear_undir (int, const std::string&);
	void clear (int, const std::string&);
	void clear (int);
	void clear (const std::string&);
	void clear ();
}; // Graph

const unsigned char Graph::UNDIRECTED = 0;
const unsigned char Graph::FROM = 1;
const unsigned char Graph::TO = 2;


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

Graph::Graph ( bool dir, float x )
: directed(dir), no_relationship(x)
{}

Graph::Graph ( const Graph& g )
: directed(g.directed)
{
	data = g.data;
}

Graph::~Graph () {}


// OPERATORS ////////////////////////////////////////////////////////

Graph& Graph::operator= ( const Graph& g )
{
	if ( this == &g ) return *this;
	directed = g.directed;
	data = g.data;
	return *this;
}

bool Graph::operator< ( const Graph& g ) const
{
	return data < g.data;
}


// OPERATIONS ///////////////////////////////////////////////////////

// Get the number of vertices represented in the graph.
size_t Graph::size () const
{
	return data.size();
}

// Get a set of neighbor IDs.
std::set<int> Graph::nbrs ( int i, unsigned char dir,
	bool limit_key, const std::string& key ) const
{
	std::set<int> out;
	
	std::map<int, std::map<int, std::map<std::string,
		std::pair<bool, float> > > >::const_iterator it = data.find(i);
	if ( it == data.end() ) return out;
	
	// vertex
	const std::map<int, std::map<
			std::string, std::pair<bool, float> > >& V =
		it->second;
	std::map<int, std::map<std::string,
			std::pair<bool, float> > >::const_iterator jt = V.begin();
	for ( ; jt != it->second.end(); ++jt ) {
		// neighbor
		int j = jt->first;
		const std::map<std::string, std::pair<bool, float> >& N =
			jt->second;
		
		if ( limit_key && N.find(key) == N.end() ) continue;
		
		// relationships
		std::map<std::string,
				std::pair<bool, float> >::const_iterator kt = N.begin();
		for ( ; kt != N.end(); ++kt ) {
			const std::string& k = kt->first;
			const std::pair<bool, float>& R = kt->second;
			if ( (!limit_key || k == key)
					&& (dir == UNDIRECTED
					|| (dir == FROM && R.first)
					|| (dir == TO && (!R.first
					|| data.find(j)->second.find(i)->second.find(k)->second.first))) )
				out.insert(jt->first);
		}
	}
	return out;
}

std::set<int> Graph::nbrs ( int i, const std::string& key ) const
{
	return nbrs(i, UNDIRECTED, true, key);
}

std::set<int> Graph::nbrs ( int i ) const
{
	return nbrs(i, UNDIRECTED, false, "");
}

std::set<int> Graph::nbrs_to ( int i, const std::string& key ) const
{
	return nbrs(i, TO, true, key);
}

std::set<int> Graph::nbrs_to ( int i ) const
{
	return nbrs(i, TO, false, "");
}

std::set<int> Graph::nbrs_from ( int i , const std::string& key ) const
{
	return nbrs(i, FROM, true, key);
}

std::set<int> Graph::nbrs_from ( int i ) const
{
	return nbrs(i, FROM, false, "");
}

// Returns a set of object IDs.
std::set<int> Graph::vertices () const
{
	std::set<int> out;
	std::map<int, std::map<int, std::map<std::string,
		std::pair<bool, float> > > >::const_iterator it =
		data.begin();
	for ( ; it != data.end(); ++it ) out.insert(it->first);
	return out;
}

// Returns all of the keys in the graph.
std::set<std::string> Graph::keys () const
{
	std::set<std::string> out;
	
	// vertices
	std::map<int, std::map<int, std::map<std::string,
			std::pair<bool, float> > > >::const_iterator it =
		data.begin();
	for ( ; it != data.end(); ++it ) {
		const std::map<int, std::map<std::string,
			std::pair<bool, float> > >& V = it->second;
		
		// neighbors
		std::map<int, std::map<std::string,
			std::pair<bool, float> > >::const_iterator jt = V.begin();
		for ( ; jt != V.end(); ++jt ) {
			const std::map<std::string, std::pair<bool, float> >& N =
				jt->second;
			
			// relationships
			std::map<std::string,
				std::pair<bool, float> >::const_iterator kt = N.begin();
			for ( ; kt != N.end(); ++kt ) out.insert(kt->first);
		}
	}
	
	return out;
}

// Returns all of the keys associated with the vertex.
std::set<std::string> Graph::keys ( int i ) const
{
	std::set<std::string> out;
	
	// vertex
	std::map<int, std::map<int, std::map<std::string,
			std::pair<bool, float> > > >::const_iterator it =
		data.find(i);
	const std::map<int, std::map<std::string,
		std::pair<bool, float> > >& V = it->second;
	
	// neighbors
	std::map<int, std::map<std::string,
		std::pair<bool, float> > >::const_iterator jt = V.begin();
	for ( ; jt != V.end(); ++jt ) {
		const std::map<std::string, std::pair<bool, float> >& N =
			jt->second;
		
		// relationships
		std::map<std::string,
			std::pair<bool, float> >::const_iterator kt = N.begin();
		for ( ; kt != N.end(); ++kt ) out.insert(kt->first);
	}
	
	return out;
}

// Returns a set of the relationship's keys or properties.
std::set<std::string> Graph::keys ( int i, int j ) const
{
	std::set<std::string> out;
	
	// vertex
	std::map<int, std::map<int, std::map<std::string,
			std::pair<bool, float> > > >::const_iterator t_i =
		data.find(i);
	if ( t_i == data.end() ) return out;
	const std::map<int, std::map<std::string,
		std::pair<bool, float> > >& V = t_i->second;
	
	// neighbor
	std::map<int, std::map<std::string,
			std::pair<bool, float> > >::const_iterator t_j =
		V.find(j);
	if ( t_j == V.end() ) return out;
	const std::map<std::string, std::pair<bool, float> >& N =
		t_j->second;
	
	// relationships
	std::map<std::string, std::pair<bool, float> >::const_iterator
		t_k = N.begin();
	for ( ; t_k != N.end(); ++t_k ) out.insert(t_k->first);
	
	return out;
}

// Returns true if the relationship exists for the given key.
bool Graph::contains (
	int i, int j, const std::string& key, bool undir ) const
{
	// vertex
	std::map<int, std::map<int, std::map<std::string,
			std::pair<bool, float> > > >::const_iterator t_i =
		data.find(i);
	if ( t_i == data.end() ) return false;
	const std::map<int, std::map<std::string,
			std::pair<bool, float> > >& V = t_i->second;
	
	// neighbor
	std::map<int, std::map<std::string,
			std::pair<bool, float> > >::const_iterator t_j =
		V.find(j);
	if ( t_j == V.end() ) return false;
	const std::map<std::string, std::pair<bool, float> >& N =
		t_j->second;
	
	// relationship
	std::map<std::string,
			std::pair<bool, float> >::const_iterator t_k =
		N.find(key);
	if ( t_k == N.end() ) return false;
	const std::pair<bool, float>& R = t_k->second;
	
	return undir || R.first;
}

bool Graph::contains_dir (
	int i, int j, const std::string& key ) const
{
	return contains(i, j, key, false);
}

bool Graph::contains_undir (
	int i, int j, const std::string& key ) const
{
	return contains(i, j, key, true);
}

bool Graph::contains ( int i, int j, const std::string& key ) const
{
	return contains(i, j, key, !directed);
}

// Returns true if a relationship exists between the given vertices.
bool Graph::contains_dir ( int i, int j ) const
{
	// vertex
	std::map<int, std::map<int, std::map<std::string,
			std::pair<bool, float> > > >::const_iterator it =
		data.find(i);
	if ( it == data.end() ) return false;
	const std::map<int, std::map<std::string,
			std::pair<bool, float> > >& V = it->second;
	
	// neighbor
	std::map<int, std::map<std::string,
			std::pair<bool, float> > >::const_iterator jt = V.find(j);
	if ( jt == V.end() ) return false;
	const std::map<std::string, std::pair<bool, float> >& N =
		jt->second;
	
	// relationships
	bool flag = false;
	std::map<std::string,
		std::pair<bool, float> >::const_iterator kt = N.begin();
	for ( ; !flag && kt != N.end(); ++kt ) flag = kt->second.first;
	
	return flag;
}

bool Graph::contains_undir ( int i, int j ) const
{
	// vertex
	std::map<int, std::map<int, std::map<std::string,
			std::pair<bool, float> > > >::const_iterator it =
		data.find(i);
	if ( it == data.end() ) return false;
	const std::map<int, std::map<std::string,
			std::pair<bool, float> > >& V = it->second;
	
	// neighbor
	std::map<int, std::map<std::string,
			std::pair<bool, float> > >::const_iterator jt = V.find(j);
	return jt != V.end() && jt->second.size() > 0;
}

bool Graph::contains ( int i, int j ) const
{
	if ( directed ) return contains_dir(i, j);
	return contains_undir(i, j);
}

// Returns true if the given vertex exists. If specifying directed
// (undir=false), only returns true if the vertex has an outgoing
// 'from' relationship.
bool Graph::contains_dir ( int i ) const
{
	// vertex
	std::map<int, std::map<int, std::map<std::string,
			std::pair<bool, float> > > >::const_iterator it =
		data.find(i);
	if ( it == data.end() ) return false;
	const std::map<int, std::map<std::string,
			std::pair<bool, float> > >& V = it->second;
	
	// neighbors
	bool flag = false;
	std::map<int, std::map<std::string,
		std::pair<bool, float> > >::const_iterator jt = V.begin();
	for ( ; !flag && jt != V.end(); ++jt ) {
		const std::map<std::string, std::pair<bool, float> >& N =
			jt->second;
		
		// relationships
		std::map<std::string,
			std::pair<bool, float> >::const_iterator kt = N.begin();
		for ( ; !flag && kt != N.end(); ++kt ) flag = kt->second.first;
	}
	
	return flag;
}

bool Graph::contains_undir ( int i ) const
{
	// vertex
	std::map<int, std::map<int, std::map<std::string,
			std::pair<bool, float> > > >::const_iterator it =
		data.find(i);
	return it != data.end() && it->second.size() > 0;
}

bool Graph::contains ( int i ) const
{
	if ( directed ) return contains_dir(i);
	return contains_undir(i);
}

// Returns the value of the relationship. If the relationship does
// not exist, returns the no_relationship value.
float Graph::get ( int i, int j, const std::string& key ) const
{
	// vertex
	std::map<int, std::map<int, std::map<std::string,
			std::pair<bool, float> > > >::const_iterator t_i =
		data.find(i);
	if ( t_i == data.end() ) return no_relationship;
	const std::map<int, std::map<std::string,
			std::pair<bool, float> > >& V = t_i->second;
	
	// neighbor
	std::map<int, std::map<std::string,
			std::pair<bool, float> > >::const_iterator t_j =
		V.find(j);
	if ( t_j == V.end() ) return no_relationship;
	const std::map<std::string, std::pair<bool, float> >& N =
		t_j->second;
	
	// relationship
	std::map<std::string,
			std::pair<bool, float> >::const_iterator t_k =
		N.find(key);
	if ( t_k == N.end() ) return no_relationship;
	const std::pair<bool, float>& R = t_k->second;
	
	if ( R.first == false ) return -1 * R.second;
	return R.second;
}

float Graph::get ( int i, int j ) const
{
	return get(i, j, "");
}

// Set the value of the given relationship. If it does not exist,
// creates it. If it already exists, overwrites it.
void Graph::update (
	int i, int j, const std::string& key, bool outward, float x )
{
	if ( data.find(i) == data.end() ) {
		data[i] = std::map<int, std::map<std::string,
			std::pair<bool, float> > > ();
	}
	if ( data[i].find(j) == data[i].end() ) {
		data[i][j] = std::map<std::string,
			std::pair<bool, float> > ();
	}
	data[i][j][key] = std::pair<bool, float> (outward, x);
}

// Undirected if undir == true
void Graph::set (
	int i, int j, const std::string& key, bool undir, float x )
{
	// check for no-relationship value
	if ( fabs(x - no_relationship) < 0.0000001 ) {
		if ( undir ) clear_undir(i, j, key);
		else clear_dir(i, j, key);
		return;
	}
	
	// create or update relationships
	update(i, j, key, true, x);
	if ( undir ) update(j, i, key, true, x);
	else if ( !contains(j, i, key) || !data[j][i][key].first )
		update(j, i, key, false, x);
}

void Graph::set (
	int i, int j, const std::string& key, float x )
{
	set(i, j, key, !directed, x);
}
void Graph::set_dir (
	int i, int j, const std::string& key, float x )
{
	set(i, j, key, false, x);
}
void Graph::set_undir (
	int i, int j, const std::string& key, float x )
{
	set(i, j, key, true, x);
}

void Graph::set ( int i, int j, float x )
{
	set(i, j, "", !directed, x);
}
void Graph::set_dir ( int i, int j, float x )
{
	set(i, j, "", false, x);
}
void Graph::set_undir ( int i, int j, float x )
{
	set(i, j, "", true, x);
}

// Remove the relationship(s) between the given vertices.
void Graph::clear_dir ( int i, int j )
{
	if ( !contains(i,j) ) return;
	
	std::map<std::string, std::pair<bool, float> >::iterator it;
	for ( it = data[i][j].begin(); it != data[i][j].end(); ) {
		if ( !it->second.first )
			++it;
		else if ( !data[j][i][it->first].first ) {
			data[j][i].erase(it->first);
			data[i][j].erase(it++);
		}
		else {
			it->second.first = false;
			it->second.second = data[j][i][it->first].second;
			++it;
		}
	}
	
	if ( data[i][j].size() == 0 ) data[i].erase(j);
	if ( data[j][i].size() == 0 ) data[j].erase(i);
	if ( data[i].size() == 0 ) data.erase(i);
	if ( data[j].size() == 0 ) data.erase(j);
}

void Graph::clear_undir ( int i, int j )
{
	if ( !contains_undir(i,j) ) return;
	
	data[i].erase(j);
	data[j].erase(i);
	if ( data[i].size() == 0 ) data.erase(i);
	if ( data[j].size() == 0 ) data.erase(j);
}

void Graph::clear ( int i, int j )
{
	if ( directed ) clear_dir(i,j);
	else clear_undir(i,j);
}

void Graph::clear (
	int i, int j, const std::string& key, bool undir )
{
	if ( !contains_undir(i,j,key) ) return;
	
	if ( undir ) {
		data[i][j].erase(key);
		data[j][i].erase(key);
	}
	
	else if ( !data[i][j][key].first ) {}
	else if ( !data[j][i][key].first ) {
		data[j][i].erase(key);
		data[i][j].erase(key);
	}
	else {
		data[i][j][key].first = false;
		data[i][j][key].second = data[j][i][key].second;
	}
	
	if ( data[i][j].size() == 0 ) data[i].erase(j);
	if ( data[j][i].size() == 0 ) data[j].erase(i);
	if ( data[i].size() == 0 ) data.erase(i);
	if ( data[j].size() == 0 ) data.erase(j);
}

void Graph::clear_dir ( int i, int j, const std::string& key )
{
	clear(i, j, key, false);
}

void Graph::clear_undir ( int i, int j, const std::string& key )
{
	clear(i, j, key, true);
}

void Graph::clear ( int i, int j, const std::string& key )
{
	clear(i, j, key, !directed);
}

// Remove relationships from vertex.
void Graph::clear_dir ( int i, const std::string& key )
{
	std::set<int> N = nbrs(i);
	std::set<int>::iterator it = N.begin();
	for ( ; it != N.end(); ++it ) clear_dir(i, *it, key);
}

void Graph::clear_undir ( int i, const std::string& key )
{
	std::set<int> N = nbrs(i, UNDIRECTED, true, key);
	std::set<int>::iterator it = N.begin();
	for ( ; it != N.end(); ++it ) clear_undir(i, *it, key);
}

void Graph::clear ( int i, const std::string& key )
{
	if ( directed ) clear_dir(i, key);
	else clear_undir(i, key);
}

// Remove vertex.
void Graph::clear ( int i )
{
	if ( !contains_undir(i) ) return;
	
	std::set<int> N = nbrs(i, UNDIRECTED, false, "");
	data.erase(i);
	std::set<int>::iterator it = N.begin();
	for ( ; it != N.end(); ++it ) {
		if ( data.find(*it) == data.end() ) continue;
		data[*it].erase(i);
		if ( data[*it].size() == 0 ) data.erase(*it);
	}
}

// Remove key.
void Graph::clear ( const std::string& key )
{
	// vertices
	std::map<int, std::map<int, std::map<
			std::string, std::pair<bool, float> > > >::iterator it =
		data.begin();
	for ( ; it != data.end(); ) {
		std::map<int, std::map<std::string,
			std::pair<bool, float> > >& V = it->second;
		
		// neighbors
		std::map<int, std::map<std::string,
			std::pair<bool, float> > >::iterator jt = V.begin();
		for ( ; jt != V.end(); ) {
			std::map<std::string, std::pair<bool, float> >& N =
				jt->second;
			
			// relationships
			std::map<std::string, std::pair<bool, float> >::iterator kt =
				N.begin();
			for ( ; kt != N.end(); ) {
				if ( kt->first == key ) N.erase(kt++);
				else ++kt;
			}
			
			if ( N.size() == 0 ) V.erase(jt++);
			else ++jt;
		}
		
		if ( V.size() == 0 ) data.erase(it++);
		else ++it;
	}
}

// Remove all relationships.
void Graph::clear ()
{
	data.clear();
}

} // namespace bygis

#endif // YOUNG_GIS_GRAPH_20221111