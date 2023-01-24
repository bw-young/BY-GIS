/////////////////////////////////////////////////////////////////////
// Structures for handling variables in the general program or     //
// application environment.                                        //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 05/20/2022 - Brennan Young                                      //
// - created.                                                      //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_ENVIRONMENT_20220520
#define YOUNG_ENVIRONMENT_20220520

#include<map>
#include<string>

#include "../BYstdlib/ProgCounter.h"
#include "CompWindow.hpp"

namespace bygis { // Brennan Young GIS namespace

class Environment {
public:
	std::string name;
	
	bystd::ProgCounter appProgress;
	bystd::ProgCounter progress;
	
	CompWindow minWin;
	CompWindow maxWin;
	
	std::map<std::string, double> values;
	
	// constructors / destructors
	Environment (const std::string& nm="Application");
	Environment (const Environment&);
	~Environment ();
	
	// operators
	Environment& operator= (const Environment&);
}; // Environment


// CONSTRUCTORS / DESTRUCTORS ///////////////////////////////////////

Environment::Environment ( const std::string& nm )
: name(nm)
{
	appProgress.title = "elapsed time";
	appProgress.mode = bystd::ProgCounter::TIME;
	appProgress.interval = 1;
	
	progress.mode = bystd::ProgCounter::PERC
		| bystd::ProgCounter::SUBP
		| bystd::ProgCounter::TIME;
	progress.interval = 1;
	progress.level = 1;
	progress.reportLevel = 0; // doesn't print anything by default
	
	minWin = CompWindow(100);
	maxWin = CompWindow(5000);
}

Environment::Environment ( const Environment& env )
{
	name = env.name;
	appProgress = env.appProgress;
	progress = env.progress;
	minWin = env.minWin;
	maxWin = env.maxWin;
	values = env.values;
}

Environment::~Environment ()
{}


// OPERATORS ////////////////////////////////////////////////////////

Environment& Environment::operator= ( const Environment& env )
{
	if ( &env == this ) return *this;
	
	name = env.name;
	appProgress = env.appProgress;
	progress = env.progress;
	minWin = env.minWin;
	maxWin = env.maxWin;
	values = env.values;
	
	return *this;
}

} // bygis namespace

#endif // YOUNG_ENVIRONMENT_20220520