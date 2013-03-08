//*-- Author :    Ole Hansen  08-March-2013
//
// dbconvert.cxx
//
// Utility to convert a lobsolgem-style database file to the format expected
// by the SoLID::GEMTracker/GEMPlane classes.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "TMath.h"
#include "VarDef.h"
#include <cstdlib>
#include <algorithm>

#define DEBUG

using namespace std;

static
string find_key( ifstream& inp, const string& key )
{
  static const string empty("");
  string line;
  string::size_type keylen = key.size();
  inp.seekg(0); // could probably be more efficient, but it's fast enough
  while( getline(inp,line) ) {
    if( line.size() <= keylen )
      continue;
    if( line.compare(0,keylen,key) == 0 ) {
      if( keylen < line.size() ) {
	string::size_type pos = line.find_first_not_of(" \t=", keylen);
	if( pos != string::npos )
	  return line.substr(pos);
      }
      break;
    }
  }
  return empty;
}

static
int load_db( ifstream& inp, DBRequest* request, const string& prefix )
{
  DBRequest* item = request;
  while( item->name ) {
    ostringstream sn(prefix, ios_base::ate);
    sn << item->name;
    const string& key = sn.str();
    string val = find_key(inp,key);
    if( !val.empty() ) {
      istringstream sv(val);
      sv >> *((double*)item->var);
      if( !sv ) {
	cerr << "Error converting key/value = " << key << "/" << val << endl;
	return 1;
      }
#ifdef DEBUG
      cout << "Read: " << key << " = " << *((double*)item->var) << endl;
#endif
    } else {
      cerr << "key \"" << key << "\" not found" << endl;
      return 2;
    }
    ++item;
  }
  return 0;
}

#ifdef DEBUG
static
void print_req( const DBRequest* req )
{
  const DBRequest* item = req;
  while( item->name ) {
    cout << " " << item->name << " = " << *((double*)item->var) << endl;
    ++item;
  }
}
#endif

int main( int /*argc*/, char** /*argv*/ )
{
  const int nsect = 30;
  const int nplanes = 4;
  const double offs[4] = { 2, 4, 0, 0 };
  const vector<double> offsets( offs, offs+4 );
  const string infile  = "/data/data/SoLID/tracking_1/db_gemc.dat";
  const string outfile = "db_solid.tracker.dat";
  const string prefix = "gemc.";
  const string out_prefix = "solid.tracker.";
  const char* proj_names[] = { "u", "v" };

  struct GlobalValue_t {
    const char* key;
    const char* sourcekey;
    double value;
    bool set;
    //    GlobalValue() : key(0), sourcekey(0), value(0.0), set(false) {}
  };

  GlobalValue_t per_plane_vals[] = {
    { "rmin", "r0" },
    { "rmax", "r1" },
    { "z",    "z0" },
    { "dphi", "dphi" },
    { "dz", "depth" },
    { "phioff", "phioff" },
    { "nstrips", "nstrips" },
    { "strip.pos", "strip.pos" },
    { "strip.pitchs", "strip.pitch" },
    { 0 }
  };
  GlobalValue_t per_sector_vals[] = {
    { "phi", "phi0" },
    { 0 },
  };
  GlobalValue_t per_tracker_vals[] = {
    { "dphi", "dphi" },
    { "dz", "depth" },
    { "nstrips", "nstrips" },
    { "strip.pos", "strip.pos" },
    { "strip.pitchs", "strip.pitch" },
    { 0 },
  };

  ifstream inp(infile.c_str());
  if( !inp ) {
    cerr << "Error opening " << infile << endl;
    exit(1);
  }
  for( int ip = 0; ip < nplanes; ++ip ) {
    for( int is = 0; is < nsect; ++is ) {
      struct ValueSet {
	double rmin, rmax, phi0, dphi, z0, depth;
	double xangle, xpitch, yangle, ypitch;
      } vals;
      int idx = is + nsect*ip; // linear index of this plane/sector combo
      ostringstream sector_prefix(prefix, ios_base::ate);
      sector_prefix << "gem" << idx+1 << ".";
      DBRequest request[] = {
	{ "r0",    &vals.rmin },
	{ "r1",    &vals.rmax },
	{ "phi0",  &vals.phi0 },
	{ "dphi",  &vals.dphi },
	{ "z0",    &vals.z0 },
	{ "depth", &vals.depth },
	{ 0 }
      };
      int err = load_db( inp, request, sector_prefix.str() );
      if( err )
	exit(2);

      if( vals.rmin <= 0 or vals.rmax <= 0 or vals.rmax <= vals.rmin ) {
	cerr << "Invalid radii r0 = " << vals.rmin
	     << ", r1 = " << vals.rmax << endl;
	exit(3);
      }

      ostringstream plane_prefix(sector_prefix.str(), ios_base::ate);
      plane_prefix << "gem" << idx+1; // sic, the same thing again
      DBRequest plane_request[] = {
	{ "x.stripangle", &vals.xangle },
	{ "x.pitch",      &vals.xpitch },
	{ "y.stripangle", &vals.yangle },
	{ "y.pitch",      &vals.ypitch },
	{ 0 }
      };
      err = load_db( inp, plane_request, plane_prefix.str() );
      if( err )
	exit(2);

      if( vals.xpitch <= 0 or vals.ypitch <= 0 ) {
	cerr << "Illegal strip pitch xpitch = " << vals.xpitch
	     << ", ypitch = " << vals.ypitch << endl;
	exit(3);
      }

#ifdef DEBUG
      cout << sector_prefix.str() << endl;
      print_req( request );
      print_req( plane_request );
#endif

      ostringstream out_sector_prefix( out_prefix, ios_base::ate );
      out_sector_prefix << is << ".";
    }
  }
  inp.close();


  return 0;
}
