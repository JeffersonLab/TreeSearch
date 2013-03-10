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
#include "TRotation.h"
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
  const double phi_offset[4] = { 2.0, 4.0, 0.0, 0.0 };
  const string infile  = "/data/data/SoLID/tracking_1/db_gemc.dat";
  const string outfile = "db_solid.tracker.dat";
  const string prefix = "gemc.";
  const string out_prefix = "solid.tracker.";
  const char* proj_names[] = { "u", "v" };

  // Parameters for one GEM chamber (two readout coordinates)
  struct ValueSet_t {
    // Values read from source file
    double rmin, rmax, phi0, dphi, z, dz;
    double xangle, xpitch, yangle, ypitch;
    // Computed/converted quantities
    double phi, phioff;
    double uangle, ustart, upitch;  // angle is that of the axis!
    double vangle, vstart, vpitch;
    int nu, nv;                     // number of u/v strips
  } vals;

  struct GlobalValue_t {
    const char* key;
    double* var;
    double value;
    bool set;
    //    GlobalValue() : key(0), sourcekey(0), value(0.0), set(false) {}
  };

  GlobalValue_t per_plane_vals[] = {
    { "rmin",   &vals.rmin },
    { "rmax",   &vals.rmax },
    { "z",      &vals.z },
    { "dz",     &vals.dz },
    { "phioff", &vals.phioff },
    // { "nstrips", &uvals.nstrips },
    // { "strip.pos", &uvals.start },
    // { "strip.pitchs", &uvals.pitch },
    { 0 }
  };
  GlobalValue_t per_sector_vals[] = {
    { "phi",    &vals.phi },
    { 0 },
  };
  GlobalValue_t per_tracker_vals[] = {
    { "dphi",   &vals.dphi },
    { "dz",     &vals.dz },
    // { "nstrips", "nstrips" },
    // { "strip.pos", "strip.pos" },
    // { "strip.pitchs", "strip.pitch" },
    { 0 },
  };

  ifstream inp(infile.c_str());
  if( !inp ) {
    cerr << "Error opening " << infile << endl;
    exit(1);
  }
  for( int ip = 0; ip < nplanes; ++ip ) {
    for( int is = 0; is < nsect; ++is ) {
      int idx = is + nsect*ip; // linear index of this plane/sector combo
      ostringstream sector_prefix(prefix, ios_base::ate);
      sector_prefix << "gem" << idx+1 << ".";
      DBRequest request[] = {
	{ "r0",    &vals.rmin },
	{ "r1",    &vals.rmax },
	{ "phi0",  &vals.phi0 },
	{ "dphi",  &vals.dphi },
	{ "z0",    &vals.z },
	{ "depth", &vals.dz },
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
      if( vals.dphi < 0 or vals.dphi >= 90.0 ) {
	cerr << "Invalid opening angle dphi = " << vals.dphi << endl;
	exit(3);
      }
      if( vals.dz < 0 ) {
	cerr << "Invalid  z = " << vals.z << endl;
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

      // Convert parameters from libsolgem conventions to ours
      double phi2 = 0.5*vals.dphi;  // half opening angle
      vals.phioff = phi_offset[ip];
      vals.phi    = vals.phi0 + phi2 - vals.phioff;

      // Calculate strip start positions as in TSolGEMPlane::ReadGeometry
      double torad = TMath::DegToRad(), phi2rad = phi2 * torad;
      double xs = 0.5 * ( vals.rmax - vals.rmin * TMath::Cos(phi2rad) );
      double ys = vals.rmax * TMath::Sin(phi2rad);
#ifdef DEBUG
      cout << " xs/ys = " << xs << "/" << ys << endl;
#endif      
      TRotation plane_to_xstrip, plane_to_ystrip;
      plane_to_xstrip.RotateZ(-vals.xangle*torad);
      plane_to_ystrip.RotateZ(-vals.yangle*torad);
      TVector3 TR(xs,ys), BR(xs,-ys), TL(-xs,ys), BL(-xs,-ys);
      TVector3 C[4] = { TR, BR, TL, BL };
      int sminx = 1e9, sminy = 1e9, smaxx = -1e9, smaxy = -1e9;
      for( int i = 0; i < 4; ++i ) {
#ifdef DEBUG
	cout << " i = " <<  i << " "; C[i].Print();
#endif
	TVector3 vx = plane_to_xstrip * C[i];
	TVector3 vy = plane_to_ystrip * C[i];
	int sx = (int) (vx.X() / vals.xpitch);
	int sy = (int) (vy.X() / vals.ypitch);
	if( sx < sminx ) sminx = sx;
	if( sx > smaxx ) smaxx = sx;
	if( sy < sminy ) sminy = sy;
	if( sy > smaxy ) smaxy = sy;
      }
      vals.nu = smaxx - sminx + 1;
      vals.nv = smaxy - sminy + 1;
      vals.upitch = vals.xpitch;
      vals.vpitch = vals.ypitch;
      vals.ustart = -vals.nu * vals.upitch * 0.5;
      vals.vstart = -vals.nv * vals.vpitch * 0.5;
      // Correct angles for possible offset
      vals.uangle = vals.xangle - vals.phioff;
      vals.vangle = vals.yangle - vals.phioff;

      // Find common values
      // TODO

#ifdef DEBUG
      // Print results
      cout << sector_prefix.str() << endl;
      print_req( request );
      print_req( plane_request );
      cout << " phi/offset = " << vals.phi << "/" << vals.phioff << endl;

      cout << " u/v n/ang/start/pitch = "
	   << vals.nu << "/" << vals.uangle << "/" << vals.ustart << "/" << vals.upitch
	   << "  "
	   << vals.nv << "/" << vals.vangle << "/" << vals.vstart << "/" << vals.vpitch
	   << endl;

#endif

      // Write output
      ostringstream out_sector_prefix( out_prefix, ios_base::ate );
      out_sector_prefix << is << ".";
    }
  }
  inp.close();


  return 0;
}
