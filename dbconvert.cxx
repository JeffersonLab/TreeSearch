//*-- Author :    Ole Hansen  08-March-2013
//
// dbconvert.cxx
//
// Utility to convert a libsolgem-style database file to the format expected
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
#include <cctype>
#include <algorithm>
#include <functional>
#include <map>
#include <utility>

//#include "TDatime.h"
#include <unistd.h>    // for basename()

using namespace std;

#define ALL(c) (c).begin(), (c).end()

#define INFILE_DEFAULT "/data/data/SoLID/tracking_1/db_gemc.dat"
#define OUTFILE_DEFAULT "db_solid.tracker.dat"

// Command line parameters
static bool do_debug = false;
static string infile = INFILE_DEFAULT;
static string outfile = OUTFILE_DEFAULT;
static const char* prgname = "";

//-----------------------------------------------------------------------------
static string find_key( ifstream& inp, const string& key )
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

//-----------------------------------------------------------------------------
static int load_db( ifstream& inp, DBRequest* request, const string& prefix )
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
      if( do_debug )
	cout << "Read: " << key << " = " << *((double*)item->var) << endl;
    } else {
      cerr << "key \"" << key << "\" not found" << endl;
      return 2;
    }
    ++item;
  }
  return 0;
}

//-----------------------------------------------------------------------------
static
void print_req( const DBRequest* req )
{
  const DBRequest* item = req;
  while( item->name ) {
    cout << " " << item->name << " = " << *((double*)item->var) << endl;
    ++item;
  }
}

//-----------------------------------------------------------------------------
// Parameters for one GEM chamber (two readout coordinates)
struct ValueSet_t {
  // Values read from source file
  double rmin, rmax, phi0, dphi, z, dz;
  double xangle, xpitch, yangle, ypitch;
  // Computed/converted quantities
  double phi, phioff;
  double angle[2], start[2], pitch[2];  // angle is that of the axis!
  int nstrips[2];                       // number of u/v strips
  int iplane, isector;                  // plane index (1-4), sector number
};

//-----------------------------------------------------------------------------
// Helper functors for STL algos
struct BySectorThenPlane
  : public std::binary_function< const ValueSet_t&, const ValueSet_t&, bool >
{
  bool operator() ( const ValueSet_t& L, const ValueSet_t& R ) const
  { 
    return (L.isector != R.isector) ?
      (L.isector < R.isector) : (L.iplane < R.iplane);
  }
};

//-----------------------------------------------------------------------------
void usage()
{
  // Print usage message and exit
  cerr << "Usage: " << prgname << "[-hd] [-o outfile] [infile]" << endl;
  cerr << " Convert libsolgem database <infile> to TreeSearch-SoLID database"
       << " <outfile>" << endl;
  cerr << " -h: Print this help message" << endl;
  cerr << " -d: Output extensive debug information" << endl;
  cerr << " -o <outfile>: Write output to <outfile>. Default: "
       << OUTFILE_DEFAULT << endl;
  cerr << " <infile>: Read input from <infile>. Default: "
       << INFILE_DEFAULT << endl;
  exit(255);
}

//-----------------------------------------------------------------------------
void getargs( int argc, const char** argv )
{
  // Get command line parameters

  while (argc-- > 1) {
    const char *opt = *++argv;
    if (*opt == '-') {
      while (*++opt != '\0') {
	switch (*opt) {
	case 'h':
	  usage();
	  break;
	case 'd':
	  do_debug = true;
	  break;
	case 'o':
	  if (!*++opt) {
	    if (argc-- < 1)
	      usage();
	    opt = *++argv;
	  }
	  outfile = opt;
	  opt = "?";
	  break;
	default:
	  usage();
	}
      }
    } else {
      infile = *argv;
    }
  }
}

//-----------------------------------------------------------------------------
int main( int argc, const char** argv )
{
  const int nsect = 30;
  const int nplanes = 4, nproj = 2;
  const double phi_offset[4] = { 2.0, 4.0, 0.0, 0.0 };
  const string prefix = "gemc.";
  const string out_prefix = "solid.tracker.";
  //  const string out_prefix = "${DET}.";
  const string allsect_prefix = out_prefix + "${allsectors}.";
  const char* proj_name[2] = { "u", "v" };
  const string dashes =
    "#-----------------------------------------------------------------";

  // Parse command line
  prgname = basename(argv[0]);
  getargs(argc,argv);


  vector<ValueSet_t> values;

  // TODO: use this later to store common values
  // struct GlobalValue_t {
  //   const char* key;
  //   double* var;
  //   double value;
  //   bool set;
  //   //    GlobalValue() : key(0), sourcekey(0), value(0.0), set(false) {}
  // };

  // GlobalValue_t per_plane_vals[] = {
  //   { "rmin",   &vals.rmin },
  //   { "rmax",   &vals.rmax },
  //   { "z",      &vals.z },
  //   { "dz",     &vals.dz },
  //   { "phioff", &vals.phioff },
  //   // { "nstrips", &uvals.nstrips },
  //   // { "strip.pos", &uvals.start },
  //   // { "strip.pitchs", &uvals.pitch },
  //   { 0 }
  // };
  // GlobalValue_t per_sector_vals[] = {
  //   { "phi",    &vals.phi },
  //   { 0 },
  // };
  // GlobalValue_t per_tracker_vals[] = {
  //   { "dphi",   &vals.dphi },
  //   { "dz",     &vals.dz },
  //   // { "nstrips", "nstrips" },
  //   // { "strip.pos", "strip.pos" },
  //   // { "strip.pitchs", "strip.pitch" },
  //   { 0 },
  // };

  //==== Read input ====
  ifstream inp(infile.c_str());
  if( !inp ) {
    cerr << "Error opening " << infile << endl;
    exit(1);
  }

  int max_nstrips = 0;
  map<int,int> nstrip_map;

  for( int ip = 0; ip < nplanes; ++ip ) {
    for( int is = 0; is < nsect; ++is ) {
      ValueSet_t vals;
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

      // Calculate strip start positions in the same way as in
      // TSolGEMPlane::ReadGeometry
      double torad = TMath::DegToRad(), phi2rad = phi2 * torad;
      double xs = 0.5 * ( vals.rmax - vals.rmin * TMath::Cos(phi2rad) );
      double ys = vals.rmax * TMath::Sin(phi2rad);
      if( do_debug )
	cout << " xs/ys = " << xs << "/" << ys << endl;
      TRotation plane_to_xstrip, plane_to_ystrip;
      plane_to_xstrip.RotateZ(-vals.xangle*torad);
      plane_to_ystrip.RotateZ(-vals.yangle*torad);
      TVector3 TR(xs,ys), BR(xs,-ys), TL(-xs,ys), BL(-xs,-ys);
      TVector3 C[4] = { TR, BR, TL, BL };
      int sminx = 1e9, sminy = 1e9, smaxx = -1e9, smaxy = -1e9;
      for( int i = 0; i < 4; ++i ) {
	if( do_debug ) {
	  cout << " i = " <<  i << " ";
	  C[i].Print();
	}
	TVector3 vx = plane_to_xstrip * C[i];
	TVector3 vy = plane_to_ystrip * C[i];
	int sx = (int) (vx.X() / vals.xpitch);
	int sy = (int) (vy.X() / vals.ypitch);
	sminx = min(sx,sminx);
	smaxx = max(sx,smaxx);
	sminy = min(sy,sminy);
	smaxy = max(sy,smaxy);
      }
      vals.nstrips[0] = smaxx - sminx + 1;
      vals.nstrips[1] = smaxy - sminy + 1;
      vals.pitch[0] = vals.xpitch;
      vals.pitch[1] = vals.ypitch;
      for( int ij = 0; ij < nproj; ij++ ) {
	vals.start[ij] = -vals.nstrips[ij] * vals.pitch[ij] * 0.5;
	// For the tracking code, the strip positions should be the _center_
	// of the strips, not the lower edge
	vals.start[ij] += 0.5 * vals.pitch[ij];
      }
      // Correct angles for offset. The results are the projection axis angles.
      // For a given projection, they must be the same for all planes.
      vals.angle[0] = vals.xangle + vals.phioff;
      vals.angle[1] = vals.yangle + vals.phioff;

      // Save the number of strips of each readout for later use when writing
      // the detector maps
      for( int ij = 0; ij < nproj; ij++ ) {
	max_nstrips = max(max_nstrips,vals.nstrips[ij]);
	int jx = ij + nproj*( ip + nplanes*is );
	pair<map<int,int>::iterator,bool> itx =
	  nstrip_map.insert( make_pair(jx,vals.nstrips[ij]) );
	if( !itx.second ) {
	  cerr << "Duplicate index " << jx << " for sector/plane/proj = "
	       << is+1 << "/" << ip+1 << "/" << ij << endl;
	  cerr << "Bug - should never happen." << endl;
	  exit(6);
	}
      }

      // Remember who we are. Counting starts at 1, as the database is
      // meant for humans ...
      vals.iplane  = ip+1;
      vals.isector = is+1;

      // Save results
      values.push_back( vals );

      if( do_debug ) {
	// Display results for debugging
	//      cout << sector_prefix.str() << endl;
	print_req( request );
	print_req( plane_request );
	cout << " phi/offset = " << vals.phi << "/" << vals.phioff << endl;

	cout << " " << proj_name[0] << "/" << proj_name[1] << " n/ang/start/pitch = ";
	for( int i=0; i < 2; ++i ) {
	  cout << vals.nstrips[i] << "/" << vals.angle[i] << "/" << vals.start[i]
	       << "/" << vals.pitch[i];
	  if( i == 0 )
	    cout << "  ";
	}
	cout << endl;
      } // do_debug

    } // all sectors
  }   // all planes

  inp.close();

  // Find common values
  // TODO - this is a bit more ambitious than what I have time for right now

  // Check values for consistency (i.e., the values that MUST be common)
  vector<double> proj_angle(nproj,-1e10);
  vector<double> sect_phi(nsect,-1e10);
  for( vector<ValueSet_t>::size_type i = 0; i < values.size(); ++i ) {
    ValueSet_t& v = values[i];
    if( i == 0 ) {
      proj_angle[0] = v.angle[0];
      proj_angle[1] = v.angle[1];
    } else if( proj_angle[0] != v.angle[0] or proj_angle[1] != v.angle[1] ) {
      cerr << "Error: inconsistent projection angles = "
	   << v.angle[0] << "/" << v.angle[1] << " "
	   << "in sector/plane/index = " << v.isector << "/" << v.iplane-1
	   << "/" << (v.iplane-1)*nsect + v.isector
	   << endl
	   << "Expected " << proj_angle[0] << "/" << proj_angle[1] << endl;
      exit(4);
    }
    int is = v.isector-1;
    if( sect_phi[is] < -1e9 ) {
      sect_phi[is] = v.phi;
    } else if( sect_phi[is] != v.phi ) {
      cerr << "Error: inconsistent sector angle = " << v.phi
	   << " in sector " << is+1 << " at plane " << v.iplane
	   << endl
	   << "Expected " << sect_phi[is] << endl;
      exit(4);
    }      
  }

  //==== Write output ====
  ofstream outp( outfile.c_str(), ios_base::out|ios_base::trunc );
  if( !outp ) {
    cerr << " Error opening output file " << outfile << endl;
    exit(5);
  }

  // Header
  //TDatime now;
  outp << "# -*- mode: Text -*-" << endl
       << "#" << endl
       << "# Database for TreeSearch-SoLID" << endl
       << "#" << endl
       << "# Converted from " << infile << endl
       // << "# on " << now.AsString() << " by " << basename(argv[0]) << endl
       << endl;

  // Plane configuration
  outp << "# Plane configuration. One string of the all plane names." << endl;
  outp << endl;
  outp << allsect_prefix << "planeconfig = ";
  for( int ip = 0; ip < nplanes; ++ip ) {
    for( int ij = 0; ij < nproj; ++ij ) {
      outp << proj_name[ij] << ip+1;
      if( ip+1 != nplanes or ij+1 != nproj )
	outp << " ";
      else
	outp << endl;
    }
  }
  outp << endl;

  // Strip angles
  outp << "# Strip angles: angles of the normal to the strips, pointing" << endl;
  outp << "# along the direction of increasing strip number." << endl;
  outp << endl;
  for( int ij = 0; ij < nproj; ++ij ) {
    outp << allsect_prefix << proj_name[ij] << ".angle = "
	 << proj_angle[ij] << endl;
  }
  outp << endl;

  // Crate map. Lots of modules here
  int model = 6425;    // Dummy model for virtual APV25
  int nchan = 1280;    // Number of channels per module (arbitrary)
  int MAXSLOT = 27;    // Max slots per crate (from THaCrateMap.h)
  int modules_per_readout = max_nstrips/nchan+1;
  int modules_per_chamber = 2*modules_per_readout; // Modules needed per chamber
  int chambers_per_crate = (MAXSLOT/modules_per_chamber/nplanes)*nplanes;
  int slot_hi = chambers_per_crate*modules_per_chamber-1;
  if( do_debug ) {
    cout << "Crate map: modules_per_readout = " << modules_per_readout << endl;
    cout << "           modules_per_chamber = " << modules_per_chamber << endl;
    cout << "           chambers_per_crate  = " << chambers_per_crate << endl;
    cout << "           sectors_per_crate   = " << chambers_per_crate/nplanes
	 << endl;
  }
  const string spc = "                   ";
  outp << "# \"Crate map\". Specifies the overall DAQ module configuration." << endl;
  outp << "# The map can be common to all sectors. It's just a lookup table" << endl;
  outp << "# for (crate,slot) -> (model,nchan)" << endl;
  outp << "#" << endl;
  outp << "# Each row is:     crate slot_lo slot_hi model# nchan" << endl;
  outp << "#" << endl;
  outp << allsect_prefix << "cratemap = \\" << endl;

  // Remember that each sector must have its own virtual hardware.
  // However, all sectors can share the same crate map - it's just a catalog
  // of which modules are in which slot, regardless of whether or not we use
  // the slot.
  int maxcrates = TMath::CeilNint( (double)nplanes*nsect / 
				   (double)chambers_per_crate );
  for( int ic = 0; ic < maxcrates; ++ ic ) {
    outp << spc << ic << "    ";
    if( ic < 10 )  outp << " ";
    outp << 0 << "       " << slot_hi << "      ";
    if( slot_hi < 10 ) outp << " ";
    outp << model << "   " << nchan;
    if( ic+1 != maxcrates ) outp << " \\";
    outp << endl;
  }
  outp << endl;

  // Per-plane detector maps. Boyoboy
  outp << "# GEM detector maps" << endl;
  outp << "# " << nsect << " sectors * " << nplanes << " planes * "
       << nproj << " readout coordinates = " << nsect*nplanes*nproj
       << " detectors" << endl;
  outp << "#" << endl;
  for( int is = 0; is < nsect; ++is ) {
    for( int ip = 0; ip < nplanes; ++ip ) {
      for( int ij = 0; ij < nproj; ++ij ) {
	outp << out_prefix << is+1 << "." << proj_name[ij] << ip+1
	     << ".detmap = ";
	if( modules_per_readout > 1 )
	  outp << "\\" << endl;
	// Find the actual number of channels used by this particular
	// readout. We don't strictly need this, but setting a limit here
	// provides extra protection against decoder bugs.
	int jx = ij + nproj*( ip + nplanes*is );
	map<int,int>::iterator it = nstrip_map.find(jx);
	if( it == nstrip_map.end() ) {
	  cerr << "Error retrieving nstrips for sector/plane/proj = "
	       << is+1 << "/" << ip+1 << "/" << ij << endl;
	  cerr << "Bug - should never happen." << endl;
	  exit(6);
	}
	int the_nstrips = (*it).second;
	// Write detector map for this readout
	for( int im = 0; im < modules_per_readout; ++im ) {
	  int ix = im + modules_per_readout*( ij + nproj*( ip + nplanes*is ));
	  int cr = ix / (chambers_per_crate*modules_per_chamber);
	  int sl = ix - cr*chambers_per_crate*modules_per_chamber;
	  int lo = im*nchan;
	  int hi = min((im+1)*nchan,the_nstrips)-1;
	  if( modules_per_readout > 1 ) outp << spc;
	  outp << cr << " ";
	  if( cr < 10 )  outp << " ";
	  outp << sl << " ";
	  if( sl < 10 ) outp << " ";
	  outp << lo << " " << hi;
	  if( 2*(im+1) < modules_per_chamber )
	    outp << " \\";
	  outp << endl;
	}
      }
    }
  }
  outp << endl;

  // Phi angles of sectors
  outp << dashes << endl;
  outp << "#  Nominal phi angles of the sectors" << endl;
  outp << "#  Individual planes within each sector may have an angle offset" << endl;
  outp << dashes << endl;
  outp << endl;
  for( int is = 0; is < nsect; ++is ) {
    outp << out_prefix << is+1 << ".phi = " << sect_phi[is] << endl;
  }
  outp << endl;

  // Tracker configuration
  outp << dashes << endl;
  outp << "# Tracker feature configuration" << endl;
  outp << dashes << endl;
  outp << endl;

  outp << allsect_prefix << "search_depth = " << 10 << endl;
  // For decoder checkout/debugging
  outp << allsect_prefix << "maxthreads = " << 1 << endl;
  outp << allsect_prefix << "disable_tracking = " << 1 << endl;
  outp << endl;

  outp << dashes << endl;
  outp << "#  Global reconstruction parameters" << endl;
  outp << dashes << endl;
  outp << endl;

  outp << allsect_prefix << "do_adc_cut = " << 1 << endl;
  outp << allsect_prefix << "do_noise = " << 0 << endl;
  outp << allsect_prefix << "3d_ampcorr_maxmiss = " << 1 << endl;
  outp << allsect_prefix << "3d_ampcorr_nsigma = " << 0.18 << endl;
  outp << allsect_prefix << "3d_chi2_conflevel = " << 1e-6 << endl;
  outp << allsect_prefix << "3d_maxmiss = " << 2 << endl;
  outp << endl;

  // Projection parameters
  outp << dashes << endl;
  outp << "#  Global projection parameters" << endl;
  outp << dashes << endl;
  outp << endl;

  outp << allsect_prefix << "maxmiss = " << 1 << endl;  
  outp << allsect_prefix << "maxpat = " << 1000 << endl;  
  outp << allsect_prefix << "chi2_conflevel = " << 1e-6 << endl;  
  outp << endl;

  // Default plane parameters
  outp << dashes << endl;
  outp << "#  Global plane parameters" << endl;
  outp << dashes << endl;
  outp << endl;

  outp << out_prefix << "xp.res = " << 4e-5 << endl;  
  outp << out_prefix << "maxclustsiz = " << 4 << endl;  
  outp << out_prefix << "adc.min = " << 1 << endl;  
  outp << out_prefix << "split.frac = " << 0.1 << endl;  
  outp << out_prefix << "maxhits = " << 1000 << endl;  
  outp << out_prefix << "maxsamp = " << 3 << endl;  
  outp << out_prefix << "adc.sigma = " << 0.2 << endl;  
  outp << endl;

  // Per-plane parameters
  outp << dashes << endl;
  outp << "#   Plane-specific data" << endl;
  outp << "#   Detector maps are above" << endl;
  outp << dashes << endl;
  outp << endl;

  sort( ALL(values), BySectorThenPlane() );

  for( vector<ValueSet_t>::size_type i = 0; i < values.size(); ++i ) {
    ValueSet_t& v = values[i];
    for( int ij = 0; ij < nproj; ++ij ) {
      ostringstream the_plane_prefix( out_prefix, ios_base::ate );
      the_plane_prefix << v.isector << "." 
		       << proj_name[ij] << v.iplane << ".";
      const string& pfx = the_plane_prefix.str();

      outp << pfx << "nstrips = " << v.nstrips[ij] << endl;
      outp << pfx << "description = " << proj_name[ij] << v.iplane
	   << " plane in sector " << v.isector << endl;
      // The !ij trick only works when nproj == 2
      outp << pfx << "partner = " << proj_name[!ij] << v.iplane << endl;
      outp << pfx << "strip.pos = " << v.start[ij] << endl;
      outp << pfx << "strip.pitch = " << v.pitch[ij] << endl;
      // TODO: these may be redundant - test for defaults per plane or tracker
      outp << pfx << "rmin = " << v.rmin << endl;
      outp << pfx << "rmax = " << v.rmax << endl;
      outp << pfx << "z = " << v.z << endl;
      outp << pfx << "dz = " << v.dz << endl;
      outp << pfx << "dphi = " << v.dphi << endl;
      if( v.phioff != 0.0 )
	outp << pfx << "phioff = " << v.phioff << endl;
      outp << endl;
    }
  }

  return 0;
}
