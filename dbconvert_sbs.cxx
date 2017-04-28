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
#include <iomanip>
#include <cassert>

#include "TDatime.h"
#include <unistd.h>    // for basename()

using namespace std;

#define ALL(c) (c).begin(), (c).end()

// These will be changed as this file will be dedicated to SBS DB
// TODO: configure all parameters with only the detector suffix as an input.
//#define DET_SUFFIX_DEFAULT  "fpp"
#define INFILE_DEFAULT "db_g4sbs_fpp.dat"
#define OUTFILE_DEFAULT "db_sbs_fpp.tracker.dat"

// Command line parameter defaults
static bool do_debug = false, do_dummies = false;
static string infile = INFILE_DEFAULT;//"db_g4sbs_"+DET_SUFFIX_DEFAULT+".dat";
static string outfile = OUTFILE_DEFAULT;//"db_sbs_"+DET_SUFFIX_DEFAULT+".tracker.dat";
static const char* prgname = "";
static const string spc = "                   ";

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
  double dmag, d0, xoff, dx, dy, thetaH, thetaV, depth;
  double xangle, xpitch, yangle, ypitch;
  // Computed/converted quantities
  double angle[2], start[2], pitch[2];  // angle is that of the axis!
  int nstrips[2];                       // number of u/v strips
  int iplane, isector;                  // plane index, sector number
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
  cerr << "Usage: " << prgname << "[-hdm] [-o outfile] [infile]" << endl;
  cerr << " Convert libsolgem database <infile> to TreeSearch-SoLID database"
       << " <outfile>" << endl;
  cerr << " -h: Print this help message" << endl;
  cerr << " -d: Output extensive debug information" << endl;
  cerr << " -m: Add extra dummy mode planes to emulate calorimeter" << endl;
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
	case 'm':
	  do_dummies = true;
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
static inline
void write_module( ostream& outp, int ic, int slot_hi, int model, int nchan )
{
  outp << spc << ic << "    ";
  if( ic < 10 )  outp << " ";
  outp << 0 << "       " << slot_hi << "      ";
  if( slot_hi < 10 ) outp << " ";
  outp << model << "   " << nchan;
}

//-----------------------------------------------------------------------------
static inline
void write_cslh( ostream& outp, int cr, int sl, int lo, int hi )
{
  outp << cr << " ";
  if( cr < 10 )  outp << " ";
  outp << sl << " ";
  if( sl < 10 ) outp << " ";
  outp << lo << " " << hi;
}

//-----------------------------------------------------------------------------
int main( int argc, const char** argv )
{
  // Parse command line
  prgname = basename(argv[0]);
  getargs(argc,argv);
  
  int crate_offset_v = 0;
  int nsect_v = 4;
  int nplanes_v = 10;
  double d0_dummy_v = 6.8;
  string prefix_v = "g4sbs_fpp.";
  string out_prefix_v = "sbs_fpp.tracker.";
  if(do_debug){
    cout << "Compare " << infile.c_str() << ": " << endl 
	 << " - with db_g4sbs_bbgem.dat : " << strcmp(infile.c_str(), "db_g4sbs_bbgem.dat") << endl
	 << " - with db_g4sbs_ft.dat : " << strcmp(infile.c_str(), "db_g4sbs_ft.dat") << endl
	 << " - with db_g4sbs_fpp.dat : " << strcmp(infile.c_str(), "db_g4sbs_fpp.dat") << endl;
  }
  bool validfile = false;
  if(strcmp(infile.c_str(), "db_g4sbs_ft.dat")==0){
    crate_offset_v = 4;
    nsect_v = 3;
    nplanes_v = 6;
    d0_dummy_v = 1.6;
    prefix_v = "g4sbs_ft.";
    out_prefix_v = "sbs_ft.tracker.";
    validfile = true;
  }
  if(strcmp(infile.c_str(), "db_g4sbs_bbgem.dat")==0){
    crate_offset_v = 0;
    nsect_v = 1;
    nplanes_v = 5;
    d0_dummy_v = 1.853320; 
    prefix_v = "g4sbs_bbgem.";
    out_prefix_v = "sbs_bbgem.tracker.";
    validfile = true;
  }
  if(strcmp(infile.c_str(), "db_g4sbs_fpp.dat")==0){
    validfile = true;
  }
  if(!validfile){
    cout << "Invalid input file: exit !" << endl;
    exit(-1);
  }
  
  // Important parameters
  const int crate_dummy = 6;
  const int crate_offset = crate_offset_v;
  const int nsect = nsect_v;
  const int nplanes = nplanes_v;
  const double d0_dummy = d0_dummy_v;
  const int nproj = 2;
  const string prefix = prefix_v;
  const string out_prefix = out_prefix_v;
  const string allsect_prefix = out_prefix + "${allsectors}.";
  const char* proj_name[2] = { "x", "y" };
  const string dashes =
    "#-----------------------------------------------------------------";


  vector<ValueSet_t> values;

  //==== Read input ====
  ifstream inp(infile.c_str());
  if( !inp ) {
    cerr << "Error opening " << infile << endl;
    exit(1);
  }

  int max_nstrips = 0;
  map<int,int> nstrip_map;

  int nplanes_eff = nplanes;
  if( do_dummies ) ++nplanes_eff;

  for( int is = 0; is < nsect; ++is ) {
    for( int ip = 0; ip < nplanes_eff; ++ip ) {
      ValueSet_t vals;
      DBRequest request[] = {
	{"dmag",        &vals.dmag   },
	{"d0",          &vals.d0     },
	{"xoffset",     &vals.xoff   },
	{"dx",          &vals.dx     },
	{"dy",          &vals.dy     },
	{"thetaH",      &vals.thetaH },
	{"thetaV",      &vals.thetaV },
	{"depth",       &vals.depth  },
	{ 0 }
      };
      DBRequest plane_request[] = {
	{ "x.stripangle", &vals.xangle },
	{ "x.pitch",      &vals.xpitch },
	{ "y.stripangle", &vals.yangle },
	{ "y.pitch",      &vals.ypitch },
	{ 0 }
      };
      if( ip < nplanes ) {
	// Regular GEM planes
	int idx = is*nplanes + ip; // linear index of this plane/sector combo
	ostringstream sector_prefix(prefix, ios_base::ate);
	sector_prefix << "gem" << idx << ".";
	
	if( do_debug )
	  cout << "plane_index/sector/plane (input) " 
	       << idx << "/" << is << "/" << ip << endl;
	
	int err = load_db( inp, request, sector_prefix.str() );
	if( err )
	  exit(2);

	if( vals.dmag < 0 or vals.d0 < 0) {
	  cerr << "Invalid  dmag = " << vals.dmag 
	       << ", d0 = " << vals.d0 << endl;
	  exit(3);
	}
	if( vals.dx <= 0 or vals.dy <= 0 ) {
	  cerr << "Invalid dimensions dx = " << vals.dx
	       << ", dy = " << vals.dy << endl;
	  exit(3);
	}
	if( vals.depth < 0 ) {
	  cerr << "Invalid  depth = " << vals.depth << endl;
	  exit(3);
	}
	
	ostringstream plane_prefix(sector_prefix.str(), ios_base::ate);
	plane_prefix << "gem" << idx; // sic, the same thing again
	err = load_db( inp, plane_request, plane_prefix.str() );
	if( err )
	  exit(2);

	if( vals.xpitch <= 0 or vals.ypitch <= 0 ) {
	  cerr << "Illegal strip pitch xpitch = " << vals.xpitch
	       << ", ypitch = " << vals.ypitch << endl;
	  exit(3);
	}
	
      }else {
	// Dummy planes
	assert( do_dummies );
	if(do_debug)cout << "Values size " << values.size() 
			 << " dummy plane index " << (is+1)*nplanes_eff-1
			 << " previous plane index " << (is+1)*nplanes_eff-2
			 << endl;
	assert( values.size() >= nplanes*(is+1) );
	// Use some values from the last GEM plane for the dummy plane as well
	const ValueSet_t& last_gem = values[(is+1)*nplanes_eff-2];
	vals.dmag = last_gem.dmag;
	vals.d0 = d0_dummy;
	vals.xoff = last_gem.xoff;
	vals.dx = last_gem.dx;
	vals.dy = last_gem.dy;
	vals.xangle = last_gem.xangle;
	vals.xpitch = 10.*last_gem.xpitch;
	vals.yangle = last_gem.yangle;
	vals.ypitch = 10.*last_gem.ypitch;
	if( do_debug ) {
	  cout << "Dummy plane in sector " << is+1 << endl;
	}
      }

      // Remember who we are. Counting starts at 1, as the database is
      // meant for humans ...
      vals.isector = is+1;
      vals.iplane  = ip+1;
      
      Int_t Nsx = (Int_t) (vals.dx / vals.xpitch);
      //The following 3 lines are to avoid to understimate the number of strips 
      //(hence the active area) because of stupid rounding issues at the 10^-12 level.
      if( (double)round(vals.dx / vals.xpitch)- (vals.dx / vals.xpitch) < 1.0e-10 ){
	Nsx = round(vals.dx / vals.xpitch);
      }
      Int_t Nsy = (Int_t) (vals.dy / vals.ypitch);
      //The following 3 lines are to avoid to understimate the number of strips 
      //(hence the active area) because of stupid rounding issues at the 10^-12 level.
      if( (double)round(vals.dy / vals.ypitch)- (vals.dy / vals.ypitch) < 1.0e-10 ){
	Nsy = round(vals.dy / vals.ypitch);
      }
      
      double xs = -Nsx/2*vals.xpitch;
      double ys = -Nsy/2*vals.ypitch;
      if( do_debug )
	cout << " xs/ys = " << xs << "/" << ys << endl;
      
      vals.nstrips[0] = Nsx;
      vals.nstrips[1] = Nsy;
      vals.start[0] = xs;
      vals.start[1] = ys;
      vals.pitch[0] = vals.xpitch;
      vals.pitch[1] = vals.ypitch;
      vals.angle[0] = vals.xangle;
      vals.angle[1] = vals.yangle;
 
      if( ip < nplanes ) {
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
      }
      
      // Save results
      values.push_back( vals );

      if( do_debug ) {
	// Display results for debugging
	//      cout << sector_prefix.str() << endl;
	print_req( request );
	print_req( plane_request );
	cout << " d0/xoffset = " << vals.d0 << "/" << vals.xoff << endl;

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
  vector<double> sect_dmag(nsect,-1e10);
  vector<double> sect_xoff(nsect,-1e10);
  vector<double> sect_thetaH(nsect,-1e10);
  vector<double> sect_thetaV(nsect,-1e10);
  for( vector<ValueSet_t>::size_type i = 0; i < values.size(); ++i ) {
    ValueSet_t& v = values[i];
    if( i == 0 ) {
      proj_angle[0] = v.angle[0];
      proj_angle[1] = v.angle[1];
    } else if( proj_angle[0] != v.angle[0] or proj_angle[1] != v.angle[1] ) {
      cerr << "Error: inconsistent projection angles = "
	   << v.angle[0] << "/" << v.angle[1] << " "
	   << " in sector/plane/index = " << v.isector-1 << "/" << v.iplane-1
	   << "/" << (v.isector-1)*nplanes + v.iplane-1
	   << endl
	   << "Expected " << proj_angle[0] << "/" << proj_angle[1] << endl;
      exit(4);
    }
    int is = v.isector-1;
    if( sect_dmag[is] < -1e9 ) {
      sect_dmag[is] = v.dmag;
    }else if( sect_dmag[is] != v.dmag ) {
      cerr << "Error: inconsistent sector dmag = " << v.dmag
	   << " in sector/plane/index = " << v.isector-1 << "/" << v.iplane-1
	   << "/" << (v.isector-1)*nplanes + v.iplane-1
	   << endl
    	   << "Expected " << sect_dmag[is] << endl;
      exit(4);
    }
    if( sect_xoff[is] < -1e9 ) {
      sect_xoff[is] = v.xoff;
    }else if( sect_xoff[is] != v.xoff ) {
      cerr << "Error: inconsistent sector xoff = " << v.xoff
	   << " in sector/plane/index = " << v.isector-1 << "/" << v.iplane-1
	   << "/" << (v.isector-1)*nplanes + v.iplane-1
	   << endl
	   << "Expected " << sect_xoff[is] << endl;
      exit(4);
    }
    if( sect_thetaH[is] < -1e9 ) {
      sect_thetaH[is] = v.thetaH;
    }else if( sect_thetaH[is] != v.thetaH ) {
      cerr << "Error: inconsistent sector thetaH = " << v.thetaH
	   << " in sector/plane/index = " << v.isector-1 << "/" << v.iplane-1
	   << "/" << (v.isector-1)*nplanes + v.iplane-1
	   << endl
    	   << "Expected " << sect_thetaH[is] << endl;
      exit(4);
    }
    if( sect_thetaV[is] < -1e9 ) {
      sect_thetaV[is] = v.thetaV;
    }else if( sect_thetaV[is] != v.thetaV ) {
      cerr << "Error: inconsistent sector thetaV = " << v.thetaV
	   << " in sector/plane/index = " << v.isector-1 << "/" << v.iplane-1
	   << "/" << (v.isector-1)*nplanes + v.iplane-1
	   << endl
    	   << "Expected " << sect_thetaV[is] << endl;
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

  // This database is for replay of Monte Carlo data, so say so
  outp << allsect_prefix << "MCdata = 1" << endl;
  outp << endl;

  // Plane configuration
  outp << "# Plane configuration. One string of the all plane names." << endl;
  outp << endl;
  outp << allsect_prefix << "planeconfig = ";
  for( int ip = 0; ip < nplanes_eff; ++ip ) {
    for( int ij = 0; ij < nproj; ++ij ) {
      outp << proj_name[ij];
      if( ip < nplanes )
	outp << ip+1;
      else
	outp << "d";
      if( ip+1 != nplanes_eff or ij+1 != nproj )
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
  int nchan = 2048;    // Number of channels per module (arbitrary) 
  // -> We put the real stuff: 1 MPD = 16*128 = 2048 channels
  int MAXSLOT = 30;    // Max slots per crate (from THaCrateMap.h)
  int modules_per_readout = max_nstrips/nchan+1;
  int modules_per_chamber = nproj*modules_per_readout; // Modules needed per chamber
  int chambers_per_crate = (MAXSLOT/modules_per_chamber/nplanes)*nplanes;
  int slot_hi = chambers_per_crate*modules_per_chamber-1;
  if( do_debug ) {
    cout << "Crate map: modules_per_readout = " << modules_per_readout << endl;
    cout << "           modules_per_chamber = " << modules_per_chamber << endl;
    cout << "           chambers_per_crate  = " << chambers_per_crate << endl;
    cout << "           sectors_per_crate   = " << chambers_per_crate/nplanes
	 << endl;
  }
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
  for( int ic = crate_offset; ic < crate_offset+maxcrates; ++ ic ) {
    write_module( outp, ic, slot_hi, model, nchan );
    if( ic+1 != maxcrates or do_dummies ) outp << " \\";
    outp << endl;
  }
  if( do_dummies ) {
    write_module( outp, crate_dummy, nproj-1, model, crate_offset+nsect );
    outp << endl;
  }
  outp << endl;
  
  if(do_debug)
    cout << "// Per-plane detector maps. Boyoboy " << endl;
  
  int crate_g = crate_offset;
  int slot_g = 0;
  
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
	  if(lo<=hi){
	    if( modules_per_readout > 1 ) outp << spc;
	    //write_cslh( outp, crate_offset+cr, sl, lo, hi );
	    write_cslh( outp, crate_offset+cr, sl, lo%nchan, hi%nchan );
	    //write_cslh( outp, crate_g, slot_g, lo, hi );
	    slot_g++;
	    if(slot_g==30){
	      crate_g++;
	      slot_g = 0;
	    }
	  }
	  if( 2*(im+1) < modules_per_chamber && (im+1)*nchan < the_nstrips){
	    outp << " \\";
	  }
	  if(lo<=hi || 2*(im+1)>modules_per_chamber)
	    outp << endl;
	}
      }
    }
  }
  outp << endl;

  if(do_debug)
    cout << "do dummies " << endl;

  if( do_dummies) {
    outp << "# Dummy GEM planes recording emulated calorimeter hits" << endl;
    outp << "#" << endl;
    for( int is = 0; is < nsect; ++is ) {
      for( int ij = 0; ij < nproj; ++ij ) {
	// E.g.: solid.tracker.1.ud.detmap =
	outp << out_prefix << is+1 << "." << proj_name[ij] << "d.detmap = ";
	// The crate for the dummy calorimeter planes follows right after those
	// of the regular GEM planes. u and v each get one module (slot).
	// The channel number in each slot is the sector number. Hits in each
	// channel represent measured coordinates.
	write_cslh( outp, crate_dummy, ij, crate_offset+is, crate_offset+is );
	//FIXME: it is a little bit "ad-hoc" for the moment. Do this better when I have time...
	outp << endl;
      }
    }
    outp << endl;
  }

  if(do_debug)
    cout << "# X offsets " << endl;

  // Phi angles of sectors
  outp << dashes << endl;
  outp << "#  X offset of the sectors (in transport coordinates)" << endl;
  outp << "#  and common position parameters (D_magnet, theta_H, theta_V)" << endl;
  outp << dashes << endl;
  outp << endl;
  for( int is = 0; is < nsect; ++is ) {
    outp << out_prefix << is+1 << ".dmag = " << sect_dmag[is] << endl;
    outp << out_prefix << is+1 << ".xoff = " << sect_xoff[is] << endl;
    outp << out_prefix << is+1 << ".thetaH = " << sect_thetaH[is] << endl;
    outp << out_prefix << is+1 << ".thetaV = " << sect_thetaV[is] << endl;
  }
  outp << endl;
  
  if(do_debug)
    cout << "# Tracker feature configuration " << endl;

  // Tracker configuration
  outp << dashes << endl;
  outp << "# Tracker feature configuration" << endl;
  outp << dashes << endl;
  outp << endl;

  outp << allsect_prefix << "search_depth = " << 10 << endl;
  // For decoder checkout/debugging
  outp << allsect_prefix << "maxthreads = " << 1 << endl;
  outp << allsect_prefix << "disable_tracking = " << 0 << endl;
  outp << allsect_prefix << "disable_finetrack = " << 1 << endl;
  outp << allsect_prefix << "3d_disable_chi2 = " << 0 << endl;
  outp << endl;

  if(do_debug)
    cout << "#  Global reconstruction parameters" << endl;

  outp << dashes << endl;
  outp << "#  Global reconstruction parameters" << endl;
  outp << dashes << endl;
  outp << endl;

  outp << allsect_prefix << "3d_ampcorr_maxmiss = " << 3 << endl;
  outp << allsect_prefix << "3d_ampcorr_nsigma = " << 0.18 << endl;
  outp << allsect_prefix << "3d_chi2_conflevel = " << 1e-6 << endl;
  outp << allsect_prefix << "3d_maxmiss = " << 2 << endl;
  outp << allsect_prefix << "proj_to_z0 = " << 0 << endl;
  outp << endl;

  if(do_debug)
    cout << "#  Global projection parameters" << endl;
  
  // Projection parameters
  outp << dashes << endl;
  outp << "#  Global projection parameters" << endl;
  outp << dashes << endl;
  outp << endl;

  outp << allsect_prefix << "maxmiss = " << 2 << endl;
  outp << allsect_prefix << "maxpat = " << 100000 << endl;
  outp << allsect_prefix << "chi2_conflevel = " << 1e-6 << endl;
  outp << allsect_prefix << "disable_chi2 = " << 0 << endl;
  outp << endl;

  if(do_debug)
    cout << "#  Global plane parameters" << endl;
  
  // Default plane parameters
  outp << dashes << endl;
  outp << "#  Global plane parameters" << endl;
  outp << dashes << endl;
  outp << endl;

  outp << out_prefix << "xp.res = " << 9e-5 << endl;
  outp << out_prefix << "maxclustsiz = " << 4 << endl;
  outp << out_prefix << "adc.min = " << 500 << endl;
  outp << out_prefix << "split.frac = " << 0.1 << endl;
  outp << out_prefix << "maxhits = " << 1000 << endl;
  outp << out_prefix << "maxsamp = " << 3 << endl;
  outp << out_prefix << "adc.sigma = " << 0.2 << endl;
  outp << out_prefix << "do_noise = " << 0 << endl;
  outp << out_prefix << "check_pulse_shape = " << 1 << endl;
  outp << out_prefix << "do_histos = " << 0 << endl;
  outp << endl;

  if(do_debug)
    cout << "#   Plane-specific data" << endl;

  // Per-plane parameters
  outp << dashes << endl;
  outp << "#   Plane-specific data" << endl;
  outp << "#   Detector maps are above" << endl;
  outp << dashes << endl;
  outp << endl;

  if(do_debug)
    cout << "sort ()" << endl;
  
  sort( ALL(values), BySectorThenPlane() );

  for( vector<ValueSet_t>::size_type i = 0; i < values.size(); ++i ) {
    ValueSet_t& v = values[i];
    bool dummy = (v.iplane == nplanes+1);
    assert( not dummy or do_dummies );

    for( int ij = 0; ij < nproj; ++ij ) {
      ostringstream the_plane_prefix( out_prefix, ios_base::ate );
      the_plane_prefix << v.isector << "." << proj_name[ij];
      if( not dummy )
	the_plane_prefix << v.iplane << ".";
      else
	the_plane_prefix << "d.";
      const string& pfx = the_plane_prefix.str();

      outp << pfx << "description = " << proj_name[ij];
      if( not dummy )
	outp << v.iplane;
      else
	outp << " dummy";
      outp << " plane in sector " << v.isector << endl;
      if( dummy ) {
	// Extra parameters for dummy planes
	outp << pfx << "dummy = "  << 1 << endl;
	outp << pfx << "xp.res = " << 1e-2 << endl;    // 1cm resolution
      }
      outp << pfx << "nstrips = " << v.nstrips[ij] << endl;
      // The !ij trick only works when nproj == 2
      outp << pfx << "partner = " << proj_name[!ij];
      if( not dummy )
	outp << v.iplane;
      else
	outp << "d";
      outp << endl;
      outp << pfx << "strip.pos = " << v.start[ij] << endl;
      outp << pfx << "strip.pitch = " << v.pitch[ij] << endl;
      // TODO: these may be redundant - test for defaults per plane or tracker
      outp << pfx << "d0 = " << setprecision(7) << v.d0 << endl;
      outp << pfx << "dx = " << v.dx << endl;
      outp << pfx << "dy = " << v.dy << endl;
      outp << pfx << "dz = " << v.depth << endl;
      outp << endl;
    }
  }
  if(do_debug)
    cout << "return" << endl;

  return 0;
}
