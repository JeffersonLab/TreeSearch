//*-- Author :    Ole Hansen, Jefferson Lab   31-Jan-2013

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Tracker                                                       //
//                                                                           //
// Abstract base class for set of tracker planes analyzed using TreeSearch.  //
// TreeSearch is the recusrsive template matching algorithm described in     //
// DellOrso et al., NIM A 287, 436 (1990).                                   //
//                                                                           //
// Examples of actual implementations are TreeSearch::MWDC for sets of       //
// horizontal drift chambers and TreeSearch::GEMTracker for sets of GEM      //
// tracker planes.                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Tracker.h"
#include "Plane.h"
#include "Projection.h"
#include "Road.h"
#include "Helper.h"

#include "THaDetMap.h"
#include "THaTrack.h"

#include "TString.h"
#include "TMath.h"
#include "THashTable.h"
#include "TVector2.h"
#include "TDecompChol.h"
#include "TSystem.h"
#include "TThread.h"
#include "TCondition.h"
#include "TMutex.h"
#include "TBits.h"
#include "TClass.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <map>
#include <string>
#include <stdexcept>
#include <cstring>   // for memset

#ifdef TESTCODE
#include "TStopwatch.h"
#endif
#ifdef MCDATA
 #include "Hit.h"
#include "GEMHit.h"
#endif

using namespace std;
using namespace Podd;

typedef string::size_type ssiz_t;

namespace {

// Helper classes for describing a DAQ hardware module and using it with a
// THashTable. CSpair can be used to find elements in the THashTable w/o
// creating large dummy objects.
//TODO: replace in favor or retrieving the information from db_cratemap
class CSpair : public TObject {
public:
  CSpair( UShort_t crate, UShort_t slot ) : fCrate(crate), fSlot(slot) {}
  virtual ULong_t Hash() const {
    UInt_t cs = (static_cast<UInt_t>(fCrate)<<16) + static_cast<UInt_t>(fSlot);
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,16,0)
    return TString::Hash( &cs, sizeof(cs) );
#else
    return TMath::Hash( &cs, sizeof(cs) );
#endif
  }
  virtual Bool_t IsEqual( const TObject* obj ) const {
    const CSpair* m;
    if( !obj or !(m = dynamic_cast<const CSpair*>(obj)) ) return kFALSE;
    return ( fCrate == m->fCrate and fSlot == m->fSlot );
  }
  UShort_t  fCrate;
  UShort_t  fSlot;
};
class DAQmodule : public CSpair {
public:
  DAQmodule( UShort_t crate, UShort_t slot, UInt_t model, UInt_t nchan )
    : CSpair(crate, slot), fModel(model), fNchan(nchan),
      fHasResolution(false) {}
  DAQmodule( UShort_t crate, UShort_t slot, UInt_t model, UInt_t nchan,
	     Double_t res )
    : CSpair(crate, slot), fModel(model), fNchan(nchan),
      fResolution(res*1e-12) /* NB: resolution in ps! */,
      fHasResolution(true) {}
  virtual ~DAQmodule() {}
  virtual void Copy( TObject& obj ) const {
    TObject::Copy(obj);
    assert( dynamic_cast<DAQmodule*>(&obj) );
    DAQmodule* m = static_cast<DAQmodule*>(&obj);
    m->fCrate = fCrate; m->fSlot = fSlot; m->fModel = fModel;
    m->fNchan = fNchan; m->fResolution = fResolution;
    m->fHasResolution = fHasResolution;
  }
  virtual void Print( Option_t* ) const {
    cout << "DAQmodule: "
	 << " crate = " << fCrate
	 << " slot = "  << fSlot
	 << " model = " << fModel
	 << " nchan = " << fNchan;
    if( fHasResolution )
      cout << " res = "   << fResolution*1e12 << " ps";
    cout << endl;
  }
  UInt_t    fModel;
  UInt_t    fNchan;
  Double_t  fResolution;
  bool      fHasResolution;
};

///////////////////////////////////////////////////////////////////////////////

} // end namespace

namespace TreeSearch {

typedef vector<Plane*>::size_type vrsiz_t;
typedef vector<Plane*>::iterator  vriter_t;
typedef vector<Projection*>::size_type vpsiz_t;
typedef vector<Projection*>::iterator  vpiter_t;
typedef vector<vector<Int_t> >::iterator vviter_t;

// Global constant indicating invalid/uninitialized data
const Double_t kBig = 1e38;

#define ALL(c) (c).begin(), (c).end()

// Default value for minimum difference between all projection angles
static const Double_t kMinProjAngleDiff = 5.0 * TMath::DegToRad();

#ifdef MCDATA
// Reconstruction status bit numbers, for evaluating tracking with MC data
enum EReconBits {
  kHitsFound = 0       // Sufficient number of hits found
};
#endif

//_____________________________________________________________________________
// Structures we want to put into STL containers

// A projection angle along with a pointer to the corresponding projection,
// sortable by angle
struct ProjAngle_t {
  ProjAngle_t( const Projection* p ) : m_proj(p) {
    assert(p);
    m_angle = p->GetAngle();
    assert( TMath::Abs(m_angle) <= TMath::Pi() ); // ensured in Projection::SetAngle
  }
  Double_t angle( bool flip = false ) const {
    if( flip )
      return m_angle - TMath::Sign(TMath::Pi(), m_angle);
    else
      return m_angle;
  }
  bool operator<( const ProjAngle_t & rhs ) const {
    return (m_angle < rhs.m_angle);
  }
  Double_t m_angle;
  const Projection* m_proj;
};

//_____________________________________________________________________________
// Support classes for data-processing threads

struct TrackThread {
  Projection*  proj;    // Projection to be processed
  Int_t*       status;  // Return status (shared)
  UInt_t*      running; // Bitfield indicating threads to wait for (shared)
  TMutex*      start_m; // Mutex for start condition
  TMutex*      done_m;  // Mutex for done condition
  TCondition*  start;   // Start condition of tracking threads
  TCondition*  done;    // Condition indicating all tracking threads done
  TThread*     thread;  // The actual thread running with these arguments
  TrackThread()
    : proj(0), status(0), running(0), start(0), done(0), thread(0) {}
  ~TrackThread() {
    if( thread ) {
      TThread::Delete(thread);
      delete thread;
    }
  }
};

//_____________________________________________________________________________
class ThreadCtrl {
public:
  ThreadCtrl( const vector<Projection*>& vp )
    : fTrackStatus(0), fTrackToDo(0),
      fTrackStartM(new TMutex), fTrackDoneM(new TMutex),
      fTrackStart(new TCondition(fTrackStartM)),
      fTrackDone(new TCondition(fTrackDoneM))
  {
    if( !vp.empty() ) {
      fTrack.reserve( vp.size() );
      fTrackStartM->Lock();
      for( vpsiz_t k = 0; k < vp.size(); ++k ) {
	TThread* t = AddTrackThread( vp[k] );
	t->Run();
      }
      fTrackStartM->UnLock();
    }
  }

  //___________________________________________________________________________
  ~ThreadCtrl()
  {
    // Terminate all running tracking threads
    if( !fTrack.empty() ) {
      fTrackDoneM->Lock();
      fTrackStartM->Lock();
      fTrackToDo = 0;
      for( vector<TrackThread>::iterator it = fTrack.begin(); it !=
	     fTrack.end(); ++it ) {
	TThread* th = (*it).thread;
	if( th and th->GetState() == TThread::kRunningState )
	  SETBIT(fTrackToDo, (*it).proj->GetType());
      }
      // MSB is termination flag
      SETBIT(fTrackToDo, kThreadTerminateBit);
      fTrackStart->Broadcast();
      fTrackStartM->UnLock();
      while( true ) {
#ifndef NDEBUG
	Int_t ret =
#endif
	  fTrackDone->Wait();
	assert( ret == 0 );
	if( fTrackToDo == BIT(kThreadTerminateBit) )
	  break;
      }
      fTrackDoneM->UnLock();
    }
    // At this point, all threads should have released their condition
    // waits and mutexes, so clean up is safe
    delete fTrackDone;
    delete fTrackStart;
    delete fTrackDoneM;
    delete fTrackStartM;
  }

  //___________________________________________________________________________
  Int_t Run( UInt_t maxthreads )
  {
    // Run our threads, at most maxthreads at a time
    if( maxthreads == 0 or fTrack.empty() )
      return 0;
    if( maxthreads > kThreadTerminateBit )
      maxthreads = kThreadTerminateBit;
    fTrackDoneM->Lock();
    fTrackStatus = 0;
    vector<TrackThread>::iterator it = fTrack.begin();
    while( it != fTrack.end() ) {
      fTrackStartM->Lock();
      fTrackToDo = 0;
      UInt_t n = 0;
      while( n < maxthreads and it != fTrack.end() ) {
	SETBIT(fTrackToDo, (*it).proj->GetType());
	++n;
	++it;
      }
      // Start tracking within the threads
      fTrackStart->Broadcast();
      fTrackStartM->UnLock();
      while( true ) {
	// Wait for end of processing
#ifndef NDEBUG
	Int_t ret =
#endif
	  fTrackDone->Wait();
	assert( ret == 0 );
	if( fTrackToDo == 0 )
	  break;
      }
    }
    fTrackDoneM->UnLock();
    return fTrackStatus;
  }

  //___________________________________________________________________________
  static void DoTrack( void* ptr )
  {
    TrackThread* arg = reinterpret_cast<TrackThread*>(ptr);

    //  TThread::SetCancelOn();
    bool terminate = false;
    arg->start_m->Lock();
    while( !terminate ) {

      // Wait for start condition
      while( true ) {
#ifndef NDEBUG
	Int_t ret =
#endif
	  arg->start->Wait(); // unlocks arg->start_m
	assert( ret == 0 );
	terminate = TESTBIT(*arg->running, kThreadTerminateBit);
	if( TESTBIT(*arg->running, arg->proj->GetType()) or terminate )
	  break;
      }
      arg->start_m->UnLock();

      Int_t nrd = 0;
      if( !terminate )
	// Process this event
	nrd = arg->proj->Track();

      // Ensure we're the only one modifying/testing thread
      // control data (arg->running)
      arg->done_m->Lock();
      // Set error flag, if necessary
      if( nrd < 0 ) {
	*arg->status = 1;
      }
      // Clear the bit for this projection in the status bitfield
      CLRBIT( *arg->running, arg->proj->GetType() );
      // If all bits are zero, all threads have finished processing,
      // in which case we wake up the main thread
      if( (*arg->running & ~BIT(kThreadTerminateBit)) == 0 )
	arg->done->Signal();

      if( !terminate )
	// Ensure that we enter Wait() before the main thread can send the
	// next Broadcast(). This must come before unlocking arg->done_m,
	// or else we have a race condition in the main thread.
	arg->start_m->Lock();

      arg->done_m->UnLock();
    }
  }

  static UInt_t GetMaxThreads() { return kThreadTerminateBit; }

private:
  static const UInt_t kThreadTerminateBit = 8*sizeof(UInt_t)-1; // = 31

  TThread* AddTrackThread( Projection* proj ) {
    assert( proj );
    fTrack.push_back( TrackThread() );
    TrackThread* t = &fTrack.back();
    t->proj    = proj;
    t->status  = &fTrackStatus;
    t->running = &fTrackToDo;
    t->start_m = fTrackStartM;
    t->done_m  = fTrackDoneM;
    t->start   = fTrackStart;
    t->done    = fTrackDone;
    string tn;
    tn = "trk_"; tn.append( proj->GetName() );
    t->thread  = new TThread( tn.c_str(), DoTrack, (void*)t );
    return t->thread;
  }
  vector<TrackThread>  fTrack;        // Tracking thread descriptors
  Int_t                fTrackStatus;  // Common status variable
  UInt_t               fTrackToDo;    // Bitfield of threads to wait for
  TMutex*              fTrackStartM;  // Mutex for start condition
  TMutex*              fTrackDoneM;   // Mutex for done condition
  TCondition*          fTrackStart;   // Start condition
  TCondition*          fTrackDone;    // Finish condition
};

//====================== Tracker class ========================================

//_____________________________________________________________________________
Tracker::Tracker( const char* name, const char* desc, THaApparatus* app )
  : THaTrackingDetector(name,desc,app), fCrateMap(0),
    fMinProjAngleDiff(kMinProjAngleDiff), fIsRotated(false),
    fAllPartnered(false), fMaxThreads(1), fThreads(0),
    fMinReqProj(3), f3dMatchvalScalefact(1), f3dMatchCut(0),
    fMinNdof(1), fTrkStat(kTrackOK),
    fNcombos(0), fN3dFits(0), fEvNum(0),
    t_track(0), t_3dmatch(0), t_3dfit(0), t_coarse(0)
#ifdef MCDATA
  , fMCDecoder(0), fMCPointUpdater(0), fChecked(false)
#endif
  , fDBmaxmiss(-1), fDBconf_level(1e-9)
{
  // Constructor

  SetBit(kProjTrackToZ0);
}

//_____________________________________________________________________________
Tracker::~Tracker()
{
  // Destructor. Delete objects & subdetectors and unregister variables
  if (fIsSetup)
    RemoveVariables();

  delete fThreads;
  if( fMaxThreads > 1 )
    gSystem->Unload("libThread");

  DeleteContainer( fPlanes );
  DeleteContainer( fProj );

#ifdef MCDATA
  delete fMCPointUpdater;
#endif
}

//_____________________________________________________________________________
Int_t Tracker::Begin( THaRunBase* run )
{
#ifdef TESTCODE
  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    fPlanes[iplane]->Begin(run);
#endif
  return 0;
}

//_____________________________________________________________________________
void Tracker::Clear( Option_t* opt )
{
  // Clear event-by-event data, including those of the planes and projections
  THaTrackingDetector::Clear(opt);

  // Clear the planes and projections, but only if we're not called from Init()
  if( !opt or *opt != 'I' ) {
    for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
      fPlanes[iplane]->Clear(opt);

    for( vpsiz_t k = 0; k < fProj.size(); ++k )
      fProj[k]->Clear(opt);
  }
#ifdef TESTCODE
  size_t nbytes = (char*)&t_coarse - (char*)&fNcombos + sizeof(t_coarse);
  memset( &fNcombos, 0, nbytes );
#endif

  // Clear tracking status
  fTrkStat = kTrackOK;

  f3dMatchBits.clear();
#ifdef MCDATA
  fMCHitBits.clear();
#endif
}

//_____________________________________________________________________________
Int_t Tracker::Decode( const THaEvData& evdata )
{
  // Decode all planes and fill hitpatterns per projection

#ifdef MCDATA
  const char* const here = "Tracker::Decode";

  struct MCHitCount {
    UInt_t min;
    UInt_t nfound;
  } mchitcount[kTypeEnd-kTypeBegin];

  Bool_t mcdata = TestBit(kMCdata);
  if( mcdata )
    memset( mchitcount, 0, (kTypeEnd-kTypeBegin)*sizeof(MCHitCount) );

  if( !fChecked ) {
    if( mcdata ) {
      fMCDecoder = dynamic_cast<const Podd::SimDecoder*>(&evdata);
      if( fMCDecoder == 0 ) {
	Error( Here(here), "MCdata flag set, but decoder is not a SimDecoder. "
	       "Fix database or replay configuration." );
	throw bad_config("Attempt to decode MCdata without a SimDecoder");
      }
    }
    fChecked = true;
  }
#endif
#ifdef TESTCODE
  // Save current event number for diagnostics
  fEvNum = evdata.GetEvNum();
#endif

#ifdef VERBOSE
  if( fDebug > 0 ) {
    cout << "============== Event ";
#ifdef TESTCODE
    cout << "# " << fEvNum << " ";
#endif
    cout << "==============" << endl;
  }
#endif

  // Decode the planes, then fill the hitpatterns in the projections
  //TODO: multithread?
  for( vpiter_t it = fProj.begin(); it != fProj.end(); ++it ) {
    Projection* theProj = *it;
     Int_t nhits = theProj->Decode( evdata );
    //  cout<<theProj->GetName()<<": "<<nhits<<endl;//getchar();
#ifdef MCDATA
    if( mcdata )
      mchitcount[theProj->GetType()].min = theProj->GetMinFitPlanes();
#endif
    // Sanity cut on overfull planes. nhits < 0 indicates overflow
    if( nhits < 0 ) {
      fTrkStat = kTooManyRawHits;
      continue;
    }
    // Fill the hitpattern if doing tracking
    if( TestBit(kDoCoarse) )
      theProj->FillHitpattern();
  }

#ifdef MCDATA
  const TSBSSimDecoder* simdata = 0;
  assert((&evdata) );
  assert( dynamic_cast<const TSBSSimDecoder*>(&evdata) );
  simdata = static_cast<const TSBSSimDecoder*>(&evdata);
  bool outputInfoForStdCode = false;
  std::vector<std::vector<Double_t>> mchits;// = simdata->GetAllMCHits();//get all mchits from signal files, excluding mchits from background files  
  if(outputInfoForStdCode){
    

    std::ofstream ofTrack;
    ofTrack.open("trackingInput.txt",ios::app);
  
    if( simdata->GetNMCTracks() > 0 ) {
      assert( dynamic_cast<TSBSSimTrack*>(simdata->GetMCTrack(0)) );
      TSBSSimTrack* trk = static_cast<TSBSSimTrack*>( simdata->GetMCTrack(0) );
      //cout<<trk->VX()<<"  "<<trk->VY()<<"  "<<trk->VZ()<<"    in mm"<<endl;
      //cout<<trk->fMomentum.X()<<"  "<<trk->fMomentum.Y()<<"  "<<trk->fMomentum.Z()<<"    in MeV"<<endl;

      int hitType = 2; //mc tracks
      double track_dx = trk->fMomentum.X()/trk->fMomentum.Z();
      double track_dy = trk->fMomentum.Y()/trk->fMomentum.Z();
      ofTrack<<hitType<<" "
	     <<trk->VX()<<" "
	     <<trk->VY()<<" "
	     <<trk->VZ()<<" "
	     <<track_dx<<" "
	     <<track_dy<<" "
	     <<trk->vertex_target.X()<<" "
	     <<trk->vertex_target.Y()<<" "
	     <<trk->vertex_target.Z()<<" "
	     <<trk->momentum_target.X()<<" "
	     <<trk->momentum_target.Y()<<" "
	     <<trk->momentum_target.Z()<<" "
	     <<sqrt( trk->fMomentum.X()*trk->fMomentum.X() + trk->fMomentum.Y()*trk->fMomentum.Y() + trk->fMomentum.Z()*trk->fMomentum.Z() )<<" "
	     <<endl;
    
    
      /*  cout<<" track projection on 5 planes in x: "
	  <<trk->VX() + track_dx * 50<<"  "
	  <<trk->VX() + track_dx * 200<<"  "
	  <<trk->VX() + track_dx * 350<<"  "
	  <<trk->VX() + track_dx * 500<<"  "
	  <<trk->VX() + track_dx * 1530<<"  "
	  <<endl;
	  cout<<" slope: "<<track_dx<<endl;
      */
    
      assert(trk);
    }
  
    for(auto& cluster:mchits){
      int hitType = 0;//type = 0 for MC hits, 1 for reconstructed
      int planeid = cluster[3];
      int moduleid= cluster[4];
      double charge = cluster[2];
      double time = cluster[5];
      double xPos = cluster[0];
      double yPos = cluster[1];
      ofTrack<<hitType<<" "
	     <<planeid<<" "
	     <<moduleid<<" "
	     <<charge<<" "
	     <<time<<" "
	     <<xPos<<" "
	     <<yPos<<" "
	     <<endl;
      //  cout<<planeid<<"  :  "<<xPos*1e3<<endl;getchar();
    }
    for(int planeid=0;planeid<5;planeid++)
      {
	Projection *projX = fProj[0];
	Projection *projY = fProj[1];
	Plane* plx = projX->GetPlane(planeid);
	Plane* ply = projY->GetPlane(planeid);
	TIterator* it;
	MCGEMHit* phit = 0;
	int hitType = 1;
	int projection;
	projection = 0;
	it = plx->GetHits()->MakeIterator();
	while( (phit = static_cast<MCGEMHit*>(it->Next())) ) {
	  ofTrack<<hitType<<" "
		 <<planeid<<" "
		 <<phit->GetModule()<<" "
		 <<projection<<" "
		 <<phit->GetADCmax()<<" "
		 <<phit->GetPos()<<" "
		 <<phit->Getprim_ratio()<<" "
		 <<phit->GetPeaktime()<<" "
		 <<endl;
	}
	projection = 1;
	it = ply->GetHits()->MakeIterator();
	while( (phit = static_cast<MCGEMHit*>(it->Next())) ) {
	  ofTrack<<hitType<<" "
		 <<planeid<<" "
		 <<phit->GetModule()<<" "
		 <<projection<<" "
		 <<phit->GetADCmax()<<" "
		 <<phit->GetPos()<<" "
		 <<phit->Getprim_ratio()<<" "
		 <<phit->GetPeaktime()<<" "
		 <<endl;
	}
      }
    int evtEndType = -1;
    ofTrack<<evtEndType<<endl;
    ofTrack.close();
  }
  bool checkflag = true, checkflag0 = true;
  std::map<int,int> mp;
  //if mchits gives more than 1 primary hit in any plane, drop this event

  for(auto& cluster:mchits){
    int planeid = cluster[3];
    if(mp.find(planeid)==mp.end()){
      mp[planeid] = 1;
    }else{
      checkflag = false;
      checkflag0 = false;
      break;
    }
  }
  
  //at most 1 hit exists in all plane
  if(0&&checkflag){
    std::ofstream outfile0;
    outfile0.open("cleanCluster_nhits.txt",ios::app);
    for(auto& cluster:mchits){
     
      int planeid = cluster[3];
      int moduleid= cluster[4];
      double& xpos = cluster[0];
      double& ypos = cluster[1];
      //double& charge = cluster[2];
      Projection *projX = fProj[0];
      Projection *projY = fProj[1];
      Plane* plx = projX->GetPlane(planeid);
      Plane* ply = projY->GetPlane(planeid);

      TIterator* it;
      MCGEMHit* phit = 0;
      MCGEMHit* phitToWrite = 0;
      Double_t max_prim_ratio = 0.0;
      //int nx=0,ny=0;
      int nbkgdPassTimeCut=0;

      // hitNo++;//
      it = plx->GetHits()->MakeIterator();
      while( (phit = static_cast<MCGEMHit*>(it->Next())) ) {
	Double_t tmp = phit->Getprim_ratio();
	if(phit->GetPeaktime()>30 && phit->GetPeaktime()<90)
	  nbkgdPassTimeCut++;
	if(tmp>max_prim_ratio){
	  max_prim_ratio = phit->Getprim_ratio();
	  phitToWrite = phit;
	}
      }
      
      if(max_prim_ratio>0) {
	//cout<<"plane: "<< planeid<<" moduleid: "<<moduleid<<" mcpos: "<<xpos<<"  pos: "<<phitToWrite->GetPos()<<" diff: "<<xpos - phitToWrite->GetPos()<<endl;getchar();
	outfile0<<planeid<<" "
		<<moduleid<<" "
		<<1e6*xpos<<" "
		<<1e6*phitToWrite->GetPos()<<" "
		<<phitToWrite->Getprim_ratio()<<" "
		<<phitToWrite->Getp_over_total()<<" "
		<<phitToWrite->Getb_over_total()<<" "
		<<phitToWrite->GetType()<<" "
		<<phitToWrite->GetSize()<<" "
		<<phitToWrite->GetPeaktime()<<" "
		<<phitToWrite->Getpos_sigma()<<" "
		<<phitToWrite->Gettime_sigma()<<" "
		<<phitToWrite->Getadc_sigma()<<" "
		<<phitToWrite->Getnb_bg()<<" "
		<<nbkgdPassTimeCut<<" "
		<<plx->GetNhits()<<" "
		<<endl;
      }else{
	outfile0<<planeid<<" "
		<<moduleid<<" "
		<<1e6*xpos<<" "
		<<endl;
      }

      // hitNo++;
      nbkgdPassTimeCut = 0;
      max_prim_ratio = 0.0;
      it = ply->GetHits()->MakeIterator();
      while( (phit = static_cast<MCGEMHit*>(it->Next())) ) {
	if(phit->GetPeaktime()>30 && phit->GetPeaktime()<90)
	  nbkgdPassTimeCut++;
	if(phit->Getprim_ratio()>max_prim_ratio){
	  max_prim_ratio = phit->Getprim_ratio();
	  phitToWrite = phit;
	}
      }
    
      if(max_prim_ratio>0) {
	outfile0<<planeid<<" "
		<<moduleid<<" "
		<<1e6*ypos<<" "
		<<1e6*phitToWrite->GetPos()<<" "
		<<phitToWrite->Getprim_ratio()<<" "
		<<phitToWrite->Getp_over_total()<<" "
		<<phitToWrite->Getb_over_total()<<" "
		<<phitToWrite->GetType()<<" "
		<<phitToWrite->GetSize()<<" "
		<<phitToWrite->GetPeaktime()<<" "
		<<phitToWrite->Getpos_sigma()<<" "
		<<phitToWrite->Gettime_sigma()<<" "
		<<phitToWrite->Getadc_sigma()<<" "
		<<phitToWrite->Getnb_bg()<<" "
		<<nbkgdPassTimeCut<<" "
		<<ply->GetNhits()<<" "
		<<endl;
      }else{
	outfile0<<planeid<<" "
		<<moduleid<<" "
		<<1e6*xpos<<" "
		<<endl;
      }

      delete it;
      //outfile0<<nx<<" "<<ny<<endl;
      //getchar();
      if(plx->GetNhits()!=1 || ply->GetNhits()!=1){
	checkflag0 = false;
      }

    }
    outfile0.close();
  }
  if(0&&checkflag0){// all planes have only 1 hits in both data and mc
    std::ofstream outfile;
    outfile.open("cleanCluster.txt",ios::app);
    for(auto& cluster:mchits){
      int planeid = cluster[3];
      double& xpos = cluster[0];
      double& ypos = cluster[1];
      double& charge = cluster[2];
     
      Projection *projX = fProj[0];
      Projection *projY = fProj[1];
      Plane* plx = projX->GetPlane(planeid);
      Plane* ply = projY->GetPlane(planeid);
      //    if(plx->GetNhits()!=1 || ply->GetNhits()!=1){
	
      //   }
      TIterator* it;
      GEMHit* phit = 0;
      //cout<<plx->GetNhits()<<" <-x:y-> "<<ply->GetNhits()<<endl;
      it = plx->GetHits()->MakeIterator();
      while( (phit = static_cast<GEMHit*>(it->Next())) ) {
	outfile<<planeid<<" "<<1e6*(phit->GetPos()-xpos)<<" "<<charge/phit->GetADCsum()
	       <<" "<<phit->GetPeaktime();
      }
      it = ply->GetHits()->MakeIterator();
      while( (phit = static_cast<GEMHit*>(it->Next())) ) {
	outfile<<" "<<1e6*(phit->GetPos()-ypos)<<" "<<charge/phit->GetADCsum()
	       <<" "<<phit->GetPeaktime()<<endl;
      }
      delete it;
    }
    outfile.close();
    
  }
  


#endif

#ifdef MCDATA
  // For MCdata, check which MC track points were detected
  if( mcdata ) {
    if( fMCDecoder->GetNMCTracks() > 0 ) {
      // TODO: support multiple tracks
      assert( dynamic_cast<MCTrack*>(fMCDecoder->GetMCTrack(0)) );
      MCTrack* trk = static_cast<MCTrack*>( fMCDecoder->GetMCTrack(0) );
      assert(trk);
      trk->fNHitsFound = 0;
      for( Int_t i = 0; i < fMCDecoder->GetNMCPoints(); ++i ) {
	MCTrackPoint* pt = fMCDecoder->GetMCPoint(i);
	//assert( pt->fMCTrack == trk->fNumber ); // temporary
	if( pt->fMCTrack != trk->fNumber )
	  continue;
	FindHitForMCPoint( pt, fMCPointUpdater );
	if( TESTBIT(pt->fStatus, MCTrackPoint::kCorrectFound) ) {
	  trk->fNHitsFound++;
	  mchitcount[pt->fType].nfound++;
	}
      }
      // Is number of found hits sufficient for the 3D fit?
      if( trk->fNHitsFound-4 >= fMinNdof ) {
	// If so, check if there are enough hits for each 2D fit, too
	bool good = !fProj.empty();
	for( vpiter_t it = fProj.begin(); it != fProj.end(); ++it ) {
	  const MCHitCount& mcnt = mchitcount[(*it)->GetType()];
	  assert( mcnt.min > 0 );
	  if( mcnt.nfound < mcnt.min ) {
	    good = false;
	    break;
	  }
	}
	if( good )
	  SETBIT( trk->fReconFlags, kHitsFound );
	else
	  CLRBIT( trk->fReconFlags, kHitsFound );
	
	// Fit the MC points. Results are stored in fMCFitPar of the MC track
	FitMCPoints( trk );
      }
    }
    else
      Warning( Here(here), "No MC track in event %d?", evdata.GetEvNum() );
  }
#endif

  return 0;
}

#ifdef MCDATA
//_____________________________________________________________________________
void Tracker::MCPointUpdater::UpdateHit( MCTrackPoint* pt, Hit* hit,
					 Double_t x ) const
{
  // Standard updater of MCTrackPoint hit residual. x is the track point
  // position converted to the coordinate system of the plane that the hit
  // belongs to.

  assert( pt );
  if( hit )
    pt->fHitResid = hit->GetPos() - x;
}

//_____________________________________________________________________________
void Tracker::MCPointUpdater::UpdateTrack( MCTrackPoint* pt, Plane* pl,
					   THaTrack* trk, Double_t x ) const
{
  // Standard updater of MCTrackPoint tracks residual. x is the track point
  // position converted to the coordinate system of the plane that the hit
  // belongs to.

  assert( pt );
  if( trk ) {
    Double_t cosa   = pl->GetProjection()->GetCosAngle();
    Double_t sina   = pl->GetProjection()->GetSinAngle();
    Double_t slope  = trk->GetDTheta()*cosa + trk->GetDPhi()*sina;
    Double_t trkpos = trk->GetDX()*cosa + trk->GetDY()*sina + slope*pl->GetZ();
    pt->fTrackResid = trkpos - x;
    SETBIT(pt->fStatus, MCTrackPoint::kUsedInTrack);
  }
  else
    CLRBIT(pt->fStatus, MCTrackPoint::kUsedInTrack);
}

//_____________________________________________________________________________
Hit* Tracker::FindHitForMCPoint( MCTrackPoint* pt,
				 MCPointUpdater* Updater ) const
{
  // Find the best-matching hit for the given MC track point

  typedef Rpvec_t::const_iterator citer_t;
  assert( pt );
  assert( pt->fMCTrack > 0 );

  // Convert the MC point to the tracker frame
  TVector3 point = pt->fMCPoint - fOrigin;
  // Find the plane for this point. Since the plane numbering can be ambiguous,
  // we search by z-position. Certain trackers have more than one plane with the
  // same z, which are then distinguished by the plane type (u,v,x etc.)
  pair<citer_t,citer_t> plane_range =
    equal_range( ALL(fPlanes), point.Z(), Plane::ZIsNear(0.001) );
  
  Hit* hit = 0;
  Double_t x = kBig;
  for( citer_t it = plane_range.first; it != plane_range.second; ++it ) {
    const Plane* pl = *it;
    if( pl->GetType() != pt->fType )
      continue;
    
    // Calculate position of the point in its projection coordinate system
    Double_t cosa = pl->GetProjection()->GetCosAngle();
    Double_t sina = pl->GetProjection()->GetSinAngle();
    x = point.X()*cosa + point.Y()*sina;

    pair<Int_t,Int_t> hit_range = pl->FindHitsInRange( x - pt->fgWindowSize,
						       x + pt->fgWindowSize );
    Int_t nfound = hit_range.second - hit_range.first;
    pt->fNFound = nfound;
    if( nfound > 0 )
      SETBIT(pt->fStatus, MCTrackPoint::kFound);

    // Find two types of reconstructed hits closest to the point's position:
    // 1) any hit
    // 2) a hit originating at least in part from the point's MC track
    Double_t anyHitPos, mcHitPos;
    Hit* anyHit = pl->FindNearestHitAndPos( x, anyHitPos );
    Hit* mcHit  = pl->FindNearestMCHitAndPos( x, pt->fMCTrack, mcHitPos );
    // cout<<x<<" "<<anyHitPos<<" "<<mcHitPos<<endl;
    // Sanity checks
    assert( not mcHit or anyHit ); // mcHit implies anyHit
    assert( nfound == 0 or  // nfound > 0 implies anyHit and pos within window
	    (anyHit != 0 and TMath::Abs(x-anyHitPos) <= pt->fgWindowSize) );
    assert( mcHit == 0 or   // mcHit, if any, cannot be closer than anyHit
	    TMath::Abs(x-anyHitPos) <= TMath::Abs(x-mcHitPos) );

    // Flag too-distant MC hits. In this case, nfound may or may not be zero
    if( mcHit and TMath::Abs(x-mcHitPos) > pt->fgWindowSize )
      SETBIT(pt->fStatus, MCTrackPoint::kCorrectFarAway);

    // Record what we found in the status flags of the point
    if( mcHit ) {
      SETBIT(pt->fStatus, MCTrackPoint::kCorrectFound);
      if( anyHit == mcHit )
	SETBIT(pt->fStatus, MCTrackPoint::kCorrectClosest);
      MCHitInfo* mcinfo = dynamic_cast<Podd::MCHitInfo*>(mcHit);
      assert( mcinfo );
      if( mcinfo->fContam > 0 )
	SETBIT(pt->fStatus, MCTrackPoint::kContaminated);
    }
    if( anyHit ) {
      MCHitInfo* mcinfo = dynamic_cast<Podd::MCHitInfo*>(anyHit);
      assert( mcinfo );
      Int_t mctrack = mcinfo->fMCTrack;
      if( mctrack > 0 and mctrack != pt->fMCTrack )
	SETBIT(pt->fStatus, MCTrackPoint::kWrongTrack);
    }
    hit = anyHit;
    break;  // plane of correct type found
  }
  if( Updater )
    Updater->UpdateHit( pt, hit, x );
  return hit;
}

//_____________________________________________________________________________
THaTrack* Tracker::FindTrackForMCPoint( MCTrackPoint* /* pt */,
					TClonesArray& /* tracks */,
					MCPointUpdater* /* Updater */ ) const
{
  // Find the best-matching track for the given MC track point

  //TODO
  return 0;
}

//_____________________________________________________________________________
Int_t Tracker::FitMCPoints( MCTrack* mctrk ) const
{
  // Fit the MC points to the expected track shape. This version performs a
  // straight-line fit. Results are stored in fMCFitPar of the MC track.
  //
  // These results can be used to get the baseline for the achievable
  // resolutions, free of all reconstruction effects, but purely a
  // measure of the effect of materials and possibly, the fitting algorithm.

  assert( mctrk );
  assert( not fProj.empty() );

  // Lookup table for projections by type
  const Projection* projections[ kTypeEnd-kTypeBegin ];
  memset( projections, 0, (kTypeEnd-kTypeBegin)*sizeof(Projection*) );
  for( vpsiz_t ip = 0; ip < fProj.size(); ++ip ) {
    const Projection* proj = fProj[ip];
    projections[ proj->GetType() ] = proj;
  }

  // Mimic the fit procedure used in FitTrack, except for the weighting
  TMatrixDSym AtA(4);
  TVectorD Aty(4);
  
  // Fill fit matrixes
  Int_t npoints = 0;
  for( Int_t i = 0; i < fMCDecoder->GetNMCPoints(); ++i ) {
    MCTrackPoint* pt = fMCDecoder->GetMCPoint(i);
    if( pt->fMCTrack != mctrk->fNumber )
      continue;
    const Projection* proj = projections[ pt->fType ];
    assert( proj );
    TVector3 hitpos = pt->fMCPoint - fOrigin;
    if( IsRotated() )
      hitpos *= fInvRot;
    Double_t cosa = proj->GetCosAngle();
    Double_t sina = proj->GetSinAngle();
    Double_t x = hitpos.X()*cosa + hitpos.Y()*sina;
    Double_t z = hitpos.Z();
    Double_t Ai[4] = { cosa, cosa*z, sina, sina*z };
    for( int j = 0; j<4; ++j ) {
      for( int k = j; k<4; ++k ) {
	AtA(j,k) += Ai[j] * Ai[k];
      }
      Aty(j) += Ai[j] * x;
    }
    ++npoints;
  }
  assert( npoints > 4 );  // assured by caller

  for( int j = 0; j<4; ++j )
    for( int k = j+1; k<4; ++k )
      AtA(k,j) = AtA(j,k);

  // Solve the normal equations
  TDecompChol chol(AtA);
  Bool_t ok = chol.Solve(Aty);
  if( !ok ) return -1;

  // Results are here, in order x, x', y, y'
  assert( Aty.GetNrows() == 4 );
  Double_t* coef = Aty.GetMatrixArray();

  // Calculate chi2
  //FIXME: too much code duplication here somehow
  Double_t chi2 = 0;
  for( Int_t i = 0; i < fMCDecoder->GetNMCPoints(); ++i ) {
    MCTrackPoint* pt = fMCDecoder->GetMCPoint(i);
    if( pt->fMCTrack != mctrk->fNumber )
      continue;
    const Projection* proj = projections[ pt->fType ];
    TVector3 hitpos = pt->fMCPoint - fOrigin;
    if( IsRotated() )
      hitpos *= fInvRot;
    Double_t cosa = proj->GetCosAngle();
    Double_t sina = proj->GetSinAngle();
    Double_t x = hitpos.X()*cosa + hitpos.Y()*sina;
    Double_t z = hitpos.Z();
    Double_t slope = coef[1]*cosa + coef[3]*sina;
    Double_t pos   = coef[0]*cosa + coef[2]*sina + slope*z;
    Double_t diff  = pos - x;
    chi2 += diff*diff;
  }

  // Convert to lab frame
  TVector3 pos( coef[0], coef[2], 0.0 );
  TVector3 dir( coef[1], coef[3], 1.0 );
  if( fIsRotated ) {
    dir *= fRotation;
    dir *= 1.0/dir.Z();
    pos *= fRotation;
  }
  pos += fOrigin;
  if( TestBit(kProjTrackToZ0) ) {
    pos -= pos.Z() * dir;
  }
  coef[0] = pos.X();
  coef[1] = dir.X();
  coef[2] = pos.Y();
  coef[3] = dir.Y();

  // Save results with the MC track
  assert( sizeof(mctrk->fMCFitPar)/sizeof(mctrk->fMCFitPar[0]) >= 6 );
  Double_t* res = mctrk->fMCFitPar;
  memcpy( res, coef, 4*sizeof(Double_t) );
  res[4] = chi2;
  res[5] = npoints;

  return npoints;
}

#endif // MCDATA

//_____________________________________________________________________________
Int_t Tracker::End( THaRunBase* run )
{
#ifdef TESTCODE
  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    fPlanes[iplane]->End(run);
#endif
  return 0;
}

//_____________________________________________________________________________
class AddFitCoord
{
public:
  AddFitCoord() {}
  void operator() ( Road* rd, Road::Point* p, const vector<Double_t>& coef )
  {
    const Projection* proj = rd->GetProjection();
    Double_t slope2d  = rd->GetSlope();
    Double_t trkpos2d = rd->GetPos() + p->z * slope2d;
    Double_t cosa     = proj->GetCosAngle();
    Double_t sina     = proj->GetSinAngle();
    Double_t slope    = coef[1]*cosa + coef[3]*sina;
    Double_t pos      = coef[0]*cosa + coef[2]*sina + slope * p->z;
    p->hit->GetPlane()->AddFitCoord
      ( FitCoord( p->hit, rd, p->x, trkpos2d, slope2d, pos, slope ));
  }
};

//_____________________________________________________________________________
#ifdef VERBOSE
class PrintFitPoint
{
public:
  PrintFitPoint() {}
  void operator() ( Road*, Road::Point* p, const vector<Double_t>& )
  {
    cout << p->hit->GetPlane()->GetName() << " "
	 << "z = " << p->z << " x = " << p->x << endl;
  }
};
#endif

//_____________________________________________________________________________
class FillFitMatrix
{
public:
  FillFitMatrix( TMatrixDSym& AtA, TVectorD& Aty )
    : fAtA(AtA), fAty(Aty), fNpoints(0)
  {
    assert( fAtA.GetNrows() == 4 && fAty.GetNrows() == 4 );
  }

  void operator() ( Road* rd, Road::Point* p, const vector<Double_t>& )
  {
    const Projection* proj = rd->GetProjection();
    Double_t cosa = proj->GetCosAngle();
    Double_t sina = proj->GetSinAngle();
    Double_t Ai[4] = { cosa, cosa * p->z, sina, sina * p->z };
    Double_t s2 = 1.0/(p->res()*p->res());
    for( int j = 0; j<4; ++j ) {
      for( int k = j; k<4; ++k ) {
	fAtA(j,k) += Ai[j] * s2 * Ai[k];
      }
      fAty(j) += Ai[j] * s2 * p->x;
    }
    ++fNpoints;
  }
  Int_t GetNpoints() const { return fNpoints; }
private:
  TMatrixDSym& fAtA;
  TVectorD&    fAty;
  Int_t        fNpoints;
};

//_____________________________________________________________________________
class CalcChisquare
{
public:
  CalcChisquare() : fChi2(0) {}
  void operator() ( Road* rd, Road::Point* p, const vector<Double_t>& coef )
  {
    const Projection* proj = rd->GetProjection();
    Double_t cosa  = proj->GetCosAngle();
    Double_t sina  = proj->GetSinAngle();
    Double_t slope = coef[1]*cosa + coef[3]*sina;
    Double_t pos   = coef[0]*cosa + coef[2]*sina + slope*p->z;
    Double_t diff  = (pos - p->x) / p->res();
    fChi2 += diff*diff;
  }
  void     Clear() { fChi2 = 0; }
  Double_t GetResult() const { return fChi2; }
private:
  Double_t fChi2;
};

//_____________________________________________________________________________
// Set bit for the unique number (fDefinedNum) of the plane of the given point
class SetPlaneBit
{
public:
  SetPlaneBit() : fBits(0) {}
  void operator() ( Road*, Road::Point* p, const vector<Double_t>& )
  {
    fBits |= 1U << p->hit->GetPlane()->GetDefinedNum();
  }
  void     Clear() { fBits = 0; }
  UInt_t   GetBits() const { return fBits; }
private:
  UInt_t fBits;
};

#ifdef MCDATA
//_____________________________________________________________________________
// Set bit for the unique number (fDefinedNum) of the plane of the given point
// if the point corresponds to a true MC hit (contaminated or not)
class SetMCPlaneBit
{
public:
  SetMCPlaneBit() : fBits(0) {}
  void operator() ( Road*, Road::Point* p, const vector<Double_t>& )
  {
    Hit* hit = p->hit;
    MCHitInfo* mcinfo = dynamic_cast<Podd::MCHitInfo*>(hit);
    assert(mcinfo); // may only call this with MC-generated data
    if( mcinfo->fMCTrack > 0 )
      fBits |= 1U << hit->GetPlane()->GetDefinedNum();
  }
  void     Clear() { fBits = 0; }
  UInt_t   GetBits() const { return fBits; }
private:
  UInt_t fBits;
};
#endif

//_____________________________________________________________________________
template< typename Action >
Action Tracker::ForAllTrackPoints( const Rvec_t& roads,
				   const vector<Double_t>& coef, Action action )
{
  // Apply action to all points from roads and the track given by coef

  for( Rvec_t::const_iterator it = roads.begin(); it != roads.end(); ++it ) {
    Road* rd = *it;
    const Road::Pvec_t& points = rd->GetPoints();
    for( Road::Pvec_t::const_iterator it2 = points.begin();
	 it2 != points.end(); ++it2 ) {
      Road::Point* p = *it2;
      action( rd, p, coef );
    }
  }
  return action;
}

//_____________________________________________________________________________
void Tracker::FitErrPrint( Int_t err ) const
{
  static const char* const here = "FitTrack";

#ifdef TESTCODE
  Error( Here(here), "Event %d: Failure fitting 3D track, err = %d. "
	 "Should never happen. Call expert.", fEvNum, err );
#else
  Error( Here(here), "Failure fitting 3D track, err = %d. "
	 "Should never happen. Call expert.", err );
#endif
}

//_____________________________________________________________________________
Int_t Tracker::FitTrack( const Rvec_t& roads, vector<Double_t>& coef,
			 Double_t& chi2, TMatrixDSym* coef_covar ) const
{
  // Linear minimization routine to fit a physical straight line track through
  // the hits registered in the different projection planes.
  //
  // This is a much streamlined version of ROOT's TLinearFitter that solves
  // the normal equations with weights, (At W A) b = (At W) y, where AWb = Wy,
  // using Cholesky decomposition, TDecompChol. The model used is
  //   y_i = P_i * T_i
  //       = ( x + z_i * mx, y + z_i * my) * ( cos(a_i), sin(a_i) )
  // where
  //   y_i is the measured coordinate in the i-th plane at z_i, i=0...nplanes-1
  //   P_i is the physical track intersection point with the z_i plane
  //   T_i is the axis unit vector of the i-th plane
  //   a_i is the angle of the coordinate axis of the i-th plane w.r.t. x
  //   b = (x,mx,y,my) are the track parameters to be fitted, origin x,y and
  //       slopes mx,my.
  // Hence the i-th row of A is
  // A_i,0..3 = ( cos(a_i), z_i*cos(a_i), sin(a_i), z_i*sin(a_i) )
  //
  // "roads" contains a set of Roads that successfully combine in 3-d, one
  // Road* per projection. Each road, in turn, contains a set of (y_i,z_i)
  // coordinates, at most one per plane of the projection type, that
  // give the best 2D track fit within the road.
  //
  // "coef" (=b) are the fitted track parameters, x, x'(=mx), y, y'(=my).
  // The reference system is the z=0 plane (usually the first chamber).
  //
  // "chi2" is the chi2 of the fit (not normalized).
  //
  // "coef_covar" is the covariance matrix of the fit results. Since it is
  // expensive to calculate, it is only filled if the argument is supplied.
  //
  // The return value is the number of degrees of freedom of the fit, i.e.
  // npoints-4 > 0, or negative if too few points or matrix inversion error

  // Fill the (At W A) matrix and (At W y) vector with the measured points
  TMatrixDSym AtA(4);
  TVectorD Aty(4);
  FillFitMatrix f = ForAllTrackPoints(roads, coef, FillFitMatrix(AtA, Aty));

  Int_t npoints = f.GetNpoints();
  assert( npoints > 4 );
  if( npoints <=4 ) return -1; // Meaningful fit not possible

  // FillFitMatrix only fills the upper triangle for efficiency
  //TODO: make a function of FillFitMatrix?
  for( int j = 0; j<4; ++j )
    for( int k = j+1; k<4; ++k )
      AtA(k,j) = AtA(j,k);

  // Invert the characteristic matrix and solve the normal equations.
  // As in ROOT's TLinearFitter, we use a Cholesky decomposition
  // to do this (since AtA is symmetric and positive-definite).
  // For more speed but less accuracy, one could use TMatrixDSymCramerInv.
  TDecompChol chol(AtA);
  Bool_t ok = chol.Solve(Aty);
  assert(ok);
  if( !ok ) return -2; //Urgh, decomposition failed. Should never happen

  // Copy results to output vector in order x, x', y, y'
  assert( Aty.GetNrows() == 4 );
  coef.assign( Aty.GetMatrixArray(), Aty.GetMatrixArray()+Aty.GetNrows() );

#ifdef VERBOSE
  if( fDebug > 2 ) {
    cout << "Points in 3D fit:" << endl;
    ForAllTrackPoints( roads, coef, PrintFitPoint() );
  }
#endif

  // Calculate chi2
  chi2 = ForAllTrackPoints(roads, coef, CalcChisquare()).GetResult();

  // Calculate covariance matrix of the parameters, if requested
  if( coef_covar ) {
    if( coef_covar->GetNrows() != 4 )
      coef_covar->ResizeTo(4,4);
    Bool_t ok = chol.Invert(*coef_covar);
    assert(ok);
    if( !ok ) return -3; // Urgh, inversion failed. Should never happen
  }

#ifdef VERBOSE
  if( fDebug > 1 )
    cout << "3D fit:  x/y = " << coef[0] << "/" << coef[2] << " "
	 << "mx/my = " << coef[1] << "/" << coef[3] << " "
	 << "ndof = " << npoints-4 << " rchi2 = " << chi2/(double)(npoints-4)
	 << endl;
#endif

  return npoints-4;
}

//_____________________________________________________________________________
Int_t Tracker::NewTrackCalc( Int_t , THaTrack*, const TVector3&,
			     const TVector3&, const FitRes_t& )
{
  // Hook for additional calculations for a newly generated track

  return 0;
}

//_____________________________________________________________________________
THaTrack* Tracker::NewTrack( TClonesArray& tracks, const FitRes_t& fit_par )
{
  // Make new track with given parameters. Correct for coordinate offset
  // and rotation of this detector. Used by CoarseTrack to generate all tracks.

  // The following coordinates ("detector" coordinates, e.g. "tr.d_x") are
  // in the Tracker frame. fOrigin of the planes is already being taken into
  // account when fitting tracks since the hit positions are corrected for it
  // (see Plane::GetStart() and Plane::GetZ())

  Double_t d_x  = fit_par.coef[0];
  Double_t d_xp = fit_par.coef[1];
  Double_t d_y  = fit_par.coef[2];
  Double_t d_yp = fit_par.coef[3];
 

  // Correct the track coordinates for the position offset and rotation of the
  // Tracker system, resulting in "global" reference coordinates, e.g. "tr.x".
  // For a focusing spectrometer, the "global" system usually is the "focal
  // plane" system, e.g. the detector stack frame.
  // fOrigin is specified in the global frame. It is never rotated.
  //
  // NB: The track origin is projected to z=0 of the "global" system, i.e. the
  // one enclosing the tracker
  //
  // Target/vertex coordinate reconstruction is done by the parent spectrometer
  // apparatus.

  assert( fIsRotated == !fRotation.IsIdentity() );
  Double_t xp, yp, x, y;
  TVector3 pos( d_x, d_y, 0.0 );
  TVector3 dir( d_xp, d_yp, 1.0 );
  
  if( fIsRotated ) {
    // Rotate the TRANSPORT direction vector to global frame
    dir *= fRotation;
    // TODO: use different type of track without this limitation
    if( TMath::Abs(dir.Theta()-TMath::PiOver2()) < 1e-2 ) {
      Error( Here("NewTrack"), "Limitation: Track (nearly) perpendicular "
	     "to z-axis not supported by TRANSPORT formalism in THaTrack "
	     "class. Skipping this track." );
      pos.Print();
      dir.Print();
      return 0;
    }
    // Normalize to z=1 to get new TRANSPORT direction vector
    dir *= 1.0/dir.Z();
    assert( TMath::Abs(dir.Z()-1.0) < 1e-6 );
    // Rotate track origin to global frame
    pos *= fRotation;
  }
  xp = dir.X();
  yp = dir.Y();
  pos += fOrigin;
  if( TestBit(kProjTrackToZ0) ) {
    // If configured, project along the track direction to z=0 in global frame
    pos -= pos.Z() * dir;
    assert( TMath::Abs(pos.Z()) < 1e-6 );
  }
  x = pos.X();
  y = pos.Y();
 
  THaTrack* newTrack = AddTrack( tracks, x, y, xp, yp );
  assert( newTrack );
  newTrack->SetD( d_x, d_y, d_xp, d_yp );
  newTrack->SetChi2( fit_par.chi2, fit_par.ndof );
  // Index of currect track
  Int_t idx = tracks.GetLast();
  assert( idx >= 0 );
  //TODO: make a TrackID?

  // Save the bitpattern of planes whose hits were used in the track fit
  UInt_t hitbits = ForAllTrackPoints( *fit_par.roads, fit_par.coef,
				      SetPlaneBit() ).GetBits();
  newTrack->SetFlag( hitbits );
  
  // UInt_t matchbits = ForAllTrackPoints( *fit_par.roads, fit_par.coef,
  // 					  SetPlaneBit() ).GetBits();
  // f3dMatchBits.resize(idx+1);
  // f3dMatchBits[idx] = matchbits;
  
#ifdef MCDATA
  // For MC data, determine which of the hits that were used for fitting this
  // track are true MC hits, i.e. hits left by the MC track we are trying to
  // reconstruct.
  if( TestBit(kMCdata) ) {
    hitbits = ForAllTrackPoints( *fit_par.roads, fit_par.coef,
				 SetMCPlaneBit() ).GetBits();
    // Use resize instead of push_back since there could be gaps in the
    // sequence of track indexes if other trackers also fill the tracks array
    assert( static_cast<vrsiz_t>(idx) >= fMCHitBits.size() );
    fMCHitBits.resize(idx+1);
    fMCHitBits[idx] = hitbits;
  }
#endif

  // For any planes in calibration mode, save the hits closest to this track.
  // Minimizing the hit residuals is the standard procedure for alignment
  // calibration.
  ForAllTrackPoints( *fit_par.roads, fit_par.coef, AddFitCoord() );
  for( Rpvec_t::const_iterator it = fCalibPlanes.begin(); it !=
	 fCalibPlanes.end(); ++it ) {
    Plane* pl = *it;
    pl->RecordNearestHits( newTrack, *fit_par.roads );
  }

  // Do additional calculations for the new track (used by derived classes)
  NewTrackCalc( idx, newTrack, pos, dir, fit_par );

  return newTrack;
}

//_____________________________________________________________________________
class CheckTypes : public unary_function<Road*,void>
{
  // Function object for testing the projection occupancy of roads.
  // Applied to a road, it saves the road's plane type (u,v,x..).
  // Testing the object returns true if all types specified in the 'req'
  // argument to the constructor have been seen in all the roads tested.
public:
  CheckTypes( UInt_t req ) : fReq(req), fActive(0) {}
  void operator() ( const Road* rd )
  { fActive |= 1U << rd->GetProjection()->GetType(); }
  operator bool()       const { return (fActive & fReq) == fReq; }
  bool     operator!()  const { return (fActive & fReq) != fReq; }
  UInt_t   GetTypes()   const { return fActive; }
private:
  UInt_t fReq;
  UInt_t fActive;
};

//_____________________________________________________________________________
class AnySharedHits : public unary_function<Road*,bool>
{
public:
  AnySharedHits() {}
  bool operator() ( const Road* rd )
  {
    const Hset_t& test_hits = rd->GetHits();
//     cout << "--- testing: " << endl;
//     PrintHits(test_hits);
    for( Hset_t::const_iterator it = test_hits.begin(); it != test_hits.end();
	 ++it ) {
      if( fHits.find(*it) != fHits.end() )
	return true;
    }
    return false;
  }
  const AnySharedHits& use( const Rset_t& tuple )
  {
    fHits.clear();
    for( Rset_t::const_iterator it = tuple.begin(); it != tuple.end(); ++it ) {
      const Road* rd = *it;
      fHits.insert( ALL(rd->GetHits()) );
    }
//     PrintHits(fHits);
    return *this;
  }
private:
  Hset_t fHits;
};

//_____________________________________________________________________________
class Tracker::TrackFitWeight
{
  // Object for sorting fits. The smaller a track's TrackFitWeight, the
  // "better" a track is considered to be.
public:
  TrackFitWeight( const Tracker::FitRes_t& fit_par ) :
    fNdof(fit_par.ndof), fChi2(fit_par.chi2) { assert(fNdof); }

  bool operator<( const TrackFitWeight& rhs ) const
  {
    // The "best" tracks have the largest number of hits and the smallest chi2
    //TODO: devalue tracks with very large chi2?
    return (fNdof != rhs.fNdof) ?
      (fNdof > rhs.fNdof) : (fChi2 < rhs.fChi2);
  }
  operator double() const { return fChi2/(double)fNdof; }
private:
  UInt_t   fNdof;
  Double_t fChi2;
};

//_____________________________________________________________________________
template< typename Container, typename Weight, typename TestF,
	  typename QuitF > Double_t
OptimalN( const Container& choices, const multimap<Weight,Container>& weights,
	  vector<Container>& picks, TestF testf, QuitF quitf )
{
  // This is the second-level de-ghosting algorithm, operating on fitted
  // 3D tracks.  It selects the best set of tracks if multiple roads
  // combinations are present in one or more projection.
  //
  // Requirements:
  //  "Container": STL container (set or sorted vector) supporting begin(),
  //               insert(), end(), swap(), and containing "Element"s sorted
  //               by Element's operator<.  Here, Element = Road*
  //  "Weight":    sorting functor, supporting operator< and operator double()
  //  "TestF":     Functor for eliminating additional leftover elements,
  //               must support:
  //                  operator() (const Element&): test this element,
  //                                         eliminate if true
  //                  use(const Container&): test against elements in this
  //                                         container
  //               (use() is in lieu of a constructor because the object
  //                is constructed by the caller (possibly with parameters)).
  //  "QuitF":     Functor for terminating search early. Applied to all
  //               leftover elements. Must support
  //                  operator() (const Element&)
  //                  operator bool(): End search if false after running on
  //                                   all elements

  //TODO: how efficient? O(N^2)? (N iterations on include)

  picks.clear();
  Container choices_left(choices);

  //TODO: try minimizing chi2 sum
  double wsum = 0;
  typename multimap<Weight,Container>::const_iterator it = weights.begin();
  for( ; it != weights.end(); ++it ) {
    const Container& tuple = (*it).second;
    if( includes(ALL(choices_left), ALL(tuple)) ) {
      // Pick next set of still-available elements in order of ascending weight
      picks.push_back( tuple );
      // Sum of weights of chosen sets, converted to double
      wsum += (*it).first;
      // Remove chosen elements from the remaining choices
      Container choices_less_tuple;
      set_difference( ALL(choices_left), ALL(tuple),
		      inserter(choices_less_tuple, choices_less_tuple.end()) );
      // Test each remaining element against the chosen set of elements
      // and remove it if it compares true (=remove elems with common traits)
      Container new_choices_left;
      remove_copy_if( ALL(choices_less_tuple),
		      inserter(new_choices_left, new_choices_left.end()),
		      testf.use(tuple) );
      // Quit if the quit function, applied to each element still remaining
      // now, is no longer true
      if( !for_each(ALL(new_choices_left), quitf) )
	break;
      // Update the set of remaining elements
      choices_left.swap( new_choices_left );
    }
  }

  return wsum;
}

//_____________________________________________________________________________
void Tracker::Add3dMatch( const Rvec_t& selected, Double_t matchval,
			  list< pair<Double_t,Rvec_t> >& combos_found,
			  Rset_t& unique_found ) const
{
  // Save road combination with good matchvalue

  combos_found.push_back( make_pair(matchval,selected) );

  // Since not all of the roads in 'roads' may make a match,
  // keep track of unique roads found for each projection type

  unique_found.insert( ALL(selected) );

#ifdef VERBOSE
  if( fDebug > 3 )
    cout << "ACCEPTED" << endl;
#endif
};

//_____________________________________________________________________________
UInt_t Tracker::MatchRoadsGeneric( vector<Rvec_t>& roads, const UInt_t ncombos,
				   list< pair<Double_t,Rvec_t> >& combos_found,
				   Rset_t& unique_found )
{
  // General MatchRoad algorithm:
  //  - find all front and back intersections [ n(n-1)/2 each ]
  //  - compute weighted center of gravity of intersection points
  //  - sum dist^2 of points to center of gravity -> matchval
  // Requires at least 3 projections

  vector<Rvec_t>::size_type nproj = roads.size();
  assert( nproj >= 3 );

#ifdef VERBOSE
  if( fDebug > 0 )
    cout << "generic algo):";
#endif

  // Vector holding a combination of roads to test. One road from each
  // projection
  Rvec_t selected;

  UInt_t nfound = 0;
  Double_t zback = fPlanes.back()->GetZ();

  vector<TVector2> fxpts, bxpts;
  fxpts.reserve( nproj*(nproj-1)/2 );
  bxpts.reserve( nproj*(nproj-1)/2 );

  for( UInt_t i = 0; i < ncombos; ++i ) {
    Double_t matchval = 0.0;

    NthCombination( i, roads, selected );
    assert( selected.size() == nproj );

    fxpts.clear();
    bxpts.clear();
    TVector2 fctr, bctr;
    for( Rvec_t::iterator it1 = selected.begin(); it1 != selected.end();
	 ++it1 ) {
      for( Rvec_t::iterator it2 = it1+1; it2 != selected.end(); ++it2 ) {
	//TODO: weigh with uncertainties of coordinates?
	fxpts.push_back( (*it1)->Intersect(*it2, 0.0) );
	bxpts.push_back( (*it1)->Intersect(*it2, zback) );
	fctr += fxpts.back();
	bctr += bxpts.back();
#ifdef VERBOSE
	if( fDebug > 3 ) {
	  cout << (*it1)->GetProjection()->GetName()
	       << (*it2)->GetProjection()->GetName()
	       << " front(" << fxpts.size() << ") = ";
	  fxpts.back().Print();
	  cout << (*it1)->GetProjection()->GetName()
	       << (*it2)->GetProjection()->GetName()
	       << " back (" << bxpts.size() << ") = ";
	  bxpts.back().Print();
	}
#endif
      }
    }
    assert( fxpts.size() <= nproj*(nproj-1)/2 );
    assert( bxpts.size() == fxpts.size() );
    fctr /= static_cast<Double_t>( fxpts.size() );
    if( !fPlanes.front()->Contains(fctr) )
      continue;
    bctr /= static_cast<Double_t>( fxpts.size() );
    if( !fPlanes.back()->Contains(bctr) )
      continue;
    for( vector<TVector2>::size_type k = 0; k < fxpts.size(); ++k ) {
      matchval += (fxpts[k]-fctr).Mod2() + (bxpts[k]-bctr).Mod2();
    }
#ifdef VERBOSE
    if( fDebug > 3 ) {
      cout << "fctr = "; fctr.Print();
      cout << "bctr = "; bctr.Print();
      cout << "matchval = " << matchval << endl;
    }
#endif
    // We could just connect fctr and bctr here to get an approximate
    // 3D track. But the linear minimization in FitTrack is the right
    // way to do this.

    if( matchval < f3dMatchCut ) {
      ++nfound;
      Add3dMatch( selected, matchval, combos_found, unique_found );
    }
  } //for(ncombos)

  return nfound;
}

//_____________________________________________________________________________
// Functor for sorting a vector< vector<Road*> > according to the order
// of projections types defined by the supplied lookup table.
// Used for generalizing MatchRoadsFast3D for arbitrary combinations of
// projections.

class ByProjTypeMap
  : public std::binary_function< Rvec_t, Rvec_t, bool >
{
public:
  ByProjTypeMap( const vec_uint_t& lookup_table ) : fLookup(lookup_table) {}
  bool operator() ( const Rvec_t& a, const Rvec_t& b ) const
  {
    assert( !a.empty() and !b.empty() );
    EProjType ta = a.front()->GetProjection()->GetType();
    EProjType tb = b.front()->GetProjection()->GetType();
    assert( ta != tb );
    assert( fLookup[ta] < 3 and fLookup[tb] < 3 ); // else bug in Tracker::Init
    return ( fLookup[ta] < fLookup[tb] );
  }
private:
  const vec_uint_t& fLookup; // Lookup index proj type -> sort order index
};

//_____________________________________________________________________________
UInt_t Tracker::MatchRoadsFast3D( vector<Rvec_t>& roads, UInt_t /* ncombos */,
				  list< pair<Double_t,Rvec_t> >& combos_found,
				  Rset_t& unique_found )
{
  // Fast MatchRoad algorithm for the special case n==3 and symmetric angles
  // of planes 0 and 1:
  //  - intersect first two proj in front and back
  //  - calc perp distances to 3rd, add in quadrature -> matchval
  // Requires exactly 3 projections
  // Note that the elements of the Rvec_t in the output 'combos_found' are not
  // necessarily in order of ascending projection type.

  vector<Rvec_t>::size_type nproj = roads.size();
  assert( nproj == 3 and fProj.size() == nproj );
  assert( f3dIdx.empty() or
	  (Int_t)f3dIdx.size() > (Int_t)fProj.back()->GetType() );

#ifdef VERBOSE
  if( fDebug > 0 )
    cout << "fast algo):";
#endif

  // Vector holding a combination of roads to test. One road from each
  // projection
  Rvec_t selected( nproj, 0 );

  // Put the road lists in the order of the projection symmetry axes.
  // After sorting, roads[0] and roads[1] correspond to the two projections
  // with symmetric angles around the symmetry axis, represented by roads[2].
  // In the comments below, indices 0, 1 and 2 are referred to as u, v and x,
  // respectively, even though the actual projection types may be different.
  if( !f3dIdx.empty() )
    sort( ALL(roads), ByProjTypeMap(f3dIdx) );

  // Fetch coefficients for coordinate transformations
  const Projection* ip = roads[0].front()->GetProjection();
  Double_t su = ip->GetSinAngle();
  Double_t cu = ip->GetCosAngle();
  assert( ip != roads[1].front()->GetProjection() );
  ip = roads[1].front()->GetProjection();
  Double_t sv = ip->GetSinAngle();
  Double_t cv = ip->GetCosAngle();
  Double_t inv_denom = 1.0/(sv*cu-su*cv);
  // Only the scaled coefficients are needed (cf. Road::Intersect)
  su *= inv_denom; cu *= inv_denom; sv *= inv_denom; cv *= inv_denom;
  assert( ip != roads[2].front()->GetProjection() );
  ip = roads[2].front()->GetProjection();
  assert( ip != roads[0].front()->GetProjection() );
  // Components of the 3rd projection's axis
  Double_t xax_x = ip->GetAxis().X();
  Double_t xax_y = ip->GetAxis().Y();

  // The selected roads from each of the three projections
  Road* tuple[3];
  Plane *front_plane = fPlanes.front(), *back_plane = fPlanes.back();

  // For fast access to the relevant position range, sort the 3rd projection
  // by ascending front position
  sort( ALL(roads[2]), Road::PosIsLess() );
  Road::PosIsNear pos_near( TMath::Sqrt(f3dMatchCut) );

  UInt_t nfound = 0;
  Double_t zback = fPlanes.back()->GetZ();
  Double_t matchval = 0.0;
  // Number of roads in u/v projections
  UInt_t nrd0 = roads[0].size(), nrd1 = roads[1].size();
  // ird0 and ird1 (below) are the indices of the currently selected u/v pair
  UInt_t ird0 = 0;
  // Time-critical loop, may run O(1e5) times per event with noisy input
  while( ird0 != nrd0 ) {
    tuple[0] = roads[0][ird0];
    Double_t uf = tuple[0]->GetPos();
    Double_t ub = tuple[0]->GetPos(zback);
    Double_t usf = uf * sv;
    Double_t ucf = uf * cv;
    Double_t usb = ub * sv;
    Double_t ucb = ub * cv;
    UInt_t ird1 = 0;
    while( ird1 != nrd1 ) {
      tuple[1] = roads[1][ird1];
      Double_t v = tuple[1]->GetPos();
      Double_t xf = usf - v * su;
      Double_t yf = v * cu - ucf;
      if( front_plane->Contains(xf,yf) ) {
	v = tuple[1]->GetPos(zback);
	Double_t xb = usb - v * su;
	Double_t yb = v * cu - ucb;
	if( back_plane->Contains(xb,yb) ) {
	  Double_t pf = xf*xax_x + yf*xax_y; // front x from u/v
	  Double_t pb = xb*xax_x + yb*xax_y; // back x from u/v
	  // Find range of roads in 3rd projection near this front x
	  pair<Rvec_t::iterator,Rvec_t::iterator> range =
	    equal_range( ALL(roads[2]), pf, pos_near );
	  // Test the candidate x-roads for complete matches
	  for( Rvec_t::iterator it=range.first; it != range.second; ++it ) {
	    tuple[2] = *it;
	    Double_t d1 = tuple[2]->GetPos()      - pf;
	    Double_t d2 = tuple[2]->GetPos(zback) - pb;
	    //TODO; weigh with uncertainties?
	    matchval = d1*d1 + d2*d2;
#ifdef VERBOSE
	    if( fDebug > 3 ) {
	      if( matchval < f3dMatchCut || fDebug > 4 ) {
		cout << tuple[0]->GetProjection()->GetName()
		     << tuple[1]->GetProjection()->GetName()
		     << " front = " << "(" << xf << "," << yf << ")" << endl;
		cout << "front " << tuple[2]->GetProjection()->GetName()
		     << " = ("
		     << tuple[2]->GetPos() * xax_x << ","
		     << tuple[2]->GetPos() * xax_y << ")" << endl;
		cout << "front dist = " << d1 << endl;
		cout << tuple[0]->GetProjection()->GetName()
		     << tuple[1]->GetProjection()->GetName()
		     << " back = " << "(" << xb << "," << yb << ")" << endl;
		cout << "back " << tuple[2]->GetProjection()->GetName()
		     << " =  ("
		     << tuple[2]->GetPos(zback) * xax_x << ","
		     << tuple[2]->GetPos(zback) * xax_y << ")" << endl;
		cout << "back dist = " << d2 << endl;
		cout << "matchval = " << matchval*f3dMatchvalScalefact
		     << endl;
	      }
	    }
#endif
	    // Check if match, if so then add
	    if( matchval < f3dMatchCut ) {
	      ++nfound;
	      selected.assign( tuple, tuple+3 );
	      Add3dMatch( selected, matchval, combos_found, unique_found );
	    }
	  }
	}
      }
      ++ird1;
    }
    ++ird0;
  }

  return nfound;
};

//_____________________________________________________________________________
UInt_t Tracker::MatchRoads( vector<Rvec_t>& roads,
			    list< pair<Double_t,Rvec_t> >& combos_found,
			    Rset_t& unique_found )
{
  // Match roads from different projections
  // The input vector 'roads' may be altered unpredictably.
  // Output in 'combos_found' and 'unique_found'

  // The number of projections that we work with
  vector<Rvec_t>::size_type nproj = roads.size();

  combos_found.clear();
  unique_found.clear();

  // Number of all possible combinations of the input roads
  // May overflow for extremely busy events
  UInt_t ncombos;
  bool inrange = true;
  try {
    ncombos = accumulate( ALL(roads), (UInt_t)1, SizeMul<Rvec_t>() );
  }
  catch( overflow_error ) {
    ncombos = 0;
    inrange = false;
  }

#ifdef VERBOSE
  if( fDebug > 0 ) {
    if( inrange )
      cout << "Matching ";
    else
      cout << "Too many combinations trying to match ";
    for( vector<Rvec_t>::size_type i = 0; i < nproj; ++i ) {
      cout << roads[i].size();
      if( i+1 < nproj ) cout << "x";
    }
    cout << " track projections in 3D";
    if( inrange ) {
      cout << "(" << ncombos << " combination";
      if( ncombos != 1 ) cout << "s";
      cout << ", ";
    } else {
      cout << ". Giving up.";
      cout << endl;
    }
  }
#endif
#ifdef TESTCODE
  fNcombos = ncombos;
#endif

  if( ncombos == 0 ) {
    if( inrange )
      fTrkStat = kNoRoadCombos;   // bug?
    else
      fTrkStat = kTooManyRoadCombos;
    return 0;
  }

  UInt_t nfound = MatchRoadsImpl( roads, ncombos, combos_found, unique_found );

#ifdef VERBOSE
  if( fDebug > 0 ) {
    cout << nfound << " match";
    if( nfound != 1 )
      cout << "es";
    cout << " found" << endl;
  }
#endif
  return nfound;
}

//_____________________________________________________________________________
inline
Bool_t Tracker::PassTrackCuts( const FitRes_t& fit_par ) const
{
  // Test results of 3D track fit (in fit_par) against hard cuts from
  // database

  if( fit_par.ndof < fMinNdof ){
    //cout<<"ndof failed: "<<fit_par.ndof<<endl;
    return false;
  }

  if( TestBit(kDoChi2) ) {
    pdbl_t chi2_interval = GetChisqLimits(fit_par.ndof);
    if( fit_par.chi2 < chi2_interval.first or
	fit_par.chi2 > chi2_interval.second )
      return false;
  }

  return true;
}

//_____________________________________________________________________________
Int_t Tracker::CoarseTrack( TClonesArray& tracks )
{
  // Find tracks from the hitpatterns, using the coarse hit drift times
  // uncorrected for track slope, timing offset, fringe field effects etc.

  //  static const char* const here = "CoarseTrack";

  if( fTrkStat != 0 )
    return -1;

  if( !TestBit(kDoCoarse) ) {
    fTrkStat = kNoTrackingRequested;
    return 0;
  }

#ifdef TESTCODE
  TStopwatch timer, timer_tot;
#endif

  Int_t err = 0;
  if( fMaxThreads > 1 ) {
    err = fThreads->Run( fMaxThreads );
  } else {
    // Single-threaded execution: Track() each projection in turn
    for( vpiter_t it = fProj.begin(); it != fProj.end(); ++it ) {
      Int_t ret = (*it)->Track();
      if( ret < 0 )
	err = 1;
    }
  }
  // Abort on error (e.g. too many patterns)
  if( err != 0 ) {
    fTrkStat = kProjTrackError;
    return -1;
  }
  // Copy pointers to roads from each projection into local 2D vector
  // (projections, roads).
  // NB: subsequent code assumes that the first index runs in order
  // of ascending projection type, so this can't easily be done in the
  // tracking threads.
  vector<Rvec_t>::size_type nproj = 0;
  vector<Rvec_t> roads;
  roads.reserve( fProj.size() );
  UInt_t found_types = 0;
  for( vpiter_t it = fProj.begin(); it != fProj.end(); ++it ) {
    Projection* proj = *it;
    EProjType type = proj->GetType();

    Int_t nrd = proj->GetNgoodRoads();
    // cout<<nrd<<" goodroads "<<endl;
    if( nrd > 0 ) {
      // Count number of projections with at least one road
      ++nproj;
      SETBIT( found_types, type );
      // Add pointers to the good roads to the appropriate vector
      roads.push_back( Rvec_t() );
      roads.back().reserve(nrd);
      for( UInt_t i = 0; i < proj->GetNroads(); ++i ) {
	Road* rd = proj->GetRoad(i);
	assert(rd);
	if( rd->IsGood() )
	  roads.back().push_back(rd);
      }
    }
  }
  assert( roads.size() == nproj );
#ifdef TESTCODE
  t_track = 1e6*timer.RealTime();
  timer.Start();
#endif

  int numberOfTracks = 0;
  // Combine track projections to 3D tracks
  if( nproj >= fMinReqProj ) {
    // Vector holding the results (vectors of roads with good matchval)
    list< pair<Double_t,Rvec_t> > road_combos;
    // Set of the unique roads occurring in the road_combos elements
    Rset_t unique_found;




    // Find matching combinations of roads
    UInt_t nfits = MatchRoads( roads, road_combos, unique_found );
    //cout<<"##@@ nfits: "<<nfits<<endl;
    

#ifdef TESTCODE
    t_3dmatch = 1e6*timer.RealTime();
    timer.Start();
#endif
    // Fit each set of matched roads using linear least squares, yielding
    // the 3D track parameters, x, x'(=mx), y, y'(=my)
    FitRes_t fit_par;
    fit_par.coef.reserve(4);
    if( nfits == 1 ) {
      // If there is only one combo (typical case), life is simple:
      fit_par.matchval    = road_combos.front().first;
      Rvec_t& these_roads = road_combos.front().second;
      fit_par.roads       = &these_roads;
      fit_par.ndof = FitTrack( these_roads, fit_par.coef, fit_par.chi2 );
      if( fit_par.ndof > 0 ) {
	if( PassTrackCuts(fit_par) ){
	  numberOfTracks++;
	  NewTrack( tracks, fit_par );
	}
	else
	  fTrkStat = kFailedTrackCuts;
      }
      else
	FitErrPrint( fit_par.ndof );
    }
    else if( nfits > 1 ) {
      // For multiple road combinations, find the set of tracks with the
      // lowest chi2s that uses each road at most once
      typedef map<Rset_t, FitRes_t> FitResMap_t;
      FitResMap_t fit_results;
      multimap< TrackFitWeight, Rset_t > fit_chi2;
      // Fit all combinations and sort the results by ascending chi2
      for( list< pair<Double_t,Rvec_t> >::iterator it = road_combos.begin();
	   it != road_combos.end(); ++it ) {
	fit_par.matchval    = (*it).first;
	Rvec_t& these_roads = (*it).second;
	fit_par.roads       = &these_roads;
	fit_par.ndof = FitTrack( these_roads, fit_par.coef, fit_par.chi2 );
	if( fit_par.ndof > 0 ) {
	  if( PassTrackCuts(fit_par) ) {
	    Rset_t road_tuple( ALL(these_roads) );
	    pair<FitResMap_t::iterator,bool> ins =
	      fit_results.insert( make_pair(road_tuple, fit_par) );
	    assert( ins.second );
	    fit_chi2.insert( make_pair( TrackFitWeight(fit_par), road_tuple) );
	  }
	} else
	  FitErrPrint( fit_par.ndof );
      }

#ifdef VERBOSE
      if( fDebug > 2 ) {
	cout << "Track candidates:" << endl;
	for( multimap<TrackFitWeight,Rset_t>::iterator it = fit_chi2.begin();
	     it != fit_chi2.end(); ++it ) {
	  FitResMap_t::iterator found = fit_results.find((*it).second);
	  FitRes_t& r = (*found).second;
	  cout 	 << "ndof = " << r.ndof
		 << " rchi2 = " << r.chi2/(double)r.ndof
		 << " x/y = " << r.coef[0] << "/" << r.coef[2]
		 << " mx/my = " << r.coef[1] << "/" << r.coef[3]
		 << endl;
	}
      }
#endif
      if( !fit_chi2.empty() ) {
	// Select "optimal" set of roads, minimizing sum of chi2s
	vector<Rset_t> best_roads;
	OptimalN( unique_found, fit_chi2, best_roads,
		  AnySharedHits(), CheckTypes(found_types) );

	if( best_roads.empty() )
	  fTrkStat = kFailedOptimalN;


	// Now each selected road tuple corresponds to a new track
	for( vector<Rset_t>::iterator it = best_roads.begin(); it !=
	       best_roads.end(); ++it ) {
	  // Retrieve the fit results for this tuple
	  FitResMap_t::iterator found = fit_results.find(*it);
	  assert( found != fit_results.end() );
	  numberOfTracks++;
	  NewTrack( tracks, (*found).second );
	}
      }
      else {
	fTrkStat = kFailedTrackCuts;
      }
    }
    else {
      fTrkStat = kFailed3DMatch;
    } //if(nfits)

#ifdef TESTCODE
    fN3dFits = nfits;
    t_3dfit = 1e6*timer.RealTime();
#endif
#ifdef VERBOSE
    if( fDebug > 0 ) {
      Int_t ntr = tracks.GetLast()+1;
      cout << ntr << " track";
      if( ntr != 1 ) cout << "s";
      if( nfits > 1 ) cout << " after chi2 optimization";
      if( fDebug > 1 and ntr > 0 ) {
	cout << ":" << endl;
	for( Int_t i = 0; i < ntr; ++i ) {
	  THaTrack* tr = (THaTrack*)tracks.UncheckedAt(i);
	  cout << "3D track:  x/y = " << tr->GetX() << "/" << tr->GetY()
	       << " mx/my = " << tr->GetTheta() << "/" << tr->GetPhi()
	       << " ndof = "  << tr->GetNDoF()
	       << " rchi2 = " << tr->GetChi2()/(double)tr->GetNDoF()
	       << endl;
	}
      } else
	cout << endl;
    }
#endif
  }
  else {
    fTrkStat = kTooFewProj;
#ifdef VERBOSE
    if( fDebug > 0 ) {
      cout << "No tracks, number of projections with roads = " << nproj;
      if( nproj > 0 ) {
	// Write out the names of the active projections
	cout << " (";
	UInt_t ndone = 0;
	for( EProjType k = kTypeBegin; k < kTypeEnd; ++k ) {
	  if( TESTBIT(found_types,k) ) {
	    cout << kProjParam[k].name;
	    if( ++ndone != nproj ) cout << ",";
	  }
	}
	assert( ndone == nproj );
	cout << ")";
      }
      cout << ", >=" << fMinReqProj << " required" << endl;
    }
#endif
  } //if(nproj>=fMinReqProj)

#ifdef TESTCODE
    t_coarse = 1e6*timer_tot.RealTime();
#endif
  // Quit here to let detectors CoarseProcess() the approximate tracks,
  // so that they can determine the corrections that we need when we
  // continue in FineTrack
  //cout<<"corseTrack ENds: number of tracks: "<<numberOfTracks<<endl;
  
  return 0;
}

//_____________________________________________________________________________
Int_t Tracker::FineTrack( TClonesArray& /* tracks */ )
{
  // Second-level tracking, applying fine corrections.
  //
  // - Correct the hit drift distances using timing, track slope,
  //   track momentum, magnet fringe field, etc., as applicable.
  // - Re-collect hits of roads and re-fit roads
  // - Re-combine track projections in 3D and re-fit 3D track

  if( !TestBit(kDoFine) )
    return 0;

  if( fTrkStat != 0 )
    return -1;

  // TODO

  return 0;;
}

//_____________________________________________________________________________
Int_t Tracker::DefineVariables( EMode mode )
{
  // Initialize global variables

  if( mode == kDefine and fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list

  RVarDef vars[] = {
    { "trkstat", "Track reconstruction status",  "fTrkStat" },
    { 0 },
  };
  DefineVarsFromList( vars, mode );

  // Define tracking-related variables only if doing tracking
  if( TestBit(kDoCoarse) ) {
    RVarDef vars_tracking[] = {
#ifdef TESTCODE
      { "ncombos",    "Number of road combinations",        "fNcombos" },
      { "nfits",      "Number of 3D track fits done",       "fN3dFits" },
      { "t_track",    "Time in 1st stage tracking (us)",    "t_track" },
      { "t_3dmatch",  "Time in MatchRoads (us)",            "t_3dmatch" },
      { "t_3dfit",    "Time fitting/selecting 3D tracks (us)", "t_3dfit" },
      { "t_coarse",   "Total time in CoarseTrack (us)",     "t_coarse" },
#endif
      { 0 }
    };
    DefineVarsFromList( vars_tracking, mode );
  }
  return 0;
}

//_____________________________________________________________________________
Projection* Tracker::MakeProjection( EProjType type, const char* name,
				     Double_t angle,
				     THaDetectorBase* parent ) const
{
  // Create an object of the projection class used by this implementation.
  // The type of projection determines basically the method with which the
  // Hitpattern will be filled.

  return new Projection( type, name, angle, parent );
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus Tracker::Init( const TDatime& date )
{
  // Initialize GEM. Calls standard Init(), then initializes subdetectors.

  static const char* const here = "Init";

  // The "cratemap" is only needed during Init of Tracker and Plane
  assert( fCrateMap == 0 );
  fCrateMap = new THashTable(100);
  fCrateMap->SetOwner();

  // Initialize ourselves. This calls our ReadDatabase() and DefineVariables()
  EStatus status = THaTrackingDetector::Init(date);

  // Initialize the planes
  UInt_t nplanes = 0;  // number of active planes
  if( status == kOK ) {
    vrsiz_t iplane = 0;
    for( ; iplane < fPlanes.size(); ++iplane ) {
      try {
	status = fPlanes[iplane]->Init(date);
	if( status )
	  break;
      }
      catch( exception& e ) {
	Error( Here(here), "Caught %s when initializing plane \"%s\". "
	       "Call expert.", e.what(), fPlanes[iplane]->GetName() );
	status = kInitError;
	break;
      }
      if( !fPlanes[iplane]->IsDummy() ) {
	++nplanes;
      }
    }
  }
  delete fCrateMap; fCrateMap = 0;
  if( status )
    return fStatus = status;

  // Sort planes by increasing z-position
  sort( ALL(fPlanes), Plane::ZIsLess() );
  // FIXED: this works now independently of the number of planes per tracker
  // because we allow an element a to be compared with itself.

  // Calculate track fitting parameters. This needs to be done here because
  // we need to know the number of active (non-dummy) planes, which is not
  // available until after Plane initialization.
  if( TestBit(kDoCoarse) ) {   // only needed if doing tracking
    if( nplanes < 5 ) {
      Error( Here(here), "Insufficient number of planes = %u. Need at least 5. "
	     "Fix database.", nplanes );
      return kInitError;
    }
    if( fDBmaxmiss+5 > static_cast<Int_t>(nplanes) ) {
      Error( Here(here), "3d_maxmiss = %d (number of allowed planes with "
	     "missing hits) too large. Must be %u or less. Fix database.",
	     fDBmaxmiss, nplanes-5 );
      return kInitError;
    }
    // Convert maximum number of missing hits to ndof of fits
    fMinNdof = ( fDBmaxmiss >= 0 ) ? nplanes - 4 - fDBmaxmiss : 1;

    // Determine Chi2 confidence interval limits for the selected CL and the
    // possible degrees of freedom of the 3D track fit
    if( TestBit(kDoChi2) ) {
      if( fDBconf_level < 0.0 || fDBconf_level > 1.0 ) {
	Error( Here(here), "Illegal fit confidence level = %lf. "
	       "Must be 0-1. Fix database.", fDBconf_level );
	return kInitError;
      }
      fChisqLimits.clear();
      fChisqLimits.resize( nplanes-3, make_pair<Double_t,Double_t>(0,0) );
      for( vec_pdbl_t::size_type dof = fMinNdof; dof < fChisqLimits.size();
	   ++dof ) {
	fChisqLimits[dof].first  = TMath::ChisquareQuantile( fDBconf_level, dof );
	fChisqLimits[dof].second =
	  TMath::ChisquareQuantile( 1.0-fDBconf_level, dof );
      }
    }
  }

  // Find partners for the planes, if appropriate. Details depend on the type
  // of detector.
  status = PartnerPlanes();
  if( status )
    return fStatus = status;
  cout<<"going to "<<endl;
  // Set up the projections based on which plane types are defined
  assert( fProj.empty() );
  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    Plane* thePlane = fPlanes[iplane];
    assert( !thePlane->GetProjection() );
    EProjType type = thePlane->GetType();
    cout<<fProj.size()<<" plane type: "<<type<<endl;
    vpiter_t it = find_if( ALL(fProj), Projection::TypeEquals(type) );
    Projection* proj;
    if( it != fProj.end() ) {
      proj = *it;
    } else {
      try {
	proj = MakeProjection( type,
			       kProjParam[type].name,
			       kProjParam[type].angle*TMath::DegToRad(),
			       this
			       );
      }
      catch( bad_alloc ) {
	Error( Here(here), "Out of memory creating projection \"%s\". "
	       "Call expert.", kProjParam[type].name );
      }
      catch( Projection::bad_angle ) {
	// Error message is printed by Projection constructor in this case
	delete proj; proj = 0;
      }
      catch( exception& e ) {
	Error( Here(here), "Caught %s when creating projection \"%s\". "
	       "Call expert.", e.what(), kProjParam[type].name );
	delete proj; proj = 0;
      }
      if( proj && proj->IsZombie() ) {
	delete proj; proj = 0;
      }
      if( !proj )
	return fStatus = kInitError;

      proj->SetDebug( fDebug ); // Projections inherit parent's debug level
      fProj.push_back( proj );
    }
    if( thePlane->IsDummy() )
      proj->AddDummyPlane( thePlane );
    else
      proj->AddPlane( thePlane );
  }

  if( fDebug > 0 )
    Info( Here(here), "Set up %u projections", (UInt_t)fProj.size() );

  // Sort projections by ascending EProjType
  sort( ALL(fProj), Projection::ByType() );
  
  // Initialize the projections. This will read the database and set
  // the projections' angle, width and maxslope
  try {
    for( vpiter_t it = fProj.begin(); it != fProj.end(); ++it ) { 
      status = (*it)->Init(date);  
      if( status )
	return fStatus = status;
    }
  }
  catch(...) {
    Error( Here(here), "Failed to initialize projections. Call expert." );
    return fStatus = kInitError;
  }

  // Sanity check of the projection angles
  for( vpiter_t it = fProj.begin(); it != fProj.end(); ++it ) {
    Projection* theProj = *it;
    // Ensure minimum angle differences
    vpiter_t jt(it);
    ++jt;
    while( jt != fProj.end() ) {
      Projection* otherProj = *jt;
      Double_t dp = theProj->GetAxis().DeltaPhi( otherProj->GetAxis() );
      if( TMath::Abs(dp) < fMinProjAngleDiff or
	  TMath::Abs(dp-TMath::Pi()) < fMinProjAngleDiff ) {
	Error( Here(here), "Projection misconfiguration: \"%s\"-angle (%6.2lf) "
	       "and \"%s\"-angle (%6.2lf) are too colinear. Review geometry.",
	       kProjParam[theProj->GetType()].name,
	       theProj->GetAngle()*TMath::RadToDeg(),
	       kProjParam[otherProj->GetType()].name,
	       otherProj->GetAngle()*TMath::RadToDeg() );
	return fStatus = kInitError;
      }
      ++jt;
    }
  }

  // Check if we can use the simplified 3D matching algorithm
  //TODO: generalize to 4 projections with symmetry?
  vector<ProjAngle_t>::size_type nproj = fProj.size();
  if( nproj == 3 ) {
    // Algorithm:
    // (1) Find the minimum angle spanned by the three axes. We need to
    //     check which combination of axis flips by 180 degrees gives the
    //     smallest value. (Remember, only the axis orientation matters, not
    //     the sign, so axes at angle and angle+/-180 are equivalent.)
    // (2) Sort the projections by (flipped) angle.
    // (3) Check if the max/min angles are symmetric wrt the center one.
    vector<ProjAngle_t> proj_angle( ALL(fProj) );
    assert( proj_angle.size() == nproj );
    Double_t smin = kBig;
    UInt_t imin = kMaxUInt, ncombos = 1U<<(nproj-1);
    for( UInt_t ix = 0; ix < ncombos; ++ix ) {
      Double_t amin = kBig, amax = -kBig;
      for( vector<ProjAngle_t>::size_type j=0; j<nproj; ++j ) {
	// The set bits in ix indicate a trial flip
	amin = TMath::Min( amin, proj_angle[j].angle(TESTBIT(ix,j)) );
	amax = TMath::Max( amax, proj_angle[j].angle(TESTBIT(ix,j)) );
      }
      assert( amin + 2.0*fMinProjAngleDiff <= amax );
      Double_t span = amax - amin;
      if( span < smin ) {
	smin = span;
	imin = ix;
      }
    }
    assert( imin < ncombos );
    // Actually flip the angles to the arrangement with smallest span
    // so we can sort by angle. Note that the last element is never
    // flipped, without loss of generality
    assert( !TESTBIT(imin,nproj-1) );
    for( vector<ProjAngle_t>::size_type j=0; j<nproj-1; ++j ) {
      if( TESTBIT(imin,j) )
	proj_angle[j].m_angle = proj_angle[j].angle(true);
    }
    sort( ALL(proj_angle) );    // Re-sort by (flipped) angle
    // The abs(angle) of the min/max projection axes wrt to the center one
    // must be (nearly) the same for the fast_3d algorithm to apply
    Double_t d1 = proj_angle.back().m_angle - proj_angle[1].m_angle;
    Double_t d2 = proj_angle[1].m_angle - proj_angle.front().m_angle;
    assert( d1 > 0 and d2 > 0 ); // else not sorted correctly
    if( TMath::Abs(d1-d2) < 0.5*TMath::DegToRad() ) {
      SetBit(k3dFastMatch);
      Double_t tan = TMath::Tan( 0.5 * (TMath::Pi()-(d1+d2)) );
      // The scale factor converts the fast_3d matchvalue to the one computed
      // by the general algorithm
      f3dMatchvalScalefact = 2.0 * (1.0/3.0 + tan*tan);
      // Avoid scaling for every event
      f3dMatchCut /= f3dMatchvalScalefact;
      // Keep a lookup table of the sorted projection order, with the
      // symmetry axis last
      EProjType
	t0 = proj_angle.front().m_proj->GetType(),
	t1 = proj_angle.back().m_proj->GetType(),
	t2 = proj_angle[1].m_proj->GetType();
      if( t0 < t1 and t1 < t2 ) { // already in the right order
	f3dIdx.clear();
      } else {
	f3dIdx.assign( kTypeEnd-kTypeBegin, kMaxUInt );
	f3dIdx[t0] = 0;
	f3dIdx[t1] = 1;
	f3dIdx[t2] = 2;
      }
      if( fDebug > 0 )
	Info( Here(here), "Enabled fast 3D projection matching" );
    }
  }

  // If threading requested, load thread library and start up threads
  if( fMaxThreads > 1 ) {
    delete fThreads; fThreads = 0;
    if( gSystem->Load("libThread") >= 0 ) {
      fThreads = new ThreadCtrl( fProj );
    } else {
      // Error loading library
      Warning( Here(here), "Error loading thread library. Falling back to "
	       "single-threaded processing." );
      fMaxThreads = 1;
    }
  }

  // Keep a simple flag for the rotation status for efficiency.
  fIsRotated = !fRotation.IsIdentity();

#ifdef MCDATA
  // Set up handler for MC data
  if( TestBit(kMCdata) ) {
    delete fMCPointUpdater;
    fMCPointUpdater = new MCPointUpdater;
  }
#endif

  return fStatus = kOK;
}

//_____________________________________________________________________________
Int_t Tracker::ReadDatabase( const TDatime& date )
{
  // Read generic Tracker keys from database

  static const char* const here = "ReadDatabase";

  fIsInit = kFALSE;
  // Delete existing configuration (in case we are re-initializing)
  DeleteContainer( fProj );
  DeleteContainer( fPlanes );
  fCalibPlanes.clear();

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  // Read fOrigin (detector position) and fSize. fOrigin is the position of
  // the Tracker relative to some superior coordinate system
  // (typically the spectrometer detector stack reference frame).
  // fOrigin will be added to all tracks generated; if fOrigin.Z() is not
  // zero, tracks will be projected into the z=0 plane.
  fOrigin.SetXYZ(0,0,0);
  fRotation.SetToIdentity();
  fInvRot.SetToIdentity();
  Int_t err = ReadGeometry( file, date );
  if( err ) {
    fclose(file);
    return err;
  }

  // Putting this container on the stack may cause strange stack overflows!
  vector<vector<Int_t> > *cmap = new vector<vector<Int_t> >;

  string planeconfig, calibconfig;
  f3dMatchCut = 1e-4;
  Int_t event_display = 0, disable_tracking = 0,
    disable_finetrack = 0, disable_chi2 = 0, proj_to_z0 = 1;
#ifdef MCDATA
  Int_t mc_data = 0;
#endif
  Int_t maxthreads = -1;
  fDBmaxmiss = -1;
  fDBconf_level = 1e-9;
  ResetBit( k3dFastMatch ); // Set in Init()
  assert( GetCrateMapDBcols() >= 5 );
  DBRequest request[] = {
    { "planeconfig",       &planeconfig,       kString },
    { "cratemap",          cmap,               kIntM,   GetCrateMapDBcols() },
#ifdef MCDATA
    { "MCdata",            &mc_data,           kInt,    0, 1 },
#endif
    { "3d_matchcut",       &f3dMatchCut,       kDouble, 0, 1 },
    { "event_display",     &event_display,     kInt,    0, 1 },
    { "disable_tracking",  &disable_tracking,  kInt,    0, 1 },
    { "disable_finetrack", &disable_finetrack, kInt,    0, 1 },
    { "proj_to_z0",        &proj_to_z0,        kInt,    0, 1 },
    { "calibrate",         &calibconfig,       kString, 0, 1 },
    { "3d_maxmiss",        &fDBmaxmiss,        kInt,    0, 1 },
    { "3d_chi2_conflevel", &fDBconf_level,     kDouble, 0, 1 },
    { "3d_disable_chi2",   &disable_chi2,      kInt,    0, 1 },
    { "maxthreads",        &maxthreads,        kInt,    0, 1 },
    { 0 }
  };

  Int_t status = kInitError;
  err = LoadDB( file, date, request, fPrefix );
  fclose(file);

  if( !err ) {
    if( cmap->empty() ) {
      Error(Here(here), "No cratemap defined. Set \"cratemap\" in database.");
    } else {
      // Build the list of crate map elements
      for( vviter_t it = cmap->begin(); it != cmap->end(); ++it ) {
	vector<int>& row = *it;
	assert( row.size() == GetCrateMapDBcols() ); // failure indicates bug in LoadDB
	for( Int_t slot = row[1]; slot <= row[2]; ++slot ) {
	  DAQmodule* m = 0;
	  if( GetCrateMapDBcols() < 6 )
	    m = new DAQmodule( row[0], slot, row[3], row[4] );
	  else
	    m = new DAQmodule( row[0], slot, row[3], row[4], row[5] );
	  DAQmodule* found = static_cast<DAQmodule*>(fCrateMap->FindObject(m));
	  if( found ) {
	    m->Copy(*found);  // Later entries override earlier ones
	    delete m;
	  }
	  else
	    fCrateMap->Add(m);
	}
      }
      status = kOK;
    }
  }
  delete cmap; cmap = 0;
  if( status != kOK )
    return status;

  // Set common analysis control flags
#ifdef MCDATA
  SetBit( kMCdata,        mc_data );
#endif
  SetBit( kEventDisplay,  event_display );
  SetBit( kDoCoarse,      !disable_tracking );
  SetBit( kDoFine,        !(disable_tracking or disable_finetrack) );
  SetBit( kDoChi2,        !disable_chi2 );
  SetBit( kProjTrackToZ0, proj_to_z0 );

  cout << endl;
  if( fDebug > 0 ) {
#ifdef MCDATA
    Info( Here(here), "Tracker flags mcdata/evdisp/coarse/fine/projz0 = "
	  "%d/%d/%d/%d/%d", TestBit(kMCdata), TestBit(kEventDisplay),
	  TestBit(kDoCoarse), TestBit(kDoFine), TestBit(kProjTrackToZ0) );
#else
    Info( Here(here), "Tracker flags evdisp/coarse/fine/projz0 = "
	  "%d/%d/%d/%d", TestBit(kEventDisplay), TestBit(kDoCoarse),
	  TestBit(kDoFine), TestBit(kProjTrackToZ0) );
#endif
  }

  vector<string> planes = vsplit(planeconfig);
  if( planes.empty() ) {
    Error( Here(here), "No planes defined. Set \"planeconfig\" in database." );
    return kInitError;
  }
  vector<string> calibplanes = vsplit(calibconfig);

  // Set up the wire/readout planes
  for( ssiz_t i=0; i<planes.size(); ++i ) {
    assert( !planes[i].empty() );
    const char* name = planes[i].c_str();
    vriter_t it = find_if( ALL(fPlanes), Plane::NameEquals( name ) );
    if( it != fPlanes.end() ) {
      Error( Here(here), "Duplicate plane name: %s. Fix database.", name );
      return kInitError;
    }
    Plane* newplane = 0;
    try { newplane = MakePlane( name, name, this ); }
    catch( ... ) { newplane = 0; }
    if( !newplane or newplane->IsZombie() ) {
      // Urgh. Something is very bad
      Error( Here(here), "Error creating readout plane %s. Call expert.",
	     name );
      delete newplane;
      return kInitError;
    }
    
    fPlanes.push_back( newplane );
    newplane->SetDebug( fDebug );
    newplane->SetDefinedNum( i );

    // Set calibration mode if requested for any planes
    if( !calibplanes.empty() ) {
      vector<string>::iterator its =
	find_if( ALL(calibplanes), bind2nd(equal_to<string>(), planes[i]) );
      if( its != calibplanes.end() ) {
	Info( Here(here), "Plane %s in calibration mode", name );
	newplane->EnableCalibration();
	fCalibPlanes.push_back( newplane );
	calibplanes.erase( its );
      }
    }
  }

  // Warn if any requested calibration plane(s) do not exist
  if( !calibplanes.empty() ) {
    string s("plane");
    if( calibplanes.size() > 1 )
      s.append("s");
    for( vector<string>::size_type i = 0; i < calibplanes.size(); ++i ) {
      s.append(" ");
      s.append(calibplanes[i]);
    }
    Warning( Here(here), "Requested calibration for undefined %s. "
	     "Error in database?", s.c_str() );
  }

  if( fDebug > 0 )
    Info( Here(here), "Loaded %u planes", static_cast<UInt_t>(fPlanes.size()) );

  // Determine maximum number of threads to run. maxthreads from the database
  // has priority. maxthreads = 0 or negative indicates that the number of
  // CPUs/cores of the current host should be used. If not available, use 1.
  // To ensure single-threaded processing, set maxthreads = 1 in the database.
  // The number of threads is capped by the maximum that the ThreadCtrl class
  // supports - currently 31 due to the size of UInt_t, which is more than
  // enough since we never run more threads than the number of projections.
  bool warn = false;
  if( maxthreads > 0 )
    fMaxThreads = maxthreads;
  else {
    SysInfo_t sysifo;
    gSystem->GetSysInfo( &sysifo );
    if( sysifo.fCpus > 0 )
      fMaxThreads = sysifo.fCpus;
    else {
      warn = true;
      fMaxThreads = 1;
    }
  }
  // Limit number of threads to what ThreadCtrl can support (currently = 31)
  fMaxThreads = min(fMaxThreads, ThreadCtrl::GetMaxThreads());
  if( warn )
    Warning( Here(here), "Cannot determine number of CPU cores. "
	     "Falling back to single-threaded processing." );
  else if( fDebug > 0 )
    Info( Here(here), "Enabled up to %u threads", fMaxThreads );

  fIsInit = kTRUE;
  return kOK;
}

//_____________________________________________________________________________
void Tracker::Print(const Option_t* opt) const
{
  THaTrackingDetector::Print(opt);

  Int_t verbose = 0;
  if( opt ) {
    TString opt_s(opt);
    opt_s.ToLower();
    verbose = opt_s.CountChar('v');
  }

  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane ) {
    fPlanes[iplane]->Print(opt);
    if( iplane+1 == fPlanes.size() )
      cout << endl;
  }

  for( vpsiz_t k = 0; k < fProj.size(); ++k ) {
    fProj[k]->Print(opt);
    if( verbose > 0 and k+1 != fProj.size() )
      cout << endl;
  }

}


//_____________________________________________________________________________
void Tracker::SetDebug( Int_t level )
{
  // Set debug level of this detector, including all planes (subdetectors)
  // and projections

  THaTrackingDetector::SetDebug( level );

  for( vrsiz_t iplane = 0; iplane < fPlanes.size(); ++iplane )
    fPlanes[iplane]->SetDebug( level );

  for( vpsiz_t k = 0; k < fProj.size(); ++k )
    fProj[k]->SetDebug( level );
}

//_____________________________________________________________________________
void Tracker::EnableEventDisplay( Bool_t b )
{
  // Enable event display support. Can only be called before initialization.

  if( !fIsInit ) {
    Error( Here("EnableEventDisplay"), "Cannot enable/disable event display "
	   "support after initialization." );
    return;
  }

  SetBit( kEventDisplay, b );
}

//_____________________________________________________________________________
inline
static DAQmodule* FindDAQmodule( UShort_t crate, UShort_t slot,
				 const THashTable* table )
{
  assert(table);
  CSpair m( crate, slot );
  return static_cast<DAQmodule*>( table->FindObject(&m) );
}

//_____________________________________________________________________________
UInt_t Tracker::LoadDAQmodel( THaDetMap::Module* mod ) const
{
  // Update detector map module 'mod' with the model number from the cratemap
  DAQmodule* found = FindDAQmodule( mod->crate, mod->slot, fCrateMap );
  UInt_t num = found ? found->fModel : 0;
  mod->SetModel( num );
  return num;
}

//_____________________________________________________________________________
Double_t Tracker::LoadDAQresolution( THaDetMap::Module* mod ) const
{
  // Update detector map module 'mod' with the resolution from the cratemap
  DAQmodule* found = FindDAQmodule( mod->crate, mod->slot, fCrateMap );
  Double_t res = found ? found->fResolution : 0.0;
  mod->SetResolution( res );
  return res;
}

//_____________________________________________________________________________
UInt_t Tracker::GetDAQnchan( THaDetMap::Module* mod ) const
{
  // Return number of channels for detector map module 'mod' from cratemap
  DAQmodule* found = FindDAQmodule( mod->crate, mod->slot, fCrateMap );
  return found ? found->fNchan : 0;
}

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

ClassImp(TreeSearch::Tracker)
