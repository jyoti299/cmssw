#ifndef TrackingTools_TransientTrack_TrackingTransientTrack_h
#define TrackingTools_TransientTrack_TrackingTransientTrack_h
#include <tuple>
#include <atomic>
#include <cmath>
/**
   * Concrete implementation of the TransientTrack for a reco::Track
   */

#include "TrackingTools/TransientTrack/interface/BasicTransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "DataFormats/VertexReco/interface/MtdtimeHostCollection.h"
using MTDsoaElements = MtdtimeSoAConstView; //MtdtimeHostCollection;
using namespace std;
namespace reco {
  class TrackTransientTrack : public Track, public BasicTransientTrack {
  public:
    // constructor from persistent track
    TrackTransientTrack();
    TrackTransientTrack(const Track& tk, const MagneticField* field);
    TrackTransientTrack(const Track& tk, const double time, const double dtime, const MagneticField* field);

    TrackTransientTrack(const TrackRef& tk, const MagneticField* field);
    TrackTransientTrack(const TrackRef& tk, const double time, const double dtime, const MagneticField* field);
    
    TrackTransientTrack(const Track& tk, const MTDsoaElements& soaelement, const double time,
                        const double dtime,
                        const MagneticField* field);
                       

    TrackTransientTrack(const TrackRef& tk,
                        const MagneticField* field,
                        const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);
    TrackTransientTrack(const TrackRef& tk,
                        const double time,
                        const double dtime,
                        const MagneticField* field,
                        const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);

    TrackTransientTrack(const Track& tk,
                        const MagneticField* field,
                        const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);
    TrackTransientTrack(const Track& tk,
                        const double time,
                        const double dtime,
                        const MagneticField* field,
                        const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);

    TrackTransientTrack(const TrackTransientTrack& tt);

    TrackTransientTrack& operator=(const TrackTransientTrack& tt);

    void setTrackingGeometry(const edm::ESHandle<GlobalTrackingGeometry>&) override;
 //   std::tuple<MTDsoaElements&, int32_t> getMTDtimeElements(const MTDsoaElements&) override; 
    //auto  getMTDtimeElements(const MTDsoaElements&); 
    void setBeamSpot(const reco::BeamSpot& beamSpot) override;
    FreeTrajectoryState initialFreeState() const override { return initialFTS; }

    TrajectoryStateOnSurface outermostMeasurementState() const override;

    TrajectoryStateOnSurface innermostMeasurementState() const override;

    TrajectoryStateClosestToPoint trajectoryStateClosestToPoint(const GlobalPoint& point) const override {
      return builder(initialFTS, point);
    }

    /**
    * The TSOS at any point. The initial state will be used for the propagation.
    */
    TrajectoryStateOnSurface stateOnSurface(const GlobalPoint& point) const override;

    TrajectoryStateClosestToPoint impactPointTSCP() const override;

    TrajectoryStateOnSurface impactPointState() const override;

    bool impactPointStateAvailable() const override { return (m_TSOS.load() == kSet ? true : false); }

    /**
   * access to original persistent track
   */
    TrackRef persistentTrackRef() const { return tkr_; }

    TrackBaseRef trackBaseRef() const override { return TrackBaseRef(tkr_); }

    TrackCharge charge() const override { return Track::charge(); }
    const MagneticField* field() const override { return theField; }
    //const MTDsoaElements& soa() const override {return soavalue; } 
    const Track& track() const override { return *this; }
   
    TrajectoryStateClosestToBeamLine stateAtBeamLine() const override;
    double timeExt() const override { return (hasTime ? timeExt_ : std::numeric_limits<double>::quiet_NaN()); }
    double dtErrorExt() const override { return (hasTime ? dtErrorExt_ : std::numeric_limits<double>::quiet_NaN()); }

  private:
    TrackRef tkr_;
    bool hasTime;
    double timeExt_, dtErrorExt_;
    const MagneticField* theField;
    FreeTrajectoryState initialFTS;
    // mutable member data, those should be treated very carefully to guarantee
    // thread safeness of the code by using atomic thread-safe helpers, see below
    mutable TrajectoryStateOnSurface initialTSOS;
    mutable TrajectoryStateClosestToPoint initialTSCP;
    mutable TrajectoryStateClosestToBeamLine trajectoryStateClosestToBeamLine;
    // thread-safe helpers to guarantee proper update of mutable member data
    mutable std::atomic<char> m_TSOS;
    mutable std::atomic<char> m_TSCP;
    mutable std::atomic<char> m_SCTBL;
    //MTDsoaElements& soavalue;
   // auto mtdTimeElements_;
    //std::tuple<MTDsoaElements&, int32_t> mtdTimeElements_;
    TSCPBuilderNoMaterial builder;
    edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    reco::BeamSpot theBeamSpot;
    // to be used to setup thread states of class mutables
    enum CacheStates { kUnset, kSetting, kSet };
  };

}  // namespace reco

#endif
