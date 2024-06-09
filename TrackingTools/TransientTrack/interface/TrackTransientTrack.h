#ifndef TrackingTools_TransientTrack_TrackingTransientTrack_h
#define TrackingTools_TransientTrack_TrackingTransientTrack_h

#include <atomic>

/**
   * Concrete implementation of the TransientTrack for a reco::Track
   */

#include "TrackingTools/TransientTrack/interface/BasicTransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"

namespace reco {

  class TrackTransientTrack : public Track, public BasicTransientTrack {
  public:
    // constructor from persistent track
    TrackTransientTrack();
    TrackTransientTrack(const Track& tk, const MagneticField* field);
    TrackTransientTrack(const Track& tk, const double time, const double dtime, const MagneticField* field);

    TrackTransientTrack(const TrackRef& tk, const MagneticField* field);
    TrackTransientTrack(const TrackRef& tk, const double time, const double dtime, const MagneticField* field);

    TrackTransientTrack(const TrackRef& tk,
                        const MagneticField* field,
                        const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);
    TrackTransientTrack(const TrackRef& tk,
                        const double time,
                        const double dtime,
                        const MagneticField* field,
                        const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);
    TrackTransientTrack(const TrackRef& tk,
		        const float trkAssoc,
                        const double time,
                        const double dtime,
                       /* float mva, float pathlength, float btlchi2,  float btltimechi2, float etlchi2, float etltimechi2,float time_pi, float time_k, float time_p, float sigma_time_pi, float sigma_time_k, float sigma_time_p,*/
                        const float mva, const float pathlength, const float btlchi2, const float btltimechi2, const float etlchi2, const float etltimechi2, const float time_pi, const float time_k, const float time_p, const float sigma_time_pi, const float sigma_time_k, const float sigma_time_p,
                        const MagneticField* field,
                        const edm::ESHandle<GlobalTrackingGeometry>& tg);
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

    const Track& track() const override { return *this; }

    TrajectoryStateClosestToBeamLine stateAtBeamLine() const override;

    double timeExt() const override { return (hasTime ? timeExt_ : std::numeric_limits<double>::quiet_NaN()); }
    double dtErrorExt() const override { return (hasTime ? dtErrorExt_ : std::numeric_limits<double>::quiet_NaN()); }
    float trackAsocMTD() const { return trkAssoc_; }
    float MVAquality() const { return mva_; }
    float pathLength() const { return pathlength_; }
    float btlMatch_chi2() const { return btlchi2_; }
    float btlMatchTime_chi2() const { return btltimechi2_; }
    float etlMatch_chi2() const { return etlchi2_; }
    float etlMatchTime_chi2() const { return etltimechi2_; }
    float trackTime_pi() const { return time_pi_; }
    float trackTime_k() const { return time_k_; }
    float trackTime_p() const { return time_p_; }
    float sigma_time_pi() const {return sigma_time_pi_; }
    float sigma_time_k() const {return sigma_time_k_; }
    float sigma_time_p() const {return sigma_time_p_; }
  private:
    TrackRef tkr_;
    bool hasTime;
    float trkAssoc_;
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

    TSCPBuilderNoMaterial builder;
    edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    float mva_;
    float pathlength_;
    float btlchi2_, btltimechi2_, etlchi2_, etltimechi2_,time_pi_,time_k_,time_p_,sigma_time_pi_, sigma_time_k_, sigma_time_p_;
    reco::BeamSpot theBeamSpot;

    // to be used to setup thread states of class mutables
    enum CacheStates { kUnset, kSetting, kSet };
  };

}  // namespace reco

#endif
