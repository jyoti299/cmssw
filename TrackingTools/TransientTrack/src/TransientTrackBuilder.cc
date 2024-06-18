#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/Common/interface/Handle.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTS.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

using namespace reco;
using namespace std;
using namespace edm;

TransientTrack TransientTrackBuilder::build(const Track* t) const {
  return TransientTrack(*t, theField, theTrackingGeometry);
}

TransientTrack TransientTrackBuilder::build(const Track& t) const {
  return TransientTrack(t, theField, theTrackingGeometry);
}

TransientTrack TransientTrackBuilder::build(const GsfTrack* t) const {
  return TransientTrack(new GsfTransientTrack(*t, theField, theTrackingGeometry));
}

TransientTrack TransientTrackBuilder::build(const GsfTrack& t) const {
  return TransientTrack(new GsfTransientTrack(t, theField, theTrackingGeometry));
}

TransientTrack TransientTrackBuilder::build(const CandidatePtr* t) const {
  reco::PFCandidatePtr tryPF(*t);
  edm::Ptr<pat::PackedCandidate> tryPacked(*t);
  if (tryPF.get() != nullptr && tryPF->isTimeValid()) {
    return TransientTrack(*t, tryPF->time(), tryPF->timeError(), theField, theTrackingGeometry);
  } else if (tryPacked.get() != nullptr && tryPacked->timeError() > 0.f) {
    return TransientTrack(*t, (double)tryPacked->time(), (double)tryPacked->timeError(), theField, theTrackingGeometry);
  }
  return TransientTrack(*t, theField, theTrackingGeometry);
}

TransientTrack TransientTrackBuilder::build(const CandidatePtr& t) const { return this->build(&t); }

TransientTrack TransientTrackBuilder::build(const TrackRef* t) const {
  return TransientTrack(*t, theField, theTrackingGeometry);
}

TransientTrack TransientTrackBuilder::build(const TrackRef& t) const {
  return TransientTrack(t, theField, theTrackingGeometry);
}

TransientTrack TransientTrackBuilder::build(const GsfTrackRef* t) const {
  return TransientTrack(new GsfTransientTrack(*t, theField, theTrackingGeometry));
}

TransientTrack TransientTrackBuilder::build(const GsfTrackRef& t) const {
  return TransientTrack(new GsfTransientTrack(t, theField, theTrackingGeometry));
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::TrackCollection>& trkColl) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    ttVect.push_back(TransientTrack(TrackRef(trkColl, i), theField, theTrackingGeometry));
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::GsfTrackCollection>& trkColl) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    ttVect.push_back(TransientTrack(new GsfTransientTrack(GsfTrackRef(trkColl, i), theField, theTrackingGeometry)));
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<edm::View<Track> >& trkColl) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    const Track* trk = &(*trkColl)[i];
    const GsfTrack* gsfTrack = dynamic_cast<const GsfTrack*>(trk);
    if (gsfTrack) {
      ttVect.push_back(TransientTrack(
          new GsfTransientTrack(RefToBase<Track>(trkColl, i).castTo<GsfTrackRef>(), theField, theTrackingGeometry)));
    } else {  // gsf
      ttVect.push_back(TransientTrack(RefToBase<Track>(trkColl, i).castTo<TrackRef>(), theField, theTrackingGeometry));
    }
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::TrackCollection>& trkColl,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    TrackRef ref(trkColl, i);
    double time = trackTimes[ref];
    double timeReso = trackTimeResos[ref];
    timeReso = (timeReso > 1e-6 ? timeReso
                                : defaultInvalidTrackTimeReso);  // make the error much larger than the BS time width
    if (edm::isNotFinite(time)) {
      time = 0.0;
      timeReso = defaultInvalidTrackTimeReso;
    }
    ttVect.push_back(TransientTrack(ref, time, timeReso, theField, theTrackingGeometry));
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::GsfTrackCollection>& trkColl,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    GsfTrackRef ref(trkColl, i);
    double time = trackTimes[ref];
    double timeReso = trackTimeResos[ref];
    timeReso = (timeReso > 1e-6 ? timeReso
                                : defaultInvalidTrackTimeReso);  // make the error much larger than the BS time width
    if (edm::isNotFinite(time)) {
      time = 0.0;
      timeReso = defaultInvalidTrackTimeReso;
    }
    ttVect.push_back(TransientTrack(new GsfTransientTrack(ref, time, timeReso, theField, theTrackingGeometry)));
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<edm::View<Track> >& trkColl,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    const Track* trk = &(*trkColl)[i];
    const GsfTrack* gsfTrack = dynamic_cast<const GsfTrack*>(trk);
    if (gsfTrack) {
      GsfTrackRef ref = RefToBase<Track>(trkColl, i).castTo<GsfTrackRef>();
      double time = trackTimes[ref];
      double timeReso = trackTimeResos[ref];
      timeReso = (timeReso > 1e-6 ? timeReso
                                  : defaultInvalidTrackTimeReso);  // make the error much larger than the BS time width
      if (edm::isNotFinite(time)) {
        time = 0.0;
        timeReso = defaultInvalidTrackTimeReso;
      }
      ttVect.push_back(TransientTrack(new GsfTransientTrack(
          RefToBase<Track>(trkColl, i).castTo<GsfTrackRef>(), time, timeReso, theField, theTrackingGeometry)));
    } else {  // gsf
      TrackRef ref = RefToBase<Track>(trkColl, i).castTo<TrackRef>();
      double time = trackTimes[ref];
      double timeReso = trackTimeResos[ref];
      timeReso = (timeReso > 1e-6 ? timeReso
                                  : defaultInvalidTrackTimeReso);  // make the error much larger than the BS time width
      if (edm::isNotFinite(time)) {
        time = 0.0;
        timeReso = defaultInvalidTrackTimeReso;
      }
      ttVect.push_back(TransientTrack(
          RefToBase<Track>(trkColl, i).castTo<TrackRef>(), time, timeReso, theField, theTrackingGeometry));
    }
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::TrackCollection>& trkColl,
                                                    const reco::BeamSpot& beamSpot) const {
  vector<TransientTrack> ttVect = build(trkColl);
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::GsfTrackCollection>& trkColl,
                                                    const reco::BeamSpot& beamSpot) const {
  vector<TransientTrack> ttVect = build(trkColl);
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<edm::View<Track> >& trkColl,
                                                    const reco::BeamSpot& beamSpot) const {
  vector<TransientTrack> ttVect = build(trkColl);
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::TrackCollection>& trkColl,
                                                    const reco::BeamSpot& beamSpot,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect = build(trkColl, trackTimes, trackTimeResos);
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::GsfTrackCollection>& trkColl,
                                                    const reco::BeamSpot& beamSpot,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect = build(trkColl, trackTimes, trackTimeResos);
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<edm::View<Track> >& trkColl,
                                                    const reco::BeamSpot& beamSpot,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect = build(trkColl, trackTimes, trackTimeResos);
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}

TransientTrack TransientTrackBuilder::build(const FreeTrajectoryState& fts) const {
  return TransientTrack(new TransientTrackFromFTS(fts));
}


vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::TrackCollection>& trkColl,
						    const reco::BeamSpot& beamSpot,
                                                    const edm::ValueMap<float>& MTDtimes,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    TrackRef ref(trkColl, i);
    double time = trackTimes[ref];
    double timeReso = trackTimeResos[ref];
    int trackAsoc_mtd = 1;
    float MTDtime_mtd = 2.0;
    float MTDtimeErr_mtd = 3.0;
    float MVAquality_mtd = 3.0;
    float pathLength_mtd = 3.0;
    float btlMatch_chi2_mtd = 3.0;
    float btlMatchTime_chi2_mtd = 3.0;
    float etlMatch_chi2_mtd = 3.0;
    float etlMatchTime_chi2_mtd = 3.0;
    float trackTime_pi_mtd = 3.0;
    float trackTime_k_mtd = 3.0;
    float trackTime_p_mtd = 3.0;
    float track_sigmaTime_pi_mtd = 3.0;
    float track_sigmaTime_k_mtd = 3.0;
    float track_sigmaTime_p_mtd = 3.0;

/*    int trackAsocMTD = (*soa).const_view().trackAsocMTD()[i];
    float MTDtime = (*soa).const_view().time()[i];
    float MTDtimeErr = (*soa).const_view().timeErr()[i];
    float MVAquality = (*soa).const_view().MVAquality()[i];
    float pathLength = (*soa).const_view().pathLength()[i];
    float btlMatch_chi2 = (*soa).const_view().btlMatch_chi2()[i];
    float btlMatchTime_chi2 = (*soa).const_view().btlMatchTime_chi2()[i];
    float etlMatch_chi2 = (*soa).const_view().etlMatch_chi2()[i];
    float etlMatchTime_chi2 = (*soa).const_view().etlMatchTime_chi2()[i];
    float trackTime_pi = (*soa).const_view().trackTime_pi()[i];
    float trackTime_k = (*soa).const_view().trackTime_k()[i];
    float trackTime_p = (*soa).const_view().trackTime_p()[i];
    float track_sigmaTime_pi = (*soa).const_view().track_sigmaTime_pi()[i];
    float track_sigmaTime_k = (*soa).const_view().track_sigmaTime_k()[i];
    float track_sigmaTime_p = (*soa).const_view().track_sigmaTime_p()[i];
  */  
    timeReso = (timeReso > 1e-6 ? timeReso
                                : defaultInvalidTrackTimeReso);  // make the error much larger than the BS time width
    if (edm::isNotFinite(time)) {
      time = 0.0;
      timeReso = defaultInvalidTrackTimeReso;
    }
    ttVect.push_back(TransientTrack(ref, time, timeReso, theField,  theTrackingGeometry, trackAsoc_mtd, MTDtime_mtd, MTDtimeErr_mtd, MVAquality_mtd, pathLength_mtd, btlMatch_chi2_mtd, btlMatchTime_chi2_mtd, etlMatch_chi2_mtd, etlMatchTime_chi2_mtd, trackTime_pi_mtd, trackTime_k_mtd, trackTime_p_mtd, track_sigmaTime_pi_mtd, track_sigmaTime_k_mtd, track_sigmaTime_p_mtd));
  }
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}
