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
                                                    const edm::Handle<MtdtimeHostCollection>& soa,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    TrackRef ref(trkColl, i);
    float trackAsocMTD = (*soa).const_view().trackAsocMTD()[i];
    double time = (*soa).const_view().time()[i];
    double timeReso = (*soa).const_view().timeErr()[i];
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
    /*
    const float trackAsocMTD = (*soa).const_view().trackAsocMTD()[i]; 
    const double time = (*soa).const_view().time()[i];
    const double timeReso = (*soa).const_view().timeErr()[i];
    const float MVAquality = (*soa).const_view().MVAquality()[i];
    const float pathLength = (*soa).const_view().pathLength()[i];
    const float btlMatch_chi2 = (*soa).const_view().btlMatch_chi2()[i];
    const float btlMatchTime_chi2 = (*soa).const_view().btlMatchTime_chi2()[i];
    const float etlMatch_chi2 = (*soa).const_view().etlMatch_chi2()[i];
    const float etlMatchTime_chi2 = (*soa).const_view().etlMatchTime_chi2()[i];
    const float trackTime_pi = (*soa).const_view().trackTime_pi()[i];
    const float trackTime_k = (*soa).const_view().trackTime_k()[i];
    const float trackTime_p = (*soa).const_view().trackTime_p()[i];
    const float track_sigmaTime_pi = (*soa).const_view().track_sigmaTime_pi()[i];
    const float track_sigmaTime_k = (*soa).const_view().track_sigmaTime_k()[i];
    const float track_sigmaTime_p = (*soa).const_view().track_sigmaTime_p()[i];*/
 /*   timeReso = (timeReso > 1e-6 ? timeReso
                                : defaultInvalidTrackTimeReso);  // make the error much larger than the BS time width
    if (edm::isNotFinite(time)) {
      time = 0.0;
      timeReso = defaultInvalidTrackTimeReso;
    }*/
    //ttVect.push_back(TransientTrack(ref, time, timeReso, theField,  theTrackingGeometry));
    ttVect.push_back(TransientTrack(ref, trackAsocMTD, time, timeReso, MVAquality, pathLength, btlMatch_chi2, btlMatchTime_chi2, etlMatch_chi2, etlMatchTime_chi2, trackTime_pi, trackTime_k, trackTime_p, track_sigmaTime_pi, track_sigmaTime_k, track_sigmaTime_p, theField,  theTrackingGeometry));
  }
  return ttVect;
}
