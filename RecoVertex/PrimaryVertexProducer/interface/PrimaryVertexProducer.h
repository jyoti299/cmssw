// -*- C++ -*-
//
// Package:    PrimaryVertexProducer
// Class:      PrimaryVertexProducer
//
/**\class PrimaryVertexProducer PrimaryVertexProducer.cc RecoVertex/PrimaryVertexProducer/src/PrimaryVertexProducer.cc

 Description: steers tracker primary vertex reconstruction and storage

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Pascal Vanlaer
//         Created:  Tue Feb 28 11:06:34 CET 2006
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFindingBase.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ_vect.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZT_vect.h"

#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"
#include "RecoVertex/PrimaryVertexProducer/interface/HITrackFilterForPVFinding.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GapClusterizerInZ.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexFitterBase.h"
#include "RecoVertex/PrimaryVertexProducer/interface/SequentialPrimaryVertexFitterAdapter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/AdaptiveChisquarePrimaryVertexFitter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/WeightedMeanFitter.h"

#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include <algorithm>
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoVertex/VertexTools/interface/VertexCompatibleWithBeam.h"
#include "DataFormats/Common/interface/ValueMap.h"
// vertex timing
#include "RecoVertex/PrimaryVertexProducer/interface/VertexTimeAlgorithmBase.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexTimeAlgorithmFromTracksPID.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexTimeAlgorithmLegacy4D.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "TracksGraph.h"
#include "DataFormats/VertexReco/interface/MtdtimeHostCollection.h"
#include <iostream>
#include <vector>
//
// class declaration
//
using namespace cms::Ort;

class PrimaryVertexProducer : public edm::stream::EDProducer<> {
public:
  PrimaryVertexProducer(const edm::ParameterSet&, cms::Ort::ONNXRuntime const* onnxRuntime = nullptr);
  ~PrimaryVertexProducer() override;

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  std::unique_ptr<TrackGraph>  produce_tracks_graph(const std::vector<reco::TransientTrack>& transientTrack);
  // access to config
  edm::ParameterSet config() const { return theConfig; }

private:
  // ----------member data ---------------------------
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;

  TrackFilterForPVFindingBase* theTrackFilter;
  TrackClusterizerInZ* theTrackClusterizer;

  // vtx fitting algorithms
  struct algo {
    PrimaryVertexFitterBase* pv_fitter;
    VertexCompatibleWithBeam* vertexSelector;
    std::string label;
    bool useBeamConstraint;
    double minNdof;
    VertexTimeAlgorithmBase* pv_time_estimator;
  };

  std::vector<algo> algorithms;

  edm::ParameterSet theConfig;
  bool fVerbose;

  bool fRecoveryIteration;
  edm::EDGetTokenT<reco::VertexCollection> recoveryVtxToken;

  edm::EDGetTokenT<reco::BeamSpot> bsToken;
  edm::EDGetTokenT<reco::TrackCollection> trkToken;
  edm::EDGetTokenT<MtdtimeHostCollection> inputTimingToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > trkTimesToken;
  edm::EDGetTokenT<edm::ValueMap<float> > trkTimeResosToken;
  edm::EDGetTokenT<edm::ValueMap<float> > trackMTDTimeQualityToken;

  bool useTransientTrackTime_;
  bool useMVASelection_;
  edm::ValueMap<float> trackMTDTimeQualities_;
  edm::ValueMap<float> trackTimes_;
  double minTrackTimeQuality_;

  std::vector<float> features;
  std::vector<float> edge_features;

  std::vector<float> node_degrees;
  std::vector<float> degree_centr;
  std::string nnVersion_;       // Version identifier of the NN (either CNN or a GNN, to choose which inputs to use)
  double nnWorkingPoint_;       // Working point for neural network (above this score, consider the t
  cms::Ort::ONNXRuntime const* onnxRuntime_;
};
