#ifndef VertexProducerbyGNN_VertexProducerbyGNN_h
#define VertexProducerbyGNN_VertexProducerbyGNN_h

/**_________________________________________________________________
   class:   VertexProducerbyGNN.h

________________________________________________________________**/

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "CondFormats/DataRecord/interface/BeamSpotObjectsRcd.h"
#include "CondFormats/BeamSpotObjects/interface/BeamSpotObjects.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
class VertexProducerbyGNN : public edm::stream::EDProducer<> {
public:
  typedef std::vector<edm::ParameterSet> Parameters;

  /// constructor
  explicit VertexProducerbyGNN(const edm::ParameterSet& conf, cms::Ort::ONNXRuntime const* onnxRuntime = nullptr);
  /// destructor
  ~VertexProducerbyGNN() override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet &);
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  std::unique_ptr<TracksGraph> produce_tracks_graph(const std::vector<reco::TransientTrack>& transientTrack);

private:
  edm::ESGetToken<BeamSpotObjects, BeamSpotObjectsRcd> m_beamToken;

  TrackFilterForPVFindingBase* theTrackFilter;
  TrackClusterizerInZ* theTrackClusterizer;
  edm::EDGetTokenT<reco::BeamSpot> bsToken;
  edm::EDGetTokenT<reco::TrackCollection> trkToken;
  edm::EDGetTokenT<edm::ValueMap<float> > trkTimesToken;
  edm::EDGetTokenT<edm::ValueMap<float> > trkTimeResosToken;
  edm::EDGetTokenT<edm::ValueMap<float> > trackMTDTimeQualityToken;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;


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
};

#endif
