
#include "RecoVertex/PrimaryVertexProducer/interface/VertexProducerbyGNN.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/IOVSyncValue.h"

#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"
//

enum class DNNVersion {
    GNN_V1,
    MLP_EDGE_FEATURES,
    MLP_NO_EDGE_FEATURES
};

// Define constants for DNN versions
std::unordered_map<std::string, DNNVersion> DNN_VERSION_MAP = {
    {"gnn_v1", DNNVersion::GNN_V1},
    {"mlp_edge_features", DNNVersion::MLP_EDGE_FEATURES},
    {"mlp_no_edge_features", DNNVersion::MLP_NO_EDGE_FEATURES}
};

// Define a map to store DNN input configurations
std::unordered_map<DNNVersion, std::pair<std::vector<std::string>, int>> dnnInputConfigurations = {
    {DNNVersion::GNN_V1, {{"features", "edge_index", "edge_features"}, FEATURE_SHAPE_GNN_V1}},
    {DNNVersion::MLP_EDGE_FEATURES, {{"features", "edge_index", "edge_features"}, FEATURE_SHAPE_MLP}},
    {DNNVersion::MLP_NO_EDGE_FEATURES, {{"features", "edge_index"}, FEATURE_SHAPE_MLP}}
};


VertexProducerbyGNN::VertexProducerbyGNN(const edm::ParameterSet& conf, cms::Ort::ONNXRuntime const* onnxRuntime)
  : theTTBToken(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))), theConfig(conf) {
  trkToken = consumes<reco::TrackCollection>(conf.getParameter<edm::InputTag>("TrackLabel"));
  bsToken = consumes<reco::BeamSpot>(conf.getParameter<edm::InputTag>("beamSpotLabel"));
  inputTimingToken_ = consumes<MtdtimeHostCollection>(conf.getParameter<edm::InputTag>("timingSoA"));
  trackBuilderToken_ = esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))} {
  useTransientTrackTime_ = false;

  // select and configure the track selection
  //fVerbose=conf.getUntrackedParameter<bool>("verbose", false);
  m_beamToken = esConsumes<BeamSpotObjects, BeamSpotObjectsRcd>();

  produces<reco::BeamSpot>();
}

VertexProducerbyGNN::~VertexProducerbyGNN() {}


/// ******* work on this function to produce a graph or in produce method. 
std::unique_ptr<TracksGraph> VertexProducerbyGNN::produce_tracks_graph(const std::vector<reco::TransientTrack>& transientTracks) {
	// need work on this, have to take reference from this https://github.com/cms-sw/cmssw/blob/34d92a71fb188acad9af249df6a9bc48436432e3/DataFormats/Math/interface/Graph.h
    auto tracksGraph = std::make_unique<TracksGraph>();
    for (const auto& track : transientTracks) {
        tracksGraph->addNode(track);
    }
}

/// *******************
/// ------------ method called to produce the data  ------------
void VertexProducerbyGNN::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
     using namespace edm;

     DNNVersion versionEnum = DNN_VERSION_MAP[nnVersion_];

    // Access the input configuration based on the DNN version
    std::vector<std::string> input_names;
    int shapeFeatures;

    auto dnnConfig = dnnInputConfigurations.find(versionEnum);
    if (dnnConfig != dnnInputConfigurations.end()) {
        // Valid DNN version found, retrieve input names and shape features
        input_names = dnnConfig->second.first;
        shapeFeatures = dnnConfig->second.second;
    } else {
        // Invalid DNN version
        edm::LogError("Vertex Producer") << "Architecture not defined: " << nnVersion_ ;
        return;
    }

    std::cout << "Input names:";
    for (const auto& name : input_names) {
        std::cout << " " << name;
    }
    std::cout << std::endl;
    std::cout << "Shape of features: " << shapeFeatures << std::endl;

    // Array of data to be filled as a network input. Should be a float array of flattened values.
    cms::Ort::FloatArrays data;
    // Network input shapes.
    std::vector<std::vector<int64_t>> input_shapes;


  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(bsToken, recoBeamSpotHandle);
  if (recoBeamSpotHandle.isValid()) {
    beamSpot = *recoBeamSpotHandle;
  } else {
    edm::LogError("UnusableBeamSpot") << "No beam spot available from EventSetup";
  }

  edm::Handle<reco::TrackCollection> t_tks;
  iEvent.getByToken(trkToken, tks);
  const auto& theB = &iSetup.getData(theTTBToken);
  std::vector<reco::TransientTrack> tks;
 
  edm::Handle<MtdtimeHostCollection> inputTiming_h;
  auto const &inputTimingView = iEvent.get(inputTimingToken_).const_view();
  auto const maxsize = inputTimingView.metadata().size();

  tks = (*theB).build(t_tks, inputTiming_h, trackTimes_, trackTimeResos_);

   std::vector<float> edges_src;
   std::vector<float> edges_dst;
  // edge features are computed from node features only

 
    for (size_t i = 0; i < tks.size() - 1; ++i) {
        const reco::TransientTrack& ttrack1 = tks[i];
        const reco::TransientTrack& ttrack2 = tks[i + 1];

        float trackMTD = ttrack1.trackAsocMTD();
        float mtdtime = ttrack1.time();
        float mtdtimeErr = ttrack1.timeErr();		
        float mvaquality = ttrack1.MVAquality();
	float pathlength = ttrack1.pathLength();       	
        float btlMatchchi2 = ttrack1.btlMatch_chi2();
	float btlMatchTimechi2 = ttrack1.btlMatchTime_chi2();
        float etlMatchchi2 = ttrack1.etlMatch_chi2();
        float etlMatchTimechi2 = ttrack1.etlMatchTime_chi2();
	float Time_pi_hypo = ttrack1.trackTime_pi(); 
        float Time_k_hypo = ttrack1.trackTime_k();
        float Time_p_hypo = ttrack1.trackTime_p();
        float sigmaTime_pi = ttrack1.track_sigmaTime_pi(); 
	float sigmaTime_k = ttrack1.track_sigmaTime_k();
        float sigmaTime_p = ttrack1.track_sigmaTime_p();

        features.push_back(trackMTD);
	features.push_back(mtdtime);
        features.push_back(mtdtimeErr);
        features.push_back(mvaquality);
        features.push_back(pathlength);
        features.push_back(btlMatchchi2);
        features.push_back(btlMatchTimechi2);		
        features.push_back(etlMatchchi2);
        features.push_back(etlMatchTimechi2);
 	features.push_back(Time_pi_hypo);
	features.push_back(Time_k_hypo);
        features.push_back(Time_p_hypo);
	features.push_back(sigmaTime_pi);
	features.push_back(sigmaTime_k);
	features.push_back(sigmaTime_p);

        double pt_diff = std::abs(ttrack1.pt() - ttrack2.pt());
        double zpca_diff = std::abs(ttrack1.vz() - ttrack2.vz());
        double dz_sign = std::abs(ttrack1.vz() - ttrack2.vz())/sqrt(((ttrack1.dzError[i])*(ttrack1.dzError[i]))+((ttrack2.dzError[i])*(ttrack2.dzError[i])));
        float time_pi_diff = std::abs(ttrack1.trackTime_pi() -  ttrack2.trackTime_pi()); 
        float time_k_diff = std::abs(ttrack1.trackTime_k() -  ttrack2.trackTime_k());
        float time_p_diff = std::abs(ttrack1.trackTime_p() -  ttrack2.trackTime_p());
        float time_pi_k_diff = std::abs(ttrack1.trackTime_pi() - ttrack2.trackTime_k());
        float time_pi_p_diff = std::abs(ttrack1.trackTime_pi() - ttrack2.trackTime_p());
        float time_k_pi_diff = std::abs(ttrack1.trackTime_k() - ttrack2.trackTime_pi());
        float time_k_p_diff = std::abs(ttrack1.trackTime_k() - ttrack2.trackTime_p());
        float time_p_pi_diff = std::abs(ttrack1.trackTime_p() - ttrack2.trackTime_pi());
        float time_p_k_diff = std::abs(ttrack1.trackTime_p() - ttrack2.trackTime_k());

	
	edges_src.push_back(ttrack2); // this should be neighbour of the track to create an edge between the tracks.. I hope this should be filled into it. 
        edges_dst.push_back(ttrack1);


    edge_features.push_back(pt_diff);
    edge_features.push_back(zpca_diff);
    edge_features.push_back(dz_sign);		    
    edge_features.push_back(time_pi_diff);
    edge_features.push_back(time_k_diff);
    edge_features.push_back(time_p_diff);
    edge_features.push_back(time_pi_k_diff);
    edge_features.push_back(time_pi_p_diff);
    edge_features.push_back(time_k_pi_diff);
    edge_features.push_back(time_k_pi_diff);
    edge_features.push_back(time_k_p_diff);
    edge_features.push_back(time_p_pi_diff);
    edge_features.push_back(time_p_k_diff);
  }
     
  // check what is with edges_src
    auto numEdges = static_cast<int>(edges_src.size());

    input_shapes.push_back({1, N, shapeFeatures});
    data.emplace_back(features);

    input_shapes.push_back({1, 2, numEdges});

    data.emplace_back(edges_src);
    for (auto &dst : edges_dst)
    {
        data.back().push_back(dst);
    }

// run prediction 
    std::vector<float> edge_predictions = onnxRuntime_->run(input_names, data, input_shapes)[0];

}

void VertexProducerbyGNN::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.add<edm::FileInPath>("onnxModelPath", edm::FileInPath("RecoVertex/PrimaryVertexProducer/data/model.onnx"))->setComment("Path to GNN (as ONNX model)");
        desc.add<edm::InputTag>("timingSoA", edm::InputTag("mtdSoA"));
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexProducerbyGNN);
