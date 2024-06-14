#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
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

const int FEATURE_SHAPE_GNN_V1 = 19;
const int FEATURE_SHAPE_MLP = 30;
const int NUM_EDGE_FEATURES = 12;

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

PrimaryVertexProducer::PrimaryVertexProducer(const edm::ParameterSet& conf, cms::Ort::ONNXRuntime const* onnxRuntime)
    : theTTBToken(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))), theConfig(conf) {
  fVerbose = conf.getUntrackedParameter<bool>("verbose", false);
  useMVASelection_ = conf.getParameter<bool>("useMVACut");

  trkToken = consumes<reco::TrackCollection>(conf.getParameter<edm::InputTag>("TrackLabel"));
  bsToken = consumes<reco::BeamSpot>(conf.getParameter<edm::InputTag>("beamSpotLabel"));
  useTransientTrackTime_ = false;

  // select and configure the track selection
  std::string trackSelectionAlgorithm =
      conf.getParameter<edm::ParameterSet>("TkFilterParameters").getParameter<std::string>("algorithm");
  if (trackSelectionAlgorithm == "filter") {
    theTrackFilter = new TrackFilterForPVFinding(conf.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else if (trackSelectionAlgorithm == "filterWithThreshold") {
    theTrackFilter = new HITrackFilterForPVFinding(conf.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else {
    throw VertexException("PrimaryVertexProducer: unknown track selection algorithm: " + trackSelectionAlgorithm);
  }

  // select and configure the track clusterizer
  std::string clusteringAlgorithm =
      conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<std::string>("algorithm");
  if (clusteringAlgorithm == "gap") {
    theTrackClusterizer = new GapClusterizerInZ(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkGapClusParameters"));
  } else if (clusteringAlgorithm == "DA_vect") {
    theTrackClusterizer = new DAClusterizerInZ_vect(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
  } else if (clusteringAlgorithm == "DA2D_vect") {
    theTrackClusterizer = new DAClusterizerInZT_vect(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
    useTransientTrackTime_ = true;
  } else {
    throw VertexException("PrimaryVertexProducer: unknown clustering algorithm: " + clusteringAlgorithm);
  }

  if (useTransientTrackTime_) {
    trkTimesToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("TrackTimesLabel"));
    trkTimeResosToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("TrackTimeResosLabel"));
    trackMTDTimeQualityToken =
        consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("trackMTDTimeQualityVMapTag"));
    minTrackTimeQuality_ = conf.getParameter<double>("minTrackTimeQuality");
    inputTimingToken_ = consumes<MtdtimeHostCollection>(conf.getParameter<edm::InputTag>("timingSoA"));
  }

  // select and configure the vertex fitters
  std::vector<edm::ParameterSet> vertexCollections =
      conf.getParameter<std::vector<edm::ParameterSet> >("vertexCollections");

  for (std::vector<edm::ParameterSet>::const_iterator algoconf = vertexCollections.begin();
       algoconf != vertexCollections.end();
       algoconf++) {
    algo algorithm;

    algorithm.label = algoconf->getParameter<std::string>("label");

    // configure the fitter and selector
    std::string fitterAlgorithm = algoconf->getParameter<std::string>("algorithm");
    if (fitterAlgorithm == "KalmanVertexFitter") {
      algorithm.pv_fitter = new SequentialPrimaryVertexFitterAdapter(new KalmanVertexFitter());
    } else if (fitterAlgorithm == "AdaptiveVertexFitter") {
      auto fitter = new AdaptiveVertexFitter(GeometricAnnealing(algoconf->getParameter<double>("chi2cutoff")));
      algorithm.pv_fitter = new SequentialPrimaryVertexFitterAdapter(fitter);
    } else if (fitterAlgorithm.empty()) {
      algorithm.pv_fitter = nullptr;
    } else if (fitterAlgorithm == "AdaptiveChisquareVertexFitter") {
      algorithm.pv_fitter = new AdaptiveChisquarePrimaryVertexFitter(algoconf->getParameter<double>("chi2cutoff"),
                                                                     algoconf->getParameter<double>("zcutoff"),
                                                                     algoconf->getParameter<double>("mintrkweight"),
                                                                     false);
    } else if (fitterAlgorithm == "MultiPrimaryVertexFitter") {
      algorithm.pv_fitter = new AdaptiveChisquarePrimaryVertexFitter(algoconf->getParameter<double>("chi2cutoff"),
                                                                     algoconf->getParameter<double>("zcutoff"),
                                                                     algoconf->getParameter<double>("mintrkweight"),
                                                                     true);
    } else if (fitterAlgorithm == "WeightedMeanFitter") {
      algorithm.pv_fitter = new WeightedMeanPrimaryVertexEstimator();
    } else {
      throw VertexException("PrimaryVertexProducer: unknown algorithm: " + fitterAlgorithm);
    }
    algorithm.minNdof = algoconf->getParameter<double>("minNdof");
    algorithm.useBeamConstraint = algoconf->getParameter<bool>("useBeamConstraint");
    algorithm.vertexSelector =
        new VertexCompatibleWithBeam(VertexDistanceXY(), algoconf->getParameter<double>("maxDistanceToBeam"));

    // configure separate vertex time reconstruction if applicable
    // note that the vertex time could, in principle, also come from the clusterizer or the vertex fit

    const auto& pv_time_conf = algoconf->getParameter<edm::ParameterSet>("vertexTimeParameters");
    const std::string vertexTimeAlgorithm = pv_time_conf.getParameter<std::string>("algorithm");
    edm::ConsumesCollector&& collector = consumesCollector();

    if (vertexTimeAlgorithm.empty()) {
      algorithm.pv_time_estimator = nullptr;
    } else if (vertexTimeAlgorithm == "legacy4D") {
      useTransientTrackTime_ = true;
      algorithm.pv_time_estimator =
          new VertexTimeAlgorithmLegacy4D(pv_time_conf.getParameter<edm::ParameterSet>("legacy4D"), collector);
    } else if (vertexTimeAlgorithm == "fromTracksPID") {
      algorithm.pv_time_estimator = new VertexTimeAlgorithmFromTracksPID(
          pv_time_conf.getParameter<edm::ParameterSet>("fromTracksPID"), collector);
    } else {
      edm::LogWarning("MisConfiguration") << "unknown vertexTimeParameters.algorithm" << vertexTimeAlgorithm;
    }
    algorithms.push_back(algorithm);

    produces<reco::VertexCollection>(algorithm.label);
  }

  //check if this is a recovery iteration
  fRecoveryIteration = conf.getParameter<bool>("isRecoveryIteration");
  if (fRecoveryIteration) {
    if (algorithms.empty()) {
      throw VertexException("PrimaryVertexProducer: No algorithm specified. ");
    } else if (algorithms.size() > 1) {
      throw VertexException(
          "PrimaryVertexProducer: Running in Recovery mode and more than one algorithm specified.  Please "
          "only one algorithm.");
    }
    recoveryVtxToken = consumes<reco::VertexCollection>(conf.getParameter<edm::InputTag>("recoveryVtxCollection"));
  }
}

PrimaryVertexProducer::~PrimaryVertexProducer() {
  if (theTrackFilter)
    delete theTrackFilter;
  if (theTrackClusterizer)
    delete theTrackClusterizer;
  for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
    if (algorithm->pv_fitter)
      delete algorithm->pv_fitter;
    if (algorithm->pv_time_estimator)
      delete algorithm->pv_time_estimator;
    if (algorithm->vertexSelector)
      delete algorithm->vertexSelector;
  }
}
// can shift this function to TracksGraph.h as well
std::unique_ptr<TrackGraph>  PrimaryVertexProducer::produce_tracks_graph(const std::vector<reco::TransientTrack>& transientTracks) {
    std::vector<Node> allNodes;
    for (size_t i = 0; i < transientTracks.size(); ++i) {
              Node tNode(i);
          allNodes.push_back(tNode);
    }
     auto resultGraph = std::make_unique<TrackGraph>(allNodes);

    return resultGraph;
}
void PrimaryVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // get the BeamSpot, it will always be needed, even when not used as a constraint
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(bsToken, recoBeamSpotHandle);
  if (recoBeamSpotHandle.isValid()) {
    beamSpot = *recoBeamSpotHandle;
  } else {
    edm::LogError("UnusableBeamSpot") << "No beam spot available from EventSetup";
  }

  bool validBS = true;
  VertexState beamVertexState(beamSpot);
  if ((beamVertexState.error().cxx() <= 0.) || (beamVertexState.error().cyy() <= 0.) ||
      (beamVertexState.error().czz() <= 0.)) {
    edm::LogError("UnusableBeamSpot") << "Beamspot with invalid errors " << beamVertexState.error().matrix();
    validBS = false;
  }

  //if this is a recovery iteration, check if we already have a valid PV
  if (fRecoveryIteration) {
    auto const& oldVertices = iEvent.get(recoveryVtxToken);
    //look for the first valid (not-BeamSpot) vertex
    for (auto const& old : oldVertices) {
      if (!(old.isFake())) {
        //found a valid vertex, write the first one to the collection and return
        //otherwise continue with regular vertexing procedure
        auto result = std::make_unique<reco::VertexCollection>();
        result->push_back(old);
        iEvent.put(std::move(result), algorithms.begin()->label);
        return;
      }
    }
  }

  // get RECO tracks from the event
  // `tks` can be used as a ptr to a reco::TrackCollection
  edm::Handle<reco::TrackCollection> tks;
  iEvent.getByToken(trkToken, tks);

  // mechanism to put the beamspot if the track collection is empty
  if (!tks.isValid()) {
    for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
      auto result = std::make_unique<reco::VertexCollection>();
      reco::VertexCollection& vColl = (*result);

      GlobalError bse(beamSpot.rotatedCovariance3D());
      if ((bse.cxx() <= 0.) || (bse.cyy() <= 0.) || (bse.czz() <= 0.)) {
        AlgebraicSymMatrix33 we;
        we(0, 0) = 10000;
        we(1, 1) = 10000;
        we(2, 2) = 10000;
        vColl.push_back(reco::Vertex(beamSpot.position(), we, 0., 0., 0));
        if (fVerbose) {
          std::cout << "RecoVertex/PrimaryVertexProducer: "
                    << "Beamspot with invalid errors " << bse.matrix() << std::endl;
          std::cout << "Will put Vertex derived from dummy-fake BeamSpot into Event.\n";
        }
      } else {
        vColl.push_back(reco::Vertex(beamSpot.position(), beamSpot.rotatedCovariance3D(), 0., 0., 0));
        if (fVerbose) {
          std::cout << "RecoVertex/PrimaryVertexProducer: "
                    << " will put Vertex derived from BeamSpot into Event.\n";
        }
      }
      iEvent.put(std::move(result), algorithm->label);
    }

    return;  // early return
  }

  for (auto& algo : algorithms) {
    if (algo.pv_time_estimator) {
      algo.pv_time_estimator->setEvent(iEvent, iSetup);
    }
  }
   // for GNN version
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
    std::vector<float> edges_src;
    std::vector<float> edges_dst;

  // interface RECO tracks to vertex reconstruction
  const auto& theB = &iSetup.getData(theTTBToken);
  std::vector<reco::TransientTrack> t_tks;
  edm::Handle<MtdtimeHostCollection> inputTiming_h;


  if (useTransientTrackTime_) {
    auto const& trackTimeResos_ = iEvent.get(trkTimeResosToken);
    auto trackTimes_ = iEvent.get(trkTimesToken);
    iEvent.getByToken(inputTimingToken_, inputTiming_h);
    auto const &inputTimingView = (*inputTiming_h).const_view();
    std::cout <<" if we are in here "<<std::endl;
    if (useMVASelection_) {
      trackMTDTimeQualities_ = iEvent.get(trackMTDTimeQualityToken);

      for (unsigned int i = 0; i < (*tks).size(); i++) {
        const reco::TrackRef ref(tks, i);
        auto const trkTimeQuality = trackMTDTimeQualities_[ref];
        if (trkTimeQuality < minTrackTimeQuality_) {
          trackTimes_[ref] = std::numeric_limits<double>::max();
        }
      }
      t_tks = (*theB).build(tks, beamSpot, trackTimes_, trackTimeResos_);
    } else {
      t_tks = (*theB).build(tks, beamSpot, trackTimes_, trackTimeResos_);
     // t_tks = (*theB).build(tks, inputTiming_h, trackTimes_, trackTimeResos_);
    }
  } else {
    t_tks = (*theB).build(tks, beamSpot);
  }

  if (useTransientTrackTime_) {
    long int N = t_tks.size();
    std::cout<<" How many N "<<N<<std::endl;
   for (size_t i = 0; i < t_tks.size() - 1; ++i) {
        const reco::TransientTrack& ttrack1 = t_tks[i];
        float trackMTD = ttrack1.trackAsocMTD();
        float mtdtime = ttrack1.timeExt();
        float mtdtimeErr = ttrack1.dtErrorExt();
        float mvaquality = ttrack1.MVAquality();
        float pathlength = ttrack1.pathLength();
        float btlMatchchi2 = ttrack1.btlMatch_chi2();
        float btlMatchTimechi2 = ttrack1.btlMatchTime_chi2();
        float etlMatchchi2 = ttrack1.etlMatch_chi2();
        float etlMatchTimechi2 = ttrack1.etlMatchTime_chi2();
        float Time_pi_hypo = ttrack1.trackTime_pi();
        float Time_k_hypo = ttrack1.trackTime_k();
        float Time_p_hypo = ttrack1.trackTime_p();
        float sigmaTime_pi = ttrack1.sigma_time_pi();
        float sigmaTime_k = ttrack1.sigma_time_k();
        float sigmaTime_p = ttrack1.sigma_time_p();

       // std::cout<<"$$$$$$$$$$$$$$$ Here 1 ##########"<<mvaquality<<" or "<<Time_k_hypo<<" mva in here "<<std::endl;

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
   for (size_t j = i + 1; j < t_tks.size(); ++j) {
        const reco::TransientTrack& ttrack2 = t_tks[j];


        double neighbour_cut = std::abs(ttrack1.track().vz() - ttrack2.track().vz());
        if (neighbour_cut < 0.3 ) {
        //std::cout<<" Here 4************ "<<j<<std::endl;
        edges_src.push_back(static_cast<float>(j)); // this should be neighbour of the track to create an edge between the tracks.. I hope this should be filled into it. 
        edges_dst.push_back(static_cast<float>(i));

	double pt_diff = std::abs(ttrack1.track().pt() - ttrack2.track().pt());
        double zpca_diff = std::abs(ttrack1.track().vz() - ttrack2.track().vz());
        double dz_sign = std::abs(ttrack1.track().vz() - ttrack2.track().vz())/sqrt(((ttrack1.track().dzError())*(ttrack1.track().dzError()))+((ttrack2.track().dzError())*(ttrack2.track().dzError())));
        float time_pi_diff = std::abs(ttrack1.trackTime_pi() -  ttrack2.trackTime_pi());
        float time_k_diff = std::abs(ttrack1.trackTime_k() -  ttrack2.trackTime_k());
        float time_p_diff = std::abs(ttrack1.trackTime_p() -  ttrack2.trackTime_p());
        float time_pi_k_diff = std::abs(ttrack1.trackTime_pi() - ttrack2.trackTime_k());
        float time_pi_p_diff = std::abs(ttrack1.trackTime_pi() - ttrack2.trackTime_p());
        float time_k_pi_diff = std::abs(ttrack1.trackTime_k() - ttrack2.trackTime_pi());
        float time_k_p_diff = std::abs(ttrack1.trackTime_k() - ttrack2.trackTime_p());
        float time_p_pi_diff = std::abs(ttrack1.trackTime_p() - ttrack2.trackTime_pi());
        float time_p_k_diff = std::abs(ttrack1.trackTime_p() - ttrack2.trackTime_k());

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
	 //std::cout<<" Here 5 ##########"<<time_pi_diff<<std::endl;
    }
   }
  }
  std::cout<<" Here 6"<<std::endl; 
   auto numEdges = static_cast<int>(edges_src.size());

    input_shapes.push_back({1, N, shapeFeatures});
    data.emplace_back(features);
    input_shapes.push_back({1, 2, numEdges});
std::cout<<" Here 70000000000000000"<<std::endl;
    data.emplace_back(edges_src);
    for (auto &dst : edges_dst)
    {
            data.back().push_back(dst);
    }
    std::cout<<" Here 80000000000000000"<<std::endl;
    input_shapes.push_back({1, numEdges, NUM_EDGE_FEATURES});
    data.emplace_back(edge_features);

std::cout<<" Here 90000000000000000"<<std::endl;
// run prediction
   std::vector<float> edge_predictions = onnxRuntime_->run(input_names, data, input_shapes)[0];
 std::cout<<" &&&&&&&&&&&&&&& Something happening here in after the EdgePredictions*********************"<<std::endl;
 }

  // select tracks
  std::vector<reco::TransientTrack>&& seltks = theTrackFilter->select(t_tks);

  // clusterize tracks in Z
  std::vector<TransientVertex>&& clusters = theTrackClusterizer->vertices(seltks);
 std::cout<<" ARe we here in primaryVertexProducer "<<std::endl;
  if (fVerbose) {
    edm::LogPrint("PrimaryVertexProducer")
        << "Clustering returned " << clusters.size() << " clusters from " << seltks.size() << " selected tracks";
  }

  // vertex fits
  for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
    auto result = std::make_unique<reco::VertexCollection>();
    reco::VertexCollection& vColl = (*result);
    std::vector<TransientVertex> pvs;
 std::cout<<" Are we here 2222 in primaryVertexProducer "<<std::endl;
    if (algorithm->pv_fitter == nullptr) {
      pvs = clusters;
    } else {
      pvs = algorithm->pv_fitter->fit(seltks, clusters, beamSpot, algorithm->useBeamConstraint);
    }

    if (algorithm->pv_time_estimator != nullptr) {
      algorithm->pv_time_estimator->fill_vertex_times(pvs);
    }

    // sort vertices by pt**2  vertex
    if (pvs.size() > 1) {
      sort(pvs.begin(), pvs.end(), VertexHigherPtSquared());
    }

    // select and convert transient vertices to (reco) vertices
    for (std::vector<TransientVertex>::const_iterator iv = pvs.begin(); iv != pvs.end(); iv++) {
      if (iv->isValid() && (iv->degreesOfFreedom() >= algorithm->minNdof)) {
        reco::Vertex v = *iv;
        if (!validBS || ((*(algorithm->vertexSelector))(v, beamVertexState))) {
          vColl.push_back(v);
        }
      }
    }

    if (fVerbose) {
      edm::LogPrint("PrimaryVertexProducer") << "PrimaryVertexProducer \"" << algorithm->label << "\" contains "
                                             << pvs.size() << " reco::Vertex candidates";
    }

    if (clusters.size() > 2 && clusters.size() > 2 * pvs.size()) {
      edm::LogWarning("PrimaryVertexProducer")
          << "More than 50% of candidate vertices lost (" << pvs.size() << " out of " << clusters.size() << ")";
    }

    if (pvs.empty() && seltks.size() > 5) {
      edm::LogWarning("PrimaryVertexProducer")
          << "No vertex found with " << seltks.size() << " tracks and " << clusters.size() << " vertex candidates";
    }

    if (vColl.empty()) {
      GlobalError bse(beamSpot.rotatedCovariance3D());
      if ((bse.cxx() <= 0.) || (bse.cyy() <= 0.) || (bse.czz() <= 0.)) {
        AlgebraicSymMatrix33 we;
        we(0, 0) = 10000;
        we(1, 1) = 10000;
        we(2, 2) = 10000;
        vColl.push_back(reco::Vertex(beamSpot.position(), we, 0., 0., 0));
        edm::LogWarning("PrimaryVertexProducer") << "Zero recostructed vertices, will put reco::Vertex derived from "
                                                    "dummy/fake BeamSpot into Event, BeamSpot has invalid errors: "
                                                 << bse.matrix();
      } else {
        vColl.push_back(reco::Vertex(beamSpot.position(), beamSpot.rotatedCovariance3D(), 0., 0., 0));
        if (fVerbose) {
          edm::LogWarning("PrimaryVertexProducer")
              << "Zero recostructed vertices, will put reco::Vertex derived from BeamSpot into Event.";
        }
      }
    }
std::cout<<" ARe we here 333 in primaryVertexProducer "<<std::endl;
    if (fVerbose) {
      int ivtx = 0;
      for (reco::VertexCollection::const_iterator v = vColl.begin(); v != vColl.end(); ++v) {
        edm::LogPrint("PrimaryVertexProducer")
            << "recvtx " << std::setw(3) << std::fixed << ivtx++ << " #trk " << std::setw(3) << v->tracksSize()
            << " chi2 " << std::setw(5) << std::setprecision(1) << v->chi2() << " ndof " << std::setw(5)
            << std::setprecision(1) << v->ndof() << " x " << std::setw(7) << std::setprecision(4) << v->position().x()
            << " dx " << std::setw(6) << std::setprecision(4) << v->xError() << " y " << std::setw(7)
            << std::setprecision(4) << v->position().y() << " dy " << std::setw(6) << std::setprecision(4)
            << v->yError() << " z " << std::setw(8) << std::setprecision(4) << v->position().z() << " dz "
            << std::setw(6) << std::setprecision(4) << v->zError();
        if (v->tError() > 0) {
          edm::LogPrint("PrimaryVertexProducer") << " t " << std::setw(6) << std::setprecision(3) << v->t() << " dt "
                                                 << std::setw(6) << std::setprecision(3) << v->tError();
        }
      }
    }

    iEvent.put(std::move(result), algorithm->label);
  }
}

void PrimaryVertexProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription psd_pv_time;
  {
    edm::ParameterSetDescription psd1;
    VertexTimeAlgorithmFromTracksPID::fillPSetDescription(psd1);
    psd_pv_time.add<edm::ParameterSetDescription>("fromTracksPID", psd1);

    edm::ParameterSetDescription psd2;
    VertexTimeAlgorithmLegacy4D::fillPSetDescription(psd2);
    psd_pv_time.add<edm::ParameterSetDescription>("legacy4D", psd2);
  }
  psd_pv_time.add<std::string>("algorithm", "");  // default = none

  // vertex collections
  edm::ParameterSetDescription desc;
  {
    edm::ParameterSetDescription vpsd1;
    vpsd1.add<double>("maxDistanceToBeam", 1.0);
    vpsd1.add<std::string>("algorithm", "AdaptiveVertexFitter");
    vpsd1.add<bool>("useBeamConstraint", false);
    vpsd1.add<std::string>("label", "");
    vpsd1.add<double>("chi2cutoff", 2.5);
    vpsd1.add<double>("zcutoff", 1.0);
    vpsd1.add<double>("mintrkweight", 0.0);
    vpsd1.add<double>("minNdof", 0.0);
    vpsd1.add<edm::ParameterSetDescription>("vertexTimeParameters", psd_pv_time);

    // two default values : with- and without beam constraint
    std::vector<edm::ParameterSet> temp1;
    temp1.reserve(2);
    {
      edm::ParameterSet temp2;
      temp2.addParameter<double>("maxDistanceToBeam", 1.0);
      temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
      temp2.addParameter<bool>("useBeamConstraint", false);
      temp2.addParameter<std::string>("label", "");
      temp2.addParameter<double>("chi2cutoff", 2.5);
      temp2.addParameter<double>("zcutoff", 1.0);
      temp2.addParameter<double>("mintrkweight", 0.);
      temp2.addParameter<double>("minNdof", 0.0);
      edm::ParameterSet temp_vertexTime;
      temp_vertexTime.addParameter<std::string>("algorithm", "");
      temp2.addParameter<edm::ParameterSet>("vertexTimeParameters", temp_vertexTime);
      temp1.push_back(temp2);
    }
    {
      edm::ParameterSet temp2;
      temp2.addParameter<double>("maxDistanceToBeam", 1.0);
      temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
      temp2.addParameter<bool>("useBeamConstraint", true);
      temp2.addParameter<std::string>("label", "WithBS");
      temp2.addParameter<double>("chi2cutoff", 2.5);
      temp2.addParameter<double>("zcutoff", 1.0);
      temp2.addParameter<double>("mintrkweight", 0.);
      temp2.addParameter<double>("minNdof", 2.0);
      edm::ParameterSet temp_vertexTime;
      temp_vertexTime.addParameter<std::string>("algorithm", "");
      temp2.addParameter<edm::ParameterSet>("vertexTimeParameters", temp_vertexTime);
      temp1.push_back(temp2);
    }
    desc.addVPSet("vertexCollections", vpsd1, temp1);
  }
  desc.addUntracked<bool>("verbose", false);
  {
    edm::ParameterSetDescription psd0;
    HITrackFilterForPVFinding::fillPSetDescription(psd0);  // extension of TrackFilterForPVFinding
    desc.add<edm::ParameterSetDescription>("TkFilterParameters", psd0);
  }
  desc.add<edm::InputTag>("beamSpotLabel", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("TrackLabel", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("TrackTimeResosLabel", edm::InputTag("dummy_default"));                         // 4D only
  desc.add<edm::InputTag>("TrackTimesLabel", edm::InputTag("dummy_default"));                             // 4D only
  desc.add<edm::InputTag>("trackMTDTimeQualityVMapTag", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));  // 4D only
  desc.add<edm::FileInPath>("onnxModelPath", edm::FileInPath("RecoVertex/PrimaryVertexProducer/data/model.onnx"))->setComment("Path to GNN (as ONNX model)");
  desc.add<edm::InputTag>("timingSoA", edm::InputTag("dummy_default"));
  {
    edm::ParameterSetDescription psd0;
    {
      edm::ParameterSetDescription psd1;
      DAClusterizerInZ_vect::fillPSetDescription(psd1);

      edm::ParameterSetDescription psd2;
      DAClusterizerInZT_vect::fillPSetDescription(psd2);

      edm::ParameterSetDescription psd3;
      GapClusterizerInZ::fillPSetDescription(psd3);

      psd0.ifValue(
          edm::ParameterDescription<std::string>("algorithm", "DA_vect", true),
          "DA_vect" >> edm::ParameterDescription<edm::ParameterSetDescription>("TkDAClusParameters", psd1, true) or
              "DA2D_vect" >>
                  edm::ParameterDescription<edm::ParameterSetDescription>("TkDAClusParameters", psd2, true) or
              "gap" >> edm::ParameterDescription<edm::ParameterSetDescription>("TkGapClusParameters", psd3, true));
    }
    desc.add<edm::ParameterSetDescription>("TkClusParameters", psd0);
  }

  desc.add<bool>("isRecoveryIteration", false);
  desc.add<edm::InputTag>("recoveryVtxCollection", {""});
  desc.add<bool>("useMVACut", false);
  desc.add<double>("minTrackTimeQuality", 0.8);
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PrimaryVertexProducer);
