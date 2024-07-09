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
#include "vdt/vdtMath.h"
const int FEATURE_SHAPE_GNN_V1 = 21;
const int FEATURE_SHAPE_MLP = 30;
const int NUM_EDGE_FEATURES = 12;
double d0CutOff_ = 4.0;  // these variables are in the config file for DAclustering
double vertexSize_ = 0.06;

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

//PrimaryVertexProducer::PrimaryVertexProducer(const edm::ParameterSet& conf, cms::Ort::ONNXRuntime const *onnxRuntime)
PrimaryVertexProducer::PrimaryVertexProducer(const edm::ParameterSet& conf, const ONNXRuntime *cache) 
   : theTTBToken(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))), theConfig(conf) {
  fVerbose = conf.getUntrackedParameter<bool>("verbose", false);
  useMVASelection_ = conf.getParameter<bool>("useMVACut");
  nnVersion_ = conf.getParameter<std::string>("nnVersion");
  nnWorkingPoint_ = conf.getParameter<double>("nnWorkingPoint");
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
   /*  } else if (clusteringAlgorithm == "DA2D_vect") {
       theTrackClusterizer = new ClusterizerinGNN(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));*/
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
    trkMTDAssocToken = consumes<edm::ValueMap<int> >(conf.getParameter<edm::InputTag>("trackAssocSrc")); 
    MTDtimeToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("tmtdSrc"));
    sigmaMTDtimeToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("sigmatmtdSrc"));
    pathLengthToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("pathmtd"));
    btlMatchChi2Token = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("btlMatchChi2Src"));
    btlMatchTime_Chi2Token = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("btlMatchTimeChi2Src"));
    etlMatchChi2Token = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("etlMatchChi2Src"));
    etlMatchTime_Chi2Token = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("etlMatchTimeChi2Src"
));
    trkTimePiToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("tofPi"));
    trkTimeKToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("tofK"));
    trkTimePToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("tofP"));
    sigmaTrkTimePiToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("sigmatofpiSrc"));
    sigmaTrkTimeKToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("sigmatofkSrc"));
    sigmaTrkTimePToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("sigmatofpSrc"));
    npixBarrelToken = consumes<edm::ValueMap<int> >(conf.getParameter<edm::InputTag>("npixBarrelSrc"));
    npixEndcapToken = consumes<edm::ValueMap<int> >(conf.getParameter<edm::InputTag>("npixEndcapSrc"));
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

std::unique_ptr<ONNXRuntime> PrimaryVertexProducer::initializeGlobalCache(const edm::ParameterSet &conf) {
  return std::make_unique<ONNXRuntime>(conf.getParameter<edm::FileInPath>("onnxModelPath").fullPath());
}

void PrimaryVertexProducer::globalEndJob(const ONNXRuntime *cache) {}

namespace {
  inline double local_exp(double const& inp) { return vdt::fast_exp(inp); }
}

// can shift this function to TracksGraph.h as well
std::unique_ptr<TrackGraph>  PrimaryVertexProducer::produce_tracks_graph(const std::vector<reco::TransientTrack>& transientTracks) {
    std::vector<Node> allNodes;
//    for (size_t i = 0; i < transientTracks.size(); ++i) {
  for (size_t i = 0; i < 5000 ; ++i) {  
    Node tNode(i);
    const reco::TransientTrack& ttrack_node = transientTracks[i];
  //  for (size_t j = i + 1; j < transientTracks.size(); ++j) {
    for (size_t j = i + 1; j < 5000; ++j) {    
    const reco::TransientTrack& ttrack = transientTracks[j];
        double neighbour_cut = std::abs(ttrack_node.track().vz() - ttrack.track().vz());
        if (neighbour_cut < 0.3 ) {
          tNode.addInner(j);
    }
  }
    allNodes.push_back(tNode);
  }
     auto resultGraph = std::make_unique<TrackGraph>(allNodes);

    return resultGraph;
}
void PrimaryVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // get the BeamSpot, it will always be needed, even when not used as a constraint

  seltks.clear();
  //outTrack.clear();
  linkedTrackIdToInputTrackId.clear();
  resultTracks.clear();

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

    // Network input shapes.
    std::vector<float> features;
    std::vector<float> edge_features;
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
    auto const& trackMTDAssoc_ = iEvent.get(trkMTDAssocToken);
    auto trackMTDTimes_ = iEvent.get(MTDtimeToken);
    auto trackMTDTimesRes_ = iEvent.get(sigmaMTDtimeToken);
    auto pathlength_ = iEvent.get(pathLengthToken);
    auto btlMatch_ = iEvent.get(btlMatchChi2Token);
    auto btlMatchTime_ = iEvent.get(btlMatchTime_Chi2Token);
    auto etlMatch_ = iEvent.get(etlMatchChi2Token);   
    auto etlMatchTime_ = iEvent.get(etlMatchTime_Chi2Token);
    auto trkPiTime_ = iEvent.get(trkTimePiToken);
    auto trkKTime_ = iEvent.get(trkTimeKToken);
    auto trkPTime_ = iEvent.get(trkTimePToken);
    auto sigmatrkPiTime_ = iEvent.get(sigmaTrkTimePiToken);
    auto sigmatrkKTime_ = iEvent.get(sigmaTrkTimeKToken);
    auto sigmatrkPTime_ = iEvent.get(sigmaTrkTimePToken);
    auto npixbarrel_ = iEvent.get(npixBarrelToken); 
    auto npixendcap_ = iEvent.get(npixEndcapToken);
    auto  trackMTDTimeQualities_ = iEvent.get(trackMTDTimeQualityToken); 


    if (useMVASelection_) {
   //   trackMTDTimeQualities_ = iEvent.get(trackMTDTimeQualityToken);

      for (unsigned int i = 0; i < (*tks).size(); i++) {
        const reco::TrackRef ref(tks, i);
        auto const trkTimeQuality = trackMTDTimeQualities_[ref];
        if (trkTimeQuality < minTrackTimeQuality_) {
          trackTimes_[ref] = std::numeric_limits<double>::max();
        }
      }
      t_tks = (*theB).build(tks, beamSpot, trackTimes_, trackTimeResos_);
    } else {
      t_tks = (*theB).build(tks, beamSpot, trackTimes_, trackTimeResos_, trackMTDAssoc_, trackMTDTimes_, trackMTDTimesRes_, trackMTDTimeQualities_, pathlength_, btlMatch_, btlMatchTime_, etlMatch_, etlMatchTime_, trkPiTime_, trkKTime_, trkPTime_, sigmatrkPiTime_, sigmatrkKTime_, sigmatrkPTime_, npixbarrel_, npixendcap_);
    }
  } else {
    t_tks = (*theB).build(tks, beamSpot);
  }

std::vector<reco::TransientTrack>&& filterTracks = theTrackFilter->select(t_tks);

  long int N = 5000; //t_tks.size();
  if (useTransientTrackTime_) {
    auto trackgraph = produce_tracks_graph(filterTracks);
    TrackGraph *trkgrp = trackgraph.get();
   for (size_t i =0; i < 5000; ++i){ 
 //  for (size_t i = 0; i < filterTracks.size() ; ++i) {
        const reco::TransientTrack& ttrack1 = filterTracks[i];
        float trackpt = ttrack1.track().pt();
        float tracketa = ttrack1.track().eta();
        float trackphi = ttrack1.track().phi();
        float trackvz = ttrack1.track().vz();
        float track_dzError = ttrack1.track().dzError();
        float track_ndof = ttrack1.track().ndof();
        float track_chi2 = ttrack1.track().chi2();
        float mtdtime = ttrack1.MTDtime();
        float mtdtimeErr = ttrack1.MTDtimeErr();
        float mvaquality = ttrack1.MVAquality();
        float pathlength = ttrack1.pathLength();
        float btlMatchchi2 = ttrack1.btlMatch_chi2();
        float btlMatchTimechi2 = ttrack1.btlMatchTime_chi2();
        float etlMatchchi2 = ttrack1.etlMatch_chi2();
        float etlMatchTimechi2 = ttrack1.etlMatchTime_chi2();
        float Time_pi_hypo = ttrack1.trackTime_pi();
        float Time_k_hypo = ttrack1.trackTime_k();
        float Time_p_hypo = ttrack1.trackTime_p();
        float npixBarrel = ttrack1.nPixBarrel();
        float npixEndcap = ttrack1.nPixEndcap();
        float sigmaTime_p = ttrack1.sigma_time_p();
        //float trackMTD = ttrack1.trackAsocMTD();
       // float sigmaTime_pi = ttrack1.sigma_time_pi();
       // float sigmaTime_k = ttrack1.sigma_time_k();

//Variables we have in model : pt, eta, phi, z_pca, t_pi, t_k, t_p, mva, btl*4, path lenght, npixBarrel, npix Ebdcap, mtdtime, sigma mtdtime, dz, sigma tof_p, trk_ndof, trk_chi2
        features.push_back(trackpt);
        features.push_back(tracketa);
        features.push_back(trackphi);
        features.push_back(trackvz);
        features.push_back(Time_pi_hypo);
        features.push_back(Time_k_hypo);
        features.push_back(Time_p_hypo);
        features.push_back(mvaquality);
        features.push_back(btlMatchchi2);
        features.push_back(btlMatchTimechi2);
        features.push_back(etlMatchchi2);
        features.push_back(etlMatchTimechi2);
        features.push_back(pathlength);
        features.push_back(npixBarrel);
        features.push_back(npixEndcap);
        features.push_back(mtdtime);
        features.push_back(track_dzError);
        features.push_back(mtdtimeErr);
        features.push_back(sigmaTime_p);
        features.push_back(track_ndof);
        features.push_back(track_chi2);


for (size_t j = i + 1; j < 5000; ++j){
  //  for (auto &j : trkgrp->getNode(i).getInner()){
        const reco::TransientTrack& ttrack2 = filterTracks[j];


        double neighbour_cut = std::abs(ttrack1.track().vz() - ttrack2.track().vz());
        if (neighbour_cut < 0.3 ) {
        edges_src.push_back(static_cast<float>(j)); // this should be neighbour of the track 
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
         edge_features.push_back(time_k_p_diff);
         edge_features.push_back(time_p_k_diff);
         edge_features.push_back(time_p_pi_diff);
    }
   }
  }
   auto numEdges = static_cast<int>(edges_src.size());
   if (numEdges < 1) {
       edm::LogPrint("PrimaryVertexProducer")  << "No edges for the event - no linking is done." ;
        linkedTrackIdToInputTrackId.resize(N);
        std::vector<unsigned int> linkedTracks;

        for (int track_id = 0; track_id < N; track_id++)
        {
            linkedTrackIdToInputTrackId[track_id].push_back(track_id);
	    //seltks.push_back(filterTracks[track_id]);
            linkedTracks.push_back(seltks.size());
        }
    }
    data.clear();

    input_shapes.push_back({1, N, shapeFeatures});
    data.emplace_back(features);
    input_shapes.push_back({1, 2, numEdges});
    data.emplace_back(edges_src);
    for (auto &dst : edges_dst)
    {
            data.back().push_back(dst);
    }
    input_shapes.push_back({1, numEdges, NUM_EDGE_FEATURES});
    data.emplace_back(edge_features);
 
// run prediction
 if (numEdges >=1 ) {
 auto edge_predictions = globalCache()->run(input_names, data, input_shapes)[0];

  // Update graph weights
    for (int i = 0; i < static_cast<int>(edge_predictions.size()); i++){
   
   if  (fVerbose) { 
	    edm::LogPrint("PrimaryVertexProducer") << "Network output for edge " << data[1][i] << "-" << data[1][numEdges + i] << " is: " << edge_predictions[i];
       }
            trkgrp->setEdgeWeight(data[1][i], data[1][numEdges + i], edge_predictions[i]);
    }

    auto connected_components = trkgrp->findSubComponents(nnWorkingPoint_);
    edm::LogPrint("PrimaryVertexProducer") <<" connected_components.size() "<<connected_components.size();
    int id = 0;
   
    std::vector<TrackCollection> outTracks; 
    linkedTrackIdToInputTrackId.resize(connected_components.size());
    for (auto &component : connected_components)
    {
        std::vector<unsigned int> linkedTracks;
	TrackCollection outTrack;
	std::vector<reco::TransientTrack> tracks;
        for (auto &track_id : component)
        {
            linkedTrackIdToInputTrackId[id].push_back(track_id);
            outTrack.addItem(filterTracks[track_id]);
	    tracks.push_back(filterTracks[track_id]);
            //seltks[id].push_back(filterTracks[track_id]);
        }
        id++;
	outTracks.push_back(outTrack);
	seltks.push_back(tracks);
    }	
//Lets clusterize here as in DA
//
    const unsigned int nv = connected_components.size(); // this nv is number of components whose model output is above a given threshold, and we can cluster these. 
    const unsigned int nt = seltks.size(); // this is total number of tracks 
    
// In DA, the z-position is typically computed as a weighted mean of the z-positions of the associated tracks,
    std::vector<double> zvtx(connected_components.size());
    std::vector<double> tvtx(connected_components.size()); 
    std::vector<double> dz2;
    std::vector<double> zpca;
    std::vector<double> t_tkwt;
    std::vector<double> t_dz2;
 for (unsigned int k =0; k < nv; k++) {
    double sumwz = 0;
    double sumwt = 0;
    double sumw_z = 0;
    double sumw_t = 0;
    const auto &tracks = seltks[k];   
    dz2.resize(tracks.size());
    zpca.resize(tracks.size());
    t_tkwt.resize(tracks.size());
    t_dz2.resize(tracks.size());
     for (size_t index = 0; index < tracks.size(); ++index) 
    {
     const auto &track = tracks[index];
      if (d0CutOff_ > 0) {
      Measurement1D atIP = track.stateAtBeamLine().transverseImpactParameter();  // error contains beamspot
      t_tkwt[index] = 1. / (1. + local_exp(std::pow(atIP.value() / atIP.error(), 2) -
                                    std::pow(d0CutOff_, 2)));  // reduce weight for high ip tracks
      if (edm::isNotFinite(t_tkwt[index]) || t_tkwt[index] < std::numeric_limits<double>::epsilon()) {
        edm::LogWarning("PrimaryVertexProducer") << "rejected track t_tkwt " << t_tkwt[index];
        continue;  // usually is > 0.99
      }
    }
      auto const& t_mom = track.stateAtBeamLine().trackStateAtPCA().momentum();
      reco::BeamSpot beamspot = track.stateAtBeamLine().beamSpot();
      t_dz2[index] = std::pow(track.track().dzError(), 2)  // track errror
                   + (std::pow(beamspot.BeamWidthX() * t_mom.x(), 2) + std::pow(beamspot.BeamWidthY() * t_mom.y(), 2)) *
                         std::pow(t_mom.z(), 2) / std::pow(t_mom.perp2(), 2)  // beam spot width
                   + std::pow(vertexSize_, 2);  // intrinsic vertex size, safer for outliers and short lived decays
      t_dz2[index] = 1. / t_dz2[index];
      if (edm::isNotFinite(t_dz2[index]) || t_dz2[index] < std::numeric_limits<double>::min()) {
      edm::LogWarning("PrimaryVertexProducer") << "rejected track t_dz2 " << t_dz2[index];
      continue;
    }

      zpca[index] = track.track().vz(); 

      double w_z = t_tkwt[index] * t_dz2[index];
      double w_t = t_tkwt[index] * t_dz2[index];
      sumwz +=  w_z * zpca[index];
      sumwt += w_t * zpca[index];
      sumw_z += w_z;
      sumw_t += w_t;
    }
    zvtx[k] = sumwz / sumw_z;
 }
    GlobalError dummyError(0.01, 0, 0.01, 0., 0., 0.01);
    for (unsigned int k = 0; k < nv; k++) {

     if (!linkedTrackIdToInputTrackId[k].empty()) {
      GlobalPoint pos(0, 0, zvtx[k]);
      std::vector<reco::TransientTrack> vertexTracks;
      for (const auto& track : outTracks[k].tracks) {
        // Here, 'track' is each TransientTrack in the k-th TrackCollection
        vertexTracks.push_back(track);
    }

      TransientVertex v(pos, dummyError, vertexTracks, 0);
      clusters_time.push_back(v);
    }
  }
     
  }// loop for numEdges 
} //TransientTrack time loop 
  

   //std::vector<TransientVertex> clusters; 
  // clusterize tracks in Z
   if (useTransientTrackTime_) {
	   clusters = clusters_time;
  }  else {
	clusters = theTrackClusterizer->vertices(filterTracks);
  }
  if (fVerbose) {
    edm::LogPrint("PrimaryVertexProducer")
        << " bool useTransientTrackTime_ "<<useTransientTrackTime_<<" Clustering returned " << clusters.size() << " clusters from " << filterTracks.size() << " selected tracks";
  }

  // vertex fits
  for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
    auto result = std::make_unique<reco::VertexCollection>();
    reco::VertexCollection& vColl = (*result);
    std::vector<TransientVertex> pvs;
    if (algorithm->pv_fitter == nullptr) {
      pvs = clusters;
    } else {
      if (useTransientTrackTime_) {
	      for (const auto& tracks : seltks) {
            // Fit using tracks (vector of reco::TransientTrack)
            pvs = algorithm->pv_fitter->fit(tracks, clusters, beamSpot, algorithm->useBeamConstraint);
        }
     // pvs = algorithm->pv_fitter->fit(seltks, clusters, beamSpot, algorithm->useBeamConstraint);
      }
       else {
          pvs = algorithm->pv_fitter->fit(filterTracks, clusters, beamSpot, algorithm->useBeamConstraint);
      }
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

    if (pvs.empty() && filterTracks.size() > 5) {
      edm::LogWarning("PrimaryVertexProducer")
          << "No vertex found with " << filterTracks.size() << " tracks and " << clusters.size() << " vertex candidates";
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
  desc.add<edm::FileInPath>("onnxModelPath", edm::FileInPath("RecoVertex/PrimaryVertexProducer/data/model_v6version24June.onnx"))->setComment("Path to GNN (as ONNX model)");
  desc.add<edm::InputTag>("timingSoA", edm::InputTag("mtdSoA"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"));
  desc.add<edm::InputTag>("tmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmatmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("pathmtd", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("btlMatchChi2Src", edm::InputTag("trackExtenderWithMTD", "btlMatchChi2"));
  desc.add<edm::InputTag>("btlMatchTimeChi2Src", edm::InputTag("trackExtenderWithMTD", "btlMatchTimeChi2"));
  desc.add<edm::InputTag>("etlMatchChi2Src", edm::InputTag("trackExtenderWithMTD", "etlMatchChi2"));
  desc.add<edm::InputTag>("etlMatchTimeChi2Src", edm::InputTag("trackExtenderWithMTD", "etlMatchTimeChi2"));
  desc.add<edm::InputTag>("tofPi", edm::InputTag("trackExtenderWithMTD:generalTrackTofPi"));
  desc.add<edm::InputTag>("tofK", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"));
  desc.add<edm::InputTag>("tofP", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"));
  desc.add<edm::InputTag>("sigmatofpiSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"));
  desc.add<edm::InputTag>("sigmatofkSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"));
  desc.add<edm::InputTag>("sigmatofpSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"));
  desc.add<edm::InputTag>("npixBarrelSrc", edm::InputTag("trackExtenderWithMTD", "npixBarrel"));
  desc.add<edm::InputTag>("npixEndcapSrc", edm::InputTag("trackExtenderWithMTD", "npixEndcap"));
  desc.add<std::string>("nnVersion", "gnn_v1") // gnn_v1, mlp_no_edge_features, mlp_edge_features
        ->setComment("GNN version tag.");
   desc.add<double>("nnWorkingPoint", 0.80)
        ->setComment("Working point of the GNN (in [0, 1]). GNN score above WP will attempt to supercluster.");
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
