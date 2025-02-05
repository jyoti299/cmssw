#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GNNClusterizer.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ_vect.h"
#include "RecoVertex/PrimaryVertexProducer/interface/ObjectCondensationClustering.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DBSCANClusterizer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "vdt/vdtMath.h"
#include <iostream>

#include <cmath>
#include <cassert>
#include <limits>
#include <iomanip>
#include <map>
#include <vector>
#include <math.h>
#include <queue>

#define cputime
#ifdef cputime
#include <chrono>
typedef std::chrono::duration<int, std::micro> microseconds_type;
#endif

using namespace std;
enum class ClusterAlgo {
    Object,
    DBSCAN
};
std::unordered_map<std::string, ClusterAlgo> algoMap = {
    {"ObjectCond", ClusterAlgo::Object},
    {"dbscan", ClusterAlgo::DBSCAN}
};

GNNClusterizer::GNNClusterizer(const edm::ParameterSet& conf, const ONNXRuntime* onnxRuntime)
    : onnxRuntime_(onnxRuntime),
      nnVersion_(conf.getParameter<std::string>("nnVersion")),
      nnWorkingPoint_(conf.getParameter<double>("nnWorkingPoint")),
      AlgoVersion_(conf.getParameter<std::string>("AlgoVersion")),	
      verbose_(conf.getUntrackedParameter<bool>("verbose", true)),
      zSep(conf.getParameter<double>("zSeparation")),
      d0CutOff_(conf.getParameter<double>("d0CutOff")),
      t_beta_(conf.getParameter<double>("t_beta")),
      t_d_(conf.getParameter<double>("t_d")),
      min_cluster_size_(conf.getParameter<int>("min_cluster_size")),	
      vertexSize_(conf.getParameter<double>("vertexSize")) {}

void GNNClusterizer::globalEndJob(const ONNXRuntime* cache) {}

namespace {
  inline double local_exp(double const& inp) { return vdt::fast_exp(inp); }
}  // namespace

std::unique_ptr<ONNXRuntime> GNNClusterizer::initializeGlobalCache(const edm::ParameterSet& conf) {
  return std::make_unique<ONNXRuntime>(conf.getParameter<edm::FileInPath>("onnxModelPath").fullPath());
}

std::vector<TransientVertex> GNNClusterizer::vertices(const std::vector<reco::TransientTrack>& tracks) const {
  cms::Ort::FloatArrays input_values;     // Stores float inputs for ONNX
  std::vector<TransientVertex> clusters;  // Output cluster vertices
  int N = tracks.size();                  // Number of tracks
  int shapeFeatures = 1;                  // Only 'vz' as feature
  std::string clusterAlgo_;

  auto it = algoMap.find(AlgoVersion_);
     if (it != algoMap.end()) {
      clusterAlgo_ = it->first;
    } else {
      throw std::invalid_argument("Invalid algorithm name: " + AlgoVersion_);
    }

  if (N == 0) {
    return clusters;  // Return empty if no tracks
  }

  // Prepare input data: features ('x'), uncertainties ('dz'), and batch indices ('batch')
  std::vector<float> features;       // 'x'
  std::vector<float> uncertainties;  // 'dz'
  std::vector<float> batch;          // 'batch' as int64_t

  features.reserve(N * shapeFeatures);
  uncertainties.reserve(N);
  batch.reserve(N);

  for (size_t i = 0; i < tracks.size(); ++i) {
    const reco::TransientTrack& ttrack = tracks[i];
    float track_vz = ttrack.track().vz();            // Feature: vz
    float track_dzError = ttrack.track().dzError();  // Uncertainty: dzError

    features.push_back(track_vz);            // Add vz to features
    uncertainties.push_back(track_dzError);  // Add dzError to uncertainties
    batch.push_back(0);                      // Batch index (int64_t)
  }

  // Define input names consistent with ONNX model
  std::vector<std::string> input_names = {"x", "batch", "dz"};

  // Add float inputs to input_values_float
  input_values.emplace_back(features);  // 'x'   // 'dz'

  // Add int64_t input to input_values_int
  input_values.emplace_back(batch);          // 'batch'
  input_values.emplace_back(uncertainties);  // 'dz'
  // Define input dimensions
  std::vector<std::vector<long int>> input_dims;
  input_dims.push_back({1, N, shapeFeatures});  // 'x': [1, N, shapeFeatures]
  input_dims.push_back({1, N});                 // 'batch': [1, N]
  input_dims.push_back({1, N});                 // 'dz': [1, N]

  // Define batch size (not used directly here)
  //int64_t batch_size = 1;

  // Define output names consistent with ONNX model
  std::vector<std::string> output_names = {"beta", "embeddings"};

#ifdef cputime
  std::chrono::duration<int, std::micro> tcpu_inference(0), tcpu_out(0);
  auto start_inference = std::chrono::high_resolution_clock::now();
#endif

  // Run ONNX inference
  auto output_values = onnxRuntime_->run(input_names, input_values, input_dims);

#ifdef cputime
  auto stop_inference = std::chrono::high_resolution_clock::now();
  tcpu_inference = std::chrono::duration_cast<std::chrono::microseconds>(stop_inference - start_inference);
  edm::LogInfo("PrimaryVertexProducer") << "###TIME inference " << tcpu_inference;
#endif

  std::vector<float>& beta_predictions = output_values[0];
  std::vector<float>& embeddings_flat = output_values[1];

  if (verbose_) {
    std::cout << "beta_predictions: [";
    for (size_t i = 0; i < std::min(size_t(10), beta_predictions.size()); ++i) {
      std::cout << beta_predictions[i];
      if (i != beta_predictions.size() - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "]" << std::endl;

    std::cout << "embeddings_flat: [";
    for (size_t i = 0; i < std::min(size_t(10), embeddings_flat.size()); ++i) {
      std::cout << embeddings_flat[i];
      if (i != embeddings_flat.size() - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "]" << std::endl;
  }

  // Determine embedding dimension
  int embedding_dim = embeddings_flat.size() / N;
  std::cout << "Emb dim:" << embedding_dim << std::endl;
  // Reshape embeddings_flat into a 2D vector
  std::vector<std::vector<double>> embedding_vectors(N, std::vector<double>(embedding_dim));

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < embedding_dim; ++j) {
      embedding_vectors[i][j] = embeddings_flat[i * embedding_dim + j];
    }
  }
/*
  float beta_threshold = 0.000;  //nnWorkingPoint_; // Use nnWorkingPoint_ as the beta threshold
  std::vector<int> selected_nodes;

  for (int i = 0; i < N; ++i) {
    if (beta_predictions[i] > beta_threshold) {
      selected_nodes.push_back(i);
    }
  }

  if (selected_nodes.empty()) {
    return clusters;
  }
  */
  std::vector<float> selected_beta_values;
  std::vector<std::vector<double>> selected_embeddings;
  for (int i = 0; i < N; ++i) { 
        selected_embeddings.push_back(embedding_vectors[i]);
        selected_beta_values.push_back(beta_predictions[i]);
    }
  
#ifdef cputime
  auto start_out = std::chrono::high_resolution_clock::now();
#endif
  // DBSCAN 
 /* std::vector<std::vector<double>> dbscan_points;
  for (int idx : selected_nodes) {
    dbscan_points.push_back(embedding_vectors[idx]);
  }*/

  double eps = 0.00025;
  std::vector<int> labels;
  if (clusterAlgo_ == "dbscan") {
     DBSCANClusterizer dbscan(eps, min_cluster_size_);
     labels = dbscan.run_clustering(embedding_vectors);
  } else {
     ObjectCondensationClustering clustering(t_beta_, t_d_, min_cluster_size_);
     std::cout<<" t_beta_ "<<t_beta_<<" t_d_ "<<t_d_<<std::endl;
     labels = clustering.cluster(embedding_vectors,beta_predictions);
  }
     
#ifdef cputime
  auto stop_out = std::chrono::high_resolution_clock::now();
  tcpu_out = std::chrono::duration_cast<std::chrono::microseconds>(stop_out - start_out);
  edm::LogInfo("GNNClusterizer") << "###TIME out " << tcpu_out;
#endif

  std::map<int, std::vector<int>> clusters_map;  // Map from cluster id to list of track indices

  for (size_t i = 0; i < labels.size(); ++i) {
    int label = labels[i];
  //  int idx = selected_nodes[i];  // Original index in tracks
    if (label == -1) {
      continue;
    }
    clusters_map[label].push_back(i);
  }

  GlobalError dummyError(0.01, 0, 0.01, 0., 0., 0.01);

  for (const auto& cluster_pair : clusters_map) {
    const std::vector<int>& cluster_indices = cluster_pair.second;
    std::vector<reco::TransientTrack> cluster_tracks;
    double sumwz = 0.0;
    double sumw = 0.0;
    double t_tkwt = 1.0;
    for (int idx : cluster_indices) {
      const reco::TransientTrack& track = tracks[idx];
      cluster_tracks.push_back(track);

      // double weight = 1.0 / (track.track().dzError() * track.track().dzError());
      if (d0CutOff_ > 0) {
        Measurement1D atIP = track.stateAtBeamLine().transverseImpactParameter();  // error contains beamspot
        double t_tkwt = 1. / (1. + local_exp(std::pow(atIP.value() / atIP.error(), 2) -
                                             std::pow(d0CutOff_, 2)));  // reduce weight for high ip tracks
        if (edm::isNotFinite(t_tkwt) || t_tkwt < std::numeric_limits<double>::epsilon()) {
          edm::LogWarning("GNN Clusterizer") << "rejected track t_tkwt " << t_tkwt;
          continue;  // usually is > 0.99
        }
      }  // d0 cutoff

      auto const& t_mom = track.stateAtBeamLine().trackStateAtPCA().momentum();
      reco::BeamSpot beamspot = track.stateAtBeamLine().beamSpot();
      double t_dz2 = std::pow(track.track().dzError(), 2)  // track errror
                     +
                     (std::pow(beamspot.BeamWidthX() * t_mom.x(), 2) + std::pow(beamspot.BeamWidthY() * t_mom.y(), 2)) *
                         std::pow(t_mom.z(), 2) / std::pow(t_mom.perp2(), 2)  // beam spot width
                     + std::pow(vertexSize_, 2);  // intrinsic vertex size, safer for outliers and short lived decays
      t_dz2 = 1. / t_dz2;
      if (edm::isNotFinite(t_dz2) || t_dz2 < std::numeric_limits<double>::min()) {
        edm::LogWarning("GNN Clusterizer") << "rejected track t_dz2 " << t_dz2;
        continue;
      }

      double weight = t_tkwt * t_dz2;
      double z = track.track().vz();
      sumwz += weight * z;
      sumw += weight;
    }
    double z_cluster = sumwz / sumw;
    //if (verbose_) {
    //           std::cout << "z_cluster:" << z_cluster << std::endl;
    //           }
    GlobalPoint pos(0, 0, z_cluster);
    TransientVertex vertex(pos, dummyError, cluster_tracks, 0);
    clusters.push_back(vertex);
  }
  std::cout << "cluster size:" << clusters.size() << std::endl;
  return clusters;
}

std::vector<std::vector<reco::TransientTrack>> GNNClusterizer::clusterize(
    const std::vector<reco::TransientTrack>& tracks) const {
  return std::vector<std::vector<reco::TransientTrack>>();
}

void GNNClusterizer::fillPSetDescription(edm::ParameterSetDescription& desc) {
  DAClusterizerInZ_vect::fillPSetDescription(desc);
  desc.add<double>("zSeparation", 0.001)->setComment("Epsilon value for DBSCAN clustering");
  desc.addUntracked<bool>("verbose", true);
   desc.add<std::string>("AlgoVersion", "dbscan")
        ->setComment("Clustering Algorithm for GravNet");
  desc.add<edm::FileInPath>("onnxModelPath", edm::FileInPath("RecoVertex/PrimaryVertexProducer/data/model_v4_6.onnx"))
      ->setComment("Path to the ONNX model");
   desc.add<double>("t_beta", 0.11)
        ->setComment("Beta threshold for node selection.");
    desc.add<double>("t_d", 0.0009)
        ->setComment("Eucl Distance.");
    desc.add<int>("min_cluster_size", 4)
        ->setComment("Min tracks in a cluster");
  desc.add<std::string>("nnVersion", "gravnet_v1")->setComment("GNN version tag.");
  desc.add<double>("nnWorkingPoint", 0.00)->setComment("Beta threshold for node selection.");
}
