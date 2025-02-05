

#ifndef GNNClusterizer_h
#define GNNClusterizer_h

#include "RecoVertex/PrimaryVertexProducer/interface/TrackClusterizerInZ.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

#include <memory>
#include <vector>

using namespace cms::Ort;
class GNNClusterizer final : public TrackClusterizerInZ {
public:
  GNNClusterizer(const edm::ParameterSet& conf, const ONNXRuntime* onnxRuntime);
  ~GNNClusterizer() override {}

  std::vector<TransientVertex> vertices(const std::vector<reco::TransientTrack>& tracks) const override;
  std::vector<std::vector<reco::TransientTrack>> clusterize(
      const std::vector<reco::TransientTrack>& tracks) const override;
  static void fillPSetDescription(edm::ParameterSetDescription& desc);
  static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet& conf);
  static void globalEndJob(const ONNXRuntime* cache);

private:
  const ONNXRuntime* onnxRuntime_;
  std::string nnVersion_;
  double nnWorkingPoint_;
  std::string AlgoVersion_;
  bool verbose_;
  double zSep;
  double d0CutOff_;
  double t_beta_;
  double t_d_;
  int min_cluster_size_;
  double vertexSize_;
};

class UnionFind {
public:
  UnionFind(int n);
  int find(int x);
  void unite(int x, int y);

private:
  std::vector<int> parent;
};

#endif
