#ifndef RecoVertex_PrimaryVertexProducer_interface_TrackGraph_h
#define RecoVertex_PrimaryVertexProducer_interface_TrackGraph_h

#include "DataFormats/TrackReco/interface/Track.h"
#include <unordered_set>
#include <vector>

class Node {
public:
  Node() = default;
  Node(unsigned index, float zPosition , bool isTrack = true)  : index_(index), zPosition_(zPosition), isTrack_(isTrack), alreadyVisited_(false) {}; 
  void addInner(unsigned int track_id) { innerNodes_.push_back(track_id); addNeighbour(track_id);}
  void addNeighbour(unsigned int track_id, double weight = 0.0) {
    neighboursId_.push_back(track_id);
    weights_.push_back(weight);
  }
   float getZPosition() const{
     return zPosition_;
 } 
  void updateWeight(unsigned int neighborIndex, float weight) {
        // Find the neighbor index in the list of neighbors
        auto it = std::find(neighboursId_.begin(), neighboursId_.end(), neighborIndex);
        if (it != neighboursId_.end()) {
            // Update the weight if the neighbor is found
            size_t index = std::distance(neighboursId_.begin(), it);
            weights_[index] = weight;
        }
    }
  
 float getWeight(unsigned int neighborIndex) const {
    // Find the neighbor in the list of neighbors
    auto it = std::find(neighboursId_.begin(), neighboursId_.end(), neighborIndex);
    if (it != neighboursId_.end()) {
        // Find the index of the neighbor and return the corresponding weight
        size_t index = std::distance(neighboursId_.begin(), it);
        return weights_[index];
    }
    return 0.0f; // Return 0 if no weight exists (no edge)
}
  const unsigned int getId() const { return index_; }
  std::vector<unsigned int> getNeighbours() const { return neighboursId_; }
  std::vector<unsigned int> getInner() const { return innerNodes_; }
  
  void findSubComponents(std::vector<Node>& graph, std::vector<int>& subComponent, float threshold) {
    if (!alreadyVisited_) {
      alreadyVisited_ = true;
      subComponent.push_back(index_);
      // Using a const iterator since we don't intend to modify elements
      auto neighbourIt = neighboursId_.cbegin(); // Iterator for neighboursId_
      auto weightIt = weights_.cbegin(); // Iterator for weights_
      // Loop through both vectors simultaneously
      for (; neighbourIt != neighboursId_.cend() && weightIt != weights_.cend(); ++neighbourIt, ++weightIt) {
          int neighbour = *neighbourIt;

          double weight = *weightIt;
          if (weight >= threshold){
    graph[neighbour].findSubComponents(graph, subComponent, threshold);
	    //graph[0].findSubComponents(graph, subComponent, threshold) will start the DFS from node 0, identifying all connected nodes that meet the weight threshold and storing their indices in subComponent.
          }
      }
    }
  }
  std::vector<int> regionQuery(const std::vector<Node>& graph, float eps) const {
        std::vector<int> neighbors;

        // Iterate through all the nodes in the graph
        for (size_t j : getInner()) {
            if (j == this->getId()) continue;  // Skip the current node
            
            float edgeWeight = this->getWeight(j);

            if (edgeWeight > 0) {
                float distance = 1.0f - edgeWeight; 
                if (distance <= eps) {
                    neighbors.push_back(j);  // Add this neighbor to the list
                }
            }
        }

        return neighbors;
    }

  bool isAlreadyVisited() const {
        return alreadyVisited_;
    }

    void setAlreadyVisited(bool visited) {
        alreadyVisited_ = visited;
    }

  ~Node() = default;

private:
  unsigned index_;
  float zPosition_; 
  bool isTrack_;

  std::vector<unsigned int> neighboursId_;
  std::vector<float> weights_;
  std::vector<unsigned int> innerNodes_;
  bool alreadyVisited_;

  //bool areCompatible(const std::vector<Node>& graph, const unsigned int& outerNode) { return true; };

};

class TrackGraph {
public:
  TrackGraph() = default;
  TrackGraph(std::vector<Node> &n) { nodes_ = n; };
  
  const std::vector<Node>& getNodes() const { return nodes_; }
  const Node& getNode(int i) const { return nodes_[i]; }
  
  void setEdgeWeight(unsigned int nodeIndexI, unsigned int nodeIndexJ, float weight) {
        // Check if the node indices are valid
        if (nodeIndexI >= nodes_.size() || nodeIndexJ >= nodes_.size()) {
            // Handle invalid node indices
            return;
        }
        // Update the weight for the edge between node i and node j
        nodes_[nodeIndexI].updateWeight(nodeIndexJ, weight);
        // If j is bidirectionally connected, update the weight for j -> i as well
        nodes_[nodeIndexJ].updateWeight(nodeIndexI, weight);
  }
 // DFS algorithm 
    std::vector<std::vector<int>> findSubComponents(float threshold) {
    std::vector<std::vector<int>> components;
for (auto& node : nodes_) {
            node.setAlreadyVisited(false);
        }
    
    for (auto const& node: nodes_) {
      auto const id = node.getId();
      std::vector<int> tmpSubComponents;
      nodes_[id].findSubComponents(nodes_, tmpSubComponents, threshold);
      if (!tmpSubComponents.empty()) {
        components.push_back(tmpSubComponents);
      }
    }
    return components;
  }

   // DBSCAN Neighbour of Neighbour i.e ORIGINAL DBSCAN 
    std::vector<std::vector<int>> dbscanClustering(float eps, unsigned int minPts){
        int numNodes = nodes_.size();
        std::vector<int> labels(numNodes, -1); // Initialize all labels as -1 (unvisited)
        std::vector<std::vector<int>> clusters; // To hold the clusters
        int clusterID = 0;

        for (auto const& node: nodes_) {
                auto const i = node.getId();
            if (labels[i] != -1) {  // Node already processed
              continue;
              }
            // Find neighbors based on inverse edge prediction score
            std::vector<int> neighbors = node.regionQuery(nodes_, eps); //regionQuery(i, eps);// edgePredictionScores);
            if (neighbors.size() < minPts) {
                labels[i] = -2; // Label as noise
            }

            else {
            clusters.emplace_back();
            labels[i] = clusterID;  // Mark the current node as part of the current cluster
            clusters[clusterID].push_back(i);
            std::vector<int> toProcess = neighbors;
	    size_t j = 0;
            while (j < toProcess.size()) {
                int neighbor = toProcess[j];

                if (labels[neighbor] == -1) {  // If unvisited
                    labels[neighbor] = clusterID;
                    clusters[clusterID].push_back(neighbor);

                    std::vector<int> neighborsOfNeighbor = nodes_[neighbor].regionQuery(nodes_, eps); //regionQuery(neighbor, eps);
                  if (neighborsOfNeighbor.size() >= minPts) {  // Check if the neighbor has enough points to expand
                        // Add all the neighbors of the current neighbor to the list
                       for (int nn : neighborsOfNeighbor) {
                          if (labels[nn] == -1 && std::find(toProcess.begin(), toProcess.end(), nn) == toProcess.end()) {
                          toProcess.push_back(nn);
                    }
                  }
              }
        }
       j++;
      }
	   if (clusters[clusterID].size() < minPts) {
                for (int node : clusters[clusterID]) {
                    labels[node] = -2; // Label as noise
                }
                clusters.pop_back(); // Remove the invalid cluster
            } else {
                clusterID++; // Increment cluster ID for the next cluster
            }
        }
    }
        return clusters; // Return the clusters containing the node indices
    } 


	std::vector<std::vector<int>> dbscanClustering_Improvised(float eps, unsigned int minPts){
        int numNodes = nodes_.size();
        std::vector<int> labels(numNodes, -1); // Initialize all labels as -1 (unvisited)
        std::vector<std::vector<int>> clusters; // To hold the clusters
        int clusterID = 0;

        for (auto const& node: nodes_) {
                auto const i = node.getId();
            if (labels[i] != -1) {  // Node already processed
              continue;
              }
            std::vector<int> neighbors = node.regionQuery(nodes_, eps); //regionQuery(i, eps);// edgePredictionScores);
            if (neighbors.size() < minPts) {
                labels[i] = -2; // Label as noise
            }

            else {
            clusters.emplace_back();
            labels[i] = clusterID;  // Mark the current node as part of the current cluster
            clusters[clusterID].push_back(i);
            std::vector<int> toProcess = neighbors;
            size_t j = 0;
            while (j < toProcess.size()) {
                int neighbor = toProcess[j];
                 if (labels[neighbor] == -1) {  // If unvisited
	           labels[neighbor] = clusterID;
                    clusters[clusterID].push_back(neighbor);
              }
            j++;
          }

             if (clusters[clusterID].size() < minPts) {
                for (int node : clusters[clusterID]) {
                    labels[node] = -2; // Label as noise
                }
                clusters.pop_back(); // Remove the invalid cluster
            } else {
                clusterID++; // Increment cluster ID for the next cluster
            }
        }
    }
        return clusters; // Return the clusters containing the node indices
    }

  ~TrackGraph() = default;


private:
  std::vector<Node> nodes_;

};

#endif
