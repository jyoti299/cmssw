#ifndef RecoVertex_PrimaryVertexProducer_interface_TrackGraph_h
#define RecoVertex_PrimaryVertexProducer_interface_TrackGraph_h

#include "DataFormats/TrackReco/interface/Track.h"
#include <unordered_set>
#include <vector>
#include <TH1F.h>

class Node {
public:
  Node() = default;
  //Node(unsigned index, bool isTrack = true) : index_(index), isTrack_(isTrack), alreadyVisited_{false}{};
  Node(unsigned index, float zPosition , bool isTrack = true)  : index_(index), zPosition_(zPosition),  isTrack_(isTrack), alreadyVisited_(false) {}; //JB
  void addInner(unsigned int track_id) { innerNodes_.push_back(track_id); addNeighbour(track_id);}
  void addOuter(unsigned int track_id) { outerNodes_.push_back(track_id); addNeighbour(track_id);}
  void addNeighbour(unsigned int track_id, double weight = 0.0) {
    neighboursId_.push_back(track_id);
    weights_.push_back(weight);
  }
//JB 
   float getZPosition() const {
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
  

  const unsigned int getId() const { return index_; }
  std::vector<unsigned int> getNeighbours() const { return neighboursId_; }
  std::vector<float> getWeights() const { return weights_; }
  std::vector<unsigned int> getInner() const { return innerNodes_; }
  std::vector<unsigned int> getOuter() const { return outerNodes_; }
  std::vector<double> zDiffHistogram;
  double calculateSpatialDistance(const Node& node1, const Node& node2) {
    return std::abs(node1.getZPosition() - node2.getZPosition());
}

  void findSubComponents(std::vector<Node>& graph, std::vector<unsigned int>& subComponent, float threshold, unsigned int originNodeIndex, TH1F* h_zDiff) {
    //if (!alreadyVisited_) {
         if (alreadyVisited_) return;
	  alreadyVisited_ = true;
      subComponent.push_back(index_);

      std::cout << "Subcomponent after adding Node " << index_ << ": ";
    for (auto nodeId : subComponent) {
        std::cout << nodeId << " ";
    }
    std::cout << std::endl;
 
      //subComponent.insert(index_);
      auto neighbourIt = neighboursId_.cbegin(); // Iterator for neighboursId_
      auto weightIt = weights_.cbegin(); // Iterator for weights_
      // Loop through both vectors simultaneously
      std::vector<std::string> neighborDetails;
      for (; neighbourIt != neighboursId_.cend() && weightIt != weights_.cend(); ++neighbourIt, ++weightIt) {
          int neighbour = *neighbourIt;

          double weight = *weightIt;
          if (weight >= threshold){
          double spatialDistance = calculateSpatialDistance(graph[originNodeIndex], graph[neighbour]);
           const double spatialThreshold = 0.1;
	   std::cout<<" ## Original Index "<<originNodeIndex<<" neighbour "<<neighbour<<" spatialDistance "<<spatialDistance<<std::endl;
            // Check if the spatial distance is within the allowed range
           // if (spatialDistance < spatialThreshold) {

         float zPos = graph[neighbour].getZPosition();
	 float zPosCurrent = graph[index_].getZPosition();
         double zDiff = std::abs(zPosCurrent - zPos);
         h_zDiff->Fill(zDiff); 
        neighborDetails.push_back("Node " + std::to_string(neighbour) + " (z: " + std::to_string(zPos) + "), Weight: " + std::to_string(weight));		  
        
   	if (graph[index_].hasEdgeTo(neighbour)) {
		         // subComponent.insert(neighbour);
                    graph[neighbour].findSubComponents(graph, subComponent, threshold, originNodeIndex, h_zDiff);
              
	    //graph[0].findSubComponents(graph, subComponent, threshold) will start the DFS from node 0, identifying all connected nodes that meet the weight threshold and storing their indices in subComponent.
          //}
      }
    }
  }
      /*
if (!neighborDetails.empty()) {
    std::cout << "Neighbors of Node " << index_ << ": ";
    for (const auto& detail : neighborDetails) {
        std::cout << detail << ", ";
    }
    std::cout << std::endl;
}
*/
  }
 //}
    void setAlreadyVisited(bool visited) {
        alreadyVisited_ = visited;
    }
    bool hasAlreadyVisited() const {
        return alreadyVisited_;
    }
    bool hasEdgeTo(unsigned int otherId) const {
        return std::find(neighboursId_.begin(), neighboursId_.end(), otherId) != neighboursId_.end();
    }
  ~Node() = default;

private:
  unsigned index_;
  float zPosition_;    //JB    // Vector of features including z position  
  bool isTrack_;

  std::vector<unsigned int> neighboursId_;
  std::vector<float> weights_;
  std::vector<unsigned int> innerNodes_;
  std::vector<unsigned int> outerNodes_;
  bool alreadyVisited_;


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

    std::vector<std::vector<unsigned int>> findSubComponents(float threshold, TH1F* h_zDiff) {
    std::vector<std::vector<unsigned int>> components;
    for (auto& node : nodes_) {
            node.setAlreadyVisited(false);
        }
    
    for (auto const& node: nodes_) {
      auto const id = node.getId();
      std::vector<unsigned int> tmpSubComponents;
      
      //std::unordered_set<unsigned int> subComponentSet;
       if (!node.hasAlreadyVisited()) {  // Only start DFS if the node hasn't been visited
      nodes_[id].findSubComponents(nodes_, tmpSubComponents, threshold, id, h_zDiff);
      } 
      //std::vector<unsigned int> tmpSubComponents(subComponentSet.begin(), subComponentSet.end());

      if (!tmpSubComponents.empty()) {
        components.push_back(tmpSubComponents);
	//JB
       std::cout << "Subcomponent found: ";
            for (unsigned int nodeId : tmpSubComponents) {
                float zPos = nodes_[nodeId].getZPosition();  // Retrieve the z position
                std::cout << "Node " << nodeId << " (z: " << zPos << "), ";
            }
                        std::cout << std::endl;
//JB till here

      }
    }
    return components;
  }

   ~TrackGraph() = default;

private:
  std::vector<Node> nodes_;

};

#endif
