#ifndef RecoVertex_PrimaryVertexProducer_interface_GraphProducer_h
#define RecoVertex_PrimaryVertexProducer_interface_GraphProducer_h
#include "DataFormats/VertexReco/interface/MtdtimeHostCollection.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include <vector>
#include <map>
 
std::unique_ptr<TracksGraph> produce_tracks_graph(const std::vector<reco::TransientTrack>& transientTrack) {
   
    track_pt = transientTrack.track.pt();
    std::vector<Node> allNodes;
    reco::TransientTrack tt = trackBuilder->build(track, trkAssoc, time, dtime, mva, pathlength, btlchi2, btltimechi2, etlchi2, etltimechi2, time_pi, time_k, time_p, sigma_time_pi, sigma_time_k, sigma_time_p, field, tg) ;

    for() i



	    /// should get handles for transient track and builder 
	    //edm::ESHandle<TransientTrackBuilder> trackBuilder; i
            // reco::TransientTrack tt = trackBuilder->build(track, all the variables we need);
