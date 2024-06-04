#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/MtdtimeHostCollection.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/TrackReco/interface/TrackBase.h"

using namespace edm;

class MTDSoAProducer : public edm::stream::EDProducer<> {
public:
  MTDSoAProducer(const ParameterSet& pset);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event& ev, const edm::EventSetup& es) final;

private:
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> trackAssocToken_;
 // edm::EDGetTokenT<edm::ValueMap<float>> t0Token_;
 // edm::EDGetTokenT<edm::ValueMap<float>> sigmat0Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> betaToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> pathToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> MVAQualityToken_;
  edm::EDGetTokenT<edm::ValueMap<GlobalPoint>> posInMtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> momentumWithMTDToken_;
 /* edm::EDGetTokenT<edm::ValueMap<float>> probPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probPToken_;*/
  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchTimeChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> etlMatchChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> etlMatchTimeChi2Token_;  
  edm::EDGetTokenT<edm::ValueMap<float>> tofPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofpiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofkToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofpToken_;
};

MTDSoAProducer::MTDSoAProducer(const ParameterSet& iConfig)
    : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracksSrc"))),
      trackAssocToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackAssocSrc"))),
     /* t0Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"))),
      sigmat0Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"))),*/
      tmtdToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtdSrc"))),
      sigmatmtdToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatmtdSrc"))),
      betaToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("betamtd"))),
      pathToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pathmtd"))),
      MVAQualityToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("mvaquality"))),
      posInMtdToken_(consumes<edm::ValueMap<GlobalPoint>>(iConfig.getParameter<edm::InputTag>("posmtd"))),
      momentumWithMTDToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("momentum"))),
      /*probPiToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probPi"))),
      probKToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probK"))),
      probPToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probP"))),*/
      btlMatchChi2Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchChi2Src"))),
      btlMatchTimeChi2Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchTimeChi2Src"))),
      etlMatchChi2Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchChi2Src"))),
      etlMatchTimeChi2Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchTimeChi2Src"))),	
      tofPiToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofPi"))),
      tofKToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofK"))),
      tofPToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofP"))),	
      sigmatofpiToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofpiSrc"))),
      sigmatofkToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofkSrc"))),
      sigmatofpToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofpSrc"))){
  produces<MtdtimeHostCollection>();
}

// Configuration descriptions
void MTDSoAProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracksSrc", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"));
//  desc.add<edm::InputTag>("t0Src", edm::InputTag("tofPID:t0"));
//  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("tmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmatmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("betamtd", edm::InputTag("trackExtenderWithMTD:generalTrackBeta"));
  desc.add<edm::InputTag>("pathmtd", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("mvaquality", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<edm::InputTag>("posmtd", edm::InputTag("trackExtenderWithMTD:generalTrackmtdpos"));
  desc.add<edm::InputTag>("momentum", edm::InputTag("trackExtenderWithMTD:generalTrackp"));
 /* desc.add<edm::InputTag>("probPi", edm::InputTag("tofPID:probPi"));
  desc.add<edm::InputTag>("probK", edm::InputTag("tofPID:probK"));
  desc.add<edm::InputTag>("probP", edm::InputTag("tofPID:probP"));*/
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
  descriptions.add("mtdSoAProducer", desc);
}

void MTDSoAProducer::produce(edm::Event& ev, const edm::EventSetup& es) {
  edm::Handle<reco::TrackCollection> tracksH;
  ev.getByToken(tracksToken_, tracksH);
  const auto& tracks = *tracksH;

  const auto& trackAssoc = ev.get(trackAssocToken_);

  //const auto& t0 = ev.get(t0Token_);
  //const auto& sigmat0 = ev.get(sigmat0Token_);

  const auto& tmtd = ev.get(tmtdToken_);
  const auto& sigmatmtd = ev.get(sigmatmtdToken_);

  const auto& beta = ev.get(betaToken_);
  const auto& path = ev.get(pathToken_);
  const auto& MVAquality = ev.get(MVAQualityToken_);
  const auto& posInMTD = ev.get(posInMtdToken_);
  const auto& momentum = ev.get(momentumWithMTDToken_);
  /*const auto& probPi = ev.get(probPiToken_);
  const auto& probK = ev.get(probKToken_);
  const auto& probP = ev.get(probPToken_);*/
  const auto& btlchi2 = ev.get(btlMatchChi2Token_);
  const auto& btltimechi2 = ev.get(btlMatchTimeChi2Token_);	  
  const auto& etlchi2  = ev.get(etlMatchChi2Token_);
  const auto& etltimechi2 = ev.get(etlMatchTimeChi2Token_);
  const auto& tofPi = ev.get(tofPiToken_);
  const auto& tofK = ev.get(tofKToken_);
  const auto& tofP = ev.get(tofPToken_);
  const auto& sigmatofPi = ev.get(sigmatofpiToken_);
  const auto& sigmatofK = ev.get(sigmatofkToken_);
  const auto& sigmatofP = ev.get(sigmatofpToken_);


  auto MtdInfo = std::make_unique<MtdtimeHostCollection>(tracks.size(), cms::alpakatools::host());

  auto& MtdInfoView = MtdInfo->view();
  for (unsigned int iTrack = 0; iTrack < tracks.size(); ++iTrack) {
    const reco::TrackRef trackref(tracksH, iTrack);

    if (trackAssoc[trackref] == -1) {
      MtdInfoView.trackAsocMTD()[iTrack] = -1;
      MtdInfoView.IndxTrackAsocMTD() = 999999999; 
      //MtdInfoView.time0()[iTrack] = 0.f;
      //MtdInfoView.time0Err()[iTrack] = -1.f;
      MtdInfoView.time()[iTrack] = 0.f;
      MtdInfoView.timeErr()[iTrack] = -1.f;
      MtdInfoView.MVAquality()[iTrack] = 0.f;
      MtdInfoView.pathLength()[iTrack] = 0.f;
      MtdInfoView.beta()[iTrack] = 0.f;
      MtdInfoView.posInMTD_x()[iTrack] = 0.f;
      MtdInfoView.posInMTD_y()[iTrack] = 0.f;
      MtdInfoView.posInMTD_z()[iTrack] = 0.f;
      MtdInfoView.momentumWithMTD()[iTrack] = 0.f;
      /*MtdInfoView.probPi()[iTrack] = 0.f;
      MtdInfoView.probK()[iTrack] = 0.f;
      MtdInfoView.probP()[iTrack] = 0.f;*/
      MtdInfoView.btlMatch_chi2()[iTrack] = 0.f;
      MtdInfoView.btlMatchTime_chi2()[iTrack] = 0.f;	 
      MtdInfoView.etlMatch_chi2()[iTrack] = 0.f;
      MtdInfoView.etlMatchTime_chi2()[iTrack] = 0.f;   
      MtdInfoView.trackTime_pi()[iTrack] = 0.f;
      MtdInfoView.trackTime_k()[iTrack] = 0.f;
      MtdInfoView.trackTime_p()[iTrack] = 0.f;  
      MtdInfoView.track_sigmaTime_pi()[iTrack] = 0.f;
      MtdInfoView.track_sigmaTime_k()[iTrack] = 0.f;
      MtdInfoView.track_sigmaTime_p()[iTrack] = 0.f;
      continue;
    }

    MtdInfoView.trackAsocMTD()[iTrack] = trackAssoc[trackref];
    MtdInfoView.IndxTrackAsocMTD() = tracks.size();
    //MtdInfoView.time0()[iTrack] = t0[trackref];
    //MtdInfoView.time0Err()[iTrack] = sigmat0[trackref];
    MtdInfoView.time()[iTrack] = tmtd[trackref];
    MtdInfoView.timeErr()[iTrack] = sigmatmtd[trackref];
    MtdInfoView.MVAquality()[iTrack] = MVAquality[trackref];
    MtdInfoView.pathLength()[iTrack] = path[trackref];
    MtdInfoView.beta()[iTrack] = beta[trackref];
    MtdInfoView.posInMTD_x()[iTrack] = posInMTD[trackref].x();
    MtdInfoView.posInMTD_y()[iTrack] = posInMTD[trackref].y();
    MtdInfoView.posInMTD_z()[iTrack] = posInMTD[trackref].z();
    MtdInfoView.momentumWithMTD()[iTrack] = momentum[trackref];
    /*MtdInfoView.probPi()[iTrack] = probPi[trackref];
    MtdInfoView.probK()[iTrack] = probK[trackref];
    MtdInfoView.probP()[iTrack] = probP[trackref];*/
    MtdInfoView.btlMatch_chi2()[iTrack] = btlchi2[trackref];
    MtdInfoView.btlMatchTime_chi2()[iTrack] = btltimechi2[trackref];
    MtdInfoView.etlMatch_chi2()[iTrack] =  etlchi2[trackref];
    MtdInfoView.etlMatchTime_chi2()[iTrack] = etltimechi2[trackref];
    MtdInfoView.trackTime_pi()[iTrack] = tofPi[trackref];
    MtdInfoView.trackTime_k()[iTrack] = tofK[trackref];
    MtdInfoView.trackTime_p()[iTrack] = tofP[trackref];
    MtdInfoView.track_sigmaTime_pi()[iTrack] = sigmatofPi[trackref];
    MtdInfoView.track_sigmaTime_k()[iTrack] = sigmatofK[trackref];
    MtdInfoView.track_sigmaTime_p()[iTrack] = sigmatofP[trackref];
  }

  ev.put(std::move(MtdInfo));
}

//define this as a plug-in
#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(MTDSoAProducer);
