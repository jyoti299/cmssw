#ifndef PFAnalyzer_H
#define PFAnalyzer_H

/** \class JetMETAnalyzer
 *
 *  DQM PF candidate analysis monitoring
 *
 *  \author J. Roloff - Brown University
 *
 */

#include <memory>
#include <fstream>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"


#include "DataFormats/JetReco/interface/Jet.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include <map>
#include <string>


class PFAnalyzer : public DQMEDAnalyzer {
public:
  /// Constructor
  PFAnalyzer(const edm::ParameterSet&);

  /// Destructor
  ~PFAnalyzer() override;

  /// Inizialize parameters for histo binning
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  /// Get the analysis
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  /// Initialize run-based parameters
  void dqmBeginRun(const edm::Run&, const edm::EventSetup&) override;

private:

  // A map between an observable name and a function that obtains that observable from a  PFCandidate.
  // This allows us to construct more complicated observables easily, and have it more configurable
  // in the config file.
  std::map<std::string, std::function<double (const reco::PFCandidate)> > m_funcMap;
  std::map<std::string, std::function<double (const reco::PFJet)> > m_jetFuncMap;

  // Book MonitorElements
  void bookMESetSelection(std::string, DQMStore::IBooker&);

  bool passesEventSelection(const edm::Event& iEvent);

  int getPFBin(const reco::PFCandidate pfCand);
  int getJetBin(const reco::PFJet jetCand);


  int getBinNumber(double binVal, std::vector<double> bins);
  int getBinNumbers(std::vector<double> binVal, std::vector<std::vector<double> > bins);
  std::vector<double> getBinList(std::string binString);

  std::vector<std::string> getAllSuffixes(std::vector<std::string> observables, std::vector<std::vector<double> > binnings);
  std::string stringWithDecimals(int bin, std::vector<double> bins);




  std::string getSuffix(std::vector<int> binList,  std::vector<std::string> observables, std::vector<std::vector<double> > binnings);


  // Various functions designed to get information from a PF canddidate
  static double getPt(const reco::PFCandidate pfCand){ return pfCand.pt();}
  static double getEta(const reco::PFCandidate pfCand){ return pfCand.eta();}
  static double getPhi(const reco::PFCandidate pfCand){ return pfCand.phi();}



  static double getTime(const reco::PFCandidate pfCand){ return pfCand.time();}

  static double getHcalEnergy_depth1(const reco::PFCandidate pfCand){ return pfCand.hcalDepthEnergyFraction(1);}
  static double getHcalEnergy_depth2(const reco::PFCandidate pfCand){ return pfCand.hcalDepthEnergyFraction(2);}
  static double getHcalEnergy_depth3(const reco::PFCandidate pfCand){ return pfCand.hcalDepthEnergyFraction(3);}
  static double getHcalEnergy_depth4(const reco::PFCandidate pfCand){ return pfCand.hcalDepthEnergyFraction(4);}
  static double getHcalEnergy_depth5(const reco::PFCandidate pfCand){ return pfCand.hcalDepthEnergyFraction(5);}
  static double getHcalEnergy_depth6(const reco::PFCandidate pfCand){ return pfCand.hcalDepthEnergyFraction(6);}
  static double getHcalEnergy_depth7(const reco::PFCandidate pfCand){ return pfCand.hcalDepthEnergyFraction(7);}

  static double getEcalEnergy(const reco::PFCandidate pfCand){ return pfCand.ecalEnergy();}
  static double getRawEcalEnergy(const reco::PFCandidate pfCand){ return pfCand.rawEcalEnergy();}
  static double getHcalEnergy(const reco::PFCandidate pfCand){ return pfCand.hcalEnergy();}
  static double getRawHcalEnergy(const reco::PFCandidate pfCand){ return pfCand.rawHcalEnergy();}
  static double getHOEnergy(const reco::PFCandidate pfCand){ return pfCand.hoEnergy();}
  static double getRawHOEnergy(const reco::PFCandidate pfCand){ return pfCand.rawHoEnergy();}

  static double getMVAIsolated(const reco::PFCandidate pfCand){ return pfCand.mva_Isolated();}
  static double getMVAEPi(const reco::PFCandidate pfCand){ return pfCand.mva_e_pi();}
  static double getMVAEMu(const reco::PFCandidate pfCand){ return pfCand.mva_e_mu();}
  static double getMVAPiMu(const reco::PFCandidate pfCand){ return pfCand.mva_pi_mu();}
  static double getMVANothingGamma(const reco::PFCandidate pfCand){ return pfCand.mva_nothing_gamma();}
  static double getMVANothingNH(const reco::PFCandidate pfCand){ return pfCand.mva_nothing_nh();}
  static double getMVAGammaNH(const reco::PFCandidate pfCand){ return pfCand.mva_gamma_nh();}
  
  static double getDNNESigIsolated(const reco::PFCandidate pfCand){ return pfCand.dnn_e_sigIsolated();}
  static double getDNNESigNonIsolated(const reco::PFCandidate pfCand){ return pfCand.dnn_e_sigNonIsolated();}
  static double getDNNEBkgNonIsolated(const reco::PFCandidate pfCand){ return pfCand.dnn_e_bkgNonIsolated();}
  static double getDNNEBkgTauIsolated(const reco::PFCandidate pfCand){ return pfCand.dnn_e_bkgTau();}
  static double getDNNEBkgPhotonIsolated(const reco::PFCandidate pfCand){ return pfCand.dnn_e_bkgPhoton();}

  static double getECalEFrac(const reco::PFCandidate pfCand){ return pfCand.ecalEnergy() / pfCand.energy();}
  static double getHCalEFrac(const reco::PFCandidate pfCand){ return pfCand.hcalEnergy() / pfCand.energy();}
  static double getTrackPt(const reco::PFCandidate pfCand){
                           if(pfCand.trackRef().isNonnull()) return (pfCand.trackRef())->pt();
                           return 0;
                          }

  static double getEoverP(const reco::PFCandidate pfCand){ 
    double energy = 0;
    int maxElement = pfCand.elementsInBlocks().size();
    for (int e = 0; e < maxElement; ++e) {
      // Get elements from block
      reco::PFBlockRef blockRef = pfCand.elementsInBlocks()[e].first;
      const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
      for (unsigned iEle = 0; iEle < elements.size(); iEle++) {
        if (elements[iEle].index() == pfCand.elementsInBlocks()[e].second) {
          if (elements[iEle].type() == reco::PFBlockElement::HCAL || elements[iEle].type() == reco::PFBlockElement::ECAL) {  // Element is HB or HE
            reco::PFClusterRef clusterref = elements[iEle].clusterRef();
            reco::PFCluster cluster = *clusterref;
            energy += cluster.energy();
          }
        }
      }
    }
    return energy / pfCand.p();
  }

  static double getHCalEnergy(const reco::PFCandidate pfCand){ 
    double energy = 0;
    int maxElement = pfCand.elementsInBlocks().size();
    for (int e = 0; e < maxElement; ++e) {
      // Get elements from block
      reco::PFBlockRef blockRef = pfCand.elementsInBlocks()[e].first;
      const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
      for (unsigned iEle = 0; iEle < elements.size(); iEle++) {
        if (elements[iEle].index() == pfCand.elementsInBlocks()[e].second) {
          if (elements[iEle].type() == reco::PFBlockElement::HCAL) {  // Element is HB or HE
            // Get cluster and hits
            reco::PFClusterRef clusterref = elements[iEle].clusterRef();
            reco::PFCluster cluster = *clusterref;
            //std::vector<std::pair<DetId, float>> hitsAndFracs = cluster.hitsAndFractions();
            energy += cluster.energy();
          }
        }
      }
    }
    return energy;
  }


  static double getECalEnergy(const reco::PFCandidate pfCand){
    double energy = 0;
    int maxElement = pfCand.elementsInBlocks().size();
    for (int e = 0; e < maxElement; ++e) {
      // Get elements from block
      reco::PFBlockRef blockRef = pfCand.elementsInBlocks()[e].first;
      const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
      for (unsigned iEle = 0; iEle < elements.size(); iEle++) {
        if (elements[iEle].index() == pfCand.elementsInBlocks()[e].second) {
          if (elements[iEle].type() == reco::PFBlockElement::ECAL) {  // Element is HB or HE
            // Get cluster and hits
            reco::PFClusterRef clusterref = elements[iEle].clusterRef();
            // TODO make sure this is correct -- I assume the cluster energy can't be shared in bad ways, but I'm not sure
            // When we don't have isolated tracks, this will be a bit useless, since the energy is shared across multiple tracks
            reco::PFCluster cluster = *clusterref;
            //std::vector<std::pair<DetId, float>> hitsAndFracs = cluster.hitsAndFractions();
            energy += cluster.energy();
          }
        }
      }
    }
    return energy;
  }

 
  static double getNTracksInBlock(const reco::PFCandidate pfCand){
    // We need this function to return a double, even though this is an integer value
    double nTrack = 0;
    int maxElement = pfCand.elementsInBlocks().size();
    for (int e = 0; e < maxElement; ++e) {
      // Get elements from block
      reco::PFBlockRef blockRef = pfCand.elementsInBlocks()[e].first;
      const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
      for (unsigned iEle = 0; iEle < elements.size(); iEle++) {
        if (elements[iEle].index() == pfCand.elementsInBlocks()[e].second) {
          if (elements[iEle].type() == reco::PFBlockElement::TRACK) {  // Element is HB or HE
            nTrack+=1;
          }
        }
      }
    }
    return nTrack;
  }




  static double getJetPt(const reco::PFJet jet){ return jet.pt();}

  edm::EDGetTokenT<reco::PFCandidateCollection> thePfCandidateCollection_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken_;
  edm::EDGetTokenT<reco::PFJetCollection> pfJetsToken_;


  std::vector<std::string> m_allSuffixes;
  std::vector<std::string> m_allJetSuffixes;

  // The directory where the output is stored
  std::string m_directory;

  // All of the histograms, stored as a map between the histogram name and the histogram
  std::map<std::string, MonitorElement*> map_of_MEs;

  //check later if we need only one set of parameters
  edm::ParameterSet parameters_;

  typedef std::vector<std::string> vstring;
  // Information on which observables to make histograms for.
  // In the config file, this should come as a comma-separated list of 
  // the observable name, the number of bins for the histogram, and
  // the lowest and highest values for the histogram.
  // The observable name should have an entry in m_funcMap to define how
  // it can be retrieved from a PFCandidate.
  vstring m_observables;
  vstring m_neutralPFObservables;
  vstring m_chargedPFObservables;

  // Information on what cuts should be applied to PFCandidates that are
  // being monitored. In the config file, this should come as a comma-separated list of 
  // the observable name, and the lowest and highest values for the histogram.
  // The observable name should have an entry in m_funcMap to define how
  // it can be retrieved from a PFCandidate.
  vstring m_cutList;
  //std::vector<double> m_cutMins;
  //std::vector<double> m_cutMaxes;
  std::vector<std::vector<double> > m_binList;

  // Information on what cuts should be applied to PFJets, in the case that we 
  // match PFCs to jets.In the config file, this should come as a comma-separated list of 
  // the observable name, and the lowest and highest values for the histogram.
  // The observable name should have an entry in m_jetFuncMap to define how
  // it can be retrieved from a PFJet.
  vstring m_jetCutList;
  //std::vector<double> m_jetCutMins;
  //std::vector<double> m_jetCutMaxes;
  std::vector<std::vector<double> > m_jetBinList;


  // The dR radius used to match PFCs to jets.
  // Making this configurable is useful in case you want to look at the core of a jet.
  double m_matchingRadius;


};
#endif
