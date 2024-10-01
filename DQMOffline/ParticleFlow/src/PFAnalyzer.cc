/** \class PFAnalyzer
 *
 *  DQM ParticleFlow analysis monitoring
 *
 *  \author J. Roloff - Brown University
 *
 */

#include "DQMOffline/ParticleFlow/interface/PFAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"

#include <string>
#include <cmath>

// ***********************************************************
PFAnalyzer::PFAnalyzer(const edm::ParameterSet& pSet){
  m_directory = "ParticleFlow";
  parameters_ = pSet.getParameter<edm::ParameterSet>("pfAnalysis");

  thePfCandidateCollection_ = consumes<reco::PFCandidateCollection>(pSet.getParameter<edm::InputTag>("pfCandidates"));
  pfJetsToken_ = consumes<reco::PFJetCollection>(pSet.getParameter<edm::InputTag>("pfJetCollection"));

  // TODO we might want to define some observables that are only for PFCs in jets,
  // and that might take a jet as an input to the function as well
  m_observables = parameters_.getParameter<vstring>("observables");

  // List of cuts applied to PFCs that we want to plot
  m_cutList = parameters_.getParameter<vstring>("cutList");
  // List of jet cuts that we apply for the case of plotting PFCs in jets
  m_jetCutList = parameters_.getParameter<vstring>("jetCutList");

  // Link observable strings to the static functions defined in the header file
  // Many of these are quite trivial, but this enables a simple way to include a 
  // variety of observables on-the-fly. 
  m_funcMap["pt"] = getPt;
  m_funcMap["eta"] = getEta;
  m_funcMap["phi"] = getPhi;



  m_funcMap["HCalE_depth1"] = getHcalEnergy_depth1;
  m_funcMap["HCalE_depth2"] = getHcalEnergy_depth2;
  m_funcMap["HCalE_depth3"] = getHcalEnergy_depth3;
  m_funcMap["HCalE_depth4"] = getHcalEnergy_depth4;
  m_funcMap["HCalE_depth5"] = getHcalEnergy_depth5;
  m_funcMap["HCalE_depth6"] = getHcalEnergy_depth6;
  m_funcMap["HCalE_depth7"] = getHcalEnergy_depth7;

  m_funcMap["ECal_E"] = getEcalEnergy;
  m_funcMap["RawECal_E"] = getRawEcalEnergy;
  m_funcMap["HCal_E"] = getHcalEnergy;
  m_funcMap["RawHCal_E"] = getRawHcalEnergy;
  m_funcMap["HO_E"] = getHOEnergy;
  m_funcMap["RawHO_E"] = getRawHOEnergy;

  m_funcMap["MVAIsolated"] = getMVAIsolated;
  m_funcMap["MVAEPi"] = getMVAEPi;
  m_funcMap["MVAEMu"] = getMVAEMu;
  m_funcMap["MVAPiMu"] = getMVAPiMu;
  m_funcMap["MVANothingGamma"] = getMVANothingGamma;
  m_funcMap["MVANothingNH"] = getMVANothingNH;
  m_funcMap["MVAGammaNH"] = getMVAGammaNH;

  m_funcMap["DNNESigIsolated"] = getDNNESigIsolated;
  m_funcMap["DNNESigNonIsolated"] = getDNNESigNonIsolated;
  m_funcMap["DNNEBkgNonIsolated"] = getDNNEBkgNonIsolated;
  m_funcMap["DNNEBkgTauIsolated"] = getDNNEBkgTauIsolated;
  m_funcMap["DNNEBkgPhotonIsolated"] = getDNNEBkgPhotonIsolated;

  m_funcMap["hcalE"] = getHCalEnergy;
  m_funcMap["eOverP"] = getEoverP;
  m_funcMap["nTrkInBlock"] = getNTracksInBlock;

  // Link jet observables to static functions in the header file.
  // This is very similar to m_funcMap, but for jets instead.
  m_jetFuncMap["pt"] = getJetPt;

  // Convert the cutList strings into real cuts that can be applied
  // The format should be three comma separated values
  // with the first number being the name of the observable
  // (corresponding to a key in m_funcMap),
  // the second being the minimum value, and the last being the max.
  for(unsigned int i=0; i<m_cutList.size(); i++){
    size_t pos = m_cutList[i].find(";");
    std::string observableName = m_cutList[i].substr(0, pos);
    m_cutList[i].erase(0, pos + 1);

    m_binList.push_back(getBinList(m_cutList[i]));
    m_cutList[i] = observableName;
  }

  // Convert the jetCutList strings into real cuts that can be applied
  // The format should be three comma separated values,
  // with the first number being the name of the observable
  // (corresponding to a key in m_jetFuncMap),
  // the second being the minimum value, and the last being the max value.
  for(unsigned int i=0; i<m_jetCutList.size(); i++){
    size_t pos = m_jetCutList[i].find(";");
    std::string observableName = m_jetCutList[i].substr(0, pos);
    m_jetCutList[i].erase(0, pos + 1);

    m_jetBinList.push_back(getBinList(m_jetCutList[i]));
    m_jetCutList[i] = observableName;
  }
}

// ***********************************************************
PFAnalyzer::~PFAnalyzer() {
  LogTrace("PFAnalyzer") << "[PFAnalyzer] Saving the histos";
}

// ***********************************************************
void PFAnalyzer::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& iRun, edm::EventSetup const&) {
  ibooker.setCurrentFolder(m_directory);

  m_allSuffixes = getAllSuffixes(m_cutList, m_binList);
  m_allJetSuffixes = getAllSuffixes(m_jetCutList, m_jetBinList);
  


  // TODO: Make it possible to use an arbitrary list of bins instead of evenly space bins?
  //
  // Books a histogram for each histogram in the config file.
  // The format for the observables should be four comma separated values,
  // with the first being the observable name (corresponding to one of 
  // the keys in m_funcMap), the second being the number of bins,
  // and the last two being the min and max value for the histogram respectively.
  for(unsigned int i=0; i<m_observables.size(); i++){
    size_t pos = m_observables[i].find(";");
    std::string observableName = m_observables[i].substr(0, pos);
    m_observables[i].erase(0, pos + 1);
    
    pos = m_observables[i].find(";");
    std::string axisString = m_observables[i].substr(0, pos);
    m_observables[i].erase(0, pos + 1);

    pos = m_observables[i].find(";");
    int nBins = atoi(m_observables[i].substr(0, pos).c_str());
    m_observables[i].erase(0, pos + 1);
    
    pos = m_observables[i].find(";");
    float binMin = atof(m_observables[i].substr(0, pos).c_str());
    m_observables[i].erase(0, pos + 1);
    
    float binMax = atof(m_observables[i].c_str());
    m_observables[i] = observableName;
    
    for(unsigned int j=0; j<m_allSuffixes.size(); j++){


      // For each observable, we make a couple histograms based on a few generic categorizations.
      // In all cases, the PFCs that go into these histograms must pass the PFC selection from m_cutList.

      // All PFCs
      std::string histName = Form("allPFC_%s%s", observableName.c_str(), m_allSuffixes[j].c_str());
      MonitorElement* mHist = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
      map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHist));

      // Neutral hadron PFCs
      histName = Form("neutralHadPFC_%s%s", observableName.c_str(), m_allSuffixes[j].c_str());
      MonitorElement* mHistNeutral = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
      map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistNeutral));

      // Charged hadron PFCs
      histName = Form("chargedHadPFC_%s%s", observableName.c_str(), m_allSuffixes[j].c_str());
      MonitorElement* mHistCharged = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
      map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistCharged));

      // Electron PFCs
      histName = Form("electronPFC_%s%s", observableName.c_str(), m_allSuffixes[j].c_str());
      MonitorElement* mHistElectron = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
      map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistElectron));

      // Muon PFCs
      histName = Form("muonPFC_%s%s", observableName.c_str(), m_allSuffixes[j].c_str());
      MonitorElement* mHistMuon = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
      map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistMuon));

      // Gamma PFCs
      histName = Form("gammaPFC_%s%s", observableName.c_str(), m_allSuffixes[j].c_str());
      MonitorElement* mHistGamma = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
      map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistGamma));
  
      // Had HF PFCs
      histName = Form("hadHFPFC_%s%s", observableName.c_str(), m_allSuffixes[j].c_str());
      MonitorElement* mHistHadHF = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
      map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistHadHF));
  
      // EM HF PFCs
      histName = Form("emHFPFC_%s%s", observableName.c_str(), m_allSuffixes[j].c_str());
      MonitorElement* mHistEMHF = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
      map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistEMHF));



      for(unsigned int k=0; k<m_allJetSuffixes.size(); k++){
        // These histograms are for PFCs passing the basic selection, and which are matched to jets
        // that pass the jet selection
        histName = Form("allPFC_jetMatched_%s%s_jetCuts%s", observableName.c_str(), m_allSuffixes[j].c_str(), m_allJetSuffixes[k].c_str());
        MonitorElement* mHistInJet = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
        map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistInJet));
  
        // This histogram is all neutral PFCs that pass the basic selection
        histName = Form("neutralHadPFC_jetMatched_%s%s_jetCuts%s", observableName.c_str(), m_allSuffixes[j].c_str(), m_allJetSuffixes[k].c_str());
        MonitorElement* mHistNeutralInJet = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
        map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistNeutralInJet));
  
        // This histogram is all neutral PFCs that pass the basic selection
        histName = Form("chargedHadPFC_jetMatched_%s%s_jetCuts%s", observableName.c_str(), m_allSuffixes[j].c_str(), m_allJetSuffixes[k].c_str());
        MonitorElement* mHistChargedInJet = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
        map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistChargedInJet));
  
        // This histogram is electron PFCs that pass the basic selection
        histName = Form("electronPFC_jetMatched_%s%s_jetCuts%s", observableName.c_str(), m_allSuffixes[j].c_str(), m_allJetSuffixes[k].c_str());
        MonitorElement* mHistElectronInJet = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
        map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistElectronInJet));
  
        // This histogram is electron PFCs that pass the basic selection
        histName = Form("muonPFC_jetMatched_%s%s_jetCuts%s", observableName.c_str(), m_allSuffixes[j].c_str(), m_allJetSuffixes[k].c_str());
        MonitorElement* mHistMuonInJet = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
        map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistMuonInJet));
  
        // This histogram is electron PFCs that pass the basic selection
        histName = Form("gammaPFC_jetMatched_%s%s_jetCuts%s", observableName.c_str(), m_allSuffixes[j].c_str(), m_allJetSuffixes[k].c_str());
        MonitorElement* mHistGammaInJet = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
        map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistGammaInJet));
  
        // This histogram is electron PFCs that pass the basic selection
        histName = Form("hadHFPFC_jetMatched_%s%s_jetCuts%s", observableName.c_str(), m_allSuffixes[j].c_str(), m_allJetSuffixes[k].c_str());
        MonitorElement* mHistHadHFInJet = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
        map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistHadHFInJet));
        
        // This histogram is electron PFCs that pass the basic selection
        histName = Form("emHFPFC_jetMatched_%s%s_jetCuts%s", observableName.c_str(), m_allSuffixes[j].c_str(), m_allJetSuffixes[k].c_str());
        MonitorElement* mHistEMHFInJet = ibooker.book1D(histName, Form(";%s;", axisString.c_str()), nBins, binMin, binMax);
        map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistEMHFInJet));
      }
    }    
  }
}

void PFAnalyzer::bookMESetSelection(std::string DirName, DQMStore::IBooker& ibooker) {
  ibooker.setCurrentFolder(DirName);
}

// ***********************************************************
void PFAnalyzer::dqmBeginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
}


bool PFAnalyzer::passesEventSelection(const edm::Event& iEvent){
  return true;
}


// How many significant digits do we need to save for the values to be distinct?
std::string PFAnalyzer::stringWithDecimals(int bin, std::vector<double> bins){
  double diff = bins[bin+1] - bins[bin];
  double sigFigs = log10(diff);

  // We only want to save as many significant digits as we need to.
  // Currently, we might lose some information, so we should think about 
  // if we want to specify more digits
  if(sigFigs >=1){
    return Form("%.0f_%.0f", bins[bin], bins[bin+1]);
  }

  int nDecimals = int(-1*sigFigs) + 1;
  // We do not want to use decimals since these can mess up histogram retrieval in some cases.
  // Instead, we use a '_' to indicate the decimal.
  double newDigit =  (bins[bin] -  int(bins[bin]))*pow(10, nDecimals);
  double newDigit2 =  (bins[bin+1] -  int(bins[bin+1]))*pow(10,nDecimals);
  std::string signStringLow = "";
  std::string signStringHigh = "";
  if(bins[bin]<0) signStringLow = "m";
  if(bins[bin+1]<0) signStringHigh = "m";
  return Form("%s%.0fp%.0f__%s%.0fp%.0f", signStringLow.c_str(), bins[bin], newDigit, signStringHigh.c_str(), bins[bin+1], newDigit2);
  

}


std::vector<double> PFAnalyzer::getBinList(std::string binString){
  std::vector<double> binList;

  while(binString.find(";")!= std::string::npos){
    size_t pos = binString.find(";");
    binList.push_back( atof(binString.substr(0, pos).c_str()));
    binString.erase(0, pos + 1);
  }
  binList.push_back( atof(binString.c_str()));
  //std::cout << binList.size() << std::endl;

  if(binList.size()==3){
    int nBins = int(binList[0]);
    double minVal = binList[1];
    double maxVal = binList[2];
    binList.clear();
    
    for(int i=0; i<=nBins; i++){
      binList.push_back( minVal + i* (maxVal - minVal) / nBins);
    }
  }


  return binList;
}


std::vector<std::string> PFAnalyzer::getAllSuffixes(std::vector<std::string> observables, std::vector<std::vector<double> > binnings){
  int nTotalBins = 1;
  std::vector<int> nBins;
  for(unsigned int i=0; i<binnings.size(); i++){
    nTotalBins = (binnings[i].size()-1)*nTotalBins;
    nBins.push_back(binnings[i].size()-1);
    std::cout << (binnings[i].size()-1) << "\t" ;
  }
  std::cout << nTotalBins << std::endl;
  

  std::vector<std::vector<int> > binList;
  
  for(int i=0; i<nTotalBins; i++){
    binList.push_back(std::vector<int>());
  }

  int factor = nTotalBins;
  int otherFactor = 1;
  for(unsigned int i=0; i<binnings.size(); i++){
    factor = factor / nBins[i];
    
    for(int j=0; j<nBins[i]; j++){
      for(int k=0; k<factor; k++){
        for(int m=0; m<otherFactor; m++){
          binList[m*otherFactor + j*factor + k].push_back(j);
        }
      }
    }
    otherFactor = otherFactor * nBins[i];
  }


  std::vector<std::string> allSuffixes;
  for(int i=0; i<nTotalBins; i++){
    allSuffixes.push_back(getSuffix(binList[i], observables, binnings));
  }

  return allSuffixes;
}

// Get a unique string corresponding to the selection cuts
std::string PFAnalyzer::getSuffix(std::vector<int> binList, std::vector<std::string> observables, std::vector<std::vector<double> > binnings){
  // TODO: we might not want to add a suffix if something is just a cut (e.g. just one bin)
  std::string suffix = "";
  for(unsigned int i=0; i<binList.size(); i++){
    // TODO we need to determine how many decimal places to save, and save this in a root-compatible way
    // (using '.' does not work nicely in root always)
    if(binList[i] < 0) return "";
    std::string digitString = stringWithDecimals(binList[i], binnings[i]);

    //suffix = Form("%s_%s_%.0f_%.0f", suffix.c_str(), observables[i].c_str(), binnings[i][binList[i]], binnings[i][binList[i]+1]);
    suffix = Form("%s_%s_%s", suffix.c_str(), observables[i].c_str(), digitString.c_str());
  }

  return suffix;
}

int PFAnalyzer::getBinNumber(double binVal, std::vector<double> bins){
  if(binVal < bins[0]) return -1;
  for(unsigned int i=0; i<bins.size(); i++){
    if(binVal<bins[i]) return i-1;
  }

  return -1;
}

int PFAnalyzer::getBinNumbers(std::vector<double> binVal, std::vector<std::vector<double> > bins){
  std::vector<int> cbins;
  std::vector<int> nBins;
  for(unsigned int i=0; i<binVal.size(); i++){
    int cbin = getBinNumber(binVal[i], bins[i]);
    if(cbin <0) return -1;
    nBins.push_back(bins[i].size()-1);
    cbins.push_back(cbin);
  }

  int bin = 0;
  int factor = 1;
  for(unsigned int i=0; i<binVal.size(); i++){
    bin += cbins[i]*factor;
    factor = factor*nBins[i];
  }
  //bin += 1;

  return bin;
}

int PFAnalyzer::getPFBin(const reco::PFCandidate pfCand){
  std::vector<double> binVals;
  for(unsigned int i=0; i<m_cutList.size(); i++){
    binVals.push_back(m_funcMap[m_cutList[i]](pfCand));
  }

  return getBinNumbers(binVals, m_binList);
}

int PFAnalyzer::getJetBin(const reco::PFJet jetCand){
  std::vector<double> binVals;
  for(unsigned int i=0; i<m_jetCutList.size(); i++){
    binVals.push_back(m_jetFuncMap[m_jetCutList[i]](jetCand));
  }
  
  return getBinNumbers(binVals, m_jetBinList);
}


// ***********************************************************
void PFAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Retrieve the PFCs
  edm::Handle<reco::PFCandidateCollection> pfCollection;
  iEvent.getByToken(thePfCandidateCollection_, pfCollection);
  if (!pfCollection.isValid()) {
    edm::LogError("PFAnalyzer") << "invalid collection: PF candidate \n";
    return;
  }

  edm::Handle<reco::PFJetCollection> pfJets;
  iEvent.getByToken(pfJetsToken_, pfJets);
  if (!pfJets.isValid()) {
    edm::LogError("PFAnalyzer") << "invalid collection: PF jets \n";
    return;
  }

  // TODO Perform some event selection
  // Probably we want to define a few different options for how the selection will work
  // Currently it is just a dummy function
  if(!passesEventSelection(iEvent)) return;


  for (reco::PFCandidateCollection::const_iterator recoPF = pfCollection->begin(); recoPF != pfCollection->end(); ++recoPF) {
    int binNumber = getPFBin(*recoPF);
    if(binNumber<0) continue;
    if(binNumber >= int(m_allSuffixes.size())) {
      continue;
    }
    std::string binString = m_allSuffixes[binNumber];

    // Eventually, we might want the hist name to include the cuts that we are applying, 
    // so I am keepking it as a separate string for now, even though it is redundant.
    // Make plots of all observables
    for(unsigned int i=0; i<m_observables.size(); i++){
      std::string histName = Form("%s%s", m_observables[i].c_str(), binString.c_str());
      map_of_MEs[m_directory + "/allPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));

      switch(recoPF->particleId()){
        case reco::PFCandidate::ParticleType::h :
          map_of_MEs[m_directory + "/chargedHadPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        case reco::PFCandidate::ParticleType::h0 :
          map_of_MEs[m_directory + "/neutralHadPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        case reco::PFCandidate::ParticleType::e :
          map_of_MEs[m_directory + "/electronPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        case reco::PFCandidate::ParticleType::mu :
          map_of_MEs[m_directory + "/muonPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        case reco::PFCandidate::ParticleType::gamma :
          map_of_MEs[m_directory + "/gammaPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        case reco::PFCandidate::ParticleType::h_HF :
          map_of_MEs[m_directory + "/hadHFPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        case reco::PFCandidate::ParticleType::egamma_HF :
          map_of_MEs[m_directory + "/emHFPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        default:
           break;
      }
    }
  }

  std::cout << std::endl;

  // Make plots of all observables, this time for PF candidates within jets
  for (reco::PFJetCollection::const_iterator cjet = pfJets->begin(); cjet != pfJets->end(); ++cjet) {

    int jetBinNumber = getJetBin(*cjet);
    if(jetBinNumber<0) continue;
    std::string jetBinString = m_allJetSuffixes[jetBinNumber];

    std::vector<reco::PFCandidatePtr> pfConstits = cjet->getPFConstituents();

    // TODO the binning doesn't actually work right yet
    for(auto recoPF:pfConstits){
      int binNumber = getPFBin(*recoPF);
      if(binNumber<0) continue;
      if(binNumber >= int(m_allSuffixes.size())) {
        continue;
      }
      std::string binString = m_allSuffixes[binNumber];

      for(unsigned int i=0; i<m_observables.size(); i++){
        std::string histName = Form("%s%s_jetCuts%s", m_observables[i].c_str(), binString.c_str(), jetBinString.c_str());
        map_of_MEs[m_directory + "/allPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));

        switch(recoPF->particleId()){
          case reco::PFCandidate::ParticleType::h :
            map_of_MEs[m_directory + "/chargedHadPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
            break;
          case reco::PFCandidate::ParticleType::h0 :
            map_of_MEs[m_directory + "/neutralHadPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
             break;
          case reco::PFCandidate::ParticleType::e :
            map_of_MEs[m_directory + "/electronPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
            break;
          case reco::PFCandidate::ParticleType::mu :
            map_of_MEs[m_directory + "/muonPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
            break;
          case reco::PFCandidate::ParticleType::gamma :
            map_of_MEs[m_directory + "/gammaPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
            break;
          case reco::PFCandidate::ParticleType::h_HF :
            map_of_MEs[m_directory + "/hadHFPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
            break;
          case reco::PFCandidate::ParticleType::egamma_HF :
            map_of_MEs[m_directory + "/emHFPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
            break;
          default:
             break;
        }
      }
    }
  }

}
