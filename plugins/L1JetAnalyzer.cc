/// -*- C++ -*-
//
// Package:    L1Trigger/L1JetAnalyzer
// Class:      L1JetAnalyzer
//
/**\class L1JetAnalyzer L1JetAnalyzer.cc L1Trigger/L1JetAnalyzer/plugins/L1JetAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Elliott Kauffman
//         Created:  Tue, 13 Jun 2023 09:23:09 GMT
//
//

// system include files
#include <memory>
#include <iostream>
#include <fstream>

// user include files
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "JetMETCorrections/Type1MET/interface/PFJetMETcorrInputProducerT.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "PhysicsTools/PatUtils/interface/PATJetCorrExtractor.h"


#include "TLorentzVector.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

using namespace l1extra;
using namespace std;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


class L1JetAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit L1JetAnalyzer(const edm::ParameterSet&);
  ~L1JetAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  const edm::EDGetTokenT<l1t::JetBxCollection> jetBXCollectionToken_; // l1 jets
  edm::EDGetTokenT< BXVector<GlobalAlgBlk> > gtAlgBlkToken;
  edm::Handle< BXVector<GlobalAlgBlk> > gtAlgBlkHandle;
  edm::EDGetTokenT<vector<pat::Jet>> slimmedJetsToken_; // reco jets to match to l1 jets
  edm::EDGetTokenT<vector<pat::Muon>> slimmedMuonsToken_; // needed for HLT_IsoMu20
  edm::EDGetTokenT<vector<pat::MET>> slimmedMETToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsORIGToken_; // need to require HLT_IsoMu20 in order to match reco jets
  const edm::ESGetToken<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd> l1GtMenuToken_;

  edm::ESGetToken<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd> L1TUtmTriggerMenuEventToken;
  const L1TUtmTriggerMenu* l1GtMenu;
  const std::map<std::string, L1TUtmAlgorithm>* algorithmMap;

  // triggers
  bool L1_SingleJet180;
  bool L1_HTT280er;
  bool L1_ETMHF90;
  int idx_L1_SingleJet180;
  int idx_L1_HTT280er;
  int idx_L1_ETMHF90;

  // histograms
  TH1F *h_SingleJet180_den;
  TH1F *h_SingleJet180_num;
  TH1F *h_SingleJet180_L130_num;
  TH1F *h_SingleJet180_L160_num;
  TH1F *h_SingleJet180_L1120_num;
  TH1F *h_SingleJet180_L1180_num;
  TH1F *h_SingleJet180_central_den;
  TH1F *h_SingleJet180_central_num;
  TH1F *h_SingleJet180_forward_den;
  TH1F *h_SingleJet180_forward_num;
  TH1F *h_L1_HTT280er_den;
  TH1F *h_L1_HTT280er_num;
  TH1F *h_ETMHF90_den;
  TH1F *h_ETMHF90_num;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constructors and destructor
//
L1JetAnalyzer::L1JetAnalyzer(const edm::ParameterSet& iConfig):
  jetBXCollectionToken_(consumes<l1t::JetBxCollection>(edm::InputTag("caloStage2Digis","Jet","RECO"))),
  gtAlgBlkToken( consumes< BXVector<GlobalAlgBlk> >(edm::InputTag("gtStage2Digis","","RECO")) ),
  l1GtMenuToken_(esConsumes<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd>())
{
  
  // tokens
  L1TUtmTriggerMenuEventToken = consumesCollector().esConsumes<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd>();
  slimmedJetsToken_ = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJetsPuppi") );
  slimmedMuonsToken_ = consumes< std::vector<pat::Muon> >(edm::InputTag("slimmedMuons") );
  slimmedMETToken_ = consumes< std::vector<pat::MET> >(edm::InputTag("slimmedMETsPuppi") );

  trgresultsORIGToken_ = consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT") );

  edm::Service<TFileService> fs;
  
  h_SingleJet180_den= fs->make<TH1F>("h_SingleJet180_den","",40,0,400);
  h_SingleJet180_num= fs->make<TH1F>("h_SingleJet180_num","",40,0,400);
  h_SingleJet180_L130_num = fs->make<TH1F>("h_SingleJet180_L130_num","",40,0,400);
  h_SingleJet180_L160_num = fs->make<TH1F>("h_SingleJet180_L160_num","",40,0,400);
  h_SingleJet180_L1120_num = fs->make<TH1F>("h_SingleJet180_L1120_num","",40,0,400);
  h_SingleJet180_L1180_num = fs->make<TH1F>("h_SingleJet180_L1180_num","",40,0,400);
  h_SingleJet180_central_den= fs->make<TH1F>("h_SingleJet180_central_den","",40,0,400);
  h_SingleJet180_central_num= fs->make<TH1F>("h_SingleJet180_central_num","",40,0,400);
  h_SingleJet180_forward_den= fs->make<TH1F>("h_SingleJet180_forward_den","",40,0,400);
  h_SingleJet180_forward_num= fs->make<TH1F>("h_SingleJet180_forward_num","",40,0,400);
  h_L1_HTT280er_den= fs->make<TH1F>("h_L1_HTT280er_den","",40,0,600);
  h_L1_HTT280er_num= fs->make<TH1F>("h_L1_HTT280er_num","",40,0,600);
  h_ETMHF90_den= fs->make<TH1F>("h_ETMHF90_den","",40,0,400);
  h_ETMHF90_num= fs->make<TH1F>("h_ETMHF90_num","",40,0,400);

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

L1JetAnalyzer::~L1JetAnalyzer() {
}

//
// member functions
//

// AKpuppi jetID
template <typename RecoJet>
bool puppiJetID(const RecoJet& jet) {
  bool tmp = true;
  if (std::abs(jet.eta())<2.6) {
    tmp &= jet.neutralHadronEnergyFraction() < 0.9;
    tmp &= jet.muonEnergyFraction() < 0.8;
    tmp &= jet.chargedEmEnergyFraction() < 0.8;
    tmp &= jet.chargedMultiplicity() > 0;
    tmp &= jet.chargedHadronEnergyFraction() > 0.01;
    tmp &= (jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1;
    tmp &= jet.neutralEmEnergyFraction() < 0.9;
  }
  if (std::abs(jet.eta())>2.6 && std::abs(jet.eta())<=2.7){
    tmp &= jet.chargedEmEnergyFraction() < 0.8;
    tmp &= jet.neutralEmEnergyFraction() < 0.99;
    tmp &= jet.muonEnergyFraction() < 0.8;
    tmp &= jet.neutralHadronEnergyFraction() < 0.9;
  }
  if (std::abs(jet.eta())>2.7 && std::abs(jet.eta())<= 3.0){
    tmp &= jet.neutralEmEnergyFraction() < 0.99;
  }
  if (std::abs(jet.eta())>3.0){
    tmp &= jet.neutralEmEnergyFraction() < 0.9;
    tmp &= jet.neutralMultiplicity() > 2;
  }
  return tmp;
}

// deltaR method (written using https://github.com/cms-sw/cmssw/blob/master/DataFormats/Math/interface/deltaR.h)
float deltaR(float p1, float p2, float e1, float e2) {

  auto dp = std::abs(p1 - p2);
  if (dp > M_PI) dp -= 2 * M_PI;
  return sqrt((e1 - e2) * (e1 - e2) + dp * dp);
}


template <typename L1JetCollection>
bool checkMatchBX(float recoPhi, float recoEta, const L1JetCollection& l1jetcollection, int bx, TLorentzVector* l1_jet ) {

  bool match = false;
  float temp = 99.0;
  float matchedPt = 0.0;
  float matchedEta = 0.0;
  float matchedPhi = 0.0;
  float matchedE = 0.0;

  for (auto it = l1jetcollection.begin(bx); it!=l1jetcollection.end(bx); it++){
    // check if deltaR is lower than temp
    if(deltaR(it->phi(), recoPhi, it->eta(), recoEta)<temp) {
      temp = deltaR(it->phi(), recoPhi, it->eta(), recoEta);
      matchedPt = it->pt();
      matchedEta = it->eta();
      matchedPhi = it->phi();
      matchedE = it->energy();
    }
  }

  if(temp < 0.4) {
    match = true;
    l1_jet->SetPtEtaPhiE(matchedPt, matchedEta, matchedPhi, matchedE);
  }

  return match;
}


// ------------ method called for each event  ------------
void L1JetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::ESHandle<L1TUtmTriggerMenu> menu;
  menu = iSetup.getHandle(l1GtMenuToken_);

  iEvent.getByToken(gtAlgBlkToken, gtAlgBlkHandle);

  // get index of L1 triggers
  int idx_L1_SingleJet180 = -1;
  int idx_L1_HTT280er = -1;
  int idx_L1_ETMHF90 = -1;
  for (auto const &keyval : menu->getAlgorithmMap()) {
    
    std::string const &name = keyval.second.getName();
    //std::cout<<name<<"    ";
    unsigned int indx = keyval.second.getIndex();
    //std::cout<<indx<<std::endl;
    if(name.find("L1_SingleJet180")!=string::npos) {
      idx_L1_SingleJet180 = indx;
      //std::cout<<name<<std::endl;
    }
    if(name.find("L1_HTT280er")!=string::npos) {
      idx_L1_HTT280er = indx;
      //std::cout<<name<<std::endl;
    }
    if(name.find("L1_ETMHF90")!=string::npos) {
      idx_L1_ETMHF90 = indx;
      //std::cout<<name<<std::endl;
    }
  }
  
  // get l1 trigger results
  L1_SingleJet180 = gtAlgBlkHandle->begin(0)->getAlgoDecisionFinal(idx_L1_SingleJet180);
  L1_HTT280er = gtAlgBlkHandle->begin(0)->getAlgoDecisionFinal(idx_L1_HTT280er);
  L1_ETMHF90 = gtAlgBlkHandle->begin(0)->getAlgoDecisionFinal(idx_L1_ETMHF90);
  std::cout << "L1_SingleJet180 = " << L1_SingleJet180 << std::endl;
  std::cout << "L1_HTT280er = " <<L1_HTT280er << std::endl;
  std::cout << "L1_ETMHF90 = " << L1_ETMHF90<< std::endl;

  // get reco jets and muons
  edm::Handle< std::vector<pat::Jet> > slimmedJets;
  iEvent.getByToken(slimmedJetsToken_,slimmedJets ); 
  edm::Handle< std::vector<pat::Muon> > slimmedMuons;
  iEvent.getByToken(slimmedMuonsToken_,slimmedMuons );

  // get met collection
  edm::Handle< std::vector<pat::MET> > slimmedMET;
  iEvent.getByToken(slimmedMETToken_, slimmedMET );


  /*auto menuRcd = iSetup.get<L1TUtmTriggerMenuRcd>();
  l1GtMenu = &menuRcd.get(L1TUtmTriggerMenuEventToken);
  algorithmMap = &(l1GtMenu->getAlgorithmMap());
  if(gtAlgBlkHandle.isValid()){
    std::vector<GlobalAlgBlk>::const_iterator algBlk = gtAlgBlkHandle->begin(0);
    if(algBlk != gtAlgBlkHandle->end(0)){
      for (std::map<std::string, L1TUtmAlgorithm>::const_iterator itAlgo = algorithmMap->begin(); itAlgo != algorithmMap->end(); itAlgo++) {
        std::string algName = itAlgo->first;
        int algBit = itAlgo->second.getIndex();
        std::cout<<algName<<"\t"<<algBit<<std::endl;
      }
    }
  }*/

  //get HLT_IsoMu20 result
  bool HLT_IsoMu20(false); 
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsORIGToken_, trigResults);
  if( !trigResults.failedToGet() ) {
    int N_Triggers = trigResults->size();
    const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      if (trigResults.product()->accept(i_Trig)) {
        TString TrigPath =trigName.triggerName(i_Trig);
        if(TrigPath.Index("HLT_IsoMu20_v") >=0) HLT_IsoMu20=true; 
      }
    }
  }
  std::cout<<"HLT_IsoMu20 = "<<HLT_IsoMu20<<endl;
 
  bool isMediumMuon = false;
  bool allLooseMuon = true;
  float muon_phi = 0.0;
  float muon_eta = 0.0;
  // iterate through muons for event selection
  for(long unsigned int i = 0; i<(*slimmedMuons).size(); i++){
    
    if (muon::isMediumMuon((*slimmedMuons)[i]) && (*slimmedMuons)[i].pt()>21 && (!isMediumMuon)) {
      isMediumMuon = true;
      muon_phi = (*slimmedMuons)[i].phi();
      muon_eta = (*slimmedMuons)[i].eta();
    }
    if (!muon::isLooseMuon((*slimmedMuons)[i]) && ((*slimmedMuons)[i].pt() < 10)) {
      allLooseMuon = false;
    }
  }

  std::cout<<"isMediumMuon = "<<isMediumMuon<<endl;
  std::cout<<"allLooseMuon = "<<allLooseMuon<<endl;

  if(HLT_IsoMu20 && isMediumMuon && allLooseMuon){

    float leadingjetpt = 0; float leadingjeteta = 99; float leadingjetphi = 0; float offlineHT = 0;

    // iterate through reco jets
    for(long unsigned int i = 0; i<(*slimmedJets).size(); i++){
   
      bool currentJetId = puppiJetID((*slimmedJets)[i]);
 
      if(((*slimmedJets)[i].pt()<30) || !(currentJetId)) continue;

      // get leading jet pt and eta
      float current_pt = (*slimmedJets)[i].pt();
      if(current_pt>leadingjetpt) {
        leadingjetpt = current_pt;
	leadingjeteta = (*slimmedJets)[i].eta();
        leadingjetphi = (*slimmedJets)[i].phi();
      }

      // add to offline HT
      if(abs((*slimmedJets)[i].eta()) < 2.4) {
        offlineHT += (*slimmedJets)[i].pt();
      }
 
    }

    // require deltaR>0.4 between muon and leading jet
    if (deltaR(leadingjetphi, muon_phi, leadingjeteta, muon_eta) > 0.4) {

      // fill histograms
      h_SingleJet180_den->Fill(leadingjetpt); 
      h_L1_HTT280er_den->Fill(offlineHT);
      if(L1_SingleJet180) h_SingleJet180_num->Fill(leadingjetpt); 
      if(L1_HTT280er) h_L1_HTT280er_num->Fill(offlineHT);

      if(abs(leadingjeteta) < 3.0) {
        h_SingleJet180_central_den->Fill(leadingjetpt);
        if(L1_SingleJet180) h_SingleJet180_central_num->Fill(leadingjetpt);
      }

      if(abs(leadingjeteta) >= 3.0) {
        h_SingleJet180_forward_den->Fill(leadingjetpt);
        if(L1_SingleJet180) h_SingleJet180_forward_num->Fill(leadingjetpt);
      }

      float PuppiMetPt = (float) (*slimmedMET)[0].pt();
      h_ETMHF90_den->Fill(PuppiMetPt);
      if(L1_ETMHF90) h_ETMHF90_num->Fill(PuppiMetPt);
      
      // get l1 jets
      edm::Handle<l1t::JetBxCollection> jetColl;
      iEvent.getByToken(jetBXCollectionToken_, jetColl);
      l1t::JetBxCollection jets;
      jets = (*jetColl.product());

      // find match
      TLorentzVector l1_jet_bx0;
      bool match_bx0 = checkMatchBX(leadingjetphi, leadingjeteta, jets, 0, &l1_jet_bx0);
      if (match_bx0) {
        if (l1_jet_bx0.Pt() > 30) h_SingleJet180_L130_num->Fill(leadingjetpt);
	if (l1_jet_bx0.Pt() > 60) h_SingleJet180_L160_num->Fill(leadingjetpt);
	if (l1_jet_bx0.Pt() > 120) h_SingleJet180_L1120_num->Fill(leadingjetpt);
	if (l1_jet_bx0.Pt() > 180) h_SingleJet180_L1180_num->Fill(leadingjetpt);
      }

    }
  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void L1JetAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void L1JetAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1JetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1JetAnalyzer);
