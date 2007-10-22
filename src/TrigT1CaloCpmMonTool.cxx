// ********************************************************************
//
// NAME:     TrigT1CaloCpmMonTool.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Peter Faulkner
//           
//
// ********************************************************************

#include <sstream>
#include <utility>

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "StoreGate/StoreGateSvc.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigT1Calo/CMMCPHits.h"
#include "TrigT1Calo/CPMHits.h"
#include "TrigT1Calo/CPMTower.h"
#include "TrigT1Calo/CPMRoI.h"
#include "TrigT1Calo/DataError.h"
#include "TrigT1Calo/TriggerTower.h"
#include "TrigT1Calo/TriggerTowerKey.h"
#include "TrigT1Interfaces/CoordinateRange.h"
#include "TrigT1Interfaces/CPRoIDecoder.h"
#include "TrigT1Interfaces/TrigT1CaloDefs.h"

#include "TrigT1CaloMonitoring/TrigT1CaloCpmMonTool.h"


/*---------------------------------------------------------*/
TrigT1CaloCpmMonTool::TrigT1CaloCpmMonTool(const std::string & type, 
				           const std::string & name,
				           const IInterface* parent)
  : ManagedMonitorToolBase(type, name, parent),
    m_storeGate("StoreGateSvc", name), m_monGroup(0)
/*---------------------------------------------------------*/
{
  declareInterface<IMonitorToolBase>(this); 

  declareProperty("CPMTowerLocation",
                 m_cpmTowerLocation  = LVL1::TrigT1CaloDefs::CPMTowerLocation);
  declareProperty("CPMHitsLocation",
                 m_cpmHitsLocation   = LVL1::TrigT1CaloDefs::CPMHitsLocation);
  declareProperty("CMMCPHitsLocation",
                 m_cmmCpHitsLocation = LVL1::TrigT1CaloDefs::CMMCPHitsLocation);
  declareProperty("CPMRoILocation",
                 m_cpmRoiLocation    = LVL1::TrigT1CaloDefs::CPMRoILocation);
  declareProperty("TriggerTowerLocation",
                 m_triggerTowerLocation =
		                   LVL1::TrigT1CaloDefs::TriggerTowerLocation);

  declareProperty("RootDirectory", m_rootDir = "L1Calo");
  declareProperty("SingleDirectory", m_oneDir = false);
  declareProperty("PhiUnits", m_phiUnits = "channels",
                  "Phi Units: radians, degrees or channels");

}

/*---------------------------------------------------------*/
TrigT1CaloCpmMonTool::~TrigT1CaloCpmMonTool()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloCpmMonTool:: initialize()
/*---------------------------------------------------------*/
{
  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;

  MsgStream log(msgSvc(), name());
  
  sc = m_storeGate.retrieve();
  if( sc.isFailure() ) {
    log << MSG::ERROR << "Unable to locate Service StoreGateSvc" << endreq;
    return sc;
  }

  // Phi units
  const double twoPi = 2.*M_PI;
  if      (m_phiUnits == "radians")  m_phiMax = twoPi;
  else if (m_phiUnits == "degrees")  m_phiMax = 360.;
  else if (m_phiUnits == "channels") m_phiMax = 64.;
  else {
    log << MSG::ERROR << "Invalid PhiUnits: " << m_phiUnits
        << ", using radians" << endreq;
    m_phiMax = twoPi;
  }
  m_phiScale = m_phiMax/twoPi;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloCpmMonTool::bookHistograms(bool isNewEventsBlock,
                                           bool isNewLumiBlock, bool isNewRun)
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());

  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
  	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  std::string pprDir("4_CP_PPr");
  std::string cpmDir("5_CP_CPM");
  std::string cmmDir("6_CP_CMM");
  std::string cpErrDir("05_Errors_CP");
  std::string cpmErrDir("05_Errors_CPM");
  std::string cmmErrDir("06_Errors_CMM");

  if ( isNewEventsBlock|| isNewLumiBlock) { }

  if ( isNewRun ) {	
  
  if (m_oneDir) newGroup(cpmDir, shift, eventsBlock );

  //  Timeslice checks

  newGroup(cpmDir + "_slices", shift, eventsBlock );

  for (int crate = 0; crate < 4; ++crate) {
    std::ostringstream cnum;
    cnum << crate;
    std::string name("CPM_slices_crate_" + cnum.str());
    std::string title("CPM Timeslices Crate " + cnum.str());
    TH2D* hist = book2D(name, title + ";Number of Slices;Triggered Slice",
			6, -0.5, 5.5, 5, -0.5, 4.5);
    m_v_CPM_slices.push_back(hist);
    name = "PP_CP_slice_crate_" + cnum.str();
    title = "PPr/CPM Tower Slice Match Crate " + cnum.str();
    hist = book2D(name, title + ";PPr Slice;CPM Slice",
                        5, -0.5, 4.5, 5, -0.5, 4.5); 
    m_v_PP_CP_slice.push_back(hist);
  }

  newGroup(cmmDir + "_slices", shift, eventsBlock );

  for (int crate = 0; crate < 4; ++crate) {
    std::ostringstream cnum;
    cnum << crate;
    std::string name("CMM_slices_crate_" + cnum.str());
    std::string title("CMM Timeslices Crate " + cnum.str());
    TH2D* hist = book2D(name, title + ";Number of Slices;Triggered Slice",
			6, -0.5, 5.5, 5, -0.5, 4.5);
    m_v_CMM_slices.push_back(hist);
    name = "CP_CM_slice_crate_" + cnum.str();
    title = "CPM/CMM Hits Slice Match Crate " + cnum.str();
    hist = book2D(name, title + ";CPM Slice;CMM Slice",
                        5, -0.5, 4.5, 5, -0.5, 4.5); 
    m_v_CP_CM_slice.push_back(hist);
  }

  //  CPM Tower - Trigger Tower comparison Histos

  newGroup(pprDir + "_Towers", expert, eventsBlock );

  m_h_TT_Em_Et = book1D("TT_EM_Et","Trigger Tower EM Et Noise",20,0,20);
  m_h_TT_Had_Et = book1D("TT_HAD_Et","Trigger Tower HAD Et Noise",20,0,20);
  m_h_TT_Em_Et_s = book1D("TT_EM_Et_s","Trigger Tower EM Et Signal",235,20,255);
  m_h_TT_Had_Et_s = book1D("TT_HAD_Et_s","Trigger Tower HAD Et Signal",235,20,255);
  m_h_TT_Em_eta = book1D("TT_EM_eta","Trigger Tower EM eta",50,-2.5,2.5);
  m_h_TT_Had_eta = book1D("TT_HAD_eta","Trigger Tower HAD eta",50,-2.5,2.5);
  m_h_TT_Em_phi = book1D("TT_EM_phi","Trigger Tower EM phi ",64,0,m_phiMax);
  m_h_TT_Had_phi = book1D("TT_HAD_phi","Trigger Tower HAD phi ",64,0,m_phiMax);
  m_h_TT_Em_eta_phi = book2D("TT_EM_eta_phi",
         "Trigger Tower EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TT_Had_eta_phi = book2D("TT_HAD_eta_phi",
         "Trigger Tower HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TT_Em_eta_phi_w = book2D("TT_EM_eta_phi_w",
      "Trigger Tower EM eta/phi weighted;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TT_Had_eta_phi_w = book2D("TT_HAD_eta_phi_w",
      "Trigger Tower HAD eta/phi weighted;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);

  newGroup(cpmDir + "_Towers", shift, eventsBlock );

  m_h_CT_Em_Et = book1D("CT_EM_Et","CPM Tower EM Et Noise",20,0,20);
  m_h_CT_Had_Et = book1D("CT_HAD_Et","CPM Tower HAD Et Noise",20,0,20);
  m_h_CT_Em_Et_s = book1D("CT_EM_Et_s","CPM Tower EM Et Signal",235,20,255);
  m_h_CT_Had_Et_s = book1D("CT_HAD_Et_s","CPM Tower HAD Et Signal",235,20,255);
  m_h_CT_Em_eta = book1D("CT_EM_eta","CPM Tower EM eta",50,-2.5,2.5);
  m_h_CT_Had_eta = book1D("CT_HAD_eta","CPM Tower HAD eta",50,-2.5,2.5);
  m_h_CT_Em_phi = book1D("CT_EM_phi","CPM Tower EM phi ",64,0,m_phiMax);
  m_h_CT_Had_phi = book1D("CT_HAD_phi","CPM Tower HAD phi ",64,0,m_phiMax);
  m_h_CT_Em_eta_phi = book2D("CT_EM_eta_phi",
         "CPM Tower EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CT_Had_eta_phi = book2D("CT_HAD_eta_phi",
         "CPM Tower HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CT_Em_eta_phi_w = book2D("CT_EM_eta_phi_w",
         "CPM Tower EM eta/phi weighted;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CT_Had_eta_phi_w = book2D("CT_HAD_eta_phi_w",
         "CPM Tower HAD eta/phi weighted;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);

  //  CPM Tower error bits

  newGroup(cpmErrDir + "_Towers", shift, eventsBlock );

  m_h_CT_Em_parity = book2D("CT_EM_parity",
            "CPM Tower EM Parity Errors;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CT_Had_parity = book2D("CT_HAD_parity",
            "CPM Tower HAD Parity Errors;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CT_Em_link = book2D("CT_EM_link",
            "CPM Tower EM Link Down Errors;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CT_Had_link = book2D("CT_HAD_link",
         "CPM Tower HAD Link Down Errors;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CT_status = book1D("CT_status", "CPM Sub-status bits", 8, 0., 8.);
  setStatusLabels(m_h_CT_status);
  m_h_CT_status_eta_phi = book2D("CT_status_eta_phi",
            "CPM Sub-status hit-map;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);

  //  Multiple slice CPM Tower plots

  //....

  //  Trigger Tower/CPM Tower event by event comparison plots

  newGroup(cpmDir + "_input", shift, eventsBlock );

  m_h_TTeqCT_Em_eta_phi = book2D("TTeqCT_EM_eta_phi",
    "Trigger/CPM Tower match EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TTneCT_Em_eta_phi = book2D("TTneCT_EM_eta_phi",
    "Trigger/CPM Tower mismatch EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TTnoCT_Em_eta_phi = book2D("TTnoCT_EM_eta_phi",
    "Trigger Tower/no CPM Tower EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CTnoTT_Em_eta_phi = book2D("CTnoTT_EM_eta_phi",
    "CPM Tower/no Trigger Tower EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TTeqCT_Had_eta_phi = book2D("TTeqCT_HAD_eta_phi",
    "Trigger/CPM Tower match HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TTneCT_Had_eta_phi = book2D("TTneCT_HAD_eta_phi",
   "Trigger/CPM Tower mismatch HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TTnoCT_Had_eta_phi = book2D("TTnoCT_HAD_eta_phi",
   "Trigger Tower/no CPM Tower HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CTnoTT_Had_eta_phi = book2D("CTnoTT_HAD_eta_phi",
   "CPM Tower/no Trigger Tower HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);

  //  CPM RoIs

  newGroup(cpmDir + "_RoIs", shift, eventsBlock );

  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("RoI_Thresh_" + cnum.str());
    std::string title("RoI Threshold " + cnum.str() + ";Crate/CPM");
    TH1D* hist = book1D(name, title, 56, 0, 56);
    setThresholdLabels(hist);
    m_v_RoI_thresholds.push_back(hist);
  }
  const double halfPhiBin = m_phiMax/128.;
  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("RoI_2D_Thresh_" + cnum.str());
    std::string title("RoI eta/phi Threshold " + cnum.str() + ";eta;phi");
    m_v_RoI_2D_thresholds.push_back(book2D(name, title, 51, -2.55, 2.55,
                                    65, -halfPhiBin, m_phiMax+halfPhiBin));
  }

  newGroup(cpmErrDir + "_RoIs", shift, eventsBlock );

  m_h_RoI_Parity = book2D("CPM_RoI_Parity",
            "CPM RoI Parity Errors;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);

  //  CPM Hits
  
  newGroup(cpmDir + "_Hits", shift, eventsBlock );

  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("CH_Thresh_" + cnum.str());
    std::string title("CPM Hits Threshold " + cnum.str() + ";Crate/CPM");
    TH1D* hist = book1D(name, title, 56, 0, 56);
    setThresholdLabels(hist);
    m_v_thresholds.push_back(hist);
  }

  //  Multiple slice CPM Hit plots

  //....

  //  CMM-CP Hits

  newGroup(cmmDir + "_Hits", shift, eventsBlock );

  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("CM_Thresh_" + cnum.str());
    std::string title("CMM-CP Hits Threshold " + cnum.str() + ";Crate/CPM");
    TH1D* hist = book1D(name, title, 56, 0, 56);
    setThresholdLabels(hist);
    m_v_CMM_thresholds.push_back(hist);
  }
  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("CM_T_Thresh_" + cnum.str());
    std::string title("CMM-CP Hits Totals Threshold "
                                                  + cnum.str() + ";By Crate");
    TH1D* hist = book1D(name, title, 20, 0, 20);
    int bin = 0;
    for (int crate = 0; crate < 4; ++crate) {
      hist->GetXaxis()->SetBinLabel(++bin, "Remote");
      hist->GetXaxis()->SetBinLabel(++bin, "Remote");
      hist->GetXaxis()->SetBinLabel(++bin, "Remote");
      hist->GetXaxis()->SetBinLabel(++bin, "Local");
      hist->GetXaxis()->SetBinLabel(++bin, "Total");
    }
    m_v_CMM_T_thresholds.push_back(hist);
  }

  //  CMM error bits

  newGroup(cmmErrDir, shift, eventsBlock );

  m_h_CMM_L_parity = book1D("CMM_L_parity",
                 "CMM Parity Errors Em/Tau (Left);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CMM_L_parity);
  m_h_CMM_R_parity = book1D("CMM_R_parity",
                 "CMM Parity Errors Em (Right);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CMM_R_parity);
  m_h_CMM_status = book1D("CMM_status", "CMM Sub-status bits", 8, 0., 8.);
  setStatusLabels(m_h_CMM_status);
  m_h_CMM_status_loc = book2D("CMM_status_loc",
                              "CMM Sub-status location;;Crate/Left-Right",
                              8, 0., 8., 8, 0., 8.);
  setStatusLabels(m_h_CMM_status_loc);
  setCmmLocLabels(m_h_CMM_status_loc);

  //  Multiple slice CMM Hit plots

  //....

  //  CPM/CMM Hits event by event comparisons

  newGroup(cmmDir + "_input", shift, eventsBlock );

  m_h_CPMeqCMM_hits1 = book1D("CPMeqCMM_hits1",
                 "CPM-CMM Hits match Em/Tau (Left);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPMeqCMM_hits1);
  m_h_CPMeqCMM_hits0 = book1D("CPMeqCMM_hits0",
                 "CPM-CMM Hits match Em (Right);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPMeqCMM_hits0);
  m_h_CPMneCMM_hits1 = book1D("CPMneCMM_hits1",
                 "CPM-CMM Hits mismatch Em/Tau (Left);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPMneCMM_hits1);
  m_h_CPMneCMM_hits0 = book1D("CPMneCMM_hits0",
                 "CPM-CMM Hits mismatch Em (Right);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPMneCMM_hits0);
  m_h_CPMnoCMM_hits1 = book1D("CPMnoCMM_hits1",
                 "CPM /no CMM Hits Em/Tau (Left);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPMnoCMM_hits1);
  m_h_CPMnoCMM_hits0 = book1D("CPMnoCMM_hits0",
                 "CPM /no CMM Hits Em (Right);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPMnoCMM_hits0);
  m_h_CMMnoCPM_hits1 = book1D("CMMnoCPM_hits1",
                 "CMM /no CPM Hits Em/Tau (Left);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CMMnoCPM_hits1);
  m_h_CMMnoCPM_hits0 = book1D("CMMnoCPM_hits0",
                 "CMM /no CPM Hits Em (Right);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CMMnoCPM_hits0);

  //  Error Summary

  newGroup(cpErrDir + "_summary", shift, eventsBlock );

  m_h_CP_errors = book1D("CP_Error_Summary", "CP Error Summary", 9, 0, 9);

  delete m_monGroup;
  m_monGroup = 0;

  } // end if (isNewEventsBlock ...

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloCpmMonTool::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  

  //Retrieve Trigger Towers from SG
  const TriggerTowerCollection* triggerTowerTES = 0; 
  StatusCode sc = m_storeGate->retrieve(triggerTowerTES,
                                                     m_triggerTowerLocation); 
  if( sc.isFailure()  ||  !triggerTowerTES ) {
    log << MSG::ERROR<< "No Trigger Tower container found"<< endreq; 
  }

  //Retrieve CPM Towers from SG
  const CpmTowerCollection* cpmTowerTES = 0; 
  sc = m_storeGate->retrieve(cpmTowerTES, m_cpmTowerLocation); 
  if( sc.isFailure()  ||  !cpmTowerTES ) {
    log << MSG::ERROR<< "No CPM Tower container found"<< endreq; 
  }
  
  //Retrieve CPM RoIs from SG
  const CpmRoiCollection* cpmRoiTES = 0;
  sc = m_storeGate->retrieve( cpmRoiTES, m_cpmRoiLocation);
  if( sc.isFailure()  ||  !cpmRoiTES ) {
    log << MSG::ERROR << "No CPM RoIs container found"<< endreq; 
  }
  
  //Retrieve CPM Hits from SG
  const CpmHitsCollection* cpmHitsTES = 0;
  sc = m_storeGate->retrieve( cpmHitsTES, m_cpmHitsLocation);
  if( sc.isFailure()  ||  !cpmHitsTES ) {
    log << MSG::ERROR << "No CPM Hits container found"<< endreq; 
  }
  
  //Retrieve CMM-CP Hits from SG
  const CmmCpHitsCollection* cmmCpHitsTES = 0;
  sc = m_storeGate->retrieve( cmmCpHitsTES, m_cmmCpHitsLocation);
  if( sc.isFailure()  ||  !cmmCpHitsTES ) {
    log << MSG::ERROR << "No CMM-CP Hits container found"<< endreq; 
  }

  //=============================================
  //   CPM Tower - Trigger Tower comparison plots
  //=============================================

  // Maps for one-one comparisons
  TriggerTowerMap ttMap;
  CpmTowerMap     cpMap;
  LVL1::TriggerTowerKey towerKey;
  unsigned int maxKey = 0;

  // Global plots

  if (triggerTowerTES) {
    TriggerTowerCollection::const_iterator ttIterator =
                                                      triggerTowerTES->begin(); 
    TriggerTowerCollection::const_iterator ttIteratorEnd =
                                                      triggerTowerTES->end(); 
    for (; ttIterator != ttIteratorEnd; ++ttIterator) {
      const int    em  = (*ttIterator)->emEnergy();
      const int    had = (*ttIterator)->hadEnergy();
      const double eta = (*ttIterator)->eta();
      const double phi = (*ttIterator)->phi();
      const double phiMod = phi * m_phiScale;
      if (em && eta > -2.5 && eta < 2.5) {
        m_h_TT_Em_Et->Fill(em, 1.);
        if (em > 10) m_h_TT_Em_Et_s->Fill(em, 1.);
        m_h_TT_Em_eta->Fill(eta, 1.);
        m_h_TT_Em_phi->Fill(phiMod, 1.);
        m_h_TT_Em_eta_phi->Fill(eta, phiMod, 1.);
        m_h_TT_Em_eta_phi_w->Fill(eta, phiMod, em);
      }
      if (had && eta > -2.5 && eta < 2.5) {
        m_h_TT_Had_Et->Fill(had, 1.);
        if (had > 10) m_h_TT_Had_Et_s->Fill(had, 1.);
        m_h_TT_Had_eta->Fill(eta, 1.);
        m_h_TT_Had_phi->Fill(phiMod, 1.);
        m_h_TT_Had_eta_phi->Fill(eta, phiMod, 1.);
        m_h_TT_Had_eta_phi_w->Fill(eta, phiMod, had);
      }
      if (eta > -2.5 && eta < 2.5) {
	const unsigned int key = towerKey.ttKey(phi, eta);
	if (key > maxKey) maxKey = key;
        ttMap.insert(std::make_pair(key, *ttIterator));
      }
    }
  }

  if (cpmTowerTES) {
    CpmTowerCollection::const_iterator ctIterator    = cpmTowerTES->begin(); 
    CpmTowerCollection::const_iterator ctIteratorEnd = cpmTowerTES->end(); 

    for (; ctIterator != ctIteratorEnd; ++ctIterator) {
      const int    em  = (*ctIterator)->emEnergy();
      const int    had = (*ctIterator)->hadEnergy();
      const double eta = (*ctIterator)->eta();
      const double phi = (*ctIterator)->phi();
      const double phiMod = phi * m_phiScale;
      const int    crate  = static_cast<int>(phi/(M_PI/2.));
      const int    peak = (*ctIterator)->peak();
      const int    slices = ((*ctIterator)->emEnergyVec()).size();
      m_v_CPM_slices[crate]->Fill(slices, peak, 1.);
      if (em) {
        m_h_CT_Em_Et->Fill(em, 1.);
        if (em > 10) m_h_CT_Em_Et_s->Fill(em, 1.);
        m_h_CT_Em_eta->Fill(eta, 1.);
        m_h_CT_Em_phi->Fill(phiMod, 1.);
        m_h_CT_Em_eta_phi->Fill(eta, phiMod, 1.);
        m_h_CT_Em_eta_phi_w->Fill(eta, phiMod, em);
      }
      if (had) {
        m_h_CT_Had_Et->Fill(had, 1.);
        if (had > 10) m_h_CT_Had_Et_s->Fill(had, 1.);
        m_h_CT_Had_eta->Fill(eta, 1.);
        m_h_CT_Had_phi->Fill(phiMod, 1.);
        m_h_CT_Had_eta_phi->Fill(eta, phiMod, 1.);
        m_h_CT_Had_eta_phi_w->Fill(eta, phiMod, had);
      }
      // Errors
      const LVL1::DataError emError((*ctIterator)->emError());
      const LVL1::DataError hadError((*ctIterator)->hadError());
      m_h_CT_Em_parity->Fill(eta, phiMod, emError.get(LVL1::DataError::Parity));
      m_h_CT_Had_parity->Fill(eta, phiMod,
                                   hadError.get(LVL1::DataError::Parity));
      m_h_CT_Em_link->Fill(eta, phiMod, emError.get(LVL1::DataError::LinkDown));
      m_h_CT_Had_link->Fill(eta, phiMod,
                                 hadError.get(LVL1::DataError::LinkDown));
      // Sub-status errors
      const int status = emError.error() >> LVL1::DataError::GLinkParity;
      for (int bit = 0; bit < 8; ++bit) {
        m_h_CT_status->Fill(bit, (status >> bit) & 0x1);
      }
      m_h_CT_status_eta_phi->Fill(eta, phiMod, status != 0);

      const unsigned int key = towerKey.ttKey(phi, eta);
      if (key > maxKey) maxKey = key;
      cpMap.insert(std::make_pair(key, *ctIterator));
    }
  }
  ++maxKey;

  // One-to-one tower comparison

  TriggerTowerMap::const_iterator ttMapIter    = ttMap.begin();
  TriggerTowerMap::const_iterator ttMapIterEnd = ttMap.end();
  CpmTowerMap::const_iterator     cpMapIter    = cpMap.begin();
  CpmTowerMap::const_iterator     cpMapIterEnd = cpMap.end();

  while (ttMapIter != ttMapIterEnd || cpMapIter != cpMapIterEnd) {

    unsigned int ttKey = maxKey;
    unsigned int cpKey = maxKey;
    int ttEm  = 0;
    int ttHad = 0;
    int cpEm  = 0;
    int cpHad = 0;
    int ttEmErr  = 0;
    int ttHadErr = 0;
    int cpEmErr  = 0;
    int cpHadErr = 0;
    double eta = 0.;
    double phi = 0.;

    if (ttMapIter != ttMapIterEnd) ttKey = ttMapIter->first;
    if (cpMapIter != cpMapIterEnd) cpKey = cpMapIter->first;

    if ((cpMapIter == cpMapIterEnd) || (cpKey > ttKey)) {

      // TriggerTower but no CPMTower

      const LVL1::TriggerTower* tt = ttMapIter->second;
      ttEm  = tt->emEnergy();
      ttHad = tt->hadEnergy();
      ttEmErr  = tt->emError();
      ttHadErr = tt->hadError();
      eta = tt->eta();
      phi = tt->phi();
      ++ttMapIter;

    } else if ((ttMapIter == ttMapIterEnd) || (ttKey > cpKey)) {

      // CPMTower but no TriggerTower

      const LVL1::CPMTower* cp = cpMapIter->second;
      cpEm  = cp->emEnergy();
      cpHad = cp->hadEnergy();
      cpEmErr  = cp->emError();
      cpHadErr = cp->hadError();
      eta = cp->eta();
      phi = cp->phi();
      ++cpMapIter;

    } else {

      // Have both

      const LVL1::TriggerTower* tt = ttMapIter->second;
      const LVL1::CPMTower*     cp = cpMapIter->second;
      ttEm  = tt->emEnergy();
      ttHad = tt->hadEnergy();
      cpEm  = cp->emEnergy();
      cpHad = cp->hadEnergy();
      ttEmErr  = tt->emError();
      ttHadErr = tt->hadError();
      cpEmErr  = cp->emError();
      cpHadErr = cp->hadError();
      eta = tt->eta();
      phi = tt->phi();
      // Slice match
      const std::vector<int>& emLut(tt->emLUT());
      const std::vector<int>& hadLut(tt->hadLUT());
      const std::vector<int>& emVec(cp->emEnergyVec());
      const std::vector<int>& hadVec(cp->hadEnergyVec());
      const int sliceEmLut = emLut.size();
      const int sliceHadLut = hadLut.size();
      const int sliceEmVec = emVec.size();
      const int sliceHadVec = hadVec.size();
      const int crate = static_cast<int>(phi/(M_PI/2.));
      for (int slice = 0; slice < sliceEmLut; ++slice) {
        if (emLut[slice] > 0) {
	  for (int slice2 = 0; slice2 < sliceEmVec; ++slice2) {
	    if (emLut[slice] == emVec[slice2]) {
	      m_v_PP_CP_slice[crate]->Fill(slice, slice2, 1.);
            }
          }
        }
      }
      for (int slice = 0; slice < sliceHadLut; ++slice) {
        if (hadLut[slice] > 0) {
	  for (int slice2 = 0; slice2 < sliceHadVec; ++slice2) {
	    if (hadLut[slice] == hadVec[slice2]) {
	      m_v_PP_CP_slice[crate]->Fill(slice, slice2, 1.);
            }
          }
        }
      }
      ++ttMapIter;
      ++cpMapIter;
    }
    const double phiMod = phi * m_phiScale;
    // Ignore data with known errors
    //if ( !ttEmErr && !cpEmErr ) {
      m_h_TTeqCT_Em_eta_phi->Fill(eta, phiMod, ttEm && cpEm && ttEm == cpEm);
      m_h_TTneCT_Em_eta_phi->Fill(eta, phiMod, ttEm && cpEm && ttEm != cpEm);
      m_h_TTnoCT_Em_eta_phi->Fill(eta, phiMod, ttEm  && !cpEm);
      m_h_CTnoTT_Em_eta_phi->Fill(eta, phiMod, !ttEm && cpEm);
    //}
    //if ( !ttHadErr && !cpHadErr ) {
      m_h_TTeqCT_Had_eta_phi->Fill(eta, phiMod,
                                        ttHad && cpHad && ttHad == cpHad);
      m_h_TTneCT_Had_eta_phi->Fill(eta, phiMod,
                                        ttHad && cpHad && ttHad != cpHad);
      m_h_TTnoCT_Had_eta_phi->Fill(eta, phiMod, ttHad && !cpHad);
      m_h_CTnoTT_Had_eta_phi->Fill(eta, phiMod, !ttHad && cpHad);
    //}
    //if ((!ttEmErr && !cpEmErr && ttEm && cpEm && ttEm != cpEm) ||
    //    (!ttHadErr && !cpHadErr && ttHad && cpHad && ttHad != cpHad)) {
    if ((ttEm && cpEm && ttEm != cpEm) || (ttHad && cpHad && ttHad != cpHad)) {
      log << MSG::DEBUG
          << "Trigger/CPM Tower mismatch, eta/phi/ttEm/ttHad/cpEm/cpHad: "
	  << eta << "/" << phi << "/" << ttEm << "/" << ttHad << "/"
	  << cpEm << "/" << cpHad << endreq;
    }
  }

  //=============================================
  //  CPM RoIs
  //=============================================

  if (cpmRoiTES) {
    LVL1::CPRoIDecoder decoder;
    CpmRoiCollection::const_iterator crIterator    = cpmRoiTES->begin(); 
    CpmRoiCollection::const_iterator crIteratorEnd = cpmRoiTES->end(); 
    for (; crIterator != crIteratorEnd; ++crIterator) {
      const int hits  = (*crIterator)->hits();
      const LVL1::CoordinateRange coord(
                                decoder.coordinate((*crIterator)->roiWord()));
      const double eta = coord.eta();
      const double phi = coord.phi();
      const double phiMod = phi * m_phiScale;
      const int crate = (*crIterator)->crate();
      const int cpm   = (*crIterator)->cpm();
      //const int chip  = (*crIterator)->chip();
      //const int loc   = (*crIterator)->location();
      const int bin1  = crate * 14 + cpm - 1;
      //const int bin2  = chip * 8 + loc;
      std::vector<TH1D*>::const_iterator hist1 = m_v_RoI_thresholds.begin();
      std::vector<TH2D*>::const_iterator hist2 = m_v_RoI_2D_thresholds.begin();
      for (int thresh = 0; thresh < 16; ++thresh) {
        const int hit = (hits >> thresh) & 0x1;
        if (hit) {
          (*hist1)->Fill(bin1, 1.);
	  (*hist2)->Fill(eta, phiMod, 1.);
        }
        ++hist1;
        ++hist2;
      }
      if ((*crIterator)->error()) m_h_RoI_Parity->Fill(eta, phiMod, 1.);
    }
  }

  //=============================================
  //  CPM Hits
  //=============================================

  CpmHitsMap cpmMap;
  maxKey = 0;

  if (cpmHitsTES) {
    CpmHitsCollection::const_iterator chIterator    = cpmHitsTES->begin(); 
    CpmHitsCollection::const_iterator chIteratorEnd = cpmHitsTES->end(); 
    for (; chIterator != chIteratorEnd; ++chIterator) {
      const unsigned int hits0 = (*chIterator)->HitWord0();
      const unsigned int hits1 = (*chIterator)->HitWord1();
      const int crate = (*chIterator)->crate();
      const int cpm   = (*chIterator)->module();
      const int bin   = crate * 14 + cpm - 1;
      const int peak   = (*chIterator)->peak();
      const int slices = ((*chIterator)->HitsVec0()).size();
      m_v_CPM_slices[crate]->Fill(slices, peak, 1.);
      std::vector<TH1D*>::const_iterator hist = m_v_thresholds.begin();
      for (int thresh = 0; thresh < 8; ++thresh) {
        const int hit = (hits0 >> thresh*3) & 0x7;
        if (hit) (*hist)->Fill(bin, hit);
        ++hist;
      }
      for (int thresh = 0; thresh < 8; ++thresh) {
        const int hit = (hits1 >> thresh*3) & 0x7;
        if (hit) (*hist)->Fill(bin, hit);
        ++hist;
      }
      const unsigned int key = crate * 100 + cpm;
      if (key > maxKey) maxKey = key;
      cpmMap.insert(std::make_pair(key, *chIterator));
    }
  }

  //=============================================
  //  CMM-CP Hits
  //=============================================

  CmmCpHitsMap cmmMap;

  if (cmmCpHitsTES) {
    CmmCpHitsCollection::const_iterator cmIterator    = cmmCpHitsTES->begin(); 
    CmmCpHitsCollection::const_iterator cmIteratorEnd = cmmCpHitsTES->end(); 
    for (; cmIterator != cmIteratorEnd; ++cmIterator) {
      const unsigned int hits0 = (*cmIterator)->HitWord0();
      const unsigned int hits1 = (*cmIterator)->HitWord1();
      const int crate  = (*cmIterator)->crate();
      const int dataId = (*cmIterator)->dataID();
      const int bin    = (dataId < 15) ? crate * 14 + dataId - 1
  		                       : crate * 5  + dataId - 15; 
      const int peak   = (*cmIterator)->peak();
      const int slices = ((*cmIterator)->HitsVec0()).size();
      m_v_CMM_slices[crate]->Fill(slices, peak, 1.);
      std::vector<TH1D*>::const_iterator hist1 = m_v_CMM_thresholds.begin();
      std::vector<TH1D*>::const_iterator hist2 = m_v_CMM_T_thresholds.begin();
      for (int thresh = 0; thresh < 8; ++thresh) {
        const int hit = (hits0 >> thresh*3) & 0x7;
        if (dataId < 15) {
          if (hit) (*hist1)->Fill(bin, hit);
          ++hist1;
        } else {
          if (hit) (*hist2)->Fill(bin, hit);
	  ++hist2;
        }
      }
      for (int thresh = 0; thresh < 8; ++thresh) {
        const int hit = (hits1 >> thresh*3) & 0x7;
        if (dataId < 15) {
          if (hit) (*hist1)->Fill(bin, hit);
          ++hist1;
        } else {
          if (hit) (*hist2)->Fill(bin, hit);
	  ++hist2;
        }
      }
      // Errors
      const LVL1::DataError hit0Err((*cmIterator)->Error0());
      const LVL1::DataError hit1Err((*cmIterator)->Error1());
      if (dataId < 15) {
        m_h_CMM_R_parity->Fill(bin, hit0Err.get(LVL1::DataError::Parity));
        m_h_CMM_L_parity->Fill(bin, hit1Err.get(LVL1::DataError::Parity));
      }
      // Sub-status errors
      const int status0 = hit0Err.error() >> LVL1::DataError::GLinkParity;
      const int status1 = hit1Err.error() >> LVL1::DataError::GLinkParity;
      for (int bit = 0; bit < 8; ++bit) {
        m_h_CMM_status->Fill(bit, (status0 >> bit) & 0x1);
        m_h_CMM_status->Fill(bit, (status1 >> bit) & 0x1);
	m_h_CMM_status_loc->Fill(bit, 2*crate, (status1 >> bit) & 0x1);
	m_h_CMM_status_loc->Fill(bit, 2*crate + 1, (status0 >> bit) & 0x1);
      }

      if (dataId < 15) {
        const unsigned int key = crate * 100 + dataId;
        if (key > maxKey) maxKey = key;
        cmmMap.insert(std::make_pair(key, *cmIterator));
      }
    }
  }
  ++maxKey;

  // One-to-one CPM-CMM hit comparison

  CpmHitsMap::const_iterator   cpmMapIter    = cpmMap.begin();
  CpmHitsMap::const_iterator   cpmMapIterEnd = cpmMap.end();
  CmmCpHitsMap::const_iterator cmmMapIter    = cmmMap.begin();
  CmmCpHitsMap::const_iterator cmmMapIterEnd = cmmMap.end();

  while (cpmMapIter != cpmMapIterEnd || cmmMapIter != cmmMapIterEnd) {

    unsigned int cpmKey = maxKey;
    unsigned int cmmKey = maxKey;
    unsigned int cpmHits0 = 0;
    unsigned int cpmHits1 = 0;
    unsigned int cmmHits0 = 0;
    unsigned int cmmHits1 = 0;
    int crate = 0;
    int cpm   = 0;

    if (cpmMapIter != cpmMapIterEnd) cpmKey = cpmMapIter->first;
    if (cmmMapIter != cmmMapIterEnd) cmmKey = cmmMapIter->first;

    if ((cmmMapIter == cmmMapIterEnd) || (cmmKey > cpmKey)) {

      // CPM Hits but no CMM Hits

      const LVL1::CPMHits* cpmh = cpmMapIter->second;
      cpmHits0 = cpmh->HitWord0();
      cpmHits1 = cpmh->HitWord1();
      crate    = cpmh->crate();
      cpm      = cpmh->module();
      ++cpmMapIter;

    } else if ((cpmMapIter == cpmMapIterEnd) || (cpmKey > cmmKey)) {

      // CMM Hits but no CPM Hits

      const LVL1::CMMCPHits* cmmh = cmmMapIter->second;
      cmmHits0 = cmmh->HitWord0();
      cmmHits1 = cmmh->HitWord1();
      crate    = cmmh->crate();
      cpm      = cmmh->dataID();
      ++cmmMapIter;

    } else {

      // Have both

      const LVL1::CPMHits*   cpmh = cpmMapIter->second;
      const LVL1::CMMCPHits* cmmh = cmmMapIter->second;
      cpmHits0 = cpmh->HitWord0();
      cpmHits1 = cpmh->HitWord1();
      cmmHits0 = cmmh->HitWord0();
      cmmHits1 = cmmh->HitWord1();
      crate    = cpmh->crate();
      cpm      = cpmh->module();
      // Slice match
      const std::vector<unsigned int>& cpmVec0(cpmh->HitsVec0());
      const std::vector<unsigned int>& cpmVec1(cpmh->HitsVec1());
      const std::vector<unsigned int>& cmmVec0(cmmh->HitsVec0());
      const std::vector<unsigned int>& cmmVec1(cmmh->HitsVec1());
      const int sliceCpmVec0 = cpmVec0.size();
      const int sliceCpmVec1 = cpmVec1.size();
      const int sliceCmmVec0 = cmmVec0.size();
      const int sliceCmmVec1 = cmmVec1.size();
      if (sliceCpmVec0 == sliceCpmVec1 && sliceCmmVec0 == sliceCmmVec1) {
        for (int slice = 0; slice < sliceCpmVec0; ++slice) {
          if (cpmVec0[slice] > 0 || cpmVec1[slice] > 0) {
	    for (int slice2 = 0; slice2 < sliceCmmVec0; ++slice2) {
	      if (cpmVec0[slice] == cmmVec0[slice2] &&
	          cpmVec1[slice] == cmmVec1[slice2]) {
	        m_v_CP_CM_slice[crate]->Fill(slice, slice2, 1.);
	      }
            }
          }
        }
      }
      ++cpmMapIter;
      ++cmmMapIter;
    }
    const int bin = crate * 14 + cpm - 1;
    m_h_CPMeqCMM_hits0->Fill(bin, cpmHits0 && cmmHits0 && cpmHits0 == cmmHits0);
    m_h_CPMeqCMM_hits1->Fill(bin, cpmHits1 && cmmHits1 && cpmHits1 == cmmHits1);
    m_h_CPMneCMM_hits0->Fill(bin, cpmHits0 && cmmHits0 && cpmHits0 != cmmHits0);
    m_h_CPMneCMM_hits1->Fill(bin, cpmHits1 && cmmHits1 && cpmHits1 != cmmHits1);
    m_h_CPMnoCMM_hits0->Fill(bin, cpmHits0 && !cmmHits0);
    m_h_CPMnoCMM_hits1->Fill(bin, cpmHits1 && !cmmHits1);
    m_h_CMMnoCPM_hits0->Fill(bin, !cpmHits0 && cmmHits0);
    m_h_CMMnoCPM_hits1->Fill(bin, !cpmHits1 && cmmHits1);
    if ((cpmHits0 && cmmHits0 && cpmHits0 != cmmHits0) ||
        (cpmHits1 && cmmHits1 && cpmHits1 != cmmHits1)) {
      log << MSG::DEBUG << "CPM/CMM Hits mismatch, crate/cpm/hitsCPM/hitsCMM: "
          << crate << "/" << cpm << "/";
      const int nthresh = 8;
      const int bits = 3;
      const int mask = 0x7;
      for (int thresh = 0; thresh < nthresh; ++thresh) {
        const int hit = (cpmHits0 >> thresh*bits) & mask;
	log << MSG::DEBUG << hit;
	if (thresh < nthresh - 1) log << MSG::DEBUG << ",";
	else log << MSG::DEBUG << ";";
      }
      for (int thresh = 0; thresh < nthresh; ++thresh) {
        const int hit = (cpmHits1 >> thresh*bits) & mask;
	log << MSG::DEBUG << hit;
	if (thresh < nthresh - 1) log << MSG::DEBUG << ",";
	else log << MSG::DEBUG << "/";
      }
      for (int thresh = 0; thresh < nthresh; ++thresh) {
        const int hit = (cmmHits0 >> thresh*bits) & mask;
	log << MSG::DEBUG << hit;
	if (thresh < nthresh - 1) log << MSG::DEBUG << ",";
	else log << MSG::DEBUG << ";";
      }
      for (int thresh = 0; thresh < nthresh; ++thresh) {
        const int hit = (cmmHits1 >> thresh*bits) & mask;
	log << MSG::DEBUG << hit;
	if (thresh < nthresh - 1) log << MSG::DEBUG << ",";
	else log << MSG::DEBUG << endreq;
      }
    }
  }

  return StatusCode::SUCCESS;

}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloCpmMonTool::procHistograms(bool isEndOfEventsBlock,
                                  bool isEndOfLumiBlock, bool isEndOfRun)
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());

  if (isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun) {
    // Fill Error summary hist
    int bin = 0;
    m_h_CP_errors->Fill(bin++, m_h_TTneCT_Em_eta_phi->GetEffectiveEntries() +
                               m_h_TTnoCT_Em_eta_phi->GetEffectiveEntries() +
			       m_h_CTnoTT_Em_eta_phi->GetEffectiveEntries() +
			       m_h_TTneCT_Had_eta_phi->GetEffectiveEntries() +
			       m_h_TTnoCT_Had_eta_phi->GetEffectiveEntries() +
			       m_h_CTnoTT_Had_eta_phi->GetEffectiveEntries()
			                                                > 0.);
    m_h_CP_errors->Fill(bin++, m_h_CT_Em_parity->GetEffectiveEntries() +
                               m_h_CT_Had_parity->GetEffectiveEntries() > 0.);
    m_h_CP_errors->Fill(bin++, m_h_CT_Em_link->GetEffectiveEntries() +
                               m_h_CT_Had_link->GetEffectiveEntries() > 0.);
    m_h_CP_errors->Fill(bin++, m_h_CT_status->GetEffectiveEntries() > 0.);
    m_h_CP_errors->Fill(bin++, m_h_RoI_Parity->GetEffectiveEntries() > 0.);
    m_h_CP_errors->Fill(bin++, m_h_CPMneCMM_hits0->GetEffectiveEntries() +
                               m_h_CPMneCMM_hits1->GetEffectiveEntries() +
			       m_h_CPMnoCMM_hits0->GetEffectiveEntries() +
			       m_h_CPMnoCMM_hits1->GetEffectiveEntries() +
			       m_h_CMMnoCPM_hits0->GetEffectiveEntries() +
			       m_h_CMMnoCPM_hits1->GetEffectiveEntries() > 0.);
    m_h_CP_errors->Fill(bin++, m_h_CMM_R_parity->GetEffectiveEntries() +
                               m_h_CMM_L_parity->GetEffectiveEntries() > 0.);
    m_h_CP_errors->Fill(bin++, m_h_CMM_status->GetEffectiveEntries() > 0.);
    int sumTe = 0;
    for (int thresh = 0; thresh < 16; ++thresh) {
      sumTe += (m_v_CMM_T_thresholds[thresh]->GetBinContent(4) !=
                m_v_CMM_T_thresholds[thresh]->GetBinContent(16)) +
               (m_v_CMM_T_thresholds[thresh]->GetBinContent(9) !=
	        m_v_CMM_T_thresholds[thresh]->GetBinContent(17)) +
               (m_v_CMM_T_thresholds[thresh]->GetBinContent(14) !=
	        m_v_CMM_T_thresholds[thresh]->GetBinContent(18));
    }
    m_h_CP_errors->Fill(bin, sumTe > 0);
    bin = 1;
    m_h_CP_errors->GetXaxis()->SetBinLabel(bin++, "PPr-CPM");
    m_h_CP_errors->GetXaxis()->SetBinLabel(bin++, "CPM parity");
    m_h_CP_errors->GetXaxis()->SetBinLabel(bin++, "CPM link");
    m_h_CP_errors->GetXaxis()->SetBinLabel(bin++, "CPM status");
    m_h_CP_errors->GetXaxis()->SetBinLabel(bin++, "RoI parity");
    m_h_CP_errors->GetXaxis()->SetBinLabel(bin++, "CPM-CMM");
    m_h_CP_errors->GetXaxis()->SetBinLabel(bin++, "CMM parity");
    m_h_CP_errors->GetXaxis()->SetBinLabel(bin++, "CMM status");
    m_h_CP_errors->GetXaxis()->SetBinLabel(bin,   "crate-sys");
  }

  return StatusCode::SUCCESS;
}

TH1D* TrigT1CaloCpmMonTool::book1D(const std::string& name,
                                   const std::string& title,
                                   int nx, double xmin, double xmax)
{
  TH1D *hist = new TH1D(TString(name), TString(title), nx, xmin, xmax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    MsgStream log(msgSvc(), this->name());
    log << MSG::WARNING << "Could not register histogram : " 
	<< name << endreq;
  }
  
  return hist;
}

TH2D* TrigT1CaloCpmMonTool::book2D(const std::string& name,
                                   const std::string& title,
                                   int nx, double xmin, double xmax,  
	                           int ny, double ymin, double ymax)
{		
  TH2D *hist = new TH2D(TString(name), TString(title), nx, xmin, xmax,
                                                       ny, ymin, ymax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    MsgStream log(msgSvc(), this->name());
    log << MSG::WARNING << "Could not register histogram : " 
	<< name << endreq;
  }
  
  return hist;
}

void TrigT1CaloCpmMonTool::newGroup(const std::string& system,
                                    LevelOfDetail_t level, Interval_t interval)
{
  if (!m_oneDir || !m_monGroup) {
    delete m_monGroup;
    std::string dir(m_rootDir + "/" + system);
    m_monGroup = new MonGroup(this, dir, level, interval);
  }
}

void TrigT1CaloCpmMonTool::setThresholdLabels(TH1* hist)
{
  hist->GetXaxis()->SetBinLabel(1, "0/1");
  hist->GetXaxis()->SetBinLabel(8, "0/8");
  hist->GetXaxis()->SetBinLabel(15, "1/1");
  hist->GetXaxis()->SetBinLabel(22, "1/8");
  hist->GetXaxis()->SetBinLabel(29, "2/1");
  hist->GetXaxis()->SetBinLabel(36, "2/8");
  hist->GetXaxis()->SetBinLabel(43, "3/1");
  hist->GetXaxis()->SetBinLabel(50, "3/8");
  hist->GetXaxis()->SetLabelSize(0.06);
  hist->GetXaxis()->SetTitleOffset(1.25);
  hist->GetXaxis()->SetNdivisions(-1404); // why doesn't this work?
}

void TrigT1CaloCpmMonTool::setStatusLabels(TH1* hist)
{
  const LVL1::DataError err(0); // should have made bitName static
  for (int bit = 0; bit < 8; ++bit) {
    hist->GetXaxis()->SetBinLabel(bit + 1,
                   (err.bitName(bit + LVL1::DataError::GLinkParity)).c_str()); 
  }
  hist->GetXaxis()->SetLabelSize(0.06);
}

void TrigT1CaloCpmMonTool::setCmmLocLabels(TH2* hist)
{
  hist->GetYaxis()->SetBinLabel(1, "0/L");
  hist->GetYaxis()->SetBinLabel(2, "0/R");
  hist->GetYaxis()->SetBinLabel(3, "1/L");
  hist->GetYaxis()->SetBinLabel(4, "1/R");
  hist->GetYaxis()->SetBinLabel(5, "2/L");
  hist->GetYaxis()->SetBinLabel(6, "2/R");
  hist->GetYaxis()->SetBinLabel(7, "3/L");
  hist->GetYaxis()->SetBinLabel(8, "3/R");
  hist->GetYaxis()->SetLabelSize(0.06);
}
