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
  : ManagedMonitorToolBase(type, name, parent)
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

  //ROOT File directory
  //declareProperty("RootFileDirectory", m_dir = "CPM" );
  declareProperty("HistogramPrefix", m_prefix = "" );
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
  
  sc = service( "StoreGateSvc", m_StoreGate);
  if( sc.isFailure() ) {
    log << MSG::ERROR << "Unable to locate Service StoreGateSvc" << endreq;
    return sc;
  }

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

  std::string rootDir("L1Calo/CPM/expert");
  if (m_prefix != "") rootDir += "/" + m_prefix;
  MonGroup monExpert ( this, rootDir, expert, eventsBlock );
  rootDir = "L1Calo/CPM/shift";
  if (m_prefix != "") rootDir += "/" + m_prefix;
  MonGroup monShift ( this, rootDir, shift, eventsBlock );

  if ( isNewEventsBlock || isNewLumiBlock ) {	

  //  Timeslice checks

  m_monGroup = &monShift;

  m_h_CPM_slices = book2D("CPM_slices",
                          "CPM Timeslices;Number of Slices;Triggered Slice",
			  6, -0.5, 5.5, 6, -0.5, 5.5);
  m_h_CMM_slices = book2D("CMM_slices",
                          "CMM Timeslices;Number of Slices;Triggered Slice",
			  6, -0.5, 5.5, 6, -0.5, 5.5);

  //  CPM Tower - Trigger Tower comparison Histos

  m_monGroup = &monExpert;

  m_h_TT_Em_Et = book1D("TT_EM_Et","Trigger Tower EM Et Noise",10,0,10);
  m_h_TT_Had_Et = book1D("TT_HAD_Et","Trigger Tower HAD Et Noise",10,0,10);
  m_h_TT_Em_Et_s = book1D("TT_EM_Et_s","Trigger Tower EM Et Signal",245,10,255);
  m_h_TT_Had_Et_s = book1D("TT_HAD_Et_s","Trigger Tower HAD Et Signal",245,10,255);
  m_h_TT_Em_eta = book1D("TT_EM_eta","Trigger Tower EM eta",50,-2.5,2.5);
  m_h_TT_Had_eta = book1D("TT_HAD_eta","Trigger Tower HAD eta",50,-2.5,2.5);
  m_h_TT_Em_phi = book1D("TT_EM_phi","Trigger Tower EM phi ",64,0,2*M_PI);
  m_h_TT_Had_phi = book1D("TT_HAD_phi","Trigger Tower HAD phi ",64,0,2*M_PI);
  m_h_TT_Em_eta_phi = book2D("TT_EM_eta_phi",
         "Trigger Tower EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_TT_Had_eta_phi = book2D("TT_HAD_eta_phi",
         "Trigger Tower HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_TT_Em_eta_phi_w = book2D("TT_EM_eta_phi_w",
         "Trigger Tower EM eta/phi weighted;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_TT_Had_eta_phi_w = book2D("TT_HAD_eta_phi_w",
         "Trigger Tower HAD eta/phi weighted;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);

  m_monGroup = &monShift;

  m_h_CT_Em_Et = book1D("CT_EM_Et","CPM Tower EM Et Noise",10,0,10);
  m_h_CT_Had_Et = book1D("CT_HAD_Et","CPM Tower HAD Et Noise",10,0,10);
  m_h_CT_Em_Et_s = book1D("CT_EM_Et_s","CPM Tower EM Et Signal",245,10,255);
  m_h_CT_Had_Et_s = book1D("CT_HAD_Et_s","CPM Tower HAD Et Signal",245,10,255);
  m_h_CT_Em_eta = book1D("CT_EM_eta","CPM Tower EM eta",50,-2.5,2.5);
  m_h_CT_Had_eta = book1D("CT_HAD_eta","CPM Tower HAD eta",50,-2.5,2.5);
  m_h_CT_Em_phi = book1D("CT_EM_phi","CPM Tower EM phi ",64,0,2*M_PI);
  m_h_CT_Had_phi = book1D("CT_HAD_phi","CPM Tower HAD phi ",64,0,2*M_PI);
  m_h_CT_Em_eta_phi = book2D("CT_EM_eta_phi",
         "CPM Tower EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_CT_Had_eta_phi = book2D("CT_HAD_eta_phi",
         "CPM Tower HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_CT_Em_eta_phi_w = book2D("CT_EM_eta_phi_w",
         "CPM Tower EM eta/phi weighted;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_CT_Had_eta_phi_w = book2D("CT_HAD_eta_phi_w",
         "CPM Tower HAD eta/phi weighted;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);

  m_h_TTCT_Em_eta_phi = book2D("TTCT_EM_eta_phi",
    "Trigger/CPM Tower mismatch EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_TTCT_Had_eta_phi = book2D("TTCT_HAD_eta_phi",
    "Trigger/CPM Tower mismatch HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);

  m_h_CT_Em_parity = book2D("CT_EM_parity",
            "CPM Tower EM Parity Errors;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_CT_Had_parity = book2D("CT_HAD_parity",
            "CPM Tower HAD Parity Errors;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_CT_Em_link = book2D("CT_EM_link",
            "CPM Tower EM Link Down Errors;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_CT_Had_link = book2D("CT_HAD_link",
            "CPM Tower HAD Link Down Errors;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_CT_status = book1D("CT_status", "CPM Sub-status bits", 8, 0., 8.);
  setStatusLabels(m_h_CT_status);
  m_h_CT_status_eta_phi = book2D("CT_status_eta_phi",
            "CPM Sub-status hit-map;eta;phi", 50,-2.5,2.5,64,0,2*M_PI);

  //  CPM Hits

  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("CH_Thresh_" + cnum.str());
    std::string title("CPM Hits Threshold " + cnum.str() + ";Crate/CPM");
    TH1D* hist = book1D(name, title, 56, 0, 56);
    setThresholdLabels(hist);
    m_v_thresholds.push_back(hist);
  }

  //  CMM-CP Hits

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
  m_h_CPM_CMM_hits0 = book1D("CPM_CMM_hits0",
                 "CPM-CMM Hits mismatch (Right);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPM_CMM_hits0);
  m_h_CPM_CMM_hits1 = book1D("CPM_CMM_hits1",
                 "CPM-CMM Hits mismatch (Left);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPM_CMM_hits1);
  m_h_CMM_R_parity = book1D("CMM_R_parity",
                 "CMM Parity Errors (Right);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CMM_R_parity);
  m_h_CMM_L_parity = book1D("CMM_L_parity",
                 "CMM Parity Errors (Left);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CMM_L_parity);
  m_h_CMM_status = book1D("CMM_status", "CMM Sub-status bits", 8, 0., 8.);
  setStatusLabels(m_h_CMM_status);

  //  CPM RoIs

  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("RoI_Thresh_" + cnum.str());
    std::string title("RoI Threshold " + cnum.str() + ";Crate/CPM");
    TH1D* hist = book1D(name, title, 56, 0, 56);
    setThresholdLabels(hist);
    m_v_RoI_thresholds.push_back(hist);
  }
  const double halfPhiBin = 2.*M_PI/128.;
  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("RoI_2D_Thresh_" + cnum.str());
    std::string title("RoI eta/phi Threshold " + cnum.str() + ";eta;phi");
    m_v_RoI_2D_thresholds.push_back(book2D(name, title, 51, -2.55, 2.55,
                                    65, -halfPhiBin, 2*M_PI+halfPhiBin));
  }

  } // end if (isNewEventsBlock ...

  if( isNewRun ) { }

  //SetBookStatus(true);
  
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloCpmMonTool::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  

  //Retrieve Trigger Towers from SG
  const TriggerTowerCollection* triggerTowerTES = 0; 
  StatusCode sc = m_StoreGate->retrieve(triggerTowerTES,
                                                     m_triggerTowerLocation); 
  if( sc.isFailure()  ||  !triggerTowerTES ) {
    log << MSG::ERROR<< "No Trigger Tower container found"<< endreq; 
  }

  //Retrieve CPM Towers from SG
  const CpmTowerCollection* cpmTowerTES = 0; 
  sc = m_StoreGate->retrieve(cpmTowerTES, m_cpmTowerLocation); 
  if( sc.isFailure()  ||  !cpmTowerTES ) {
    log << MSG::ERROR<< "No CPM Tower container found"<< endreq; 
  }
  
  //Retrieve CPM Hits from SG
  const CpmHitsCollection* cpmHitsTES = 0;
  sc = m_StoreGate->retrieve( cpmHitsTES, m_cpmHitsLocation);
  if( sc.isFailure()  ||  !cpmHitsTES ) {
    log << MSG::ERROR << "No CPM Hits container found"<< endreq; 
  }
  
  //Retrieve CMM-CP Hits from SG
  const CmmCpHitsCollection* cmmCpHitsTES = 0;
  sc = m_StoreGate->retrieve( cmmCpHitsTES, m_cmmCpHitsLocation);
  if( sc.isFailure()  ||  !cmmCpHitsTES ) {
    log << MSG::ERROR << "No CMM-CP Hits container found"<< endreq; 
  }
  
  //Retrieve CPM RoIs from SG
  const CpmRoiCollection* cpmRoiTES = 0;
  sc = m_StoreGate->retrieve( cpmRoiTES, m_cpmRoiLocation);
  if( sc.isFailure()  ||  !cpmRoiTES ) {
    log << MSG::ERROR << "No CPM RoIs container found"<< endreq; 
  }

  //=============================================
  //   CPM Tower - Trigger Tower comparison plots
  //=============================================

  // Maps for one-one comparisons
  TriggerTowerMap ttMap;
  CpmTowerMap     cpMap;
  LVL1::TriggerTowerKey towerKey;

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
      if (em && eta > -2.5 && eta < 2.5) {
        m_h_TT_Em_Et->Fill(em, 1.);
        if (em > 10) m_h_TT_Em_Et_s->Fill(em, 1.);
        m_h_TT_Em_eta->Fill(eta, 1.);
        m_h_TT_Em_phi->Fill(phi, 1.);
        m_h_TT_Em_eta_phi->Fill(eta, phi, 1.);
        m_h_TT_Em_eta_phi_w->Fill(eta, phi, em);
      }
      if (had && eta > -2.5 && eta < 2.5) {
        m_h_TT_Had_Et->Fill(had, 1.);
        if (had > 10) m_h_TT_Had_Et_s->Fill(had, 1.);
        m_h_TT_Had_eta->Fill(eta, 1.);
        m_h_TT_Had_phi->Fill(phi, 1.);
        m_h_TT_Had_eta_phi->Fill(eta, phi, 1.);
        m_h_TT_Had_eta_phi_w->Fill(eta, phi, had);
      }
      if (eta > -2.5 && eta < 2.5) {
        ttMap.insert(std::make_pair(towerKey.ttKey(phi, eta), *ttIterator));
      }
    }
  }

  if (cpmTowerTES) {
    CpmTowerCollection::const_iterator ctIterator    = cpmTowerTES->begin(); 
    CpmTowerCollection::const_iterator ctIteratorEnd = cpmTowerTES->end(); 

    for (; ctIterator != ctIteratorEnd; ++ctIterator) {
      const int    peak = (*ctIterator)->peak();
      const int    slices = ((*ctIterator)->emEnergyVec()).size();
      m_h_CPM_slices->Fill(slices, peak, 1.);
      const int    em  = (*ctIterator)->emEnergy();
      const int    had = (*ctIterator)->hadEnergy();
      const double eta = (*ctIterator)->eta();
      const double phi = (*ctIterator)->phi();
      if (em) {
        m_h_CT_Em_Et->Fill(em, 1.);
        if (em > 10) m_h_CT_Em_Et_s->Fill(em, 1.);
        m_h_CT_Em_eta->Fill(eta, 1.);
        m_h_CT_Em_phi->Fill(phi, 1.);
        m_h_CT_Em_eta_phi->Fill(eta, phi, 1.);
        m_h_CT_Em_eta_phi_w->Fill(eta, phi, em);
      }
      if (had) {
        m_h_CT_Had_Et->Fill(had, 1.);
        if (had > 10) m_h_CT_Had_Et_s->Fill(had, 1.);
        m_h_CT_Had_eta->Fill(eta, 1.);
        m_h_CT_Had_phi->Fill(phi, 1.);
        m_h_CT_Had_eta_phi->Fill(eta, phi, 1.);
        m_h_CT_Had_eta_phi_w->Fill(eta, phi, had);
      }
      // Errors
      const LVL1::DataError emError((*ctIterator)->emError());
      const LVL1::DataError hadError((*ctIterator)->hadError());
      m_h_CT_Em_parity->Fill(eta, phi, emError.get(LVL1::DataError::Parity));
      m_h_CT_Had_parity->Fill(eta, phi, hadError.get(LVL1::DataError::Parity));
      m_h_CT_Em_link->Fill(eta, phi, emError.get(LVL1::DataError::LinkDown));
      m_h_CT_Had_link->Fill(eta, phi, hadError.get(LVL1::DataError::LinkDown));
      // Sub-status errors
      const int status = emError.error() >> LVL1::DataError::GLinkParity;
      for (int bit = 0; bit < 8; ++bit) {
        m_h_CT_status->Fill(bit, (status >> bit) & 0x1);
      }
      m_h_CT_status_eta_phi->Fill(eta, phi, status != 0);

      cpMap.insert(std::make_pair(towerKey.ttKey(phi, eta), *ctIterator));
    }
  }

  // One-to-one tower comparison

  TriggerTowerMap::const_iterator ttMapIter    = ttMap.begin();
  TriggerTowerMap::const_iterator ttMapIterEnd = ttMap.end();
  CpmTowerMap::const_iterator     cpMapIter    = cpMap.begin();
  CpmTowerMap::const_iterator     cpMapIterEnd = cpMap.end();

  while (ttMapIter != ttMapIterEnd || cpMapIter != cpMapIterEnd) {

    unsigned int ttKey = 0xffffffff;
    unsigned int cpKey = 0xffffffff;
    int ttEm  = 0;
    int ttHad = 0;
    int cpEm  = 0;
    int cpHad = 0;
    double eta = 0.;
    double phi = 0.;

    if (ttMapIter != ttMapIterEnd) ttKey = ttMapIter->first;
    if (cpMapIter != cpMapIterEnd) cpKey = cpMapIter->first;

    if ((cpMapIter == cpMapIterEnd) || (cpKey > ttKey)) {

      // TriggerTower but no CPMTower

      const LVL1::TriggerTower* tt = ttMapIter->second;
      ttEm  = tt->emEnergy();
      ttHad = tt->hadEnergy();
      eta = tt->eta();
      phi = tt->phi();
      ++ttMapIter;

    } else if ((ttMapIter == ttMapIterEnd) || (ttKey > cpKey)) {

      // CPMTower but no TriggerTower

      const LVL1::CPMTower* cp = cpMapIter->second;
      cpEm  = cp->emEnergy();
      cpHad = cp->hadEnergy();
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
      eta = tt->eta();
      phi = tt->phi();
      ++ttMapIter;
      ++cpMapIter;
    }
    m_h_TTCT_Em_eta_phi->Fill(eta, phi, ttEm != cpEm);
    m_h_TTCT_Had_eta_phi->Fill(eta, phi, ttHad != cpHad);
    if ((ttEm != cpEm) || (ttHad != cpHad)) {
      log << MSG::DEBUG
          << "Trigger/CPM Tower mismatch, eta/phi/ttEm/ttHad/cpEm/cpHad: "
	  << eta << "/" << phi << "/" << ttEm << "/" << ttHad << "/"
	  << cpEm << "/" << cpHad << endreq;
    }
  }

  //=============================================
  //  CPM Hits
  //=============================================

  CpmHitsMap cpmMap;

  if (cpmHitsTES) {
    CpmHitsCollection::const_iterator chIterator    = cpmHitsTES->begin(); 
    CpmHitsCollection::const_iterator chIteratorEnd = cpmHitsTES->end(); 
    for (; chIterator != chIteratorEnd; ++chIterator) {
      const int peak   = (*chIterator)->peak();
      const int slices = ((*chIterator)->HitsVec0()).size();
      m_h_CPM_slices->Fill(slices, peak, 1.);
      const unsigned int hits0 = (*chIterator)->HitWord0();
      const unsigned int hits1 = (*chIterator)->HitWord1();
      const int crate = (*chIterator)->crate();
      const int cpm   = (*chIterator)->module();
      const int bin   = crate * 14 + cpm - 1;
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
      const int key = crate * 100 + cpm;
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
      const int peak   = (*cmIterator)->peak();
      const int slices = ((*cmIterator)->HitsVec0()).size();
      m_h_CMM_slices->Fill(slices, peak, 1.);
      const unsigned int hits0 = (*cmIterator)->HitWord0();
      const unsigned int hits1 = (*cmIterator)->HitWord1();
      const int crate  = (*cmIterator)->crate();
      const int dataId = (*cmIterator)->dataID();
      const int bin    = (dataId < 15) ? crate * 14 + dataId - 1
  		                       : crate * 5  + dataId - 15; 
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
      }

      if (dataId < 15) {
        const int key = crate * 100 + dataId;
        cmmMap.insert(std::make_pair(key, *cmIterator));
      }
    }
  }

  // One-to-one CPM-CMM hit comparison

  CpmHitsMap::const_iterator   cpmMapIter    = cpmMap.begin();
  CpmHitsMap::const_iterator   cpmMapIterEnd = cpmMap.end();
  CmmCpHitsMap::const_iterator cmmMapIter    = cmmMap.begin();
  CmmCpHitsMap::const_iterator cmmMapIterEnd = cmmMap.end();

  while (cpmMapIter != cpmMapIterEnd || cmmMapIter != cmmMapIterEnd) {

    int cpmKey = 999;
    int cmmKey = 999;
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
      ++cpmMapIter;
      ++cmmMapIter;
    }
    const int bin = crate * 14 + cpm - 1;
    m_h_CPM_CMM_hits0->Fill(bin, cpmHits0 != cmmHits0);
    m_h_CPM_CMM_hits1->Fill(bin, cpmHits1 != cmmHits1);
    if (cpmHits0 != cmmHits0 || cpmHits1 != cmmHits1) {
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
	  (*hist2)->Fill(eta, phi, 1.);
        }
        ++hist1;
        ++hist2;
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

  if (isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun) { }

  return StatusCode::SUCCESS;
}

TH1D* TrigT1CaloCpmMonTool::book1D(std::string nam, std::string tit,
                                   int nx, double xmin, double xmax)
{
  const std::string newName = (m_prefix == "") ? nam : m_prefix + "_" + nam;
  const std::string newTitle = (m_prefix == "") ? tit : m_prefix + " " + tit;
  TH1D *hist = new TH1D(TString(newName), TString(newTitle), nx, xmin, xmax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    MsgStream log(msgSvc(), name());
    log << MSG::WARNING << "Could not register histogram : " 
	<< newName << endreq;
  }
  
  return hist;
}

TH2D* TrigT1CaloCpmMonTool::book2D(std::string nam, std::string tit,
                                   int nx, double xmin, double xmax,  
	                           int ny, double ymin, double ymax)
{		
  const std::string newName = (m_prefix == "") ? nam : m_prefix + "_" + nam;
  const std::string newTitle = (m_prefix == "") ? tit : m_prefix + " " + tit;
  TH2D *hist = new TH2D(TString(newName), TString(newTitle), nx, xmin, xmax,
                                                             ny, ymin, ymax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    MsgStream log(msgSvc(), name());
    log << MSG::WARNING << "Could not register histogram : " 
	<< newName << endreq;
  }
  
  return hist;
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
