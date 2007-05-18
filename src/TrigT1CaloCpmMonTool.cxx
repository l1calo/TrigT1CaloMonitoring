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
#include "TrigT1Calo/TriggerTower.h"
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

  MonGroup monExpert ( this, "CPM/expert", expert, eventsBlock );

  if ( isNewEventsBlock || isNewLumiBlock ) {	

  m_monGroup = &monExpert;
  
  //  CPM Tower - Trigger Tower comparison Histos

  m_h_TT_Em_Et = book1D("TT_EM_Et","Trigger Tower EM Et Noise",10,0,10);
  m_h_TT_Had_Et = book1D("TT_HAD_Et","Trigger Tower HAD Et Noise",10,0,10);
  m_h_CT_Em_Et = book1D("CT_EM_Et","CPM Tower EM Et Noise",10,0,10);
  m_h_CT_Had_Et = book1D("CT_HAD_Et","CPM Tower HAD Et Noise",10,0,10);
  m_h_TT_Em_Et_s = book1D("TT_EM_Et_s","Trigger Tower EM Et Signal",245,10,255);
  m_h_TT_Had_Et_s = book1D("TT_HAD_Et_s","Trigger Tower HAD Et Signal",245,10,255);
  m_h_CT_Em_Et_s = book1D("CT_EM_Et_s","CPM Tower EM Et Signal",245,10,255);
  m_h_CT_Had_Et_s = book1D("CT_HAD_Et_s","CPM Tower HAD Et Signal",245,10,255);
  m_h_TT_Em_eta = book1D("TT_EM_eta","Trigger Tower EM eta",50,-2.5,2.5);
  m_h_TT_Had_eta = book1D("TT_HAD_eta","Trigger Tower HAD eta",50,-2.5,2.5);
  m_h_CT_Em_eta = book1D("CT_EM_eta","CPM Tower EM eta",50,-2.5,2.5);
  m_h_CT_Had_eta = book1D("CT_HAD_eta","CPM Tower HAD eta",50,-2.5,2.5);
  m_h_TT_Em_phi = book1D("TT_EM_phi","Trigger Tower EM phi ",64,0,2*M_PI);
  m_h_TT_Had_phi = book1D("TT_HAD_phi","Trigger Tower HAD phi ",64,0,2*M_PI);
  m_h_CT_Em_phi = book1D("CT_EM_phi","CPM Tower EM phi ",64,0,2*M_PI);
  m_h_CT_Had_phi = book1D("CT_HAD_phi","CPM Tower HAD phi ",64,0,2*M_PI);
  m_h_TT_Em_eta_phi = book2D("TT_EM_eta_phi","Trigger Tower EM eta/phi",
                                                50,-2.5,2.5,64,0,2*M_PI);
  m_h_TT_Had_eta_phi = book2D("TT_HAD_eta_phi","Trigger Tower HAD eta/phi",
                                                50,-2.5,2.5,64,0,2*M_PI);
  m_h_CT_Em_eta_phi = book2D("CT_EM_eta_phi","CPM Tower EM eta/phi",
                                                50,-2.5,2.5,64,0,2*M_PI);
  m_h_CT_Had_eta_phi = book2D("CT_HAD_eta_phi","CPM Tower HAD eta/phi",
                                                50,-2.5,2.5,64,0,2*M_PI);
  m_h_TT_Em_eta_phi_w = book2D("TT_EM_eta_phi_w",
               "Trigger Tower EM eta/phi weighted", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_TT_Had_eta_phi_w = book2D("TT_HAD_eta_phi_w",
               "Trigger Tower HAD eta/phi weighted", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_CT_Em_eta_phi_w = book2D("CT_EM_eta_phi_w",
                    "CPM Tower EM eta/phi weighted", 50,-2.5,2.5,64,0,2*M_PI);
  m_h_CT_Had_eta_phi_w = book2D("CT_HAD_eta_phi_w",
                    "CPM Tower HAD eta/phi weighted", 50,-2.5,2.5,64,0,2*M_PI);

  //  CPM Hits

  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("CH_Thresh_" + cnum.str());
    std::string title("CPM Hits Threshold " + cnum.str());
    m_v_thresholds.push_back(book1D(name, title, 56, 0, 56));
  }

  //  CMM-CP Hits

  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("CM_Thresh_" + cnum.str());
    std::string title("CMM-CP Hits Threshold " + cnum.str());
    m_v_CMM_thresholds.push_back(book1D(name, title, 56, 0, 56));
  }
  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("CM_T_Thresh_" + cnum.str());
    std::string title("CMM-CP Hits Totals Threshold " + cnum.str());
    m_v_CMM_T_thresholds.push_back(book1D(name, title, 20, 0, 20));
  }

  //  CPM RoIs

  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("RoI_Thresh_" + cnum.str());
    std::string title("RoI Threshold " + cnum.str());
    m_v_RoI_thresholds.push_back(book1D(name, title, 56, 0, 56));
  }
  const double halfPhiBin = 2.*M_PI/128.;
  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("RoI_2D_Thresh_" + cnum.str());
    std::string title("RoI eta/phi Threshold " + cnum.str());
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
    }
  }

  //=============================================
  //  CPM Hits
  //=============================================

  if (cpmHitsTES) {
    CpmHitsCollection::const_iterator chIterator    = cpmHitsTES->begin(); 
    CpmHitsCollection::const_iterator chIteratorEnd = cpmHitsTES->end(); 
    for (; chIterator != chIteratorEnd; ++chIterator) {
      const unsigned int hits0 = (*chIterator)->HitWord0();
      const unsigned int hits1 = (*chIterator)->HitWord1();
      const int bin = (*chIterator)->crate() * 14 + (*chIterator)->module() - 1;
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
    }
  }

  //=============================================
  //  CMM-CP Hits
  //=============================================

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
