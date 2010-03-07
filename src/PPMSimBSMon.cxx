// ********************************************************************
//
// NAME:     PPMSimBSMon.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHORS:  Peter Faulkner
//           Sky French          
//
// ********************************************************************

#include <sstream>
#include <utility>
#include <cmath>
#include "SGTools/StlVectorClids.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TString.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IToolSvc.h"
#include "StoreGate/StoreGateSvc.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"

#include "TrigT1CaloUtils/CoordToHardware.h"
#include "TrigT1CaloEvent/TriggerTower.h"
#include "TrigT1CaloUtils/TriggerTowerKey.h"
#include "TrigT1CaloToolInterfaces/IL1TriggerTowerTool.h"
#include "TrigT1Interfaces/Coordinate.h"
#include "TrigT1Interfaces/CoordinateRange.h"
#include "TrigT1Interfaces/TrigT1CaloDefs.h"

#include "TrigT1CaloMonitoring/PPMSimBSMon.h"

#include "Identifier/HWIdentifier.h"

/*---------------------------------------------------------*/
PPMSimBSMon::PPMSimBSMon(const std::string & type, 
			 const std::string & name,
			 const IInterface* parent)
  : ManagedMonitorToolBase(type, name, parent),
    m_storeGate("StoreGateSvc", name),
    m_ttTool("LVL1::L1TriggerTowerTool/L1TriggerTowerTool"),
    m_log(msgSvc(), name), m_debug(false),
    m_monGroup(0), m_events(0)
/*---------------------------------------------------------*/
{
  declareInterface<IMonitorToolBase>(this); 

  declareProperty("TriggerTowerLocation",
                 m_triggerTowerLocation =
		  LVL1::TrigT1CaloDefs::TriggerTowerLocation);
  
  declareProperty("RootDirectory", m_rootDir = "L1Calo");

  declareProperty("EventSamples", m_eventSamples = 10,
                  "Number of Error Event Number Samples");

}

/*---------------------------------------------------------*/
PPMSimBSMon::~PPMSimBSMon()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode PPMSimBSMon:: initialize()
/*---------------------------------------------------------*/
{
  m_log.setLevel(outputLevel());
  m_debug = outputLevel() <= MSG::DEBUG;

  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;

  sc = m_storeGate.retrieve();
  if( sc.isFailure() ) {
    m_log << MSG::ERROR << "Unable to locate Service StoreGateSvc" << endreq;
    return sc;
  }

  sc = m_ttTool.retrieve();
  if( sc.isFailure() ) {
    m_log << MSG::ERROR << "Unable to locate Tool L1TriggerTowerTool" << endreq;
    return sc;
  }

  ///========
  IToolSvc* toolSvc;
  StatusCode status = service( "ToolSvc",toolSvc  );

  if(status.isSuccess()) 
    {
      IAlgTool *algtool;
      sc = toolSvc->retrieveTool("L1CaloTTIdTools", algtool);
      m_log<<MSG::DEBUG<<"L1CaloTTIdTools retrieved"<<endreq;
      if (sc!=StatusCode::SUCCESS) 
	{
	  m_log << MSG::ERROR << " Cannot get L1CaloTTIdTools !" << endreq;
	  return sc;
	}
      m_l1CaloTTIdTools = dynamic_cast<L1CaloTTIdTools*> (algtool);
    } 
  else 
    {
      return StatusCode::FAILURE;
    }

  ISvcLocator* svcLoc = Gaudi::svcLocator( );
  toolSvc = 0; // Pointer to Tool Service
  sc = svcLoc->service( "ToolSvc",toolSvc  );
  if(sc.isSuccess())
    {
      sc = toolSvc->retrieveTool("CaloTriggerTowerService",m_ttSvc);
      if(sc.isFailure())
        {
          m_log << MSG::ERROR << "Could not retrieve CaloTriggerTowerService Tool" << endreq;
          return StatusCode::FAILURE;
        }
    }
  else
    {
      m_log << MSG::ERROR << "Could not retrieve ToolSvc" << endreq;
      return StatusCode::FAILURE;
    }
  
  // Get a pointer to DetectorStore services
  sc = service("DetectorStore", m_detStore);
  if (sc.isFailure())
    {
      m_log << MSG::ERROR << "Cannot access DetectorStore" << endreq;
      return StatusCode::FAILURE;
    }
  
  // Retrieve the CaloIdManager from the detector store
  sc = m_detStore->retrieve(m_caloMgr);
  if (sc.isFailure())
    {
      m_log << MSG::ERROR << "Unable to retrieve CaloIdManager from DetectorStore" << endreq;
      return StatusCode::FAILURE;
    }

  // Use the CaloIdManager to get a pointer to an instance of the CaloLVL1_ID helper
  m_lvl1Helper = m_caloMgr->getLVL1_ID();
  if (!m_lvl1Helper)
    {
      m_log << MSG::ERROR << "Could not access CaloLVL1_ID helper" << endreq;
      return StatusCode::FAILURE;
    }

  // Use the CaloIdManager to get a pointer to an instance of the TTOnlineID helper
  m_l1ttonlineHelper = m_caloMgr->getTTOnlineID();
  if (!m_l1ttonlineHelper )
    {
      m_log << MSG::ERROR << "Could not access TTOnlineID helper" << endreq;
      return StatusCode::FAILURE;
    }

  

  //===========

  // Phi units
  const double twoPi = 2.*M_PI;
  m_phiMax = 64;
  m_phiScale = m_phiMax/twoPi;

  m_log << MSG::INFO << "PPMSimBSMon initialised" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode PPMSimBSMon:: finalize()
/*---------------------------------------------------------*/
{
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode PPMSimBSMon::bookHistograms(bool isNewEventsBlock,
                                       bool isNewLumiBlock, bool isNewRun)
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "bookHistograms entered" << endreq;

  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
  	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  if ( isNewEventsBlock || isNewLumiBlock ) { }
  
  if( isNewRun ) {

  std::string dir(m_rootDir + "/PPM/Errors/Data_Simulation");
  MonGroup monShift ( this, dir, shift, run );
  MonGroup monExpert( this, dir, expert, run );
  MonGroup monPPM   ( this, dir + "/PPMLUTSim", expert, run );
  MonGroup monEvent ( this, dir + "/MismatchEventNumbers", expert, run, "", "eventSample" );

  // LUT

  m_monGroup = &monPPM;

  m_h_ppm_em_2d_etaPhi_tt_lut_SimEqData = bookEtaPhi("ppm_em_2d_etaPhi_tt_lut_SimEqData",
						     "PPM LUT EM Data/Simulation Non-zero Matches;eta;phi");
  m_h_ppm_em_2d_etaPhi_tt_lut_SimNeData = bookEtaPhi("ppm_em_2d_etaPhi_tt_lut_SimNeData",
						     "PPM LUT EM Data/Simulation Non-zero Mismatches;eta;phi");
  m_h_ppm_em_2d_etaPhi_tt_lut_SimNoData = bookEtaPhi("ppm_em_2d_etaPhi_tt_lut_SimNoData",
						     "PPM LUT EM Simulation but no Data;eta;phi");
  m_h_ppm_em_2d_etaPhi_tt_lut_DataNoSim = bookEtaPhi("ppm_em_2d_etaPhi_tt_lut_DataNoSim",
						     "PPM LUT EM Data but no Simulation;eta;phi");
  m_h_ppm_had_2d_etaPhi_tt_lut_SimEqData = bookEtaPhi("ppm_had_2d_etaPhi_tt_lut_SimEqData",
						      "PPM LUT HAD Data/Simulation Non-zero Matches;eta;phi");
  m_h_ppm_had_2d_etaPhi_tt_lut_SimNeData = bookEtaPhi("ppm_had_2d_etaPhi_tt_lut_SimNeData",
						      "PPM LUT HAD Data/Simulation Non-zero Mismatches;eta;phi");
  m_h_ppm_had_2d_etaPhi_tt_lut_SimNoData = bookEtaPhi("ppm_had_2d_etaPhi_tt_lut_SimNoData",
						      "PPM LUT HAD Simulation but no Data;eta;phi");
  m_h_ppm_had_2d_etaPhi_tt_lut_DataNoSim = bookEtaPhi("ppm_had_2d_etaPhi_tt_lut_DataNoSim",
						      "PPM LUT HAD Data but no Simulation;eta;phi");

  // Mismatch Histograms

  m_monGroup = &monEvent;

  m_sampleCounts.clear();
  m_sampleCounts.resize(8*32,0);

  m_h_ppm_2d_LUT_MismatchEvents_cr0cr1 = book2I("ppm_2d_LUT_MismatchEvents_cr0cr1","PPM LUT Mismatch Event Numbers;Sample;Crate/Module",m_eventSamples,0,m_eventSamples,32,0,32);
  m_h_ppm_2d_LUT_MismatchEvents_cr2cr3 = book2I("ppm_2d_LUT_MismatchEvents_cr2cr3","PPM LUT Mismatch Event Numbers;Sample;Crate/Module",m_eventSamples,0,m_eventSamples,32,0,32);
  m_h_ppm_2d_LUT_MismatchEvents_cr4cr5 = book2I("ppm_2d_LUT_MismatchEvents_cr4cr5","PPM LUT Mismatch Event Numbers;Sample;Crate/Module",m_eventSamples,0,m_eventSamples,32,0,32);
  m_h_ppm_2d_LUT_MismatchEvents_cr6cr7 = book2I("ppm_2d_LUT_MismatchEvents_cr6cr7","PPM LUT Mismatch Event Numbers;Sample;Crate/Module",m_eventSamples,0,m_eventSamples,32,0,32);

  setLabelsCM(m_h_ppm_2d_LUT_MismatchEvents_cr0cr1,false,0);
  setLabelsCM(m_h_ppm_2d_LUT_MismatchEvents_cr2cr3,false,2);
  setLabelsCM(m_h_ppm_2d_LUT_MismatchEvents_cr4cr5,false,4);
  setLabelsCM(m_h_ppm_2d_LUT_MismatchEvents_cr6cr7,false,6);

  } // end if (isNewRun ...

  m_log << MSG::DEBUG << "Leaving bookHistograms" << endreq;
  
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode PPMSimBSMon::fillHistograms()
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "fillHistograms entered" << endreq;

  // NB. Collection retrieves wrapped in m_storeGate->contains<..>(..)
  // are for those not expected to be on ESD. They should be on bytestream.
  
  StatusCode sc;

  //Retrieve Trigger Towers from SG
  const TriggerTowerCollection* triggerTowerTES = 0; 
  sc = m_storeGate->retrieve(triggerTowerTES, m_triggerTowerLocation); 
  if( sc.isFailure()  ||  !triggerTowerTES ) {
    m_log << MSG::DEBUG<< "No Trigger Tower container found"<< endreq; 
  }
  ++m_events;

  // Compare LUT simulated from FADC with LUT from data

  if (triggerTowerTES) {
    simulateAndCompare(triggerTowerTES);
  }
 
  m_log << MSG::DEBUG << "Leaving fillHistograms" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode PPMSimBSMon::procHistograms(bool isEndOfEventsBlock,
                                  bool isEndOfLumiBlock, bool isEndOfRun)
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "procHistograms entered" << endreq;

  if (isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun) {
  }

  return StatusCode::SUCCESS;
}

TH1F* PPMSimBSMon::book1F(std::string name, std::string title,
                                   int nx, double xmin, double xmax)
{
  TH1F *hist = new TH1F(TString(name), TString(title), nx, xmin, xmax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    m_log << MSG::WARNING << "Could not register histogram : " 
 	  << name << endreq;
  }
  hist->SetStats(kFALSE);
  
  return hist;
}

TH2F* PPMSimBSMon::bookEtaPhi(std::string name, std::string title)
{
  TH2F* hist = 0;
  const double phiBin = M_PI/32.;
  const double halfPhiBin = M_PI/64.;
  const int nxbins = 66;

  const double xbins[nxbins+1] = {-4.9,-4.475,-4.050,-3.625,-3.2,-3.1,-2.9,-2.7,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.7,2.9,3.1,3.2,3.625,4.050,4.475,4.9};
  std::string newTitle = title + ";eta";
  hist = book2F(name, newTitle, nxbins, xbins, 64,0.,64);
  for (int chan = 0; chan < 64; chan +=2 ) {
    const double rad = chan*phiBin + halfPhiBin;
    std::ostringstream cnum;
    cnum << chan << "/"
        << std::setiosflags(std::ios::fixed | std::ios::showpoint)
        << std::setprecision(2) << rad;
    hist->GetYaxis()->SetBinLabel(chan+1, cnum.str().c_str());
  }
  hist->GetYaxis()->SetTitleOffset(1.3);

  return hist;

}

TH2F* PPMSimBSMon::book2F(std::string name, std::string title,
			  int nx, const double* xbins,
			  int ny, double ymin, double ymax)
{		
  TH2F *hist = new TH2F(TString(name), TString(title), nx, xbins,
                                                       ny, ymin, ymax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    m_log << MSG::WARNING << "Could not register histogram : " 
	  << name << endreq;
  }
  hist->SetOption("colz");
  hist->SetStats(kFALSE);
  
  return hist;
}

TH2I* PPMSimBSMon::book2I(std::string name, std::string title,
                                   int nx, double xmin, double xmax,  
	                           int ny, double ymin, double ymax)
{		
  TH2I *hist = new TH2I(TString(name), TString(title), nx, xmin, xmax,
                                                       ny, ymin, ymax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    m_log << MSG::WARNING << "Could not register histogram : " 
	  << name << endreq;
  }
  hist->SetOption("text");
  hist->SetStats(kFALSE);
  
  return hist;
}

void PPMSimBSMon::simulateAndCompare(const TriggerTowerCollection* ttIn)
{
  m_log << MSG::DEBUG << "Simulate LUT data from FADC data" << endreq;

  StatusCode sc = m_ttTool->retrieveConditions();
  if (sc.isFailure()) return;

  int nCrates = 8;
  ErrorVector crateError(nCrates);
  ErrorVector moduleError(nCrates);
  
  m_ttTool->setDebug(false);
  m_eventNumber = -1;
  TriggerTowerCollection::const_iterator iter  = ttIn->begin();
  TriggerTowerCollection::const_iterator iterE = ttIn->end();

  for (; iter != iterE; ++iter) {
    
    Identifier EmTowerId,HadTowerId;
    
    const LVL1::TriggerTower* tt = *iter;
 
    std::vector<int> emLut;
    std::vector<int> emBcidR;
    std::vector<int> emBcidD;

    int detside = m_l1CaloTTIdTools->pos_neg_z(tt->eta());
    int detregion = m_l1CaloTTIdTools->regionIndex(tt->eta());
    int eta_index = m_l1CaloTTIdTools->etaIndex(tt->eta());
    int phi_index = m_l1CaloTTIdTools->phiIndex(tt->eta(),tt->phi());
    
    EmTowerId = m_lvl1Helper->tower_id(detside, 0, detregion,eta_index,phi_index );
    HadTowerId = m_lvl1Helper->tower_id(detside, 1, detregion,eta_index,phi_index );
    
    int had_crate = -1;
    int had_module = -1;
    int em_crate = -1;
    int em_module = -1;
    
    try
      {
    HWIdentifier ttOnlId = m_ttSvc->createTTChannelID(HadTowerId);
    had_crate     = m_l1ttonlineHelper->crate(ttOnlId);
    had_module    = m_l1ttonlineHelper->module(ttOnlId);
      }
        
    catch(CaloID_Exception& except) {
      m_log << MSG::ERROR<< "CaloID_Exception" << (std::string) except << endreq;
    }
    
    try
      {
    HWIdentifier ttOnlId = m_ttSvc->createTTChannelID(EmTowerId);
    em_crate     = m_l1ttonlineHelper->crate(ttOnlId);
    em_module    = m_l1ttonlineHelper->module(ttOnlId);
      }

    catch(CaloID_Exception& except) {
      m_log << MSG::ERROR<< "CaloID_Exception" << (std::string) except << endreq;
     }

    m_ttTool->process(tt->emADC(), tt->eta(), tt->phi(), 0, emLut, emBcidR, emBcidD);
    const int emPeak = tt->emADCPeak();
    std::vector<int> emLut1;
    int emSlices = (tt->emADC()).size();
    if (emSlices < 7 || emBcidD[emPeak]) emLut1.push_back(emLut[emPeak]);
    else                 emLut1.push_back(0);
    std::vector<int> emBcidR1;
    emBcidR1.push_back(emBcidR[emPeak]);
    if (m_debug && emLut1[0] != tt->emEnergy() && (emSlices>=7 || tt->emEnergy()!=0)) { // mismatch - repeat with debug on
      std::vector<int> emLut2; 
      std::vector<int> emBcidR2;
      std::vector<int> emBcidD2;
      m_ttTool->setDebug(true);
      m_ttTool->process(tt->emADC(), tt->eta(), tt->phi(), 0, emLut2, emBcidR2, emBcidD2);
      m_ttTool->setDebug(false);
    }
    
    std::vector<int> hadLut;
    std::vector<int> hadBcidR;
    std::vector<int> hadBcidD;
    m_ttTool->process(tt->hadADC(), tt->eta(), tt->phi(), 1, hadLut, hadBcidR, hadBcidD);
    const int hadPeak = tt->hadADCPeak();
    std::vector<int> hadLut1;
    int hadSlices = (tt->hadADC()).size();
    if (hadSlices < 7 || hadBcidD[hadPeak]) hadLut1.push_back(hadLut[hadPeak]);
    else                   hadLut1.push_back(0);
    std::vector<int> hadBcidR1;
    hadBcidR1.push_back(hadBcidR[hadPeak]);
    if (m_debug && hadLut1[0] != tt->hadEnergy() && (hadSlices>=7 || tt->hadEnergy()!=0)) {
      std::vector<int> hadLut2;
      std::vector<int> hadBcidR2;
      std::vector<int> hadBcidD2;
      m_ttTool->setDebug(true);
      m_ttTool->process(tt->hadADC(), tt->eta(), tt->phi(), 1, hadLut2, hadBcidR2, hadBcidD2);
      m_ttTool->setDebug(false);
    }
    
    const int simEm  = emLut1[0];
    const int simHad = hadLut1[0];
    const int datEm  = tt->emEnergy();
    const int datHad = tt->hadEnergy();
    if (!simEm && !simHad && !datEm && !datHad) continue;
 
    const double eta = tt->eta();
    const double phi = tt->phi();
    const double phiMod = phi * m_phiScale;
    
    //  Fill in error plots
    
    int em_mismatch = 0;
    int had_mismatch = 0;
    
    TH2F* hist1 = 0;
    if (simEm && simEm == datEm) { // non-zero match
      hist1 = m_h_ppm_em_2d_etaPhi_tt_lut_SimEqData;
    } else if (simEm != datEm) {  // mis-match
      em_mismatch = 1;
      if (simEm && datEm) {       // non-zero mis-match
        hist1 = m_h_ppm_em_2d_etaPhi_tt_lut_SimNeData;
      } else if (!datEm) {        // no data
	if(emSlices>=7) {
	  hist1 = m_h_ppm_em_2d_etaPhi_tt_lut_SimNoData;
	} else em_mismatch = 0;
      } else {                    // no sim
	hist1 = m_h_ppm_em_2d_etaPhi_tt_lut_DataNoSim;
      }
      if (m_debug) {
        m_log << MSG::DEBUG << " EMTowerMismatch eta/phi/sim/dat: "
              << eta << "/" << phi << "/" << simEm << "/" << datEm << endreq;
      }
    }
    
    if (hist1) hist1->Fill(eta, phiMod);
    
    if(em_mismatch==1) {
      crateError[em_crate]=1;
      if (!((moduleError[em_crate]>>em_module)&0x1)) {
	fillEventSample(em_crate,em_module);
	moduleError[em_crate] |= (1 << em_module);
      }
    }
    
    hist1 = 0;
    if (simHad && simHad == datHad) { // non-zero match
      hist1 = m_h_ppm_had_2d_etaPhi_tt_lut_SimEqData;
    } else if (simHad != datHad) {   // mis-match
      had_mismatch = 1;
      if (simHad && datHad) {        // non-zero mis-match
        hist1 = m_h_ppm_had_2d_etaPhi_tt_lut_SimNeData;
      } else if (!datHad) {          // no data
	if(hadSlices>=7) {
	  hist1 = m_h_ppm_had_2d_etaPhi_tt_lut_SimNoData;
	} else had_mismatch = 0;
      } else {                       // no sim
	hist1 = m_h_ppm_had_2d_etaPhi_tt_lut_DataNoSim;
      }
      if (m_debug) {
        m_log << MSG::DEBUG << " HadTowerMismatch eta/phi/sim/dat: "
              << eta << "/" << phi << "/" << simHad << "/" << datHad << endreq;
      }
    }

    if (hist1) hist1->Fill(eta, phiMod);
      
    if(had_mismatch==1) {
      crateError[had_crate] = 1;
      if (!((moduleError[had_crate]>>had_module)&0x1)) {
	fillEventSample(had_crate,had_module);
	moduleError[had_crate] |= (1 << had_module);
      }
    }
  
  }    
    
  ErrorVector* save = new ErrorVector(crateError);
  sc = m_storeGate->record(save, "L1CaloPPMMismatchVector");
  if (sc.isFailure()) {
    m_log << MSG::ERROR << "Error recording PPM mismatch vector in TES "
	  << endreq;
  }
  
  m_ttTool->setDebug(true);
  
}

void PPMSimBSMon::setLabelsCM(TH2* hist, bool xAxis, int first)
{
  const int nPPMmodulespercrate = 16;
  for (int crate = 0; crate < 2; ++crate) {
    for (int module = 0; module < nPPMmodulespercrate; module += 2) {
      std::ostringstream cnum;
      cnum << crate+first << "/" << module;
      if (xAxis) {
        hist->GetXaxis()->SetBinLabel(crate*nPPMmodulespercrate + module + 1,
                                                        cnum.str().c_str());
      } else {
        hist->GetYaxis()->SetBinLabel(crate*nPPMmodulespercrate + module + 1,
                                                        cnum.str().c_str());
      }
    }
  }
}

void PPMSimBSMon::fillEventSample(int crate, int module)
{
  
  if (m_eventNumber < 0) {
    m_eventNumber = 0;
    const EventInfo* evInfo = 0;
    StatusCode sc = m_storeGate->retrieve(evInfo);
    if (sc.isFailure()) {
      m_log << MSG::DEBUG << "No EventInfo found" << endreq;
    } else {
      const EventID* evID = evInfo->event_ID();
      if (evID) m_eventNumber = evID->event_number();
    }
  }
  if (m_eventNumber > 0) {
    int offset = 0;
    int y = 0;
    int x = 0;
    int count = 0;
    if(crate %2 == 0) offset = 0;
    else offset = 16;
    y = (module-5) + offset;
    count = crate*16 + (module-5);
    x = m_sampleCounts[count];
    if (x < m_eventSamples) {
      if(crate==0 || crate==1) {
	TH2I* hist = m_h_ppm_2d_LUT_MismatchEvents_cr0cr1;
	if (hist) hist->Fill(x,y, m_eventNumber);
      }
      if(crate==2 || crate==3) {
	TH2I* hist = m_h_ppm_2d_LUT_MismatchEvents_cr2cr3;
	if (hist) hist->Fill(x,y, m_eventNumber);
      }
      if(crate==4 || crate==5) {
	TH2I* hist = m_h_ppm_2d_LUT_MismatchEvents_cr4cr5;
	if (hist) hist->Fill(x,y, m_eventNumber);
      }
      if(crate==6 || crate==7) {
	TH2I* hist = m_h_ppm_2d_LUT_MismatchEvents_cr6cr7;
	if (hist) hist->Fill(x,y, m_eventNumber);
      }
      ++m_sampleCounts[count];
    }
  }

}


