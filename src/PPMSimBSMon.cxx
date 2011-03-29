// ********************************************************************
//
// NAME:     PPMSimBSMon.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHORS:  Peter Faulkner
//           Sky French          
//
// ********************************************************************

#include <cmath>

#include "LWHists/TH2F_LW.h"
#include "LWHists/TH2I_LW.h"
#include "LWHists/TProfile2D_LW.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/StatusCode.h"
#include "SGTools/StlVectorClids.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigT1CaloEvent/TriggerTower.h"
#include "TrigT1CaloUtils/TriggerTowerKey.h"
#include "TrigT1CaloToolInterfaces/IL1TriggerTowerTool.h"
#include "TrigT1Interfaces/TrigT1CaloDefs.h"
#include "TrigT1CaloCondSvc/L1CaloCondSvc.h"
#include "TrigT1CaloCalibConditions/L1CaloCoolChannelId.h"
#include "TrigT1CaloCalibConditions/L1CaloPprLutContainer.h"
#include "TrigT1CaloCalibConditions/L1CaloPprLut.h"
#include "TrigT1CaloMonitoring/TrigT1CaloLWHistogramTool.h"

#include "TrigT1CaloMonitoring/PPMSimBSMon.h"

/*---------------------------------------------------------*/
PPMSimBSMon::PPMSimBSMon(const std::string & type, 
			 const std::string & name,
			 const IInterface* parent)
  : ManagedMonitorToolBase(type, name, parent),
    m_l1CondSvc("L1CaloCondSvc", name),
    m_ttTool("LVL1::L1TriggerTowerTool/L1TriggerTowerTool"), 
    m_histTool("TrigT1CaloLWHistogramTool"),
    m_LutContainer(0), m_debug(false), m_events(0),
    m_histBooked(false),
    m_h_ppm_em_2d_etaPhi_tt_lut_SimEqData(0),
    m_h_ppm_em_2d_etaPhi_tt_lut_SimNeData(0),
    m_h_ppm_em_2d_etaPhi_tt_lut_SimNoData(0),
    m_h_ppm_em_2d_etaPhi_tt_lut_DataNoSim(0),
    m_h_ppm_had_2d_etaPhi_tt_lut_SimEqData(0),
    m_h_ppm_had_2d_etaPhi_tt_lut_SimNeData(0),
    m_h_ppm_had_2d_etaPhi_tt_lut_SimNoData(0),
    m_h_ppm_had_2d_etaPhi_tt_lut_DataNoSim(0),
    m_h_ppm_em_2d_etaPhi_tt_ped_runavg(0),
    m_h_ppm_had_2d_etaPhi_tt_ped_runavg(0),
    m_h_ppm_em_2d_etaPhi_tt_ped_worstavg(0),
    m_h_ppm_had_2d_etaPhi_tt_ped_worstavg(0),
    m_h_ppm_em_2d_etaPhi_tt_ped_runrms(0),
    m_h_ppm_had_2d_etaPhi_tt_ped_runrms(0),
    m_h_ppm_em_2d_etaPhi_tt_ped_instavg(0),
    m_h_ppm_had_2d_etaPhi_tt_ped_instavg(0),
    m_h_ppm_em_2d_etaPhi_tt_ped_instrms(0),
    m_h_ppm_had_2d_etaPhi_tt_ped_instrms(0),
    m_h_ppm_em_2d_etaPhi_tt_ped_instavg_B(0),
    m_h_ppm_had_2d_etaPhi_tt_ped_instavg_B(0),
    m_h_ppm_em_2d_etaPhi_tt_ped_instrms_B(0),
    m_h_ppm_had_2d_etaPhi_tt_ped_instrms_B(0),
    m_h_ppm_2d_LUT_MismatchEvents_cr0cr1(0),
    m_h_ppm_2d_LUT_MismatchEvents_cr2cr3(0),
    m_h_ppm_2d_LUT_MismatchEvents_cr4cr5(0),
    m_h_ppm_2d_LUT_MismatchEvents_cr6cr7(0)
/*---------------------------------------------------------*/
{
  declareProperty("TriggerTowerLocation",
                 m_triggerTowerLocation =
		  LVL1::TrigT1CaloDefs::TriggerTowerLocation);
  
  declareProperty("RootDirectory", m_rootDir = "L1Calo");

  declareProperty("PedestalSampleSize", m_instantaneous=400);
  declareProperty("OnlineTest", m_onlineTest = false,
                  "Test online code when running offline");
}

/*---------------------------------------------------------*/
PPMSimBSMon::~PPMSimBSMon()
/*---------------------------------------------------------*/
{
}

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

/*---------------------------------------------------------*/
StatusCode PPMSimBSMon:: initialize()
/*---------------------------------------------------------*/
{
  msg(MSG::INFO) << "Initializing " << name() << " - package version "
                 << PACKAGE_VERSION << endreq;
  m_debug = msgLvl(MSG::DEBUG);

  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;

  sc = m_l1CondSvc.retrieve();
  if(sc.isFailure()){
    msg(MSG::WARNING) << "Could not retrieve L1CaloCondSvc" << endreq;
  }

  sc = m_ttTool.retrieve();
  if( sc.isFailure() ) {
    msg(MSG::ERROR) << "Unable to locate Tool L1TriggerTowerTool" << endreq;
    return sc;
  }

  sc = m_histTool.retrieve();
  if( sc.isFailure() ) {
    msg(MSG::ERROR) << "Unable to locate Tool TrigT1CaloLWHistogramTool"
                    << endreq;
    return sc;
  }

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
  msg(MSG::DEBUG) << "bookHistograms entered" << endreq;

  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
  	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  if ( isNewEventsBlock || isNewLumiBlock ) { }
  
  if( isNewRun ) {
    
    
    m_LutContainer = 0;
    if (m_l1CondSvc) {
      msg(MSG::DEBUG) << "Retrieving Conditions Containers " << endreq;
      
      StatusCode sc = m_l1CondSvc->retrieve(m_LutContainer);
      if (sc.isFailure()) {
	msg(MSG::WARNING) << "Could not retrieve LutContainer " << endreq;
	return sc;
      }
      msg(MSG::DEBUG) << "Retrieved LutContainer " << endreq;
    }
    
    else {
      msg(MSG::WARNING) << "Could not retrieve Conditions Containers" << endreq;
    }
    
    

  std::string dir(m_rootDir + "/PPM/Errors/Data_Simulation");
  MonGroup monShift ( this, dir, shift, run );
  MonGroup monExpert( this, dir, expert, run );
  MonGroup monPPM   ( this, dir + "/PPMLUTSim", expert, run );
  MonGroup monEvent ( this, dir + "/MismatchEventNumbers", expert, run, "", "eventSample" );
  std::string dirPed (m_rootDir + "/PPM/ADC/Pedestal");
  MonGroup monPed ( this, dirPed, expert, run);
  MonGroup monPedrms ( this, dirPed, expert, run,"","mergeRMS");

  // LUT

  m_histTool->setMonGroup(&monPPM);

  m_h_ppm_em_2d_etaPhi_tt_lut_SimEqData = m_histTool->bookPPMEmEtaVsPhi(
    "ppm_em_2d_etaPhi_tt_lut_SimEqData",
    "PPM LUT EM Data/Simulation Non-zero Matches");
  m_h_ppm_em_2d_etaPhi_tt_lut_SimNeData = m_histTool->bookPPMEmEtaVsPhi(
    "ppm_em_2d_etaPhi_tt_lut_SimNeData",
    "PPM LUT EM Data/Simulation Non-zero Mismatches");
  m_h_ppm_em_2d_etaPhi_tt_lut_SimNoData = m_histTool->bookPPMEmEtaVsPhi(
    "ppm_em_2d_etaPhi_tt_lut_SimNoData",
    "PPM LUT EM Simulation but no Data");
  m_h_ppm_em_2d_etaPhi_tt_lut_DataNoSim = m_histTool->bookPPMEmEtaVsPhi(
    "ppm_em_2d_etaPhi_tt_lut_DataNoSim",
    "PPM LUT EM Data but no Simulation");
  m_h_ppm_had_2d_etaPhi_tt_lut_SimEqData = m_histTool->bookPPMHadEtaVsPhi(
    "ppm_had_2d_etaPhi_tt_lut_SimEqData",
    "PPM LUT HAD Data/Simulation Non-zero Matches");
  m_h_ppm_had_2d_etaPhi_tt_lut_SimNeData = m_histTool->bookPPMHadEtaVsPhi(
    "ppm_had_2d_etaPhi_tt_lut_SimNeData",
    "PPM LUT HAD Data/Simulation Non-zero Mismatches");
  m_h_ppm_had_2d_etaPhi_tt_lut_SimNoData = m_histTool->bookPPMHadEtaVsPhi(
    "ppm_had_2d_etaPhi_tt_lut_SimNoData",
    "PPM LUT HAD Simulation but no Data");
  m_h_ppm_had_2d_etaPhi_tt_lut_DataNoSim = m_histTool->bookPPMHadEtaVsPhi(
    "ppm_had_2d_etaPhi_tt_lut_DataNoSim",
    "PPM LUT HAD Data but no Simulation");

  // Overall pedestal
 
  m_histTool->setMonGroup(&monPed);

  m_h_ppm_em_2d_etaPhi_tt_ped_runavg = m_histTool->bookProfilePPMEmEtaVsPhi(
    "ppm_em_2d_etaPhi_tt_ped_runavg",
    "PPM Mean Pedestal Difference EM (over run)");
  m_h_ppm_had_2d_etaPhi_tt_ped_runavg = m_histTool->bookProfilePPMHadEtaVsPhi(
    "ppm_had_2d_etaPhi_tt_ped_runavg",
    "PPM Mean Pedestal Difference Had (over run)");
  
  m_histTool->setMonGroup(&monPedrms);

  m_h_ppm_em_2d_etaPhi_tt_ped_runrms = m_histTool->bookPPMEmEtaVsPhi(
    "ppm_em_2d_etaPhi_tt_ped_runrms",
    "PPM rms Pedestal Difference EM (over run)");
  m_h_ppm_had_2d_etaPhi_tt_ped_runrms = m_histTool->bookPPMHadEtaVsPhi(
    "ppm_had_2d_etaPhi_tt_ped_runrms",
    "PPM rms Pedestal Difference Had (over run)");

  if(m_onlineTest || m_environment == AthenaMonManager::online) {

    m_histTool->setMonGroup(&monPed);

    m_h_ppm_em_2d_etaPhi_tt_ped_worstavg = m_histTool->bookPPMEmEtaVsPhi(
      "ppm_em_2d_etaPhi_tt_ped_worstavg",
      "PPM Worst Pedestal Difference EM (over run)");
    m_h_ppm_had_2d_etaPhi_tt_ped_worstavg = m_histTool->bookPPMHadEtaVsPhi(
      "ppm_had_2d_etaPhi_tt_ped_worstavg",
      "PPM Worst Pedestal Difference Had (over run)");

    m_h_ppm_em_2d_etaPhi_tt_ped_instavg = m_histTool->bookProfilePPMEmEtaVsPhi(
      "ppm_em_2d_etaPhi_tt_ped_instavg",
      "PPM Mean Pedestal Difference EM (instantaneous)");
    m_h_ppm_had_2d_etaPhi_tt_ped_instavg = m_histTool->bookProfilePPMHadEtaVsPhi(
      "ppm_had_2d_etaPhi_tt_ped_instavg",
      "PPM Mean Pedestal Difference Had (instantaneous)");
    
    m_histTool->setMonGroup(&monPedrms);

    m_h_ppm_em_2d_etaPhi_tt_ped_instrms = m_histTool->bookPPMEmEtaVsPhi(
      "ppm_em_2d_etaPhi_tt_ped_instrms",
      "PPM rms Pedestal Difference EM (instantaneous)");
    m_h_ppm_had_2d_etaPhi_tt_ped_instrms = m_histTool->bookPPMHadEtaVsPhi(
      "ppm_had_2d_etaPhi_tt_ped_instrms",
      "PPM rms Pedestal Difference Had (instantaneous)");
    
    m_histTool->unsetMonGroup();

    m_h_ppm_em_2d_etaPhi_tt_ped_instavg_B = m_histTool->bookProfilePPMEmEtaVsPhi(
      "ppm_em_2d_etaPhi_tt_ped_instavg_B",
      "PPM Mean Pedestal Difference EM (instantaneous [B])");
    m_h_ppm_had_2d_etaPhi_tt_ped_instavg_B = m_histTool->bookProfilePPMHadEtaVsPhi(
      "ppm_had_2d_etaPhi_tt_ped_instavg_B",
      "PPM Mean Pedestal Difference Had (instantaneous [B])");
   
    m_h_ppm_em_2d_etaPhi_tt_ped_instrms_B = m_histTool->bookPPMEmEtaVsPhi(
      "ppm_em_2d_etaPhi_tt_ped_instrms_B",
      "PPM rms Pedestal Difference EM (instantaneous [B])");
    m_h_ppm_had_2d_etaPhi_tt_ped_instrms_B = m_histTool->bookPPMHadEtaVsPhi(
      "ppm_had_2d_etaPhi_tt_ped_instrms_B",
      "PPM rms Pedestal Difference Had (instantaneous [B])");
  
  }
							
  // Mismatch Histograms

  m_histTool->setMonGroup(&monEvent);

  m_h_ppm_2d_LUT_MismatchEvents_cr0cr1 = m_histTool->bookPPMEventVsCrateModule(
    "ppm_2d_LUT_MismatchEvents_cr0cr1","PPM LUT Mismatch Event Numbers",0,1);
  m_h_ppm_2d_LUT_MismatchEvents_cr2cr3 = m_histTool->bookPPMEventVsCrateModule(
    "ppm_2d_LUT_MismatchEvents_cr2cr3","PPM LUT Mismatch Event Numbers",2,3);
  m_h_ppm_2d_LUT_MismatchEvents_cr4cr5 = m_histTool->bookPPMEventVsCrateModule(
    "ppm_2d_LUT_MismatchEvents_cr4cr5","PPM LUT Mismatch Event Numbers",4,5);
  m_h_ppm_2d_LUT_MismatchEvents_cr6cr7 = m_histTool->bookPPMEventVsCrateModule(
    "ppm_2d_LUT_MismatchEvents_cr6cr7","PPM LUT Mismatch Event Numbers",6,7);

  m_histTool->unsetMonGroup();
  m_histBooked = true;

  } // end if (isNewRun ...

  msg(MSG::DEBUG) << "Leaving bookHistograms" << endreq;
  
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode PPMSimBSMon::fillHistograms()
/*---------------------------------------------------------*/
{
  if (m_debug) msg(MSG::DEBUG) << "fillHistograms entered" << endreq;

  if (!m_histBooked) {
    if (debug) msg(MSG::DEBUG) << "Histogram(s) not booked" << endreq;
    return StatusCode::SUCCESS;
  }

  StatusCode sc;

  //Retrieve Trigger Towers from SG
  const TriggerTowerCollection* triggerTowerTES = 0; 
  sc = evtStore()->retrieve(triggerTowerTES, m_triggerTowerLocation); 
  if( sc.isFailure()  ||  !triggerTowerTES ) {
    if (m_debug) msg(MSG::DEBUG) << "No Trigger Tower container found"<< endreq; 
  }
  ++m_events;

  // Compare LUT simulated from FADC with LUT from data
 
  if (triggerTowerTES) {
    simulateAndCompare(triggerTowerTES);
  }


  if ((m_environment == AthenaMonManager::online || m_onlineTest) && (m_events%int(0.5*m_instantaneous)==0)) {
    double worst_value = 0.;
    double minValue    = 0.;
    double maxValue    = 0.;
    int worst_binx = 0;
    int worst_biny = 0;
    int minBinx    = 0;
    int minBiny    = 0;
    int maxBinx    = 0;
    int maxBiny    = 0;
    m_histTool->getMinMaxBin(m_h_ppm_em_2d_etaPhi_tt_ped_instavg,
                       minBinx, minBiny, maxBinx, maxBiny, minValue, maxValue);
    if(fabs(minValue) > fabs(maxValue)){
      worst_value = minValue;
      worst_binx  = minBinx;
      worst_biny  = minBiny;
    }
    else {
      worst_value = maxValue;
      worst_binx  = maxBinx;
      worst_biny  = maxBiny;
    }
    if(fabs(worst_value)>fabs(m_h_ppm_em_2d_etaPhi_tt_ped_worstavg->GetBinContent(worst_binx, worst_biny))) {
      m_histTool->setBinPPMEmEtaVsPhi(m_h_ppm_em_2d_etaPhi_tt_ped_worstavg,
                                      worst_binx, worst_biny, worst_value);
    }
    m_histTool->replaceContents(m_h_ppm_em_2d_etaPhi_tt_ped_instavg,
                                m_h_ppm_em_2d_etaPhi_tt_ped_instavg_B);
    m_h_ppm_em_2d_etaPhi_tt_ped_instavg_B->Reset();
    m_histTool->replaceContents(m_h_ppm_em_2d_etaPhi_tt_ped_instrms,
                                m_h_ppm_em_2d_etaPhi_tt_ped_instrms_B);
    m_h_ppm_em_2d_etaPhi_tt_ped_instrms_B->Reset();	    
    
    m_histTool->getMinMaxBin(m_h_ppm_had_2d_etaPhi_tt_ped_instavg,
                       minBinx, minBiny, maxBinx, maxBiny, minValue, maxValue);
    if(fabs(minValue) > fabs(maxValue)){
      worst_value = minValue;
      worst_binx  = minBinx;
      worst_biny  = minBiny;
    }
    else {
      worst_value = maxValue;
      worst_binx  = maxBinx;
      worst_biny  = maxBiny;
    }
    if(fabs(worst_value)>fabs(m_h_ppm_had_2d_etaPhi_tt_ped_worstavg->GetBinContent(worst_binx, worst_biny))) {
      m_histTool->setBinPPMHadEtaVsPhi(m_h_ppm_had_2d_etaPhi_tt_ped_worstavg,
                                       worst_binx, worst_biny, worst_value);
    }
    m_histTool->replaceContents(m_h_ppm_had_2d_etaPhi_tt_ped_instavg,
                                m_h_ppm_had_2d_etaPhi_tt_ped_instavg_B);
    m_h_ppm_had_2d_etaPhi_tt_ped_instavg_B->Reset();
    m_histTool->replaceContents(m_h_ppm_had_2d_etaPhi_tt_ped_instrms,
                                m_h_ppm_had_2d_etaPhi_tt_ped_instrms_B);
    m_h_ppm_had_2d_etaPhi_tt_ped_instrms_B->Reset();
  }
  
  if (m_debug) msg(MSG::DEBUG) << "Leaving fillHistograms" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode PPMSimBSMon::procHistograms(bool isEndOfEventsBlock,
                                  bool isEndOfLumiBlock, bool isEndOfRun)
/*---------------------------------------------------------*/
{
  msg(MSG::DEBUG) << "procHistograms entered" << endreq;

  if (isEndOfEventsBlock || isEndOfLumiBlock) {
  }

  if (isEndOfRun) {
    if(m_onlineTest || m_environment == AthenaMonManager::online) {
      LWHist::safeDelete(m_h_ppm_em_2d_etaPhi_tt_ped_instavg_B);
      LWHist::safeDelete(m_h_ppm_had_2d_etaPhi_tt_ped_instavg_B);
      LWHist::safeDelete(m_h_ppm_em_2d_etaPhi_tt_ped_instrms_B);
      LWHist::safeDelete(m_h_ppm_had_2d_etaPhi_tt_ped_instrms_B);
    }
  }

  return StatusCode::SUCCESS;
}

void PPMSimBSMon::simulateAndCompare(const TriggerTowerCollection* ttIn)
{
  if (m_debug) msg(MSG::DEBUG) << "Simulate LUT data from FADC data" << endreq;

  StatusCode sc = m_ttTool->retrieveConditions();
  if (sc.isFailure()) return;

  const int nCrates = 8;
  ErrorVector crateError(nCrates);
  ErrorVector moduleError(nCrates);
  
  m_ttTool->setDebug(false);
  TriggerTowerCollection::const_iterator iter  = ttIn->begin();
  TriggerTowerCollection::const_iterator iterE = ttIn->end();

  for (; iter != iterE; ++iter) {
    
    const LVL1::TriggerTower* tt = *iter;
    
    const L1CaloCoolChannelId em_coolId(m_ttTool->channelID(tt->eta(),
                                                            tt->phi(),0));
    const L1CaloCoolChannelId had_coolId(m_ttTool->channelID(tt->eta(),
                                                             tt->phi(),1));

    const int had_crate  = had_coolId.crate();
    const int had_module = had_coolId.module();
    const int em_crate   = em_coolId.crate();
    const int em_module  = em_coolId.module();
 
    std::vector<int> emLut;
    std::vector<int> emBcidR;
    std::vector<int> emBcidD;
    m_ttTool->process(tt->emADC(),em_coolId, emLut, emBcidR, emBcidD);
    const int emPeak = tt->emADCPeak();
    std::vector<int> emLut1;
    const int emSlices = (tt->emADC()).size();
    if (emSlices < 7 || emBcidD[emPeak]) emLut1.push_back(emLut[emPeak]);
    else                 emLut1.push_back(0);
    std::vector<int> emBcidR1;
    emBcidR1.push_back(emBcidR[emPeak]);
    if (m_debug && emLut1[0] != tt->emEnergy() && (emSlices>=7 || tt->emEnergy()!=0)) { // mismatch - repeat with debug on
      std::vector<int> emLut2; 
      std::vector<int> emBcidR2;
      std::vector<int> emBcidD2;
      m_ttTool->setDebug(true);
      m_ttTool->process(tt->emADC(),em_coolId, emLut2, emBcidR2, emBcidD2);
      m_ttTool->setDebug(false);
    }
    
    std::vector<int> hadLut;
    std::vector<int> hadBcidR;
    std::vector<int> hadBcidD;
    m_ttTool->process(tt->hadADC(),had_coolId, hadLut, hadBcidR, hadBcidD);
    const int hadPeak = tt->hadADCPeak();
    std::vector<int> hadLut1;
    const int hadSlices = (tt->hadADC()).size();
    if (hadSlices < 7 || hadBcidD[hadPeak]) hadLut1.push_back(hadLut[hadPeak]);
    else                   hadLut1.push_back(0);
    std::vector<int> hadBcidR1;
    hadBcidR1.push_back(hadBcidR[hadPeak]);
    if (m_debug && hadLut1[0] != tt->hadEnergy() && (hadSlices>=7 || tt->hadEnergy()!=0)) {
      std::vector<int> hadLut2;
      std::vector<int> hadBcidR2;
      std::vector<int> hadBcidD2;
      m_ttTool->setDebug(true);
      m_ttTool->process(tt->hadADC(),had_coolId, hadLut2, hadBcidR2, hadBcidD2);
      m_ttTool->setDebug(false);
    }
    
    const int simEm  = emLut1[0];
    const int simHad = hadLut1[0];
    const int datEm  = tt->emEnergy();
    const int datHad = tt->hadEnergy();

    int em_offset = 0;
    
    if (m_LutContainer) {
      const L1CaloPprLut* em_LUT = m_LutContainer->pprLut(em_coolId);
      if (em_LUT) {
	em_offset   = em_LUT->pedMean();
      } 
      else if (m_debug) msg(MSG::DEBUG) << "::lut: No L1CaloPprLut found" << endreq;
    } 
    else if (m_debug) msg(MSG::DEBUG) << "::lut: No LUT Container retrieved" << endreq;

    const bool em_notSignal = (datEm==0 && !(m_ttTool->disabledChannel(em_coolId)));

    int had_offset = 0;

    if (m_LutContainer) {
      const L1CaloPprLut* had_LUT = m_LutContainer->pprLut(had_coolId);
      if (had_LUT) {
	had_offset   = had_LUT->pedMean();
      } 
      else if (m_debug) msg(MSG::DEBUG) << "::lut: No L1CaloPprLut found" << endreq;
    } 
    else if (m_debug) msg(MSG::DEBUG) << "::lut: No LUT Container retrieved" << endreq;

    const bool had_notSignal = (datHad==0 && !(m_ttTool->disabledChannel(had_coolId)));
    
    const double eta = tt->eta();
    const double phi = tt->phi();
  
    int tt_accept = 1;

    for(int i=0; i<emSlices; i++) {
      if(em_notSignal) {
	if(fabs((tt->emADC()).at(i)-em_offset)>10) tt_accept = 0;
      }
    }

    if(tt_accept==1) {
      if(em_notSignal) {
	
	for(int i=0; i<emSlices; i++) {
	  
	  m_histTool->fillPPMEmEtaVsPhi(m_h_ppm_em_2d_etaPhi_tt_ped_runavg,
	                          eta, phi, ((tt->emADC()).at(i)-em_offset));

	  if (m_environment == AthenaMonManager::online || m_onlineTest) {
	    m_histTool->fillPPMEmEtaVsPhi(m_h_ppm_em_2d_etaPhi_tt_ped_instavg,
	                          eta, phi, ((tt->emADC()).at(i)-em_offset));
	    m_histTool->fillPPMEmEtaVsPhi(m_h_ppm_em_2d_etaPhi_tt_ped_instavg_B,
	                          eta, phi, ((tt->emADC()).at(i)-em_offset));
	  }
	  
          int binx = 0;
          int biny = 0;
	  m_histTool->findBinPPMEmEtaVsPhi(m_h_ppm_em_2d_etaPhi_tt_ped_runavg,
	                                                 eta, phi, binx, biny);
	  double entries = 0.;
	  double val = 0.;
	  double rms = 0.;
	  m_h_ppm_em_2d_etaPhi_tt_ped_runavg->GetBinInfo(binx, biny, entries,
	                                                             val, rms);
	  m_histTool->setBinPPMEmEtaVsPhi(m_h_ppm_em_2d_etaPhi_tt_ped_runrms,
	                                binx, biny, rms, rms/sqrt(2.*entries));
	  if (m_environment == AthenaMonManager::online || m_onlineTest) {
	    m_h_ppm_em_2d_etaPhi_tt_ped_instavg->GetBinInfo(binx, biny, entries,
	                                                             val, rms);
	    m_histTool->setBinPPMEmEtaVsPhi(m_h_ppm_em_2d_etaPhi_tt_ped_instrms,
	                                binx, biny, rms, rms/sqrt(2.*entries));
	    m_h_ppm_em_2d_etaPhi_tt_ped_instavg_B->GetBinInfo(binx, biny,
	                                                    entries, val, rms);
	    m_histTool->setBinPPMEmEtaVsPhi(m_h_ppm_em_2d_etaPhi_tt_ped_instrms_B,
	                                binx, biny, rms, rms/sqrt(2.*entries));
	    
	  }
	}	 
      }	
    }
      
    tt_accept = 1;

    for(int i=0;i<hadSlices; i++) {
      if(had_notSignal) {
	if(fabs((tt->hadADC()).at(i)-had_offset)>10) tt_accept = 0;
      }
    }

    if(tt_accept==1) {
      if(had_notSignal) {

	for(int i=0; i<hadSlices; i++) {
	  
	  m_histTool->fillPPMHadEtaVsPhi(m_h_ppm_had_2d_etaPhi_tt_ped_runavg,
	                          eta, phi, ((tt->hadADC()).at(i)-had_offset));

	  if (m_environment == AthenaMonManager::online || m_onlineTest) {
	    m_histTool->fillPPMHadEtaVsPhi(m_h_ppm_had_2d_etaPhi_tt_ped_instavg,
	                          eta, phi, ((tt->hadADC()).at(i)-had_offset));
	    m_histTool->fillPPMHadEtaVsPhi(m_h_ppm_had_2d_etaPhi_tt_ped_instavg_B,
	                          eta, phi, ((tt->hadADC()).at(i)-had_offset));
	  }
          int binx = 0;
          int biny = 0;
	  m_histTool->findBinPPMHadEtaVsPhi(m_h_ppm_had_2d_etaPhi_tt_ped_runavg,
	                                                 eta, phi, binx, biny);
	  double entries = 0.;
	  double val = 0.;
	  double rms = 0.;
	  m_h_ppm_had_2d_etaPhi_tt_ped_runavg->GetBinInfo(binx, biny, entries,
	                                                             val, rms);
	  m_histTool->setBinPPMHadEtaVsPhi(m_h_ppm_had_2d_etaPhi_tt_ped_runrms,
	                                binx, biny, rms, rms/sqrt(2.*entries));
	  if (m_environment == AthenaMonManager::online || m_onlineTest) {
	    m_h_ppm_had_2d_etaPhi_tt_ped_instavg->GetBinInfo(binx, biny,
	                                                    entries, val, rms);
	    m_histTool->setBinPPMHadEtaVsPhi(m_h_ppm_had_2d_etaPhi_tt_ped_instrms,
	                                binx, biny, rms, rms/sqrt(2.*entries));
	    m_h_ppm_had_2d_etaPhi_tt_ped_instavg_B->GetBinInfo(binx, biny,
	                                                    entries, val, rms);
	    m_histTool->setBinPPMHadEtaVsPhi(m_h_ppm_had_2d_etaPhi_tt_ped_instrms_B,
	                                binx, biny, rms, rms/sqrt(2.*entries));
	  
	  }
	}
      }
    }

    if (!simEm && !simHad && !datEm && !datHad) continue;
    
    //  Fill in error plots
    
    int em_mismatch = 0;
    int had_mismatch = 0;
    
    TH2F_LW* hist1 = 0;
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
        msg(MSG::DEBUG) << " EMTowerMismatch eta/phi/sim/dat: "
              << eta << "/" << phi << "/" << simEm << "/" << datEm << endreq;
      }
    }
    
    if (hist1) m_histTool->fillPPMEmEtaVsPhi(hist1, eta, phi);
    
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
        msg(MSG::DEBUG) << " HadTowerMismatch eta/phi/sim/dat: "
              << eta << "/" << phi << "/" << simHad << "/" << datHad << endreq;
      }
    }

    if (hist1) m_histTool->fillPPMHadEtaVsPhi(hist1, eta, phi);
      
    if(had_mismatch==1) {
      crateError[had_crate] = 1;
      if (!((moduleError[had_crate]>>had_module)&0x1)) {
	fillEventSample(had_crate,had_module);
	moduleError[had_crate] |= (1 << had_module);
      }
    }
  
  }    
    
  ErrorVector* save = new ErrorVector(crateError);
  sc = evtStore()->record(save, "L1CaloPPMMismatchVector");
  if (sc.isFailure()) {
    msg(MSG::ERROR) << "Error recording PPM mismatch vector in TES "
	            << endreq;
  }
  
  m_ttTool->setDebug(true);
  
}

void PPMSimBSMon::fillEventSample(int crate, int module)
{
  const int y = module + 16*(crate%2);
  TH2I_LW* hist = 0;
  if     (crate==0 || crate==1) hist = m_h_ppm_2d_LUT_MismatchEvents_cr0cr1;
  else if(crate==2 || crate==3) hist = m_h_ppm_2d_LUT_MismatchEvents_cr2cr3;
  else if(crate==4 || crate==5) hist = m_h_ppm_2d_LUT_MismatchEvents_cr4cr5;
  else if(crate==6 || crate==7) hist = m_h_ppm_2d_LUT_MismatchEvents_cr6cr7;
  if (hist) m_histTool->fillEventNumber(hist, y);
}


