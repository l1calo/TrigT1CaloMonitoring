// ********************************************************************
//
// NAME:     CPMSimBSMon.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Peter Faulkner
//           
//
// ********************************************************************

#include <sstream>
#include <utility>
#include <cmath>

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "StoreGate/StoreGateSvc.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigT1Calo/CPAlgorithm.h"
#include "TrigT1Calo/CMMCPHits.h"
#include "TrigT1Calo/CPMHits.h"
#include "TrigT1Calo/CPMTower.h"
#include "TrigT1Calo/CPMRoI.h"
#include "TrigT1Calo/CoordToHardware.h"
#include "TrigT1Calo/TriggerTower.h"
#include "TrigT1Calo/TriggerTowerKey.h"
#include "TrigT1CaloTools/IL1EmTauTools.h"
#include "TrigT1CaloTools/IL1CPHitsTools.h"
#include "TrigT1Interfaces/Coordinate.h"
#include "TrigT1Interfaces/CoordinateRange.h"
#include "TrigT1Interfaces/CPRoIDecoder.h"
#include "TrigT1Interfaces/TrigT1CaloDefs.h"

#include "TrigT1CaloMonitoring/CPMSimBSMon.h"


/*---------------------------------------------------------*/
CPMSimBSMon::CPMSimBSMon(const std::string & type, 
			 const std::string & name,
			 const IInterface* parent)
  : ManagedMonitorToolBase(type, name, parent),
    m_storeGate("StoreGateSvc", name),
    m_emTauTool("LVL1::L1EmTauTools/L1EmTauTools"),
    m_cpHitsTool("LVL1::L1CPHitsTools/L1CPHitsTools"),
    m_monGroup(0), m_events(0)
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
  declareProperty("CPMRoILocationRoIB",
                 m_cpmRoiLocationRoib =
		                 LVL1::TrigT1CaloDefs::CPMRoILocation+"RoIB");
  declareProperty("TriggerTowerLocation",
                 m_triggerTowerLocation =
		                 LVL1::TrigT1CaloDefs::TriggerTowerLocation);

  declareProperty("RootDirectory", m_rootDir = "L1Calo");
  declareProperty("PhiUnits", m_phiUnits = "channels",
                  "Phi Units: radians, degrees or channels");
  declareProperty("CompareWithSimulation", m_compareWithSim = true,
                  "Include the checks that need to rerun simulation");
}

/*---------------------------------------------------------*/
CPMSimBSMon::~CPMSimBSMon()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode CPMSimBSMon:: initialize()
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

  if (m_compareWithSim) {

    sc = m_emTauTool.retrieve();
    if( sc.isFailure() ) {
      log << MSG::ERROR << "Unable to locate Tool L1EmTauTools" << endreq;
      return sc;
    }

    sc = m_cpHitsTool.retrieve();
    if( sc.isFailure() ) {
      log << MSG::ERROR << "Unable to locate Tool L1CPHitsTools" << endreq;
      return sc;
    }
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
StatusCode CPMSimBSMon::bookHistograms(bool isNewEventsBlock,
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

  if ( isNewEventsBlock || isNewLumiBlock ) { }

  if( isNewRun ) {

  std::string dir(m_rootDir + "/5_CP_Sim_Check");
  MonGroup monShift ( this, dir, shift, run );
  MonGroup monExpert( this, dir, expert, run );
  MonGroup monRoIs  ( this, dir + "/RoIs", expert, run );
  MonGroup monCPM   ( this, dir + "/CPM",  expert, run );
  MonGroup monCMM   ( this, dir + "/CMM",  expert, run );

  // CPMTowers

  m_monGroup = &monCPM;

  m_h_EMTowerSIMeqDAT = book2F("EMTowerSIMeqDAT",
    "CPM Towers EM Data/Simulation Non-zero Matches;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_EMTowerSIMneDAT = book2F("EMTowerSIMneDAT",
    "CPM Towers EM Data/Simulation Non-zero Mismatches;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_EMTowerSIMnoDAT = book2F("EMTowerSIMnoDAT",
    "CPM Towers EM Simulation but no Data;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_EMTowerDATnoSIM = book2F("EMTowerDATnoSIM",
    "CPM Towers EM Data but no Simulation;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_HadTowerSIMeqDAT = book2F("HadTowerSIMeqDAT",
    "CPM Towers HAD Data/Simulation Non-zero Matches;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_HadTowerSIMneDAT = book2F("HadTowerSIMneDAT",
    "CPM Towers HAD Data/Simulation Non-zero Mismatches;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_HadTowerSIMnoDAT = book2F("HadTowerSIMnoDAT",
    "CPM Towers HAD Simulation but no Data;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_HadTowerDATnoSIM = book2F("HadTowerDATnoSIM",
    "CPM Towers HAD Data but no Simulation;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);

  if (m_compareWithSim) {

  // RoIs

  m_monGroup = &monRoIs;

  m_h_RoISIMeqDAT = book2F("RoISIMeqDAT",
     "CPM RoI Data/Simulation Non-zero Matches;Crate/Module;Chip/Local Coord",
             56, 0, 56, 64, 0, 64);
  setLabelsCMCC(m_h_RoISIMeqDAT);
  m_h_RoISIMneDAT = book2F("RoISIMneDAT",
   "CPM RoI Data/Simulation Non-zero Mismatches;Crate/Module;Chip/Local Coord",
             56, 0, 56, 64, 0, 64);
  setLabelsCMCC(m_h_RoISIMneDAT);
  m_h_RoISIMnoDAT = book2F("RoISIMnoDAT",
            "CPM RoI Simulation but no Data;Crate/Module;Chip/Local Coord",
             56, 0, 56, 64, 0, 64);
  setLabelsCMCC(m_h_RoISIMnoDAT);
  m_h_RoIDATnoSIM = book2F("RoIDATnoSIM",
            "CPM RoI Data but no Simulation;Crate/Module;Chip/Local Coord",
             56, 0, 56, 64, 0, 64);
  setLabelsCMCC(m_h_RoIDATnoSIM);
  m_h_RoIThreshSIMeqDAT = book2F("RoIThreshSIMeqDAT",
     "CPM RoI Data/Simulation Threshold Matches;Crate/Module;Threshold",
             56, 0, 56, 16, 0, 16);
  setLabelsCMT(m_h_RoIThreshSIMeqDAT);
  m_h_RoIThreshSIMneDAT = book2F("RoIThreshSIMneDAT",
     "CPM RoI Data/Simulation Threshold Mismatches;Crate/Module;Threshold",
             56, 0, 56, 16, 0, 16);
  setLabelsCMT(m_h_RoIThreshSIMneDAT);
  const double halfPhiBin = m_phiMax/128.;
  m_h_RoIEtaPhiSIMeqDAT = book2F("RoIEtaPhiSIMeqDAT",
     "CPM RoI Data/Simulation Non-zero Matches;eta;phi",
             51, -2.55, 2.55, 65, -halfPhiBin, m_phiMax+halfPhiBin);
  m_h_RoIEtaPhiSIMneDAT = book2F("RoIEtaPhiSIMneDAT",
     "CPM RoI Data/Simulation Non-zero Mismatches;eta;phi",
             51, -2.55, 2.55, 65, -halfPhiBin, m_phiMax+halfPhiBin);
  m_h_RoIEtaPhiSIMnoDAT = book2F("RoIEtaPhiSIMnoDAT",
            "CPM RoI Simulation but no Data;eta;phi",
             51, -2.55, 2.55, 65, -halfPhiBin, m_phiMax+halfPhiBin);
  m_h_RoIEtaPhiDATnoSIM = book2F("RoIEtaPhiDATnoSIM",
            "CPM RoI Data but no Simulation;eta;phi",
             51, -2.55, 2.55, 65, -halfPhiBin, m_phiMax+halfPhiBin);

  // CPMHits

  m_monGroup = &monCPM;

  m_h_CPMHitsSIMeqDAT = book2F("CPMHitsSIMeqDAT",
     "CPM Hits Data/Simulation Non-zero Matches;Module;Crate",
             14, 1, 15, 4, 0, 4);
  setLabelsMC(m_h_CPMHitsSIMeqDAT);
  m_h_CPMHitsSIMneDAT = book2F("CPMHitsSIMneDAT",
     "CPM Hits Data/Simulation Non-zero Mismatches;Module;Crate",
             14, 1, 15, 4, 0, 4);
  setLabelsMC(m_h_CPMHitsSIMneDAT);
  m_h_CPMHitsSIMnoDAT = book2F("CPMHitsSIMnoDAT",
     "CPM Hits Simulation but no Data;Module;Crate", 14, 1, 15, 4, 0, 4);
  setLabelsMC(m_h_CPMHitsSIMnoDAT);
  m_h_CPMHitsDATnoSIM = book2F("CPMHitsDATnoSIM",
     "CPM Hits Data but no Simulation;Module;Crate", 14, 1, 15, 4, 0, 4);
  setLabelsMC(m_h_CPMHitsDATnoSIM);
  m_h_CPMHitsThreshSIMeqDAT = book2F("CPMHitsThreshSIMeqDAT",
     "CPM Hits Data/Simulation Threshold Matches;Crate/Module;Threshold",
             56, 0, 56, 16, 0, 16);
  setLabelsCMT(m_h_CPMHitsThreshSIMeqDAT);
  m_h_CPMHitsThreshSIMneDAT = book2F("CPMHitsThreshSIMneDAT",
     "CPM Hits Data/Simulation Threshold Mismatches;Crate/Module;Threshold",
             56, 0, 56, 16, 0, 16);
  setLabelsCMT(m_h_CPMHitsThreshSIMneDAT);

  } // end if (m_compareWithSim)

  // CMMHits

  m_monGroup = &monCMM;

  m_h_CMMHitsSIMeqDAT = book2F("CMMHitsSIMeqDAT",
     "CMM-CP Hits Data/Simulation Non-zero Matches;Module;Crate/Left-Right",
             14, 1, 15, 8, 0, 8);
  setLabelsMCLR(m_h_CMMHitsSIMeqDAT);
  m_h_CMMHitsSIMneDAT = book2F("CMMHitsSIMneDAT",
     "CMM-CP Hits Data/Simulation Non-zero Mismatches;Module;Crate/Left-Right",
             14, 1, 15, 8, 0, 8);
  setLabelsMCLR(m_h_CMMHitsSIMneDAT);
  m_h_CMMHitsSIMnoDAT = book2F("CMMHitsSIMnoDAT",
     "CMM-CP Hits Simulation but no Data;Module;Crate/Left-Right",
             14, 1, 15, 8, 0, 8);
  setLabelsMCLR(m_h_CMMHitsSIMnoDAT);
  m_h_CMMHitsDATnoSIM = book2F("CMMHitsDATnoSIM",
     "CMM-CP Hits Data but no Simulation;Module;Crate/Left-Right",
             14, 1, 15, 8, 0, 8);
  setLabelsMCLR(m_h_CMMHitsDATnoSIM);
  m_h_CMMHitsThreshSIMeqDAT = book2F("CMMHitsThreshSIMeqDAT",
     "CMM-CP Hits Data/Simulation Threshold Matches;Crate/Module;Threshold",
             56, 0, 56, 16, 0, 16);
  setLabelsCMT(m_h_CMMHitsThreshSIMeqDAT);
  m_h_CMMHitsThreshSIMneDAT = book2F("CMMHitsThreshSIMneDAT",
     "CMM-CP Hits Data/Simulation Threshold Mismatches;Crate/Module;Threshold",
             56, 0, 56, 16, 0, 16);
  setLabelsCMT(m_h_CMMHitsThreshSIMneDAT);

  if (m_compareWithSim) {

  // Local/Remote/Total sums

  m_h_SumsSIMeqDAT = book1F("SumsSIMeqDAT",
     "CMM-CP Hit Sums Data/Simulation Non-zero Matches;Sum/Left-Right",
             16, 0, 16);
  setLabelsSLR(m_h_SumsSIMeqDAT);
  m_h_SumsSIMneDAT = book1F("SumsSIMneDAT",
     "CMM-CP Hit Sums Data/Simulation Non-zero Mismatches;Sum/Left-Right",
             16, 0, 16);
  setLabelsSLR(m_h_SumsSIMneDAT);
  m_h_SumsSIMnoDAT = book1F("SumsSIMnoDAT",
     "CMM-CP Hit Sums Simulation but no Data;Sum/Left-Right",
             16, 0, 16);
  setLabelsSLR(m_h_SumsSIMnoDAT);
  m_h_SumsDATnoSIM = book1F("SumsDATnoSIM",
     "CMM-CP Hit Sums Data but no Simulation;Sum/Left-Right",
             16, 0, 16);
  setLabelsSLR(m_h_SumsDATnoSIM);
  m_h_SumsThreshSIMeqDAT = book2F("SumsThreshSIMeqDAT",
     "CMM-CP Hit Sums Data/Simulation Threshold Matches;Sum;Threshold",
             8, 0, 8, 16, 0, 16);
  setLabelsST(m_h_SumsThreshSIMeqDAT);
  m_h_SumsThreshSIMneDAT = book2F("SumsThreshSIMneDAT",
     "CMM-CP Hit Sums Data/Simulation Threshold Mismatches;Sum;Threshold",
             8, 0, 8, 16, 0, 16);
  setLabelsST(m_h_SumsThreshSIMneDAT);

  } else {

  // Remote sums only

  m_h_SumsSIMeqDAT = book1F("RemoteSIMeqDAT",
     "CMM-CP Remote Sums Data/Simulation Non-zero Matches;Sum/Left-Right",
             6, 0, 6);
  setLabelsSRLR(m_h_SumsSIMeqDAT);
  m_h_SumsSIMneDAT = book1F("RemoteSIMneDAT",
     "CMM-CP Remote Sums Data/Simulation Non-zero Mismatches;Sum/Left-Right",
             6, 0, 6);
  setLabelsSRLR(m_h_SumsSIMneDAT);
  m_h_SumsSIMnoDAT = book1F("RemoteSIMnoDAT",
     "CMM-CP Remote Sums Simulation but no Data;Sum/Left-Right",
             6, 0, 6);
  setLabelsSRLR(m_h_SumsSIMnoDAT);
  m_h_SumsDATnoSIM = book1F("RemoteDATnoSIM",
     "CMM-CP Remote Sums Data but no Simulation;Sum/Left-Right",
             6, 0, 6);
  setLabelsSRLR(m_h_SumsDATnoSIM);
  m_h_SumsThreshSIMeqDAT = book2F("RemoteThreshSIMeqDAT",
     "CMM-CP Remote Sums Data/Simulation Threshold Matches;Sum;Threshold",
             3, 0, 3, 16, 0, 16);
  setLabelsSRT(m_h_SumsThreshSIMeqDAT);
  m_h_SumsThreshSIMneDAT = book2F("RemoteThreshSIMneDAT",
     "CMM-CP Remote Sums Data/Simulation Threshold Mismatches;Sum;Threshold",
             3, 0, 3, 16, 0, 16);
  setLabelsSRT(m_h_SumsThreshSIMneDAT);

  } // end if (m_compareWithSim) else...

  // Summary

  m_monGroup = &monExpert;

  m_h_CPeqSIM = book2F("CPeqSIMOverview",
   "CP Comparison with Simulation Overview - Events with Matches;Crate/Module",
             64, 0, 64, 8, 0, 8);
  m_h_CPeqSIM->SetStats(kFALSE);
  setLabels(m_h_CPeqSIM);

  m_h_CPneSIM = book2F("CPneSIMOverview",
   "CP Comparison with Simulation Overview - Events with Mismatches;Crate/Module",
             64, 0, 64, 8, 0, 8);
  m_h_CPneSIM->SetStats(kFALSE);
  setLabels(m_h_CPneSIM);

  m_monGroup = &monShift;

  m_h_CPneSIMSummary = book1F("CPneSIMSummary",
   "CP Comparison with Simulation Mismatch Summary for 0 Events;;Events",
    7, 0, 7);
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(1, "Towers");
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(2, "RoIs");
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(3, "CPMHits");
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(4, "CMMHits");
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(5, "Local");
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(6, "Remote");
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(7, "Total");
  m_h_CPneSIMSummary->GetXaxis()->SetLabelSize(0.045);

  m_events = 0;

  } // end if (isNewRun ...
  
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode CPMSimBSMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  

  //Retrieve Trigger Towers from SG
  const TriggerTowerCollection* triggerTowerTES = 0; 
  StatusCode sc = m_storeGate->retrieve(triggerTowerTES,
                                                     m_triggerTowerLocation); 
  if( sc.isFailure()  ||  !triggerTowerTES ) {
    log << MSG::DEBUG<< "No Trigger Tower container found"<< endreq; 
  }

  //Retrieve CPM Towers from SG
  const CpmTowerCollection* cpmTowerTES = 0; 
  sc = m_storeGate->retrieve(cpmTowerTES, m_cpmTowerLocation); 
  if( sc.isFailure()  ||  !cpmTowerTES ) {
    log << MSG::DEBUG<< "No CPM Tower container found"<< endreq; 
  }
  
  //Retrieve CPM RoIs from SG
  const CpmRoiCollection* cpmRoiTES = 0;
  sc = m_storeGate->retrieve( cpmRoiTES, m_cpmRoiLocation);
  if( sc.isFailure()  ||  !cpmRoiTES  ||  cpmRoiTES->empty() ) {
    log << MSG::DEBUG << "No DAQ CPM RoIs found, trying RoIB"
        << endreq; 
    cpmRoiTES = 0;
    sc = m_storeGate->retrieve( cpmRoiTES, m_cpmRoiLocationRoib);
    if( sc.isFailure()  ||  !cpmRoiTES ) {
      log << MSG::DEBUG << "No RoIB CPM RoIs container found"<< endreq;
    }
  }
  
  //Retrieve CPM Hits from SG
  const CpmHitsCollection* cpmHitsTES = 0;
  sc = m_storeGate->retrieve( cpmHitsTES, m_cpmHitsLocation);
  if( sc.isFailure()  ||  !cpmHitsTES ) {
    log << MSG::DEBUG << "No CPM Hits container found"<< endreq; 
  }
  
  //Retrieve CMM-CP Hits from SG
  const CmmCpHitsCollection* cmmCpHitsTES = 0;
  sc = m_storeGate->retrieve( cmmCpHitsTES, m_cmmCpHitsLocation);
  if( sc.isFailure()  ||  !cmmCpHitsTES ) {
    log << MSG::DEBUG << "No CMM-CP Hits container found"<< endreq; 
  }

  // Maps to simplify comparisons
  
  TriggerTowerMap ttMap;
  CpmTowerMap     cpMap;
  CpmRoiMap       crMap;
  CpmHitsMap      chMap;
  CmmCpHitsMap    cmMap;
  setupMap(triggerTowerTES, ttMap);
  setupMap(cpmTowerTES, cpMap);
  if (m_compareWithSim) setupMap(cpmRoiTES, crMap);
  setupMap(cpmHitsTES, chMap);
  setupMap(cmmCpHitsTES, cmMap);

  // Summary error vectors
  const int nCrates = 4;
  const int nCPMs   = 14;
  const int nCMMs   = 2;
  const int vecsizeCpm = 2 * nCrates * nCPMs;
  const int vecsizeCmm = 2 * nCrates * nCMMs;
  ErrorVector errors1(vecsizeCpm);
  ErrorVector errors2(vecsizeCpm);
  ErrorVector errors3(vecsizeCpm);
  ErrorVector errors4(vecsizeCpm);
  ErrorVector errors5(vecsizeCpm);
  ErrorVector errors6(vecsizeCmm);
  ErrorVector errors7(vecsizeCmm);
  ErrorVector errors8(vecsizeCmm);
  ErrorVector errors9(vecsizeCmm);

  // Note - Simulation steps which are simply copies of data from
  // one container to another are not actually done to save time
  // and space.  The comparisons are made with the input instead.

  // Compare Trigger Towers and CPM Towers from data

  compare(ttMap, cpMap, errors1, errors2);

  if (m_compareWithSim) {

  // Compare RoIs simulated from CPM Towers with CPM RoIs from data

  CpmRoiCollection* cpmRoiSIM = 0;
  if (cpmTowerTES) {
    cpmRoiSIM = new CpmRoiCollection;
    simulate(cpMap, cpmRoiSIM);
  }
  CpmRoiMap crSimMap;
  setupMap(cpmRoiSIM, crSimMap);
  compare(crSimMap, crMap, errors3);
  crSimMap.clear();
  delete cpmRoiSIM;

  // Compare CPM Hits simulated from CPM RoIs with CPM Hits from data

  CpmHitsCollection* cpmHitsSIM = 0;
  if (cpmRoiTES) {
    cpmHitsSIM = new CpmHitsCollection;
    simulate(cpmRoiTES, cpmHitsSIM);
  }
  CpmHitsMap chSimMap;
  setupMap(cpmHitsSIM, chSimMap);
  compare(chSimMap, chMap, errors4);
  chSimMap.clear();
  delete cpmHitsSIM;

  } // end if (m_compareWithSim)

  // Compare CPM hits with CMM Hits from data

  compare(chMap, cmMap, errors5, errors6);

  if (m_compareWithSim) {

  // Compare Local sums simulated from CMM Hits with Local sums from data

  CmmCpHitsCollection* cmmLocalSIM = 0;
  if (cmmCpHitsTES) {
    cmmLocalSIM = new CmmCpHitsCollection;
    simulate(cmmCpHitsTES, cmmLocalSIM, LVL1::CMMCPHits::LOCAL);
  }
  CmmCpHitsMap cmmLocalSimMap;
  setupMap(cmmLocalSIM, cmmLocalSimMap);
  compare(cmmLocalSimMap, cmMap, errors7, LVL1::CMMCPHits::LOCAL);
  cmmLocalSimMap.clear();
  delete cmmLocalSIM;

  } // end if (m_compareWithSim)

  // Compare Local sums with Remote sums from data

  compare(cmMap, cmMap, errors8, LVL1::CMMCPHits::REMOTE_0);

  if (m_compareWithSim) {

  // Compare Total sums simulated from Remote sums with Total sums from data

  CmmCpHitsCollection* cmmTotalSIM = 0;
  if (cmmCpHitsTES) {
    cmmTotalSIM = new CmmCpHitsCollection;
    simulate(cmmCpHitsTES, cmmTotalSIM, LVL1::CMMCPHits::TOTAL);
  }
  CmmCpHitsMap cmmTotalSimMap;
  setupMap(cmmTotalSIM, cmmTotalSimMap);
  compare(cmmTotalSimMap, cmMap, errors9, LVL1::CMMCPHits::TOTAL);
  cmmTotalSimMap.clear();
  delete cmmTotalSIM;

  } // end if (m_compareWithSim)

  // Fill summary histograms

  std::vector<int> summary(7);
  const int cpmBins = nCrates * nCPMs;
  const int cmmBins = nCrates * nCMMs;
  for (int crate = 0; crate < nCrates; ++crate) {
    for (int module = 0; module < nCPMs; ++module) {
      int xBin = crate * nCPMs + module;
      int yBin = 0;
      m_h_CPeqSIM->Fill(xBin, yBin++, errors1[xBin]);
      m_h_CPeqSIM->Fill(xBin, yBin++, errors2[xBin]);
      m_h_CPeqSIM->Fill(xBin, yBin++, errors3[xBin]);
      m_h_CPeqSIM->Fill(xBin, yBin++, errors4[xBin]);
      m_h_CPeqSIM->Fill(xBin, yBin,   errors5[xBin]);
      if (module < nCMMs) {
        xBin = crate * nCMMs + module + cpmBins;
        int loc = crate * nCMMs + module;
        m_h_CPeqSIM->Fill(xBin, yBin++, errors6[loc]);
        m_h_CPeqSIM->Fill(xBin, yBin++, errors7[loc]);
        m_h_CPeqSIM->Fill(xBin, yBin++, errors8[loc]);
        m_h_CPeqSIM->Fill(xBin, yBin,   errors9[loc]);
        xBin = crate * nCPMs + module;
      }
      yBin = 0;
      m_h_CPneSIM->Fill(xBin, yBin++, errors1[xBin+cpmBins]);
      m_h_CPneSIM->Fill(xBin, yBin++, errors2[xBin+cpmBins]);
      m_h_CPneSIM->Fill(xBin, yBin++, errors3[xBin+cpmBins]);
      m_h_CPneSIM->Fill(xBin, yBin++, errors4[xBin+cpmBins]);
      m_h_CPneSIM->Fill(xBin, yBin,   errors5[xBin+cpmBins]);
      summary[0] += errors1[xBin+cpmBins];
      summary[0] += errors2[xBin+cpmBins];
      summary[1] += errors3[xBin+cpmBins];
      summary[2] += errors4[xBin+cpmBins];
      summary[3] += errors5[xBin+cpmBins];
      if (module < nCMMs) {
        xBin = crate * nCMMs + module + cpmBins;
        int loc = crate * nCMMs + module;
        m_h_CPneSIM->Fill(xBin, yBin++, errors6[loc+cmmBins]);
        m_h_CPneSIM->Fill(xBin, yBin++, errors7[loc+cmmBins]);
        m_h_CPneSIM->Fill(xBin, yBin++, errors8[loc+cmmBins]);
        m_h_CPneSIM->Fill(xBin, yBin,   errors9[loc+cmmBins]);
	summary[3] += errors6[loc+cmmBins];
	summary[4] += errors7[loc+cmmBins];
	summary[5] += errors8[loc+cmmBins];
	summary[6] += errors9[loc+cmmBins];
      }
    }
  }
  for (int i = 0; i < 7; ++i) m_h_CPneSIMSummary->Fill(i, summary[i] > 0);
  ++m_events;
  std::ostringstream cnum;
  cnum << m_events;
  std::string title("CP Comparison with Simulation Mismatch Summary for "
                    + cnum.str() + " Events");
  m_h_CPneSIMSummary->SetTitle(TString(title));

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode CPMSimBSMon::procHistograms(bool isEndOfEventsBlock,
                                  bool isEndOfLumiBlock, bool isEndOfRun)
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());

  if (isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun) {
  }

  return StatusCode::SUCCESS;
}

TH1F* CPMSimBSMon::book1F(std::string name, std::string title,
                                   int nx, double xmin, double xmax)
{
  TH1F *hist = new TH1F(TString(name), TString(title), nx, xmin, xmax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    MsgStream log(msgSvc(), this->name());
    log << MSG::WARNING << "Could not register histogram : " 
	<< name << endreq;
  }
  
  return hist;
}

TH2F* CPMSimBSMon::book2F(std::string name, std::string title,
                                   int nx, double xmin, double xmax,  
	                           int ny, double ymin, double ymax)
{		
  TH2F *hist = new TH2F(TString(name), TString(title), nx, xmin, xmax,
                                                       ny, ymin, ymax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    MsgStream log(msgSvc(), this->name());
    log << MSG::WARNING << "Could not register histogram : " 
	<< name << endreq;
  }
  hist->SetOption("colz");
  
  return hist;
}

//  Compare Trigger Towers and CPM Towers

void CPMSimBSMon::compare(const TriggerTowerMap& ttMap,
                          const CpmTowerMap& cpMap, ErrorVector& errors,
			  ErrorVector& errors2)
{

  const int maxKey = 0x7fffffff;
  TriggerTowerMap::const_iterator ttMapIter    = ttMap.begin();
  TriggerTowerMap::const_iterator ttMapIterEnd = ttMap.end();
  CpmTowerMap::const_iterator     cpMapIter    = cpMap.begin();
  CpmTowerMap::const_iterator     cpMapIterEnd = cpMap.end();

  while (ttMapIter != ttMapIterEnd || cpMapIter != cpMapIterEnd) {

    int ttKey = maxKey;
    int cpKey = maxKey;
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

    if (!ttEm && !ttHad && !cpEm && !cpHad) continue;
    
    //  Fill in error plots

    const LVL1::Coordinate coord(phi, eta);
    LVL1::CoordToHardware converter;
    const int crate = converter.cpCrate(coord);
    const int cpm   = converter.cpModule(coord);
    const int loc   = crate * 14 + cpm - 1;
    const int cpmBins = 4 * 14;
    if (ttEm && ttEm == cpEm)    errors[loc] = 1;
    if (ttHad && ttHad == cpHad) errors2[loc] = 1;
    if (ttEm != cpEm)   errors[loc+cpmBins] = 1;
    if (ttHad != cpHad) errors2[loc+cpmBins] = 1;
    double phiMod = phi * m_phiScale;
    m_h_EMTowerSIMeqDAT->Fill(eta, phiMod, ttEm && ttEm == cpEm);
    m_h_EMTowerSIMneDAT->Fill(eta, phiMod, ttEm && cpEm && ttEm != cpEm);
    m_h_EMTowerSIMnoDAT->Fill(eta, phiMod, ttEm && !cpEm);
    m_h_EMTowerDATnoSIM->Fill(eta, phiMod, cpEm && !ttEm);
    m_h_HadTowerSIMeqDAT->Fill(eta, phiMod, ttHad && ttHad == cpHad);
    m_h_HadTowerSIMneDAT->Fill(eta, phiMod, ttHad && cpHad && ttHad != cpHad);
    m_h_HadTowerSIMnoDAT->Fill(eta, phiMod, ttHad && !cpHad);
    m_h_HadTowerDATnoSIM->Fill(eta, phiMod, cpHad && !ttHad);
  }
}

//  Compare Simulated RoIs with data

void CPMSimBSMon::compare(const CpmRoiMap& roiSimMap, const CpmRoiMap& roiMap,
                                                      ErrorVector& errors)
{
  MsgStream log(msgSvc(), name());

  const int maxKey = 0xffff;
  LVL1::CPRoIDecoder decoder;
  CpmRoiMap::const_iterator simMapIter    = roiSimMap.begin();
  CpmRoiMap::const_iterator simMapIterEnd = roiSimMap.end();
  CpmRoiMap::const_iterator datMapIter    = roiMap.begin();
  CpmRoiMap::const_iterator datMapIterEnd = roiMap.end();

  while (simMapIter != simMapIterEnd || datMapIter != datMapIterEnd) {

    int simKey = maxKey;
    int datKey = maxKey;
    unsigned int simHits = 0;
    unsigned int datHits = 0;
    const LVL1::CPMRoI* roi = 0;

    if (simMapIter != simMapIterEnd) simKey = simMapIter->first;
    if (datMapIter != datMapIterEnd) datKey = datMapIter->first;

    if ((datMapIter == datMapIterEnd) || (datKey > simKey)) {

      // Simulated RoI but no data RoI

      roi     = simMapIter->second;
      simHits = roi->hits();
      ++simMapIter;

    } else if ((simMapIter == simMapIterEnd) || (simKey > datKey)) {

      // Data RoI but no simulated RoI

      roi     = datMapIter->second;
      datHits = roi->hits();
      ++datMapIter;

    } else {

      // Have both

      const LVL1::CPMRoI* roiS = simMapIter->second;
      roi     = datMapIter->second;
      simHits = roiS->hits();
      datHits = roi ->hits();
      ++simMapIter;
      ++datMapIter;
    }

    if (!simHits && !datHits) continue;
    
    //  Fill in error plots

    const int crate = roi->crate();
    const int cpm   = roi->cpm();
    const int chip  = roi->chip();
    const int local = roi->location();
    const int locX  = crate * 14 + cpm - 1;
    const int locY  = chip * 8 + local;
    const int cpmBins = 4 * 14;
    if (simHits && simHits == datHits) errors[locX] = 1;
    if (simHits != datHits) errors[locX+cpmBins] = 1;

    m_h_RoISIMeqDAT->Fill(locX, locY, simHits && simHits == datHits);
    m_h_RoISIMneDAT->Fill(locX, locY, simHits && datHits && simHits != datHits);
    m_h_RoISIMnoDAT->Fill(locX, locY, simHits && !datHits);
    m_h_RoIDATnoSIM->Fill(locX, locY, datHits && !simHits);
    for (int thr = 0; thr < 16; ++thr) {
      m_h_RoIThreshSIMeqDAT->Fill(locX, thr, ((datHits >> thr) & 0x1) ==
                                             ((simHits >> thr) & 0x1));
      m_h_RoIThreshSIMneDAT->Fill(locX, thr, ((datHits >> thr) & 0x1) !=
                                             ((simHits >> thr) & 0x1));
    }

    const LVL1::CoordinateRange coord(decoder.coordinate(roi->roiWord()));
    const double eta = coord.eta();
    const double phi = coord.phi() * m_phiScale;
    m_h_RoIEtaPhiSIMeqDAT->Fill(eta, phi, simHits && simHits == datHits);
    m_h_RoIEtaPhiSIMneDAT->Fill(eta, phi, simHits && datHits
                                                  && simHits != datHits);
    m_h_RoIEtaPhiSIMnoDAT->Fill(eta, phi, simHits && !datHits);
    m_h_RoIEtaPhiDATnoSIM->Fill(eta, phi, datHits && !simHits);

    log << MSG::DEBUG << "DataHits/SimHits: ";
    for (int i = 15; i >= 0; --i) {
      int bit = (datHits >> i) & 0x1;
      log << MSG::DEBUG << bit;
    }
    log << MSG::DEBUG << "/";
    for (int i = 15; i >= 0; --i) {
      int bit = (simHits >> i) & 0x1;
      log << MSG::DEBUG << bit;
    }
    log << MSG::DEBUG << endreq;
  }
}

//  Compare simulated CPM Hits with data

void CPMSimBSMon::compare(const CpmHitsMap& cpmSimMap, const CpmHitsMap& cpmMap,
                                                       ErrorVector& errors)
{

  const int maxKey = 0x7fffffff;
  CpmHitsMap::const_iterator simMapIter    = cpmSimMap.begin();
  CpmHitsMap::const_iterator simMapIterEnd = cpmSimMap.end();
  CpmHitsMap::const_iterator datMapIter    = cpmMap.begin();
  CpmHitsMap::const_iterator datMapIterEnd = cpmMap.end();

  while (simMapIter != simMapIterEnd || datMapIter != datMapIterEnd) {

    int simKey = maxKey;
    int datKey = maxKey;
    unsigned int simHits0 = 0;
    unsigned int simHits1 = 0;
    unsigned int datHits0 = 0;
    unsigned int datHits1 = 0;
    int crate = 0;
    int cpm   = 0;

    if (simMapIter != simMapIterEnd) simKey = simMapIter->first;
    if (datMapIter != datMapIterEnd) datKey = datMapIter->first;

    if ((datMapIter == datMapIterEnd) || (datKey > simKey)) {

      // Simulation Hits but no data Hits

      const LVL1::CPMHits* simh = simMapIter->second;
      simHits0 = simh->HitWord0();
      simHits1 = simh->HitWord1();
      crate    = simh->crate();
      cpm      = simh->module();
      ++simMapIter;

    } else if ((simMapIter == simMapIterEnd) || (simKey > datKey)) {

      // Data Hits but no simulation Hits

      const LVL1::CPMHits* dath = datMapIter->second;
      datHits0 = dath->HitWord0();
      datHits1 = dath->HitWord1();
      crate    = dath->crate();
      cpm      = dath->module();
      ++datMapIter;

    } else {

      // Have both

      const LVL1::CPMHits* simh = simMapIter->second;
      const LVL1::CPMHits* dath = datMapIter->second;
      simHits0 = simh->HitWord0();
      simHits1 = simh->HitWord1();
      datHits0 = dath->HitWord0();
      datHits1 = dath->HitWord1();
      crate    = dath->crate();
      cpm      = dath->module();
      ++simMapIter;
      ++datMapIter;
    }

    if (!simHits0 && !simHits1 && !datHits0 && !datHits1) continue;
    
    //  Fill in error plots

    const int loc = crate * 14 + cpm - 1;
    const int cpmBins = 4 * 14;
    if ((simHits0 && simHits0 == datHits0) ||
        (simHits1 && simHits1 == datHits1)) errors[loc] = 1;
    if (simHits0 != datHits0 || simHits1 != datHits1) errors[loc+cpmBins] = 1;

    m_h_CPMHitsSIMeqDAT->Fill(cpm, crate, (simHits0 || simHits1) &&
                              simHits0 == datHits0 && simHits1 == datHits1);
    m_h_CPMHitsSIMneDAT->Fill(cpm, crate, (simHits0 || simHits1) &&
      (datHits0 || datHits1) && (simHits0 != datHits0 || simHits1 != datHits1));
    m_h_CPMHitsSIMnoDAT->Fill(cpm, crate, (simHits0 || simHits1) &&
                                          !datHits0 && !datHits1);
    m_h_CPMHitsDATnoSIM->Fill(cpm, crate, (datHits0 || datHits1) &&
                                          !simHits0 && !simHits1);

    for (int thr = 0; thr < 8; ++thr) {
      m_h_CPMHitsThreshSIMeqDAT->Fill(loc, thr,
	           ((datHits0 >> 3*thr) & 0x7) == ((simHits0 >> 3*thr) & 0x7));
      m_h_CPMHitsThreshSIMeqDAT->Fill(loc, thr + 8,
	           ((datHits1 >> 3*thr) & 0x7) == ((simHits1 >> 3*thr) & 0x7));
      m_h_CPMHitsThreshSIMneDAT->Fill(loc, thr,
	           ((datHits0 >> 3*thr) & 0x7) != ((simHits0 >> 3*thr) & 0x7));
      m_h_CPMHitsThreshSIMneDAT->Fill(loc, thr + 8,
	           ((datHits1 >> 3*thr) & 0x7) != ((simHits1 >> 3*thr) & 0x7));
    }
  }
}

//  Compare CPM Hits and CMM Hits

void CPMSimBSMon::compare(const CpmHitsMap& cpmMap, const CmmCpHitsMap& cmmMap,
                          ErrorVector& errors, ErrorVector& errors2)
{

  const int maxKey = 0x7fffffff;
  CpmHitsMap::const_iterator   cpmMapIter    = cpmMap.begin();
  CpmHitsMap::const_iterator   cpmMapIterEnd = cpmMap.end();
  CmmCpHitsMap::const_iterator cmmMapIter    = cmmMap.begin();
  CmmCpHitsMap::const_iterator cmmMapIterEnd = cmmMap.end();

  while (cpmMapIter != cpmMapIterEnd || cmmMapIter != cmmMapIterEnd) {

    int cpmKey = maxKey;
    int cmmKey = maxKey;
    unsigned int cpmHits0 = 0;
    unsigned int cpmHits1 = 0;
    unsigned int cmmHits0 = 0;
    unsigned int cmmHits1 = 0;
    int crate  = 0;
    int cpm    = 0;

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
      if (cpm > 14) continue;

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

    if (!cpmHits0 && !cpmHits1 && !cmmHits0 && !cmmHits1) continue;
    
    //  Fill in error plots

    const int loc  = crate * 14 + cpm - 1;
    const int loc2 = crate * 2;
    const int cpmBins = 4 * 14;
    const int cmmBins = 4 * 2;
    if ((cpmHits0 && cpmHits0 == cmmHits0) ||
        (cpmHits1 && cpmHits1 == cmmHits1)) errors[loc] = 1;
    if (cpmHits0 != cmmHits0 || cpmHits1 != cmmHits1) errors[loc+cpmBins] = 1;
    if (cpmHits1 && cpmHits1 == cmmHits1) errors2[loc2]   = 1;
    if (cpmHits0 && cpmHits0 == cmmHits0) errors2[loc2+1] = 1;
    if (cpmHits1 != cmmHits1) errors2[loc2+cmmBins]   = 1; // hits1==>cmm 0
    if (cpmHits0 != cmmHits0) errors2[loc2+cmmBins+1] = 1; // hits0==>cmm 1

    m_h_CMMHitsSIMeqDAT->Fill(cpm, 2*crate,   cpmHits1 && cpmHits1 == cmmHits1);
    m_h_CMMHitsSIMeqDAT->Fill(cpm, 2*crate+1, cpmHits0 && cpmHits0 == cmmHits0);
    m_h_CMMHitsSIMneDAT->Fill(cpm, 2*crate,   cpmHits1 && cmmHits1 &&
                                              cpmHits1 != cmmHits1);
    m_h_CMMHitsSIMneDAT->Fill(cpm, 2*crate+1, cpmHits0 && cmmHits0 &&
                                              cpmHits0 != cmmHits0);
    m_h_CMMHitsSIMnoDAT->Fill(cpm, 2*crate,   cpmHits1 && !cmmHits1);
    m_h_CMMHitsSIMnoDAT->Fill(cpm, 2*crate+1, cpmHits0 && !cmmHits0);
    m_h_CMMHitsDATnoSIM->Fill(cpm, 2*crate,   cmmHits1 && !cpmHits1);
    m_h_CMMHitsDATnoSIM->Fill(cpm, 2*crate+1, cmmHits0 && !cpmHits0);

    for (int thr = 0; thr < 8; ++thr) {
      m_h_CMMHitsThreshSIMeqDAT->Fill(loc, thr,
	           ((cmmHits0 >> 3*thr) & 0x7) == ((cpmHits0 >> 3*thr) & 0x7));
      m_h_CMMHitsThreshSIMeqDAT->Fill(loc, thr + 8,
	           ((cmmHits1 >> 3*thr) & 0x7) == ((cpmHits1 >> 3*thr) & 0x7));
      m_h_CMMHitsThreshSIMneDAT->Fill(loc, thr,
	           ((cmmHits0 >> 3*thr) & 0x7) != ((cpmHits0 >> 3*thr) & 0x7));
      m_h_CMMHitsThreshSIMneDAT->Fill(loc, thr + 8,
	           ((cmmHits1 >> 3*thr) & 0x7) != ((cpmHits1 >> 3*thr) & 0x7));
    }
  }
}

//  Compare Simulated CMM Hit Sums and Data CMM Hit Sums

void CPMSimBSMon::compare(const CmmCpHitsMap& cmmSimMap,
                          const CmmCpHitsMap& cmmMap,
                          ErrorVector& errors, int selection)
{

  const bool local  = (selection == LVL1::CMMCPHits::LOCAL);
  const bool remote = (selection == LVL1::CMMCPHits::REMOTE_0);
  const bool total  = (selection == LVL1::CMMCPHits::TOTAL);
  if (!local && !remote && !total) return;
  std::vector<unsigned int> hits0Sim(3);
  std::vector<unsigned int> hits1Sim(3);
  std::vector<unsigned int> hits0(3);
  std::vector<unsigned int> hits1(3);
  const int maxKey = 0x7fffffff;
  CmmCpHitsMap::const_iterator cmmSimMapIter    = cmmSimMap.begin();
  CmmCpHitsMap::const_iterator cmmSimMapIterEnd = cmmSimMap.end();
  CmmCpHitsMap::const_iterator cmmMapIter       = cmmMap.begin();
  CmmCpHitsMap::const_iterator cmmMapIterEnd    = cmmMap.end();

  while (cmmSimMapIter != cmmSimMapIterEnd || cmmMapIter != cmmMapIterEnd) {

    int cmmSimKey = maxKey;
    int cmmKey    = maxKey;
    unsigned int cmmSimHits0 = 0;
    unsigned int cmmSimHits1 = 0;
    unsigned int cmmHits0 = 0;
    unsigned int cmmHits1 = 0;
    int crate  = 0;
    int dataId = 0;

    if (cmmSimMapIter != cmmSimMapIterEnd) cmmSimKey = cmmSimMapIter->first;
    if (cmmMapIter    != cmmMapIterEnd)    cmmKey    = cmmMapIter->first;

    if ((cmmMapIter == cmmMapIterEnd) || (cmmKey > cmmSimKey)) {

      // Sim CMM Hits but no Data CMM Hits

      const LVL1::CMMCPHits* cmmS = cmmSimMapIter->second;
      ++cmmSimMapIter;
      dataId = cmmS->dataID();
      if (local  && dataId != LVL1::CMMCPHits::LOCAL) continue;
      if (remote && dataId != LVL1::CMMCPHits::LOCAL) continue;
      if (total  && dataId != LVL1::CMMCPHits::TOTAL) continue;
      cmmSimHits0 = cmmS->HitWord0();
      cmmSimHits1 = cmmS->HitWord1();
      crate       = cmmS->crate();

    } else if ((cmmSimMapIter == cmmSimMapIterEnd) || (cmmSimKey > cmmKey)) {

      // Data CMM Hits but no Sim CMM Hits

      const LVL1::CMMCPHits* cmmD = cmmMapIter->second;
      ++cmmMapIter;
      dataId   = cmmD->dataID();
      if (local  && dataId != LVL1::CMMCPHits::LOCAL)    continue;
      if (remote && dataId != LVL1::CMMCPHits::REMOTE_0 &&
                    dataId != LVL1::CMMCPHits::REMOTE_1 &&
		    dataId != LVL1::CMMCPHits::REMOTE_2) continue;
      if (total  && dataId != LVL1::CMMCPHits::TOTAL)    continue;
      cmmHits0 = cmmD->HitWord0();
      cmmHits1 = cmmD->HitWord1();
      crate    = cmmD->crate();

    } else {

      // Have both

      const LVL1::CMMCPHits* cmmS = cmmSimMapIter->second;
      const LVL1::CMMCPHits* cmmD = cmmMapIter->second;
      ++cmmSimMapIter;
      ++cmmMapIter;
      dataId   = cmmS->dataID();
      if (local  && dataId != LVL1::CMMCPHits::LOCAL)    continue;
      if (remote && dataId != LVL1::CMMCPHits::LOCAL    &&
                    dataId != LVL1::CMMCPHits::REMOTE_0 &&
                    dataId != LVL1::CMMCPHits::REMOTE_1 &&
		    dataId != LVL1::CMMCPHits::REMOTE_2) continue;
      if (total  && dataId != LVL1::CMMCPHits::TOTAL)    continue;
      cmmSimHits0 = cmmS->HitWord0();
      cmmSimHits1 = cmmS->HitWord1();
      cmmHits0    = cmmD->HitWord0();
      cmmHits1    = cmmD->HitWord1();
      crate       = cmmS->crate();
    }

    if (!cmmSimHits0 && !cmmSimHits1 && !cmmHits0 && !cmmHits1) continue;
    
    //  Fill in error plots

    if (local || total) {
      int loc = crate * 2;
      const int cmmBins = 4 * 2;
      if (cmmSimHits1 && cmmSimHits1 == cmmHits1) errors[loc]   = 1;
      if (cmmSimHits0 && cmmSimHits0 == cmmHits0) errors[loc+1] = 1;
      if (cmmSimHits1 != cmmHits1) errors[loc+cmmBins]   = 1; // hits1==>cmm 0
      if (cmmSimHits0 != cmmHits0) errors[loc+cmmBins+1] = 1; // hits0==>cmm 1
      loc = (local) ? loc : 14;
      m_h_SumsSIMeqDAT->Fill(loc,   cmmSimHits1 && cmmSimHits1 == cmmHits1);
      m_h_SumsSIMeqDAT->Fill(loc+1, cmmSimHits0 && cmmSimHits0 == cmmHits0);
      m_h_SumsSIMneDAT->Fill(loc,   cmmSimHits1 && cmmHits1 &&
                                                    cmmSimHits1 != cmmHits1);
      m_h_SumsSIMneDAT->Fill(loc+1, cmmSimHits0 && cmmHits0 &&
                                                    cmmSimHits0 != cmmHits0);
      m_h_SumsSIMnoDAT->Fill(loc,   cmmSimHits1 && !cmmHits1);
      m_h_SumsSIMnoDAT->Fill(loc+1, cmmSimHits0 && !cmmHits0);
      m_h_SumsDATnoSIM->Fill(loc,   cmmHits1 && !cmmSimHits1);
      m_h_SumsDATnoSIM->Fill(loc+1, cmmHits0 && !cmmSimHits0);
      loc /= 2;
      for (int thr = 0; thr < 8; ++thr) {
        m_h_SumsThreshSIMeqDAT->Fill(loc, thr,
	        ((cmmSimHits0 >> 3*thr) & 0x7) == ((cmmHits0 >> 3*thr) & 0x7));
        m_h_SumsThreshSIMeqDAT->Fill(loc, thr + 8,
	        ((cmmSimHits1 >> 3*thr) & 0x7) == ((cmmHits1 >> 3*thr) & 0x7));
        m_h_SumsThreshSIMneDAT->Fill(loc, thr,
	        ((cmmSimHits0 >> 3*thr) & 0x7) != ((cmmHits0 >> 3*thr) & 0x7));
        m_h_SumsThreshSIMneDAT->Fill(loc, thr + 8,
	        ((cmmSimHits1 >> 3*thr) & 0x7) != ((cmmHits1 >> 3*thr) & 0x7));
      }
    } else {
      if (dataId == LVL1::CMMCPHits::LOCAL) {
        if (crate != 3) {
	  hits0Sim[crate] = cmmSimHits0;
	  hits1Sim[crate] = cmmSimHits1;
        }
      } else if (dataId == LVL1::CMMCPHits::REMOTE_0) {
        hits0[0] = cmmHits0;
        hits1[0] = cmmHits1;
      } else if (dataId == LVL1::CMMCPHits::REMOTE_1) {
        hits0[1] = cmmHits0;
        hits1[1] = cmmHits1;
      } else {
        hits0[2] = cmmHits0;
        hits1[2] = cmmHits1;
      }
    }
  }
  if (remote) {
    for (int crate = 0; crate < 3; ++crate) {
      int loc = crate * 2;
      const int cmmBins = 4 * 2;
      if (hits1Sim[crate] && hits1Sim[crate] == hits1[crate]) errors[loc]   = 1;
      if (hits0Sim[crate] && hits0Sim[crate] == hits0[crate]) errors[loc+1] = 1;
      if (hits1Sim[crate] != hits1[crate]) errors[loc+cmmBins]   = 1;
      if (hits0Sim[crate] != hits0[crate]) errors[loc+cmmBins+1] = 1;
      loc = (m_compareWithSim) ? loc + 8 : loc;
      m_h_SumsSIMeqDAT->Fill(loc,   hits1Sim[crate] &&
                                    hits1Sim[crate] == hits1[crate]);
      m_h_SumsSIMeqDAT->Fill(loc+1, hits0Sim[crate] &&
                                    hits0Sim[crate] == hits0[crate]);
      m_h_SumsSIMneDAT->Fill(loc,   hits1Sim[crate] && hits1[crate] &&
                                    hits1Sim[crate] != hits1[crate]);
      m_h_SumsSIMneDAT->Fill(loc+1, hits0Sim[crate] && hits0[crate] &&
                                    hits0Sim[crate] != hits0[crate]);
      m_h_SumsSIMnoDAT->Fill(loc,   hits1Sim[crate] && !hits1[crate]);
      m_h_SumsSIMnoDAT->Fill(loc+1, hits0Sim[crate] && !hits0[crate]);
      m_h_SumsDATnoSIM->Fill(loc,   hits1[crate] && !hits1Sim[crate]);
      m_h_SumsDATnoSIM->Fill(loc+1, hits0[crate] && !hits0Sim[crate]);
      loc /= 2;
      for (int thr = 0; thr < 8; ++thr) {
        m_h_SumsThreshSIMeqDAT->Fill(loc, thr,
	 ((hits0Sim[crate] >> 3*thr) & 0x7) == ((hits0[crate] >> 3*thr) & 0x7));
        m_h_SumsThreshSIMeqDAT->Fill(loc, thr + 8,
	 ((hits1Sim[crate] >> 3*thr) & 0x7) == ((hits1[crate] >> 3*thr) & 0x7));
        m_h_SumsThreshSIMneDAT->Fill(loc, thr,
	 ((hits0Sim[crate] >> 3*thr) & 0x7) != ((hits0[crate] >> 3*thr) & 0x7));
        m_h_SumsThreshSIMneDAT->Fill(loc, thr + 8,
	 ((hits1Sim[crate] >> 3*thr) & 0x7) != ((hits1[crate] >> 3*thr) & 0x7));
      }
    }
  }
}

void CPMSimBSMon::setLabels(TH2* hist)
{
  setLabelsCPM(hist);
  hist->GetXaxis()->SetBinLabel(1, "CPM");
  hist->GetXaxis()->SetBinLabel(57, "CMM");
  hist->GetXaxis()->SetBinLabel(59, "1/L");
  hist->GetXaxis()->SetBinLabel(61, "2/L");
  hist->GetXaxis()->SetBinLabel(63, "3/L");
  hist->GetXaxis()->SetTitleOffset(1.25);
  int bin = 1;
  hist->GetYaxis()->SetBinLabel(bin++, "EM tt");
  hist->GetYaxis()->SetBinLabel(bin++, "Had tt");
  if (m_compareWithSim) {
    hist->GetYaxis()->SetBinLabel(bin++, "RoIs");
    hist->GetYaxis()->SetBinLabel(bin++, "CPMHits");
  } else {
    hist->GetYaxis()->SetBinLabel(bin++, "n/a");
    hist->GetYaxis()->SetBinLabel(bin++, "n/a");
  }
  hist->GetYaxis()->SetBinLabel(bin++, "CMMHits");
  if (m_compareWithSim) {
    hist->GetYaxis()->SetBinLabel(bin++, "Local");
  } else {
    hist->GetYaxis()->SetBinLabel(bin++, "n/a");
  }
  hist->GetYaxis()->SetBinLabel(bin++, "Remote");
  if (m_compareWithSim) {
    hist->GetYaxis()->SetBinLabel(bin++, "Total");
  } else {
    hist->GetYaxis()->SetBinLabel(bin++, "n/a");
  }
  hist->GetYaxis()->SetLabelSize(0.045);
}

void CPMSimBSMon::setLabelsCMCC(TH2* hist)
{
  setLabelsCPM(hist);
  for (int chip = 0; chip < 8; ++chip) {
    for (int loc = 0; loc < 8; loc += 4) {
      std::ostringstream cnum;
      cnum << chip << "/" << loc;
      hist->GetYaxis()->SetBinLabel(chip*8 + loc + 1, cnum.str().c_str());
    }
  }
}

void CPMSimBSMon::setLabelsCMT(TH2* hist)
{
  setLabelsCPM(hist);
  setLabelsT(hist);
}

void CPMSimBSMon::setLabelsT(TH2* hist)
{
  for (int thresh = 0; thresh < 16; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    hist->GetYaxis()->SetBinLabel(thresh + 1, cnum.str().c_str());
  }
  hist->GetYaxis()->SetLabelSize(0.05);
}

void CPMSimBSMon::setLabelsCPM(TH2* hist)
{
  const int nCPMs = 14;
  for (int crate = 0; crate < 4; ++crate) {
    for (int module = 1; module <= 14; module += 7) {
      std::ostringstream cnum;
      cnum << crate << "/" << module;
      hist->GetXaxis()->SetBinLabel(crate*nCPMs + module, cnum.str().c_str());
    }
  }
}

void CPMSimBSMon::setLabelsMC(TH2* hist)
{
  for (int module = 1; module <= 14; ++module) {
    std::ostringstream cnum;
    cnum << module;
    hist->GetXaxis()->SetBinLabel(module, cnum.str().c_str());
  }
  hist->GetXaxis()->SetLabelSize(0.05);
  for (int crate = 0; crate < 4; ++crate) {
    std::ostringstream cnum;
    cnum << crate;
    hist->GetYaxis()->SetBinLabel(crate + 1, cnum.str().c_str());
  }
  hist->GetYaxis()->SetLabelSize(0.05);
}

void CPMSimBSMon::setLabelsMCLR(TH2* hist)
{
  for (int module = 1; module <= 14; ++module) {
    std::ostringstream cnum;
    cnum << module;
    hist->GetXaxis()->SetBinLabel(module, cnum.str().c_str());
  }
  hist->GetXaxis()->SetLabelSize(0.05);
  for (int crate = 0; crate < 4; ++crate) {
    for (int cmm = 0; cmm < 2; ++cmm) {
      std::ostringstream cnum;
      if (cmm == 0) cnum << crate << "/L";
      else          cnum << crate << "/R";
      hist->GetYaxis()->SetBinLabel(crate*2 + cmm + 1, cnum.str().c_str());
    }
  }
  hist->GetYaxis()->SetLabelSize(0.05);
}

void CPMSimBSMon::setLabelsSLR(TH1* hist)
{
  hist->GetXaxis()->SetBinLabel(1, "L0/L");
  hist->GetXaxis()->SetBinLabel(2, "L0/R");
  hist->GetXaxis()->SetBinLabel(3, "L1/L");
  hist->GetXaxis()->SetBinLabel(4, "L1/R");
  hist->GetXaxis()->SetBinLabel(5, "L2/L");
  hist->GetXaxis()->SetBinLabel(6, "L2/R");
  hist->GetXaxis()->SetBinLabel(7, "L3/L");
  hist->GetXaxis()->SetBinLabel(8, "L3/R");
  hist->GetXaxis()->SetBinLabel(9, "R0/L");
  hist->GetXaxis()->SetBinLabel(10, "R0/R");
  hist->GetXaxis()->SetBinLabel(11, "R1/L");
  hist->GetXaxis()->SetBinLabel(12, "R1/R");
  hist->GetXaxis()->SetBinLabel(13, "R2/L");
  hist->GetXaxis()->SetBinLabel(14, "R2/R");
  hist->GetXaxis()->SetBinLabel(15, "T/L");
  hist->GetXaxis()->SetBinLabel(16, "T/R");
}

void CPMSimBSMon::setLabelsST(TH2* hist)
{
  hist->GetXaxis()->SetBinLabel(1, "L0");
  hist->GetXaxis()->SetBinLabel(2, "L1");
  hist->GetXaxis()->SetBinLabel(3, "L2");
  hist->GetXaxis()->SetBinLabel(4, "L3");
  hist->GetXaxis()->SetBinLabel(5, "R0");
  hist->GetXaxis()->SetBinLabel(6, "R1");
  hist->GetXaxis()->SetBinLabel(7, "R2");
  hist->GetXaxis()->SetBinLabel(8, "T");
  hist->GetXaxis()->SetLabelSize(0.05);
  setLabelsT(hist);
}

void CPMSimBSMon::setLabelsSRLR(TH1* hist)
{
  hist->GetXaxis()->SetBinLabel(1, "R0/L");
  hist->GetXaxis()->SetBinLabel(2, "R0/R");
  hist->GetXaxis()->SetBinLabel(3, "R1/L");
  hist->GetXaxis()->SetBinLabel(4, "R1/R");
  hist->GetXaxis()->SetBinLabel(5, "R2/L");
  hist->GetXaxis()->SetBinLabel(6, "R2/R");
  hist->GetXaxis()->SetLabelSize(0.05);
}

void CPMSimBSMon::setLabelsSRT(TH2* hist)
{
  hist->GetXaxis()->SetBinLabel(1, "L0");
  hist->GetXaxis()->SetBinLabel(2, "L1");
  hist->GetXaxis()->SetBinLabel(3, "L2");
  hist->GetXaxis()->SetLabelSize(0.05);
  setLabelsT(hist);
}

void CPMSimBSMon::setupMap(const TriggerTowerCollection* coll,
                                 TriggerTowerMap& map)
{
  if (coll) {
    LVL1::TriggerTowerKey towerKey;
    TriggerTowerCollection::const_iterator pos  = coll->begin();
    TriggerTowerCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const double eta = (*pos)->eta();
      if (eta > -2.5 && eta < 2.5) {
        const double phi = (*pos)->phi();
        const int key = towerKey.ttKey(phi, eta);
        map.insert(std::make_pair(key, *pos));
      }
    }
  }
}

void CPMSimBSMon::setupMap(const CpmTowerCollection* coll, CpmTowerMap& map)
{
  if (coll) {
    LVL1::TriggerTowerKey towerKey;
    CpmTowerCollection::const_iterator pos  = coll->begin();
    CpmTowerCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const double eta = (*pos)->eta();
      const double phi = (*pos)->phi();
      const int key = towerKey.ttKey(phi, eta);
      map.insert(std::make_pair(key, *pos));
    }
  }
}

void CPMSimBSMon::setupMap(const CpmRoiCollection* coll, CpmRoiMap& map)
{
  if (coll) {
    CpmRoiCollection::const_iterator pos  = coll->begin();
    CpmRoiCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const int crate = (*pos)->crate();
      const int cpm   = (*pos)->cpm();
      const int chip  = (*pos)->chip();
      const int loc   = (*pos)->location();
      const int key   = (((((crate << 4) | cpm) << 3) | chip) << 3) | loc;
      map.insert(std::make_pair(key, *pos));
    }
  }
}

void CPMSimBSMon::setupMap(const CpmHitsCollection* coll, CpmHitsMap& map)
{
  if (coll) {
    CpmHitsCollection::const_iterator pos  = coll->begin();
    CpmHitsCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const int crate = (*pos)->crate();
      const int cpm   = (*pos)->module();
      const int key   = crate * 100 + cpm;
      map.insert(std::make_pair(key, *pos));
    }
  }
}

void CPMSimBSMon::setupMap(const CmmCpHitsCollection* coll, CmmCpHitsMap& map)
{
  if (coll) {
    CmmCpHitsCollection::const_iterator pos  = coll->begin();
    CmmCpHitsCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const int crate  = (*pos)->crate();
      const int dataId = (*pos)->dataID();
      const int key  = crate * 100 + dataId;
      map.insert(std::make_pair(key, *pos));
    }
  }
}

void CPMSimBSMon::simulate(const CpmTowerMap towers,
                                 CpmRoiCollection* rois)
{
  InternalRoiCollection* intRois = new InternalRoiCollection;
  m_emTauTool->findRoIs(&towers, intRois);
  m_cpHitsTool->formCPMRoI(intRois, rois);
  delete intRois;
}

void CPMSimBSMon::simulate(const CpmRoiCollection* rois,
                                 CpmHitsCollection* hits)
{
  m_cpHitsTool->formCPMHits(rois, hits);
}

void CPMSimBSMon::simulate(const CmmCpHitsCollection* hitsIn,
                                 CmmCpHitsCollection* hitsOut, int selection)
{
  if (selection == LVL1::CMMCPHits::LOCAL) {
    m_cpHitsTool->formCMMCPHitsCrate(hitsIn, hitsOut);
  } else if (selection == LVL1::CMMCPHits::TOTAL) {
    m_cpHitsTool->formCMMCPHitsSystem(hitsIn, hitsOut);
  }
}
