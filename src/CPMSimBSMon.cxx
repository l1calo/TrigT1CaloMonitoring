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
    m_log(msgSvc(), name),
    m_monGroup(0), m_events(0)
/*---------------------------------------------------------*/
{
  declareInterface<IMonitorToolBase>(this); 

  declareProperty("CPMTowerLocation",
                 m_cpmTowerLocation  = LVL1::TrigT1CaloDefs::CPMTowerLocation);
  declareProperty("CPMTowerLocationOverlap",
                 m_cpmTowerLocationOverlap  =
		           LVL1::TrigT1CaloDefs::CPMTowerLocation+"Overlap");
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
                  "Include the checks that run simulation tools");
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
  m_log.setLevel(outputLevel());

  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;

  sc = m_storeGate.retrieve();
  if( sc.isFailure() ) {
    m_log << MSG::ERROR << "Unable to locate Service StoreGateSvc" << endreq;
    return sc;
  }

  if (m_compareWithSim) {

    sc = m_emTauTool.retrieve();
    if( sc.isFailure() ) {
      m_log << MSG::ERROR << "Unable to locate Tool L1EmTauTools" << endreq;
      return sc;
    }

    sc = m_cpHitsTool.retrieve();
    if( sc.isFailure() ) {
      m_log << MSG::ERROR << "Unable to locate Tool L1CPHitsTools" << endreq;
      return sc;
    }
  }

  // Phi units
  const double twoPi = 2.*M_PI;
  if      (m_phiUnits == "radians")  m_phiMax = twoPi;
  else if (m_phiUnits == "degrees")  m_phiMax = 360.;
  else if (m_phiUnits == "channels") m_phiMax = 64.;
  else {
    m_log << MSG::ERROR << "Invalid PhiUnits: " << m_phiUnits
          << ", using radians" << endreq;
    m_phiMax = twoPi;
  }
  m_phiScale = m_phiMax/twoPi;

  m_log << MSG::INFO << "CPMSimBSMon initialised" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode CPMSimBSMon::bookHistograms(bool isNewEventsBlock,
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

  std::string dir(m_rootDir + "/5_CP_Sim_Check");
  MonGroup monShift ( this, dir, shift, run );
  MonGroup monExpert( this, dir, expert, run );
  MonGroup monRoIs  ( this, dir + "/RoIs", expert, run );
  MonGroup monCPM   ( this, dir + "/CPM",  expert, run );
  MonGroup monCMM   ( this, dir + "/CMM",  expert, run );

  // CPMTowers

  m_monGroup = &monCPM;

  m_h_EMTowerSIMeqDAT = book2F("EMTowerSIMeqDAT",
    "Core CPM Towers EM Data/Simulation Non-zero Matches;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_EMTowerSIMneDAT = book2F("EMTowerSIMneDAT",
    "Core CPM Towers EM Data/Simulation Non-zero Mismatches;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_EMTowerSIMnoDAT = book2F("EMTowerSIMnoDAT",
    "Core CPM Towers EM Simulation but no Data;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_EMTowerDATnoSIM = book2F("EMTowerDATnoSIM",
    "Core CPM Towers EM Data but no Simulation;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_HadTowerSIMeqDAT = book2F("HadTowerSIMeqDAT",
    "Core CPM Towers HAD Data/Simulation Non-zero Matches;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_HadTowerSIMneDAT = book2F("HadTowerSIMneDAT",
    "Core CPM Towers HAD Data/Simulation Non-zero Mismatches;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_HadTowerSIMnoDAT = book2F("HadTowerSIMnoDAT",
    "Core CPM Towers HAD Simulation but no Data;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_HadTowerDATnoSIM = book2F("HadTowerDATnoSIM",
    "Core CPM Towers HAD Data but no Simulation;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_EMTowerOvSIMeqDAT = book2F("EMTowerOvSIMeqDAT",
    "Overlap CPM Towers EM Data/Simulation Non-zero Matches;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_EMTowerOvSIMneDAT = book2F("EMTowerOvSIMneDAT",
    "Overlap CPM Towers EM Data/Simulation Non-zero Mismatches;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_EMTowerOvSIMnoDAT = book2F("EMTowerOvSIMnoDAT",
    "Overlap CPM Towers EM Simulation but no Data;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_EMTowerOvDATnoSIM = book2F("EMTowerOvDATnoSIM",
    "Overlap CPM Towers EM Data but no Simulation;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_HadTowerOvSIMeqDAT = book2F("HadTowerOvSIMeqDAT",
    "Overlap CPM Towers HAD Data/Simulation Non-zero Matches;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_HadTowerOvSIMneDAT = book2F("HadTowerOvSIMneDAT",
    "Overlap CPM Towers HAD Data/Simulation Non-zero Mismatches;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_HadTowerOvSIMnoDAT = book2F("HadTowerOvSIMnoDAT",
    "Overlap CPM Towers HAD Simulation but no Data;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_HadTowerOvDATnoSIM = book2F("HadTowerOvDATnoSIM",
    "Overlap CPM Towers HAD Data but no Simulation;eta;phi",
             50, -2.5, 2.5, 64, 0, m_phiMax);
  m_h_FpgaTowerSIMeqDAT = book2F("FpgaTowerSIMeqDAT",
    "CPM Towers Data/Simulation Non-zero Matches by FPGA;Crate/Module;Serialiser FPGA",
             56, 0, 56, 20, 0, 20);
  setLabelsCPMFP(m_h_FpgaTowerSIMeqDAT);
  m_h_FpgaTowerSIMneDAT = book2F("FpgaTowerSIMneDAT",
    "CPM Towers Data/Simulation Non-zero Mismatches by FPGA;Crate/Module;Serialiser FPGA",
             56, 0, 56, 20, 0, 20);
  setLabelsCPMFP(m_h_FpgaTowerSIMneDAT);
  m_h_FpgaTowerSIMnoDAT = book2F("FpgaTowerSIMnoDAT",
    "CPM Towers Simulation but no Data by FPGA;Crate/Module;Serialiser FPGA",
             56, 0, 56, 20, 0, 20);
  setLabelsCPMFP(m_h_FpgaTowerSIMnoDAT);
  m_h_FpgaTowerDATnoSIM = book2F("FpgaTowerDATnoSIM",
    "CPM Towers Data but no Simulation by FPGA;Crate/Module;Serialiser FPGA",
             56, 0, 56, 20, 0, 20);
  setLabelsCPMFP(m_h_FpgaTowerDATnoSIM);

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
             64, 0, 64, NumberOfSummaryBins, 0, NumberOfSummaryBins);
  m_h_CPeqSIM->SetStats(kFALSE);
  setLabels(m_h_CPeqSIM);

  m_h_CPneSIM = book2F("CPneSIMOverview",
   "CP Comparison with Simulation Overview - Events with Mismatches;Crate/Module",
             64, 0, 64, NumberOfSummaryBins, 0, NumberOfSummaryBins);
  m_h_CPneSIM->SetStats(kFALSE);
  setLabels(m_h_CPneSIM);

  m_monGroup = &monShift;

  m_h_CPneSIMSummary = book1F("CPneSIMSummary",
   "CP Comparison with Simulation Mismatch Summary for 0 Events;;Events",
    NumberOfSummaryBins, 0, NumberOfSummaryBins);
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(1+EMTowerMismatch,   "EM tt");
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(1+HadTowerMismatch,  "Had tt");
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(1+RoIMismatch,       "RoIs");
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(1+CPMHitsMismatch,   "CPMHits");
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(1+CMMHitsMismatch,   "CMMHits");
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(1+LocalSumMismatch,  "Local");
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(1+RemoteSumMismatch, "Remote");
  m_h_CPneSIMSummary->GetXaxis()->SetBinLabel(1+TotalSumMismatch,  "Total");
  m_h_CPneSIMSummary->GetXaxis()->SetLabelSize(0.045);

  m_events = 0;

  } // end if (isNewRun ...

  m_log << MSG::DEBUG << "Leaving bookHistograms" << endreq;
  
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode CPMSimBSMon::fillHistograms()
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "fillHistograms entered" << endreq;

  //Retrieve Trigger Towers from SG
  const TriggerTowerCollection* triggerTowerTES = 0; 
  StatusCode sc = m_storeGate->retrieve(triggerTowerTES,
                                                     m_triggerTowerLocation); 
  if( sc.isFailure()  ||  !triggerTowerTES ) {
    m_log << MSG::DEBUG<< "No Trigger Tower container found"<< endreq; 
  }

  //Retrieve Core CPM Towers from SG
  const CpmTowerCollection* cpmTowerTES = 0; 
  sc = m_storeGate->retrieve(cpmTowerTES, m_cpmTowerLocation); 
  if( sc.isFailure()  ||  !cpmTowerTES ) {
    m_log << MSG::DEBUG<< "No Core CPM Tower container found"<< endreq; 
  }

  //Retrieve Overlap CPM Towers from SG
  const CpmTowerCollection* cpmTowerOvTES = 0; 
  sc = m_storeGate->retrieve(cpmTowerOvTES, m_cpmTowerLocationOverlap); 
  if( sc.isFailure()  ||  !cpmTowerOvTES ) {
    m_log << MSG::DEBUG<< "No Overlap CPM Tower container found"<< endreq; 
  }
  
  //Retrieve CPM RoIs from SG
  const CpmRoiCollection* cpmRoiTES = 0;
  sc = m_storeGate->retrieve( cpmRoiTES, m_cpmRoiLocation);
  if( sc.isFailure()  ||  !cpmRoiTES  ||  cpmRoiTES->empty() ) {
    m_log << MSG::DEBUG << "No DAQ CPM RoIs found, trying RoIB"
          << endreq; 
    cpmRoiTES = 0;
    sc = m_storeGate->retrieve( cpmRoiTES, m_cpmRoiLocationRoib);
    if( sc.isFailure()  ||  !cpmRoiTES ) {
      m_log << MSG::DEBUG << "No RoIB CPM RoIs container found"<< endreq;
    }
  }
  
  //Retrieve CPM Hits from SG
  const CpmHitsCollection* cpmHitsTES = 0;
  sc = m_storeGate->retrieve( cpmHitsTES, m_cpmHitsLocation);
  if( sc.isFailure()  ||  !cpmHitsTES ) {
    m_log << MSG::DEBUG << "No CPM Hits container found"<< endreq; 
  }
  
  //Retrieve CMM-CP Hits from SG
  const CmmCpHitsCollection* cmmCpHitsTES = 0;
  sc = m_storeGate->retrieve( cmmCpHitsTES, m_cmmCpHitsLocation);
  if( sc.isFailure()  ||  !cmmCpHitsTES ) {
    m_log << MSG::DEBUG << "No CMM-CP Hits container found"<< endreq; 
  }

  // Maps to simplify comparisons
  
  TriggerTowerMap ttMap;
  CpmTowerMap     cpMap;
  CpmTowerMap     ovMap;
  CpmRoiMap       crMap;
  CpmHitsMap      chMap;
  CmmCpHitsMap    cmMap;
  setupMap(triggerTowerTES, ttMap);
  setupMap(cpmTowerTES, cpMap);
  setupMap(cpmTowerOvTES, ovMap);
  if (m_compareWithSim) setupMap(cpmRoiTES, crMap);
  setupMap(cpmHitsTES, chMap);
  setupMap(cmmCpHitsTES, cmMap);

  // Vectors for error overview bits;
  const int nCrates = 4;
  const int nCPMs   = 14;
  const int nCMMs   = 2;
  const int vecsizeCpm = 2 * nCrates * nCPMs;
  const int vecsizeCmm = 2 * nCrates * nCMMs;
  ErrorVector errorsCPM(vecsizeCpm);
  ErrorVector errorsCMM(vecsizeCmm);

  // Note - Simulation steps which are simply copies of data from
  // one container to another are not actually done to save time
  // and space.  The comparisons are made with the input instead.

  // Compare Trigger Towers and CPM Towers from data

  bool overlap = false;
  compare(ttMap, cpMap, errorsCPM, overlap);
  overlap = true;
  compare(ttMap, ovMap, errorsCPM, overlap);

  if (m_compareWithSim) {

  // Compare RoIs simulated from CPM Towers with CPM RoIs from data

  CpmRoiCollection* cpmRoiSIM = 0;
  if (cpmTowerTES || cpmTowerOvTES) {
    cpmRoiSIM = new CpmRoiCollection;
    simulate(cpMap, ovMap, cpmRoiSIM);
  }
  CpmRoiMap crSimMap;
  setupMap(cpmRoiSIM, crSimMap);
  compare(crSimMap, crMap, errorsCPM);
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
  compare(chSimMap, chMap, errorsCPM);
  chSimMap.clear();
  delete cpmHitsSIM;

  } // end if (m_compareWithSim)

  // Compare CPM hits with CMM Hits from data

  compare(chMap, cmMap, errorsCPM, errorsCMM);

  if (m_compareWithSim) {

  // Compare Local sums simulated from CMM Hits with Local sums from data

  CmmCpHitsCollection* cmmLocalSIM = 0;
  if (cmmCpHitsTES) {
    cmmLocalSIM = new CmmCpHitsCollection;
    simulate(cmmCpHitsTES, cmmLocalSIM, LVL1::CMMCPHits::LOCAL);
  }
  CmmCpHitsMap cmmLocalSimMap;
  setupMap(cmmLocalSIM, cmmLocalSimMap);
  compare(cmmLocalSimMap, cmMap, errorsCMM, LVL1::CMMCPHits::LOCAL);
  cmmLocalSimMap.clear();
  delete cmmLocalSIM;

  } // end if (m_compareWithSim)

  // Compare Local sums with Remote sums from data

  compare(cmMap, cmMap, errorsCMM, LVL1::CMMCPHits::REMOTE_0);

  if (m_compareWithSim) {

  // Compare Total sums simulated from Remote sums with Total sums from data

  CmmCpHitsCollection* cmmTotalSIM = 0;
  if (cmmCpHitsTES) {
    cmmTotalSIM = new CmmCpHitsCollection;
    simulate(cmmCpHitsTES, cmmTotalSIM, LVL1::CMMCPHits::TOTAL);
  }
  CmmCpHitsMap cmmTotalSimMap;
  setupMap(cmmTotalSIM, cmmTotalSimMap);
  compare(cmmTotalSimMap, cmMap, errorsCMM, LVL1::CMMCPHits::TOTAL);
  cmmTotalSimMap.clear();
  delete cmmTotalSIM;

  } // end if (m_compareWithSim)

  // Update error summary plots

  const int cpmBins = nCrates * nCPMs;
  const int cmmBins = nCrates * nCMMs;
  for (int err = 0; err < NumberOfSummaryBins; ++err) {
    int error = 0;
    for (int loc = 0; loc < cpmBins; ++loc) {
      if ((errorsCPM[loc] >> err) & 0x1) {
        m_h_CPeqSIM->Fill(loc, err, 1.);
      }
      if ((errorsCPM[loc + cpmBins] >> err) & 0x1) {
        m_h_CPneSIM->Fill(loc, err, 1.);
	error = 1;
      }
      if (loc < cmmBins) {
        if ((errorsCMM[loc] >> err) & 0x1) {
          m_h_CPeqSIM->Fill(loc+cpmBins, err, 1.);
        }
        if ((errorsCMM[loc + cmmBins] >> err) & 0x1) {
          m_h_CPneSIM->Fill(loc+cpmBins, err, 1.);
	  error = 1;
        }
      }
    }
    m_h_CPneSIMSummary->Fill(err, error);
  }
  ++m_events;
  std::ostringstream cnum;
  cnum << m_events;
  std::string title("CP Comparison with Simulation Mismatch Summary for "
                    + cnum.str() + " Events");
  m_h_CPneSIMSummary->SetTitle(TString(title));

  m_log << MSG::DEBUG << "Leaving fillHistograms" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode CPMSimBSMon::procHistograms(bool isEndOfEventsBlock,
                                  bool isEndOfLumiBlock, bool isEndOfRun)
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "procHistograms entered" << endreq;

  if (isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun) {
  }

  return StatusCode::SUCCESS;
}

TH1F* CPMSimBSMon::book1F(std::string name, std::string title,
                                   int nx, double xmin, double xmax)
{
  TH1F *hist = new TH1F(TString(name), TString(title), nx, xmin, xmax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    m_log << MSG::WARNING << "Could not register histogram : " 
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
    m_log << MSG::WARNING << "Could not register histogram : " 
	  << name << endreq;
  }
  hist->SetOption("colz");
  hist->SetStats(kFALSE);
  
  return hist;
}

//  Compare Trigger Towers and CPM Towers

void CPMSimBSMon::compare(const TriggerTowerMap& ttMap,
                          const CpmTowerMap& cpMap, ErrorVector& errors,
			  bool overlap)
{
  m_log << MSG::DEBUG << "Compare Trigger Towers and CPM Towers" << endreq;

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
    const int crate = (overlap) ? converter.cpCrateOverlap(coord)
                                : converter.cpCrate(coord);
    const int cpm   = (overlap) ? converter.cpModuleOverlap(coord)
                                : converter.cpModule(coord);
    if (crate > 3 || cpm > 14) continue;
    const int loc = crate * 14 + cpm - 1;
    const int cpmBins = 4 * 14;
    const int bitEm  = (1 << EMTowerMismatch);
    const int bitHad = (1 << HadTowerMismatch);
    if (ttEm && ttEm == cpEm)    errors[loc] |= bitEm;
    if (ttHad && ttHad == cpHad) errors[loc] |= bitHad;
    if (ttEm != cpEm)   errors[loc+cpmBins]  |= bitEm;
    if (ttHad != cpHad) errors[loc+cpmBins]  |= bitHad;
    const double phiMod = phi * m_phiScale;
    if (overlap) {
      m_h_EMTowerOvSIMeqDAT->Fill(eta, phiMod, ttEm && ttEm == cpEm);
      m_h_EMTowerOvSIMneDAT->Fill(eta, phiMod, ttEm && cpEm && ttEm != cpEm);
      m_h_EMTowerOvSIMnoDAT->Fill(eta, phiMod, ttEm && !cpEm);
      m_h_EMTowerOvDATnoSIM->Fill(eta, phiMod, cpEm && !ttEm);
      m_h_HadTowerOvSIMeqDAT->Fill(eta, phiMod, ttHad && ttHad == cpHad);
      m_h_HadTowerOvSIMneDAT->Fill(eta, phiMod, ttHad && cpHad
                                                      && ttHad != cpHad);
      m_h_HadTowerOvSIMnoDAT->Fill(eta, phiMod, ttHad && !cpHad);
      m_h_HadTowerOvDATnoSIM->Fill(eta, phiMod, cpHad && !ttHad);
    } else {
      m_h_EMTowerSIMeqDAT->Fill(eta, phiMod, ttEm && ttEm == cpEm);
      m_h_EMTowerSIMneDAT->Fill(eta, phiMod, ttEm && cpEm && ttEm != cpEm);
      m_h_EMTowerSIMnoDAT->Fill(eta, phiMod, ttEm && !cpEm);
      m_h_EMTowerDATnoSIM->Fill(eta, phiMod, cpEm && !ttEm);
      m_h_HadTowerSIMeqDAT->Fill(eta, phiMod, ttHad && ttHad == cpHad);
      m_h_HadTowerSIMneDAT->Fill(eta, phiMod, ttHad && cpHad && ttHad != cpHad);
      m_h_HadTowerSIMnoDAT->Fill(eta, phiMod, ttHad && !cpHad);
      m_h_HadTowerDATnoSIM->Fill(eta, phiMod, cpHad && !ttHad);
    }
    if (overlap) {
      const double twoPi    = 2.*M_PI;
      const double piByFour = M_PI/4.;
      if (phi > 7.*piByFour)   phi -= twoPi;
      else if (phi < piByFour) phi += twoPi;
    }
    const int loc2 = fpga(crate, phi);
    m_h_FpgaTowerSIMeqDAT->Fill(loc, loc2,   ttEm && ttEm == cpEm);
    m_h_FpgaTowerSIMeqDAT->Fill(loc, loc2+1, ttHad && ttHad == cpHad);
    m_h_FpgaTowerSIMneDAT->Fill(loc, loc2,   ttEm && cpEm && ttEm != cpEm);
    m_h_FpgaTowerSIMneDAT->Fill(loc, loc2+1, ttHad && cpHad && ttHad != cpHad);
    m_h_FpgaTowerSIMnoDAT->Fill(loc, loc2,   ttEm && !cpEm);
    m_h_FpgaTowerSIMnoDAT->Fill(loc, loc2+1, ttHad && !cpHad);
    m_h_FpgaTowerDATnoSIM->Fill(loc, loc2,   cpEm && !ttEm);
    m_h_FpgaTowerDATnoSIM->Fill(loc, loc2+1, cpHad && !ttHad);
  }
}

//  Compare Simulated RoIs with data

void CPMSimBSMon::compare(const CpmRoiMap& roiSimMap, const CpmRoiMap& roiMap,
                                                      ErrorVector& errors)
{
  m_log << MSG::DEBUG << "Compare Simulated RoIs with data" << endreq;

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
    const int bit = (1 << RoIMismatch);
    if (simHits && simHits == datHits) errors[locX] |= bit;
    if (simHits != datHits) errors[locX+cpmBins]    |= bit;

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

    m_log << MSG::DEBUG << "DataHits/SimHits: ";
    for (int i = 15; i >= 0; --i) {
      int hit = (datHits >> i) & 0x1;
      m_log << MSG::DEBUG << hit;
    }
    m_log << MSG::DEBUG << "/";
    for (int i = 15; i >= 0; --i) {
      int hit = (simHits >> i) & 0x1;
      m_log << MSG::DEBUG << hit;
    }
    m_log << MSG::DEBUG << endreq;
  }
}

//  Compare simulated CPM Hits with data

void CPMSimBSMon::compare(const CpmHitsMap& cpmSimMap, const CpmHitsMap& cpmMap,
                                                       ErrorVector& errors)
{
  m_log << MSG::DEBUG << "Compare simulated CPM Hits with data" << endreq;

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
    const int bit = (1 << CPMHitsMismatch);
    if ((simHits0 && simHits0 == datHits0) ||
        (simHits1 && simHits1 == datHits1)) errors[loc] |= bit;
    if (simHits0 != datHits0 || simHits1 != datHits1)
                                            errors[loc+cpmBins] |= bit;

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
                          ErrorVector& errorsCPM, ErrorVector& errorsCMM)
{
  m_log << MSG::DEBUG << "Compare CPM Hits and CMM Hits" << endreq;

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
    const int bit = (1 << CMMHitsMismatch);
    if ((cpmHits0 && cpmHits0 == cmmHits0) ||
        (cpmHits1 && cpmHits1 == cmmHits1)) errorsCPM[loc] |= bit;
    if (cpmHits0 != cmmHits0 || cpmHits1 != cmmHits1)
                                            errorsCPM[loc+cpmBins] |= bit;
    if (cpmHits1 && cpmHits1 == cmmHits1) errorsCMM[loc2]   |= bit;
    if (cpmHits0 && cpmHits0 == cmmHits0) errorsCMM[loc2+1] |= bit;
    if (cpmHits1 != cmmHits1) errorsCMM[loc2+cmmBins] |= bit; // hits1==>cmm 0
    if (cpmHits0 != cmmHits0)
                            errorsCMM[loc2+cmmBins+1] |= bit; // hits0==>cmm 1

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
  m_log << MSG::DEBUG << "Compare Simulated CMM Hit Sums and Data CMM Hit Sums"
        << endreq;

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
      const int bit = (local) ? (1 << LocalSumMismatch)
                              : (1 << TotalSumMismatch);
      if (cmmSimHits1 && cmmSimHits1 == cmmHits1) errors[loc]   |= bit;
      if (cmmSimHits0 && cmmSimHits0 == cmmHits0) errors[loc+1] |= bit;
      if (cmmSimHits1 != cmmHits1) errors[loc+cmmBins]   |= bit; // cmm 0
      if (cmmSimHits0 != cmmHits0) errors[loc+cmmBins+1] |= bit; // cmm 1
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
      const int bit = (1 << RemoteSumMismatch);
      if (hits1Sim[crate] && hits1Sim[crate] == hits1[crate])
                                                         errors[loc] |= bit;
      if (hits0Sim[crate] && hits0Sim[crate] == hits0[crate])
                                                         errors[loc+1] |= bit;
      if (hits1Sim[crate] != hits1[crate]) errors[loc+cmmBins]   |= bit;
      if (hits0Sim[crate] != hits0[crate]) errors[loc+cmmBins+1] |= bit;
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
  // Simulation steps in red (#color[2]) depend on Trigger Menu
  hist->GetYaxis()->SetBinLabel(1+EMTowerMismatch,  "EM tt");
  hist->GetYaxis()->SetBinLabel(1+HadTowerMismatch, "Had tt");
  if (m_compareWithSim) {
    hist->GetYaxis()->SetBinLabel(1+RoIMismatch,     "#color[2]{RoIs}");
    hist->GetYaxis()->SetBinLabel(1+CPMHitsMismatch, "CPMHits");
  } else {
    hist->GetYaxis()->SetBinLabel(1+RoIMismatch,     "n/a");
    hist->GetYaxis()->SetBinLabel(1+CPMHitsMismatch, "n/a");
  }
  hist->GetYaxis()->SetBinLabel(1+CMMHitsMismatch, "CMMHits");
  if (m_compareWithSim) {
    hist->GetYaxis()->SetBinLabel(1+LocalSumMismatch, "Local");
  } else {
    hist->GetYaxis()->SetBinLabel(1+LocalSumMismatch, "n/a");
  }
  hist->GetYaxis()->SetBinLabel(1+RemoteSumMismatch, "Remote");
  if (m_compareWithSim) {
    hist->GetYaxis()->SetBinLabel(1+TotalSumMismatch, "Total");
  } else {
    hist->GetYaxis()->SetBinLabel(1+TotalSumMismatch, "n/a");
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
  setLabelsYNUM(hist, 0, 15);
}

void CPMSimBSMon::setLabelsCPMFP(TH2* hist)
{
  setLabelsCPM(hist);
  setLabelsYNUM(hist, 0, 19);
}

void CPMSimBSMon::setLabelsYNUM(TH2* hist, int beg, int end)
{
  int bin = 1;
  for (int val = beg; val <= end; ++val) {
    std::ostringstream cnum;
    cnum << val;
    hist->GetYaxis()->SetBinLabel(bin++, cnum.str().c_str());
  }
  hist->GetYaxis()->SetLabelSize(0.05);
}

void CPMSimBSMon::setLabelsXNUM(TH2* hist, int beg, int end)
{
  int bin = 1;
  for (int val = beg; val <= end; ++val) {
    std::ostringstream cnum;
    cnum << val;
    hist->GetXaxis()->SetBinLabel(bin++, cnum.str().c_str());
  }
  hist->GetXaxis()->SetLabelSize(0.05);
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
  setLabelsXNUM(hist, 1, 14);
  setLabelsYNUM(hist, 0, 3);
}

void CPMSimBSMon::setLabelsMCLR(TH2* hist)
{
  setLabelsXNUM(hist, 1, 14);
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
  setLabelsYNUM(hist, 0, 15);
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
  setLabelsYNUM(hist, 0, 15);
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

void CPMSimBSMon::simulate(const CpmTowerMap towers, const CpmTowerMap towersOv,
                                 CpmRoiCollection* rois)
{
  m_log << MSG::DEBUG << "Simulate CPM RoIs from CPM Towers" << endreq;

  // Process a crate at a time to use overlap data
  const int ncrates = 4;
  std::vector<CpmTowerMap> crateMaps(ncrates);
  LVL1::CoordToHardware converter;
  CpmTowerMap::const_iterator iter  = towers.begin();
  CpmTowerMap::const_iterator iterE = towers.end();
  for (; iter != iterE; ++iter) {
    LVL1::CPMTower* tt = iter->second;
    const LVL1::Coordinate coord(tt->phi(), tt->eta());
    const int crate = converter.cpCrate(coord);
    if (crate >= ncrates) continue;
    crateMaps[crate].insert(std::make_pair(iter->first, tt));
  }
  iter  = towersOv.begin();
  iterE = towersOv.end();
  for (; iter != iterE; ++iter) {
    LVL1::CPMTower* tt = iter->second;
    const LVL1::Coordinate coord(tt->phi(), tt->eta());
    const int crate = converter.cpCrateOverlap(coord);
    if (crate >= ncrates) continue;
    crateMaps[crate].insert(std::make_pair(iter->first, tt));
  }
  for (int crate = 0; crate < ncrates; ++crate) {
    InternalRoiCollection* intRois = new InternalRoiCollection;
    m_emTauTool->findRoIs(&crateMaps[crate], intRois);
    InternalRoiCollection::iterator roiIter  = intRois->begin();
    InternalRoiCollection::iterator roiIterE = intRois->end();
    for (; roiIter != roiIterE; ++roiIter) {
      LVL1::CPMRoI* roi = new LVL1::CPMRoI((*roiIter)->RoIWord());
      if (roi->crate() == crate) rois->push_back(roi);
      else delete roi;
    }
    delete intRois;
  }
}

void CPMSimBSMon::simulate(const CpmRoiCollection* rois,
                                 CpmHitsCollection* hits)
{
  m_log << MSG::DEBUG << "Simulate CPM Hits from CPM RoIs" << endreq;

  m_cpHitsTool->formCPMHits(rois, hits);
}

void CPMSimBSMon::simulate(const CmmCpHitsCollection* hitsIn,
                                 CmmCpHitsCollection* hitsOut, int selection)
{
  m_log << MSG::DEBUG << "Simulate CMM Hit sums from CMM Hits" << endreq;

  if (selection == LVL1::CMMCPHits::LOCAL) {
    m_cpHitsTool->formCMMCPHitsCrate(hitsIn, hitsOut);
  } else if (selection == LVL1::CMMCPHits::TOTAL) {
    m_cpHitsTool->formCMMCPHitsSystem(hitsIn, hitsOut);
  }
}

// Return EM FPGA for given crate/phi

int CPMSimBSMon::fpga(int crate, double phi)
{
  const double phiGran = M_PI/32.;
  const double phiBase = M_PI/2. * double(crate);
  const int phiBin = int(floor((phi - phiBase) / phiGran)) + 2;
  return 2 * (phiBin/2);
}
