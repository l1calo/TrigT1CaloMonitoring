// ********************************************************************
//
// NAME:     CPMSimBSMon.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Peter Faulkner
//           
//
// ********************************************************************

#include <iomanip>
#include <sstream>
#include <utility>
#include <cmath>

#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "GaudiKernel/ITHistSvc.h"
#include "StoreGate/StoreGateSvc.h"
#include "SGTools/StlVectorClids.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigT1CaloEvent/CMMCPHits.h"
#include "TrigT1CaloEvent/CPMHits.h"
#include "TrigT1CaloEvent/CPMTower.h"
#include "TrigT1CaloEvent/CPMRoI.h"
#include "TrigT1CaloEvent/RODHeader.h"
#include "TrigT1CaloEvent/TriggerTower.h"
#include "TrigT1CaloUtils/CoordToHardware.h"
#include "TrigT1CaloUtils/CPAlgorithm.h"
#include "TrigT1CaloUtils/DataError.h"
#include "TrigT1CaloUtils/TriggerTowerKey.h"
#include "TrigT1CaloToolInterfaces/IL1EmTauTools.h"
#include "TrigT1CaloToolInterfaces/IL1CPHitsTools.h"
#include "TrigT1Interfaces/Coordinate.h"
#include "TrigT1Interfaces/CoordinateRange.h"
#include "TrigT1Interfaces/CPRoIDecoder.h"
#include "TrigT1Interfaces/TrigT1CaloDefs.h"

#include "TrigT1CaloMonitoring/CPMSimBSMon.h"
#include "TrigT1CaloMonitoring/TrigT1CaloMonErrorTool.h"
#include "TrigT1CaloMonitoring/TrigT1CaloHistogramTool.h"

/*---------------------------------------------------------*/
CPMSimBSMon::CPMSimBSMon(const std::string & type, 
			 const std::string & name,
			 const IInterface* parent)
  : ManagedMonitorToolBase(type, name, parent),
    m_storeGate("StoreGateSvc", name),
    m_emTauTool("LVL1::L1EmTauTools/L1EmTauTools"),
    m_cpHitsTool("LVL1::L1CPHitsTools/L1CPHitsTools"),
    m_errorTool("TrigT1CaloMonErrorTool"),
    m_histTool("TrigT1CaloHistogramTool"),
    m_log(msgSvc(), name), m_debug(false), m_events(0)
/*---------------------------------------------------------*/
{
  declareInterface<IMonitorToolBase>(this); 

  declareProperty("EmTauTool", m_emTauTool);
  declareProperty("CPHitsTool", m_cpHitsTool);

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
  declareProperty("TriggerTowerLocation",
                 m_triggerTowerLocation =
		                 LVL1::TrigT1CaloDefs::TriggerTowerLocation);
  declareProperty("RodHeaderLocation",
                 m_rodHeaderLocation = "RODHeaders");

  declareProperty("RootDirectory", m_rootDir = "L1Calo");

  declareProperty("CompareWithSimulation", m_compareWithSim = true,
                  "Include the checks that run simulation tools");
  declareProperty("CompareTriggerTowers", m_compareTriggerTowers = true,
                  "Include TriggerTower/CPMTower comparison");
  declareProperty("RoIThresholds", m_roiThresh,
                  "RoI thresholds to compare, default 0-15");
  declareProperty("IgnoreTowersEM", m_ignoreTowersEm,
                  "EM TriggerTowers with known problems");
  declareProperty("IgnoreTowersHad", m_ignoreTowersHad,
                  "Had TriggerTowers with known problems");
}

/*---------------------------------------------------------*/
CPMSimBSMon::~CPMSimBSMon()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode CPMSimBSMon::initialize()
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

  sc = m_errorTool.retrieve();
  if( sc.isFailure() ) {
    m_log << MSG::ERROR << "Unable to locate Tool TrigT1CaloMonErrorTool"
                        << endreq;
    return sc;
  }

  sc = m_histTool.retrieve();
  if( sc.isFailure() ) {
    m_log << MSG::ERROR << "Unable to locate Tool TrigT1CaloHistogramTool"
                        << endreq;
    return sc;
  }

  // RoI thresholds mask
  if (m_roiThresh.empty()) m_roiMask = 0xffff;
  else {
    m_roiMask = 0;
    std::vector<int>::const_iterator pos  = m_roiThresh.begin();
    std::vector<int>::const_iterator posE = m_roiThresh.end();
    for (; pos != posE; ++pos) {
      const int thresh = *pos;
      if (thresh >= 0 && thresh < 16) m_roiMask |= (1 << thresh);
    }
    m_log << MSG::INFO << "RoI comparison mask: " << MSG::hex
          << m_roiMask << MSG::dec << endreq;
  }

  m_log << MSG::INFO << "CPMSimBSMon initialised" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode CPMSimBSMon::finalize()
/*---------------------------------------------------------*/
{
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

  std::string dir1(m_rootDir + "/CPM/Errors/Transmission_Simulation");
  MonGroup monShift( this, dir1, shift, run );
  MonGroup monExpert( this, dir1, expert, run );
  MonGroup monCPMin( this, dir1 + "/PPM2CPMTowers", expert, run );
  MonGroup monRoIs( this, dir1 + "/Towers2RoIs", expert, run );
  MonGroup monCPMout( this, dir1 + "/RoIs2Hits", expert, run);
  MonGroup monEvent1( this, dir1 + "/MismatchEventNumbers", expert, run, "",
                                                            "eventSample" );
  std::string dir2(m_rootDir + "/CPM_CMM/Errors/Transmission_Simulation");
  MonGroup monCMMin( this, dir2 + "/CPM2CMMHits", expert, run );
  MonGroup monCMMout( this, dir2 + "/Hits2Sums", expert, run );
  MonGroup monEvent2( this, dir2 + "/MismatchEventNumbers", expert, run, "",
                                                            "eventSample" );

  // CPMTowers

  m_histTool->setMonGroup(&monCPMin);

  if (m_compareTriggerTowers) {

  m_h_EMTowerSIMeqDAT = m_histTool->bookCPMEtaVsPhi(
    "cpm_em_2d_etaPhi_tt_PpmEqCore",
    "EM Core CPM Tower/PPM Tower Non-zero Matches");
  m_h_EMTowerSIMneDAT = m_histTool->bookCPMEtaVsPhi(
    "cpm_em_2d_etaPhi_tt_PpmNeCore",
    "EM Core CPM Tower/PPM Tower Non-zero Mismatches");
  m_h_EMTowerSIMnoDAT = m_histTool->bookCPMEtaVsPhi(
    "cpm_em_2d_etaPhi_tt_PpmNoCore",
    "EM PPM Towers but no Core CPM Towers");
  m_h_EMTowerDATnoSIM = m_histTool->bookCPMEtaVsPhi(
    "cpm_em_2d_etaPhi_tt_CoreNoPpm",
    "EM Core CPM Towers but no PPM Towers");
  m_h_HadTowerSIMeqDAT = m_histTool->bookCPMEtaVsPhi(
    "cpm_had_2d_etaPhi_tt_PpmEqCore",
    "HAD Core CPM Tower/PPM Tower Non-zero Matches");
  m_h_HadTowerSIMneDAT = m_histTool->bookCPMEtaVsPhi(
    "cpm_had_2d_etaPhi_tt_PpmNeCore",
    "HAD Core CPM Tower/PPM Tower Non-zero Mismatches");
  m_h_HadTowerSIMnoDAT = m_histTool->bookCPMEtaVsPhi(
    "cpm_had_2d_etaPhi_tt_PpmNoCore",
    "HAD PPM Towers but no Core CPM Towers");
  m_h_HadTowerDATnoSIM = m_histTool->bookCPMEtaVsPhi(
    "cpm_had_2d_etaPhi_tt_CoreNoPpm",
    "HAD Core CPM Towers but no PPM Towers");
  m_h_EMTowerOvSIMeqDAT = m_histTool->bookCPMEtaVsPhi(
    "cpm_em_2d_etaPhi_tt_PpmEqOverlap",
    "EM Overlap CPM Tower/PPM Tower Non-zero Matches");
  m_h_EMTowerOvSIMneDAT = m_histTool->bookCPMEtaVsPhi(
    "cpm_em_2d_etaPhi_tt_PpmNeOverlap",
    "EM Overlap CPM Tower/PPM Tower Non-zero Mismatches");
  m_h_EMTowerOvSIMnoDAT = m_histTool->bookCPMEtaVsPhi(
    "cpm_em_2d_etaPhi_tt_PpmNoOverlap",
    "EM PPM Towers but no Overlap CPM Towers");
  m_h_EMTowerOvDATnoSIM = m_histTool->bookCPMEtaVsPhi(
    "cpm_em_2d_etaPhi_tt_OverlapNoPpm",
    "EM Overlap CPM Towers but no PPM Towers");
  m_h_HadTowerOvSIMeqDAT = m_histTool->bookCPMEtaVsPhi(
    "cpm_had_2d_etaPhi_tt_PpmEqOverlap",
    "HAD Overlap CPM Tower/PPM Tower Non-zero Matches");
  m_h_HadTowerOvSIMneDAT = m_histTool->bookCPMEtaVsPhi(
    "cpm_had_2d_etaPhi_tt_PpmNeOverlap",
    "HAD Overlap CPM Tower/PPM Tower Non-zero Mismatches");
  m_h_HadTowerOvSIMnoDAT = m_histTool->bookCPMEtaVsPhi(
    "cpm_had_2d_etaPhi_tt_PpmNoOverlap",
    "HAD PPM Towers but no Overlap CPM Towers");
  m_h_HadTowerOvDATnoSIM = m_histTool->bookCPMEtaVsPhi(
    "cpm_had_2d_etaPhi_tt_OverlapNoPpm",
    "HAD Overlap CPM Towers but no PPM Towers");

  m_h_FpgaTowerSIMeqDAT = m_histTool->bookCPMCrateModuleVsFPGA(
    "cpm_2d_tt_PpmEqCpmFpga", "CPM Tower/PPM Tower Non-zero Matches by FPGA");
  m_h_FpgaTowerSIMneDAT = m_histTool->bookCPMCrateModuleVsFPGA(
   "cpm_2d_tt_PpmNeCpmFpga", "CPM Tower/PPM Tower Non-zero Mismatches by FPGA");
  m_h_FpgaTowerSIMnoDAT = m_histTool->bookCPMCrateModuleVsFPGA(
    "cpm_2d_tt_PpmNoCpmFpga", "PPM Towers but no CPM Towers by FPGA");
  m_h_FpgaTowerDATnoSIM = m_histTool->bookCPMCrateModuleVsFPGA(
    "cpm_2d_tt_CpmNoPpmFpga", "CPM Towers but no PPM Towers by FPGA");

  if ( !m_ignoreTowersEm.empty() || !m_ignoreTowersHad.empty() ) {
    m_h_IgnoreTowersEM = m_histTool->bookCPMEtaVsPhi(
      "cpm_em_2d_etaPhi_tt_IgnoredTowers",
      "EM Tower Mismatches in Ignore List");
    m_h_IgnoreTowersHad = m_histTool->bookCPMEtaVsPhi(
      "cpm_had_2d_etaPhi_tt_IgnoredTowers",
      "Had Tower Mismatches in Ignore List");
  }

  } // end if (m_compareTriggerTowers)

  if (m_compareWithSim) {

  // RoIs

  m_histTool->setMonGroup(&monRoIs);

  m_h_RoISIMeqDAT = m_histTool->bookCPMCrateModuleVsChipLocalCoord(
    "cpm_2d_roi_SimEqData", "CPM RoI Data/Simulation Non-zero Matches");
  m_h_RoISIMneDAT = m_histTool->bookCPMCrateModuleVsChipLocalCoord(
    "cpm_2d_roi_SimNeData", "CPM RoI Data/Simulation Non-zero Mismatches");
  m_h_RoISIMnoDAT = m_histTool->bookCPMCrateModuleVsChipLocalCoord(
    "cpm_2d_roi_SimNoData", "CPM RoI Simulation but no Data");
  m_h_RoIDATnoSIM = m_histTool->bookCPMCrateModuleVsChipLocalCoord(
    "cpm_2d_roi_DataNoSim", "CPM RoI Data but no Simulation");

  m_h_RoIThreshSIMeqDAT = m_histTool->bookCPMCrateModuleVsThreshold(
    "cpm_2d_roi_ThreshSimEqData", "CPM RoI Data/Simulation Threshold Matches");
  m_h_RoIThreshSIMneDAT = m_histTool->bookCPMCrateModuleVsThreshold(
    "cpm_2d_roi_ThreshSimNeData",
    "CPM RoI Data/Simulation Threshold Mismatches");

  m_h_RoIEtaPhiSIMeqDAT = m_histTool->bookCPMRoIEtaVsPhi(
    "cpm_2d_etaPhi_roi_SimEqData", "CPM RoI Data/Simulation Non-zero Matches");
  m_h_RoIEtaPhiSIMneDAT = m_histTool->bookCPMRoIEtaVsPhi(
    "cpm_2d_etaPhi_roi_SimNeData",
    "CPM RoI Data/Simulation Non-zero Mismatches");
  m_h_RoIEtaPhiSIMnoDAT = m_histTool->bookCPMRoIEtaVsPhi(
    "cpm_2d_etaPhi_roi_SimNoData", "CPM RoI Simulation but no Data");
  m_h_RoIEtaPhiDATnoSIM = m_histTool->bookCPMRoIEtaVsPhi(
    "cpm_2d_etaPhi_roi_DataNoSim", "CPM RoI Data but no Simulation");

  // CPMHits

  m_histTool->setMonGroup(&monCPMout);

  m_h_CPMHitsSIMeqDAT = m_histTool->bookCPMModuleVsCrate(
    "cpm_2d_thresh_SimEqData", "CPM Hits Data/Simulation Non-zero Matches");
  m_h_CPMHitsSIMneDAT = m_histTool->bookCPMModuleVsCrate(
    "cpm_2d_thresh_SimNeData", "CPM Hits Data/Simulation Non-zero Mismatches");
  m_h_CPMHitsSIMnoDAT = m_histTool->bookCPMModuleVsCrate(
    "cpm_2d_thresh_SimNoData", "CPM Hits Simulation but no Data");
  m_h_CPMHitsDATnoSIM = m_histTool->bookCPMModuleVsCrate(
    "cpm_2d_thresh_DataNoSim", "CPM Hits Data but no Simulation");

  m_h_CPMHitsThreshSIMeqDAT = m_histTool->bookCPMCrateModuleVsThreshold(
    "cpm_2d_thresh_ThreshSimEqData",
    "CPM Hits Data/Simulation Threshold Matches");
  m_h_CPMHitsThreshSIMneDAT = m_histTool->bookCPMCrateModuleVsThreshold(
    "cpm_2d_thresh_ThreshSimNeData",
    "CPM Hits Data/Simulation Threshold Mismatches");

  } // end if (m_compareWithSim)

  // CMMHits

  m_histTool->setMonGroup(&monCMMin);

  m_h_CMMHitsSIMeqDAT = m_histTool->bookCPMModuleVsCrateCMM(
    "cmm_2d_thresh_CpmEqCmm", "CMM Hits/CPM Hits Non-zero Matches");
  m_h_CMMHitsSIMneDAT = m_histTool->bookCPMModuleVsCrateCMM(
    "cmm_2d_thresh_CpmNeCmm", "CMM Hits/CPM Hits Non-zero Mismatches");
  m_h_CMMHitsSIMnoDAT = m_histTool->bookCPMModuleVsCrateCMM(
    "cmm_2d_thresh_CpmNoCmm", "CPM Hits but no CMM Hits");
  m_h_CMMHitsDATnoSIM = m_histTool->bookCPMModuleVsCrateCMM(
    "cmm_2d_thresh_CmmNoCpm", "CMM Hits but no CPM Hits");

  m_h_CMMHitsThreshSIMeqDAT = m_histTool->bookCPMCrateModuleVsThreshold(
    "cmm_2d_thresh_ThreshCpmEqCmm", "CMM Hits/CPM Hits Threshold Matches");
  m_h_CMMHitsThreshSIMneDAT = m_histTool->bookCPMCrateModuleVsThreshold(
    "cmm_2d_thresh_ThreshCpmNeCmm", "CMM Hits/CPM Hits Threshold Mismatches");

  m_histTool->setMonGroup(&monCMMout);

  // Local/Remote/Total sums

  m_h_SumsSIMeqDAT = m_histTool->bookCPMSumCMM("cmm_1d_thresh_SumsSimEqData",
    "CMM Hit Sums Data/Simulation Non-zero Matches");
  m_h_SumsSIMneDAT = m_histTool->bookCPMSumCMM("cmm_1d_thresh_SumsSimNeData",
    "CMM Hit Sums Data/Simulation Non-zero Mismatches");
  m_h_SumsSIMnoDAT = m_histTool->bookCPMSumCMM("cmm_1d_thresh_SumsSimNoData",
    "CMM Hit Sums Simulation but no Data");
  m_h_SumsDATnoSIM = m_histTool->bookCPMSumCMM("cmm_1d_thresh_SumsDataNoSim",
    "CMM Hit Sums Data but no Simulation");
  m_h_SumsThreshSIMeqDAT = m_histTool->bookCPMSumVsThreshold(
    "cmm_2d_thresh_SumsThreshSimEqData",
    "CMM Hit Sums Data/Simulation Threshold Matches");
  m_h_SumsThreshSIMneDAT = m_histTool->bookCPMSumVsThreshold(
    "cmm_2d_thresh_SumsThreshSimNeData",
    "CMM Hit Sums Data/Simulation Threshold Mismatches");

  // Summary

  m_histTool->setMonGroup(&monExpert);

  m_h_CPeqSIM = m_histTool->book2F("cpm_2d_SimEqDataOverview",
   "CP Transmission/Comparison with Simulation Overview - Events with Matches",
             64, 0, 64, NumberOfSummaryBins, 0, NumberOfSummaryBins);
  m_histTool->cpmCMMCrateModule(m_h_CPeqSIM);
  setLabels(m_h_CPeqSIM, false);

  m_h_CPneSIM = m_histTool->book2F("cpm_2d_SimNeDataOverview",
   "CP Transmission/Comparison with Simulation Overview - Events with Mismatches",
             64, 0, 64, NumberOfSummaryBins, 0, NumberOfSummaryBins);
  m_histTool->cpmCMMCrateModule(m_h_CPneSIM);
  setLabels(m_h_CPneSIM, false);

  m_histTool->setMonGroup(&monShift);

  m_h_CPneSIMSummary = m_histTool->book1F("cpm_1d_SimNeDataSummary",
   "CP Transmission/Comparison with Simulation Mismatch Summary;;Events",
    NumberOfSummaryBins, 0, NumberOfSummaryBins);
  setLabels(m_h_CPneSIMSummary);

  // Mismatch Event Number Samples

  m_histTool->setMonGroup(&monEvent1);

  TH2I* hist = 0;
  m_sampleHists.clear();
  m_sampleHists.resize(6, hist);
  if (m_compareTriggerTowers) {
    hist = m_histTool->bookCPMEventVsCrateModule("cpm_em_2d_tt_MismatchEvents",
                                     "CPM Towers EM Mismatch Event Numbers");
    m_sampleHists[0] = hist;
    hist = m_histTool->bookCPMEventVsCrateModule("cpm_had_2d_tt_MismatchEvents",
                                     "CPM Towers Had Mismatch Event Numbers");
    m_sampleHists[1] = hist;
  }
  if (m_compareWithSim) {
    hist = m_histTool->bookCPMEventVsCrateModule("cpm_2d_roi_MismatchEvents",
                                     "CPM RoIs Mismatch Event Numbers");
    m_sampleHists[2] = hist;
    hist = m_histTool->bookCPMEventVsCrateModule("cpm_2d_thresh_MismatchEvents",
                                     "CPM Hits Mismatch Event Numbers");
    m_sampleHists[3] = hist;
  }

  m_histTool->setMonGroup(&monEvent2);

  hist = m_histTool->bookCPMEventVsCrateModule("cmm_2d_thresh_MismatchEvents",
                                     "CMM Hits Mismatch Event Numbers");
  m_sampleHists[4] = hist;
  hist = m_histTool->bookEventNumbers("cmm_2d_thresh_SumsMismatchEvents",
                     "CMM Hit Sums Mismatch Event Numbers", 24, 0., 24.);
  m_histTool->cpmCrateCMM(hist, 0, false);
  m_histTool->cpmCrateCMM(hist, 8, false);
  m_histTool->cpmCrateCMM(hist, 16, false);
  hist->GetYaxis()->SetBinLabel(1, "Module 0/0");
  hist->GetYaxis()->SetBinLabel(9, "Local 0/0");
  hist->GetYaxis()->SetBinLabel(17, "Remote 0/0");
  hist->GetYaxis()->SetBinLabel(23, "Total 3/0");
  hist->GetYaxis()->SetTitle("");
  m_sampleHists[5] = hist;

  m_histTool->unsetMonGroup();

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
  
  // Skip events believed to be corrupt

  if (m_errorTool->corrupt()) {
    m_log << MSG::DEBUG << "Skipping corrupt event" << endreq;
    return StatusCode::SUCCESS;
  }

  // NB. Collection retrieves wrapped in m_storeGate->contains<..>(..)
  // are for those not expected to be on ESD. They should be on bytestream.

  StatusCode sc;

  //Retrieve Trigger Towers from SG
  const TriggerTowerCollection* triggerTowerTES = 0; 
  if (m_compareTriggerTowers) {
    sc = m_storeGate->retrieve(triggerTowerTES, m_triggerTowerLocation); 
    if( sc.isFailure()  ||  !triggerTowerTES ) {
      m_log << MSG::DEBUG<< "No Trigger Tower container found"<< endreq; 
    }
  }

  //Retrieve Core and Overlap CPM Towers from SG
  const CpmTowerCollection* cpmTowerTES = 0; 
  const CpmTowerCollection* cpmTowerOvTES = 0; 
  if (m_compareTriggerTowers || m_compareWithSim) {
    sc = m_storeGate->retrieve(cpmTowerTES, m_cpmTowerLocation); 
    if( sc.isFailure()  ||  !cpmTowerTES ) {
      m_log << MSG::DEBUG<< "No Core CPM Tower container found"<< endreq; 
    }
    if (m_storeGate->contains<CpmTowerCollection>(m_cpmTowerLocationOverlap)) {
      sc = m_storeGate->retrieve(cpmTowerOvTES, m_cpmTowerLocationOverlap); 
    } else sc = StatusCode::FAILURE;
    if( sc.isFailure()  ||  !cpmTowerOvTES ) {
      m_log << MSG::DEBUG<< "No Overlap CPM Tower container found"<< endreq; 
    }
  }
  m_overlapPresent = cpmTowerOvTES != 0;
  
  //Retrieve CPM RoIs from SG
  const CpmRoiCollection* cpmRoiTES = 0;
  const RodHeaderCollection* rodTES = 0;
  if (m_compareWithSim) {
    sc = m_storeGate->retrieve( cpmRoiTES, m_cpmRoiLocation);
    if( sc.isFailure()  ||  !cpmRoiTES ) {
      m_log << MSG::DEBUG << "No DAQ CPM RoIs container found"<< endreq;
    } else {
      if (m_storeGate->contains<RodHeaderCollection>(m_rodHeaderLocation)) {
        sc = m_storeGate->retrieve( rodTES, m_rodHeaderLocation);
      } else sc = StatusCode::FAILURE;
      if( sc.isFailure()  ||  !rodTES ) {
        m_log << MSG::DEBUG << "No ROD Header container found"<< endreq;
      }
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
  ++m_events;

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
  setupMap(cpmRoiTES, crMap);
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

  if (m_compareTriggerTowers) {

    // Compare Trigger Towers and CPM Towers from data

    bool overlap = false;
    compare(ttMap, cpMap, errorsCPM, overlap);
    if (m_overlapPresent) {
      overlap = true;
      compare(ttMap, ovMap, errorsCPM, overlap);
    }

  }

  if (m_compareWithSim) {

    // Compare RoIs simulated from CPM Towers with CPM RoIs from data

    CpmRoiCollection* cpmRoiSIM = 0;
    if (cpmTowerTES || cpmTowerOvTES) {
      cpmRoiSIM = new CpmRoiCollection;
      simulate(cpMap, ovMap, cpmRoiSIM);
    }
    CpmRoiMap crSimMap;
    setupMap(cpmRoiSIM, crSimMap);
    compare(crSimMap, crMap, rodTES, errorsCPM);
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

  }

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

  }

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

  }

  // Update error summary plots

  ErrorVector crateErr(nCrates);
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
	crateErr[loc/nCPMs] |= (1 << err);
	if (m_sampleHists[err]) {
	  m_histTool->fillEventNumber(m_sampleHists[err], loc);
        }
      }
      if (loc < cmmBins) {
        if ((errorsCMM[loc] >> err) & 0x1) {
          m_h_CPeqSIM->Fill(loc+cpmBins, err, 1.);
        }
        if ((errorsCMM[loc + cmmBins] >> err) & 0x1) {
          m_h_CPneSIM->Fill(loc+cpmBins, err, 1.);
	  error = 1;
	  crateErr[loc/nCMMs] |= (1 << err);
	  int offset = 0;
	  if (err == LocalSumMismatch) offset = 8;
	  else if (err == RemoteSumMismatch || err == TotalSumMismatch)
	                               offset = 16;
	  if (m_sampleHists[5]) {
	    m_histTool->fillEventNumber(m_sampleHists[5], loc+offset);
	  }
        }
      }
    }
    if (error) m_h_CPneSIMSummary->Fill(err);
  }

  // Save error vector for global summary

  ErrorVector* save = new ErrorVector(crateErr);
  sc = m_storeGate->record(save, "L1CaloCPMMismatchVector");
  if (sc != StatusCode::SUCCESS) {
    m_log << MSG::ERROR << "Error recording CPM mismatch vector in TES "
          << endreq;
    return sc;
  }

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

//  Compare Trigger Towers and CPM Towers

void CPMSimBSMon::compare(const TriggerTowerMap& ttMap,
                          const CpmTowerMap& cpMap, ErrorVector& errors,
			  bool overlap)
{
  m_log << MSG::DEBUG << "Compare Trigger Towers and CPM Towers" << endreq;

  const int nCrates = 4;
  const int nCPMs   = 14;
  const int maxKey = 0x7fffffff;
  LVL1::CoordToHardware converter;
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
    int key = 0;

    if (ttMapIter != ttMapIterEnd) ttKey = ttMapIter->first;
    if (cpMapIter != cpMapIterEnd) cpKey = cpMapIter->first;

    if ((cpMapIter == cpMapIterEnd) || (cpKey > ttKey)) {

      // TriggerTower but no CPMTower

      const LVL1::TriggerTower* tt = ttMapIter->second;
      ++ttMapIter;
      eta = tt->eta();
      phi = tt->phi();
      if (overlap) { // skip non-overlap TTs
        const LVL1::Coordinate coord(phi, eta);
	const int crate = converter.cpCrateOverlap(coord);
        if (crate >= nCrates) continue;
      }
      ttEm  = tt->emEnergy();
      ttHad = tt->hadEnergy();
      key = ttKey;

    } else if ((ttMapIter == ttMapIterEnd) || (ttKey > cpKey)) {

      // CPMTower but no TriggerTower

      const LVL1::CPMTower* cp = cpMapIter->second;
      ++cpMapIter;
      eta = cp->eta();
      phi = cp->phi();
      cpEm  = cp->emEnergy();
      cpHad = cp->hadEnergy();
      key = cpKey;

    } else {

      // Have both

      const LVL1::TriggerTower* tt = ttMapIter->second;
      const LVL1::CPMTower*     cp = cpMapIter->second;
      ++ttMapIter;
      ++cpMapIter;
      eta = tt->eta();
      phi = tt->phi();
      ttEm  = tt->emEnergy();
      ttHad = tt->hadEnergy();
      cpEm  = cp->emEnergy();
      cpHad = cp->hadEnergy();
      key = ttKey;
    }
    
    // Check for known bad towers

    if (ttEm != cpEm && ignoreTower(0, key)) {
      m_histTool->fillCPMEtaVsPhi(m_h_IgnoreTowersEM, eta, phi);
      ttEm = 0;
      cpEm = 0;
    }
    if (ttHad != cpHad && ignoreTower(1, key)) {
      m_histTool->fillCPMEtaVsPhi(m_h_IgnoreTowersHad, eta, phi);
      ttHad = 0;
      cpHad = 0;
    }

    if (!ttEm && !ttHad && !cpEm && !cpHad) continue;
    
    //  Fill in error plots

    const LVL1::Coordinate coord(phi, eta);
    const int crate = (overlap) ? converter.cpCrateOverlap(coord)
                                : converter.cpCrate(coord);
    const int cpm   = (overlap) ? converter.cpModuleOverlap(coord)
                                : converter.cpModule(coord);
    if (crate >= nCrates || cpm > nCPMs) continue;
    const int loc = crate * nCPMs + cpm - 1;
    const int cpmBins = nCrates * nCPMs;
    const int bitEm  = (1 << EMTowerMismatch);
    const int bitHad = (1 << HadTowerMismatch);
    double phiFPGA = phi;
    if (overlap) {
      const double twoPi    = 2.*M_PI;
      const double piByFour = M_PI/4.;
      if (phi > 7.*piByFour)   phiFPGA -= twoPi;
      else if (phi < piByFour) phiFPGA += twoPi;
    }
    const int loc2 = fpga(crate, phiFPGA);

    TH2F* hist1 = 0;
    TH2F* hist2 = 0;
    if (ttEm && ttEm == cpEm) { // non-zero match
      errors[loc] |= bitEm;
      hist1 = (overlap) ? m_h_EMTowerOvSIMeqDAT : m_h_EMTowerSIMeqDAT;
      hist2 = m_h_FpgaTowerSIMeqDAT;
    } else if (ttEm != cpEm) {  // mis-match
      errors[loc+cpmBins] |= bitEm;
      if (ttEm && cpEm) {       // non-zero mis-match
        hist1 = (overlap) ? m_h_EMTowerOvSIMneDAT : m_h_EMTowerSIMneDAT;
        hist2 = m_h_FpgaTowerSIMneDAT;
      } else if (!cpEm) {       // no cp
	hist1 = (overlap) ? m_h_EMTowerOvSIMnoDAT : m_h_EMTowerSIMnoDAT;
	hist2 = m_h_FpgaTowerSIMnoDAT;
      } else {                  // no tt
	hist1 = (overlap) ? m_h_EMTowerOvDATnoSIM : m_h_EMTowerDATnoSIM;
	hist2 = m_h_FpgaTowerDATnoSIM;
      }
      if (m_debug) {
        m_log << MSG::DEBUG << " EMTowerMismatch key/eta/phi/crate/cpm/tt/cp: "
              << key << "/" << eta << "/" << phi << "/" << crate << "/"
	      << cpm << "/" << ttEm << "/" << cpEm << endreq;
      }
    }
    if (hist1) m_histTool->fillCPMEtaVsPhi(hist1, eta, phi);
    if (hist2) hist2->Fill(loc, loc2);

    hist1 = 0;
    hist2 = 0;
    if (ttHad && ttHad == cpHad) { // non-zero match
      errors[loc] |= bitHad;
      hist1 = (overlap) ? m_h_HadTowerOvSIMeqDAT : m_h_HadTowerSIMeqDAT;
      hist2 = m_h_FpgaTowerSIMeqDAT;
    } else if (ttHad != cpHad) {   // mis-match
      errors[loc+cpmBins] |= bitHad;
      if (ttHad && cpHad) {        // non-zero mis-match
        hist1 = (overlap) ? m_h_HadTowerOvSIMneDAT : m_h_HadTowerSIMneDAT;
        hist2 = m_h_FpgaTowerSIMneDAT;
      } else if (!cpHad) {         // no cp
	hist1 = (overlap) ? m_h_HadTowerOvSIMnoDAT : m_h_HadTowerSIMnoDAT;
	hist2 = m_h_FpgaTowerSIMnoDAT;
      } else {                     // no tt
	hist1 = (overlap) ? m_h_HadTowerOvDATnoSIM : m_h_HadTowerDATnoSIM;
	hist2 = m_h_FpgaTowerDATnoSIM;
      }
      if (m_debug) {
        m_log << MSG::DEBUG << " HadTowerMismatch key/eta/phi/crate/cpm/tt/cp: "
              << key << "/" << eta << "/" << phi << "/" << crate << "/"
	      << cpm << "/" << ttHad << "/" << cpHad << endreq;
      }
    }
    if (hist1) m_histTool->fillCPMEtaVsPhi(hist1, eta, phi);
    if (hist2) hist2->Fill(loc, loc2+1);
  }
}

//  Compare Simulated RoIs with data

void CPMSimBSMon::compare(const CpmRoiMap& roiSimMap, const CpmRoiMap& roiMap,
                          const RodHeaderCollection* rods, ErrorVector& errors)
{
  m_log << MSG::DEBUG << "Compare Simulated RoIs with data" << endreq;

  const int nCrates = 4;
  const int nCPMs = 14;
  const int maxKey = 0xffff;
  LVL1::CPRoIDecoder decoder;
  std::vector<int> limitedRoi(nCrates);
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
    simHits &= m_roiMask;
    datHits &= m_roiMask;

    if (!simHits && !datHits) continue;

    //  Check LimitedRoISet bit

    const int crate = roi->crate();
    if (!datHits) {
      if (rods) {
	RodHeaderCollection::const_iterator rodIter  = rods->begin();
	RodHeaderCollection::const_iterator rodIterE = rods->end();
	for (; rodIter != rodIterE; ++rodIter) {
	  LVL1::RODHeader* rod = *rodIter;
	  const int rodCrate = rod->crate() - 8;
	  if (rodCrate >= 0 && rodCrate < nCrates
	      && rod->dataType() == 1 && rod->limitedRoISet()) {
	    limitedRoi[rodCrate] = 1;
	  }
	}
        rods = 0;
      }
      if (limitedRoi[crate]) continue;
    }
    
    //  Fill in error plots

    const int cpm   = roi->cpm();
    const int chip  = roi->chip();
    const int local = roi->location();
    const int locX  = crate * nCPMs + cpm - 1;
    const int locY  = chip * 8 + local;
    const int cpmBins = nCrates * nCPMs;
    const int bit = (1 << RoIMismatch);
    const LVL1::CoordinateRange coord(decoder.coordinate(roi->roiWord()));
    const double eta = coord.eta();
    const double phi = coord.phi();

    TH2F* hist1 = 0;
    TH2F* hist2 = 0;
    if (simHits == datHits) {
      errors[locX] |= bit;
      hist1 = m_h_RoISIMeqDAT;
      hist2 = m_h_RoIEtaPhiSIMeqDAT;
    } else {
      errors[locX+cpmBins] |= bit;
      if (simHits && datHits) {
        hist1 = m_h_RoISIMneDAT;
	hist2 = m_h_RoIEtaPhiSIMneDAT;
      } else if (!datHits) {
        hist1 = m_h_RoISIMnoDAT;
	hist2 = m_h_RoIEtaPhiSIMnoDAT;
      } else {
        hist1 = m_h_RoIDATnoSIM;
	hist2 = m_h_RoIEtaPhiDATnoSIM;
      }
    }
    if (hist1) hist1->Fill(locX, locY);
    if (hist2) m_histTool->fillCPMRoIEtaVsPhi(hist2, eta, phi);

    for (int thr = 0; thr < 16; ++thr) {
      if ( !((m_roiMask >> thr) & 0x1) ) continue;
      const int thrDat = (datHits >> thr) & 0x1;
      const int thrSim = (simHits >> thr) & 0x1;
      if (thrDat || thrSim) {
        if (thrDat == thrSim) m_h_RoIThreshSIMeqDAT->Fill(locX, thr);
        else                  m_h_RoIThreshSIMneDAT->Fill(locX, thr);
      }
    }

    if (m_debug) {
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

    const int nCrates = 4;
    const int nCPMs = 14;
    const int loc = crate * nCPMs + cpm - 1;
    const int cpmBins = nCrates * nCPMs;
    const int bit = (1 << CPMHitsMismatch);

    if ((simHits0 && simHits0 == datHits0) ||
        (simHits1 && simHits1 == datHits1)) errors[loc] |= bit;
    if (simHits0 != datHits0 || simHits1 != datHits1)
                                            errors[loc+cpmBins] |= bit;

    if ((simHits0 || simHits1) && simHits0 == datHits0
                               && simHits1 == datHits1) {
      m_h_CPMHitsSIMeqDAT->Fill(cpm, crate);
    }
    if ((simHits0 || simHits1) && (datHits0 || datHits1) &&
        (simHits0 != datHits0 || simHits1 != datHits1)) {
      m_h_CPMHitsSIMneDAT->Fill(cpm, crate);
    }
    if ((simHits0 || simHits1) && !datHits0 && !datHits1) {
      m_h_CPMHitsSIMnoDAT->Fill(cpm, crate);
    }
    if ((datHits0 || datHits1) && !simHits0 && !simHits1) {
      m_h_CPMHitsDATnoSIM->Fill(cpm, crate);
    }

    const int nThresh = 8;
    for (int thr = 0; thr < nThresh; ++thr) {
      const int thr2 = thr + nThresh;
      const int thrLen = 3;
      const int shift = thrLen*thr;
      const unsigned int thrMask = 0x7;
      const unsigned int d0 = (datHits0 >> shift) & thrMask;
      const unsigned int d1 = (datHits1 >> shift) & thrMask;
      const unsigned int s0 = (simHits0 >> shift) & thrMask;
      const unsigned int s1 = (simHits1 >> shift) & thrMask;
      if (d0 || s0) {
        if (d0 == s0) m_h_CPMHitsThreshSIMeqDAT->Fill(loc, thr);
        else          m_h_CPMHitsThreshSIMneDAT->Fill(loc, thr);
      }
      if (d1 || s1) {
        if (d1 == s1) m_h_CPMHitsThreshSIMeqDAT->Fill(loc, thr2);
        else          m_h_CPMHitsThreshSIMneDAT->Fill(loc, thr2);
      }
    }
  }
}

//  Compare CPM Hits and CMM Hits

void CPMSimBSMon::compare(const CpmHitsMap& cpmMap, const CmmCpHitsMap& cmmMap,
                          ErrorVector& errorsCPM, ErrorVector& errorsCMM)
{
  m_log << MSG::DEBUG << "Compare CPM Hits and CMM Hits" << endreq;

  const int nCrates = 4;
  const int nCPMs = 14;
  const int nCMMs = 2;
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
      if (cpm > nCPMs) continue;

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

    int loc  = crate * nCPMs + cpm - 1;
    const int loc2 = crate * nCMMs;
    const int cpmBins = nCrates * nCPMs;
    const int cmmBins = nCrates * nCMMs;
    const int bit = (1 << CMMHitsMismatch);
    
    if ((cpmHits0 && cpmHits0 == cmmHits0) ||
        (cpmHits1 && cpmHits1 == cmmHits1)) errorsCPM[loc] |= bit;
    if (cpmHits0 != cmmHits0 || cpmHits1 != cmmHits1)
                                            errorsCPM[loc+cpmBins] |= bit;

    TH2F* hist = 0;
    if (cpmHits1 && cpmHits1 == cmmHits1) { // hits1==>cmm 0
      errorsCMM[loc2] |= bit;
      hist = m_h_CMMHitsSIMeqDAT;
    } else if (cpmHits1 != cmmHits1) {
      errorsCMM[loc2+cmmBins] |= bit;
      if (cpmHits1 && cmmHits1) hist = m_h_CMMHitsSIMneDAT;
      else if (!cmmHits1)       hist = m_h_CMMHitsSIMnoDAT;
      else                      hist = m_h_CMMHitsDATnoSIM;
    }
    if (hist) hist->Fill(cpm, loc2);

    hist = 0;
    if (cpmHits0 && cpmHits0 == cmmHits0) { // hits0==>cmm 1
      errorsCMM[loc2+1] |= bit;
      hist = m_h_CMMHitsSIMeqDAT;
    } else if (cpmHits0 != cmmHits0) {
      errorsCMM[loc2+cmmBins+1] |= bit;
      if (cpmHits0 && cmmHits0) hist = m_h_CMMHitsSIMneDAT;
      else if (!cmmHits0)       hist = m_h_CMMHitsSIMnoDAT;
      else                      hist = m_h_CMMHitsDATnoSIM;
    }
    if (hist) hist->Fill(cpm, loc2+1);

    const int nThresh = 8;
    for (int thr = 0; thr < nThresh; ++thr) {
      const int thr2 = thr + nThresh;
      const int thrLen = 3;
      const int shift = thrLen*thr;
      const unsigned int thrMask = 0x7;
      const unsigned int d0 = (cmmHits0 >> shift) & thrMask;
      const unsigned int d1 = (cmmHits1 >> shift) & thrMask;
      const unsigned int s0 = (cpmHits0 >> shift) & thrMask;
      const unsigned int s1 = (cpmHits1 >> shift) & thrMask;
      if (d0 || s0) {
        if (d0 == s0) m_h_CMMHitsThreshSIMeqDAT->Fill(loc, thr);
        else          m_h_CMMHitsThreshSIMneDAT->Fill(loc, thr);
      }
      if (d1 || s1) {
        if (d1 == s1) m_h_CMMHitsThreshSIMeqDAT->Fill(loc, thr2);
        else          m_h_CMMHitsThreshSIMneDAT->Fill(loc, thr2);
      }
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
  const int nCrates = 4;
  const int nCMMs = 2;
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
      int loc = crate * nCMMs;
      const int cmmBins = nCrates * nCMMs;
      const int bit = (local) ? (1 << LocalSumMismatch)
                              : (1 << TotalSumMismatch);
      TH1F* hist1 = 0;
      if (cmmSimHits1 && cmmSimHits1 == cmmHits1) {
        errors[loc] |= bit;
	hist1 = m_h_SumsSIMeqDAT;
      } else if (cmmSimHits1 != cmmHits1) {
        errors[loc+cmmBins] |= bit;
	if (cmmSimHits1 && cmmHits1) hist1 = m_h_SumsSIMneDAT;
	else if (!cmmHits1)          hist1 = m_h_SumsSIMnoDAT;
	else                         hist1 = m_h_SumsDATnoSIM;
      }
      TH1F* hist0 = 0;
      if (cmmSimHits0 && cmmSimHits0 == cmmHits0) {
        errors[loc+1] |= bit;
	hist0 = m_h_SumsSIMeqDAT;
      } else if (cmmSimHits0 != cmmHits0) {
        errors[loc+cmmBins+1] |= bit;
	if (cmmSimHits0 && cmmHits0) hist0 = m_h_SumsSIMneDAT;
	else if (!cmmHits0)          hist0 = m_h_SumsSIMnoDAT;
	else                         hist0 = m_h_SumsDATnoSIM;
      }
      loc = (local) ? loc : 14;
      if (hist1) hist1->Fill(loc);
      if (hist0) hist0->Fill(loc+1);

      loc /= 2;
      const int nThresh = 8;
      for (int thr = 0; thr < nThresh; ++thr) {
        const int thr2 = thr + nThresh;
        const int thrLen = 3;
        const int shift = thrLen*thr;
        const unsigned int thrMask = 0x7;
        const unsigned int d0 = (cmmHits0 >> shift) & thrMask;
        const unsigned int d1 = (cmmHits1 >> shift) & thrMask;
        const unsigned int s0 = (cmmSimHits0 >> shift) & thrMask;
        const unsigned int s1 = (cmmSimHits1 >> shift) & thrMask;
	if (d0 || s0) {
          if (d0 == s0) m_h_SumsThreshSIMeqDAT->Fill(loc, thr);
          else          m_h_SumsThreshSIMneDAT->Fill(loc, thr);
	}
	if (d1 || s1) {
          if (d1 == s1) m_h_SumsThreshSIMeqDAT->Fill(loc, thr2);
          else          m_h_SumsThreshSIMneDAT->Fill(loc, thr2);
        }
      }
    } else {
      if (dataId == LVL1::CMMCPHits::LOCAL) {
        if (crate != nCrates-1) {
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
    for (int crate = 0; crate < nCrates-1; ++crate) {
      int loc = crate * nCMMs;
      const int cmmBins = nCrates * nCMMs;
      const int bit = (1 << RemoteSumMismatch);
      const unsigned int hd0 = hits0[crate];
      const unsigned int hd1 = hits1[crate];
      const unsigned int hs0 = hits0Sim[crate];
      const unsigned int hs1 = hits1Sim[crate];

      if (!hd0 && !hd1 && !hs0 && !hs1) continue;

      TH1F* hist1 = 0;
      if (hs1 && hs1 == hd1) {
        errors[loc] |= bit;
	hist1 = m_h_SumsSIMeqDAT;
      } else if (hs1 != hd1) {
        errors[loc+cmmBins] |= bit;
	if (hs1 && hd1) hist1 = m_h_SumsSIMneDAT;
	else if (!hd1)  hist1 = m_h_SumsSIMnoDAT;
	else            hist1 = m_h_SumsDATnoSIM;
      }
      TH1F* hist0 = 0;
      if (hs0 && hs0 == hd0) {
        errors[loc+1] |= bit;
	hist0 = m_h_SumsSIMeqDAT;
      } else if (hs0 != hd0) {
        errors[loc+cmmBins+1] |= bit;
	if (hs0 && hd0) hist0 = m_h_SumsSIMneDAT;
	else if (!hd0)  hist0 = m_h_SumsSIMnoDAT;
	else            hist0 = m_h_SumsDATnoSIM;
      }
      loc += 8;
      if (hist1) hist1->Fill(loc);
      if (hist0) hist0->Fill(loc+1);

      loc /= 2;
      const int nThresh = 8;
      for (int thr = 0; thr < nThresh; ++thr) {
        const int thr2 = thr + nThresh;
        const int thrLen = 3;
        const int shift = thrLen*thr;
        const unsigned int thrMask = 0x7;
	const unsigned int d0 = (hd0 >> shift) & thrMask;
	const unsigned int s0 = (hs0 >> shift) & thrMask;
	const unsigned int d1 = (hd1 >> shift) & thrMask;
	const unsigned int s1 = (hs1 >> shift) & thrMask;
	if (d0 || s0) {
          if (d0 == s0) m_h_SumsThreshSIMeqDAT->Fill(loc, thr);
          else          m_h_SumsThreshSIMneDAT->Fill(loc, thr);
	}
	if (d1 || s1) {
          if (d1 == s1) m_h_SumsThreshSIMeqDAT->Fill(loc, thr2);
          else          m_h_SumsThreshSIMneDAT->Fill(loc, thr2);
        }
      }
    }
  }
}

void CPMSimBSMon::setLabels(TH1* hist, bool xAxis)
{
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  // Simulation steps in red (#color[2]) depend on Trigger Menu
  axis->SetBinLabel(1+EMTowerMismatch,  "EM tt");
  axis->SetBinLabel(1+HadTowerMismatch, "Had tt");
  axis->SetBinLabel(1+RoIMismatch,      "#color[2]{RoIs}");
  axis->SetBinLabel(1+CPMHitsMismatch,  "CPMHits");
  axis->SetBinLabel(1+CMMHitsMismatch,  "CMMHits");
  axis->SetBinLabel(1+LocalSumMismatch, "#splitline{Local}{Sums}");
  axis->SetBinLabel(1+RemoteSumMismatch,"#splitline{Remote}{Sums}");
  axis->SetBinLabel(1+TotalSumMismatch, "#splitline{Total}{Sums}");
  //axis->SetLabelSize(0.045);
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
      if (eta > -2.5 && eta < 2.5 &&
                     ((*pos)->emEnergy() > 0 || (*pos)->hadEnergy() > 0)) {
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
  CpmTowerCollection* tempColl = new CpmTowerCollection;
  CpmTowerMap::const_iterator iter  = towers.begin();
  CpmTowerMap::const_iterator iterE = towers.end();
  for (; iter != iterE; ++iter) {
    LVL1::CPMTower* tt = ttCheck(iter->second, tempColl);
    const LVL1::Coordinate coord(tt->phi(), tt->eta());
    const int crate = converter.cpCrate(coord);
    if (crate >= ncrates) continue;
    crateMaps[crate].insert(std::make_pair(iter->first, tt));
  }
  // If overlap data not present take from core data
  iter  = (m_overlapPresent) ? towersOv.begin() : towers.begin();
  iterE = (m_overlapPresent) ? towersOv.end()   : towers.end();
  for (; iter != iterE; ++iter) {
    LVL1::CPMTower* tt = ttCheck(iter->second, tempColl);
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
  delete tempColl;
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

// Test for tower in ignore list

bool CPMSimBSMon::ignoreTower(int layer, int key)
{
  bool ignore = false;
  const std::vector<int>* list = (layer == 0) ? &m_ignoreTowersEm
                                              : &m_ignoreTowersHad;
  std::vector<int>::const_iterator iter  = list->begin();
  std::vector<int>::const_iterator iterE = list->end();
  for (; iter != iterE; ++iter) {
    if (key == *iter) {
      ignore = true;
      break;
    }
  }
  return ignore;
}

// Return a tower with zero energy if parity bit is set

LVL1::CPMTower* CPMSimBSMon::ttCheck(LVL1::CPMTower* tt,
                                     CpmTowerCollection* coll)
{
  const LVL1::DataError emError(tt->emError());
  const LVL1::DataError hadError(tt->hadError());
  const int emParity  = emError.get(LVL1::DataError::Parity);
  const int hadParity = hadError.get(LVL1::DataError::Parity);
  if ((emParity && tt->emEnergy()) || (hadParity && tt->hadEnergy())) {
    const int peak = tt->peak();
    std::vector<int> emEnergyVec(tt->emEnergyVec());
    std::vector<int> hadEnergyVec(tt->hadEnergyVec());
    std::vector<int> emErrorVec(tt->emErrorVec());
    std::vector<int> hadErrorVec(tt->hadErrorVec());
    if (emParity)  emEnergyVec[peak]  = 0;
    if (hadParity) hadEnergyVec[peak] = 0;
    LVL1::CPMTower* ct = new LVL1::CPMTower(tt->phi(), tt->eta(),
        emEnergyVec, emErrorVec, hadEnergyVec, hadErrorVec, peak);
    coll->push_back(ct);
    return ct;
  }
  return tt;
}
