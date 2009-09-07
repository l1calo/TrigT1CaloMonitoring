// ********************************************************************
//
// NAME:     JEPSimBSMon.cxx
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

#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"

#include "TrigT1CaloEvent/CMMJetHits.h"
#include "TrigT1CaloEvent/JEMHits.h"
#include "TrigT1CaloEvent/JetElement.h"
#include "TrigT1CaloEvent/JEMRoI.h"
#include "TrigT1CaloEvent/CMMRoI.h"
#include "TrigT1CaloEvent/JEMEtSums.h"
#include "TrigT1CaloEvent/CMMEtSums.h"
#include "TrigT1CaloEvent/RODHeader.h"
#include "TrigT1CaloEvent/TriggerTower.h"
#include "TrigT1CaloUtils/CoordToHardware.h"
#include "TrigT1CaloUtils/JetAlgorithm.h"
#include "TrigT1CaloToolInterfaces/IL1JEPHitsTools.h"
#include "TrigT1CaloToolInterfaces/IL1JetTools.h"
#include "TrigT1CaloToolInterfaces/IL1JetElementTools.h"
#include "TrigT1CaloToolInterfaces/IL1JEPEtSumsTools.h"
#include "TrigT1Interfaces/Coordinate.h"
#include "TrigT1Interfaces/JEPRoIDecoder.h"
#include "TrigT1Interfaces/TrigT1CaloDefs.h"

#include "TrigT1CaloMonitoring/JEPSimBSMon.h"


/*---------------------------------------------------------*/
JEPSimBSMon::JEPSimBSMon(const std::string & type, 
			 const std::string & name,
			 const IInterface* parent)
  : ManagedMonitorToolBase(type, name, parent),
    m_storeGate("StoreGateSvc", name),
    m_jepHitsTool("LVL1::L1JEPHitsTools/L1JEPHitsTools"),
    m_jetTool("LVL1::L1JetTools/L1JetTools"),
    m_jetElementTool("LVL1::L1JetElementTools/L1JetElementTools"),
    m_etSumsTool("LVL1::L1JEPEtSumsTools/L1JEPEtSumsTools"),
    m_log(msgSvc(), name), m_debug(false),
    m_monGroup(0), m_phiScale(16./M_PI), m_events(0)
/*---------------------------------------------------------*/
{
  declareInterface<IMonitorToolBase>(this); 

  declareProperty("JEPHitsTool", m_jepHitsTool);
  declareProperty("JetTool", m_jetTool);
  declareProperty("JetElementTool", m_jetElementTool);
  declareProperty("JEPEtSumsTool", m_etSumsTool);

  declareProperty("JetElementLocation",
                 m_jetElementLocation =
		                 LVL1::TrigT1CaloDefs::JetElementLocation);
  declareProperty("JetElementLocationOverlap",
                 m_jetElementLocationOverlap =
		           LVL1::TrigT1CaloDefs::JetElementLocation+"Overlap");
  declareProperty("JEMHitsLocation",
                 m_jemHitsLocation    = LVL1::TrigT1CaloDefs::JEMHitsLocation);
  declareProperty("CMMJetHitsLocation",
                 m_cmmJetHitsLocation =
		                 LVL1::TrigT1CaloDefs::CMMJetHitsLocation);
  declareProperty("JEMRoILocation",
                 m_jemRoiLocation     = LVL1::TrigT1CaloDefs::JEMRoILocation);
  declareProperty("CMMRoILocation",
                 m_cmmRoiLocation     = LVL1::TrigT1CaloDefs::CMMRoILocation);
  declareProperty("JEMEtSumsLocation",
                 m_jemEtSumsLocation = LVL1::TrigT1CaloDefs::JEMEtSumsLocation);
  declareProperty("CMMEtSumsLocation",
                 m_cmmEtSumsLocation = LVL1::TrigT1CaloDefs::CMMEtSumsLocation);
  declareProperty("TriggerTowerLocation",
                 m_triggerTowerLocation =
		                 LVL1::TrigT1CaloDefs::TriggerTowerLocation);
  declareProperty("RodHeaderLocation",
                 m_rodHeaderLocation = "RODHeaders");

  declareProperty("RootDirectory", m_rootDir = "L1Calo");

  declareProperty("CompareWithSimulation", m_compareWithSim = true,
                  "Include the checks that run simulation tools");
  declareProperty("EventSamples", m_eventSamples = 10,
                  "Number of Error Event Number Samples");
}

/*---------------------------------------------------------*/
JEPSimBSMon::~JEPSimBSMon()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode JEPSimBSMon:: initialize()
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

    sc = m_jetElementTool.retrieve();
    if( sc.isFailure() ) {
      m_log << MSG::ERROR << "Unable to locate Tool L1JetElementTools"
            << endreq;
      return sc;
    }

    sc = m_jetTool.retrieve();
    if( sc.isFailure() ) {
      m_log << MSG::ERROR << "Unable to locate Tool L1JetTools" << endreq;
      return sc;
    }

    sc = m_jepHitsTool.retrieve();
    if( sc.isFailure() ) {
      m_log << MSG::ERROR << "Unable to locate Tool L1JEPHitsTools" << endreq;
      return sc;
    }

    sc = m_etSumsTool.retrieve();
    if( sc.isFailure() ) {
      m_log << MSG::ERROR << "Unable to locate Tool L1JEPEtSumsTools" << endreq;
      return sc;
    }
  }

  m_log << MSG::INFO << "JEPSimBSMon initialised" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode JEPSimBSMon::bookHistograms(bool isNewEventsBlock,
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

  std::string dir1(m_rootDir + "/JEM/Errors/Transmission_Simulation");
  MonGroup monShift( this, dir1, shift, run );
  MonGroup monExpert( this, dir1, expert, run );
  MonGroup monElements( this, dir1 + "/PPM2Elements", expert, run );
  MonGroup monRoIs( this, dir1 + "/Elements2RoIs", expert, run );
  MonGroup monHits( this, dir1 + "/RoIs2Hits", expert, run );
  MonGroup monEnergy( this, dir1 + "/Elements2Energy", expert, run );
  MonGroup monEvent1( this, dir1 + "/MismatchEventNumbers", expert, run, "",
                                                            "eventSample" );
  std::string dir2(m_rootDir + "/JEM-CMM/Errors/Transmission_Simulation");
  MonGroup monHits2( this, dir2 + "/JEM2CMMHits", expert, run );
  MonGroup monHitSums( this, dir2 + "/Hits2Sums", expert, run );
  MonGroup monEnergy2( this, dir2 + "/JEM2CMMEnergy", expert, run );
  MonGroup monEnergySums( this, dir2 + "/Energy2Sums", expert, run );
  MonGroup monEvent2( this, dir2 + "/MismatchEventNumbers", expert, run, "",
                                                            "eventSample" );
  //temporary
  MonGroup monTemp( this, m_rootDir + "/3_JEP_TransmissionAndPerformance", 
                                                               shift, run );

  // JetElements

  if (m_compareWithSim) {

  m_monGroup = &monElements;

  m_h_EMEleSIMeqDAT = bookEtaPhi("jem_em_2d_etaPhi_jetEl_Sim_eq_core",
    "Core Jet Elements EM Data/Simulation Non-zero Matches");
  m_h_EMEleSIMneDAT = bookEtaPhi("jem_em_2d_etaPhi_jetEl_Sim_ne_core",
    "Core Jet Elements EM Data/Simulation Non-zero Mismatches");
  m_h_EMEleSIMnoDAT = bookEtaPhi("jem_em_2d_etaPhi_jetEl_Sim_no_core",
    "Core Jet Elements EM Simulation but no Data");
  m_h_EMEleDATnoSIM = bookEtaPhi("jem_em_2d_etaPhi_jetEl_Core_no_sim",
    "Core Jet Elements EM Data but no Simulation");
  m_h_HadEleSIMeqDAT = bookEtaPhi("jem_had_2d_etaPhi_jetEl_Sim_eq_core",
    "Core Jet Elements HAD Data/Simulation Non-zero Matches");
  m_h_HadEleSIMneDAT = bookEtaPhi("jem_had_2d_etaPhi_jetEl_Sim_ne_core",
    "Core Jet Elements HAD Data/Simulation Non-zero Mismatches");
  m_h_HadEleSIMnoDAT = bookEtaPhi("jem_had_2d_etaPhi_jetEl_Sim_no_core",
    "Core Jet Elements HAD Simulation but no Data");
  m_h_HadEleDATnoSIM = bookEtaPhi("jem_had_2d_etaPhi_jetEl_Core_no_sim",
    "Core Jet Elements HAD Data but no Simulation");
  m_h_EMEleOvSIMeqDAT = bookEtaPhi("jem_em_2d_etaPhi_jetEl_Sim_eq_overlap",
    "Overlap Jet Elements EM Data/Simulation Non-zero Matches");
  m_h_EMEleOvSIMneDAT = bookEtaPhi("jem_em_2d_etaPhi_jetEl_Sim_ne_overlap",
    "Overlap Jet Elements EM Data/Simulation Non-zero Mismatches");
  m_h_EMEleOvSIMnoDAT = bookEtaPhi("jem_em_2d_etaPhi_jetEl_Sim_no_overlap",
    "Overlap Jet Elements EM Simulation but no Data");
  m_h_EMEleOvDATnoSIM = bookEtaPhi("jem_em_2d_etaPhi_jetEl_Overlap_no_sim",
    "Overlap Jet Elements EM Data but no Simulation");
  m_h_HadEleOvSIMeqDAT = bookEtaPhi("jem_had_2d_etaPhi_jetEl_Sim_eq_overlap",
    "Overlap Jet Elements HAD Data/Simulation Non-zero Matches");
  m_h_HadEleOvSIMneDAT = bookEtaPhi("jem_had_2d_etaPhi_jetEl_Sim_ne_overlap",
    "Overlap Jet Elements HAD Data/Simulation Non-zero Mismatches");
  m_h_HadEleOvSIMnoDAT = bookEtaPhi("jem_had_2d_etaPhi_jetEl_Sim_no_overlap",
    "Overlap Jet Elements HAD Simulation but no Data");
  m_h_HadEleOvDATnoSIM = bookEtaPhi("jem_had_2d_etaPhi_jetEl_Overlap_no_sim",
    "Overlap Jet Elements HAD Data but no Simulation");

  //  RoIs

  m_monGroup = &monRoIs;

  m_h_RoISIMeqDAT = book2F("jem_2d_roi_Sim_eq_data",
     "JEM RoI Data/Simulation Non-zero Matches;Crate/Module;Frame/Local Coord",
             32, 0, 32, 32, 0, 32);
  setLabelsCMFC(m_h_RoISIMeqDAT);
  m_h_RoISIMneDAT = book2F("jem_2d_roi_Sim_ne_data",
   "JEM RoI Data/Simulation Non-zero Mismatches;Crate/Module;Frame/Local Coord",
             32, 0, 32, 32, 0, 32);
  setLabelsCMFC(m_h_RoISIMneDAT);
  m_h_RoISIMnoDAT = book2F("jem_2d_roi_Sim_no_data",
            "JEM RoI Simulation but no Data;Crate/Module;Frame/Local Coord",
	     32, 0, 32, 32, 0, 32);
  setLabelsCMFC(m_h_RoISIMnoDAT);
  m_h_RoIDATnoSIM = book2F("jem_2d_roi_Data_no_sim",
            "JEM RoI Data but no Simulation;Crate/Module;Frame/Local Coord",
	     32, 0, 32, 32, 0, 32);
  setLabelsCMFC(m_h_RoIDATnoSIM);
  m_h_RoIThreshSIMeqDAT = book2F("jem_2d_roi_Thresh_sim_eq_data",
     "JEM RoI Data/Simulation Threshold Matches;Crate/Module",
             32, 0, 32, 12, 0, 12);
  setLabelsCMT(m_h_RoIThreshSIMeqDAT);
  m_h_RoIThreshSIMneDAT = book2F("jem_2d_roi_Thresh_sim_ne_data",
     "JEM RoI Data/Simulation Threshold Mismatches;Crate/Module",
             32, 0, 32, 12, 0, 12);
  setLabelsCMT(m_h_RoIThreshSIMneDAT);
  m_h_RoIEtaPhiSIMeqDAT = bookEtaPhi("jem_2d_etaPhi_roi_Sim_eq_data",
     "JEM RoI Data/Simulation Non-zero Matches", true);
  m_h_RoIEtaPhiSIMneDAT = bookEtaPhi("jem_2d_etaPhi_roi_Sim_ne_data",
     "JEM RoI Data/Simulation Non-zero Mismatches", true);
  m_h_RoIEtaPhiSIMnoDAT = bookEtaPhi("jem_2d_etaPhi_roi_Sim_no_data",
     "JEM RoI Simulation but no Data", true);
  m_h_RoIEtaPhiDATnoSIM = bookEtaPhi("jem_2d_etaPhi_roi_Data_no_sim",
     "JEM RoI Data but no Simulation", true);

  // JEMHits

  m_monGroup = &monHits;

  m_h_JEMHitsSIMeqDAT = book2F("jem_2d_thresh_Sim_eq_data",
     "JEM Hits Data/Simulation Non-zero Matches;Module;Crate",
             16, 0, 16, 2, 0, 2);
  setLabelsMC(m_h_JEMHitsSIMeqDAT);
  m_h_JEMHitsSIMneDAT = book2F("jem_2d_thresh_Sim_ne_data",
     "JEM Hits Data/Simulation Non-zero Mismatches;Module;Crate",
             16, 0, 16, 2, 0, 2);
  setLabelsMC(m_h_JEMHitsSIMneDAT);
  m_h_JEMHitsSIMnoDAT = book2F("jem_2d_thresh_Sim_no_data",
     "JEM Hits Simulation but no Data;Module;Crate", 16, 0, 16, 2, 0, 2);
  setLabelsMC(m_h_JEMHitsSIMnoDAT);
  m_h_JEMHitsDATnoSIM = book2F("jem_2d_thresh_Data_no_sim",
     "JEM Hits Data but no Simulation;Module;Crate", 16, 0, 16, 2, 0, 2);
  setLabelsMC(m_h_JEMHitsDATnoSIM);
  m_h_JEMHitsThreshSIMeqDAT = book2F("jem_2d_thresh_Thresh_sim_eq_data",
     "JEM Hits Data/Simulation Threshold Matches;Crate/Module",
             32, 0, 32, 12, 0, 12);
  setLabelsCMT(m_h_JEMHitsThreshSIMeqDAT);
  m_h_JEMHitsThreshSIMneDAT = book2F("jem_2d_thresh_Thresh_sim_ne_data",
     "JEM Hits Data/Simulation Threshold Mismatches;Crate/Module",
             32, 0, 32, 12, 0, 12);
  setLabelsCMT(m_h_JEMHitsThreshSIMneDAT);

  } // end if (m_compareWithSim)

  // CMMHits

  m_monGroup = &monHits2;

  m_h_CMMHitsSIMeqDAT = book2F("cmm_2d_thresh_JEM_eq_CMM",
     "CMM Hits/JEM Hits Non-zero Matches;Module;Crate", 16, 0, 16, 2, 0, 2);
  setLabelsMC(m_h_CMMHitsSIMeqDAT);
  m_h_CMMHitsSIMneDAT = book2F("cmm_2d_thresh_JEM_ne_CMM",
     "CMM Hits/JEM Hits Non-zero Mismatches;Module;Crate", 16, 0, 16, 2, 0, 2);
  setLabelsMC(m_h_CMMHitsSIMneDAT);
  m_h_CMMHitsSIMnoDAT = book2F("cmm_2d_thresh_JEM_no_CMM",
     "JEM Hits but no CMM Hits;Module;Crate", 16, 0, 16, 2, 0, 2);
  setLabelsMC(m_h_CMMHitsSIMnoDAT);
  m_h_CMMHitsDATnoSIM = book2F("cmm_2d_thresh_CMM_no_JEM",
     "CMM Hits but no JEM Hits;Module;Crate", 16, 0, 16, 2, 0, 2);
  setLabelsMC(m_h_CMMHitsDATnoSIM);
  m_h_CMMHitsThreshSIMeqDAT = book2F("cmm_2d_thresh_Thresh_JEM_eq_CMM",
     "CMM Hits/JEM Hits Threshold Matches;Crate/Module",
             32, 0, 32, 12, 0, 12);
  setLabelsCMT(m_h_CMMHitsThreshSIMeqDAT);
  m_h_CMMHitsThreshSIMneDAT = book2F("cmm_2d_thresh_Thresh_JEM_ne_CMM",
     "CMM Hits/JEM Hits Threshold Mismatches;Crate/Module",
             32, 0, 32, 12, 0, 12);
  setLabelsCMT(m_h_CMMHitsThreshSIMneDAT);

  m_monGroup = &monHitSums;

  // Local/Remote/Total sums

  m_h_SumsSIMeqDAT = book1F("cmm_1d_thresh_Sums_sim_eq_data",
     "CMM Hit Sums Data/Simulation Non-zero Matches", 6, 0, 6);
  setLabelsSH(m_h_SumsSIMeqDAT);
  m_h_SumsSIMneDAT = book1F("cmm_1d_thresh_Sums_sim_ne_data",
     "CMM Hit Sums Data/Simulation Non-zero Mismatches", 6, 0, 6);
  setLabelsSH(m_h_SumsSIMneDAT);
  m_h_SumsSIMnoDAT = book1F("cmm_1d_thresh_Sums_sim_no_data",
     "CMM Hit Sums Simulation but no Data", 6, 0, 6);
  setLabelsSH(m_h_SumsSIMnoDAT);
  m_h_SumsDATnoSIM = book1F("cmm_1d_thresh_Sums_data_no_sim",
     "CMM Hit Sums Data but no Simulation", 6, 0, 6);
  setLabelsSH(m_h_SumsDATnoSIM);
  m_h_SumsThreshSIMeqDAT = book2F("cmm_2d_thresh_Sums_thresh_sim_eq_data",
     "CMM Hit Sums Data/Simulation Threshold Matches",
             6, 0, 6, 16, 0, 16);
  setLabelsSHF(m_h_SumsThreshSIMeqDAT);
  m_h_SumsThreshSIMneDAT = book2F("cmm_2d_thresh_Sums_thresh_sim_ne_data",
     "CMM Hit Sums Data/Simulation Threshold Mismatches",
             6, 0, 6, 16, 0, 16);
  setLabelsSHF(m_h_SumsThreshSIMneDAT);

  if (m_compareWithSim) {

  // JEMEtSums

  m_monGroup = &monEnergy;

  m_h_jemEtSumsSIMeqDAT = book2F("jem_2d_energy_Sim_eq_data",
    "JEM EtSums Data/Simulation Non-zero Matches;Crate/Module",
             32, 0, 32, 3, 0, 3);
  setLabelsJES(m_h_jemEtSumsSIMeqDAT);
  m_h_jemEtSumsSIMneDAT = book2F("jem_2d_energy_Sim_ne_data",
    "JEM EtSums Data/Simulation Non-zero Mismatches;Crate/Module",
             32, 0, 32, 3, 0, 3);
  setLabelsJES(m_h_jemEtSumsSIMneDAT);
  m_h_jemEtSumsSIMnoDAT = book2F("jem_2d_energy_Sim_no_data",
    "JEM EtSums Simulation but no Data;Crate/Module",
             32, 0, 32, 3, 0, 3);
  setLabelsJES(m_h_jemEtSumsSIMnoDAT);
  m_h_jemEtSumsDATnoSIM = book2F("jem_2d_energy_Data_no_sim",
    "JEM EtSums Data but no Simulation;Crate/Module",
             32, 0, 32, 3, 0, 3);
  setLabelsJES(m_h_jemEtSumsDATnoSIM);

  } // end if (m_compareWithSim)

  // CMMEtSums

  m_monGroup = &monEnergy2;

  m_h_cmmEtSumsSIMeqDAT = book2F("cmm_2d_energy_JEM_eq_CMM",
    "CMM EtSums/JEM EtSums Non-zero Matches;Crate/Module",
             32, 0, 32, 3, 0, 3);
  setLabelsJES(m_h_cmmEtSumsSIMeqDAT);
  m_h_cmmEtSumsSIMneDAT = book2F("cmm_2d_energy_JEM_ne_CMM",
    "CMM EtSums/JEM EtSums Non-zero Mismatches;Crate/Module",
             32, 0, 32, 3, 0, 3);
  setLabelsJES(m_h_cmmEtSumsSIMneDAT);
  m_h_cmmEtSumsSIMnoDAT = book2F("cmm_2d_energy_JEM_no_CMM",
    "JEM EtSums but no CMM EtSums;Crate/Module",
             32, 0, 32, 3, 0, 3);
  setLabelsJES(m_h_cmmEtSumsSIMnoDAT);
  m_h_cmmEtSumsDATnoSIM = book2F("cmm_2d_energy_CMM_no_JEM",
    "CMM EtSums but no JEM EtSums;Crate/Module",
             32, 0, 32, 3, 0, 3);
  setLabelsJES(m_h_cmmEtSumsDATnoSIM);

  m_monGroup = &monEnergySums;

  // Energy Crate/System sums

  m_h_EnSumsSIMeqDAT = book2F("cmm_2d_energy_Sums_sim_eq_data",
    "Energy Totals Data/Simulation Non-zero Matches", 5, 0, 5, 5, 0, 5);
  setLabelsEnTot(m_h_EnSumsSIMeqDAT);
  m_h_EnSumsSIMneDAT = book2F("cmm_2d_energy_Sums_sim_ne_data",
    "Energy Totals Data/Simulation Non-zero Mismatches", 5, 0, 5, 5, 0, 5);
  setLabelsEnTot(m_h_EnSumsSIMneDAT);
  m_h_EnSumsSIMnoDAT = book2F("cmm_2d_energy_Sums_sim_no_data",
    "Energy Totals Simulation but no Data", 5, 0, 5, 5, 0, 5);
  setLabelsEnTot(m_h_EnSumsSIMnoDAT);
  m_h_EnSumsDATnoSIM = book2F("cmm_2d_energy_Sums_data_no_sim",
    "Energy Totals Data but no Simulation", 5, 0, 5, 5, 0, 5);
  setLabelsEnTot(m_h_EnSumsDATnoSIM);
  m_h_EnSumsThreshSIMeqDAT = book2F("cmm_2d_energy_EtMaps_thresh_sim_eq_data",
    "Et Maps Data/Simulation Threshold Matches;;Threshold", 4, 0, 4, 8, 0, 8);
  setLabelsEnTotThr(m_h_EnSumsThreshSIMeqDAT);
  m_h_EnSumsThreshSIMneDAT = book2F("cmm_2d_energy_EtMaps_thresh_sim_ne_data",
    "Et Maps Data/Simulation Threshold Mismatches;;Threshold",
                                                            4, 0, 4, 8, 0, 8);
  setLabelsEnTotThr(m_h_EnSumsThreshSIMneDAT);

  // Summary

  m_monGroup = &monExpert;

  m_h_JEPeqSIM = book2F("jem_2d_Sim_eq_data_overview",
   "JEP Transmission/Comparison with Simulation Overview - Events with Matches;Crate/Module",
             36, 0, 36, NumberOfSummaryBins, 0, NumberOfSummaryBins);
  setLabels(m_h_JEPeqSIM);

  m_h_JEPneSIM = book2F("jem_2d_Sim_ne_data_overview",
"JEP Transmission/Comparison with Simulation Overview - Events with Mismatches;Crate/Module",
             36, 0, 36, NumberOfSummaryBins, 0, NumberOfSummaryBins);
  setLabels(m_h_JEPneSIM);

  //temporary duplicates
  m_monGroup = &monTemp;
  m_h_JEPeqSIM2 = book2F("JEP_TransCheck",
   "JEP Transmission/Comparison with Simulation Overview - Events with Matches;Crate/Module",
             36, 0, 36, NumberOfSummaryBins, 0, NumberOfSummaryBins);
  setLabels(m_h_JEPeqSIM2);

  m_h_JEPneSIM2 = book2F("JEP_Calc_Error",
"JEP Transmission/Comparison with Simulation Overview - Events with Mismatches;Crate/Module",
             36, 0, 36, NumberOfSummaryBins, 0, NumberOfSummaryBins);
  setLabels(m_h_JEPneSIM2);

  m_h_JEPeqSIM3 = book2F("hadJE_TransPerfCheck",
   "JEP Transmission/Comparison with Simulation Overview - Events with Matches;Crate/Module",
             36, 0, 36, NumberOfSummaryBins, 0, NumberOfSummaryBins);
  setLabels(m_h_JEPeqSIM3);

  m_h_JEPneSIM3 = book2F("emJE_TransPerfCheck",
"JEP Transmission/Comparison with Simulation Overview - Events with Mismatches;Crate/Module",
             36, 0, 36, NumberOfSummaryBins, 0, NumberOfSummaryBins);
  setLabels(m_h_JEPneSIM3);
  //end temporary duplicates

  m_monGroup = &monShift;

  m_h_JEPneSIMSummary = book1F("jem_1d_Sim_ne_data_summary",
   "JEP Transmission/Comparison with Simulation Mismatch Summary;;Events",
    NumberOfSummaryBins, 0, NumberOfSummaryBins);
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+EMElementMismatch,  "EM je");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+HadElementMismatch, "Had je");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+RoIMismatch,        "RoIs");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+JEMHitsMismatch,    "JEMHits");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+CMMJetHitsMismatch, "CMMHits");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+LocalJetMismatch,   "Local");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+RemoteJetMismatch,  "Remote");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+TotalJetMismatch,   "Total");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+JetEtMismatch,      "JetEt");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+JetEtRoIMismatch, "JetEt RoI");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+JEMEtSumsMismatch,  "JEMSums");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+CMMEtSumsMismatch,  "CMMSums");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+LocalEnergyMismatch, "Local");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+RemoteEnergyMismatch,"Remote");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+TotalEnergyMismatch, "Total");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+SumEtMismatch,    "SumEt");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+MissingEtMismatch,"MissingEt");
  m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+EnergyRoIMismatch,"Engy RoIs");
  if (!m_compareWithSim) { // Transmission only - grey out simulation
    m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+EMElementMismatch,
                                                        "#color[16]{EM je}");
    m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+HadElementMismatch,
                                                        "#color[16]{Had je}");
    m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+RoIMismatch,
                                                        "#color[16]{RoIs}");
    m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+JEMHitsMismatch,
                                                        "#color[16]{JEMHits}");
    m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+LocalJetMismatch,
                                                        "#color[16]{Local}");
    m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+TotalJetMismatch,
                                                        "#color[16]{Total}");
    m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+JetEtMismatch,
                                                        "#color[16]{JetEt}");
    m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+JEMEtSumsMismatch,
                                                        "#color[16]{JEMSums}");
    m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+LocalEnergyMismatch,
                                                        "#color[16]{Local}");
    m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+TotalEnergyMismatch,
                                                        "#color[16]{Total}");
    m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+SumEtMismatch,
                                                        "#color[16]{SumEt}");
    m_h_JEPneSIMSummary->GetXaxis()->SetBinLabel(1+MissingEtMismatch,
                                                      "#color[16]{MissingEt}");
  }
  m_h_JEPneSIMSummary->GetXaxis()->SetLabelSize(0.045);

  // Mismatch Event Number Samples

  m_monGroup = &monEvent1;

  m_sampleCounts.clear();
  m_sampleCounts.resize(8*32+9, 0);
  TH2I* hist = 0;
  m_sampleHists.clear();
  m_sampleHists.resize(9, hist);
  if (m_compareWithSim) {
    hist = book2I("jem_em_2d_jetEl_Mismatch_events",
           "Jet Elements EM Mismatch Event Numbers;Sample;Crate/Module",
	   m_eventSamples, 0, m_eventSamples, 32, 0, 32);
    setLabelsJMS(hist);
    m_sampleHists[0] = hist;
    hist = book2I("jem_had_2d_jetEl_Mismatch_events",
           "Jet Elements Had Mismatch Event Numbers;Sample;Crate/Module",
	   m_eventSamples, 0, m_eventSamples, 32, 0, 32);
    setLabelsJMS(hist);
    m_sampleHists[1] = hist;
    hist = book2I("jem_2d_roi_Mismatch_events",
           "JEM RoIs Mismatch Event Numbers;Sample;Crate/Module",
	   m_eventSamples, 0, m_eventSamples, 32, 0, 32);
    setLabelsJMS(hist);
    m_sampleHists[2] = hist;
    hist = book2I("jem_2d_thresh_Mismatch_events",
           "JEM Hits Mismatch Event Numbers;Sample;Crate/Module",
	   m_eventSamples, 0, m_eventSamples, 32, 0, 32);
    setLabelsJMS(hist);
    m_sampleHists[3] = hist;
    hist = book2I("jem_2d_energy_Mismatch_events",
           "JEM Energy Mismatch Event Numbers;Sample;Crate/Module",
	   m_eventSamples, 0, m_eventSamples, 32, 0, 32);
    setLabelsJMS(hist);
    m_sampleHists[5] = hist;
  }

  m_monGroup = &monEvent2;

  hist = book2I("cmm_2d_thresh_Mismatch_events",
           "CMM Hits Mismatch Event Numbers;Sample;Crate/Module",
	   m_eventSamples, 0, m_eventSamples, 32, 0, 32);
  setLabelsJMS(hist);
  m_sampleHists[4] = hist;
  hist = book2I("cmm_2d_energy_Mismatch_events",
           "CMM Energy Mismatch Event Numbers;Sample;Crate/Module",
	   m_eventSamples, 0, m_eventSamples, 32, 0, 32);
  setLabelsJMS(hist);
  m_sampleHists[6] = hist;
  hist = book2I("cmm_2d_thresh_Sums_mismatch_events",
           "CMM Hit Sums Mismatch Event Numbers;Sample",
	   m_eventSamples, 0, m_eventSamples, 8, 0, 8);
  hist->GetYaxis()->SetBinLabel(1, "Modules 0");
  hist->GetYaxis()->SetBinLabel(2, "Modules 1");
  hist->GetYaxis()->SetBinLabel(3, "Local 0");
  hist->GetYaxis()->SetBinLabel(4, "Local 1");
  hist->GetYaxis()->SetBinLabel(5, "Remote");
  hist->GetYaxis()->SetBinLabel(6, "Total");
  hist->GetYaxis()->SetBinLabel(7, "JetEt");
  hist->GetYaxis()->SetBinLabel(8, "JetEt RoI");
  if (!m_compareWithSim) {
    hist->GetYaxis()->SetBinLabel(3, "#color[16]{Local 0}");
    hist->GetYaxis()->SetBinLabel(4, "#color[16]{Local 1}");
    hist->GetYaxis()->SetBinLabel(6, "#color[16]{Total}");
    hist->GetYaxis()->SetBinLabel(7, "#color[16]{JetEt}");
  }
  if (m_eventSamples <= 10) setLabelsXNUM(hist, 1, m_eventSamples);
  m_sampleHists[7] = hist;
  hist = book2I("cmm_2d_energy_Sums_mismatch_events",
           "CMM Energy Sums Mismatch Event Numbers;Sample",
	   m_eventSamples, 0, m_eventSamples, 9, 0, 9);
  hist->GetYaxis()->SetBinLabel(1, "Modules 0");
  hist->GetYaxis()->SetBinLabel(2, "Modules 1");
  hist->GetYaxis()->SetBinLabel(3, "Local 0");
  hist->GetYaxis()->SetBinLabel(4, "Local 1");
  hist->GetYaxis()->SetBinLabel(5, "Remote");
  hist->GetYaxis()->SetBinLabel(6, "Total");
  hist->GetYaxis()->SetBinLabel(7, "SumEt");
  hist->GetYaxis()->SetBinLabel(8, "MissingEt");
  hist->GetYaxis()->SetBinLabel(9, "Engy RoIs");
  if (!m_compareWithSim) {
    hist->GetYaxis()->SetBinLabel(3, "#color[16]{Local 0}");
    hist->GetYaxis()->SetBinLabel(4, "#color[16]{Local 1}");
    hist->GetYaxis()->SetBinLabel(6, "#color[16]{Total}");
    hist->GetYaxis()->SetBinLabel(7, "#color[16]{SumEt}");
    hist->GetYaxis()->SetBinLabel(8, "#color[16]{MissingEt}");
  }
  if (m_eventSamples <= 10) setLabelsXNUM(hist, 1, m_eventSamples);
  m_sampleHists[8] = hist;

  m_events = 0;

  } // end if (isNewRun ...

  m_log << MSG::DEBUG << "Leaving bookHistograms" << endreq;
  
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode JEPSimBSMon::fillHistograms()
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "fillHistograms entered" << endreq;

  // NB. Collection retrieves wrapped in m_storeGate->contains<..>(..)
  // are for those not expected to be on ESD. They should be on bytestream.

  StatusCode sc;

  //Retrieve Trigger Towers from SG
  const TriggerTowerCollection* triggerTowerTES = 0; 
  if (m_compareWithSim) {
    sc = m_storeGate->retrieve(triggerTowerTES, m_triggerTowerLocation); 
    if( sc.isFailure()  ||  !triggerTowerTES ) {
      m_log << MSG::DEBUG<< "No Trigger Tower container found"<< endreq; 
    }
  }

  //Retrieve Core and Overlap Jet Elements from SG
  const JetElementCollection* jetElementTES = 0; 
  const JetElementCollection* jetElementOvTES = 0; 
  if (m_compareWithSim) {
    sc = m_storeGate->retrieve(jetElementTES, m_jetElementLocation); 
    if( sc.isFailure()  ||  !jetElementTES ) {
      m_log << MSG::DEBUG<< "No Core Jet Element container found"<< endreq; 
    }
    if (m_storeGate->contains<JetElementCollection>(m_jetElementLocationOverlap)) {
      sc = m_storeGate->retrieve(jetElementOvTES, m_jetElementLocationOverlap);
    } else sc = StatusCode::FAILURE;
    if( sc.isFailure()  ||  !jetElementOvTES ) {
      m_log << MSG::DEBUG<< "No Overlap Jet Element container found"<< endreq;
    }
  }
  
  //Retrieve JEM RoIs from SG
  const JemRoiCollection* jemRoiTES = 0;
  const RodHeaderCollection* rodTES = 0;
  if (m_compareWithSim) {
    sc = m_storeGate->retrieve( jemRoiTES, m_jemRoiLocation);
    if( sc.isFailure()  ||  !jemRoiTES  ||  jemRoiTES->empty() ) {
      m_log << MSG::DEBUG << "No DAQ JEM RoIs container found" << endreq; 
    } else {
      if (m_storeGate->contains<RodHeaderCollection>(m_rodHeaderLocation)) {
        sc = m_storeGate->retrieve( rodTES, m_rodHeaderLocation);
      } else sc = StatusCode::FAILURE;
      if( sc.isFailure()  ||  !rodTES ) {
        m_log << MSG::DEBUG << "No ROD Header container found"<< endreq;
      }
    }
  }
  
  //Retrieve JEM Hits from SG
  const JemHitsCollection* jemHitsTES = 0;
  sc = m_storeGate->retrieve( jemHitsTES, m_jemHitsLocation);
  if( sc.isFailure()  ||  !jemHitsTES ) {
    m_log << MSG::DEBUG << "No JEM Hits container found"<< endreq; 
  }
  
  //Retrieve CMM-Jet Hits from SG
  const CmmJetHitsCollection* cmmJetHitsTES = 0;
  sc = m_storeGate->retrieve( cmmJetHitsTES, m_cmmJetHitsLocation);
  if( sc.isFailure()  ||  !cmmJetHitsTES ) {
    m_log << MSG::DEBUG << "No CMM-Jet Hits container found"<< endreq; 
  }
  
  //Retrieve CMM RoIs from SG
  const LVL1::CMMRoI* cmmRoiTES = 0;
  sc = m_storeGate->retrieve( cmmRoiTES, m_cmmRoiLocation);
  if( sc.isFailure()  ||  !cmmRoiTES  || ( !cmmRoiTES->jetEtRoiWord() &&
         !cmmRoiTES->energyRoiWord0() && !cmmRoiTES->energyRoiWord1() &&
                                         !cmmRoiTES->energyRoiWord2())) {
    m_log << MSG::DEBUG << "No DAQ CMM RoIs found" << endreq; 
  }

  //Retrieve JEM Et Sums from SG
  const JemEtSumsCollection* jemEtSumsTES = 0;
  sc = m_storeGate->retrieve( jemEtSumsTES, m_jemEtSumsLocation);
  if( sc.isFailure()  ||  !jemEtSumsTES ) {
    m_log << MSG::DEBUG << "No JEM Et Sums container found"<< endreq;
  }

  //Retrieve CMM Et Sums from SG
  const CmmEtSumsCollection* cmmEtSumsTES = 0;
  sc = m_storeGate->retrieve( cmmEtSumsTES, m_cmmEtSumsLocation);
  if( sc.isFailure()  ||  !cmmEtSumsTES ) {
    m_log << MSG::DEBUG << "No CMM-Energy Et Sums container found"<< endreq;
  }

  // Maps to simplify comparisons
  
  JetElementMap jeMap;
  JetElementMap ovMap;
  JemRoiMap     jrMap;
  JemHitsMap    jhMap;
  CmmJetHitsMap cmMap;
  JemEtSumsMap  jsMap;
  CmmEtSumsMap  csMap;
  setupMap(jetElementTES,   jeMap);
  setupMap(jetElementOvTES, ovMap);
  setupMap(jemRoiTES,       jrMap);
  setupMap(jemHitsTES,      jhMap);
  setupMap(cmmJetHitsTES,   cmMap);
  setupMap(jemEtSumsTES,    jsMap);
  setupMap(cmmEtSumsTES,    csMap);

  // Vectors for error overview bits;
  const int nCrates = 2;
  const int nJEMs   = 16;
  const int nCMMs   = 2;
  const int vecsizeJem = 2 * nCrates * nJEMs;
  const int vecsizeCmm = 2 * nCrates * nCMMs;
  ErrorVector errorsJEM(vecsizeJem);
  ErrorVector errorsCMM(vecsizeCmm);

  if (m_compareWithSim) {

    // Compare Jet Elements simulated from Trigger Towers with Jet Elements
    // from data

    JetElementCollection* jetElementSIM = 0;
    if (triggerTowerTES) {
      jetElementSIM = new JetElementCollection;
      simulate(triggerTowerTES, jetElementSIM);
    }
    JetElementMap jeSimMap;
    setupMap(jetElementSIM, jeSimMap);
    bool overlap = false;
    compare(jeSimMap, jeMap, errorsJEM, overlap);
    if (jetElementOvTES) {
      overlap = true;
      compare(jeSimMap, ovMap, errorsJEM, overlap);
    }
    jeSimMap.clear();
    delete jetElementSIM;

    // Compare RoIs simulated from Jet Elements with JEM RoIs from data

    JemRoiCollection* jemRoiSIM = 0;
    if (jetElementTES || jetElementOvTES) {
      jemRoiSIM = new JemRoiCollection;
      simulate(jetElementTES, jetElementOvTES, jemRoiSIM);
    }
    JemRoiMap jrSimMap;
    setupMap(jemRoiSIM, jrSimMap);
    compare(jrSimMap, jrMap, rodTES, errorsJEM);
    jrSimMap.clear();
    delete jemRoiSIM;

    // Compare JEM Hits simulated from JEM RoIs with JEM Hits from data

    JemHitsCollection* jemHitsSIM = 0;
    if (jemRoiTES) {
      jemHitsSIM = new JemHitsCollection;
      simulate(jemRoiTES, jemHitsSIM);
    }
    JemHitsMap jhSimMap;
    setupMap(jemHitsSIM, jhSimMap);
    compare(jhSimMap, jhMap, errorsJEM);
    jhSimMap.clear();
    delete jemHitsSIM;

  }

  // Compare JEM hits with CMM Hits from data

  compare(jhMap, cmMap, errorsJEM, errorsCMM);

  if (m_compareWithSim) {

    // Compare Local sums simulated from CMM Hits with Local sums from data

    CmmJetHitsCollection* cmmLocalSIM = 0;
    if (cmmJetHitsTES) {
      cmmLocalSIM = new CmmJetHitsCollection;
      simulate(cmmJetHitsTES, cmmLocalSIM, LVL1::CMMJetHits::LOCAL_MAIN);
    }
    CmmJetHitsMap cmmLocalSimMap;
    setupMap(cmmLocalSIM, cmmLocalSimMap);
    compare(cmmLocalSimMap, cmMap, errorsCMM, LVL1::CMMJetHits::LOCAL_MAIN);
    cmmLocalSimMap.clear();
    delete cmmLocalSIM;

  }

  // Compare Local sums with Remote sums from data

  compare(cmMap, cmMap, errorsCMM, LVL1::CMMJetHits::REMOTE_MAIN);

  if (m_compareWithSim) {

    // Compare Total sums simulated from Remote sums with Total sums from data

    CmmJetHitsCollection* cmmTotalSIM = 0;
    if (cmmJetHitsTES) {
      cmmTotalSIM = new CmmJetHitsCollection;
      simulate(cmmJetHitsTES, cmmTotalSIM, LVL1::CMMJetHits::TOTAL_MAIN);
    }
    CmmJetHitsMap cmmTotalSimMap;
    setupMap(cmmTotalSIM, cmmTotalSimMap);
    compare(cmmTotalSimMap, cmMap, errorsCMM, LVL1::CMMJetHits::TOTAL_MAIN);
    cmmTotalSimMap.clear();
    delete cmmTotalSIM;

    // Compare JetEt Map simulated from Total sums with JetEt Map from data

    CmmJetHitsCollection* cmmJetEtSIM = 0;
    if (cmmJetHitsTES) {
      cmmJetEtSIM = new CmmJetHitsCollection;
      simulate(cmmJetHitsTES, cmmJetEtSIM, LVL1::CMMJetHits::ET_MAP);
    }
    CmmJetHitsMap cmmJetEtSimMap;
    setupMap(cmmJetEtSIM, cmmJetEtSimMap);
    compare(cmmJetEtSimMap, cmMap, errorsCMM, LVL1::CMMJetHits::ET_MAP);
    cmmJetEtSimMap.clear();
    delete cmmJetEtSIM;

  }

  // Compare JetEt Map with JetEt RoI from data

  compare(cmMap, cmmRoiTES, errorsCMM);

  if (m_compareWithSim) {

    // Compare JEMEtSums simulated from JetElements with JEMEtSums from data

    JemEtSumsCollection* jemEtSumsSIM = 0;
    if (jemEtSumsTES) {
      jemEtSumsSIM = new JemEtSumsCollection;
      simulate(jetElementTES, jemEtSumsSIM);
    }
    JemEtSumsMap jemEtSumsSimMap;
    setupMap(jemEtSumsSIM, jemEtSumsSimMap);
    compare(jemEtSumsSimMap, jsMap, errorsJEM);
    jemEtSumsSimMap.clear();
    delete jemEtSumsSIM;

  }

  // Compare JEMEtSums with CMMEtSums from data

  compare(jsMap, csMap, errorsJEM, errorsCMM);

  if (m_compareWithSim) {

    // Compare Local sums simulated from CMMEtSums with Local sums from data

    CmmEtSumsCollection* cmmEtLocalSIM = 0;
    if (cmmEtSumsTES) {
      cmmEtLocalSIM = new CmmEtSumsCollection;
      simulate(cmmEtSumsTES, cmmEtLocalSIM, LVL1::CMMEtSums::LOCAL);
    }
    CmmEtSumsMap cmmEtLocalSimMap;
    setupMap(cmmEtLocalSIM, cmmEtLocalSimMap);
    compare(cmmEtLocalSimMap, csMap, errorsCMM, LVL1::CMMEtSums::LOCAL);
    cmmEtLocalSimMap.clear();
    delete cmmEtLocalSIM;

  }

  // Compare Local Energy sums with Remote sums from data

  compare(csMap, csMap, errorsCMM, LVL1::CMMEtSums::REMOTE);

  if (m_compareWithSim) {

    // Compare Total sums simulated from Remote sums with Total sums from data

    CmmEtSumsCollection* cmmEtTotalSIM = 0;
    if (cmmEtSumsTES) {
      cmmEtTotalSIM = new CmmEtSumsCollection;
      simulate(cmmEtSumsTES, cmmEtTotalSIM, LVL1::CMMEtSums::TOTAL);
    }
    CmmEtSumsMap cmmEtTotalSimMap;
    setupMap(cmmEtTotalSIM, cmmEtTotalSimMap);
    compare(cmmEtTotalSimMap, csMap, errorsCMM, LVL1::CMMEtSums::TOTAL);
    cmmEtTotalSimMap.clear();
    delete cmmEtTotalSIM;

    // Compare Et Maps (sumEt/missingEt) simulated from Total sums
    // with Et Maps from data

    CmmEtSumsCollection* cmmSumEtSIM = 0;
    if (cmmEtSumsTES) {
      cmmSumEtSIM = new CmmEtSumsCollection;
      simulate(cmmEtSumsTES, cmmSumEtSIM, LVL1::CMMEtSums::SUM_ET_MAP);
    }
    CmmEtSumsMap cmmSumEtSimMap;
    setupMap(cmmSumEtSIM, cmmSumEtSimMap);
    compare(cmmSumEtSimMap, csMap, errorsCMM, LVL1::CMMEtSums::SUM_ET_MAP);
    cmmSumEtSimMap.clear();
    delete cmmSumEtSIM;

  }

  // Compare Total Energy sums and Et Maps with Energy RoIs from data

  compare(csMap, cmmRoiTES, errorsCMM);

  // Update error summary plots

  ErrorVector crateErr(nCrates);
  const int jemBins = nCrates * nJEMs;
  const int cmmBins = nCrates * nCMMs;
  m_eventNumber = -1;
  for (int err = 0; err < NumberOfSummaryBins; ++err) {
    int error = 0;
    for (int loc = 0; loc < jemBins; ++loc) {
      if ((errorsJEM[loc] >> err) & 0x1) {
        m_h_JEPeqSIM->Fill(loc, err, 1.);
	//temporary duplicates
        m_h_JEPeqSIM2->Fill(loc, err, 1.);
        m_h_JEPeqSIM3->Fill(loc, err, 1.);
      }
      if ((errorsJEM[loc + jemBins] >> err) & 0x1) {
        m_h_JEPneSIM->Fill(loc, err, 1.);
        m_h_JEPneSIM2->Fill(loc, err, 1.);
        m_h_JEPneSIM3->Fill(loc, err, 1.);
	error = 1;
	crateErr[loc/nJEMs] |= (1 << err);
	fillEventSample(err, loc, true);
      }
      if (loc < cmmBins) {
        if ((errorsCMM[loc] >> err) & 0x1) {
          m_h_JEPeqSIM->Fill(loc+jemBins, err, 1.);
          m_h_JEPeqSIM2->Fill(loc+jemBins, err, 1.);
          m_h_JEPeqSIM3->Fill(loc+jemBins, err, 1.);
        }
        if ((errorsCMM[loc + cmmBins] >> err) & 0x1) {
          m_h_JEPneSIM->Fill(loc+jemBins, err, 1.);
          m_h_JEPneSIM2->Fill(loc+jemBins, err, 1.);
          m_h_JEPneSIM3->Fill(loc+jemBins, err, 1.);
	  error = 1;
	  crateErr[loc/nCMMs] |= (1 << err);
	  fillEventSample(err, loc, false);
        }
      }
    }
    m_h_JEPneSIMSummary->Fill(err, error);
  }

  // Save error vector for global summary

  ErrorVector* save = new ErrorVector(crateErr);
  sc = m_storeGate->record(save, "L1CaloJEMMismatchVector");
  if (sc != StatusCode::SUCCESS) {
    m_log << MSG::ERROR << "Error recording JEM mismatch vector in TES "
          << endreq;
    return sc;
  }

  m_log << MSG::DEBUG << "Leaving fillHistograms" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode JEPSimBSMon::procHistograms(bool isEndOfEventsBlock,
                                  bool isEndOfLumiBlock, bool isEndOfRun)
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "procHistograms entered" << endreq;

  if (isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun) {
  }

  return StatusCode::SUCCESS;
}

TH1F* JEPSimBSMon::book1F(const std::string& name,
                          const std::string& title,
                          int nx, double xmin, double xmax)
{
  TH1F *hist = new TH1F(name.c_str(), title.c_str(), nx, xmin, xmax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    m_log << MSG::WARNING << "Could not register histogram : " 
  	  << name << endreq;
  }
  hist->SetStats(kFALSE);
  
  return hist;
}

TH2F* JEPSimBSMon::book2F(const std::string& name,
                          const std::string& title,
                          int nx, double xmin, double xmax,  
	                  int ny, double ymin, double ymax)
{		
  TH2F *hist = new TH2F(name.c_str(), title.c_str(), nx, xmin, xmax,
                                                     ny, ymin, ymax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    m_log << MSG::WARNING << "Could not register histogram : " 
	  << name << endreq;
  }
  hist->SetOption("colz");
  hist->SetStats(kFALSE);
  
  return hist;
}

TH2F* JEPSimBSMon::book2F(const std::string& name,
                          const std::string& title,
                          int nx, const double* xbins,
	                  int ny, double ymin, double ymax)
{		
  TH2F *hist = new TH2F(name.c_str(), title.c_str(), nx, xbins,
                                                     ny, ymin, ymax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    m_log << MSG::WARNING << "Could not register histogram : " 
	  << name << endreq;
  }
  hist->SetOption("colz");
  hist->SetStats(kFALSE);
  
  return hist;
}

TH2I* JEPSimBSMon::book2I(const std::string& name,
                          const std::string& title,
			  int nx, double xmin, double xmax,
			  int ny, double ymin, double ymax)
{
  TH2I *hist = new TH2I(name.c_str(), title.c_str(), nx, xmin, xmax,
                                                     ny, ymin, ymax);

  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    m_log << MSG::WARNING << "Could not register histogram : " 
	  << name << endreq;
  }
  hist->SetOption("text");
  hist->SetStats(kFALSE);
  
  return hist;
}

TH2F* JEPSimBSMon::bookEtaPhi(const std::string& name,
                              const std::string& title, bool isRoi)
{
  // We use phi range 0-32 so tick marks correspond to the bins
  TH2F* hist = 0;
  const double phiBin     = M_PI/16.;
  const double halfPhiBin = M_PI/32.;
  const int nxbins = 32;
  const double xbins[nxbins+1] = {-4.9,-3.2,-2.9,-2.7,-2.4,-2.2,-2.0,-1.8,-1.6,
                                  -1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,
				  0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,
				  2.7,2.9,3.2,4.9};
  const double xbinsRoi[nxbins+1] = {-4.0,-3.05,-2.8,-2.55,-2.3,-2.1,-1.9,-1.7,
                                     -1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,
				     0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,
				     2.1,2.3,2.55,2.8,3.05,4.0,4.95};
  std::string newTitle = title + ";eta";
  hist = (isRoi) ? book2F(name, newTitle, nxbins, xbinsRoi, 32, 0., 32.)
                 : book2F(name, newTitle, nxbins, xbins,    32, 0., 32.);
  for (int chan = 0; chan < 32; chan += 2 ) {
    const double rad = (isRoi) ? (chan + 1)*phiBin : chan*phiBin + halfPhiBin;
    std::ostringstream cnum;
    cnum << chan << "/"
         << std::setiosflags(std::ios::fixed | std::ios::showpoint)
	 << std::setprecision(2) << rad;
    hist->GetYaxis()->SetBinLabel(chan+1, cnum.str().c_str());
  }
  hist->GetYaxis()->SetBinLabel(32, "phi");
  return hist;
}

//  Compare Simulated JetElements with data

void JEPSimBSMon::compare(const JetElementMap& jeSimMap,
                          const JetElementMap& jeMap, ErrorVector& errors,
			  bool overlap)
{
  m_log << MSG::DEBUG << "Compare Simulated JetElements with data" << endreq;

  const int maxKey = 0x7fffffff;
  JetElementMap::const_iterator simMapIter    = jeSimMap.begin();
  JetElementMap::const_iterator simMapIterEnd = jeSimMap.end();
  JetElementMap::const_iterator datMapIter    = jeMap.begin();
  JetElementMap::const_iterator datMapIterEnd = jeMap.end();

  while (simMapIter != simMapIterEnd || datMapIter != datMapIterEnd) {

    int simKey = maxKey;
    int datKey = maxKey;
    int simEm  = 0;
    int simHad = 0;
    int datEm  = 0;
    int datHad = 0;
    double eta = 0.;
    double phi = 0.;

    if (simMapIter != simMapIterEnd) simKey = simMapIter->first;
    if (datMapIter != datMapIterEnd) datKey = datMapIter->first;

    if ((datMapIter == datMapIterEnd) || (datKey > simKey)) {

      // Simulated JetElement but no data JetElement

      const LVL1::JetElement* je = simMapIter->second;
      simEm  = je->emEnergy();
      simHad = je->hadEnergy();
      eta = je->eta();
      phi = je->phi();
      ++simMapIter;

    } else if ((simMapIter == simMapIterEnd) || (simKey > datKey)) {

      // Data JetElement but no simulated JetElement

      const LVL1::JetElement* je = datMapIter->second;
      datEm  = je->emEnergy();
      datHad = je->hadEnergy();
      eta = je->eta();
      phi = je->phi();
      ++datMapIter;

    } else {

      // Have both

      const LVL1::JetElement* jeS = simMapIter->second;
      const LVL1::JetElement* jeD = datMapIter->second;
      simEm  = jeS->emEnergy();
      simHad = jeS->hadEnergy();
      datEm  = jeD->emEnergy();
      datHad = jeD->hadEnergy();
      eta = jeD->eta();
      phi = jeD->phi();
      ++simMapIter;
      ++datMapIter;
    }

    if (!simEm && !simHad && !datEm && !datHad) continue;
    
    //  Fill in error vectors

    const LVL1::Coordinate coord(phi, eta);
    LVL1::CoordToHardware converter;
    const int crate = (overlap) ? converter.jepCrateOverlap(coord)
                                : converter.jepCrate(coord);
    const int jem   = (overlap) ? converter.jepModuleOverlap(coord)
                                : converter.jepModule(coord);
    if (crate > 1 || jem > 15) continue;
    const int loc   = crate * 16 + jem;
    const int jemBins = 2 * 16;
    const int bitEm  = (1 << EMElementMismatch);
    const int bitHad = (1 << HadElementMismatch);
    if (simEm && simEm == datEm)    errors[loc] |= bitEm;
    if (simHad && simHad == datHad) errors[loc] |= bitHad;
    if (simEm != datEm)     errors[loc+jemBins] |= bitEm;
    if (simHad != datHad)   errors[loc+jemBins] |= bitHad;
    if (simEm != datEm || simHad != datHad) {
      m_log << MSG::VERBOSE << "JE mismatch, EM data/sim: " << datEm << "/"
            << simEm << " Had data/sim: " << datHad << "/" << simHad
	    << " crate/jem: " << crate << "/" << jem
	    << " eta/phi: " << eta << "/" << phi
	    << endreq;
    }
    // Fill two phi bins for FCAL
    const bool fcal = (eta < -3.2 || eta > 3.2);
    const double phi1 = ((fcal) ? phi - M_PI/64. : phi) * m_phiScale;
    const double phi2 = (phi + M_PI/64.) * m_phiScale;
    if (overlap) {
      m_h_EMEleOvSIMeqDAT->Fill(eta, phi1, simEm && simEm == datEm);
      m_h_EMEleOvSIMneDAT->Fill(eta, phi1, simEm && datEm && simEm != datEm);
      m_h_EMEleOvSIMnoDAT->Fill(eta, phi1, simEm && !datEm);
      m_h_EMEleOvDATnoSIM->Fill(eta, phi1, datEm && !simEm);
      m_h_HadEleOvSIMeqDAT->Fill(eta, phi1, simHad && simHad == datHad);
      m_h_HadEleOvSIMneDAT->Fill(eta, phi1, simHad && datHad && simHad != datHad);
      m_h_HadEleOvSIMnoDAT->Fill(eta, phi1, simHad && !datHad);
      m_h_HadEleOvDATnoSIM->Fill(eta, phi1, datHad && !simHad);
      if (fcal) {
        m_h_EMEleOvSIMeqDAT->Fill(eta, phi2, simEm && simEm == datEm);
        m_h_EMEleOvSIMneDAT->Fill(eta, phi2, simEm && datEm && simEm != datEm);
        m_h_EMEleOvSIMnoDAT->Fill(eta, phi2, simEm && !datEm);
        m_h_EMEleOvDATnoSIM->Fill(eta, phi2, datEm && !simEm);
        m_h_HadEleOvSIMeqDAT->Fill(eta, phi2, simHad && simHad == datHad);
        m_h_HadEleOvSIMneDAT->Fill(eta, phi2, simHad && datHad && simHad != datHad);
        m_h_HadEleOvSIMnoDAT->Fill(eta, phi2, simHad && !datHad);
        m_h_HadEleOvDATnoSIM->Fill(eta, phi2, datHad && !simHad);
      }
    } else {
      m_h_EMEleSIMeqDAT->Fill(eta, phi1, simEm && simEm == datEm);
      m_h_EMEleSIMneDAT->Fill(eta, phi1, simEm && datEm && simEm != datEm);
      m_h_EMEleSIMnoDAT->Fill(eta, phi1, simEm && !datEm);
      m_h_EMEleDATnoSIM->Fill(eta, phi1, datEm && !simEm);
      m_h_HadEleSIMeqDAT->Fill(eta, phi1, simHad && simHad == datHad);
      m_h_HadEleSIMneDAT->Fill(eta, phi1, simHad && datHad && simHad != datHad);
      m_h_HadEleSIMnoDAT->Fill(eta, phi1, simHad && !datHad);
      m_h_HadEleDATnoSIM->Fill(eta, phi1, datHad && !simHad);
      if (fcal) {
        m_h_EMEleSIMeqDAT->Fill(eta, phi2, simEm && simEm == datEm);
        m_h_EMEleSIMneDAT->Fill(eta, phi2, simEm && datEm && simEm != datEm);
        m_h_EMEleSIMnoDAT->Fill(eta, phi2, simEm && !datEm);
        m_h_EMEleDATnoSIM->Fill(eta, phi2, datEm && !simEm);
        m_h_HadEleSIMeqDAT->Fill(eta, phi2, simHad && simHad == datHad);
        m_h_HadEleSIMneDAT->Fill(eta, phi2, simHad && datHad && simHad != datHad);
        m_h_HadEleSIMnoDAT->Fill(eta, phi2, simHad && !datHad);
        m_h_HadEleDATnoSIM->Fill(eta, phi2, datHad && !simHad);
      }
    }
  }
}

//  Compare Simulated RoIs with data

void JEPSimBSMon::compare(const JemRoiMap& roiSimMap,
                          const JemRoiMap& roiMap,
			  const RodHeaderCollection* rods,
                                ErrorVector& errors)
{
  m_log << MSG::DEBUG << "Compare Simulated RoIs with data" << endreq;

  const int nCrates = 2;
  const int nJEMs = 16;
  const int maxKey = 0xffff;
  LVL1::JEPRoIDecoder decoder;
  std::vector<int> limitedRoi(nCrates);
  JemRoiMap::const_iterator simMapIter    = roiSimMap.begin();
  JemRoiMap::const_iterator simMapIterEnd = roiSimMap.end();
  JemRoiMap::const_iterator datMapIter    = roiMap.begin();
  JemRoiMap::const_iterator datMapIterEnd = roiMap.end();

  while (simMapIter != simMapIterEnd || datMapIter != datMapIterEnd) {

    int simKey = maxKey;
    int datKey = maxKey;
    unsigned int simHits = 0;
    unsigned int datHits = 0;
    const LVL1::JEMRoI* roi = 0;

    if (simMapIter != simMapIterEnd) simKey = simMapIter->first;
    if (datMapIter != datMapIterEnd) datKey = datMapIter->first;

    if ((datMapIter == datMapIterEnd) || (datKey > simKey)) {

      // Simulated RoI but no data RoI

      roi = simMapIter->second;
      simHits = roi->hits();
      ++simMapIter;
      m_log << MSG::VERBOSE
            << "Sim  RoI crate/jem/frame/loc/fwd/error/hits: "
            << roi->crate() << "/" << roi->jem() << "/" << roi->frame() << "/"
	    << roi->location() << "/" << roi->forward() << "/"
	    << roi->error() << "/" << MSG::hex << roi->hits() << MSG::dec
	    << endreq;

    } else if ((simMapIter == simMapIterEnd) || (simKey > datKey)) {

      // Data RoI but no simulated RoI

      roi = datMapIter->second;
      datHits = roi->hits();
      ++datMapIter;
      m_log << MSG::VERBOSE
            << "Data RoI crate/jem/frame/loc/fwd/error/hits: "
            << roi->crate() << "/" << roi->jem() << "/" << roi->frame() << "/"
	    << roi->location() << "/" << roi->forward() << "/"
	    << roi->error() << "/" << MSG::hex << roi->hits() << MSG::dec
	    << endreq;

    } else {

      // Have both

      const LVL1::JEMRoI* roiS = simMapIter->second;
      roi     = datMapIter->second;
      simHits = roiS->hits();
      datHits = roi->hits();
      ++simMapIter;
      ++datMapIter;
      m_log << MSG::VERBOSE
            << "Sim  RoI crate/jem/frame/loc/fwd/error/hits: "
            << roiS->crate() << "/" << roiS->jem() << "/" << roiS->frame() << "/"
	    << roiS->location() << "/" << roiS->forward() << "/"
	    << roiS->error() << "/" << MSG::hex << roiS->hits() << MSG::dec
	    << endreq;
      m_log << MSG::VERBOSE
            << "Data RoI crate/jem/frame/loc/fwd/error/hits: "
            << roi->crate() << "/" << roi->jem() << "/" << roi->frame() << "/"
	    << roi->location() << "/" << roi->forward() << "/"
	    << roi->error() << "/" << MSG::hex << roi->hits() << MSG::dec
	    << endreq;
    }

    if (!simHits && !datHits) continue;

    //  Check LimitedRoISet bit

    const int crate = roi->crate();
    if (!datHits) {
      if (rods) {
        RodHeaderCollection::const_iterator rodIter  = rods->begin();
	RodHeaderCollection::const_iterator rodIterE = rods->end();
	for (; rodIter != rodIterE; ++rodIter) {
	  LVL1::RODHeader* rod = *rodIter;
	  const int rodCrate = rod->crate() - 12;
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

    const int jem   = roi->jem();
    const int frame = roi->frame();
    const int local = roi->location();
    const int forward = roi->forward();
    const int locX  = crate * nJEMs + jem;
    const int locY  = frame * 4 + local;
    const int jemBins = nCrates * nJEMs;
    const int bit = (1 << RoIMismatch);
    const LVL1::CoordinateRange coord(decoder.coordinate(roi->roiWord()));
    double eta = coord.eta();
    // Distinguish right forward columns 3 and 4 for checking purposes
    if (forward && eta > 0.0 && frame > 3) eta = (local%2) ? 4.05 : 3.9;
    const double phi = coord.phi() * m_phiScale - 0.5;

    TH2F* hist1 = 0;
    TH2F* hist2 = 0;
    if (simHits == datHits) {
      errors[locX] |= bit;
      hist1 = m_h_RoISIMeqDAT;
      hist2 = m_h_RoIEtaPhiSIMeqDAT;
    } else {
      errors[locX+jemBins] |= bit;
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
    if (hist2) hist2->Fill(eta, phi);

    const int nThresh = (forward) ? 4 : 8;
    const int offset  = (forward) ? 8 : 0;
    for (int thr = 0; thr < nThresh; ++thr) {
      const int thrDat = (datHits >> thr) & 0x1;
      const int thrSim = (simHits >> thr) & 0x1;
      if (thrDat || thrSim) {
        if (thrDat == thrSim) m_h_RoIThreshSIMeqDAT->Fill(locX, thr+offset);
        else                  m_h_RoIThreshSIMneDAT->Fill(locX, thr+offset);
      }
    }

    if (simHits != datHits) {
      m_log << MSG::VERBOSE << " RoI Mismatch Crate/JEM DataHits/SimHits: ";
      m_log << MSG::VERBOSE << crate << "/" << jem << " ";
      for (int i = 7; i >= 0; --i) {
        int hit = (datHits >> i) & 0x1;
        m_log << MSG::VERBOSE << hit;
      }
      m_log << MSG::VERBOSE << "/";
      for (int i = 7; i >= 0; --i) {
        int hit = (simHits >> i) & 0x1;
        m_log << MSG::VERBOSE << hit;
      }
      m_log << MSG::VERBOSE << endreq;
    }
  }
}

//  Compare simulated JEM Hits with data

void JEPSimBSMon::compare(const JemHitsMap& jemSimMap,
                          const JemHitsMap& jemMap,
                                ErrorVector& errors)
{
  m_log << MSG::DEBUG << "Compare simulated JEM Hits with data" << endreq;

  const int maxKey = 0x7fffffff;
  JemHitsMap::const_iterator simMapIter    = jemSimMap.begin();
  JemHitsMap::const_iterator simMapIterEnd = jemSimMap.end();
  JemHitsMap::const_iterator datMapIter    = jemMap.begin();
  JemHitsMap::const_iterator datMapIterEnd = jemMap.end();

  while (simMapIter != simMapIterEnd || datMapIter != datMapIterEnd) {

    int simKey = maxKey;
    int datKey = maxKey;
    unsigned int simHits = 0;
    unsigned int datHits = 0;
    int crate = 0;
    int jem   = 0;

    if (simMapIter != simMapIterEnd) simKey = simMapIter->first;
    if (datMapIter != datMapIterEnd) datKey = datMapIter->first;

    if ((datMapIter == datMapIterEnd) || (datKey > simKey)) {

      // Simulation Hits but no data Hits

      const LVL1::JEMHits* simh = simMapIter->second;
      simHits = simh->JetHits();
      crate   = simh->crate();
      jem     = simh->module();
      ++simMapIter;

    } else if ((simMapIter == simMapIterEnd) || (simKey > datKey)) {

      // Data Hits but no simulation Hits

      const LVL1::JEMHits* dath = datMapIter->second;
      datHits = dath->JetHits();
      crate   = dath->crate();
      jem     = dath->module();
      ++datMapIter;

    } else {

      // Have both

      const LVL1::JEMHits* simh = simMapIter->second;
      const LVL1::JEMHits* dath = datMapIter->second;
      simHits = simh->JetHits();
      datHits = dath->JetHits();
      crate   = dath->crate();
      jem     = dath->module();
      ++simMapIter;
      ++datMapIter;
    }

    if (!simHits && !datHits) continue;
    
    //  Fill in error plots

    const int loc = crate * 16 + jem;
    const int jemBins = 2 * 16;
    const int bit = (1 << JEMHitsMismatch);
    TH2F* hist = 0;
    if (simHits == datHits) {
      errors[loc] |= bit;
      hist = m_h_JEMHitsSIMeqDAT;
    } else {
      errors[loc+jemBins] |= bit;
      if (simHits && datHits) hist = m_h_JEMHitsSIMneDAT;
      else if (simHits)       hist = m_h_JEMHitsSIMnoDAT;
      else                    hist = m_h_JEMHitsDATnoSIM;
    }
    if (hist) hist->Fill(jem, crate);

    int nThresh = 8;
    int thrLen  = 3;
    unsigned int thrMask = 0x7;
    if (jem == 0 || jem == 7 || jem == 8 || jem == 15) {
      nThresh = 12;
      thrLen  = 2;
      thrMask = 0x3;
    }
    for (int thr = 0; thr < nThresh; ++thr) {
      const int shift = thrLen*thr;
      const unsigned int d0 = (datHits >> shift) & thrMask;
      const unsigned int s0 = (simHits >> shift) & thrMask;
      if (d0 || s0) {
        if (d0 == s0) m_h_JEMHitsThreshSIMeqDAT->Fill(loc, thr);
	else          m_h_JEMHitsThreshSIMneDAT->Fill(loc, thr);
      }
    }
  }
}

//  Compare JEM Hits and CMM Hits

void JEPSimBSMon::compare(const JemHitsMap& jemMap,
                          const CmmJetHitsMap& cmmMap,
                                ErrorVector& errorsJEM,
				ErrorVector& errorsCMM)
{
  m_log << MSG::DEBUG << "Compare JEM Hits and CMM Hits" << endreq;

  const int maxKey = 0x7fffffff;
  JemHitsMap::const_iterator    jemMapIter    = jemMap.begin();
  JemHitsMap::const_iterator    jemMapIterEnd = jemMap.end();
  CmmJetHitsMap::const_iterator cmmMapIter    = cmmMap.begin();
  CmmJetHitsMap::const_iterator cmmMapIterEnd = cmmMap.end();

  while (jemMapIter != jemMapIterEnd || cmmMapIter != cmmMapIterEnd) {

    int jemKey = maxKey;
    int cmmKey = maxKey;
    unsigned int jemHits = 0;
    unsigned int cmmHits = 0;
    int crate  = 0;
    int jem    = 0;

    if (jemMapIter != jemMapIterEnd) jemKey = jemMapIter->first;
    if (cmmMapIter != cmmMapIterEnd) cmmKey = cmmMapIter->first;

    if ((cmmMapIter == cmmMapIterEnd) || (cmmKey > jemKey)) {

      // JEM Hits but no CMM Hits

      const LVL1::JEMHits* jemh = jemMapIter->second;
      jemHits = jemh->JetHits();
      crate   = jemh->crate();
      jem     = jemh->module();
      ++jemMapIter;

    } else if ((jemMapIter == jemMapIterEnd) || (jemKey > cmmKey)) {

      // CMM Hits but no JEM Hits

      const LVL1::CMMJetHits* cmmh = cmmMapIter->second;
      cmmHits = cmmh->Hits();
      crate   = cmmh->crate();
      jem     = cmmh->dataID();
      ++cmmMapIter;
      if (jem > 15) continue;

    } else {

      // Have both

      const LVL1::JEMHits*    jemh = jemMapIter->second;
      const LVL1::CMMJetHits* cmmh = cmmMapIter->second;
      jemHits = jemh->JetHits();
      cmmHits = cmmh->Hits();
      crate   = jemh->crate();
      jem     = jemh->module();
      ++jemMapIter;
      ++cmmMapIter;
    }

    if (!jemHits && !cmmHits) continue;
    
    //  Fill in error plots

    const int loc  = crate * 16 + jem;
    const int loc2 = crate * 2 + 1; // CMM-Jet is right ==> cmm 1
    const int jemBins = 2 * 16;
    const int cmmBins = 2 * 2;
    const int bit = (1 << CMMJetHitsMismatch);
    TH2F* hist = 0;
    if (jemHits == cmmHits) {
      errorsJEM[loc]  |= bit;
      errorsCMM[loc2] |= bit;
      hist = m_h_CMMHitsSIMeqDAT;
    } else {
      errorsJEM[loc+jemBins]  |= bit;
      errorsCMM[loc2+cmmBins] |= bit;
      if (jemHits && cmmHits) hist = m_h_CMMHitsSIMneDAT;
      else if (jemHits)       hist = m_h_CMMHitsSIMnoDAT;
      else                    hist = m_h_CMMHitsDATnoSIM;
    }
    if (hist) hist->Fill(jem, crate);

    int nThresh = 8;
    int thrLen  = 3;
    unsigned int thrMask = 0x7;
    if (jem == 0 || jem == 7 || jem == 8 || jem == 15) {
      nThresh = 12;
      thrLen  = 2;
      thrMask = 0x3;
    }
    for (int thr = 0; thr < nThresh; ++thr) {
      const int shift = thrLen*thr;
      const unsigned int d0 = (cmmHits >> shift) & thrMask;
      const unsigned int s0 = (jemHits >> shift) & thrMask;
      if (d0 || s0) {
        if (d0 == s0) m_h_CMMHitsThreshSIMeqDAT->Fill(loc, thr);
	else          m_h_CMMHitsThreshSIMneDAT->Fill(loc, thr);
      }
    }
  }
}

//  Compare Simulated CMM Hit Sums and Data CMM Hit Sums

void JEPSimBSMon::compare(const CmmJetHitsMap& cmmSimMap,
                          const CmmJetHitsMap& cmmMap,
                                ErrorVector& errors, int selection)
{
  m_log << MSG::DEBUG << "Compare Simulated CMM Hit Sums and Data CMM Hit Sums"
        << endreq;

  const bool local  = (selection == LVL1::CMMJetHits::LOCAL_MAIN);
  const bool remote = (selection == LVL1::CMMJetHits::REMOTE_MAIN);
  const bool total  = (selection == LVL1::CMMJetHits::TOTAL_MAIN);
  const bool etmap  = (selection == LVL1::CMMJetHits::ET_MAP);
  if (!local && !remote && !total && !etmap) return;
  unsigned int hitsSimMain = 0;
  unsigned int hitsSimFwd  = 0;
  unsigned int hitsDatMain = 0;
  unsigned int hitsDatFwd  = 0;
  const int maxKey = 0x7fffffff;
  CmmJetHitsMap::const_iterator cmmSimMapIter    = cmmSimMap.begin();
  CmmJetHitsMap::const_iterator cmmSimMapIterEnd = cmmSimMap.end();
  CmmJetHitsMap::const_iterator cmmMapIter       = cmmMap.begin();
  CmmJetHitsMap::const_iterator cmmMapIterEnd    = cmmMap.end();

  while (cmmSimMapIter != cmmSimMapIterEnd || cmmMapIter != cmmMapIterEnd) {

    int cmmSimKey = maxKey;
    int cmmKey    = maxKey;
    unsigned int cmmSimHits = 0;
    unsigned int cmmHits = 0;
    int crate  = 0;
    int dataId = 0;

    if (cmmSimMapIter != cmmSimMapIterEnd) cmmSimKey = cmmSimMapIter->first;
    if (cmmMapIter    != cmmMapIterEnd)    cmmKey    = cmmMapIter->first;

    if ((cmmMapIter == cmmMapIterEnd) || (cmmKey > cmmSimKey)) {

      // Sim CMM Hits but no Data CMM Hits

      const LVL1::CMMJetHits* cmmS = cmmSimMapIter->second;
      ++cmmSimMapIter;
      dataId = cmmS->dataID();
      if (local  && dataId != LVL1::CMMJetHits::LOCAL_MAIN &&
                    dataId != LVL1::CMMJetHits::LOCAL_FORWARD) continue;
      if (remote && dataId != LVL1::CMMJetHits::LOCAL_MAIN &&
                    dataId != LVL1::CMMJetHits::LOCAL_FORWARD) continue;
      if (total  && dataId != LVL1::CMMJetHits::TOTAL_MAIN &&
                    dataId != LVL1::CMMJetHits::TOTAL_FORWARD) continue;
      if (etmap  && dataId != LVL1::CMMJetHits::ET_MAP)        continue;
      cmmSimHits = cmmS->Hits();
      crate      = cmmS->crate();

    } else if ((cmmSimMapIter == cmmSimMapIterEnd) || (cmmSimKey > cmmKey)) {

      // Data CMM Hits but no Sim CMM Hits

      const LVL1::CMMJetHits* cmmD = cmmMapIter->second;
      ++cmmMapIter;
      dataId   = cmmD->dataID();
      if (local  && dataId != LVL1::CMMJetHits::LOCAL_MAIN &&
                    dataId != LVL1::CMMJetHits::LOCAL_FORWARD)  continue;
      if (remote && dataId != LVL1::CMMJetHits::REMOTE_MAIN &&
		    dataId != LVL1::CMMJetHits::REMOTE_FORWARD) continue;
      if (total  && dataId != LVL1::CMMJetHits::TOTAL_MAIN &&
                    dataId != LVL1::CMMJetHits::TOTAL_FORWARD)  continue;
      if (etmap  && dataId != LVL1::CMMJetHits::ET_MAP)         continue;
      cmmHits = cmmD->Hits();
      crate   = cmmD->crate();

    } else {

      // Have both

      const LVL1::CMMJetHits* cmmS = cmmSimMapIter->second;
      const LVL1::CMMJetHits* cmmD = cmmMapIter->second;
      ++cmmSimMapIter;
      ++cmmMapIter;
      dataId   = cmmS->dataID();
      if (local  && dataId != LVL1::CMMJetHits::LOCAL_MAIN    &&
                    dataId != LVL1::CMMJetHits::LOCAL_FORWARD)  continue;
      if (remote && dataId != LVL1::CMMJetHits::LOCAL_MAIN    &&
                    dataId != LVL1::CMMJetHits::LOCAL_FORWARD &&
                    dataId != LVL1::CMMJetHits::REMOTE_MAIN   &&
		    dataId != LVL1::CMMJetHits::REMOTE_FORWARD) continue;
      if (total  && dataId != LVL1::CMMJetHits::TOTAL_MAIN    &&
                    dataId != LVL1::CMMJetHits::TOTAL_FORWARD)  continue;
      if (etmap  && dataId != LVL1::CMMJetHits::ET_MAP)         continue;
      cmmSimHits = cmmS->Hits();
      cmmHits    = cmmD->Hits();
      crate      = cmmS->crate();
    }

    if (!cmmSimHits && !cmmHits) continue;
    
    //  Fill in error plots

    if (local || total || etmap) {
      const int loc = crate * 2 + 1;
      const int cmmBins = 2 * 2;
      const int bit = (local) ? (1 << LocalJetMismatch)
                              : (total) ? (1 << TotalJetMismatch)
			                : (1 << JetEtMismatch);
      TH1F* hist = 0;
      if (cmmSimHits == cmmHits) {
        errors[loc] |= bit;
	hist = m_h_SumsSIMeqDAT;
      } else {
        errors[loc+cmmBins] |= bit;
	if (cmmSimHits && cmmHits) hist = m_h_SumsSIMneDAT;
	else if (!cmmHits)         hist = m_h_SumsSIMnoDAT;
	else                       hist = m_h_SumsDATnoSIM;
      }
      const int loc1 = (total) ? 3 : (etmap) ? 4 : crate;
      if (hist) hist->Fill(loc1);

      int nThresh = 8;
      int thrLen  = 3;
      unsigned int thrMask = 0x7;
      if (dataId == LVL1::CMMJetHits::LOCAL_FORWARD ||
          dataId == LVL1::CMMJetHits::TOTAL_FORWARD) {
        thrLen  = 2;
        thrMask = 0x3;
      } else if (dataId == LVL1::CMMJetHits::ET_MAP) {
        nThresh = 4;
	thrLen  = 1;
	thrMask = 0x1;
      }
      for (int thr = 0; thr < nThresh; ++thr) {
        const int loc2 = (thrLen == 2) ? thr+8 : thr;
        const int shift = thrLen*thr;
        const unsigned int d0 = (cmmHits >> shift) & thrMask;
        const unsigned int s0 = (cmmSimHits >> shift) & thrMask;
        if (d0 || s0) {
          if (d0 == s0) m_h_SumsThreshSIMeqDAT->Fill(loc1, loc2);
	  else          m_h_SumsThreshSIMneDAT->Fill(loc1, loc2);
        }
      }
    } else {
      if (dataId == LVL1::CMMJetHits::LOCAL_MAIN) {
        if (crate == 0) hitsSimMain = cmmSimHits;
      } else if (dataId == LVL1::CMMJetHits::LOCAL_FORWARD) {
        if (crate == 0) hitsSimFwd  = cmmSimHits;
      } else if (dataId == LVL1::CMMJetHits::REMOTE_MAIN) {
                        hitsDatMain = cmmHits;
      } else            hitsDatFwd  = cmmHits;
    }
  }
  if (remote) {
    const int crate = 1;
    const int loc = crate * 2 + 1;
    const int cmmBins = 2 * 2;
    const int bit = (1 << RemoteJetMismatch);

    TH1F* hist = 0;
    if (hitsSimMain && hitsSimMain == hitsDatMain) {
      errors[loc] |= bit;
      hist = m_h_SumsSIMeqDAT;
    } else if (hitsSimMain != hitsDatMain) {
      errors[loc+cmmBins] |= bit;
      if (hitsSimMain && hitsDatMain) hist = m_h_SumsSIMneDAT;
      else if (!hitsDatMain)          hist = m_h_SumsSIMnoDAT;
      else                            hist = m_h_SumsDATnoSIM;
    }
    const int loc1 = 2;
    if (hist) hist->Fill(loc1);
    hist = 0;
    if (hitsSimFwd && hitsSimFwd == hitsDatFwd) {
      errors[loc] |= bit;
      hist = m_h_SumsSIMeqDAT;
    } else if (hitsSimFwd != hitsDatFwd) {
      errors[loc+cmmBins] |= bit;
      if (hitsSimFwd && hitsDatFwd) hist = m_h_SumsSIMneDAT;
      else if (!hitsDatFwd)         hist = m_h_SumsSIMnoDAT;
      else                          hist = m_h_SumsDATnoSIM;
    }
    if (hist) hist->Fill(loc1);

    const int nThresh = 8;
    const int thrLen  = 3;
    const unsigned int thrMask = 0x7;
    const int thrLen2 = 2;
    const unsigned int thrMask2 = 0x3;
    for (int thr = 0; thr < nThresh; ++thr) {
      int shift = thrLen*thr;
      const unsigned int d0 = (hitsDatMain >> shift) & thrMask;
      const unsigned int s0 = (hitsSimMain >> shift) & thrMask;
      if (d0 || s0) {
        if (d0 == s0) m_h_SumsThreshSIMeqDAT->Fill(loc1, thr);
        else          m_h_SumsThreshSIMneDAT->Fill(loc1, thr);
      }
      shift = thrLen2*thr;
      const int thr2 = thr + nThresh;
      const unsigned int d1 = (hitsDatFwd >> shift) & thrMask2;
      const unsigned int s1 = (hitsSimFwd >> shift) & thrMask2;
      if (d1 || s1) {
        if (d1 == s1) m_h_SumsThreshSIMeqDAT->Fill(loc1, thr2);
        else          m_h_SumsThreshSIMneDAT->Fill(loc1, thr2);
      }
    }
  }
}

// Compare JetEt Map with JetEt RoI from data

void JEPSimBSMon::compare(const CmmJetHitsMap& cmmMap,
                          const LVL1::CMMRoI* cmmRoi,
                                ErrorVector& errors)
{
  m_log << MSG::DEBUG << "Compare JetEt Map with JetEt RoI from data"
        << endreq;

  int etMap = 0;
  int etRoi = 0;
  const int key = 100 + LVL1::CMMJetHits::ET_MAP;
  CmmJetHitsMap::const_iterator iter = cmmMap.find(key);
  if (iter != cmmMap.end()) {
    const LVL1::CMMJetHits* hits = iter->second;
    etMap = hits->Hits();
  }
  if (cmmRoi) etRoi = cmmRoi->jetEtHits();
  if (etMap || etRoi) {
    const int crate = 1;
    const int loc = crate * 2 + 1;
    const int cmmBins = 2 * 2;
    const int bit = (1 << JetEtRoIMismatch);
    TH1F* hist = 0;
    if (etMap == etRoi) {
      errors[loc] |= bit;
      hist = m_h_SumsSIMeqDAT;
    } else {
      errors[loc+cmmBins] |= bit;
      if (etMap && etRoi) hist = m_h_SumsSIMneDAT;
      else if (!etRoi)    hist = m_h_SumsSIMnoDAT;
      else                hist = m_h_SumsDATnoSIM;
    }
    const int loc1 = 5;
    if (hist) hist->Fill(loc1);

    const int nThresh = 4;
    const int thrLen  = 1;
    const unsigned int thrMask = 0x1;
    for (int thr = 0; thr < nThresh; ++thr) {
      int shift = thrLen*thr;
      const unsigned int d0 = (etRoi >> shift) & thrMask;
      const unsigned int s0 = (etMap >> shift) & thrMask;
      if (d0 || s0) {
        if (d0 == s0) m_h_SumsThreshSIMeqDAT->Fill(loc1, thr);
        else          m_h_SumsThreshSIMneDAT->Fill(loc1, thr);
      }
    }
  }
}

//  Compare simulated JEM Et Sums with data

void JEPSimBSMon::compare(const JemEtSumsMap& jemSimMap,
                          const JemEtSumsMap& jemMap,
                                ErrorVector& errors)
{
  m_log << MSG::DEBUG << "Compare simulated JEM Et Sums with data" << endreq;

  const int maxKey = 0x7fffffff;
  JemEtSumsMap::const_iterator simMapIter    = jemSimMap.begin();
  JemEtSumsMap::const_iterator simMapIterEnd = jemSimMap.end();
  JemEtSumsMap::const_iterator datMapIter    = jemMap.begin();
  JemEtSumsMap::const_iterator datMapIterEnd = jemMap.end();

  while (simMapIter != simMapIterEnd || datMapIter != datMapIterEnd) {

    int simKey = maxKey;
    int datKey = maxKey;
    unsigned int simEt = 0;
    unsigned int simEx = 0;
    unsigned int simEy = 0;
    unsigned int datEt = 0;
    unsigned int datEx = 0;
    unsigned int datEy = 0;
    int crate = 0;
    int jem   = 0;

    if (simMapIter != simMapIterEnd) simKey = simMapIter->first;
    if (datMapIter != datMapIterEnd) datKey = datMapIter->first;

    if ((datMapIter == datMapIterEnd) || (datKey > simKey)) {

      // Simulation EtSums but no data EtSums

      const LVL1::JEMEtSums* sime = simMapIter->second;
      simEt = sime->Et();
      simEx = sime->Ex();
      simEy = sime->Ey();
      crate = sime->crate();
      jem   = sime->module();
      ++simMapIter;

    } else if ((simMapIter == simMapIterEnd) || (simKey > datKey)) {

      // Data EtSums but no simulation EtSums

      const LVL1::JEMEtSums* date = datMapIter->second;
      datEt = date->Et();
      datEx = date->Ex();
      datEy = date->Ey();
      crate = date->crate();
      jem   = date->module();
      ++datMapIter;

    } else {

      // Have both

      const LVL1::JEMEtSums* sime = simMapIter->second;
      const LVL1::JEMEtSums* date = datMapIter->second;
      simEt = sime->Et();
      simEx = sime->Ex();
      simEy = sime->Ey();
      datEt = date->Et();
      datEx = date->Ex();
      datEy = date->Ey();
      crate = date->crate();
      jem   = date->module();
      ++simMapIter;
      ++datMapIter;
    }

    if (!simEt && !simEx && !simEy && !datEt && !datEx && !datEy) continue;
    
    //  Fill in error vector

    const int loc = crate * 16 + jem;
    const int jemBins = 2 * 16;
    const int bit = (1 << JEMEtSumsMismatch);
    if (simEt == datEt && simEx == datEx && simEy == datEy) errors[loc] |= bit;
    if (simEt != datEt || simEx != datEx || simEy != datEy)
                                                    errors[loc+jemBins] |= bit;
    if (simEt != datEt || simEx != datEx || simEy != datEy) {
      m_log << MSG::VERBOSE
            << "EtSums Mismatch Crate/JEM, Data Et/Ex/Ey, Sim Et/Ex/Ey: "
            << crate << "/" << jem << ", "
            << datEt << "/" << datEx << "/" << datEy << ", "
	    << simEt << "/" << simEx << "/" << simEy << endreq;
    }
    TH2F* hist = 0;
    if (simEx && simEx == datEx) hist = m_h_jemEtSumsSIMeqDAT;
    else if (simEx != datEx) {
      if (simEx && datEx) hist = m_h_jemEtSumsSIMneDAT;
      else if (!datEx)    hist = m_h_jemEtSumsSIMnoDAT;
      else                hist = m_h_jemEtSumsDATnoSIM;
    }
    if (hist) hist->Fill(loc, 0);
    hist = 0;
    if (simEy && simEy == datEy) hist = m_h_jemEtSumsSIMeqDAT;
    else if (simEy != datEy) {
      if (simEy && datEy) hist = m_h_jemEtSumsSIMneDAT;
      else if (!datEy)    hist = m_h_jemEtSumsSIMnoDAT;
      else                hist = m_h_jemEtSumsDATnoSIM;
    }
    if (hist) hist->Fill(loc, 1);
    hist = 0;
    if (simEt && simEt == datEt) hist = m_h_jemEtSumsSIMeqDAT;
    else if (simEt != datEt) {
      if (simEt && datEt) hist = m_h_jemEtSumsSIMneDAT;
      else if (!datEt)    hist = m_h_jemEtSumsSIMnoDAT;
      else                hist = m_h_jemEtSumsDATnoSIM;
    }
    if (hist) hist->Fill(loc, 2);
  }
}

//  Compare JEM EtSums and CMM EtSums

void JEPSimBSMon::compare(const JemEtSumsMap& jemMap,
                          const CmmEtSumsMap& cmmMap,
                                ErrorVector& errorsJEM,
				ErrorVector& errorsCMM)
{
  m_log << MSG::DEBUG << "Compare JEM EtSums and CMM EtSums" << endreq;

  const int maxKey = 0x7fffffff;
  JemEtSumsMap::const_iterator jemMapIter    = jemMap.begin();
  JemEtSumsMap::const_iterator jemMapIterEnd = jemMap.end();
  CmmEtSumsMap::const_iterator cmmMapIter    = cmmMap.begin();
  CmmEtSumsMap::const_iterator cmmMapIterEnd = cmmMap.end();

  while (jemMapIter != jemMapIterEnd || cmmMapIter != cmmMapIterEnd) {

    int jemKey = maxKey;
    int cmmKey = maxKey;
    unsigned int jemEt = 0;
    unsigned int jemEx = 0;
    unsigned int jemEy = 0;
    unsigned int cmmEt = 0;
    unsigned int cmmEx = 0;
    unsigned int cmmEy = 0;
    int crate  = 0;
    int jem    = 0;

    if (jemMapIter != jemMapIterEnd) jemKey = jemMapIter->first;
    if (cmmMapIter != cmmMapIterEnd) cmmKey = cmmMapIter->first;

    if ((cmmMapIter == cmmMapIterEnd) || (cmmKey > jemKey)) {

      // JEM EtSums but no CMM EtSums

      const LVL1::JEMEtSums* jeme = jemMapIter->second;
      jemEt = jeme->Et();
      jemEx = jeme->Ex();
      jemEy = jeme->Ey();
      crate = jeme->crate();
      jem   = jeme->module();
      ++jemMapIter;

    } else if ((jemMapIter == jemMapIterEnd) || (jemKey > cmmKey)) {

      // CMM EtSums but no JEM EtSums

      const LVL1::CMMEtSums* cmme = cmmMapIter->second;
      cmmEt = cmme->Et();
      cmmEx = cmme->Ex();
      cmmEy = cmme->Ey();
      crate = cmme->crate();
      jem   = cmme->dataID();
      ++cmmMapIter;
      if (jem > 15) continue;

    } else {

      // Have both

      const LVL1::JEMEtSums* jeme = jemMapIter->second;
      const LVL1::CMMEtSums* cmme = cmmMapIter->second;
      jemEt = jeme->Et();
      jemEx = jeme->Ex();
      jemEy = jeme->Ey();
      cmmEt = cmme->Et();
      cmmEx = cmme->Ex();
      cmmEy = cmme->Ey();
      crate = jeme->crate();
      jem   = jeme->module();
      ++jemMapIter;
      ++cmmMapIter;
    }

    if (!jemEt && !jemEx && !jemEy && !cmmEt && !cmmEx && !cmmEy) continue;
    
    //  Fill in error vectors

    const int loc  = crate * 16 + jem;
    const int loc2 = crate * 2; // CMM-Energy is left ==> cmm 0
    const int jemBins = 2 * 16;
    const int cmmBins = 2 * 2;
    const int bit = (1 << CMMEtSumsMismatch);
    if (jemEt == cmmEt && jemEx == cmmEx && jemEy == cmmEy) {
      errorsJEM[loc]  |= bit;
      errorsCMM[loc2] |= bit;
    } else {
      errorsJEM[loc+jemBins]  |= bit;
      errorsCMM[loc2+cmmBins] |= bit;
    }
    TH2F* hist = 0;
    if (jemEx && jemEx == cmmEx) hist = m_h_cmmEtSumsSIMeqDAT;
    else if (jemEx != cmmEx) {
      if (jemEx && cmmEx) hist = m_h_cmmEtSumsSIMneDAT;
      else if (!cmmEx)    hist = m_h_cmmEtSumsSIMnoDAT;
      else                hist = m_h_cmmEtSumsDATnoSIM;
    }
    if (hist) hist->Fill(loc, 0);
    hist = 0;
    if (jemEy && jemEy == cmmEy) hist = m_h_cmmEtSumsSIMeqDAT;
    else if (jemEy != cmmEy) {
      if (jemEy && cmmEy) hist = m_h_cmmEtSumsSIMneDAT;
      else if (!cmmEy)    hist = m_h_cmmEtSumsSIMnoDAT;
      else                hist = m_h_cmmEtSumsDATnoSIM;
    }
    if (hist) hist->Fill(loc, 1);
    hist = 0;
    if (jemEt && jemEt == cmmEt) hist = m_h_cmmEtSumsSIMeqDAT;
    else if (jemEt != cmmEt) {
      if (jemEt && cmmEt) hist = m_h_cmmEtSumsSIMneDAT;
      else if (!cmmEt)    hist = m_h_cmmEtSumsSIMnoDAT;
      else                hist = m_h_cmmEtSumsDATnoSIM;
    }
    if (hist) hist->Fill(loc, 2);
  }
}

//  Compare Simulated CMM EtSums and Data CMM EtSums

void JEPSimBSMon::compare(const CmmEtSumsMap& cmmSimMap,
                          const CmmEtSumsMap& cmmMap,
                                ErrorVector& errors, int selection)
{
  m_log << MSG::DEBUG << "Compare Simulated CMM EtSums and Data CMM EtSums"
        << endreq;

  const bool local  = (selection == LVL1::CMMEtSums::LOCAL);
  const bool remote = (selection == LVL1::CMMEtSums::REMOTE);
  const bool total  = (selection == LVL1::CMMEtSums::TOTAL);
  const bool etmaps = (selection == LVL1::CMMEtSums::SUM_ET_MAP);
  if (!local && !remote && !total && !etmaps) return;
  unsigned int localEt  = 0;
  unsigned int localEx  = 0;
  unsigned int localEy  = 0;
  unsigned int remoteEt = 0;
  unsigned int remoteEx = 0;
  unsigned int remoteEy = 0;
  const int maxKey = 0x7fffffff;
  CmmEtSumsMap::const_iterator cmmSimMapIter    = cmmSimMap.begin();
  CmmEtSumsMap::const_iterator cmmSimMapIterEnd = cmmSimMap.end();
  CmmEtSumsMap::const_iterator cmmMapIter       = cmmMap.begin();
  CmmEtSumsMap::const_iterator cmmMapIterEnd    = cmmMap.end();

  while (cmmSimMapIter != cmmSimMapIterEnd || cmmMapIter != cmmMapIterEnd) {

    int cmmSimKey = maxKey;
    int cmmKey    = maxKey;
    unsigned int cmmSimEt = 0;
    unsigned int cmmSimEx = 0;
    unsigned int cmmSimEy = 0;
    unsigned int cmmEt = 0;
    unsigned int cmmEx = 0;
    unsigned int cmmEy = 0;
    int crate  = 0;
    int dataId = 0;

    if (cmmSimMapIter != cmmSimMapIterEnd) cmmSimKey = cmmSimMapIter->first;
    if (cmmMapIter    != cmmMapIterEnd)    cmmKey    = cmmMapIter->first;

    if ((cmmMapIter == cmmMapIterEnd) || (cmmKey > cmmSimKey)) {

      // Sim CMM EtSums but no Data CMM EtSums

      const LVL1::CMMEtSums* cmmS = cmmSimMapIter->second;
      ++cmmSimMapIter;
      dataId = cmmS->dataID();
      if (local  && dataId != LVL1::CMMEtSums::LOCAL) continue;
      if (remote && dataId != LVL1::CMMEtSums::LOCAL) continue;
      if (total  && dataId != LVL1::CMMEtSums::TOTAL) continue;
      if (etmaps && dataId != LVL1::CMMEtSums::SUM_ET_MAP &&
                    dataId != LVL1::CMMEtSums::MISSING_ET_MAP) continue;
      cmmSimEt = cmmS->Et();
      cmmSimEx = cmmS->Ex();
      cmmSimEy = cmmS->Ey();
      if (!etmaps) { // include overflow bit in test
        cmmSimEt |= ((cmmS->EtError() & 0x1) << 15);
        cmmSimEx |= ((cmmS->ExError() & 0x1) << 15);
        cmmSimEy |= ((cmmS->EyError() & 0x1) << 15);
      }
      crate    = cmmS->crate();

    } else if ((cmmSimMapIter == cmmSimMapIterEnd) || (cmmSimKey > cmmKey)) {

      // Data CMM EtSums but no Sim CMM EtSums

      const LVL1::CMMEtSums* cmmD = cmmMapIter->second;
      ++cmmMapIter;
      dataId   = cmmD->dataID();
      if (local  && dataId != LVL1::CMMEtSums::LOCAL)  continue;
      if (remote && dataId != LVL1::CMMEtSums::REMOTE) continue;
      if (total  && dataId != LVL1::CMMEtSums::TOTAL)  continue;
      if (etmaps && dataId != LVL1::CMMEtSums::SUM_ET_MAP &&
                    dataId != LVL1::CMMEtSums::MISSING_ET_MAP) continue;
      cmmEt = cmmD->Et();
      cmmEx = cmmD->Ex();
      cmmEy = cmmD->Ey();
      if (!etmaps) {
        cmmEt |= ((cmmD->EtError() & 0x1) << 15);
        cmmEx |= ((cmmD->ExError() & 0x1) << 15);
        cmmEy |= ((cmmD->EyError() & 0x1) << 15);
      }
      crate = cmmD->crate();

    } else {

      // Have both

      const LVL1::CMMEtSums* cmmS = cmmSimMapIter->second;
      const LVL1::CMMEtSums* cmmD = cmmMapIter->second;
      ++cmmSimMapIter;
      ++cmmMapIter;
      dataId   = cmmS->dataID();
      if (local  && dataId != LVL1::CMMEtSums::LOCAL)  continue;
      if (remote && dataId != LVL1::CMMEtSums::LOCAL &&
                    dataId != LVL1::CMMEtSums::REMOTE) continue;
      if (total  && dataId != LVL1::CMMEtSums::TOTAL)  continue;
      if (etmaps && dataId != LVL1::CMMEtSums::SUM_ET_MAP &&
                    dataId != LVL1::CMMEtSums::MISSING_ET_MAP) continue;
      cmmSimEt = cmmS->Et();
      cmmSimEx = cmmS->Ex();
      cmmSimEy = cmmS->Ey();
      cmmEt    = cmmD->Et();
      cmmEx    = cmmD->Ex();
      cmmEy    = cmmD->Ey();
      if (!etmaps) {
        cmmSimEt |= ((cmmS->EtError() & 0x1) << 15);
        cmmSimEx |= ((cmmS->ExError() & 0x1) << 15);
        cmmSimEy |= ((cmmS->EyError() & 0x1) << 15);
        cmmEt    |= ((cmmD->EtError() & 0x1) << 15);
        cmmEx    |= ((cmmD->ExError() & 0x1) << 15);
        cmmEy    |= ((cmmD->EyError() & 0x1) << 15);
      }
      crate    = cmmS->crate();
    }

    if (!cmmSimEt && !cmmSimEx && !cmmSimEy && !cmmEt && !cmmEx && !cmmEy)
                                                                     continue;
    
    //  Fill in error vectors

    if (local || total || etmaps) {
      const int loc = crate * 2;
      const int cmmBins = 2 * 2;
      const int bit = (local)
                        ? (1 << LocalEnergyMismatch)
                        : (total)
			    ? (1 << TotalEnergyMismatch)
			    : (dataId == LVL1::CMMEtSums::SUM_ET_MAP)
			       ? (1 << SumEtMismatch)
			       : (1 << MissingEtMismatch);
      if (cmmSimEt == cmmEt && cmmSimEx == cmmEx && cmmSimEy == cmmEy) {
        errors[loc] |= bit;
      } else errors[loc+cmmBins] |= bit;
      const int loc1 = (local) ? crate : 3;
      if (local || total) {
        TH2F* hist = 0;
	if (cmmSimEx && cmmSimEx == cmmEx) hist = m_h_EnSumsSIMeqDAT;
	else if (cmmSimEx != cmmEx) {
	  if (cmmSimEx && cmmEx) hist = m_h_EnSumsSIMneDAT;
	  else if (!cmmEx)       hist = m_h_EnSumsSIMnoDAT;
	  else                   hist = m_h_EnSumsDATnoSIM;
        }
	if (hist) hist->Fill(loc1, 0);
        hist = 0;
	if (cmmSimEy && cmmSimEy == cmmEy) hist = m_h_EnSumsSIMeqDAT;
	else if (cmmSimEy != cmmEy) {
	  if (cmmSimEy && cmmEy) hist = m_h_EnSumsSIMneDAT;
	  else if (!cmmEy)       hist = m_h_EnSumsSIMnoDAT;
	  else                   hist = m_h_EnSumsDATnoSIM;
        }
	if (hist) hist->Fill(loc1, 1);
        hist = 0;
	if (cmmSimEt && cmmSimEt == cmmEt) hist = m_h_EnSumsSIMeqDAT;
	else if (cmmSimEt != cmmEt) {
	  if (cmmSimEt && cmmEt) hist = m_h_EnSumsSIMneDAT;
	  else if (!cmmEt)       hist = m_h_EnSumsSIMnoDAT;
	  else                   hist = m_h_EnSumsDATnoSIM;
        }
	if (hist) hist->Fill(loc1, 2);
      } else {
        const int loc2 = (dataId == LVL1::CMMEtSums::SUM_ET_MAP) ? 3 : 4;
        TH2F* hist = 0;
	if (cmmSimEt && cmmSimEt == cmmEt) hist = m_h_EnSumsSIMeqDAT;
	else if (cmmSimEt != cmmEt) {
	  if (cmmSimEt && cmmEt) hist = m_h_EnSumsSIMneDAT;
	  else if (!cmmEt)       hist = m_h_EnSumsSIMnoDAT;
	  else                   hist = m_h_EnSumsDATnoSIM;
        }
	if (hist) hist->Fill(loc1, loc2);
	if (cmmSimEt || cmmEt) {
	  int loc3 = 0;
	  int nThresh = 4;
	  if (dataId == LVL1::CMMEtSums::MISSING_ET_MAP) {
	    loc3 = 2;
	    nThresh = 8;
          }
	  for (int thr = 0; thr < nThresh; ++thr) {
	    const int thrDat = (cmmEt    >> thr) & 0x1;
	    const int thrSim = (cmmSimEt >> thr) & 0x1;
	    if (thrDat || thrSim) {
	      if (thrDat == thrSim) m_h_EnSumsThreshSIMeqDAT->Fill(loc3, thr);
	      else                  m_h_EnSumsThreshSIMneDAT->Fill(loc3, thr);
	    }
	  }
        }
      }
    } else {
      if (dataId == LVL1::CMMEtSums::LOCAL) {
        if (crate == 0) {
	  localEt = cmmSimEt;
	  localEx = cmmSimEx;
	  localEy = cmmSimEy;
	}
      } else {
        remoteEt = cmmEt;
        remoteEx = cmmEx;
        remoteEy = cmmEy;
      }
    }
  }
  if (remote && (localEt || localEx || localEy ||
                 remoteEt || remoteEx || remoteEy)) {
    const int crate = 1;
    const int loc = crate * 2;
    const int cmmBins = 2 * 2;
    const int bit = (1 << RemoteEnergyMismatch);
    if (localEt == remoteEt && localEx == remoteEx && localEy == remoteEy) {
      errors[loc] |= bit;
    } else errors[loc+cmmBins] |= bit;
    TH2F* hist = 0;
    if (localEx && localEx == remoteEx) hist = m_h_EnSumsSIMeqDAT;
    else if (localEx != remoteEx) {
      if (localEx && remoteEx) hist = m_h_EnSumsSIMneDAT;
      else if (!remoteEx)      hist = m_h_EnSumsSIMnoDAT;
      else                     hist = m_h_EnSumsDATnoSIM;
    }
    if (hist) hist->Fill(2, 0);
    hist = 0;
    if (localEy && localEy == remoteEy) hist = m_h_EnSumsSIMeqDAT;
    else if (localEy != remoteEy) {
      if (localEy && remoteEy) hist = m_h_EnSumsSIMneDAT;
      else if (!remoteEy)      hist = m_h_EnSumsSIMnoDAT;
      else                     hist = m_h_EnSumsDATnoSIM;
    }
    if (hist) hist->Fill(2, 1);
    hist = 0;
    if (localEt && localEt == remoteEt) hist = m_h_EnSumsSIMeqDAT;
    else if (localEt != remoteEt) {
      if (localEt && remoteEt) hist = m_h_EnSumsSIMneDAT;
      else if (!remoteEt)      hist = m_h_EnSumsSIMnoDAT;
      else                     hist = m_h_EnSumsDATnoSIM;
    }
    if (hist) hist->Fill(2, 2);
  }
}

// Compare Et Maps and Energy Totals with Energy RoIs from data

void JEPSimBSMon::compare(const CmmEtSumsMap& cmmMap,
                          const LVL1::CMMRoI* cmmRoi,
                                ErrorVector& errors)
{
  m_log << MSG::DEBUG << "Compare Et Maps and Energy Totals with RoIs from data"
        << endreq;

  int sumEtMap = 0;
  int missEtMap = 0;
  int et = 0;
  int ex = 0;
  int ey = 0;
  int sumEtRoi = 0;
  int missEtRoi = 0;
  int etRoi = 0;
  int exRoi = 0;
  int eyRoi = 0;
  int key = 100 + LVL1::CMMEtSums::SUM_ET_MAP;
  CmmEtSumsMap::const_iterator iter = cmmMap.find(key);
  if (iter != cmmMap.end()) {
    const LVL1::CMMEtSums* sums = iter->second;
    sumEtMap = sums->Et();
  }
  key = 100 + LVL1::CMMEtSums::MISSING_ET_MAP;
  iter = cmmMap.find(key);
  if (iter != cmmMap.end()) {
    const LVL1::CMMEtSums* sums = iter->second;
    missEtMap = sums->Et();
  }
  key = 100 + LVL1::CMMEtSums::TOTAL;
  iter = cmmMap.find(key);
  if (iter != cmmMap.end()) {
    const LVL1::CMMEtSums* sums = iter->second;
    et = (sums->Et() | ((sums->EtError() & 0x1) << 15));
    ex = (sums->Ex() | ((sums->ExError() & 0x1) << 15));
    ey = (sums->Ey() | ((sums->EyError() & 0x1) << 15));
  }
  if (cmmRoi) {
    sumEtRoi = cmmRoi->sumEtHits();
    missEtRoi = cmmRoi->missingEtHits();
    etRoi = (cmmRoi->et() | ((cmmRoi->etError() & 0x1) << 15));
    exRoi = (cmmRoi->ex() | ((cmmRoi->exError() & 0x1) << 15));
    eyRoi = (cmmRoi->ey() | ((cmmRoi->eyError() & 0x1) << 15));
  }
  if (sumEtMap || sumEtRoi || missEtMap || missEtRoi || et || etRoi ||
                                           ex || exRoi || ey || eyRoi) {
    const int crate = 1;
    const int loc = crate * 2;
    const int cmmBins = 2 * 2;
    const int bit = (1 << EnergyRoIMismatch);
    if (sumEtMap == sumEtRoi && missEtMap == missEtRoi &&
        et == etRoi && ex == exRoi && ey == eyRoi) errors[loc] |= bit;
    else errors[loc+cmmBins] |= bit;
    TH2F* hist = 0;
    if (ex && ex == exRoi) hist = m_h_EnSumsSIMeqDAT;
    else if (ex != exRoi) {
      if (ex && exRoi) hist = m_h_EnSumsSIMneDAT;
      else if (!exRoi) hist = m_h_EnSumsSIMnoDAT;
      else             hist = m_h_EnSumsDATnoSIM;
    }
    if (hist) hist->Fill(4, 0);
    hist = 0;
    if (ey && ey == eyRoi) hist = m_h_EnSumsSIMeqDAT;
    else if (ey != eyRoi) {
      if (ey && eyRoi) hist = m_h_EnSumsSIMneDAT;
      else if (!eyRoi) hist = m_h_EnSumsSIMnoDAT;
      else             hist = m_h_EnSumsDATnoSIM;
    }
    if (hist) hist->Fill(4, 1);
    hist = 0;
    if (et && et == etRoi) hist = m_h_EnSumsSIMeqDAT;
    else if (et != etRoi) {
      if (et && etRoi) hist = m_h_EnSumsSIMneDAT;
      else if (!etRoi) hist = m_h_EnSumsSIMnoDAT;
      else             hist = m_h_EnSumsDATnoSIM;
    }
    if (hist) hist->Fill(4, 2);
    hist = 0;
    if (sumEtMap && sumEtMap == sumEtRoi) hist = m_h_EnSumsSIMeqDAT;
    else if (sumEtMap != sumEtRoi) {
      if (sumEtMap && sumEtRoi) hist = m_h_EnSumsSIMneDAT;
      else if (!sumEtRoi)       hist = m_h_EnSumsSIMnoDAT;
      else                      hist = m_h_EnSumsDATnoSIM;
    }
    if (hist) hist->Fill(4, 3);
    hist = 0;
    if (missEtMap && missEtMap == missEtRoi) hist = m_h_EnSumsSIMeqDAT;
    else if (missEtMap != missEtRoi) {
      if (missEtMap && missEtRoi) hist = m_h_EnSumsSIMneDAT;
      else if (!missEtRoi)        hist = m_h_EnSumsSIMnoDAT;
      else                        hist = m_h_EnSumsDATnoSIM;
    }
    if (hist) hist->Fill(4, 4);
    if (sumEtMap || sumEtRoi) {
      for (int thr = 0; thr < 4; ++thr) {
        const int thrDat = (sumEtRoi >> thr) & 0x1;
        const int thrSim = (sumEtMap >> thr) & 0x1;
	if (thrDat || thrSim) {
          if (thrDat == thrSim) m_h_EnSumsThreshSIMeqDAT->Fill(1, thr);
          else                  m_h_EnSumsThreshSIMneDAT->Fill(1, thr);
        }
      }
    }
    if (missEtMap || missEtRoi) {
      for (int thr = 0; thr < 8; ++thr) {
        const int thrDat = (missEtRoi >> thr) & 0x1;
        const int thrSim = (missEtMap >> thr) & 0x1;
	if (thrDat || thrSim) {
          if (thrDat == thrSim) m_h_EnSumsThreshSIMeqDAT->Fill(3, thr);
          else                  m_h_EnSumsThreshSIMneDAT->Fill(3, thr);
        }
      }
    }
  }
}

void JEPSimBSMon::fillEventSample(int err, int loc, bool isJem)
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
    int hist = 0;
    int y    = 0;
    if (isJem) {
      hist = (err < 5) ? err : err - 5;
      y    = loc;
    } else {
      hist = (err < 10) ? 7 : 8;
      y = loc/2;
      if (err == LocalJetMismatch  || err == LocalEnergyMismatch)  y += 2;
      if (err == RemoteJetMismatch || err == RemoteEnergyMismatch) y = 4;
      if (err == TotalJetMismatch  || err == TotalEnergyMismatch)  y = 5;
      if (err == JetEtMismatch     || err == SumEtMismatch)        y = 6;
      if (err == JetEtRoIMismatch  || err == MissingEtMismatch)    y = 7;
      if (err == EnergyRoIMismatch)                                y = 8;
    }
    const int count = hist*32 + y;
    const int x = m_sampleCounts[count];
    if (x < m_eventSamples && m_sampleHists[hist]) {
      m_sampleHists[hist]->Fill(x, y, m_eventNumber);
      ++m_sampleCounts[count];
    }
  }
}

void JEPSimBSMon::setLabels(TH2* hist)
{
  hist->GetXaxis()->SetBinLabel(1, "JEM");
  hist->GetXaxis()->SetBinLabel(5, "0/4");
  hist->GetXaxis()->SetBinLabel(9, "0/8");
  hist->GetXaxis()->SetBinLabel(13, "0/12");
  hist->GetXaxis()->SetBinLabel(17, "1/0");
  hist->GetXaxis()->SetBinLabel(21, "1/4");
  hist->GetXaxis()->SetBinLabel(25, "1/8");
  hist->GetXaxis()->SetBinLabel(29, "1/12");
  hist->GetXaxis()->SetBinLabel(33, "CMM");
  hist->GetXaxis()->SetBinLabel(35, "1/L");
  hist->GetXaxis()->SetTitleOffset(1.25);
  // Simulation steps in red (#color[2]) depend on Trigger Menu
  hist->GetYaxis()->SetBinLabel(1+EMElementMismatch,    "EM je");
  hist->GetYaxis()->SetBinLabel(1+HadElementMismatch,   "Had je");
  hist->GetYaxis()->SetBinLabel(1+RoIMismatch,          "#color[2]{RoIs}");
  hist->GetYaxis()->SetBinLabel(1+JEMHitsMismatch,      "JEMHits");
  hist->GetYaxis()->SetBinLabel(1+CMMJetHitsMismatch,   "CMMHits");
  hist->GetYaxis()->SetBinLabel(1+LocalJetMismatch,     "Local");
  hist->GetYaxis()->SetBinLabel(1+RemoteJetMismatch,    "Remote");
  hist->GetYaxis()->SetBinLabel(1+TotalJetMismatch,     "Total");
  hist->GetYaxis()->SetBinLabel(1+JetEtMismatch,        "#color[2]{JetEt}");
  hist->GetYaxis()->SetBinLabel(1+JetEtRoIMismatch,     "JetEt RoI");
  hist->GetYaxis()->SetBinLabel(1+JEMEtSumsMismatch,    "JEMSums");
  hist->GetYaxis()->SetBinLabel(1+CMMEtSumsMismatch,    "CMMSums");
  hist->GetYaxis()->SetBinLabel(1+LocalEnergyMismatch,  "Local");
  hist->GetYaxis()->SetBinLabel(1+RemoteEnergyMismatch, "Remote");
  hist->GetYaxis()->SetBinLabel(1+TotalEnergyMismatch,  "Total");
  hist->GetYaxis()->SetBinLabel(1+SumEtMismatch,        "#color[2]{SumEt}");
  hist->GetYaxis()->SetBinLabel(1+MissingEtMismatch,    "#color[2]{MissingEt}");
  hist->GetYaxis()->SetBinLabel(1+EnergyRoIMismatch,    "Engy RoIs");
  if (!m_compareWithSim) {
    hist->GetYaxis()->SetBinLabel(1+EMElementMismatch,  "#color[16]{EM je}");
    hist->GetYaxis()->SetBinLabel(1+HadElementMismatch, "#color[16]{Had je}");
    hist->GetYaxis()->SetBinLabel(1+RoIMismatch,        "#color[16]{RoIs}");
    hist->GetYaxis()->SetBinLabel(1+JEMHitsMismatch,    "#color[16]{JEMHits}");
    hist->GetYaxis()->SetBinLabel(1+LocalJetMismatch,   "#color[16]{Local}");
    hist->GetYaxis()->SetBinLabel(1+TotalJetMismatch,   "#color[16]{Total}");
    hist->GetYaxis()->SetBinLabel(1+JetEtMismatch,      "#color[16]{JetEt}");
    hist->GetYaxis()->SetBinLabel(1+JEMEtSumsMismatch,  "#color[16]{JEMSums}");
    hist->GetYaxis()->SetBinLabel(1+LocalEnergyMismatch,"#color[16]{Local}");
    hist->GetYaxis()->SetBinLabel(1+TotalEnergyMismatch,"#color[16]{Total}");
    hist->GetYaxis()->SetBinLabel(1+SumEtMismatch,      "#color[16]{SumEt}");
    hist->GetYaxis()->SetBinLabel(1+MissingEtMismatch, "#color[16]{MissingEt}");
  }
  hist->GetYaxis()->SetLabelSize(0.045);
}

void JEPSimBSMon::setLabelsCMFC(TH2* hist)
{
  setLabelsJEM(hist);
  for (int frame = 0; frame < 8; ++frame) {
    for (int loc = 0; loc < 4; loc += 2) {
      std::ostringstream cnum;
      cnum << frame << "/" << loc;
      hist->GetYaxis()->SetBinLabel(frame*4 + loc + 1, cnum.str().c_str());
    }
  }
}

void JEPSimBSMon::setLabelsJMS(TH2* hist)
{
  setLabelsJEM(hist, false);
  if (m_eventSamples <= 10) setLabelsXNUM(hist, 1, m_eventSamples);
}

void JEPSimBSMon::setLabelsCMT(TH2* hist)
{
  setLabelsJEM(hist);
  setLabelsYNUM(hist, 0, 7);
  hist->GetYaxis()->SetBinLabel(1, "Main 0");
  hist->GetYaxis()->SetBinLabel(9, "Fwd 0");
  hist->GetYaxis()->SetBinLabel(10, "1");
  hist->GetYaxis()->SetBinLabel(11, "2");
  hist->GetYaxis()->SetBinLabel(12, "3");
}

void JEPSimBSMon::setLabelsYNUM(TH2* hist, int beg, int end)
{
  int bin = 1;
  for (int val = beg; val <= end; ++val) {
    std::ostringstream cnum;
    cnum << val;
    hist->GetYaxis()->SetBinLabel(bin++, cnum.str().c_str());
  }
  hist->GetYaxis()->SetLabelSize(0.05);
}

void JEPSimBSMon::setLabelsXNUM(TH2* hist, int beg, int end)
{
  int bin = 1;
  for (int val = beg; val <= end; ++val) {
    std::ostringstream cnum;
    cnum << val;
    hist->GetXaxis()->SetBinLabel(bin++, cnum.str().c_str());
  }
  hist->GetXaxis()->SetLabelSize(0.05);
}

void JEPSimBSMon::setLabelsJEM(TH2* hist, bool xAxis)
{
  const int nJEMs = 16;
  for (int crate = 0; crate < 2; ++crate) {
    for (int module = 0; module < nJEMs; module += 2) {
      std::ostringstream cnum;
      cnum << crate << "/" << module;
      if (xAxis) {
        hist->GetXaxis()->SetBinLabel(crate*nJEMs + module + 1,
                                                        cnum.str().c_str());
      } else {
        hist->GetYaxis()->SetBinLabel(crate*nJEMs + module + 1,
                                                        cnum.str().c_str());
      }
    }
  }
}

void JEPSimBSMon::setLabelsMC(TH2* hist)
{
  setLabelsXNUM(hist, 0, 15);
  setLabelsYNUM(hist, 0, 1);
}

void JEPSimBSMon::setLabelsSH(TH1* hist)
{
  hist->GetXaxis()->SetBinLabel(1, "Local0");
  hist->GetXaxis()->SetBinLabel(2, "Local1");
  hist->GetXaxis()->SetBinLabel(3, "Remote");
  hist->GetXaxis()->SetBinLabel(4, "Total");
  hist->GetXaxis()->SetBinLabel(5, "JetEt");
  hist->GetXaxis()->SetBinLabel(6, "JetEt RoI");
  if (!m_compareWithSim) {
    hist->GetXaxis()->SetBinLabel(1, "#color[16]{Local0}");
    hist->GetXaxis()->SetBinLabel(2, "#color[16]{Local1}");
    hist->GetXaxis()->SetBinLabel(4, "#color[16]{Total}");
    hist->GetXaxis()->SetBinLabel(5, "#color[16]{JetEt}");
  }
}

void JEPSimBSMon::setLabelsSHF(TH2* hist)
{
  setLabelsSH(hist);
  setLabelsYNUM(hist, 0, 7);
  hist->GetYaxis()->SetBinLabel(1, "Main 0");
  hist->GetYaxis()->SetBinLabel(9, "FwdL 0");
  hist->GetYaxis()->SetBinLabel(10, "1");
  hist->GetYaxis()->SetBinLabel(11, "2");
  hist->GetYaxis()->SetBinLabel(12, "3");
  hist->GetYaxis()->SetBinLabel(13, "FwdR 0");
  hist->GetYaxis()->SetBinLabel(14, "1");
  hist->GetYaxis()->SetBinLabel(15, "2");
  hist->GetYaxis()->SetBinLabel(16, "3");
}

void JEPSimBSMon::setLabelsJES(TH2* hist)
{
  setLabelsJEM(hist);
  hist->GetYaxis()->SetBinLabel(1, "Ex");
  hist->GetYaxis()->SetBinLabel(2, "Ey");
  hist->GetYaxis()->SetBinLabel(3, "Et");
}

void JEPSimBSMon::setLabelsEnTot(TH2* hist)
{
  hist->GetXaxis()->SetBinLabel(1, "Local0");
  hist->GetXaxis()->SetBinLabel(2, "Local1");
  hist->GetXaxis()->SetBinLabel(3, "Remote");
  hist->GetXaxis()->SetBinLabel(4, "Total");
  hist->GetXaxis()->SetBinLabel(5, "RoI");
  hist->GetYaxis()->SetBinLabel(1, "Ex");
  hist->GetYaxis()->SetBinLabel(2, "Ey");
  hist->GetYaxis()->SetBinLabel(3, "Et");
  hist->GetYaxis()->SetBinLabel(4, "SumEt");
  hist->GetYaxis()->SetBinLabel(5, "MissingEt");
  if (!m_compareWithSim) {
    hist->GetXaxis()->SetBinLabel(1, "#color[16]{Local0}");
    hist->GetXaxis()->SetBinLabel(2, "#color[16]{Local1}");
    hist->GetXaxis()->SetBinLabel(4, "#color[16]{Total}");
  }
}

void JEPSimBSMon::setLabelsEnTotThr(TH2* hist)
{
  hist->GetXaxis()->SetBinLabel(1, "SumEt");
  hist->GetXaxis()->SetBinLabel(2, "SumEt RoI");
  hist->GetXaxis()->SetBinLabel(3, "MissingEt");
  hist->GetXaxis()->SetBinLabel(4, "MissingEt RoI");
  if (!m_compareWithSim) {
    hist->GetXaxis()->SetBinLabel(1, "#color[16]{SumEt}");
    hist->GetXaxis()->SetBinLabel(3, "#color[16]{MissingEt}");
  }
  setLabelsYNUM(hist, 0, 7);
}

void JEPSimBSMon::setupMap(const JetElementCollection* coll,
                                 JetElementMap& map)
{
  if (coll) m_jetElementTool->mapJetElements(coll, &map);
}

void JEPSimBSMon::setupMap(const JemRoiCollection* coll, JemRoiMap& map)
{
  if (coll) {
    JemRoiCollection::const_iterator pos  = coll->begin();
    JemRoiCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const int crate = (*pos)->crate();
      const int jem   = (*pos)->jem();
      const int frame = (*pos)->frame();
      const int loc   = (*pos)->location();
      const int fwd   = (*pos)->forward();
      const int key   =
             (((((((crate << 4) | jem) << 3) | frame) << 3) | loc) << 1) | fwd;
      map.insert(std::make_pair(key, *pos));
    }
  }
}

void JEPSimBSMon::setupMap(const JemHitsCollection* coll, JemHitsMap& map)
{
  if (coll) {
    JemHitsCollection::const_iterator pos  = coll->begin();
    JemHitsCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const int crate = (*pos)->crate();
      const int jem   = (*pos)->module();
      const int key   = crate * 100 + jem;
      map.insert(std::make_pair(key, *pos));
    }
  }
}

void JEPSimBSMon::setupMap(const CmmJetHitsCollection* coll, CmmJetHitsMap& map)
{
  if (coll) {
    CmmJetHitsCollection::const_iterator pos  = coll->begin();
    CmmJetHitsCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const int crate  = (*pos)->crate();
      const int dataId = (*pos)->dataID();
      const int key  = crate * 100 + dataId;
      map.insert(std::make_pair(key, *pos));
    }
  }
}

void JEPSimBSMon::setupMap(const JemEtSumsCollection* coll, JemEtSumsMap& map)
{
  if (coll) {
    JemEtSumsCollection::const_iterator pos  = coll->begin();
    JemEtSumsCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const int crate = (*pos)->crate();
      const int jem   = (*pos)->module();
      const int key   = crate * 100 + jem;
      map.insert(std::make_pair(key, *pos));
    }
  }
}

void JEPSimBSMon::setupMap(const CmmEtSumsCollection* coll, CmmEtSumsMap& map)
{
  if (coll) {
    CmmEtSumsCollection::const_iterator pos  = coll->begin();
    CmmEtSumsCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const int crate  = (*pos)->crate();
      const int dataId = (*pos)->dataID();
      const int key  = crate * 100 + dataId;
      map.insert(std::make_pair(key, *pos));
    }
  }
}

void JEPSimBSMon::simulate(const TriggerTowerCollection* towers,
                                 JetElementCollection* elements)
{
  m_log << MSG::DEBUG << "Simulate Jet Elements from Trigger Towers" << endreq;

  // Make zero-suppressed collection to speed up simulation

  TriggerTowerCollection* towersZ =
                              new TriggerTowerCollection(SG::VIEW_ELEMENTS);
  TriggerTowerCollection::const_iterator pos  = towers->begin();
  TriggerTowerCollection::const_iterator posE = towers->end();
  for (; pos != posE; ++pos) {
    if ((*pos)->emEnergy() > 0 || (*pos)->hadEnergy() > 0) {
      towersZ->push_back(*pos);
    }
  }
  m_jetElementTool->makeJetElements(towersZ, elements);
  delete towersZ;
}

void JEPSimBSMon::simulate(const JetElementCollection* elements,
                           const JetElementCollection* elementsOv,
                                 JemRoiCollection* rois)
{
  m_log << MSG::DEBUG << "Simulate JEM RoIs from Jet Elements" << endreq;

  // Process a crate at a time to use overlap data
  const int ncrates = 2;
  std::vector<JetElementCollection*> crateColl;
  for (int crate = 0; crate < ncrates; ++crate) {
    crateColl.push_back(new JetElementCollection(SG::VIEW_ELEMENTS));
  }
  LVL1::CoordToHardware converter;
  JetElementCollection::const_iterator iter;
  JetElementCollection::const_iterator iterE;
  if (elements) {  // core data
    iter  = elements->begin();
    iterE = elements->end();
    for (; iter != iterE; ++iter) {
      LVL1::JetElement* je = *iter;
      const LVL1::Coordinate coord(je->phi(), je->eta());
      const int crate = converter.jepCrate(coord);
      if (crate < ncrates) crateColl[crate]->push_back(je);
    }
  }
  if (elementsOv) {  // overlap data
    iter  = elementsOv->begin();
    iterE = elementsOv->end();
    for (; iter != iterE; ++iter) {
      LVL1::JetElement* je = *iter;
      const LVL1::Coordinate coord(je->phi(), je->eta());
      const int crate = converter.jepCrateOverlap(coord);
      if (crate < ncrates) crateColl[crate]->push_back(je);
    }
  } else if (elements) {  // take overlap from core
    iter  = elements->begin();
    iterE = elements->end();
    for (; iter != iterE; ++iter) {
      LVL1::JetElement* je = *iter;
      const LVL1::Coordinate coord(je->phi(), je->eta());
      const int crate = converter.jepCrateOverlap(coord);
      if (crate < ncrates) crateColl[crate]->push_back(je);
    }
  }
  for (int crate = 0; crate < ncrates; ++crate) {
    InternalRoiCollection* intRois = new InternalRoiCollection;
    m_jetTool->findRoIs(crateColl[crate], intRois);
    InternalRoiCollection::iterator roiIter  = intRois->begin();
    InternalRoiCollection::iterator roiIterE = intRois->end();
    for (; roiIter != roiIterE; ++roiIter) {
      LVL1::JEMRoI* roi = new LVL1::JEMRoI((*roiIter)->RoIWord());
      if (roi->crate() == crate) rois->push_back(roi);
      else delete roi;
    }
    delete intRois;
    delete crateColl[crate];
  }
}

void JEPSimBSMon::simulate(const JemRoiCollection* rois,
                                 JemHitsCollection* hits)
{
  m_log << MSG::DEBUG << "Simulate JEM Hits from JEM RoIs" << endreq;

  m_jepHitsTool->formJEMHits(rois, hits);
}

void JEPSimBSMon::simulate(const CmmJetHitsCollection* hitsIn,
                                 CmmJetHitsCollection* hitsOut,
				 int selection)
{
  m_log << MSG::DEBUG << "Simulate CMM-Jet Hit sums from CMM-Jet Hits"
        << endreq;

  if (selection == LVL1::CMMJetHits::LOCAL_MAIN) {
    m_jepHitsTool->formCMMJetHitsCrate(hitsIn, hitsOut);
  } else if (selection == LVL1::CMMJetHits::TOTAL_MAIN) {
    m_jepHitsTool->formCMMJetHitsSystem(hitsIn, hitsOut);
  } else if (selection == LVL1::CMMJetHits::ET_MAP) {
    m_jepHitsTool->formCMMJetHitsEtMap(hitsIn, hitsOut);
  }
}

void JEPSimBSMon::simulate(const JetElementCollection* elements,
                                 JemEtSumsCollection* sums)
{
  m_log << MSG::DEBUG << "Simulate JEM EtSums from JetElements" << endreq;

  m_etSumsTool->formJEMEtSums(elements, sums);
}

void JEPSimBSMon::simulate(const CmmEtSumsCollection* sumsIn,
                                 CmmEtSumsCollection* sumsOut,
				 int selection)
{
  m_log << MSG::DEBUG << "Simulate CMM-Energy Total sums from CMM-Energy Sums"
        << endreq;

  if (selection == LVL1::CMMEtSums::LOCAL) {
    m_etSumsTool->formCMMEtSumsCrate(sumsIn, sumsOut);
  } else if (selection == LVL1::CMMEtSums::TOTAL) {
    m_etSumsTool->formCMMEtSumsSystem(sumsIn, sumsOut);
  } else if (selection == LVL1::CMMEtSums::SUM_ET_MAP) {
    m_etSumsTool->formCMMEtSumsEtMaps(sumsIn, sumsOut);
  }
}
