// ********************************************************************
//
// NAME:     TrigT1CaloCpmMonTool.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Peter Faulkner
//           
//
// ********************************************************************

#include <iomanip>
#include <numeric>
#include <sstream>
#include <utility>

#include "TH1F.h"
#include "TH2F.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "GaudiKernel/ITHistSvc.h"
#include "StoreGate/StoreGateSvc.h"
#include "SGTools/StlVectorClids.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigT1CaloEvent/CMMCPHits.h"
#include "TrigT1CaloEvent/CPMHits.h"
#include "TrigT1CaloEvent/CPMTower.h"
#include "TrigT1CaloEvent/CPMRoI.h"
#include "TrigT1CaloEvent/TriggerTower.h"
#include "TrigT1CaloUtils/CoordToHardware.h"
#include "TrigT1CaloUtils/DataError.h"
#include "TrigT1CaloUtils/TriggerTowerKey.h"
#include "TrigT1Interfaces/CoordinateRange.h"
#include "TrigT1Interfaces/CPRoIDecoder.h"
#include "TrigT1Interfaces/TrigT1CaloDefs.h"

#include "TrigT1CaloMonitoring/TrigT1CaloCpmMonTool.h"
#include "TrigT1CaloMonitoring/TrigT1CaloMonErrorTool.h"
#include "TrigT1CaloMonitoring/TrigT1CaloHistogramTool.h"

const int TrigT1CaloCpmMonTool::s_crates;
const int TrigT1CaloCpmMonTool::s_modules;
const int TrigT1CaloCpmMonTool::s_maxSlices;
const int TrigT1CaloCpmMonTool::s_thresholds;
const int TrigT1CaloCpmMonTool::s_threshBits;
const int TrigT1CaloCpmMonTool::s_threshMask;

/*---------------------------------------------------------*/
TrigT1CaloCpmMonTool::TrigT1CaloCpmMonTool(const std::string & type, 
				           const std::string & name,
				           const IInterface* parent)
  : ManagedMonitorToolBase(type, name, parent),
    m_storeGate("StoreGateSvc", name),
    m_errorTool("TrigT1CaloMonErrorTool"),
    m_histTool("TrigT1CaloHistogramTool"),
    m_log(msgSvc(), name), m_events(0)
/*---------------------------------------------------------*/
{
  declareInterface<IMonitorToolBase>(this); 

  declareProperty("CPMTowerLocation",
                 m_cpmTowerLocation  = LVL1::TrigT1CaloDefs::CPMTowerLocation);
  declareProperty("CPMTowerLocationOverlap",
                 m_cpmTowerLocationOverlap =
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

  declareProperty("RootDirectory", m_rootDir = "L1Calo");
  declareProperty("MaxEnergyRange", m_maxEnergyRange = 256);

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
  m_log.setLevel(outputLevel());

  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;
  
  sc = m_storeGate.retrieve();
  if( sc.isFailure() ) {
    m_log << MSG::ERROR << "Unable to locate Service StoreGateSvc" << endreq;
    return sc;
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

  m_log << MSG::INFO << "TrigT1CaloCpmMonTool initialised" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloCpmMonTool:: finalize()
/*---------------------------------------------------------*/
{
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloCpmMonTool::bookHistograms(bool isNewEventsBlock,
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

  if ( isNewRun ) {

  std::string dir1(m_rootDir + "/CPM");
  MonGroup monShift( this, dir1 + "/Errors/Hardware", shift, run );
  MonGroup monExpert( this, dir1 + "/Errors/Hardware", expert, run );
  MonGroup monDetail( this, dir1 + "/Errors/Hardware/Detail", expert, run );
  MonGroup monCPMin( this, dir1 + "/Input", expert, run );
  MonGroup monRoIs( this, dir1 + "/Output/RoI", expert, run );
  MonGroup monCPMout( this, dir1 + "/Output/Thresholds", expert, run);
  std::string dir2(m_rootDir + "/CPM_CMM");
  MonGroup monCMM( this, dir2 + "/Errors/Hardware",  expert, run );
  MonGroup monCMMin( this, dir2 + "/Input",  expert, run );
  MonGroup monCMMout( this, dir2 + "/Output",  expert, run );

  //  Timeslice checks

  m_histTool->setMonGroup(&monCPMin);

  const int xbins = s_crates*s_maxSlices;
  m_h_CPM_slices = m_histTool->book2F("cpm_2d_tt_Slices",
       "CPM Slices and Triggered Slice;Crate/Number of Slices;Triggered Slice",
        xbins, 0, xbins, s_maxSlices, 0, s_maxSlices);
  m_histTool->numberPairs(m_h_CPM_slices, 0, s_crates-1, 1, s_maxSlices);
  m_histTool->numbers(m_h_CPM_slices, 0, s_maxSlices-1, 1, 0, false);
  m_h_PP_CP_slice = m_histTool->book2F("cpm_2d_tt_SliceMatch",
        "PPM/CPM Tower Slice Match;Crate/CPM Slice;PPM Slice",
        xbins, 0, xbins, s_maxSlices, 0, s_maxSlices);
  m_histTool->numberPairs(m_h_PP_CP_slice, 0, s_crates-1, 0, s_maxSlices-1);
  m_histTool->numbers(m_h_PP_CP_slice, 0, s_maxSlices-1, 1, 0, false);

  m_histTool->setMonGroup(&monCMMin);

  m_h_CMM_slices = m_histTool->book2F("cmm_2d_thresh_Slices",
       "CMM Slices and Triggered Slice;Crate/Number of Slices;Triggered Slice",
        xbins, 0, xbins, s_maxSlices, 0, s_maxSlices);
  m_histTool->numberPairs(m_h_CMM_slices, 0, s_crates-1, 1, s_maxSlices);
  m_histTool->numbers(m_h_CMM_slices, 0, s_maxSlices-1, 1, 0, false);
  m_h_CP_CM_slice = m_histTool->book2F("cmm_2d_thresh_SliceMatch",
        "CPM/CMM Hits Slice Match;Crate/CMM Slice;CPM Slice",
	xbins, 0, xbins, s_maxSlices, 0, s_maxSlices);
  m_histTool->numberPairs(m_h_CP_CM_slice, 0, s_crates-1, 0, s_maxSlices-1);
  m_histTool->numbers(m_h_CP_CM_slice, 0, s_maxSlices-1, 1, 0, false);

  //  CPM Tower - Trigger Tower Histos

  m_histTool->setMonGroup(&monCPMin);

  m_h_TT_Em_eta_phi = m_histTool->bookCPMEtaVsPhi(
    "ppm_em_2d_etaPhi_tt_Hitmap", "PPM Trigger Tower EM eta/phi");
  m_h_TT_Had_eta_phi = m_histTool->bookCPMEtaVsPhi(
    "ppm_had_2d_etaPhi_tt_Hitmap", "PPM Trigger Tower HAD eta/phi");

  m_h_CT_Em_Et = m_histTool->book1F("cpm_em_1d_tt_Et", "CPM Tower EM Et",
                                    m_maxEnergyRange, 0, m_maxEnergyRange);
  m_h_CT_Had_Et = m_histTool->book1F("cpm_had_1d_tt_Et", "CPM Tower HAD Et",
                                    m_maxEnergyRange, 0, m_maxEnergyRange);
  m_h_CT_Em_eta = m_histTool->book1F("cpm_em_1d_tt_Eta",
                                     "CPM Tower EM eta",50,-2.5,2.5);
  m_h_CT_Had_eta = m_histTool->book1F("cpm_had_1d_tt_Eta",
                                      "CPM Tower HAD eta",50,-2.5,2.5);
  m_h_CT_Em_phi = m_histTool->book1F("cpm_em_1d_tt_Phi",
                                     "CPM Tower EM phi",64,0.,2.*M_PI);
  m_h_CT_Had_phi = m_histTool->book1F("cpm_had_1d_tt_Phi",
                                      "CPM Tower HAD phi", 64,0.,2.*M_PI);
  m_h_CT_Em_eta_phi = m_histTool->bookCPMEtaVsPhi(
    "cpm_em_2d_etaPhi_tt_Hitmap", "CPM Tower EM eta/phi");
  m_h_CT_Had_eta_phi = m_histTool->bookCPMEtaVsPhi(
    "cpm_had_2d_etaPhi_tt_Hitmap", "CPM Tower HAD eta/phi");
  m_h_CT_Em_eta_phi_w = m_histTool->bookCPMEtaVsPhi(
    "cpm_em_2d_etaPhi_tt_EtWeighted", "CPM Tower EM eta/phi weighted");
  m_h_CT_Had_eta_phi_w = m_histTool->bookCPMEtaVsPhi(
    "cpm_had_2d_etaPhi_tt_EtWeighted", "CPM Tower HAD eta/phi weighted");

  //  CPM Tower error bits

  m_histTool->setMonGroup(&monDetail);

  m_h_CT_Em_parity = m_histTool->bookCPMEtaVsPhi(
    "cpm_em_2d_etaPhi_tt_Parity", "CPM Tower EM Parity Errors");
  m_h_CT_Had_parity = m_histTool->bookCPMEtaVsPhi(
    "cpm_had_2d_etaPhi_tt_Parity", "CPM Tower HAD Parity Errors");
  m_h_CT_Em_link = m_histTool->bookCPMEtaVsPhi(
    "cpm_em_2d_etaPhi_tt_LinkDown", "CPM Tower EM Link Down Errors");
  m_h_CT_Had_link = m_histTool->bookCPMEtaVsPhi(
    "cpm_had_2d_etaPhi_tt_LinkDown", "CPM Tower HAD Link Down Errors");
  m_h_CT_status = m_histTool->bookCPMSubStatusVsCrateModule(
                               "cpm_2d_Status", "CPM Sub-status bits");

  //  CPM RoIs

  m_histTool->setMonGroup(&monRoIs);

  m_h_RoI_thresholds = m_histTool->bookCPMCrateModuleVsThreshold(
                       "cpm_2d_roi_Thresholds", "CPM RoI Thresholds");
  m_h_RoI_eta_phi = m_histTool->bookCPMRoIEtaVsPhi(
    "cpm_2d_etaPhi_roi_Hitmap", "CPM RoIs Eta-Phi Hit Map");
  m_h_RoI_Em_eta_phi = m_histTool->bookCPMRoIEtaVsPhi(
    "cpm_2d_etaPhi_roi_EmHitmap", "CPM RoIs EM Eta-Phi Hit Map");
  m_h_RoI_Tau_eta_phi = m_histTool->bookCPMRoIEtaVsPhi(
    "cpm_2d_etaPhi_roi_TauHitmap", "CPM RoIs Tau Eta-Phi Hit Map");
  m_h_RoI_Saturation = m_histTool->bookCPMRoIEtaVsPhi(
    "cpm_2d_etaPhi_roi_Saturation", "CPM RoI Tower Saturation");

  m_histTool->setMonGroup(&monDetail);

  m_h_RoI_Parity = m_histTool->bookCPMRoIEtaVsPhi("cpm_2d_etaPhi_roi_Parity",
                                                   "CPM RoI Parity Errors");

  //  CPM Hits
  
  m_histTool->setMonGroup(&monCPMout);

  m_h_CPM_thresholds = m_histTool->bookCPMCrateModuleVsThreshold(
    "cpm_2d_thresh_Weighted", "CPM Hits Thresholds Weighted");

  //  CMM-CP Hits

  m_histTool->setMonGroup(&monCMMin);

  m_h_CMM_thresholds = m_histTool->bookCPMCrateModuleVsThreshold(
    "cmm_2d_thresh_Weighted", "CMM-CP Hits Thresholds Weighted");

  m_histTool->setMonGroup(&monCMMout);

  m_h_CMM_T_thresholds = m_histTool->bookCPMSumVsThreshold(
    "cmm_2d_thresh_SumsWeighted", "CMM-CP Hit Sums Thresholds Weighted");

  //  CMM error bits

  m_histTool->setMonGroup(&monCMM);

  m_h_CMM_parity = m_histTool->book2F("cmm_2d_thresh_Parity",
                          "CMM Parity Errors;Module or Remote;Crate/CMM",
		           15, 1, 16, 8, 0, 8);
  m_histTool->numbers(m_h_CMM_parity, 1, s_modules);
  m_h_CMM_parity->GetXaxis()->SetBinLabel(15, "REM");
  m_histTool->numberPairs(m_h_CMM_parity, 0, s_crates-1, 0, 1, 1, 0, false);
  m_h_CMM_status = m_histTool->book2F("cmm_2d_Status",
                          "CMM Sub-status bits;;Crate/CMM",
                           8, 0., 8., 8, 0., 8.);
  m_histTool->subStatus(m_h_CMM_status);
  m_histTool->numberPairs(m_h_CMM_status, 0, s_crates-1, 0, 1, 1, 0, false);
  
  //  Error Overview and Summary

  m_histTool->setMonGroup(&monExpert);

  m_h_CP_overview = m_histTool->book2F("cpm_2d_ErrorOverview",
                           "CP Error Overview", 64, 0, 64,
			    NumberOfSummaryBins, 0, NumberOfSummaryBins);
  m_histTool->cpmCMMCrateModule(m_h_CP_overview);

  m_histTool->setMonGroup(&monShift);

  m_h_CP_errors = m_histTool->book1F("cpm_1d_ErrorSummary",
                         "CP Error Summary;;Events",
                          NumberOfSummaryBins, 0, NumberOfSummaryBins);

  TH1*   hist = m_h_CP_overview;
  TAxis* axis = hist->GetYaxis();
  for (int i = 0; i < 2; ++i) {
    axis->SetBinLabel(1+EMParity,  "EM parity");
    axis->SetBinLabel(1+EMLink,    "EM link");
    axis->SetBinLabel(1+HadParity, "Had parity");
    axis->SetBinLabel(1+HadLink,   "Had link");
    axis->SetBinLabel(1+CPMStatus, "CPM status");
    axis->SetBinLabel(1+RoIParity, "RoI parity");
    axis->SetBinLabel(1+CMMParity, "CMM parity");
    axis->SetBinLabel(1+CMMStatus, "CMM status");
    //axis->SetLabelSize(0.06);
    hist = m_h_CP_errors;
    axis = hist->GetXaxis();
  }

  m_histTool->unsetMonGroup();

  m_events = 0;

  } // end if (isNewRun ...

  m_log << MSG::DEBUG << "Leaving bookHistograms" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloCpmMonTool::fillHistograms()
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "fillHistograms entered" << endreq;

  // Skip events believed to be corrupt

  if (m_errorTool->corrupt()) {
    m_log << MSG::DEBUG << "Skipping corrupt event" << endreq;
    return StatusCode::SUCCESS;
  }

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
  const CpmTowerCollection* cpmTowerOverlapTES = 0; 
  if (m_storeGate->contains<CpmTowerCollection>(m_cpmTowerLocationOverlap)) {
    sc = m_storeGate->retrieve(cpmTowerOverlapTES, m_cpmTowerLocationOverlap); 
  } else sc = StatusCode::FAILURE;
  if( sc.isFailure()  ||  !cpmTowerOverlapTES ) {
    m_log << MSG::DEBUG<< "No Overlap CPM Tower container found"<< endreq; 
  }
  
  //Retrieve CPM RoIs from SG
  const CpmRoiCollection* cpmRoiTES = 0;
  sc = m_storeGate->retrieve( cpmRoiTES, m_cpmRoiLocation);
  if( sc.isFailure()  ||  !cpmRoiTES ) {
    m_log << MSG::DEBUG << "No CPM RoIs container found"<< endreq;
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

  // Vectors for error overview bits;
  std::vector<int> errorsCPM(s_crates*s_modules);
  std::vector<int> errorsCMM(s_crates*2); // L/R

  //=============================================
  //   CPM Tower - Trigger Tower comparison plots
  //=============================================

  // Maps for slice match
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
      LVL1::TriggerTower* tt = *ttIterator;
      const double eta = tt->eta();
      if (eta < -2.5 || eta > 2.5) continue;
      const std::vector<int>& emLut(tt->emLUT());
      const std::vector<int>& hadLut(tt->hadLUT());
      if (std::accumulate(emLut.begin(), emLut.end(), 0) == 0 &&
          std::accumulate(hadLut.begin(), hadLut.end(), 0) == 0) continue;
      const int    em  = tt->emEnergy();
      const int    had = tt->hadEnergy();
      const double phi = tt->phi();
      if (em)  m_histTool->fillCPMEtaVsPhi(m_h_TT_Em_eta_phi, eta, phi);
      if (had) m_histTool->fillCPMEtaVsPhi(m_h_TT_Had_eta_phi, eta, phi);
      const unsigned int key = towerKey.ttKey(phi, eta);
      if (key > maxKey) maxKey = key;
      ttMap.insert(std::make_pair(key, tt));
    }
  }

  for (int i = 0; i < 2; ++i) {
    const bool core = (i == 0);
    const CpmTowerCollection* cpmTwrTES = (core) ? cpmTowerTES
                                                 : cpmTowerOverlapTES;
    if (cpmTwrTES) {
      CpmTowerCollection::const_iterator ctIterator    = cpmTwrTES->begin(); 
      CpmTowerCollection::const_iterator ctIteratorEnd = cpmTwrTES->end(); 

      for (; ctIterator != ctIteratorEnd; ++ctIterator) {
	LVL1::CPMTower* ct = *ctIterator;
        const int    em  = ct->emEnergy();
        const int    had = ct->hadEnergy();
        const double eta = ct->eta();
        const double phi = ct->phi();
        const LVL1::Coordinate coord(phi, eta);
        LVL1::CoordToHardware converter;
        const int crate  = (core) ? converter.cpCrate(coord)
	                          : converter.cpCrateOverlap(coord);
        const int cpm    = (core) ? converter.cpModule(coord)
	                          : converter.cpModuleOverlap(coord);
        const int loc    = crate * s_modules + cpm - 1;
        const int peak   = ct->peak();
        const int slices = (ct->emEnergyVec()).size();
        m_h_CPM_slices->Fill(crate*s_maxSlices + slices - 1, peak, 1.);
        if (em && core) {
          m_h_CT_Em_Et->Fill(em, 1.);
          m_h_CT_Em_eta->Fill(eta, 1.);
          m_h_CT_Em_phi->Fill(phi, 1.);
	  m_histTool->fillCPMEtaVsPhi(m_h_CT_Em_eta_phi, eta, phi);
	  m_histTool->fillCPMEtaVsPhi(m_h_CT_Em_eta_phi_w, eta, phi, em);
        }
        if (had && core) {
          m_h_CT_Had_Et->Fill(had, 1.);
          m_h_CT_Had_eta->Fill(eta, 1.);
          m_h_CT_Had_phi->Fill(phi, 1.);
	  m_histTool->fillCPMEtaVsPhi(m_h_CT_Had_eta_phi, eta, phi);
	  m_histTool->fillCPMEtaVsPhi(m_h_CT_Had_eta_phi_w, eta, phi, had);
        }
        // Errors
	int error = ct->emError();
	if (error) {
	  const LVL1::DataError emError(error);
	  if (emError.get(LVL1::DataError::Parity)) {
	    m_histTool->fillCPMEtaVsPhi(m_h_CT_Em_parity, eta, phi);
	    errorsCPM[loc] |= (1 << EMParity);
	  }
          if (emError.get(LVL1::DataError::LinkDown)) {
	    m_histTool->fillCPMEtaVsPhi(m_h_CT_Em_link, eta, phi);
	    errorsCPM[loc] |= (1 << EMLink);
          }
	  const int status = error >> LVL1::DataError::GLinkParity;
	  if (status) {
	    for (int bit = 0; bit < 8; ++bit) {
	      if ((status >> bit) & 0x1) m_h_CT_status->Fill(bit, loc);
            }
	    errorsCPM[loc] |= (1 << CPMStatus);
          }
        }
	error = ct->hadError();
	if (error) {
	  const LVL1::DataError hadError(error);
	  if (hadError.get(LVL1::DataError::Parity)) {
	    m_histTool->fillCPMEtaVsPhi(m_h_CT_Had_parity, eta, phi);
	    errorsCPM[loc] |= (1 << HadParity);
	  }
          if (hadError.get(LVL1::DataError::LinkDown)) {
	    m_histTool->fillCPMEtaVsPhi(m_h_CT_Had_link, eta, phi);
	    errorsCPM[loc] |= (1 << HadLink);
          }
        }

        if (core) {
          const unsigned int key = towerKey.ttKey(phi, eta);
          if (key > maxKey) maxKey = key;
          cpMap.insert(std::make_pair(key, ct));
        }
      }
    }
  }

  // Slice match

  TriggerTowerMap::const_iterator ttMapIter    = ttMap.begin();
  TriggerTowerMap::const_iterator ttMapIterEnd = ttMap.end();
  CpmTowerMap::const_iterator     cpMapIter    = cpMap.begin();
  CpmTowerMap::const_iterator     cpMapIterEnd = cpMap.end();

  while (ttMapIter != ttMapIterEnd && cpMapIter != cpMapIterEnd) {

    unsigned int ttKey = ttMapIter->first;
    unsigned int cpKey = cpMapIter->first;

    if      (cpKey > ttKey) ++ttMapIter;
    else if (ttKey > cpKey) ++cpMapIter;

    else {
      const LVL1::TriggerTower* tt = ttMapIter->second;
      const LVL1::CPMTower*     cp = cpMapIter->second;
      const std::vector<int>& emLut(tt->emLUT());
      const std::vector<int>& hadLut(tt->hadLUT());
      const std::vector<int>& emVec(cp->emEnergyVec());
      const std::vector<int>& hadVec(cp->hadEnergyVec());
      const int sliceEmLut = emLut.size();
      const int sliceHadLut = hadLut.size();
      const int sliceEmVec = emVec.size();
      const int sliceHadVec = hadVec.size();
      const int crate = static_cast<int>(tt->phi()/(M_PI/2.));
      for (int slice = 0; slice < sliceEmLut; ++slice) {
        if (emLut[slice] > 0) {
	  for (int slice2 = 0; slice2 < sliceEmVec; ++slice2) {
	    if (emLut[slice] == emVec[slice2]) {
	      m_h_PP_CP_slice->Fill(crate*s_maxSlices+slice2, slice, 1.);
            }
          }
        }
      }
      for (int slice = 0; slice < sliceHadLut; ++slice) {
        if (hadLut[slice] > 0) {
	  for (int slice2 = 0; slice2 < sliceHadVec; ++slice2) {
	    if (hadLut[slice] == hadVec[slice2]) {
	      m_h_PP_CP_slice->Fill(crate*s_maxSlices+slice2, slice, 1.);
            }
          }
        }
      }
      ++ttMapIter;
      ++cpMapIter;
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
      int bin  = crate * s_modules + cpm - 1;
      for (int thresh = 0; thresh < s_thresholds; ++thresh) {
        const int hit = (hits >> thresh) & 0x1;
        if (hit) {
          m_h_RoI_thresholds->Fill(bin, thresh, 1.);
	  m_histTool->fillCPMRoIEtaVsPhi(m_h_RoI_eta_phi, eta, phi);
	  if (thresh < s_thresholds/2) {
	    m_histTool->fillCPMRoIEtaVsPhi(m_h_RoI_Em_eta_phi, eta, phi);
	  } else {
	    m_histTool->fillCPMRoIEtaVsPhi(m_h_RoI_Tau_eta_phi, eta, phi);
          }
        }
      }
      const LVL1::DataError err((*crIterator)->error());
      if (err.get(LVL1::DataError::Overflow)) {
	m_histTool->fillCPMRoIEtaVsPhi(m_h_RoI_Saturation, eta, phi);
      }
      if (err.get(LVL1::DataError::Parity)) {
	m_histTool->fillCPMRoIEtaVsPhi(m_h_RoI_Parity, eta, phi);
	errorsCPM[bin] |= (1 << RoIParity);
      }
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
      const unsigned int hits0 = (*chIterator)->HitWord0();
      const unsigned int hits1 = (*chIterator)->HitWord1();
      const int crate = (*chIterator)->crate();
      const int cpm   = (*chIterator)->module();
      const int bin   = crate * s_modules + cpm - 1;
      const int peak   = (*chIterator)->peak();
      const int slices = ((*chIterator)->HitsVec0()).size();
      m_h_CPM_slices->Fill(crate*s_maxSlices + slices -1, peak, 1.);
      for (int thresh = 0; thresh < s_thresholds/2; ++thresh) {
        int hit0 = (hits0 >> thresh*s_threshBits) & s_threshMask;
        int hit1 = (hits1 >> thresh*s_threshBits) & s_threshMask;
        if (hit0) m_h_CPM_thresholds->Fill(bin, thresh,   hit0);
        if (hit1) m_h_CPM_thresholds->Fill(bin, thresh+8, hit1);
      }
      const unsigned int key = crate * s_modules + cpm;
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
      const int peak   = (*cmIterator)->peak();
      const int slices = ((*cmIterator)->HitsVec0()).size();
      m_h_CMM_slices->Fill(crate*s_maxSlices + slices -1, peak, 1.);
      int bin = 0;
      if (dataId <= s_modules) bin = crate*s_modules + dataId - 1;
      else {
	bin = crate;
	if (dataId == LVL1::CMMCPHits::REMOTE_0) bin = s_crates;
	if (dataId == LVL1::CMMCPHits::REMOTE_1) bin = s_crates + 1;
	if (dataId == LVL1::CMMCPHits::REMOTE_2) bin = s_crates + 2;
	if (dataId == LVL1::CMMCPHits::TOTAL)    bin = s_crates + 3;
      }
      for (int thresh = 0; thresh < s_thresholds/2; ++thresh) {
        int hit0 = (hits0 >> thresh*s_threshBits) & s_threshMask;
        int hit1 = (hits1 >> thresh*s_threshBits) & s_threshMask;
        if (dataId <= s_modules) {
          if (hit0) m_h_CMM_thresholds->Fill(bin, thresh,   hit0);
          if (hit1) m_h_CMM_thresholds->Fill(bin, thresh+8, hit1);
        } else {
          if (hit0) m_h_CMM_T_thresholds->Fill(bin, thresh,   hit0);
          if (hit1) m_h_CMM_T_thresholds->Fill(bin, thresh+8, hit1);
        }
      }
      // Errors
      int error0 = (*cmIterator)->Error0();
      int error1 = (*cmIterator)->Error1();
      if (error0 || error1) {
        const LVL1::DataError hit0Err(error0);
        const LVL1::DataError hit1Err(error1);
	const int parity0 = hit0Err.get(LVL1::DataError::Parity);
	const int parity1 = hit1Err.get(LVL1::DataError::Parity);
	if (parity0 || parity1) {
          if (dataId <= s_modules) {
	    if (parity1) m_h_CMM_parity->Fill(dataId, 2*crate);
	    if (parity0) m_h_CMM_parity->Fill(dataId, 2*crate + 1);
	    errorsCPM[bin] |= (1 << CMMParity);
          } else {
  	    int remBin   = s_modules + 1;
	    int remCrate = -1;
            if      (dataId == LVL1::CMMCPHits::REMOTE_0) remCrate = 0;
            else if (dataId == LVL1::CMMCPHits::REMOTE_1) remCrate = 1;
            else if (dataId == LVL1::CMMCPHits::REMOTE_2) remCrate = 2;
	    if (remCrate >= 0) {
	      if (parity1) m_h_CMM_parity->Fill(remBin, 2*remCrate);
	      if (parity0) m_h_CMM_parity->Fill(remBin, 2*remCrate + 1);
	    }
          }
	  if (parity1) errorsCMM[crate*2]     |= (1 << CMMParity);
	  if (parity0) errorsCMM[crate*2 + 1] |= (1 << CMMParity);
        }
        // Sub-status errors
        const int status0 = error0 >> LVL1::DataError::GLinkParity;
        const int status1 = error1 >> LVL1::DataError::GLinkParity;
	if (status0 || status1) {
          for (int bit = 0; bit < 8; ++bit) {
	    if ((status1 >> bit) & 0x1) m_h_CMM_status->Fill(bit, 2*crate);
	    if ((status0 >> bit) & 0x1) m_h_CMM_status->Fill(bit, 2*crate + 1);
	  }
	  if (status1) errorsCMM[crate*2]     |= (1 << CMMStatus);
	  if (status0) errorsCMM[crate*2 + 1] |= (1 << CMMStatus);
        }
      }

      if (dataId <= s_modules) {
        const unsigned int key = crate * s_modules + dataId;
        cmmMap.insert(std::make_pair(key, *cmIterator));
      }
    }
  }

  // Slice match

  CpmHitsMap::const_iterator   cpmMapIter    = cpmMap.begin();
  CpmHitsMap::const_iterator   cpmMapIterEnd = cpmMap.end();
  CmmCpHitsMap::const_iterator cmmMapIter    = cmmMap.begin();
  CmmCpHitsMap::const_iterator cmmMapIterEnd = cmmMap.end();

  while (cpmMapIter != cpmMapIterEnd && cmmMapIter != cmmMapIterEnd) {

    unsigned int cpmKey = cpmMapIter->first;
    unsigned int cmmKey = cmmMapIter->first;

    if      (cmmKey > cpmKey) ++cpmMapIter;
    else if (cpmKey > cmmKey) ++cmmMapIter;

    else {
      const LVL1::CPMHits*   cpmh = cpmMapIter->second;
      const LVL1::CMMCPHits* cmmh = cmmMapIter->second;
      const std::vector<unsigned int>& cpmVec0(cpmh->HitsVec0());
      const std::vector<unsigned int>& cpmVec1(cpmh->HitsVec1());
      const std::vector<unsigned int>& cmmVec0(cmmh->HitsVec0());
      const std::vector<unsigned int>& cmmVec1(cmmh->HitsVec1());
      const int sliceCpmVec0 = cpmVec0.size();
      const int sliceCpmVec1 = cpmVec1.size();
      const int sliceCmmVec0 = cmmVec0.size();
      const int sliceCmmVec1 = cmmVec1.size();
      const int crate = cpmh->crate();
      if (sliceCpmVec0 == sliceCpmVec1 && sliceCmmVec0 == sliceCmmVec1) {
        for (int slice = 0; slice < sliceCpmVec0; ++slice) {
          if (cpmVec0[slice] > 0 || cpmVec1[slice] > 0) {
	    for (int slice2 = 0; slice2 < sliceCmmVec0; ++slice2) {
	      if (cpmVec0[slice] == cmmVec0[slice2] &&
	          cpmVec1[slice] == cmmVec1[slice2]) {
		m_h_CP_CM_slice->Fill(crate*s_maxSlices+slice2, slice, 1.);
	      }
            }
          }
        }
      }
      ++cpmMapIter;
      ++cmmMapIter;
    }
  }

  // Update error summary plot

  std::vector<int> crateErr(4);
  for (int err = 0; err < NumberOfSummaryBins; ++err) {
    int error = 0;
    for (int loc = 0; loc < s_crates*s_modules; ++loc) {
      if ((errorsCPM[loc] >> err) & 0x1) {
        m_h_CP_overview->Fill(loc, err, 1.);
	error = 1;
	crateErr[loc/s_modules] |= (1 << err);
      }
      if (loc < s_crates*2) {
        if ((errorsCMM[loc] >> err) & 0x1) {
          m_h_CP_overview->Fill(loc+s_crates*s_modules, err, 1.);
	  error = 1;
	  crateErr[loc/2] |= (1 << err);
        }
      }
    }
    if (error) m_h_CP_errors->Fill(err);
  }
  ++m_events;

  // Save error vector for global summary

  std::vector<int>* save = new std::vector<int>(crateErr);
  sc = m_storeGate->record(save, "L1CaloCPMErrorVector");
  if (sc != StatusCode::SUCCESS) {
    m_log << MSG::ERROR << "Error recording CPM error vector in TES " << endreq;
    return sc;
  }

  m_log << MSG::DEBUG << "Leaving fillHistograms" << endreq;

  return StatusCode::SUCCESS;

}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloCpmMonTool::procHistograms(bool isEndOfEventsBlock,
                                  bool isEndOfLumiBlock, bool isEndOfRun)
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "procHistograms entered" << endreq;

  if (isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun) {
  }

  return StatusCode::SUCCESS;
}
