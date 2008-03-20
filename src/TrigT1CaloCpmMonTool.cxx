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

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "StoreGate/StoreGateSvc.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigT1Calo/CMMCPHits.h"
#include "TrigT1Calo/CoordToHardware.h"
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

const int TrigT1CaloCpmMonTool::s_crates;
const int TrigT1CaloCpmMonTool::s_modules;
const int TrigT1CaloCpmMonTool::s_thresholds;
const int TrigT1CaloCpmMonTool::s_threshBits;
const int TrigT1CaloCpmMonTool::s_threshMask;

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
  declareProperty("CPMTowerLocationOverlap",
                 m_cpmTowerLocationOverlap =
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
  declareProperty("SingleDirectory", m_oneDir = false);
  declareProperty("PhiUnits", m_phiUnits = "channels",
                  "Phi Units: radians, degrees or channels");
  declareProperty("NoiseSignalSplit", m_noiseSignalSplit = 0);
  declareProperty("MaxEnergyRange", m_maxEnergyRange = 50);
  declareProperty( "Offline", m_Offline = 1) ;

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

  if ( isNewEventsBlock || isNewLumiBlock ) { }

  if ( isNewRun ) {
  
  if (m_oneDir) newGroup(cpmDir, shift, run );

  //  Timeslice checks

  newGroup(cpmDir + "_slices", expert, run );

  for (int crate = 0; crate < s_crates; ++crate) {
    std::ostringstream cnum;
    cnum << crate;
    std::string name("CPM_slices_crate_" + cnum.str());
    std::string title("CPM Timeslices Crate " + cnum.str());
    TH2F* hist = book2F(name, title + ";Number of Slices;Triggered Slice",
			6, -0.5, 5.5, 5, -0.5, 4.5);
    m_v_CPM_slices.push_back(hist);
    name = "PP_CP_slice_crate_" + cnum.str();
    title = "PPr/CPM Tower Slice Match Crate " + cnum.str();
    hist = book2F(name, title + ";PPr Slice;CPM Slice",
                        5, -0.5, 4.5, 5, -0.5, 4.5); 
    m_v_PP_CP_slice.push_back(hist);
  }

  newGroup(cmmDir + "_slices", expert, run );

  for (int crate = 0; crate < s_crates; ++crate) {
    std::ostringstream cnum;
    cnum << crate;
    std::string name("CMM_slices_crate_" + cnum.str());
    std::string title("CMM Timeslices Crate " + cnum.str());
    TH2F* hist = book2F(name, title + ";Number of Slices;Triggered Slice",
			6, -0.5, 5.5, 5, -0.5, 4.5);
    m_v_CMM_slices.push_back(hist);
    name = "CP_CM_slice_crate_" + cnum.str();
    title = "CPM/CMM Hits Slice Match Crate " + cnum.str();
    hist = book2F(name, title + ";CPM Slice;CMM Slice",
                        5, -0.5, 4.5, 5, -0.5, 4.5); 
    m_v_CP_CM_slice.push_back(hist);
  }

  //  CPM Tower - Trigger Tower comparison Histos

  newGroup(pprDir + "_Towers", expert, run );

  const int signalBins = m_maxEnergyRange - m_noiseSignalSplit;
  if (m_noiseSignalSplit) {
    m_h_TT_Em_Et = book1F("TT_EM_Et","Trigger Tower EM Et Noise",
                                   m_noiseSignalSplit, 0, m_noiseSignalSplit);
    m_h_TT_Had_Et = book1F("TT_HAD_Et","Trigger Tower HAD Et Noise",
                                   m_noiseSignalSplit, 0, m_noiseSignalSplit);
    m_h_TT_Em_Et_s = book1F("TT_EM_Et_s","Trigger Tower EM Et Signal",
                            signalBins, m_noiseSignalSplit, m_maxEnergyRange);
    m_h_TT_Had_Et_s = book1F("TT_HAD_Et_s","Trigger Tower HAD Et Signal",
                            signalBins, m_noiseSignalSplit, m_maxEnergyRange);
  } else {
    m_h_TT_Em_Et = book1F("TT_EM_Et","Trigger Tower EM Et",
                            signalBins, m_noiseSignalSplit, m_maxEnergyRange);
    m_h_TT_Had_Et = book1F("TT_HAD_Et","Trigger Tower HAD Et",
                            signalBins, m_noiseSignalSplit, m_maxEnergyRange);
  }
  m_h_TT_Em_eta = book1F("TT_EM_eta","Trigger Tower EM eta",50,-2.5,2.5);
  m_h_TT_Had_eta = book1F("TT_HAD_eta","Trigger Tower HAD eta",50,-2.5,2.5);
  m_h_TT_Em_phi = book1F("TT_EM_phi","Trigger Tower EM phi ",64,0,m_phiMax);
  m_h_TT_Had_phi = book1F("TT_HAD_phi","Trigger Tower HAD phi ",64,0,m_phiMax);
  m_h_TT_Em_eta_phi = book2F("TT_EM_eta_phi",
         "Trigger Tower EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TT_Had_eta_phi = book2F("TT_HAD_eta_phi",
         "Trigger Tower HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TT_Em_eta_phi_w = book2F("TT_EM_eta_phi_w",
      "Trigger Tower EM eta/phi weighted;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TT_Had_eta_phi_w = book2F("TT_HAD_eta_phi_w",
      "Trigger Tower HAD eta/phi weighted;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);

  newGroup(cpmDir + "_Towers", expert, run );

  if (m_noiseSignalSplit) {
    m_h_CT_Em_Et = book1F("CT_EM_Et","CPM Tower EM Et Noise",
                                   m_noiseSignalSplit, 0, m_noiseSignalSplit);
    m_h_CT_Had_Et = book1F("CT_HAD_Et","CPM Tower HAD Et Noise",
                                   m_noiseSignalSplit, 0, m_noiseSignalSplit);
    m_h_CT_Em_Et_s = book1F("CT_EM_Et_s","CPM Tower EM Et Signal",
                            signalBins, m_noiseSignalSplit, m_maxEnergyRange);
    m_h_CT_Had_Et_s = book1F("CT_HAD_Et_s","CPM Tower HAD Et Signal",
                            signalBins, m_noiseSignalSplit, m_maxEnergyRange);
  } else {
    m_h_CT_Em_Et = book1F("CT_EM_Et","CPM Tower EM Et",
                            signalBins, m_noiseSignalSplit, m_maxEnergyRange);
    m_h_CT_Had_Et = book1F("CT_HAD_Et","CPM Tower HAD Et",
                            signalBins, m_noiseSignalSplit, m_maxEnergyRange);
  }
  m_h_CT_Em_eta = book1F("CT_EM_eta","CPM Tower EM eta",50,-2.5,2.5);
  m_h_CT_Had_eta = book1F("CT_HAD_eta","CPM Tower HAD eta",50,-2.5,2.5);
  m_h_CT_Em_phi = book1F("CT_EM_phi","CPM Tower EM phi ",64,0,m_phiMax);
  m_h_CT_Had_phi = book1F("CT_HAD_phi","CPM Tower HAD phi ",64,0,m_phiMax);
  m_h_CT_Em_eta_phi = book2F("CT_EM_eta_phi",
         "CPM Tower EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CT_Had_eta_phi = book2F("CT_HAD_eta_phi",
         "CPM Tower HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CT_Em_eta_phi_w = book2F("CT_EM_eta_phi_w",
         "CPM Tower EM eta/phi weighted;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CT_Had_eta_phi_w = book2F("CT_HAD_eta_phi_w",
         "CPM Tower HAD eta/phi weighted;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);

  //  CPM Tower error bits

  newGroup(cpErrDir + "_CPM_parity", expert, run );

  m_h_CT_Em_parity = book2F("CT_EM_parity",
            "CPM Tower EM Parity Errors;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CT_Had_parity = book2F("CT_HAD_parity",
            "CPM Tower HAD Parity Errors;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);

  newGroup(cpErrDir + "_CPM_link", expert, run );

  m_h_CT_Em_link = book2F("CT_EM_link",
            "CPM Tower EM Link Down Errors;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CT_Had_link = book2F("CT_HAD_link",
         "CPM Tower HAD Link Down Errors;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);

  newGroup(cpErrDir + "_CPM_status", expert, run );

  m_h_CT_status = book1F("CT_status", "CPM Sub-status bits", 8, 0., 8.);
  setStatusLabels(m_h_CT_status);
  m_h_CT_status_eta_phi = book2F("CT_status_eta_phi",
            "CPM Sub-status hit-map;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);

  //  Multiple slice CPM Tower plots

  //....

  //  Trigger Tower/CPM Tower event by event comparison plots

  newGroup(cpErrDir + "_PPr_CPM", expert, run );

  m_h_TTeqCT_Em_eta_phi = book2F("TTeqCT_EM_eta_phi",
    "Trigger/CPM Tower non-zero match EM eta/phi;eta;phi",
                                                 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TTneCT_Em_eta_phi = book2F("TTneCT_EM_eta_phi",
    "Trigger/CPM Tower non-zero mismatch EM eta/phi;eta;phi",
                                                 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TTnoCT_Em_eta_phi = book2F("TTnoCT_EM_eta_phi",
    "Trigger Tower/no CPM Tower EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CTnoTT_Em_eta_phi = book2F("CTnoTT_EM_eta_phi",
    "CPM Tower/no Trigger Tower EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TTeqCT_Had_eta_phi = book2F("TTeqCT_HAD_eta_phi",
    "Trigger/CPM Tower non-zero match HAD eta/phi;eta;phi",
                                                 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TTneCT_Had_eta_phi = book2F("TTneCT_HAD_eta_phi",
    "Trigger/CPM Tower non-zero mismatch HAD eta/phi;eta;phi",
                                                 50,-2.5,2.5,64,0,m_phiMax);
  m_h_TTnoCT_Had_eta_phi = book2F("TTnoCT_HAD_eta_phi",
   "Trigger Tower/no CPM Tower HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CTnoTT_Had_eta_phi = book2F("CTnoTT_HAD_eta_phi",
   "CPM Tower/no Trigger Tower HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);

  //  CPM Tower Core/Overlap event by event comparison plots

  newGroup(cpErrDir + "_Phi_overlap", expert, run );

  m_h_CTeqCO_Em_eta_phi = book2F("CTeqCO_EM_eta_phi",
    "Core/Overlap non-zero match EM eta/phi;eta;phi",
                                                50,-2.5,2.5,64,0,m_phiMax);
  m_h_CTneCO_Em_eta_phi = book2F("CTneCO_EM_eta_phi",
    "Core/Overlap non-zero mismatch EM eta/phi;eta;phi",
                                                50,-2.5,2.5,64,0,m_phiMax);
  m_h_CTnoCO_Em_eta_phi = book2F("CTnoCO_EM_eta_phi",
    "Core/no Overlap EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_COnoCT_Em_eta_phi = book2F("COnoCT_EM_eta_phi",
    "Overlap/no Core EM eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_CTeqCO_Had_eta_phi = book2F("CTeqCO_HAD_eta_phi",
    "Core/Overlap non-zero match HAD eta/phi;eta;phi",
                                                50,-2.5,2.5,64,0,m_phiMax);
  m_h_CTneCO_Had_eta_phi = book2F("CTneCO_HAD_eta_phi",
    "Core/Overlap non-zero mismatch HAD eta/phi;eta;phi",
                                                50,-2.5,2.5,64,0,m_phiMax);
  m_h_CTnoCO_Had_eta_phi = book2F("CTnoCO_HAD_eta_phi",
    "Core/no Overlap HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);
  m_h_COnoCT_Had_eta_phi = book2F("COnoCT_HAD_eta_phi",
    "Overlap/no Core HAD eta/phi;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);

  //  CPM RoIs

  newGroup(cpmDir + "_RoIs", expert, run );

  for (int thresh = 0; thresh < s_thresholds; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("RoI_Thresh_" + cnum.str());
    std::string title("RoI Threshold " + cnum.str() + ";Crate/CPM");
    TH1F* hist = book1F(name, title, 56, 0, 56);
    setThresholdLabels(hist);
    m_v_RoI_thresholds.push_back(hist);
  }
  const double halfPhiBin = m_phiMax/128.;
  for (int thresh = 0; thresh < s_thresholds; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("RoI_2D_Thresh_" + cnum.str());
    std::string title("RoI eta/phi Threshold " + cnum.str() + ";eta;phi");
    m_v_RoI_2D_thresholds.push_back(book2F(name, title, 51, -2.55, 2.55,
                                    65, -halfPhiBin, m_phiMax+halfPhiBin));
  }

  newGroup(cpErrDir + "_RoI_parity", expert, run );

  m_h_RoI_Parity = book2F("CPM_RoI_Parity",
            "CPM RoI Parity Errors;eta;phi", 50,-2.5,2.5,64,0,m_phiMax);

  //  CPM Hits
  
  newGroup(cpmDir + "_Hits", expert, run );

  for (int thresh = 0; thresh < s_thresholds; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("CH_Thresh_" + cnum.str());
    std::string title("CPM Hits Threshold " + cnum.str() + ";Crate/CPM");
    TH1F* hist = book1F(name, title, 56, 0, 56);
    setThresholdLabels(hist);
    m_v_thresholds.push_back(hist);
  }

  //  Multiple slice CPM Hit plots

  //....

  //  CMM-CP Hits

  newGroup(cmmDir + "_Hits", expert, run );

  for (int thresh = 0; thresh < s_thresholds; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("CM_Thresh_" + cnum.str());
    std::string title("CMM-CP Hits Threshold " + cnum.str() + ";Crate/CPM");
    TH1F* hist = book1F(name, title, 56, 0, 56);
    setThresholdLabels(hist);
    m_v_CMM_thresholds.push_back(hist);
  }
  for (int thresh = 0; thresh < s_thresholds; ++thresh) {
    std::ostringstream cnum;
    cnum << thresh;
    std::string name("CM_T_Thresh_" + cnum.str());
    std::string title("CMM-CP Hits Totals Threshold "
                                                  + cnum.str() + ";By Crate");
    TH1F* hist = book1F(name, title, 20, 0, 20);
    int bin = 0;
    for (int crate = 0; crate < s_crates; ++crate) {
      hist->GetXaxis()->SetBinLabel(++bin, "Remote");
      hist->GetXaxis()->SetBinLabel(++bin, "Remote");
      hist->GetXaxis()->SetBinLabel(++bin, "Remote");
      hist->GetXaxis()->SetBinLabel(++bin, "Local");
      hist->GetXaxis()->SetBinLabel(++bin, "Total");
    }
    m_v_CMM_T_thresholds.push_back(hist);
  }

  newGroup(cpErrDir + "_crate_sys", expert, run );

  m_h_LOCeqREM = book1F("CM_LOCeqREM",
      "CMM Local-Remote Totals Non-Zero Match;Crate/Left-Right", 6, 0, 6);
  m_h_LOCneREM = book1F("CM_LOCneREM",
      "CMM Local-Remote Totals Non-Zero Mismatch;Crate/Left-Right", 6, 0, 6);
  m_h_LOCnoREM = book1F("CM_LOCnoREM",
                  "CMM Local/No Remote Totals;Crate/Left-Right", 6, 0, 6);
  m_h_REMnoLOC = book1F("CM_REMnoLOC",
                  "CMM Remote/No Local Totals;Crate/Left-Right", 6, 0, 6);
  setCmmLocRemLabels(m_h_LOCeqREM);
  setCmmLocRemLabels(m_h_LOCneREM);
  setCmmLocRemLabels(m_h_LOCnoREM);
  setCmmLocRemLabels(m_h_REMnoLOC);

  //  CMM error bits

  newGroup(cpErrDir + "_CMM_parity", expert, run );

  m_h_CMM_L_parity = book1F("CMM_L_parity",
                     "CMM Parity Errors Em/Tau (Left);Crate/CPM", 59, 0, 59);
  setThresholdLabels(m_h_CMM_L_parity);
  m_h_CMM_L_parity->GetXaxis()->SetBinLabel(57, "R0");
  m_h_CMM_R_parity = book1F("CMM_R_parity",
                     "CMM Parity Errors Em (Right);Crate/CPM", 59, 0, 59);
  setThresholdLabels(m_h_CMM_R_parity);
  m_h_CMM_R_parity->GetXaxis()->SetBinLabel(57, "R0");

  newGroup(cpErrDir + "_CMM_status", expert, run );

  m_h_CMM_status = book1F("CMM_status", "CMM Sub-status bits", 8, 0., 8.);
  setStatusLabels(m_h_CMM_status);
  m_h_CMM_status_loc = book2F("CMM_status_loc",
                              "CMM Sub-status location;;Crate/Left-Right",
                              8, 0., 8., 8, 0., 8.);
  setStatusLabels(m_h_CMM_status_loc);
  setCmmLocLabels(m_h_CMM_status_loc);

  //  Multiple slice CMM Hit plots

  //....

  //  CPM/CMM Hits event by event comparisons

  newGroup(cpErrDir + "_CPM_CMM", expert, run );

  m_h_CPMeqCMM_hits1 = book1F("CPMeqCMM_hits1",
           "CPM-CMM Hits non-zero match Em/Tau (Left);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPMeqCMM_hits1);
  m_h_CPMeqCMM_hits0 = book1F("CPMeqCMM_hits0",
           "CPM-CMM Hits non-zero match Em (Right);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPMeqCMM_hits0);
  m_h_CPMneCMM_hits1 = book1F("CPMneCMM_hits1",
           "CPM-CMM Hits non-zero mismatch Em/Tau (Left);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPMneCMM_hits1);
  m_h_CPMneCMM_hits0 = book1F("CPMneCMM_hits0",
           "CPM-CMM Hits non-zero mismatch Em (Right);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPMneCMM_hits0);
  m_h_CPMnoCMM_hits1 = book1F("CPMnoCMM_hits1",
                 "CPM /no CMM Hits Em/Tau (Left);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPMnoCMM_hits1);
  m_h_CPMnoCMM_hits0 = book1F("CPMnoCMM_hits0",
                 "CPM /no CMM Hits Em (Right);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CPMnoCMM_hits0);
  m_h_CMMnoCPM_hits1 = book1F("CMMnoCPM_hits1",
                 "CMM /no CPM Hits Em/Tau (Left);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CMMnoCPM_hits1);
  m_h_CMMnoCPM_hits0 = book1F("CMMnoCPM_hits0",
                 "CMM /no CPM Hits Em (Right);Crate/CPM", 56, 0, 56);
  setThresholdLabels(m_h_CMMnoCPM_hits0);

  //  Error Summary

  newGroup(cpErrDir + "_summary", shift, run );

  m_h_CP_errors = book1F("CP_Error_Summary", "CP Error Summary;;Events",
                          NumberOfSummaryBins, 0, NumberOfSummaryBins);
  m_h_CP_errors->GetXaxis()->SetBinLabel(1+PPrCPMTransfer,   "PPr-CPM");
  m_h_CP_errors->GetXaxis()->SetBinLabel(1+CoreOverlap,      "Phi overlap");
  m_h_CP_errors->GetXaxis()->SetBinLabel(1+CPMParity,        "CPM parity");
  m_h_CP_errors->GetXaxis()->SetBinLabel(1+CPMLink,          "CPM link");
  m_h_CP_errors->GetXaxis()->SetBinLabel(1+CPMStatus,        "CPM status");
  m_h_CP_errors->GetXaxis()->SetBinLabel(1+RoIParity,        "RoI parity");
  m_h_CP_errors->GetXaxis()->SetBinLabel(1+CPMCMMTransfer,   "CPM-CMM");
  m_h_CP_errors->GetXaxis()->SetBinLabel(1+CMMParity,        "CMM parity");
  m_h_CP_errors->GetXaxis()->SetBinLabel(1+CMMStatus,        "CMM status");
  m_h_CP_errors->GetXaxis()->SetBinLabel(1+CrateSysTransfer, "crate-sys");
  m_h_CP_errors->GetXaxis()->SetLabelSize(0.06);

  delete m_monGroup;
  m_monGroup = 0;

  } // end if (isNewRun ...

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
    log << MSG::DEBUG<< "No Trigger Tower container found"<< endreq; 
  }

  //Retrieve Core CPM Towers from SG
  const CpmTowerCollection* cpmTowerTES = 0; 
  sc = m_storeGate->retrieve(cpmTowerTES, m_cpmTowerLocation); 
  if( sc.isFailure()  ||  !cpmTowerTES ) {
    log << MSG::DEBUG<< "No Core CPM Tower container found"<< endreq; 
  }

  //Retrieve Overlap CPM Towers from SG
  const CpmTowerCollection* cpmTowerOverlapTES = 0; 
  sc = m_storeGate->retrieve(cpmTowerOverlapTES, m_cpmTowerLocationOverlap); 
  if( sc.isFailure()  ||  !cpmTowerOverlapTES ) {
    log << MSG::DEBUG<< "No Overlap CPM Tower container found"<< endreq; 
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

  // Error summary plot flags
  std::vector<int> errorPlot(NumberOfSummaryBins);

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
        if (m_noiseSignalSplit) m_h_TT_Em_Et_s->Fill(em, 1.);
        m_h_TT_Em_eta->Fill(eta, 1.);
        m_h_TT_Em_phi->Fill(phiMod, 1.);
        m_h_TT_Em_eta_phi->Fill(eta, phiMod, 1.);
        m_h_TT_Em_eta_phi_w->Fill(eta, phiMod, em);
      }
      if (had && eta > -2.5 && eta < 2.5) {
        m_h_TT_Had_Et->Fill(had, 1.);
        if (m_noiseSignalSplit) m_h_TT_Had_Et_s->Fill(had, 1.);
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
        if (m_noiseSignalSplit) m_h_CT_Em_Et_s->Fill(em, 1.);
        m_h_CT_Em_eta->Fill(eta, 1.);
        m_h_CT_Em_phi->Fill(phiMod, 1.);
        m_h_CT_Em_eta_phi->Fill(eta, phiMod, 1.);
        m_h_CT_Em_eta_phi_w->Fill(eta, phiMod, em);
      }
      if (had) {
        m_h_CT_Had_Et->Fill(had, 1.);
        if (m_noiseSignalSplit) m_h_CT_Had_Et_s->Fill(had, 1.);
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

      // Error summary flags
      if (emError.get(LVL1::DataError::Parity) ||
          hadError.get(LVL1::DataError::Parity)) errorPlot[CPMParity] = 1;
      if (emError.get(LVL1::DataError::LinkDown) ||
          hadError.get(LVL1::DataError::LinkDown)) errorPlot[CPMLink] = 1;
      if (status) errorPlot[CPMStatus] = 1;

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
    // Error summary flags
    if (ttEm != cpEm || ttHad != cpHad) errorPlot[PPrCPMTransfer] = 1;
  }

  // Core/Overlap Tower comparison

  CpmTowerMap coMap;
  if (cpmTowerOverlapTES) {
    CpmTowerCollection::const_iterator ctIterator =
                                                  cpmTowerOverlapTES->begin(); 
    CpmTowerCollection::const_iterator ctIteratorEnd =
                                                  cpmTowerOverlapTES->end(); 
    for (; ctIterator != ctIteratorEnd; ++ctIterator) {
      double phi = (*ctIterator)->phi();
      if (phi >= 2.*M_PI) phi -= 2.*M_PI;//bug TrigT1CaloByteStream-00-03-14/15
      const unsigned int key = towerKey.ttKey( phi, (*ctIterator)->eta() );
      if (key > maxKey) maxKey = key;
      coMap.insert(std::make_pair(key, *ctIterator));
    }
  }
  ++maxKey;

  LVL1::CoordToHardware conv;
  cpMapIter    = cpMap.begin();
  cpMapIterEnd = cpMap.end();
  CpmTowerMap::const_iterator coMapIter    = coMap.begin();
  CpmTowerMap::const_iterator coMapIterEnd = coMap.end();

  while (cpMapIter != cpMapIterEnd || coMapIter != coMapIterEnd) {

    unsigned int cpKey = maxKey;
    unsigned int coKey = maxKey;
    int cpEm  = 0;
    int cpHad = 0;
    int coEm  = 0;
    int coHad = 0;
    int cpEmErr  = 0;
    int cpHadErr = 0;
    int coEmErr  = 0;
    int coHadErr = 0;
    double eta = 0.;
    double phi = 0.;

    if (cpMapIter != cpMapIterEnd) cpKey = cpMapIter->first;
    if (coMapIter != coMapIterEnd) coKey = coMapIter->first;

    if ((coMapIter == coMapIterEnd) || (coKey > cpKey)) {

      // Core CPMTower but no Overlap CPMTower

      const LVL1::CPMTower* cp = cpMapIter->second;
      ++cpMapIter;
      eta = cp->eta();
      phi = cp->phi();
      // Does it correspond to overlap channel?
      const int crate = conv.cpCrateOverlap(LVL1::Coordinate(phi, eta));
      if (crate >= s_crates) continue;
      cpEm  = cp->emEnergy();
      cpHad = cp->hadEnergy();
      cpEmErr  = cp->emError();
      cpHadErr = cp->hadError();

    } else if ((cpMapIter == cpMapIterEnd) || (cpKey > coKey)) {

      // Overlap CPMTower but no Core CPMTower

      const LVL1::CPMTower* co = coMapIter->second;
      ++coMapIter;
      eta = co->eta();
      phi = co->phi();
      if (phi >= 2.*M_PI) phi -= 2.*M_PI;//bug TrigT1CaloByteStream-00-03-14/15
      const int crate = conv.cpCrateOverlap(LVL1::Coordinate(phi, eta));
      if (crate >= s_crates) {
        log << MSG::DEBUG << "Overlap mapping error, eta/phi: "
                          << eta << "/" << phi << endreq;
        continue;
      }
      coEm  = co->emEnergy();
      coHad = co->hadEnergy();
      coEmErr  = co->emError();
      coHadErr = co->hadError();

    } else {

      // Have both

      const LVL1::CPMTower* cp = cpMapIter->second;
      const LVL1::CPMTower* co = coMapIter->second;
      ++cpMapIter;
      ++coMapIter;
      eta = cp->eta();
      phi = cp->phi();
      const int crate = conv.cpCrateOverlap(LVL1::Coordinate(phi, eta));
      if (crate >= s_crates) {
        log << MSG::DEBUG << "Overlap mapping error, eta/phi: "
                          << eta << "/" << phi << endreq;
        continue;
      }
      cpEm  = cp->emEnergy();
      cpHad = cp->hadEnergy();
      coEm  = co->emEnergy();
      coHad = co->hadEnergy();
      cpEmErr  = cp->emError();
      cpHadErr = cp->hadError();
      coEmErr  = co->emError();
      coHadErr = co->hadError();

    }

    const double phiMod = phi * m_phiScale;
    // Ignore data with known errors
    //if ( !cpEmErr && !coEmErr ) {
      m_h_CTeqCO_Em_eta_phi->Fill(eta, phiMod, cpEm && coEm && cpEm == coEm);
      m_h_CTneCO_Em_eta_phi->Fill(eta, phiMod, cpEm && coEm && cpEm != coEm);
      m_h_CTnoCO_Em_eta_phi->Fill(eta, phiMod, cpEm  && !coEm);
      m_h_COnoCT_Em_eta_phi->Fill(eta, phiMod, !cpEm && coEm);
    //}
    //if ( !cpHadErr && !coHadErr ) {
      m_h_CTeqCO_Had_eta_phi->Fill(eta, phiMod,
                                        cpHad && coHad && cpHad == coHad);
      m_h_CTneCO_Had_eta_phi->Fill(eta, phiMod,
                                        cpHad && coHad && cpHad != coHad);
      m_h_CTnoCO_Had_eta_phi->Fill(eta, phiMod, cpHad && !coHad);
      m_h_COnoCT_Had_eta_phi->Fill(eta, phiMod, !cpHad && coHad);
    //}
    //if ((!cpEmErr && !coEmErr && cpEm && coEm && cpEm != coEm) ||
    //    (!cpHadErr && !coHadErr && cpHad && coHad && cpHad != coHad)) {
    if ((cpEm && coEm && cpEm != coEm) || (cpHad && coHad && cpHad != coHad)) {
      log << MSG::DEBUG
          << "Core/Overlap CPM Tower mismatch, eta/phi/cpEm/cpHad/coEm/coHad: "
	  << eta << "/" << phi << "/" << cpEm << "/" << cpHad << "/"
	  << coEm << "/" << coHad << endreq;
    }
    // Error summary flags
    if (cpEm != coEm || cpHad != coHad) errorPlot[CoreOverlap] = 1;
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
      const int bin1  = crate * s_modules + cpm - 1;
      //const int bin2  = chip * 8 + loc;
      std::vector<TH1F*>::const_iterator hist1 = m_v_RoI_thresholds.begin();
      std::vector<TH2F*>::const_iterator hist2 = m_v_RoI_2D_thresholds.begin();
      for (int thresh = 0; thresh < s_thresholds; ++thresh) {
        const int hit = (hits >> thresh) & 0x1;
        if (hit) {
          (*hist1)->Fill(bin1, 1.);
	  (*hist2)->Fill(eta, phiMod, 1.);
        }
        ++hist1;
        ++hist2;
      }
      if ((*crIterator)->error()) {
        m_h_RoI_Parity->Fill(eta, phiMod, 1.);
	errorPlot[RoIParity] = 1;
      }
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
      const int bin   = crate * s_modules + cpm - 1;
      const int peak   = (*chIterator)->peak();
      const int slices = ((*chIterator)->HitsVec0()).size();
      m_v_CPM_slices[crate]->Fill(slices, peak, 1.);
      std::vector<TH1F*>::const_iterator hist = m_v_thresholds.begin();
      for (int thresh = 0; thresh < s_thresholds/2; ++thresh) {
        const int hit = (hits0 >> thresh*s_threshBits) & s_threshMask;
        if (hit) (*hist)->Fill(bin, hit);
        ++hist;
      }
      for (int thresh = 0; thresh < s_thresholds/2; ++thresh) {
        const int hit = (hits1 >> thresh*s_threshBits) & s_threshMask;
        if (hit) (*hist)->Fill(bin, hit);
        ++hist;
      }
      const unsigned int key = crate * s_modules + cpm;
      if (key > maxKey) maxKey = key;
      cpmMap.insert(std::make_pair(key, *chIterator));
    }
  }

  //=============================================
  //  CMM-CP Hits
  //=============================================

  CmmCpHitsMap cmmMap;

  std::vector<unsigned int> hits0Local(s_crates-1);
  std::vector<unsigned int> hits1Local(s_crates-1);
  std::vector<unsigned int> hits0Remote(s_crates-1);
  std::vector<unsigned int> hits1Remote(s_crates-1);
  const int systemCrate = s_crates - 1;

  if (cmmCpHitsTES) {
    CmmCpHitsCollection::const_iterator cmIterator    = cmmCpHitsTES->begin(); 
    CmmCpHitsCollection::const_iterator cmIteratorEnd = cmmCpHitsTES->end(); 
    for (; cmIterator != cmIteratorEnd; ++cmIterator) {
      const unsigned int hits0 = (*cmIterator)->HitWord0();
      const unsigned int hits1 = (*cmIterator)->HitWord1();
      const int crate  = (*cmIterator)->crate();
      const int dataId = (*cmIterator)->dataID();
      const int bin    = (dataId <= s_modules) ? crate*s_modules + dataId - 1
  		                            : crate*5 + dataId - s_modules - 1; 
      const int peak   = (*cmIterator)->peak();
      const int slices = ((*cmIterator)->HitsVec0()).size();
      m_v_CMM_slices[crate]->Fill(slices, peak, 1.);
      std::vector<TH1F*>::const_iterator hist1 = m_v_CMM_thresholds.begin();
      std::vector<TH1F*>::const_iterator hist2 = m_v_CMM_T_thresholds.begin();
      for (int thresh = 0; thresh < s_thresholds/2; ++thresh) {
        const int hit = (hits0 >> thresh*s_threshBits) & s_threshMask;
        if (dataId <= s_modules) {
          if (hit) (*hist1)->Fill(bin, hit);
          ++hist1;
        } else {
          if (hit) (*hist2)->Fill(bin, hit);
	  ++hist2;
        }
      }
      for (int thresh = 0; thresh < s_thresholds/2; ++thresh) {
        const int hit = (hits1 >> thresh*s_threshBits) & s_threshMask;
        if (dataId <= s_modules) {
          if (hit) (*hist1)->Fill(bin, hit);
          ++hist1;
        } else {
          if (hit) (*hist2)->Fill(bin, hit);
	  ++hist2;
        }
      }
      // Save hits for cable check
      if (crate != systemCrate) {
        if (dataId == LVL1::CMMCPHits::LOCAL) {
          hits0Local[crate] = hits0;;
          hits1Local[crate] = hits1;
        }
      } else {
        if (dataId == LVL1::CMMCPHits::REMOTE_0) {
	  hits0Remote[0] = hits0;
	  hits1Remote[0] = hits1;
        } else if (dataId == LVL1::CMMCPHits::REMOTE_1) {
	  hits0Remote[1] = hits0;
	  hits1Remote[1] = hits1;
        } else if (dataId == LVL1::CMMCPHits::REMOTE_2) {
	  hits0Remote[2] = hits0;
	  hits1Remote[2] = hits1;
        }
      }
      // Errors
      const LVL1::DataError hit0Err((*cmIterator)->Error0());
      const LVL1::DataError hit1Err((*cmIterator)->Error1());
      if (dataId <= s_modules) {
        m_h_CMM_R_parity->Fill(bin, hit0Err.get(LVL1::DataError::Parity));
        m_h_CMM_L_parity->Fill(bin, hit1Err.get(LVL1::DataError::Parity));
      } else {
	int remBin = 0;
        if (dataId == LVL1::CMMCPHits::REMOTE_0) {
	  remBin = s_crates*s_modules;
        } else if (dataId == LVL1::CMMCPHits::REMOTE_1) {
	  remBin = s_crates*s_modules + 1;
        } else if (dataId == LVL1::CMMCPHits::REMOTE_2) {
	  remBin = s_crates*s_modules + 2;
        }
	if (remBin) {
	  m_h_CMM_R_parity->Fill(remBin, hit0Err.get(LVL1::DataError::Parity));
	  m_h_CMM_L_parity->Fill(remBin, hit1Err.get(LVL1::DataError::Parity));
        }
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
      // Error summary flags
      if (hit0Err.get(LVL1::DataError::Parity) ||
          hit1Err.get(LVL1::DataError::Parity)) errorPlot[CMMParity] = 1;
      if (status0 || status1) errorPlot[CMMStatus] = 1;

      if (dataId <= s_modules) {
        const unsigned int key = crate * s_modules + dataId;
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
    const int bin = crate * s_modules + cpm - 1;
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
      const int nthresh = s_thresholds/2;
      for (int thresh = 0; thresh < nthresh; ++thresh) {
        const int hit = (cpmHits0 >> thresh*s_threshBits) & s_threshMask;
	log << MSG::DEBUG << hit;
	if (thresh < nthresh - 1) log << MSG::DEBUG << ",";
	else log << MSG::DEBUG << ";";
      }
      for (int thresh = 0; thresh < nthresh; ++thresh) {
        const int hit = (cpmHits1 >> thresh*s_threshBits) & s_threshMask;
	log << MSG::DEBUG << hit;
	if (thresh < nthresh - 1) log << MSG::DEBUG << ",";
	else log << MSG::DEBUG << "/";
      }
      for (int thresh = 0; thresh < nthresh; ++thresh) {
        const int hit = (cmmHits0 >> thresh*s_threshBits) & s_threshMask;
	log << MSG::DEBUG << hit;
	if (thresh < nthresh - 1) log << MSG::DEBUG << ",";
	else log << MSG::DEBUG << ";";
      }
      for (int thresh = 0; thresh < nthresh; ++thresh) {
        const int hit = (cmmHits1 >> thresh*s_threshBits) & s_threshMask;
	log << MSG::DEBUG << hit;
	if (thresh < nthresh - 1) log << MSG::DEBUG << ",";
	else log << MSG::DEBUG << endreq;
      }
    }
    // Error summary flags
    if (cpmHits0 != cmmHits0 || cpmHits1 != cmmHits1) {
      errorPlot[CPMCMMTransfer] = 1;
    }
  }

  // Hits totals local crate -> system crate transfer check

  for (int crate = 0; crate < systemCrate; ++crate) {
    const unsigned int h0Loc = hits0Local[crate];
    const unsigned int h1Loc = hits1Local[crate];
    const unsigned int h0Rem = hits0Remote[crate];
    const unsigned int h1Rem = hits1Remote[crate];
    if (h0Loc && h0Rem && h0Loc == h0Rem) {
      m_h_LOCeqREM->Fill(2*crate + 1, 1.); //hits0 == right
    } else if (h0Loc && h0Rem && h0Loc != h0Rem) {
      m_h_LOCneREM->Fill(2*crate + 1, 1.);
    } else if (h0Loc && !h0Rem) {
      m_h_LOCnoREM->Fill(2*crate + 1, 1.);
    } else if (!h0Loc && h0Rem) {
      m_h_REMnoLOC->Fill(2*crate + 1, 1.);
    }
    if (h1Loc && h1Rem && h1Loc == h1Rem) {
      m_h_LOCeqREM->Fill(2*crate, 1.);     //hits1 == left
    } else if (h1Loc && h1Rem && h1Loc != h1Rem) {
      m_h_LOCneREM->Fill(2*crate, 1.);
    } else if (h1Loc && !h1Rem) {
      m_h_LOCnoREM->Fill(2*crate, 1.);
    } else if (!h1Loc && h1Rem) {
      m_h_REMnoLOC->Fill(2*crate, 1.);
    }
    // Error summary flags
    if (h0Loc != h0Rem || h1Loc != h1Rem) errorPlot[CrateSysTransfer] = 1;
  }

  // Update error summary plot

  for (int err = 0; err < NumberOfSummaryBins; ++err) {
    m_h_CP_errors->Fill(err, errorPlot[err]);
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
  }

  return StatusCode::SUCCESS;
}

TH1F* TrigT1CaloCpmMonTool::book1F(const std::string& name,
                                   const std::string& title,
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

TH2F* TrigT1CaloCpmMonTool::book2F(const std::string& name,
                                   const std::string& title,
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

void TrigT1CaloCpmMonTool::setCmmLocRemLabels(TH1* hist)
{
  hist->GetXaxis()->SetBinLabel(1, "0/L");
  hist->GetXaxis()->SetBinLabel(2, "0/R");
  hist->GetXaxis()->SetBinLabel(3, "1/L");
  hist->GetXaxis()->SetBinLabel(4, "1/R");
  hist->GetXaxis()->SetBinLabel(5, "2/L");
  hist->GetXaxis()->SetBinLabel(6, "2/R");
  hist->GetXaxis()->SetLabelSize(0.06);
}
