// ********************************************************************
//
// NAME:     TrigT1CaloCpmMonTool.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Peter Faulkner
//	     
//
// ********************************************************************
#ifndef TRIGT1CALOCPMMONTOOL_H
#define TRIGT1CALOCPMMONTOOL_H

#include <map>
#include <string>
#include <vector>

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ServiceHandle.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "DataModel/DataVector.h"

class TH1F;
class TH2F;
class StoreGateSvc;

namespace LVL1 {
  class CPMTower;
  class CPMHits;
  class CMMCPHits;
  class CPMRoI;
  class TriggerTower;
}

class TrigT1CaloCpmMonTool: public ManagedMonitorToolBase
{

public:
  
  TrigT1CaloCpmMonTool(const std::string & type, const std::string & name,
		       const IInterface* parent);
    

  virtual ~TrigT1CaloCpmMonTool();

  virtual StatusCode initialize();
    
  virtual StatusCode bookHistograms(bool isNewEventsBlock, bool isNewLumiBlock,
                                                           bool isNewRun);
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms(bool isEndOfEventsBlock,
                                    bool isEndOfLumiBlock, bool isEndOfRun);

private:

  enum SummaryErrors { CPMParity, CPMLink, CPMStatus, RoIParity,
                       CMMParity, CMMStatus, NumberOfSummaryBins };

  typedef DataVector<LVL1::CPMTower>     CpmTowerCollection;
  typedef DataVector<LVL1::CPMHits>      CpmHitsCollection;
  typedef DataVector<LVL1::CMMCPHits>    CmmCpHitsCollection;
  typedef DataVector<LVL1::CPMRoI>       CpmRoiCollection;
  typedef DataVector<LVL1::TriggerTower> TriggerTowerCollection;

  typedef std::map<unsigned int, LVL1::TriggerTower*> TriggerTowerMap;
  typedef std::map<unsigned int, LVL1::CPMTower*>     CpmTowerMap;
  typedef std::map<unsigned int, LVL1::CPMHits*>      CpmHitsMap;
  typedef std::map<unsigned int, LVL1::CMMCPHits*>    CmmCpHitsMap;
  
  static const int s_crates     = 4;
  static const int s_modules    = 14; // Modules numbered 1-14
  static const int s_maxSlices  = 5;
  static const int s_thresholds = 16;
  static const int s_threshBits = 3;
  static const int s_threshMask = 0x7;
  
  TH1F* book1F(const std::string& name, const std::string& title,
                                    int nx, double xmin, double xmax);
  TH2F* book2F(const std::string& name, const std::string& title,
                                    int nx, double xmin, double xmax,
                                    int ny, double ymin, double ymax);
  void  setStatusLabels(TH1* hist);
  void  setLabelsCNSTS(TH2* hist);
  void  setLabelsPSCS(TH2* hist);
  void  setLabelsCMT(TH2* hist);
  void  setLabelsT(TH2* hist);
  void  setLabelsCPM(TH2* hist);
  void  setYLabelsCPM(TH2* hist);
  void  setLabelsMCLR(TH2* hist);
  void  setLabelsCLR(TH2* hist);
  void  setLabelsST(TH2* hist);

  ServiceHandle<StoreGateSvc> m_storeGate;
  mutable MsgStream m_log;

  MonGroup* m_monGroup;

  /// Core CPM tower container StoreGate key
  std::string m_cpmTowerLocation;
  /// Overlap CPM tower container StoreGate key
  std::string m_cpmTowerLocationOverlap;
  /// CPM hits container StoreGate key
  std::string m_cpmHitsLocation;
  /// CMM-CP hits container StoreGate key
  std::string m_cmmCpHitsLocation;
  /// DAQ CPM RoI container StoreGate key
  std::string m_cpmRoiLocation;
  /// RoIB CPM RoI container StoreGate key
  std::string m_cpmRoiLocationRoib;
  /// Trigger Tower container StoreGate key
  std::string m_triggerTowerLocation;
  
  /// Root directory
  std::string m_rootDir;

  /// Phi Units for eta/phi plots
  std::string m_phiUnits;
  /// Phi maximum in wanted units
  double m_phiMax;
  /// Phi scale to convert from radians to wanted units
  double m_phiScale;
  /// Noise/signal energy split
  int m_noiseSignalSplit;
  /// Maximum energy plotted
  int m_maxEnergyRange;
  /// Number of events
  int m_events;
  /// Not used
  bool m_Offline;

  //=======================
  //   Timeslice plots
  //=======================

  TH2F* m_h_CPM_slices;
  TH2F* m_h_CMM_slices;
  std::vector<TH2F*> m_v_PP_CP_slice;
  std::vector<TH2F*> m_v_CP_CM_slice;

  //=============================================
  //   CPM Tower - Trigger Tower comparison plots
  //=============================================

  // TriggerTower plots
  TH1F* m_h_TT_Em_Et;  
  TH1F* m_h_TT_Had_Et;
  TH1F* m_h_TT_Em_Et_s;  
  TH1F* m_h_TT_Had_Et_s;
  TH1F* m_h_TT_Em_eta;
  TH1F* m_h_TT_Had_eta;
  TH1F* m_h_TT_Em_phi;
  TH1F* m_h_TT_Had_phi;
  TH2F* m_h_TT_Em_eta_phi;
  TH2F* m_h_TT_Had_eta_phi;
  TH2F* m_h_TT_Em_eta_phi_w;
  TH2F* m_h_TT_Had_eta_phi_w;
  // CPMTower plots
  TH1F* m_h_CT_Em_Et;  
  TH1F* m_h_CT_Had_Et;
  TH1F* m_h_CT_Em_Et_s;
  TH1F* m_h_CT_Had_Et_s;
  TH1F* m_h_CT_Em_eta;
  TH1F* m_h_CT_Had_eta;
  TH1F* m_h_CT_Em_phi;
  TH1F* m_h_CT_Had_phi;
  TH2F* m_h_CT_Em_eta_phi;
  TH2F* m_h_CT_Had_eta_phi;
  TH2F* m_h_CT_Em_eta_phi_w;
  TH2F* m_h_CT_Had_eta_phi_w;
  // Errors
  TH2F* m_h_CT_Em_parity;
  TH2F* m_h_CT_Had_parity;
  TH2F* m_h_CT_Em_link;
  TH2F* m_h_CT_Had_link;
  TH2F* m_h_CT_status;

  //=============================================
  //  CPM RoIs
  //=============================================

  TH2F* m_h_RoI_thresholds;
  TH2F* m_h_RoI_eta_phi;
  TH2F* m_h_RoI_Em_eta_phi;
  TH2F* m_h_RoI_Tau_eta_phi;
  // Tower saturation
  TH2F* m_h_RoI_Saturation;
  // Parity errors
  TH2F* m_h_RoI_Parity;

  //=============================================
  //  CPM Hits
  //=============================================

  TH2F* m_h_CPM_thresholds;

  //=============================================
  //  CMM-CP Hits
  //=============================================

  TH2F* m_h_CMM_thresholds;
  TH2F* m_h_CMM_T_thresholds;
  // Errors
  TH2F* m_h_CMM_parity;
  TH2F* m_h_CMM_status;

  //=============================================
  //  Error summary
  //=============================================

  TH1F* m_h_CP_errors;
  TH2F* m_h_CP_overview;

};

#endif
