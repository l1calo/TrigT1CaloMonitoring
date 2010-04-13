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
#include "GaudiKernel/ToolHandle.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "DataModel/DataVector.h"

class TH1F;
class TH2F;
class StoreGateSvc;
class TrigT1CaloMonErrorTool;
class TrigT1CaloHistogramTool;

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
  virtual StatusCode finalize();
    
  virtual StatusCode bookHistograms(bool isNewEventsBlock, bool isNewLumiBlock,
                                                           bool isNewRun);
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms(bool isEndOfEventsBlock,
                                    bool isEndOfLumiBlock, bool isEndOfRun);

private:

  enum SummaryErrors { EMParity, EMLink, HadParity, HadLink, CPMStatus,
                       RoIParity, CMMParity, CMMStatus, NumberOfSummaryBins };

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

  ServiceHandle<StoreGateSvc> m_storeGate;
  ToolHandle<TrigT1CaloMonErrorTool> m_errorTool;
  ToolHandle<TrigT1CaloHistogramTool> m_histTool;
  mutable MsgStream m_log;

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
  /// Trigger Tower container StoreGate key
  std::string m_triggerTowerLocation;
  
  /// Root directory
  std::string m_rootDir;

  /// Maximum energy plotted
  int m_maxEnergyRange;
  /// Number of events
  int m_events;

  //=======================
  //   Timeslice plots
  //=======================

  TH2F* m_h_CPM_slices;
  TH2F* m_h_CMM_slices;
  TH2F* m_h_PP_CP_slice;
  TH2F* m_h_CP_CM_slice;

  //=============================================
  //   CPM Tower - Trigger Tower plots
  //=============================================

  // TriggerTower plots
  TH2F* m_h_TT_Em_eta_phi;
  TH2F* m_h_TT_Had_eta_phi;
  // CPMTower plots
  TH1F* m_h_CT_Em_Et;  
  TH1F* m_h_CT_Had_Et;
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
