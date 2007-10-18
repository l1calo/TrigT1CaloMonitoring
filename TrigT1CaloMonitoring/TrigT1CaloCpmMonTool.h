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

#include "GaudiKernel/ServiceHandle.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "DataModel/DataVector.h"

class TH1D;
class TH2D;
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

  typedef DataVector<LVL1::CPMTower>     CpmTowerCollection;
  typedef DataVector<LVL1::CPMHits>      CpmHitsCollection;
  typedef DataVector<LVL1::CMMCPHits>    CmmCpHitsCollection;
  typedef DataVector<LVL1::CPMRoI>       CpmRoiCollection;
  typedef DataVector<LVL1::TriggerTower> TriggerTowerCollection;

  typedef std::map<unsigned int, LVL1::TriggerTower*> TriggerTowerMap;
  typedef std::map<unsigned int, LVL1::CPMTower*>     CpmTowerMap;
  typedef std::map<unsigned int, LVL1::CPMHits*>      CpmHitsMap;
  typedef std::map<unsigned int, LVL1::CMMCPHits*>    CmmCpHitsMap;
  
  
  TH1D* book1D(const std::string& name, const std::string& title,
                                    int nx, double xmin, double xmax);
  TH2D* book2D(const std::string& name, const std::string& title,
                                    int nx, double xmin, double xmax,
                                    int ny, double ymin, double ymax);
  void  newGroup(const std::string& system, LevelOfDetail_t level,
                                            Interval_t interval);
  void  setThresholdLabels(TH1* hist);
  void  setStatusLabels(TH1* hist);
  void  setCmmLocLabels(TH2* hist);

  ServiceHandle<StoreGateSvc> m_storeGate;

  MonGroup* m_monGroup;

  /// CPM tower container StoreGate key
  std::string m_cpmTowerLocation;
  /// CPM hits container StoreGate key
  std::string m_cpmHitsLocation;
  /// CMM-CP hits container StoreGate key
  std::string m_cmmCpHitsLocation;
  /// CPM RoI container StoreGate key
  std::string m_cpmRoiLocation;
  /// Trigger Tower container StoreGate key
  std::string m_triggerTowerLocation;
  
  /// Root directory
  std::string m_rootDir;
  /// Flag to put all plots in one dir
  bool m_oneDir;

  /// Phi Units for eta/phi plots
  std::string m_phiUnits;
  /// Phi maximum in wanted units
  double m_phiMax;
  /// Phi scale to convert from radians to wanted units
  double m_phiScale;

  //=======================
  //   Timeslice plots
  //=======================

  std::vector<TH2D*> m_v_CPM_slices;
  std::vector<TH2D*> m_v_PP_CP_slice;
  std::vector<TH2D*> m_v_CMM_slices;
  std::vector<TH2D*> m_v_CP_CM_slice;

  //=============================================
  //   CPM Tower - Trigger Tower comparison plots
  //=============================================

  // TriggerTower plots
  TH1D* m_h_TT_Em_Et;  
  TH1D* m_h_TT_Had_Et;
  TH1D* m_h_TT_Em_Et_s;  
  TH1D* m_h_TT_Had_Et_s;
  TH1D* m_h_TT_Em_eta;
  TH1D* m_h_TT_Had_eta;
  TH1D* m_h_TT_Em_phi;
  TH1D* m_h_TT_Had_phi;
  TH2D* m_h_TT_Em_eta_phi;
  TH2D* m_h_TT_Had_eta_phi;
  TH2D* m_h_TT_Em_eta_phi_w;
  TH2D* m_h_TT_Had_eta_phi_w;
  // CPMTower plots
  TH1D* m_h_CT_Em_Et;  
  TH1D* m_h_CT_Had_Et;
  TH1D* m_h_CT_Em_Et_s;
  TH1D* m_h_CT_Had_Et_s;
  TH1D* m_h_CT_Em_eta;
  TH1D* m_h_CT_Had_eta;
  TH1D* m_h_CT_Em_phi;
  TH1D* m_h_CT_Had_phi;
  TH2D* m_h_CT_Em_eta_phi;
  TH2D* m_h_CT_Had_eta_phi;
  TH2D* m_h_CT_Em_eta_phi_w;
  TH2D* m_h_CT_Had_eta_phi_w;
  // Errors
  TH2D* m_h_CT_Em_parity;
  TH2D* m_h_CT_Had_parity;
  TH2D* m_h_CT_Em_link;
  TH2D* m_h_CT_Had_link;
  TH1D* m_h_CT_status;
  TH2D* m_h_CT_status_eta_phi;
  // Mismatch plots
  TH2D* m_h_TTeqCT_Em_eta_phi;
  TH2D* m_h_TTneCT_Em_eta_phi;
  TH2D* m_h_TTnoCT_Em_eta_phi;
  TH2D* m_h_CTnoTT_Em_eta_phi;
  TH2D* m_h_TTeqCT_Had_eta_phi;
  TH2D* m_h_TTneCT_Had_eta_phi;
  TH2D* m_h_TTnoCT_Had_eta_phi;
  TH2D* m_h_CTnoTT_Had_eta_phi;

  //=============================================
  //  CPM RoIs
  //=============================================

  std::vector<TH1D*> m_v_RoI_thresholds;
  std::vector<TH2D*> m_v_RoI_2D_thresholds;
  // Parity errors
  TH2D* m_h_RoI_Parity;

  //=============================================
  //  CPM Hits
  //=============================================

  std::vector<TH1D*> m_v_thresholds;

  //=============================================
  //  CMM-CP Hits
  //=============================================

  std::vector<TH1D*> m_v_CMM_thresholds;
  std::vector<TH1D*> m_v_CMM_T_thresholds;
  // Errors
  TH1D* m_h_CMM_R_parity;
  TH1D* m_h_CMM_L_parity;
  TH1D* m_h_CMM_status;
  TH2D* m_h_CMM_status_loc;
  // CPM-CMM mismatch
  TH1D* m_h_CPMeqCMM_hits0;
  TH1D* m_h_CPMeqCMM_hits1;
  TH1D* m_h_CPMneCMM_hits0;
  TH1D* m_h_CPMneCMM_hits1;
  TH1D* m_h_CPMnoCMM_hits0;
  TH1D* m_h_CPMnoCMM_hits1;
  TH1D* m_h_CMMnoCPM_hits0;
  TH1D* m_h_CMMnoCPM_hits1;

  //=============================================
  //  Error summary
  //=============================================

  TH1D* m_h_CP_errors;

};

#endif
