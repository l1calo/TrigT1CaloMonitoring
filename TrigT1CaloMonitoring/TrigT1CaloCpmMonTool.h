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
  typedef std::map<int, LVL1::CPMHits*>               CpmHitsMap;
  typedef std::map<int, LVL1::CMMCPHits*>             CmmCpHitsMap;
  
  TH1D* book1D(std::string nam, std::string tit,
                                    int nx, double xmin, double xmax);
  TH2D* book2D(std::string nam, std::string tit,
                                    int nx, double xmin, double xmax,
                                    int ny, double ymin, double ymax);
  void  setThresholdLabels(TH1* hist);
  void  setStatusLabels(TH1* hist);

  StoreGateSvc* m_StoreGate;
  MonGroup* m_monGroup;
  std::string m_prefix;
  //std::string m_dir;

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

  //=======================
  //   Timeslice plots
  //=======================

  TH2D* m_h_CPM_slices;
  TH2D* m_h_CMM_slices;

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
  // Mismatch plots
  TH2D* m_h_TTCT_Em_eta_phi;
  TH2D* m_h_TTCT_Had_eta_phi;
  // Errors
  TH2D* m_h_CT_Em_parity;
  TH2D* m_h_CT_Had_parity;
  TH2D* m_h_CT_Em_link;
  TH2D* m_h_CT_Had_link;
  TH1D* m_h_CT_status;
  TH2D* m_h_CT_status_eta_phi;

  //=============================================
  //  CPM Hits
  //=============================================

  std::vector<TH1D*> m_v_thresholds;

  //=============================================
  //  CMM-CP Hits
  //=============================================

  std::vector<TH1D*> m_v_CMM_thresholds;
  std::vector<TH1D*> m_v_CMM_T_thresholds;
  // CPM-CMM mismatch
  TH1D* m_h_CPM_CMM_hits0;
  TH1D* m_h_CPM_CMM_hits1;
  // Errors
  TH1D* m_h_CMM_R_parity;
  TH1D* m_h_CMM_L_parity;
  TH1D* m_h_CMM_status;

  //=============================================
  //  CPM RoIs
  //=============================================

  std::vector<TH1D*> m_v_RoI_thresholds;
  std::vector<TH2D*> m_v_RoI_2D_thresholds;

};

#endif
