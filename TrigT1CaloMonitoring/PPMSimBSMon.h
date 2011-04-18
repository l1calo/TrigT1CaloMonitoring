// ********************************************************************
//
// NAME:     PPMSimBSMon.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Peter Faulkner
//           Sky French
//	     
//
// ********************************************************************
#ifndef PPMSIMBSMON_H
#define PPMSIMBSMON_H

#include <string>
#include <vector>

#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "DataModel/DataVector.h"

class TH2F_LW;
class TH2I_LW;
class TProfile2D_LW;

class StatusCode;
class L1CaloCondSvc;
class L1CaloPprConditionsContainer;

class TrigT1CaloLWHistogramTool;

namespace LVL1 {
  class TriggerTower;
  class IL1TriggerTowerTool;
}

class PPMSimBSMon: public ManagedMonitorToolBase
{

public:
  
  PPMSimBSMon(const std::string & type, const std::string & name,
		       const IInterface* parent);
    

  virtual ~PPMSimBSMon();

  virtual StatusCode initialize();
  virtual StatusCode finalize();  
  virtual StatusCode bookHistograms(bool isNewEventsBlock, bool isNewLumiBlock,
                                                           bool isNewRun);
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms(bool isEndOfEventsBlock,
                                    bool isEndOfLumiBlock, bool isEndOfRun);

private:

  typedef DataVector<LVL1::TriggerTower> TriggerTowerCollection;
  
  typedef std::vector<int> ErrorVector;

  void  fillEventSample(int crate, int module);

  void simulateAndCompare(const TriggerTowerCollection* ttIn);

  ServiceHandle<L1CaloCondSvc> m_l1CondSvc;

  ToolHandle<LVL1::IL1TriggerTowerTool> m_ttTool;
  ToolHandle<TrigT1CaloLWHistogramTool> m_histTool;
      
  L1CaloPprConditionsContainer* m_conditionsContainer;

  bool m_debug;
  bool m_onlineTest;

  std::string m_rootDir;

  /// Trigger Tower container StoreGate key
  std::string m_triggerTowerLocation;

  /// Number of events
  int m_events;
  /// Number of events over which to sample pedestal
  int m_instantaneous;
  /// Cut on ADC digits for re-simulation
  int m_simulationADCCut;
  /// Histograms booked flag
  bool m_histBooked;

  //=======================
  //   Match/Mismatch plots
  //=======================

  // LUT
  TH2F_LW* m_h_ppm_em_2d_etaPhi_tt_lut_SimEqData;
  TH2F_LW* m_h_ppm_em_2d_etaPhi_tt_lut_SimNeData;
  TH2F_LW* m_h_ppm_em_2d_etaPhi_tt_lut_SimNoData;
  TH2F_LW* m_h_ppm_em_2d_etaPhi_tt_lut_DataNoSim;
  TH2F_LW* m_h_ppm_had_2d_etaPhi_tt_lut_SimEqData;
  TH2F_LW* m_h_ppm_had_2d_etaPhi_tt_lut_SimNeData;
  TH2F_LW* m_h_ppm_had_2d_etaPhi_tt_lut_SimNoData;
  TH2F_LW* m_h_ppm_had_2d_etaPhi_tt_lut_DataNoSim;
  
  //Overal Pedestal
  TProfile2D_LW* m_h_ppm_em_2d_etaPhi_tt_ped_runavg;
  TProfile2D_LW* m_h_ppm_had_2d_etaPhi_tt_ped_runavg;
  TH2F_LW* m_h_ppm_em_2d_etaPhi_tt_ped_worstavg;
  TH2F_LW* m_h_ppm_had_2d_etaPhi_tt_ped_worstavg;
  TH2F_LW* m_h_ppm_em_2d_etaPhi_tt_ped_runrms;
  TH2F_LW* m_h_ppm_had_2d_etaPhi_tt_ped_runrms;
  TProfile2D_LW* m_h_ppm_em_2d_etaPhi_tt_ped_instavg;
  TProfile2D_LW* m_h_ppm_had_2d_etaPhi_tt_ped_instavg;
  TH2F_LW* m_h_ppm_em_2d_etaPhi_tt_ped_instrms;
  TH2F_LW* m_h_ppm_had_2d_etaPhi_tt_ped_instrms;
  TProfile2D_LW* m_h_ppm_em_2d_etaPhi_tt_ped_instavg_B;
  TProfile2D_LW* m_h_ppm_had_2d_etaPhi_tt_ped_instavg_B;
  TH2F_LW* m_h_ppm_em_2d_etaPhi_tt_ped_instrms_B;
  TH2F_LW* m_h_ppm_had_2d_etaPhi_tt_ped_instrms_B;
  
  // Mismatch Histograms
  TH2I_LW* m_h_ppm_2d_LUT_MismatchEvents_cr0cr1;
  TH2I_LW* m_h_ppm_2d_LUT_MismatchEvents_cr2cr3;
  TH2I_LW* m_h_ppm_2d_LUT_MismatchEvents_cr4cr5;
  TH2I_LW* m_h_ppm_2d_LUT_MismatchEvents_cr6cr7;
  
};

#endif
