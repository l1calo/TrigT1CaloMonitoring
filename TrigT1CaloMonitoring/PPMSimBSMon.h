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

#include <map>
#include <string>
#include <vector>

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

#include "CaloIdentifier/CaloIdManager.h"
#include "CaloIdentifier/CaloLVL1_ID.h"
#include "CaloTriggerTool/CaloTriggerTowerService.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "DataModel/DataVector.h"
#include "TrigT1CaloCalibTools/L1CaloTTIdTools.h"

#include "Identifier/Identifier.h"

class TH1;
class TH2;
class TH1F;
class TH2F;
class TH2I;
class StoreGateSvc;

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

  typedef std::map<int, LVL1::TriggerTower*> TriggerTowerMap;
  
  TH1F* book1F(std::string nam, std::string tit,
                                    int nx, double xmin, double xmax);
  TH2F* book2F(std::string nam, std::string tit,
                                    int nx, const double* xbins,
                                    int ny, double ymin, double ymax);
  TH2F* bookEtaPhi(std::string nam, std::string tit);
  TH2I* book2I(std::string nam, std::string tit,
                                    int nx, double xmin, double xmax,
                                    int ny, double ymin, double ymax);
  void  fillEventSample(int crate, int module);

  void simulateAndCompare(const TriggerTowerCollection* ttIn);

  ServiceHandle<StoreGateSvc> m_storeGate;

  ToolHandle<LVL1::IL1TriggerTowerTool> m_ttTool;
      
  mutable MsgStream m_log;
  bool m_debug;

  MonGroup* m_monGroup;
  std::string m_rootDir;

  /// Trigger Tower container StoreGate key
  std::string m_triggerTowerLocation;

  /// Phi Units for eta/phi plots
  std::string m_phiUnits;
  /// Phi maximum in wanted units
  double m_phiMax;
  /// Phi scale to convert from radians to wanted units
  double m_phiScale;
  /// Number of events
  int m_events;
  /// Number of error event number samples to keep
  int m_eventSamples;
  /// Current event number
  int m_eventNumber;
  /// Sample event number counts
  std::vector<int> m_sampleCounts;

  //=======================
  //   Match/Mismatch plots
  //=======================

  // LUT
  TH2F* m_h_ppm_em_2d_etaPhi_tt_lut_SimEqData;
  TH2F* m_h_ppm_em_2d_etaPhi_tt_lut_SimNeData;
  TH2F* m_h_ppm_em_2d_etaPhi_tt_lut_SimNoData;
  TH2F* m_h_ppm_em_2d_etaPhi_tt_lut_DataNoSim;
  TH2F* m_h_ppm_had_2d_etaPhi_tt_lut_SimEqData;
  TH2F* m_h_ppm_had_2d_etaPhi_tt_lut_SimNeData;
  TH2F* m_h_ppm_had_2d_etaPhi_tt_lut_SimNoData;
  TH2F* m_h_ppm_had_2d_etaPhi_tt_lut_DataNoSim;
  
  // Mismatch Histograms
  TH2I* m_h_ppm_2d_LUT_MismatchEvents_cr0cr1;
  TH2I* m_h_ppm_2d_LUT_MismatchEvents_cr2cr3;
  TH2I* m_h_ppm_2d_LUT_MismatchEvents_cr4cr5;
  TH2I* m_h_ppm_2d_LUT_MismatchEvents_cr6cr7;

  void setLabelsCM(TH2* hist, bool xAxis = true, int first = 0);
  
 private:
  
 protected:
  
  StoreGateSvc* m_detStore;
 
  const CaloIdManager* m_caloMgr;
  const CaloLVL1_ID* m_lvl1Helper;
  const L1CaloTTIdTools* m_l1CaloTTIdTools;
  
  CaloTriggerTowerService* m_ttSvc;
  const TTOnlineID* m_l1ttonlineHelper;
    
};

#endif
