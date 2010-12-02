// ********************************************************************
//
// NAME:     CPMSimBSMon.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Peter Faulkner
//	     
//
// ********************************************************************
#ifndef CPMSIMBSMON_H
#define CPMSIMBSMON_H

#include <map>
#include <string>
#include <vector>

#include "GaudiKernel/ToolHandle.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "DataModel/DataVector.h"

class LWHist;
class TH1F_LW;
class TH2F_LW;
class TH2I_LW;
class StatusCode;
class TrigT1CaloMonErrorTool;
class TrigT1CaloLWHistogramTool;

namespace LVL1 {
  class CPAlgorithm;
  class CPMTower;
  class CPMHits;
  class CMMCPHits;
  class CPMRoI;
  class RODHeader;
  class TriggerTower;
  class IL1EmTauTools;
  class IL1CPHitsTools;
}

class CPMSimBSMon: public ManagedMonitorToolBase
{

public:
  
  CPMSimBSMon(const std::string & type, const std::string & name,
		                        const IInterface* parent);
    

  virtual ~CPMSimBSMon();

  virtual StatusCode initialize();
  virtual StatusCode finalize();
    
  virtual StatusCode bookHistograms(bool isNewEventsBlock, bool isNewLumiBlock,
                                                           bool isNewRun);
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms(bool isEndOfEventsBlock,
                                    bool isEndOfLumiBlock, bool isEndOfRun);

private:

  enum SummaryErrors { EMTowerMismatch, HadTowerMismatch, RoIMismatch,
                       CPMHitsMismatch, CMMHitsMismatch, LocalSumMismatch,
		       RemoteSumMismatch, TotalSumMismatch,
		       NumberOfSummaryBins };

  typedef DataVector<LVL1::CPMTower>     CpmTowerCollection;
  typedef DataVector<LVL1::CPMHits>      CpmHitsCollection;
  typedef DataVector<LVL1::CMMCPHits>    CmmCpHitsCollection;
  typedef DataVector<LVL1::CPMRoI>       CpmRoiCollection;
  typedef DataVector<LVL1::TriggerTower> TriggerTowerCollection;
  typedef DataVector<LVL1::CPAlgorithm>  InternalRoiCollection;
  typedef DataVector<LVL1::RODHeader>    RodHeaderCollection;
  
  typedef std::vector<int> ErrorVector;

  typedef std::map<int, LVL1::TriggerTower*> TriggerTowerMap;
  typedef std::map<int, LVL1::CPMTower*>     CpmTowerMap;
  typedef std::map<int, LVL1::CPMRoI*>       CpmRoiMap;
  typedef std::map<int, LVL1::CPMHits*>      CpmHitsMap;
  typedef std::map<int, LVL1::CMMCPHits*>    CmmCpHitsMap;
  
  void  compare(const TriggerTowerMap& ttMap, const CpmTowerMap& cpMap,
                      ErrorVector& errors, bool overlap);
  void  compare(const CpmRoiMap& roiSimMap, const CpmRoiMap& roiMap,
                                                 ErrorVector& errors);
  void  compare(const CpmHitsMap& cpmSimMap, const CpmHitsMap& cpmMap,
                                             ErrorVector& errors);
  void  compare(const CpmHitsMap& cpmMap, const CmmCpHitsMap& cmmMap,
                      ErrorVector& errorsCPM, ErrorVector& errorsCMM);
  void  compare(const CmmCpHitsMap& cmmSimMap, const CmmCpHitsMap& cmmMap,
                                          ErrorVector& errors, int selection);
  void  setLabels(LWHist* hist, bool xAxis = true);
  void  setupMap(const TriggerTowerCollection* coll, TriggerTowerMap& map);
  void  setupMap(const CpmTowerCollection* coll, CpmTowerMap& map);
  void  setupMap(const CpmRoiCollection* coll, CpmRoiMap& map);
  void  setupMap(const CpmHitsCollection* coll, CpmHitsMap& map);
  void  setupMap(const CmmCpHitsCollection* coll, CmmCpHitsMap& map);
  void  simulate(const CpmTowerMap towers, const CpmTowerMap towersOv,
                       CpmRoiCollection* rois);
  void  simulate(const CpmRoiCollection* rois, CpmHitsCollection* hits);
  void  simulate(const CmmCpHitsCollection* hitsIn,
                       CmmCpHitsCollection* hitsOut, int selection);
  int   fpga(int crate, double phi);
  LVL1::CPMTower* ttCheck(LVL1::CPMTower* tt, CpmTowerCollection* coll);
  bool  limitedRoiSet(int crate);

  ToolHandle<LVL1::IL1EmTauTools>       m_emTauTool;
  ToolHandle<LVL1::IL1CPHitsTools>      m_cpHitsTool;
  ToolHandle<TrigT1CaloMonErrorTool>    m_errorTool;
  ToolHandle<TrigT1CaloLWHistogramTool> m_histTool;
  bool m_debug;

  /// Root directory
  std::string m_rootDir;

  /// CPM core tower container StoreGate key
  std::string m_cpmTowerLocation;
  /// CPM overlap tower container StoreGate key
  std::string m_cpmTowerLocationOverlap;
  /// CPM hits container StoreGate key
  std::string m_cpmHitsLocation;
  /// CMM-CP hits container StoreGate key
  std::string m_cmmCpHitsLocation;
  /// CPM RoI container StoreGate key
  std::string m_cpmRoiLocation;
  /// Trigger Tower container StoreGate key
  std::string m_triggerTowerLocation;
  /// ROD header container StoreGate key
  std::string m_rodHeaderLocation;
  /// Pointer to ROD header container
  const RodHeaderCollection* m_rodTES;

  /// CPM overlap tower container present
  bool m_overlapPresent;
  /// LimitedRoISet flags
  int m_limitedRoi;

  //=======================
  //   Match/Mismatch plots
  //=======================

  // CPM Towers
  TH2F_LW* m_h_EMTowerSIMeqDAT;
  TH2F_LW* m_h_EMTowerSIMneDAT;
  TH2F_LW* m_h_EMTowerSIMnoDAT;
  TH2F_LW* m_h_EMTowerDATnoSIM;
  TH2F_LW* m_h_HadTowerSIMeqDAT;
  TH2F_LW* m_h_HadTowerSIMneDAT;
  TH2F_LW* m_h_HadTowerSIMnoDAT;
  TH2F_LW* m_h_HadTowerDATnoSIM;
  TH2F_LW* m_h_EMTowerOvSIMeqDAT;
  TH2F_LW* m_h_EMTowerOvSIMneDAT;
  TH2F_LW* m_h_EMTowerOvSIMnoDAT;
  TH2F_LW* m_h_EMTowerOvDATnoSIM;
  TH2F_LW* m_h_HadTowerOvSIMeqDAT;
  TH2F_LW* m_h_HadTowerOvSIMneDAT;
  TH2F_LW* m_h_HadTowerOvSIMnoDAT;
  TH2F_LW* m_h_HadTowerOvDATnoSIM;
  TH2F_LW* m_h_FpgaTowerSIMeqDAT;
  TH2F_LW* m_h_FpgaTowerSIMneDAT;
  TH2F_LW* m_h_FpgaTowerSIMnoDAT;
  TH2F_LW* m_h_FpgaTowerDATnoSIM;

  // RoI
  TH2F_LW* m_h_RoISIMeqDAT;
  TH2F_LW* m_h_RoISIMneDAT;
  TH2F_LW* m_h_RoISIMnoDAT;
  TH2F_LW* m_h_RoIDATnoSIM;
  TH2F_LW* m_h_RoIThreshSIMeqDAT;
  TH2F_LW* m_h_RoIThreshSIMneDAT;
  TH2F_LW* m_h_RoIEtaPhiSIMeqDAT;
  TH2F_LW* m_h_RoIEtaPhiSIMneDAT;
  TH2F_LW* m_h_RoIEtaPhiSIMnoDAT;
  TH2F_LW* m_h_RoIEtaPhiDATnoSIM;

  // CPM Hits
  TH2F_LW* m_h_CPMHitsSIMeqDAT;
  TH2F_LW* m_h_CPMHitsSIMneDAT;
  TH2F_LW* m_h_CPMHitsSIMnoDAT;
  TH2F_LW* m_h_CPMHitsDATnoSIM;
  TH2F_LW* m_h_CPMHitsThreshSIMeqDAT;
  TH2F_LW* m_h_CPMHitsThreshSIMneDAT;

  // CMM-CP Hits
  TH2F_LW* m_h_CMMHitsSIMeqDAT;
  TH2F_LW* m_h_CMMHitsSIMneDAT;
  TH2F_LW* m_h_CMMHitsSIMnoDAT;
  TH2F_LW* m_h_CMMHitsDATnoSIM;
  TH2F_LW* m_h_CMMHitsThreshSIMeqDAT;
  TH2F_LW* m_h_CMMHitsThreshSIMneDAT;

  // CMM-CP Hit Sums
  TH1F_LW* m_h_SumsSIMeqDAT;
  TH1F_LW* m_h_SumsSIMneDAT;
  TH1F_LW* m_h_SumsSIMnoDAT;
  TH1F_LW* m_h_SumsDATnoSIM;
  TH2F_LW* m_h_SumsThreshSIMeqDAT;
  TH2F_LW* m_h_SumsThreshSIMneDAT;
  
  // Summary
  TH2F_LW* m_h_CPeqSIM;
  TH2F_LW* m_h_CPneSIM;
  TH1F_LW* m_h_CPneSIMSummary;

  // Mismatch Event Number Samples
  std::vector<TH2I_LW*> m_sampleHists;

};

#endif
