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

#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "DataModel/DataVector.h"

class TH1F;
class TH2F;
class StoreGateSvc;

namespace LVL1 {
  class CPAlgorithm;
  class CPMTower;
  class CPMHits;
  class CMMCPHits;
  class CPMRoI;
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
  typedef DataVector<LVL1::CPAlgorithm>  InternalRoiCollection;
  
  typedef std::vector<int> ErrorVector;

  typedef std::map<int, LVL1::TriggerTower*> TriggerTowerMap;
  typedef std::map<int, LVL1::CPMTower*>     CpmTowerMap;
  typedef std::map<int, LVL1::CPMRoI*>       CpmRoiMap;
  typedef std::map<int, LVL1::CPMHits*>      CpmHitsMap;
  typedef std::map<int, LVL1::CMMCPHits*>    CmmCpHitsMap;
  
  TH1F* book1F(std::string nam, std::string tit,
                                    int nx, double xmin, double xmax);
  TH2F* book2F(std::string nam, std::string tit,
                                    int nx, double xmin, double xmax,
                                    int ny, double ymin, double ymax);
  void  compare(const TriggerTowerMap& ttMap, const CpmTowerMap& cpMap,
                      ErrorVector& errors, ErrorVector& errors2);
  void  compare(const CpmRoiMap& roiSimMap, const CpmRoiMap& roiMap,
                                            ErrorVector& errors);
  void  compare(const CpmHitsMap& cpmSimMap, const CpmHitsMap& cpmMap,
                                             ErrorVector& errors);
  void  compare(const CpmHitsMap& cpmMap, const CmmCpHitsMap& cmmMap,
                      ErrorVector& errors, ErrorVector& errors2);
  void  compare(const CmmCpHitsMap& cmmSimMap, const CmmCpHitsMap& cmmMap,
                                          ErrorVector& errors, int selection);
  void  setLabels(TH2* hist);
  void  setLabelsCMCC(TH2* hist);
  void  setLabelsCMT(TH2* hist);
  void  setLabelsT(TH2* hist);
  void  setLabelsCPM(TH2* hist);
  void  setLabelsMC(TH2* hist);
  void  setLabelsMCLR(TH2* hist);
  void  setLabelsSLR(TH1* hist);
  void  setLabelsST(TH2* hist);
  void  setLabelsSRLR(TH1* hist);
  void  setLabelsSRT(TH2* hist);
  void  setupMap(const TriggerTowerCollection* coll, TriggerTowerMap& map);
  void  setupMap(const CpmTowerCollection* coll, CpmTowerMap& map);
  void  setupMap(const CpmRoiCollection* coll, CpmRoiMap& map);
  void  setupMap(const CpmHitsCollection* coll, CpmHitsMap& map);
  void  setupMap(const CmmCpHitsCollection* coll, CmmCpHitsMap& map);
  void  simulate(const CpmTowerMap towers, CpmRoiCollection* rois);
  void  simulate(const CpmRoiCollection* rois, CpmHitsCollection* hits);
  void  simulate(const CmmCpHitsCollection* hitsIn,
                       CmmCpHitsCollection* hitsOut, int selection);

  ServiceHandle<StoreGateSvc> m_storeGate;
  ToolHandle<LVL1::IL1EmTauTools> m_emTauTool;
  ToolHandle<LVL1::IL1CPHitsTools> m_cpHitsTool;

  MonGroup* m_monGroup;
  std::string m_rootDir;

  /// CPM tower container StoreGate key
  std::string m_cpmTowerLocation;
  /// CPM hits container StoreGate key
  std::string m_cpmHitsLocation;
  /// CMM-CP hits container StoreGate key
  std::string m_cmmCpHitsLocation;
  /// CPM RoI container StoreGate key
  std::string m_cpmRoiLocation;
  /// CPM RoIB container StoreGate key
  std::string m_cpmRoiLocationRoib;
  /// Trigger Tower container StoreGate key
  std::string m_triggerTowerLocation;

  /// Phi Units for eta/phi plots
  std::string m_phiUnits;
  /// Phi maximum in wanted units
  double m_phiMax;
  /// Phi scale to convert from radians to wanted units
  double m_phiScale;
  /// Simulation allowed flag
  bool m_compareWithSim;

  //=======================
  //   Match/Mismatch plots
  //=======================

  // CPM Towers
  TH2F* m_h_EMTowerSIMeqDAT;
  TH2F* m_h_EMTowerSIMneDAT;
  TH2F* m_h_EMTowerSIMnoDAT;
  TH2F* m_h_EMTowerDATnoSIM;
  TH2F* m_h_HadTowerSIMeqDAT;
  TH2F* m_h_HadTowerSIMneDAT;
  TH2F* m_h_HadTowerSIMnoDAT;
  TH2F* m_h_HadTowerDATnoSIM;

  // RoI
  TH2F* m_h_RoISIMeqDAT;
  TH2F* m_h_RoISIMneDAT;
  TH2F* m_h_RoISIMnoDAT;
  TH2F* m_h_RoIDATnoSIM;
  TH2F* m_h_RoIThreshSIMeqDAT;
  TH2F* m_h_RoIThreshSIMneDAT;
  TH2F* m_h_RoIEtaPhiSIMeqDAT;
  TH2F* m_h_RoIEtaPhiSIMneDAT;
  TH2F* m_h_RoIEtaPhiSIMnoDAT;
  TH2F* m_h_RoIEtaPhiDATnoSIM;

  // CPM Hits
  TH2F* m_h_CPMHitsSIMeqDAT;
  TH2F* m_h_CPMHitsSIMneDAT;
  TH2F* m_h_CPMHitsSIMnoDAT;
  TH2F* m_h_CPMHitsDATnoSIM;
  TH2F* m_h_CPMHitsThreshSIMeqDAT;
  TH2F* m_h_CPMHitsThreshSIMneDAT;

  // CMM-CP Hits
  TH2F* m_h_CMMHitsSIMeqDAT;
  TH2F* m_h_CMMHitsSIMneDAT;
  TH2F* m_h_CMMHitsSIMnoDAT;
  TH2F* m_h_CMMHitsDATnoSIM;
  TH2F* m_h_CMMHitsThreshSIMeqDAT;
  TH2F* m_h_CMMHitsThreshSIMneDAT;

  // CMM-CP Hit Sums
  TH1F* m_h_SumsSIMeqDAT;
  TH1F* m_h_SumsSIMneDAT;
  TH1F* m_h_SumsSIMnoDAT;
  TH1F* m_h_SumsDATnoSIM;
  TH2F* m_h_SumsThreshSIMeqDAT;
  TH2F* m_h_SumsThreshSIMneDAT;
  
  // Summary
  TH2F* m_h_CPeqSIM;
  TH2F* m_h_CPneSIM;
  TH1F* m_h_CPneSIMSummary;

};

#endif
