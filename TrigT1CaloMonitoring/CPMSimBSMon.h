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
                                              ErrorVector& errors);
  void  compare(const CpmRoiMap& roiSimMap, const CpmRoiMap& roiMap,
                                            ErrorVector& errors);
  void  compare(const CpmHitsMap& cpmSimMap, const CpmHitsMap& cpmMap,
                                             ErrorVector& errors);
  void  compare(const CpmHitsMap& cpmMap, const CmmCpHitsMap& cmmMap,
                      ErrorVector& errors, ErrorVector& errors2);
  void  compare(const CmmCpHitsMap& cmmSimMap, const CmmCpHitsMap& cmmMap,
                                          ErrorVector& errors, int selection);
  void  setLabels(TH2* hist);
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

  //=======================
  //   Match/Mismatch plots
  //=======================

  TH2F* m_h_CPeqSIM;
  TH2F* m_h_CPneSIM;

};

#endif
