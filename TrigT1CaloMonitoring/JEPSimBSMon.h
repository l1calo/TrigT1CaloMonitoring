// ********************************************************************
//
// NAME:     JEPSimBSMon.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Peter Faulkner
//	     
//
// ********************************************************************
#ifndef JEPSIMBSMON_H
#define JEPSIMBSMON_H

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
  class CMMEtSums;
  class CMMJetHits;
  class CMMRoI;
  class JEMEtSums;
  class JEMHits;
  class JEMRoI;
  class JetAlgorithm;
  class JetElement;
  class RODHeader;
  class TriggerTower;
  class IL1JEPHitsTools;
  class IL1JetElementTools;
  class IL1JetTools;
  class IL1JEPEtSumsTools;
}

class JEPSimBSMon: public ManagedMonitorToolBase
{

public:
  
  JEPSimBSMon(const std::string & type, const std::string & name,
		       const IInterface* parent);
    

  virtual ~JEPSimBSMon();

  virtual StatusCode initialize();
    
  virtual StatusCode bookHistograms(bool isNewEventsBlock, bool isNewLumiBlock,
                                                           bool isNewRun);
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms(bool isEndOfEventsBlock,
                                    bool isEndOfLumiBlock, bool isEndOfRun);

private:

  enum SummaryErrors { EMElementMismatch, HadElementMismatch, RoIMismatch,
                       JEMHitsMismatch, CMMJetHitsMismatch, LocalJetMismatch,
		       RemoteJetMismatch, TotalJetMismatch, JetEtMismatch,
		       JetEtRoIMismatch, JEMEtSumsMismatch, CMMEtSumsMismatch,
		       LocalEnergyMismatch, RemoteEnergyMismatch,
		       TotalEnergyMismatch, SumEtMismatch, MissingEtMismatch,
		       EnergyRoIMismatch, NumberOfSummaryBins };

  typedef DataVector<LVL1::JetElement>   JetElementCollection;
  typedef DataVector<LVL1::JEMHits>      JemHitsCollection;
  typedef DataVector<LVL1::CMMJetHits>   CmmJetHitsCollection;
  typedef DataVector<LVL1::JEMRoI>       JemRoiCollection;
  typedef DataVector<LVL1::TriggerTower> TriggerTowerCollection;
  typedef DataVector<LVL1::JetAlgorithm> InternalRoiCollection;
  typedef DataVector<LVL1::JEMEtSums>    JemEtSumsCollection;
  typedef DataVector<LVL1::CMMEtSums>    CmmEtSumsCollection;
  typedef DataVector<LVL1::RODHeader>    RodHeaderCollection;
  
  typedef std::vector<int> ErrorVector;

  typedef std::map<int, LVL1::JetElement*>   JetElementMap;
  typedef std::map<int, LVL1::JEMRoI*>       JemRoiMap;
  typedef std::map<int, LVL1::JEMHits*>      JemHitsMap;
  typedef std::map<int, LVL1::CMMJetHits*>   CmmJetHitsMap;
  typedef std::map<int, LVL1::JEMEtSums*>    JemEtSumsMap;
  typedef std::map<int, LVL1::CMMEtSums*>    CmmEtSumsMap;
  
  void  compare(const JetElementMap& jeSimMap, const JetElementMap& jeMap,
                      ErrorVector& errors, bool overlap);
  void  compare(const JemRoiMap& roiSimMap, const JemRoiMap& roiMap,
                                                 ErrorVector& errors);
  void  compare(const JemHitsMap& jemSimMap, const JemHitsMap& jemMap,
                                             ErrorVector& errors);
  void  compare(const JemHitsMap& jemMap, const CmmJetHitsMap& cmmMap,
                      ErrorVector& errorsJEM, ErrorVector& errorsCMM);
  void  compare(const CmmJetHitsMap& cmmSimMap, const CmmJetHitsMap& cmmMap,
                                          ErrorVector& errors, int selection);
  void  compare(const CmmJetHitsMap& cmmMap, const LVL1::CMMRoI* cmmRoi,
                                             ErrorVector& errors);
  void  compare(const JemEtSumsMap& jemSimMap, const JemEtSumsMap& jemMap,
                                               ErrorVector& errors);
  void  compare(const JemEtSumsMap& jemMap, const CmmEtSumsMap& cmmMap,
                      ErrorVector& errorsJEM, ErrorVector& errorsCMM);
  void  compare(const CmmEtSumsMap& cmmSimMap, const CmmEtSumsMap& cmmMap,
                                          ErrorVector& errors, int selection);
  void  compare(const CmmEtSumsMap& cmmMap, const LVL1::CMMRoI* cmmRoi,
                                              ErrorVector& errors);
  void  fillEventSample(int err, int loc, bool isJem);
  void  setLabels(LWHist* hist);
  void  setLabelsSH(LWHist* hist);
  void  setLabelsSHF(LWHist* hist);
  void  setLabelsEnTot(LWHist* hist);
  void  setLabelsEnTotThr(LWHist* hist);
  void  setupMap(const JetElementCollection* coll, JetElementMap& map);
  void  setupMap(const JemRoiCollection* coll, JemRoiMap& map);
  void  setupMap(const JemHitsCollection* coll, JemHitsMap& map);
  void  setupMap(const CmmJetHitsCollection* coll, CmmJetHitsMap& map);
  void  setupMap(const JemEtSumsCollection* coll, JemEtSumsMap& map);
  void  setupMap(const CmmEtSumsCollection* coll, CmmEtSumsMap& map);
  void  simulate(const TriggerTowerCollection* towers,
                       JetElementCollection* elements);
  void  simulate(const JetElementCollection* elements,
                 const JetElementCollection* elementsOv,
		       JemRoiCollection* rois);
  void  simulate(const JemRoiCollection* rois, JemHitsCollection* hits);
  void  simulate(const CmmJetHitsCollection* hitsIn,
                       CmmJetHitsCollection* hitsOut, int selection);
  void  simulate(const JetElementCollection* elements,
                       JemEtSumsCollection* sums);
  void  simulate(const CmmEtSumsCollection* sumsIn,
                       CmmEtSumsCollection* sumsOut, int selection);
  bool  limitedRoiSet(int crate);

  ToolHandle<LVL1::IL1JEPHitsTools>      m_jepHitsTool;
  ToolHandle<LVL1::IL1JetTools>          m_jetTool;
  ToolHandle<LVL1::IL1JetElementTools>   m_jetElementTool;
  ToolHandle<LVL1::IL1JEPEtSumsTools>    m_etSumsTool;
  ToolHandle<TrigT1CaloMonErrorTool>     m_errorTool;
  ToolHandle<TrigT1CaloLWHistogramTool>  m_histTool;

  bool m_debug;
  std::string m_rootDir;

  /// Core Jet Element container StoreGate key
  std::string m_jetElementLocation;
  /// Overlap Jet Element container StoreGate key
  std::string m_jetElementLocationOverlap;
  /// JEM hits container StoreGate key
  std::string m_jemHitsLocation;
  /// CMM-Jet hits container StoreGate key
  std::string m_cmmJetHitsLocation;
  /// JEM RoI container StoreGate key
  std::string m_jemRoiLocation;
  /// CMM RoI container StoreGate key
  std::string m_cmmRoiLocation;
  /// JEM Et sums container StoreGate key
  std::string m_jemEtSumsLocation;
  /// CMM Et sums container StoreGate key
  std::string m_cmmEtSumsLocation;
  /// Trigger Tower container StoreGate key
  std::string m_triggerTowerLocation;
  /// ROD header container StoreGate key
  std::string m_rodHeaderLocation;
  /// Pointer to ROD header container
  const RodHeaderCollection* m_rodTES;
  /// LimitedRoISet flags
  int m_limitedRoi;

  //=======================
  //   Match/Mismatch plots
  //=======================

  // Jet Elements
  TH2F_LW* m_h_EMEleSIMeqDAT;
  TH2F_LW* m_h_EMEleSIMneDAT;
  TH2F_LW* m_h_EMEleSIMnoDAT;
  TH2F_LW* m_h_EMEleDATnoSIM;
  TH2F_LW* m_h_HadEleSIMeqDAT;
  TH2F_LW* m_h_HadEleSIMneDAT;
  TH2F_LW* m_h_HadEleSIMnoDAT;
  TH2F_LW* m_h_HadEleDATnoSIM;
  TH2F_LW* m_h_EMEleOvSIMeqDAT;
  TH2F_LW* m_h_EMEleOvSIMneDAT;
  TH2F_LW* m_h_EMEleOvSIMnoDAT;
  TH2F_LW* m_h_EMEleOvDATnoSIM;
  TH2F_LW* m_h_HadEleOvSIMeqDAT;
  TH2F_LW* m_h_HadEleOvSIMneDAT;
  TH2F_LW* m_h_HadEleOvSIMnoDAT;
  TH2F_LW* m_h_HadEleOvDATnoSIM;

  // RoIs
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

  // JEM Hits
  TH2F_LW* m_h_JEMHitsSIMeqDAT;
  TH2F_LW* m_h_JEMHitsSIMneDAT;
  TH2F_LW* m_h_JEMHitsSIMnoDAT;
  TH2F_LW* m_h_JEMHitsDATnoSIM;
  TH2F_LW* m_h_JEMHitsThreshSIMeqDAT;
  TH2F_LW* m_h_JEMHitsThreshSIMneDAT;

  // CMM-Jet Hits
  TH2F_LW* m_h_CMMHitsSIMeqDAT;
  TH2F_LW* m_h_CMMHitsSIMneDAT;
  TH2F_LW* m_h_CMMHitsSIMnoDAT;
  TH2F_LW* m_h_CMMHitsDATnoSIM;
  TH2F_LW* m_h_CMMHitsThreshSIMeqDAT;
  TH2F_LW* m_h_CMMHitsThreshSIMneDAT;

  // CMM-Jet Hit Sums
  TH1F_LW* m_h_SumsSIMeqDAT;
  TH1F_LW* m_h_SumsSIMneDAT;
  TH1F_LW* m_h_SumsSIMnoDAT;
  TH1F_LW* m_h_SumsDATnoSIM;
  TH2F_LW* m_h_SumsThreshSIMeqDAT;
  TH2F_LW* m_h_SumsThreshSIMneDAT;

  // JEMEtSums
  TH2F_LW* m_h_jemEtSumsSIMeqDAT;
  TH2F_LW* m_h_jemEtSumsSIMneDAT;
  TH2F_LW* m_h_jemEtSumsSIMnoDAT;
  TH2F_LW* m_h_jemEtSumsDATnoSIM;

  // CMMEtSums
  TH2F_LW* m_h_cmmEtSumsSIMeqDAT;
  TH2F_LW* m_h_cmmEtSumsSIMneDAT;
  TH2F_LW* m_h_cmmEtSumsSIMnoDAT;
  TH2F_LW* m_h_cmmEtSumsDATnoSIM;

  // Energy Crate/System sums
  TH2F_LW* m_h_EnSumsSIMeqDAT;
  TH2F_LW* m_h_EnSumsSIMneDAT;
  TH2F_LW* m_h_EnSumsSIMnoDAT;
  TH2F_LW* m_h_EnSumsDATnoSIM;
  TH2F_LW* m_h_EnSumsThreshSIMeqDAT;
  TH2F_LW* m_h_EnSumsThreshSIMneDAT;

  // Summary
  TH2F_LW* m_h_JEPeqSIM;
  TH2F_LW* m_h_JEPneSIM;
  TH1F_LW* m_h_JEPneSIMSummary;

  // Mismatch Event Number Samples
  std::vector<TH2I_LW*> m_sampleHists;

};

#endif
