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

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "DataModel/DataVector.h"

class TH1F;
class TH2F;
class TH2I;
class StoreGateSvc;

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
  
  TH1F* book1F(const std::string& name, const std::string& title,
                                        int nx, double xmin, double xmax);
  TH2F* book2F(const std::string& name, const std::string& title,
                                        int nx, double xmin, double xmax,
                                        int ny, double ymin, double ymax);
  TH2F* book2F(const std::string& name, const std::string& title,
                                        int nx, const double* xbins,
                                        int ny, double ymin, double ymax);
  TH2I* book2I(const std::string& name, const std::string& title,
                                        int nx, double xmin, double xmax,
                                        int ny, double ymin, double ymax);
  TH2F* bookEtaPhi(const std::string& name, const std::string& title,
                                            bool isRoi = false);
  void  compare(const JetElementMap& jeSimMap, const JetElementMap& jeMap,
                      ErrorVector& errors, bool overlap);
  void  compare(const JemRoiMap& roiSimMap, const JemRoiMap& roiMap,
                const RodHeaderCollection* rods, ErrorVector& errors);
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
  void  setLabels(TH2* hist);
  void  setLabelsCMFC(TH2* hist);
  void  setLabelsJMS(TH2* hist);
  void  setLabelsCMT(TH2* hist);
  void  setLabelsYNUM(TH2* hist, int beg, int end);
  void  setLabelsXNUM(TH2* hist, int beg, int end);
  void  setLabelsJEM(TH2* hist, bool xAxis = true);
  void  setLabelsMC(TH2* hist);
  void  setLabelsSH(TH1* hist);
  void  setLabelsSHF(TH2* hist);
  void  setLabelsJES(TH2* hist);
  void  setLabelsEnTot(TH2* hist);
  void  setLabelsEnTotThr(TH2* hist);
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

  ServiceHandle<StoreGateSvc>          m_storeGate;
  ToolHandle<LVL1::IL1JEPHitsTools>    m_jepHitsTool;
  ToolHandle<LVL1::IL1JetTools>        m_jetTool;
  ToolHandle<LVL1::IL1JetElementTools> m_jetElementTool;
  ToolHandle<LVL1::IL1JEPEtSumsTools>  m_etSumsTool;
  mutable MsgStream m_log;
  bool m_debug;

  MonGroup* m_monGroup;
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

  /// Phi scale to convert from radians to histogram units
  double m_phiScale;
  /// Number of events
  int m_events;
  /// Simulation allowed flag
  bool m_compareWithSim;
  /// Number of error event number samples to keep
  int m_eventSamples;
  /// Current event number
  int m_eventNumber;
  /// Sample event number counts
  std::vector<int> m_sampleCounts;

  //=======================
  //   Match/Mismatch plots
  //=======================

  // Jet Elements
  TH2F* m_h_EMEleSIMeqDAT;
  TH2F* m_h_EMEleSIMneDAT;
  TH2F* m_h_EMEleSIMnoDAT;
  TH2F* m_h_EMEleDATnoSIM;
  TH2F* m_h_HadEleSIMeqDAT;
  TH2F* m_h_HadEleSIMneDAT;
  TH2F* m_h_HadEleSIMnoDAT;
  TH2F* m_h_HadEleDATnoSIM;
  TH2F* m_h_EMEleOvSIMeqDAT;
  TH2F* m_h_EMEleOvSIMneDAT;
  TH2F* m_h_EMEleOvSIMnoDAT;
  TH2F* m_h_EMEleOvDATnoSIM;
  TH2F* m_h_HadEleOvSIMeqDAT;
  TH2F* m_h_HadEleOvSIMneDAT;
  TH2F* m_h_HadEleOvSIMnoDAT;
  TH2F* m_h_HadEleOvDATnoSIM;

  // RoIs
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

  // JEM Hits
  TH2F* m_h_JEMHitsSIMeqDAT;
  TH2F* m_h_JEMHitsSIMneDAT;
  TH2F* m_h_JEMHitsSIMnoDAT;
  TH2F* m_h_JEMHitsDATnoSIM;
  TH2F* m_h_JEMHitsThreshSIMeqDAT;
  TH2F* m_h_JEMHitsThreshSIMneDAT;

  // CMM-Jet Hits
  TH2F* m_h_CMMHitsSIMeqDAT;
  TH2F* m_h_CMMHitsSIMneDAT;
  TH2F* m_h_CMMHitsSIMnoDAT;
  TH2F* m_h_CMMHitsDATnoSIM;
  TH2F* m_h_CMMHitsThreshSIMeqDAT;
  TH2F* m_h_CMMHitsThreshSIMneDAT;

  // CMM-Jet Hit Sums
  TH1F* m_h_SumsSIMeqDAT;
  TH1F* m_h_SumsSIMneDAT;
  TH1F* m_h_SumsSIMnoDAT;
  TH1F* m_h_SumsDATnoSIM;
  TH2F* m_h_SumsThreshSIMeqDAT;
  TH2F* m_h_SumsThreshSIMneDAT;

  // JEMEtSums
  TH2F* m_h_jemEtSumsSIMeqDAT;
  TH2F* m_h_jemEtSumsSIMneDAT;
  TH2F* m_h_jemEtSumsSIMnoDAT;
  TH2F* m_h_jemEtSumsDATnoSIM;

  // CMMEtSums
  TH2F* m_h_cmmEtSumsSIMeqDAT;
  TH2F* m_h_cmmEtSumsSIMneDAT;
  TH2F* m_h_cmmEtSumsSIMnoDAT;
  TH2F* m_h_cmmEtSumsDATnoSIM;

  // Energy Crate/System sums
  TH2F* m_h_EnSumsSIMeqDAT;
  TH2F* m_h_EnSumsSIMneDAT;
  TH2F* m_h_EnSumsSIMnoDAT;
  TH2F* m_h_EnSumsDATnoSIM;
  TH2F* m_h_EnSumsThreshSIMeqDAT;
  TH2F* m_h_EnSumsThreshSIMneDAT;

  // Summary
  TH2F* m_h_JEPeqSIM;
  TH2F* m_h_JEPneSIM;
  TH1F* m_h_JEPneSIMSummary;

  // Mismatch Event Number Samples
  std::vector<TH2I*> m_sampleHists;

};

#endif
