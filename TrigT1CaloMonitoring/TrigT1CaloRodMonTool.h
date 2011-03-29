// ********************************************************************
//
// NAME:     TrigT1CaloRodMonTool.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Peter Faulkner
//	     
//
// ********************************************************************
#ifndef TRIGT1CALORODMONTOOL_H
#define TRIGT1CALORODMONTOOL_H

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
  class RODHeader;
}

class TrigT1CaloRodMonTool: public ManagedMonitorToolBase
{

public:
  
  TrigT1CaloRodMonTool(const std::string & type, const std::string & name,
		       const IInterface* parent);
    

  virtual ~TrigT1CaloRodMonTool();

  virtual StatusCode initialize();
    
  virtual StatusCode bookHistograms(bool isNewEventsBlock, bool isNewLumiBlock,
                                                           bool isNewRun);
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms(bool isEndOfEventsBlock,
                                    bool isEndOfLumiBlock, bool isEndOfRun);

private:

  enum StatusBits { GLink, /*CMMParity,*/ LVDSLink, FIFOOverflow, DataTransport,
                    Timeout, BCNMismatch, TriggerType, LimitedRoI, NoFragment,
		    NumberOfStatusBins, NoPayload = LimitedRoI,
		    ROBStatusError = NumberOfStatusBins, UnpackingError };

  typedef DataVector<LVL1::RODHeader> RodHeaderCollection;
  typedef std::vector<unsigned int>   ROBErrorCollection;
  typedef std::vector<int>            ErrorVector;
  
  void setLabelsStatus(LWHist* hist, bool xAxis = true);
  void setLabelsROBStatusGen(LWHist* hist, bool xAxis = true);
  void setLabelsROBStatusSpec(LWHist* hist, bool xAxis = true);
  void setLabelsUnpacking(LWHist* hist, bool xAxis = true);

  ToolHandle<TrigT1CaloMonErrorTool>    m_errorTool;
  ToolHandle<TrigT1CaloLWHistogramTool> m_histTool;

  /// DAQ ROD header container StoreGate key
  std::string m_rodHeaderLocation;
  /// CP RoIB ROD header container StoreGate key
  std::string m_cpRoibRodHeaderLocation;
  /// JEP RoIB ROD header container StoreGate key
  std::string m_jepRoibRodHeaderLocation;
  /// ROB and Unpacking Error vector StoreGate key
  std::string m_robErrorVectorLocation;
  
  /// Root directory
  std::string m_rootDir;

  /// Accumulated payload sizes
  std::vector<double> m_sumPayloads1;
  std::vector<double> m_sumPayloads2;

  /// Number of events
  int m_events;
  /// Test online code flag
  bool m_onlineTest;
  /// Histograms booked flag
  bool m_histBooked;

  //=======================
  //   Payload plots
  //=======================

  TH1F_LW* m_h_ROD_PP;
  TH1F_LW* m_h_ROD_CP;
  TH1F_LW* m_h_ROD_JEP;
  TH1F_LW* m_h_ROD_RoI;

  //=======================
  //   Status bit plots
  //=======================

  TH2F_LW* m_h_ROD_PP_stat;
  TH2F_LW* m_h_ROD_CPJEP_stat;
  TH2F_LW* m_h_ROD_RoI_stat;
  TH2F_LW* m_h_ROD_PP_robgen;
  TH2F_LW* m_h_ROD_CPJEP_robgen;
  TH2F_LW* m_h_ROD_RoI_robgen;
  TH2F_LW* m_h_ROD_PP_robspec;
  TH2F_LW* m_h_ROD_CPJEP_robspec;
  TH2F_LW* m_h_ROD_RoI_robspec;
  TH2F_LW* m_h_ROD_PP_unp;
  TH2F_LW* m_h_ROD_CPJEP_unp;
  TH2F_LW* m_h_ROD_RoI_unp;

  //=======================
  //   Summary plots
  //=======================

  TH1F_LW* m_h_ROD_summary;
  TH1F_LW* m_h_ROB_summary;
  TH1F_LW* m_h_Unp_summary;
  TH2I_LW* m_h_ROD_events;
  TH2I_LW* m_h_ROB_events;
  TH2I_LW* m_h_Unp_events;

};

#endif
