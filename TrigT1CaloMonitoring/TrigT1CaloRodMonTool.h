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

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "DataModel/DataVector.h"

class TH1;
class TH1F;
class TH2F;
class StoreGateSvc;
class TrigT1CaloMonErrorTool;

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

  // Enums for global summary plot

  // Hardware errors
  enum PPMErrors { DataStatus, DataError, PPMSubStatus };
  enum CPMErrors { CPMEMParity, CPMEMLink, CPMHadParity, CPMHadLink, CPMStatus,
                   CPMRoIParity, CMMCPParity, CMMCPStatus };
  enum JEMErrors { JEMEMParity, JEMHadParity, JEMEMLink, JEMHadLink, JEMStatus,
                   JEMRoIParity };
  enum CMMErrors { JEMCMMStatus, JEMCMMParity };
  // Transmission/Comparison with simulation errors
  enum CPMMismatch { EMTowerMismatch, HadTowerMismatch, CPMRoIMismatch,
                     CPMHitsMismatch, CMMHitsMismatch, LocalSumMismatch,
		     RemoteSumMismatch, TotalSumMismatch };
  enum JEMMismatch { EMElementMismatch, HadElementMismatch, JEMRoIMismatch,
                     JEMHitsMismatch, CMMJetHitsMismatch, LocalJetMismatch,
		     RemoteJetMismatch, TotalJetMismatch, JetEtMismatch,
		     JetEtRoIMismatch, JEMEtSumsMismatch, CMMEtSumsMismatch,
		     LocalEnergyMismatch, RemoteEnergyMismatch,
		     TotalEnergyMismatch, SumEtMismatch, MissingEtMismatch,
		     EnergyRoIMismatch };

  enum GlobalErrors { PPMDataStatus, PPMDataError, SubStatus, Parity, LinkDown,
                      RoIParity, Transmission, Simulation, CMMSubStatus,
		      GbCMMParity, CMMTransmission, CMMSimulation,
		      RODStatus, RODMissing, ROBStatus, Unpacking,
		      NumberOfGlobalErrors };

  typedef DataVector<LVL1::RODHeader> RodHeaderCollection;
  typedef std::vector<unsigned int>   ROBErrorCollection;
  typedef std::vector<int>            ErrorVector;
  
  TH1F* book1F(const std::string& name, const std::string& title,
                                    int nx, double xmin, double xmax);
  TH2F* book2F(const std::string& name, const std::string& title,
                                    int nx, double xmin, double xmax,
                                    int ny, double ymin, double ymax);
  void setLabelsStatus(TH1* hist);
  void setLabelsCSL(TH1* hist, bool xAxis, int firstBin, int lastBin,
                                           int binIncr, int slinkIncr);
  void setLabelsROBStatusGen(TH1* hist);
  void setLabelsROBStatusSpec(TH1* hist);
  void setLabelsUnpacking(TH1* hist);

  ServiceHandle<StoreGateSvc> m_storeGate;
  ToolHandle<TrigT1CaloMonErrorTool> m_errorTool;
  mutable MsgStream m_log;

  MonGroup* m_monGroup;

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

  //=======================
  //   Payload plots
  //=======================

  TH1F* m_h_ROD_PP;
  TH1F* m_h_ROD_CP;
  TH1F* m_h_ROD_JEP;
  TH1F* m_h_ROD_RoI;

  //=======================
  //   Status bit plots
  //=======================

  TH2F* m_h_ROD_PP_stat;
  TH2F* m_h_ROD_CPJEP_stat;
  TH2F* m_h_ROD_RoI_stat;
  TH2F* m_h_ROD_PP_robgen;
  TH2F* m_h_ROD_CPJEP_robgen;
  TH2F* m_h_ROD_RoI_robgen;
  TH2F* m_h_ROD_PP_robspec;
  TH2F* m_h_ROD_CPJEP_robspec;
  TH2F* m_h_ROD_RoI_robspec;
  TH2F* m_h_ROD_PP_unp;
  TH2F* m_h_ROD_CPJEP_unp;
  TH2F* m_h_ROD_RoI_unp;

  //=======================
  //   Summary plots
  //=======================

  TH1F* m_h_ROD_summary;
  TH1F* m_h_ROB_summary;
  TH1F* m_h_Unp_summary;
  TH2F* m_h_global;

};

#endif
