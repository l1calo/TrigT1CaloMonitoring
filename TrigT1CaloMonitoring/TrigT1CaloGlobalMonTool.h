// ********************************************************************
//
// NAME:     TrigT1CaloGlobalMonTool.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Peter Faulkner
//	     
//
// ********************************************************************
#ifndef TRIGT1CALOGLOBALMONTOOL_H
#define TRIGT1CALOGLOBALMONTOOL_H

#include <string>
#include <vector>

#include "GaudiKernel/ToolHandle.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"

class TH1F;
class TH2F;
class StatusCode;
class TrigT1CaloLWHistogramTool;

class TrigT1CaloGlobalMonTool: public ManagedMonitorToolBase
{

public:
  
  TrigT1CaloGlobalMonTool(const std::string & type, const std::string & name,
  		          const IInterface* parent);
    

  virtual ~TrigT1CaloGlobalMonTool();

  virtual StatusCode initialize();
  virtual StatusCode finalize();
    
  virtual StatusCode bookHistograms(bool isNewEventsBlock, bool isNewLumiBlock,
                                                           bool isNewRun);
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms(bool isEndOfEventsBlock,
                                    bool isEndOfLumiBlock, bool isEndOfRun);

private:

  // Enums for global summary plot

  // Hardware errors
  enum PPMErrors { DataStatus, DataError, PPMSubStatus };
  enum CPMErrors { CPMEMParity, CPMEMLink, CPMHadParity, CPMHadLink, CPMStatus,
                   CPMRoIParity, CMMCPParity, CMMCPStatus };
  enum JEMErrors { JEMEMParity, JEMHadParity, JEMEMLink, JEMHadLink, JEMStatus,
                   JEMRoIParity };
  enum CMMErrors { JEMCMMJetStatus, JEMCMMEnergyStatus, JEMCMMJetParity,
                   JEMCMMEnergyParity, JEMCMMRoIParity};
  enum RODErrors { GLink, /*CMMParity,*/ LVDSLink, FIFOOverflow, DataTransport,
                   Timeout, BCNMismatch, TriggerType, NoPayload, NoFragment,
	           ROBStatusError, UnpackingError };
  // Transmission/Comparison with simulation errors
  enum PPMMismatch { LUTMismatch };
  enum CPMMismatch { EMTowerMismatch, HadTowerMismatch, CPMRoIMismatch,
                     CPMHitsMismatch, CMMHitsMismatch, LocalSumMismatch,
		     RemoteSumMismatch, TotalSumMismatch };
  enum JEMMismatch { EMElementMismatch, HadElementMismatch, JEMRoIMismatch,
                     JEMHitsMismatch, CMMJetHitsMismatch, LocalJetMismatch,
		     RemoteJetMismatch, TotalJetMismatch, JetEtMismatch,
		     JetEtRoIMismatch, JEMEtSumsMismatch, CMMEtSumsMismatch,
		     LocalEnergyMismatch, RemoteEnergyMismatch,
		     TotalEnergyMismatch, SumEtMismatch, MissingEtMismatch,
		     MissingEtSigMismatch, EnergyRoIMismatch };

  enum GlobalErrors { PPMDataStatus, PPMDataError, SubStatus, Parity, LinkDown,
                      RoIParity, Transmission, Simulation, CMMSubStatus,
		      GbCMMParity, CMMTransmission, CMMSimulation,
		      RODStatus, RODMissing, ROBStatus, Unpacking,
		      NumberOfGlobalErrors };

  typedef std::vector<int> ErrorVector;

  TH2F* bookOverview(const std::string& name, const std::string& title);

  ToolHandle<TrigT1CaloLWHistogramTool> m_histTool;

  /// Root directory
  std::string m_rootDir;

  /// Threshold histogram pre-booking flags
  bool m_cpmThresh;
  bool m_jemThresh;
  bool m_cmmThresh;

  int m_recentLumi;
  bool m_onlineTest;
  unsigned int m_lumiNo;

  //========================
  //   Global Overview plots
  //========================

  TH2F* m_h_global;
  TH2F* m_h_current;
  TH2F* m_h_lumiblocks;
  TH1F* m_h_bylumi;
  TH1F* m_h_bytime;
  std::vector<TH2F*> m_v_lumi;

};

#endif
