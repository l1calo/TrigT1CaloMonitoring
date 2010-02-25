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

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ServiceHandle.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"

class TH2F;
class StoreGateSvc;

class TrigT1CaloGlobalMonTool: public ManagedMonitorToolBase
{

public:
  
  TrigT1CaloGlobalMonTool(const std::string & type, const std::string & name,
  		          const IInterface* parent);
    

  virtual ~TrigT1CaloGlobalMonTool();

  virtual StatusCode initialize();
    
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
  enum CMMErrors { JEMCMMStatus, JEMCMMParity };
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
		     EnergyRoIMismatch };

  enum GlobalErrors { PPMDataStatus, PPMDataError, SubStatus, Parity, LinkDown,
                      RoIParity, Transmission, Simulation, CMMSubStatus,
		      GbCMMParity, CMMTransmission, CMMSimulation,
		      RODStatus, RODMissing, ROBStatus, Unpacking,
		      NumberOfGlobalErrors };

  typedef std::vector<int>            ErrorVector;

  TH2F* book2F(const std::string& name, const std::string& title,
                                    int nx, double xmin, double xmax,
                                    int ny, double ymin, double ymax);

  ServiceHandle<StoreGateSvc> m_storeGate;
  mutable MsgStream m_log;

  MonGroup* m_monGroup;
  
  /// Root directory
  std::string m_rootDir;

  //=======================
  //   Global Overview plot
  //=======================

  TH2F* m_h_global;

};

#endif
