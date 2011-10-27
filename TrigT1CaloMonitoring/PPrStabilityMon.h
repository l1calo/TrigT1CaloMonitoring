// ********************************************************************
//
// NAME:        PPrStabilityMon.h
// PACKAGE:     TrigT1CaloMonitoring
//
// Author:      Rohin T Narayan (narayan@physi.uni-heidelberg.de)
//              Universitaet Heidelberg
//
// ********************************************************************

/**
 * This class does online monitoring of "FineTime"
 * for each TriggerTower
 * 
 * The class uses Trigger/TrigT1/TrigT1CaloCalibTools/L1CaloPprFineTimePlotManager
 * tool to generate the monitoring histogram.
 * */
#ifndef TRIGT1CALOMONITORING_PPRSTABILITYMON_H
#define TRIGT1CALOMONITORING_PPRSTABILITYMON_H

#include <string>
#include <vector>
#include <map>

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "GaudiKernel/ToolHandle.h"

class TH1F_LW;
class TH2F_LW;
class TH2I_LW;
class TProfile_LW;
class TProfile2D_LW;

class StatusCode;
class EventInfo;

class TrigT1CaloMonErrorTool;
class TrigT1CaloLWHistogramTool;
class L1CaloPprFineTimePlotManager;

namespace LVL1 {
  class IL1TriggerTowerTool;
}

class PPrStabilityMon: public ManagedMonitorToolBase
{

 public:
  
  PPrStabilityMon(const std::string & type, const std::string & name,const IInterface* parent);
  virtual ~PPrStabilityMon();

  virtual StatusCode initialize();
  virtual StatusCode finalize();
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms(bool isEndofEventsBlock, bool isEndofLumiBlock, bool isEndofRun);

private:
  unsigned int m_ppmADCMinValue;
  unsigned int m_lumiBlock;
  unsigned int m_lumiBlockMax;

  ServiceHandle<StoreGateSvc>             m_storeGate;
  // Tool to retrieve bytestream errors
  ToolHandle<TrigT1CaloMonErrorTool>      m_errorTool;
  ToolHandle<TrigT1CaloLWHistogramTool>   m_histTool;
  ToolHandle<LVL1::IL1TriggerTowerTool>   m_ttTool; 
  L1CaloPprFineTimePlotManager*           m_plotManager;

  std::string m_TriggerTowerContainerName;
  std::string m_PathInRootFile;

  const EventInfo* m_evtInfo;
  unsigned int m_fineTimeCut;
};

#endif
