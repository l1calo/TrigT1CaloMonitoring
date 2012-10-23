// ********************************************************************
//
// NAME:     TrigT1CaloBSMon.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Peter Faulkner
//	     
//           Force read of all L1Calo containers for cpu accounting
//
// ********************************************************************
#ifndef TRIGT1CALOBSMON_H
#define TRIGT1CALOBSMON_H

#include "GaudiKernel/ToolHandle.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"

class StatusCode;
class TrigT1CaloMonErrorTool;

namespace LVL1 {
  class IL1CaloMonitoringCaloTool;
}

class TrigT1CaloBSMon: public ManagedMonitorToolBase
{

public:
  
  TrigT1CaloBSMon(const std::string & type, const std::string & name,
		                            const IInterface* parent);
    

  virtual ~TrigT1CaloBSMon();

  virtual StatusCode initialize();
  virtual StatusCode fillHistograms();

private:

  ToolHandle<TrigT1CaloMonErrorTool>          m_errorTool;
  ToolHandle<LVL1::IL1CaloMonitoringCaloTool> m_caloTool;

  bool m_l1calo;
  bool m_caloCells;

};

#endif
