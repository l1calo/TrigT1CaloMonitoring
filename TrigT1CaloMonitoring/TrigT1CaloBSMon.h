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

/** Force load of all L1Calo containers and/or CaloCell container.
 *
 *  Allows separate accounting of cpu time for data unpacking.
 *  For testing purposes only, not used in normal running.
 *
 *  <b>Tools Used:</b>
 *
 *  <table>
 *  <tr><th> Tool                               </th><th> Description          </th></tr>
 *  <tr><td> @c TrigT1CaloMonErrorTool          </td><td> @copydoc m_errorTool </td></tr>
 *  <tr><td> @c LVL1::IL1CaloMonitoringCaloTool </td><td> @copydoc m_caloTool  </td></tr>
 *  </table>
 *
 *  <b>JobOption Properties:</b>
 *
 *  <table>
 *  <tr><th> Property         </th><th> Description          </th></tr>
 *  <tr><td> @c LoadL1Calo    </td><td> @copydoc m_l1calo    </td></tr>
 *  <tr><td> @c LoadCaloCells </td><td> @copydoc m_caloCells </td></tr>
 *  </table>
 *
 *  @author Peter Faulkner
 *
 */

class TrigT1CaloBSMon: public ManagedMonitorToolBase
{

public:
  
  TrigT1CaloBSMon(const std::string & type, const std::string & name,
		                            const IInterface* parent);
    

  virtual ~TrigT1CaloBSMon();

  virtual StatusCode initialize();
  virtual StatusCode fillHistograms();

private:

  /// ByteStream unpacking error tool (forces L1Calo data read)
  ToolHandle<TrigT1CaloMonErrorTool>          m_errorTool;
  /// CaloCell info by TT tool (forces CaloCell read)
  ToolHandle<LVL1::IL1CaloMonitoringCaloTool> m_caloTool;

  /// Switch for L1Calo
  bool m_l1calo;
  /// Switch for CaloCells
  bool m_caloCells;

};

#endif
