// ********************************************************************
//
// NAME:     TrigT1CaloBSMon.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Peter Faulkner
//           
//
// ********************************************************************

#include <vector>

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/StatusCode.h"

#include "TrigT1CaloMonitoring/TrigT1CaloBSMon.h"
#include "TrigT1CaloMonitoringTools/TrigT1CaloMonErrorTool.h"
#include "TrigT1CaloCalibToolInterfaces/IL1CaloMonitoringCaloTool.h"

/*---------------------------------------------------------*/
TrigT1CaloBSMon::TrigT1CaloBSMon(const std::string & type, 
			 const std::string & name,
			 const IInterface* parent)
  : ManagedMonitorToolBase(type, name, parent),
    m_errorTool("TrigT1CaloMonErrorTool"),
    m_caloTool("LVL1::L1CaloMonitoringCaloTool/L1CaloMonitoringCaloTool")
/*---------------------------------------------------------*/
{
  declareProperty("LoadL1Calo", m_l1calo = false);
  declareProperty("LoadCaloCells", m_caloCells = false);
}

/*---------------------------------------------------------*/
TrigT1CaloBSMon::~TrigT1CaloBSMon()
/*---------------------------------------------------------*/
{
}

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

/*---------------------------------------------------------*/
StatusCode TrigT1CaloBSMon::initialize()
/*---------------------------------------------------------*/
{
  msg(MSG::INFO) << "Initializing " << name() << " - package version "
                 << PACKAGE_VERSION << endreq;

  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;

  if (m_l1calo) {
    sc = m_errorTool.retrieve();
    if( sc.isFailure() ) {
      msg(MSG::ERROR) << "Unable to locate Tool TrigT1CaloMonErrorTool"
                      << endreq;
      return sc;
    }
  }
  if (m_caloCells) {
    sc = m_caloTool.retrieve();
    if( sc.isFailure() ) {
      msg(MSG::ERROR) << "Unable to locate Tool L1CaloMonitoringCaloTool"
                      << endreq;
      return sc;
    }
  }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloBSMon::fillHistograms()
/*---------------------------------------------------------*/
{
  // Use error tool to force read of all L1Calo containers

  if (m_l1calo) {
    const std::vector<unsigned int>* errVecTES = 0;
    StatusCode sc = m_errorTool->retrieve(errVecTES);
    if (sc.isFailure()) return sc;
  }

  // Load CaloCells for faster data retrieval

  if (m_caloCells) {
    StatusCode sc = m_caloTool->loadCaloCells();
    if (sc.isFailure()) return sc;
  }

  return StatusCode::SUCCESS;
}
