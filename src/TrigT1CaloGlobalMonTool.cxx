// ********************************************************************
//
// NAME:     TrigT1CaloGlobalMonTool.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Peter Faulkner
//           
//
// ********************************************************************

#include <sstream>

#include "TAxis.h"
#include "TH2F.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/StatusCode.h"
#include "SGTools/StlVectorClids.h"

#include "TrigT1CaloMonitoring/TrigT1CaloGlobalMonTool.h"
#include "TrigT1CaloMonitoring/TrigT1CaloHistogramTool.h"

/*---------------------------------------------------------*/
TrigT1CaloGlobalMonTool::TrigT1CaloGlobalMonTool(const std::string & type, 
				                 const std::string & name,
				                 const IInterface* parent)
  : ManagedMonitorToolBase(type, name, parent),
    m_histTool("TrigT1CaloHistogramTool")
/*---------------------------------------------------------*/
{

  declareProperty("RootDirectory", m_rootDir = "L1Calo");

}

/*---------------------------------------------------------*/
TrigT1CaloGlobalMonTool::~TrigT1CaloGlobalMonTool()
/*---------------------------------------------------------*/
{
}

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

/*---------------------------------------------------------*/
StatusCode TrigT1CaloGlobalMonTool:: initialize()
/*---------------------------------------------------------*/
{
  msg(MSG::INFO) << "Initializing " << name() << " - package version "
                 << PACKAGE_VERSION << endreq;

  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;

  sc = m_histTool.retrieve();
  if( sc.isFailure() ) {
    msg(MSG::ERROR) << "Unable to locate Tool TrigT1CaloHistogramTool"
                    << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloGlobalMonTool::bookHistograms(bool isNewEventsBlock,
                                           bool isNewLumiBlock, bool isNewRun)
/*---------------------------------------------------------*/
{
  msg(MSG::DEBUG) << "bookHistograms entered" << endreq;

  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
  	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  if ( isNewEventsBlock || isNewLumiBlock ) { }

  if ( isNewRun ) {

  std::string dir(m_rootDir + "/Overview/Errors");
  MonGroup monGlobal( this, dir, shift, run );

  // Global Error Overview

  m_histTool->setMonGroup(&monGlobal);

  m_h_global = m_histTool->book2F("l1calo_2d_GlobalOverview",
                      "L1Calo Global Error Overview",
	              NumberOfGlobalErrors, 0, NumberOfGlobalErrors,
		      14, 0, 14);
  TAxis* axis = m_h_global->GetXaxis();
  axis->SetBinLabel(1+PPMDataStatus,   "PPMDataStatus");
  axis->SetBinLabel(1+PPMDataError,    "PPMDataError");
  axis->SetBinLabel(1+SubStatus,       "SubStatus");
  axis->SetBinLabel(1+Parity,          "Parity");
  axis->SetBinLabel(1+LinkDown,        "LinkDown");
  axis->SetBinLabel(1+RoIParity,       "RoIParity");
  axis->SetBinLabel(1+Transmission,    "Transmission");
  axis->SetBinLabel(1+Simulation,      "Simulation");
  axis->SetBinLabel(1+CMMSubStatus,    "CMMSubStatus");
  axis->SetBinLabel(1+GbCMMParity,     "CMMParity");
  axis->SetBinLabel(1+CMMTransmission, "CMMTransmission");
  axis->SetBinLabel(1+CMMSimulation,   "CMMSimulation");
  axis->SetBinLabel(1+RODStatus,       "RODStatus");
  axis->SetBinLabel(1+RODMissing,      "RODMissing");
  axis->SetBinLabel(1+ROBStatus,       "ROBStatus");
  axis->SetBinLabel(1+Unpacking,       "Unpacking");

  axis = m_h_global->GetYaxis();
  for (int crate = 0; crate < 14; ++crate) {
    int cr = crate;
    if (cr >= 12) cr -= 12;
    if (cr >= 8)  cr -= 8;
    std::string type = (crate < 8) ? "PP " : (crate < 12) ? "CP " : "JEP ";
    std::ostringstream cnum;
    cnum << type << cr;
    axis->SetBinLabel(crate+1, cnum.str().c_str());
  }

  m_histTool->unsetMonGroup();

  } // end if (isNewRun ...

  msg(MSG::DEBUG) << "Leaving bookHistograms" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloGlobalMonTool::fillHistograms()
/*---------------------------------------------------------*/
{
  msg(MSG::DEBUG) << "fillHistograms entered" << endreq;

  StatusCode sc;

  // Update Global overview plot

  const int ppmCrates = 8;
  const int cpmCrates = 4;
  const int jemCrates = 2;

  // PPM Error data
  const ErrorVector* errTES = 0; 
  if (evtStore()->contains<ErrorVector>("L1CaloPPMErrorVector")) {
    sc = evtStore()->retrieve(errTES, "L1CaloPPMErrorVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(ppmCrates)) {
    msg(MSG::DEBUG) << "No PPM error vector of expected size" << endreq;
  } else {
    for (int crate = 0; crate < ppmCrates; ++crate) {
      int err = (*errTES)[crate];
      if (err == 0) continue;
      if ((err >> DataStatus) & 0x1)   m_h_global->Fill(PPMDataStatus, crate);
      if ((err >> DataError) & 0x1)    m_h_global->Fill(PPMDataError,  crate);
      if ((err >> PPMSubStatus) & 0x1) m_h_global->Fill(SubStatus,     crate);
    }
  }

  // Spare PPM Channels Error data
  errTES = 0;
  if (evtStore()->contains<ErrorVector>("L1CaloPPMSpareErrorVector")) {
    sc = evtStore()->retrieve(errTES, "L1CaloPPMSpareErrorVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(ppmCrates)) {
    msg(MSG::DEBUG) << "No PPMSpare error vector of expected size" << endreq;
  } else {
    for (int crate = 0; crate < ppmCrates; ++crate) {
      int err = (*errTES)[crate];
      if (err == 0) continue;
      if ((err >> DataStatus) & 0x1)   m_h_global->Fill(PPMDataStatus, crate);
      if ((err >> DataError) & 0x1)    m_h_global->Fill(PPMDataError,  crate);
      if ((err >> PPMSubStatus) & 0x1) m_h_global->Fill(SubStatus,     crate);
    }
  }

  // CPM and CPM CMM Error data
  errTES = 0; 
  if (evtStore()->contains<ErrorVector>("L1CaloCPMErrorVector")) {
    sc = evtStore()->retrieve(errTES, "L1CaloCPMErrorVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(cpmCrates)) {
    msg(MSG::DEBUG) << "No CPM error vector of expected size" << endreq;
  } else {
    for (int crate = 0; crate < cpmCrates; ++crate) {
      int err = (*errTES)[crate];
      if (err == 0) continue;
      const int cr = crate + ppmCrates;
      if ((err >> CPMStatus) & 0x1) m_h_global->Fill(SubStatus, cr);
      if (((err >> CPMEMParity) & 0x1) || ((err >> CPMHadParity) & 0x1))
                                             m_h_global->Fill(Parity, cr);
      if (((err >> CPMEMLink) & 0x1) || ((err >> CPMHadLink) & 0x1))
                                             m_h_global->Fill(LinkDown, cr);
      if ((err >> CPMRoIParity) & 0x1) m_h_global->Fill(RoIParity, cr);
      if ((err >> CMMCPStatus) & 0x1)  m_h_global->Fill(CMMSubStatus, cr);
      if ((err >> CMMCPParity) & 0x1)  m_h_global->Fill(GbCMMParity, cr);
    }
  }

  // JEM Error data
  errTES = 0; 
  if (evtStore()->contains<ErrorVector>("L1CaloJEMErrorVector")) {
    sc = evtStore()->retrieve(errTES, "L1CaloJEMErrorVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(jemCrates)) {
    msg(MSG::DEBUG) << "No JEM error vector of expected size" << endreq;
  } else {
    for (int crate = 0; crate < jemCrates; ++crate) {
      int err = (*errTES)[crate];
      if (err == 0) continue;
      const int cr = crate + ppmCrates + cpmCrates;
      if ((err >> JEMStatus) & 0x1) m_h_global->Fill(SubStatus, cr);
      if (((err >> JEMEMParity) & 0x1) || ((err >> JEMHadParity) & 0x1))
                                             m_h_global->Fill(Parity, cr);
      if (((err >> JEMEMLink) & 0x1) || ((err >> JEMHadLink) & 0x1))
                                             m_h_global->Fill(LinkDown, cr);
      if ((err >> JEMRoIParity) & 0x1) m_h_global->Fill(RoIParity, cr);
    }
  }

  // JEM CMM Error data
  errTES = 0; 
  if (evtStore()->contains<ErrorVector>("L1CaloJEMCMMErrorVector")) {
    sc = evtStore()->retrieve(errTES, "L1CaloJEMCMMErrorVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(jemCrates)) {
    msg(MSG::DEBUG) << "No JEM CMM error vector of expected size" << endreq;
  } else {
    for (int crate = 0; crate < jemCrates; ++crate) {
      int err = (*errTES)[crate];
      if (err == 0) continue;
      const int cr = crate + ppmCrates + cpmCrates;
      if ((err >> JEMCMMStatus) & 0x1) m_h_global->Fill(CMMSubStatus, cr);
      if ((err >> JEMCMMParity) & 0x1) m_h_global->Fill(GbCMMParity, cr);
    }
  }

  // ROD Error data
  errTES = 0; 
  if (evtStore()->contains<ErrorVector>("L1CaloRODErrorVector")) {
    sc = evtStore()->retrieve(errTES, "L1CaloRODErrorVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(ppmCrates + cpmCrates + jemCrates)) {
    msg(MSG::DEBUG) << "No ROD error vector of expected size" << endreq;
  } else {
    for (int crate = 0; crate < ppmCrates+cpmCrates+jemCrates; ++crate) {
      int err = (*errTES)[crate];
      if (err == 0) continue;
      //if (err & 0x7f) m_h_global->Fill(RODStatus, crate);
      if (err & 0x3f) m_h_global->Fill(RODStatus, crate);
      if (((err >> NoFragment) & 0x1) || ((err >> NoPayload) & 0x1))
                      m_h_global->Fill(RODMissing, crate);
      if ((err >> ROBStatusError) & 0x1) m_h_global->Fill(ROBStatus, crate);
      if ((err >> UnpackingError) & 0x1) m_h_global->Fill(Unpacking, crate);
    }
  }

  // PPM Mismatch data
  errTES = 0; 
  if (evtStore()->contains<ErrorVector>("L1CaloPPMMismatchVector")) {
    sc = evtStore()->retrieve(errTES, "L1CaloPPMMismatchVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(ppmCrates)) {
    msg(MSG::DEBUG) << "No PPM mismatch vector of expected size" << endreq;
  } else {
    for (int crate = 0; crate < ppmCrates; ++crate) {
      int err = (*errTES)[crate];
      if (err == 0) continue;
      if (((err >> LUTMismatch) & 0x1)) m_h_global->Fill(Simulation, crate);
    }
  }

  // CPM Mismatch data
  errTES = 0; 
  if (evtStore()->contains<ErrorVector>("L1CaloCPMMismatchVector")) {
    sc = evtStore()->retrieve(errTES, "L1CaloCPMMismatchVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(cpmCrates)) {
    msg(MSG::DEBUG) << "No CPM mismatch vector of expected size" << endreq;
  } else {
    for (int crate = 0; crate < cpmCrates; ++crate) {
      int err = (*errTES)[crate];
      if (err == 0) continue;
      const int cr = crate + ppmCrates;
      if (((err >> EMTowerMismatch) & 0x1) || ((err >> HadTowerMismatch) & 0x1))
                                        m_h_global->Fill(Transmission, cr);
      if (((err >> CPMRoIMismatch) & 0x1) || ((err >> CPMHitsMismatch) & 0x1))
                                        m_h_global->Fill(Simulation, cr);
      if (((err >> CMMHitsMismatch) & 0x1) || ((err >> RemoteSumMismatch) & 0x1))
                                        m_h_global->Fill(CMMTransmission, cr);
      if (((err >> LocalSumMismatch) & 0x1) || ((err >> TotalSumMismatch) & 0x1))
                                        m_h_global->Fill(CMMSimulation, cr);
    }
  }

  // JEM Mismatch data
  errTES = 0; 
  if (evtStore()->contains<ErrorVector>("L1CaloJEMMismatchVector")) {
    sc = evtStore()->retrieve(errTES, "L1CaloJEMMismatchVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(jemCrates)) {
    msg(MSG::DEBUG) << "No JEM mismatch vector of expected size" << endreq;
  } else {
    for (int crate = 0; crate < jemCrates; ++crate) {
      int err = (*errTES)[crate];
      if (err == 0) continue;
      const int cr = crate + ppmCrates + cpmCrates;
      if (((err >> EMElementMismatch) & 0x1)  ||
          ((err >> HadElementMismatch) & 0x1) ||
          ((err >> JEMRoIMismatch) & 0x1)     ||
          ((err >> JEMHitsMismatch) & 0x1)    ||
          ((err >> JEMEtSumsMismatch) & 0x1)) m_h_global->Fill(Simulation, cr);
      if (((err >> CMMJetHitsMismatch) & 0x1)   ||
          ((err >> RemoteJetMismatch) & 0x1)    ||
	  ((err >> JetEtRoIMismatch) & 0x1)     ||
	  ((err >> CMMEtSumsMismatch) & 0x1)    ||
	  ((err >> RemoteEnergyMismatch) & 0x1) ||
	  ((err >> EnergyRoIMismatch) & 0x1))
	                              m_h_global->Fill(CMMTransmission, cr);
      if (((err >> LocalJetMismatch) & 0x1)    ||
          ((err >> TotalJetMismatch) & 0x1)    ||
	  ((err >> JetEtMismatch) & 0x1)       ||
          ((err >> LocalEnergyMismatch) & 0x1) ||
          ((err >> TotalEnergyMismatch) & 0x1) ||
	  ((err >> SumEtMismatch) & 0x1)       ||
	  ((err >> MissingEtMismatch) & 0x1))
	                                m_h_global->Fill(CMMSimulation, cr);
    }
  }

  msg(MSG::DEBUG) << "Leaving fillHistograms" << endreq;

  return StatusCode::SUCCESS;

}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloGlobalMonTool::procHistograms(bool isEndOfEventsBlock,
                                  bool isEndOfLumiBlock, bool isEndOfRun)
/*---------------------------------------------------------*/
{
  msg(MSG::DEBUG) << "procHistograms entered" << endreq;

  if (isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun) {
  }

  return StatusCode::SUCCESS;
}
