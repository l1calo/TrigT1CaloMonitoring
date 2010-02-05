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
#include <vector>

#include "TH2F.h"

#include "GaudiKernel/ITHistSvc.h"
#include "StoreGate/StoreGateSvc.h"
#include "SGTools/StlVectorClids.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigT1CaloMonitoring/TrigT1CaloGlobalMonTool.h"

/*---------------------------------------------------------*/
TrigT1CaloGlobalMonTool::TrigT1CaloGlobalMonTool(const std::string & type, 
				                 const std::string & name,
				                 const IInterface* parent)
  : ManagedMonitorToolBase(type, name, parent),
    m_storeGate("StoreGateSvc", name),
    m_log(msgSvc(), name), m_monGroup(0)
/*---------------------------------------------------------*/
{
  declareInterface<IMonitorToolBase>(this); 

  declareProperty("RootDirectory", m_rootDir = "L1Calo");

}

/*---------------------------------------------------------*/
TrigT1CaloGlobalMonTool::~TrigT1CaloGlobalMonTool()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloGlobalMonTool:: initialize()
/*---------------------------------------------------------*/
{
  m_log.setLevel(outputLevel());

  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;
  
  sc = m_storeGate.retrieve();
  if( sc.isFailure() ) {
    m_log << MSG::ERROR << "Unable to locate Service StoreGateSvc" << endreq;
    return sc;
  }

  m_log << MSG::INFO << "TrigT1CaloGlobalMonTool initialised" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloGlobalMonTool::bookHistograms(bool isNewEventsBlock,
                                           bool isNewLumiBlock, bool isNewRun)
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "bookHistograms entered" << endreq;

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

  m_monGroup = &monGlobal;

  m_h_global = book2F("l1calo_2d_GlobalOverview",
                      "L1Calo Global Error Overview;;Crate",
	              NumberOfGlobalErrors, 0, NumberOfGlobalErrors,
		      15, 0, 15);
  m_h_global->GetXaxis()->SetBinLabel(1+PPMDataStatus,   "PPMDataStatus");
  m_h_global->GetXaxis()->SetBinLabel(1+PPMDataError,    "PPMDataError");
  m_h_global->GetXaxis()->SetBinLabel(1+SubStatus,       "SubStatus");
  m_h_global->GetXaxis()->SetBinLabel(1+Parity,          "Parity");
  m_h_global->GetXaxis()->SetBinLabel(1+LinkDown,        "LinkDown");
  m_h_global->GetXaxis()->SetBinLabel(1+RoIParity,       "RoIParity");
  m_h_global->GetXaxis()->SetBinLabel(1+Transmission,    "Transmission");
  m_h_global->GetXaxis()->SetBinLabel(1+Simulation,      "Simulation");
  m_h_global->GetXaxis()->SetBinLabel(1+CMMSubStatus,    "CMMSubStatus");
  m_h_global->GetXaxis()->SetBinLabel(1+GbCMMParity,     "CMMParity");
  m_h_global->GetXaxis()->SetBinLabel(1+CMMTransmission, "CMMTransmission");
  m_h_global->GetXaxis()->SetBinLabel(1+CMMSimulation,   "CMMSimulation");
  m_h_global->GetXaxis()->SetBinLabel(1+RODStatus,       "RODStatus");
  m_h_global->GetXaxis()->SetBinLabel(1+RODMissing,      "RODMissing");
  m_h_global->GetXaxis()->SetBinLabel(1+ROBStatus,       "ROBStatus");
  m_h_global->GetXaxis()->SetBinLabel(1+Unpacking,       "Unpacking");

  for (int crate = 0; crate < 14; ++crate) {
    int cr = crate;
    if (cr >= 12) cr -= 12;
    if (cr >= 8)  cr -= 8;
    std::ostringstream cnum;
    cnum << cr;
    m_h_global->GetYaxis()->SetBinLabel(crate+1, cnum.str().c_str());
  }
  m_h_global->GetYaxis()->SetBinLabel(1, "PP 0");
  m_h_global->GetYaxis()->SetBinLabel(9, "CP 0");
  m_h_global->GetYaxis()->SetBinLabel(13, "JEP 0");

  } // end if (isNewRun ...

  m_log << MSG::DEBUG << "Leaving bookHistograms" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloGlobalMonTool::fillHistograms()
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "fillHistograms entered" << endreq;

  StatusCode sc;

  // Update Global overview plot

  const int ppmCrates = 8;
  const int cpmCrates = 4;
  const int jemCrates = 2;

  // PPM Error data
  const ErrorVector* errTES = 0; 
  if (m_storeGate->contains<ErrorVector>("L1CaloPPMErrorVector")) {
    sc = m_storeGate->retrieve(errTES, "L1CaloPPMErrorVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(ppmCrates)) {
    m_log << MSG::DEBUG << "No PPM error vector of expected size" << endreq;
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
  if (m_storeGate->contains<ErrorVector>("L1CaloCPMErrorVector")) {
    sc = m_storeGate->retrieve(errTES, "L1CaloCPMErrorVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(cpmCrates)) {
    m_log << MSG::DEBUG << "No CPM error vector of expected size" << endreq;
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
  if (m_storeGate->contains<ErrorVector>("L1CaloJEMErrorVector")) {
    sc = m_storeGate->retrieve(errTES, "L1CaloJEMErrorVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(jemCrates)) {
    m_log << MSG::DEBUG << "No JEM error vector of expected size" << endreq;
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
  if (m_storeGate->contains<ErrorVector>("L1CaloJEMCMMErrorVector")) {
    sc = m_storeGate->retrieve(errTES, "L1CaloJEMCMMErrorVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(jemCrates)) {
    m_log << MSG::DEBUG << "No JEM CMM error vector of expected size" << endreq;
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
  if (m_storeGate->contains<ErrorVector>("L1CaloRODErrorVector")) {
    sc = m_storeGate->retrieve(errTES, "L1CaloRODErrorVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(ppmCrates + cpmCrates + jemCrates)) {
    m_log << MSG::DEBUG << "No ROD error vector of expected size" << endreq;
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

  // CPM Mismatch data
  errTES = 0; 
  if (m_storeGate->contains<ErrorVector>("L1CaloCPMMismatchVector")) {
    sc = m_storeGate->retrieve(errTES, "L1CaloCPMMismatchVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(cpmCrates)) {
    m_log << MSG::DEBUG << "No CPM mismatch vector of expected size" << endreq;
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
  if (m_storeGate->contains<ErrorVector>("L1CaloJEMMismatchVector")) {
    sc = m_storeGate->retrieve(errTES, "L1CaloJEMMismatchVector"); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure() || errTES->size() != size_t(jemCrates)) {
    m_log << MSG::DEBUG << "No JEM mismatch vector of expected size" << endreq;
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

  m_log << MSG::DEBUG << "Leaving fillHistograms" << endreq;

  return StatusCode::SUCCESS;

}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloGlobalMonTool::procHistograms(bool isEndOfEventsBlock,
                                  bool isEndOfLumiBlock, bool isEndOfRun)
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "procHistograms entered" << endreq;

  if (isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun) {
  }

  return StatusCode::SUCCESS;
}

TH2F* TrigT1CaloGlobalMonTool::book2F(const std::string& name,
                                   const std::string& title,
                                   int nx, double xmin, double xmax,  
	                           int ny, double ymin, double ymax)
{		
  TH2F *hist = new TH2F(name.c_str(), title.c_str(), nx, xmin, xmax,
                                                     ny, ymin, ymax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    m_log << MSG::WARNING << "Could not register histogram : " 
	  << name << endreq;
  }
  hist->SetOption("colz");
  hist->SetStats(kFALSE);
  
  return hist;
}
