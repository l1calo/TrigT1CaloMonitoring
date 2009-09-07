// ********************************************************************
//
// NAME:     TrigT1CaloRodMonTool.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Peter Faulkner
//           
//
// ********************************************************************

#include <cmath>
#include <sstream>

#include "TH1F.h"
#include "TH2F.h"

#include "GaudiKernel/ITHistSvc.h"
#include "StoreGate/StoreGateSvc.h"
#include "SGTools/StlVectorClids.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigT1CaloEvent/RODHeader.h"

#include "TrigT1CaloMonitoring/TrigT1CaloRodMonTool.h"

/*---------------------------------------------------------*/
TrigT1CaloRodMonTool::TrigT1CaloRodMonTool(const std::string & type, 
				           const std::string & name,
				           const IInterface* parent)
  : ManagedMonitorToolBase(type, name, parent),
    m_storeGate("StoreGateSvc", name),
    m_log(msgSvc(), name), m_monGroup(0)
/*---------------------------------------------------------*/
{
  declareInterface<IMonitorToolBase>(this); 

  declareProperty("RodHeaderLocation", m_rodHeaderLocation = "RODHeaders");
  m_cpRoibRodHeaderLocation  = m_rodHeaderLocation + "CPRoIB";
  m_jepRoibRodHeaderLocation = m_rodHeaderLocation + "JEPRoIB";
  declareProperty("L1CaloErrorLocation",
                   m_robErrorVectorLocation = "L1CaloUnpackingErrors");

  declareProperty("RootDirectory", m_rootDir = "L1Calo");
  declareProperty("OnlineTest", m_onlineTest = false,
                  "Test online code when running offline");

}

/*---------------------------------------------------------*/
TrigT1CaloRodMonTool::~TrigT1CaloRodMonTool()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloRodMonTool:: initialize()
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

  m_log << MSG::INFO << "TrigT1CaloRodMonTool initialised" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloRodMonTool::bookHistograms(bool isNewEventsBlock,
                                           bool isNewLumiBlock, bool isNewRun)
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "bookHistograms entered" << endreq;

  bool online = m_onlineTest;
  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
    online = true;
  }
  	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  if ( isNewEventsBlock || isNewLumiBlock ) { }

  if ( isNewRun ) {

  std::string dir1(m_rootDir + "/ROD");
  MonGroup monShift ( this, dir1, shift, run );
  MonGroup monExpert( this, dir1, expert, run );
  MonGroup monAverage( this, dir1, expert, run, "", "weightedAverage" );
  MonGroup monROB( this, dir1 + "/ROBStatus", expert, run );
  MonGroup monUnpack( this, dir1 + "/Unpacking", expert, run );
  std::string dir2(m_rootDir + "/Overview/Errors");
  MonGroup monGlobal( this, dir2, shift, run );

  //  Payload Averages

  m_monGroup = &monAverage;

  m_sumPayloads1.clear();
  m_sumPayloads1.resize(80, 0.);
  m_sumPayloads2.clear();
  m_sumPayloads2.resize(160, 0.);
  m_events = 0;
  std::string axisTitles = (online)
         ? ";Complete Run | Recent Events        Crate/S-Link;Words per Event"
	 : ";Crate/S-Link;Words per Event";
  int nbins = (online) ? 65 : 32;
  m_h_ROD_PP = book1F("rod_1d_PP_payload",
                      "ROD PP Average Payload Size"+axisTitles,
                       nbins, 0, nbins);
  setLabelsCSL(m_h_ROD_PP, true, 1, 32, 2, 2);
  if (online) {
    m_h_ROD_PP->GetXaxis()->SetBinLabel(33, "---");
    setLabelsCSL(m_h_ROD_PP, true, 34, 65, 2, 2);
  }
  nbins = (online) ? 17 : 8;
  m_h_ROD_CP = book1F("rod_1d_CP_payload",
                      "ROD CP Average Payload Size"+axisTitles,
                       nbins, 0, nbins);
  setLabelsCSL(m_h_ROD_CP, true, 1, 8, 1, 2);
  if (online) {
    m_h_ROD_CP->GetXaxis()->SetBinLabel(9, "---");
    setLabelsCSL(m_h_ROD_CP, true, 10, 17, 1, 2);
  }
  m_h_ROD_JEP = book1F("rod_1d_JEP_payload",
                       "ROD JEP Average Payload Size"+axisTitles,
                        nbins, 0, nbins);
  setLabelsCSL(m_h_ROD_JEP, true, 1, 8, 1, 1);
  if (online) {
    m_h_ROD_JEP->GetXaxis()->SetBinLabel(9, "---");
    setLabelsCSL(m_h_ROD_JEP, true, 10, 17, 1, 1);
  }
  nbins = (online) ? 25 : 12;
  m_h_ROD_RoI = book1F("rod_1d_CPJEP_RoI_payload",
                       "ROD CP and JEP RoI Average Payload Size"+axisTitles,
                                nbins, 0, nbins);
  setLabelsCSL(m_h_ROD_RoI, true, 1, 8, 1, 2);
  m_h_ROD_RoI->GetXaxis()->SetBinLabel(1, "CP");
  setLabelsCSL(m_h_ROD_RoI, true, 9, 12, 1, 2);
  m_h_ROD_RoI->GetXaxis()->SetBinLabel(9, "JEP");
  if (online) {
    m_h_ROD_RoI->GetXaxis()->SetBinLabel(13, "---");
    setLabelsCSL(m_h_ROD_RoI, true, 14, 21, 1, 2);
    m_h_ROD_RoI->GetXaxis()->SetBinLabel(14, "CP");
    setLabelsCSL(m_h_ROD_RoI, true, 22, 25, 1, 2);
    m_h_ROD_RoI->GetXaxis()->SetBinLabel(22, "JEP");
  }

  //  Status bits

  m_monGroup = &monExpert;

  m_h_ROD_PP_stat = book2F("rod_2d_PP_status",
                           "PP ROD Status Bits and Payload Check;;Crate/S-Link",
                           NumberOfStatusBins, 0, NumberOfStatusBins,
			   32, 0, 32);
  setLabelsStatus(m_h_ROD_PP_stat);
  setLabelsCSL(m_h_ROD_PP_stat, false, 1, 32, 2, 2);
  m_h_ROD_PP_stat->GetXaxis()->SetBinLabel(1+NoPayload, "No Payload");
  m_h_ROD_CPJEP_stat = book2F("rod_2d_CPJEP_status",
                   "CP and JEP ROD Status Bits and Payload Check;;Crate/S-Link",
                              NumberOfStatusBins, 0, NumberOfStatusBins,
			      16, 0, 16);
  setLabelsStatus(m_h_ROD_CPJEP_stat);
  setLabelsCSL(m_h_ROD_CPJEP_stat, false, 1, 8, 1, 2);
  m_h_ROD_CPJEP_stat->GetYaxis()->SetBinLabel(1, "CP 0/0");
  setLabelsCSL(m_h_ROD_CPJEP_stat, false, 9, 16, 1, 1);
  m_h_ROD_CPJEP_stat->GetYaxis()->SetBinLabel(9, "JEP 0/0");
  m_h_ROD_CPJEP_stat->GetXaxis()->SetBinLabel(1+NoPayload, "No Payload");
  m_h_ROD_RoI_stat = book2F("rod_2d_CPJEP_RoI_status",
                            "CP and JEP RoI ROD Status Bits;;Crate/S-Link",
                            NumberOfStatusBins, 0, NumberOfStatusBins,
			    12, 0, 12);
  setLabelsStatus(m_h_ROD_RoI_stat);
  setLabelsCSL(m_h_ROD_RoI_stat, false, 1, 8, 1, 2);
  m_h_ROD_RoI_stat->GetYaxis()->SetBinLabel(1, "CP 0/0");
  setLabelsCSL(m_h_ROD_RoI_stat, false, 9, 12, 1, 2);
  m_h_ROD_RoI_stat->GetYaxis()->SetBinLabel(9, "JEP 0/0");

  //  ROB Status bits

  m_monGroup = &monROB;

  m_h_ROD_PP_rob = book2F("rod_2d_PP_ROB_status",
                          "PP ROB Status Bits;;Crate/S-Link",
                           19, 0, 19, 32, 0, 32);
  setLabelsROBStatus(m_h_ROD_PP_rob);
  setLabelsCSL(m_h_ROD_PP_rob, false, 1, 32, 2, 2);
  m_h_ROD_CPJEP_rob = book2F("rod_2d_CPJEP_ROB_status",
                             "CP and JEP ROB Status Bits;;Crate/S-Link",
                              19, 0, 19, 16, 0, 16);
  setLabelsROBStatus(m_h_ROD_CPJEP_rob);
  setLabelsCSL(m_h_ROD_CPJEP_rob, false, 1, 8, 1, 2);
  m_h_ROD_CPJEP_rob->GetYaxis()->SetBinLabel(1, "CP 0/0");
  setLabelsCSL(m_h_ROD_CPJEP_rob, false, 9, 16, 1, 1);
  m_h_ROD_CPJEP_rob->GetYaxis()->SetBinLabel(9, "JEP 0/0");
  m_h_ROD_RoI_rob = book2F("rod_2d_CPJEP_RoI_ROB_status",
                           "CP and JEP RoI ROB Status Bits;;Crate/S-Link",
                            19, 0, 19, 12, 0, 12);
  setLabelsROBStatus(m_h_ROD_RoI_rob);
  setLabelsCSL(m_h_ROD_RoI_rob, false, 1, 8, 1, 2);
  m_h_ROD_RoI_rob->GetYaxis()->SetBinLabel(1, "CP 0/0");
  setLabelsCSL(m_h_ROD_RoI_rob, false, 9, 12, 1, 2);
  m_h_ROD_RoI_rob->GetYaxis()->SetBinLabel(9, "JEP 0/0");

  //  Unpacking Errors

  m_monGroup = &monUnpack;

  m_h_ROD_PP_unp = book2F("rod_2d_PP_unpack",
                          "PP Bytestream Unpacking Errors;;Crate/S-Link",
                           18, 1, 19, 32, 0, 32);
  setLabelsUnpacking(m_h_ROD_PP_unp);
  setLabelsCSL(m_h_ROD_PP_unp, false, 1, 32, 2, 2);
  m_h_ROD_CPJEP_unp = book2F("rod_2d_CPJEP_unpack",
                       "CP and JEP Bytestream Unpacking Errors;;Crate/S-Link",
                        18, 1, 19, 16, 0, 16);
  setLabelsUnpacking(m_h_ROD_CPJEP_unp);
  setLabelsCSL(m_h_ROD_CPJEP_unp, false, 1, 8, 1, 2);
  m_h_ROD_CPJEP_unp->GetYaxis()->SetBinLabel(1, "CP 0/0");
  setLabelsCSL(m_h_ROD_CPJEP_unp, false, 9, 16, 1, 1);
  m_h_ROD_CPJEP_unp->GetYaxis()->SetBinLabel(9, "JEP 0/0");
  m_h_ROD_RoI_unp = book2F("rod_2d_CPJEP_RoI_unpack",
                   "CP and JEP RoI Bytestream Unpacking Errors;;Crate/S-Link",
                    18, 1, 19, 12, 0, 12);
  setLabelsUnpacking(m_h_ROD_RoI_unp);
  setLabelsCSL(m_h_ROD_RoI_unp, false, 1, 8, 1, 2);
  m_h_ROD_RoI_unp->GetYaxis()->SetBinLabel(1, "CP 0/0");
  setLabelsCSL(m_h_ROD_RoI_unp, false, 9, 12, 1, 2);
  m_h_ROD_RoI_unp->GetYaxis()->SetBinLabel(9, "JEP 0/0");

  //  Error summary

  m_monGroup = &monShift;

  m_h_ROD_summary = book1F("rod_1d_Error_summary", "ROD Error Summary;;Events",
                           NumberOfStatusBins, 0, NumberOfStatusBins);
  setLabelsStatus(m_h_ROD_summary);
  m_h_ROD_summary->GetXaxis()->SetBinLabel(1+TriggerType, "");
  m_h_ROD_summary->GetXaxis()->SetBinLabel(1+LimitedRoI,  "");
  m_h_ROD_summary->GetXaxis()->SetBinLabel(1+NoPayload, "No Payload");
  m_h_ROB_summary = book1F("rod_1d_ROB_Error_summary",
                           "ROB Status Error Summary;;Events", 19, 0, 19);
  setLabelsROBStatus(m_h_ROB_summary);
  m_h_Unp_summary = book1F("rod_1d_Unpack_Error_summary",
                           "Bytestream Unpacking Error Summary;;Events",
                           18, 1, 19);
  setLabelsUnpacking(m_h_Unp_summary);

  // Global Error Overview

  m_monGroup = &monGlobal;

  m_h_global = book2F("l1calo_2d_Global_overview",
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
StatusCode TrigT1CaloRodMonTool::fillHistograms()
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "fillHistograms entered" << endreq;

  const bool online = m_onlineTest ||
                     (m_environment == AthenaMonManager::online);

  //Retrieve DAQ ROD Headers from SG
  StatusCode sc;
  const RodHeaderCollection* rodTES = 0; 
  if (m_storeGate->contains<RodHeaderCollection>(m_rodHeaderLocation)) {
    sc = m_storeGate->retrieve(rodTES, m_rodHeaderLocation); 
  } else sc = StatusCode::FAILURE;
  if( sc.isFailure()  ||  !rodTES ) {
    m_log << MSG::DEBUG<< "No DAQ ROD Header container found"<< endreq; 
  }

  //Retrieve CP RoIB ROD Headers from SG
  const RodHeaderCollection* cpRoibTES = 0; 
  if (m_storeGate->contains<RodHeaderCollection>(m_cpRoibRodHeaderLocation)) {
    sc = m_storeGate->retrieve(cpRoibTES, m_cpRoibRodHeaderLocation); 
  } else sc = StatusCode::FAILURE;
  if( sc.isFailure()  ||  !cpRoibTES ) {
    m_log << MSG::DEBUG<< "No CP RoIB ROD Header container found"<< endreq; 
  }

  //Retrieve JEP RoIB ROD Headers from SG
  const RodHeaderCollection* jepRoibTES = 0; 
  if (m_storeGate->contains<RodHeaderCollection>(m_jepRoibRodHeaderLocation)) {
    sc = m_storeGate->retrieve(jepRoibTES, m_jepRoibRodHeaderLocation); 
  } else sc = StatusCode::FAILURE;
  if( sc.isFailure()  ||  !jepRoibTES ) {
    m_log << MSG::DEBUG<< "No JEP RoIB ROD Header container found"<< endreq; 
  }

  //Retrieve ROB and Unpacking Error vector from SG
  const ROBErrorCollection* errVecTES = 0;
  if (m_storeGate->contains<ROBErrorCollection>(m_robErrorVectorLocation)) {
    sc = m_storeGate->retrieve(errVecTES, m_robErrorVectorLocation);
  } else sc = StatusCode::FAILURE;
  if( sc.isFailure()  ||  !errVecTES ) {
    m_log << MSG::DEBUG<< "No ROB and Unpacking Error vector found"<< endreq;
  }

  //if ( !rodTES && !cpRoibTES && !jepRoibTES ) {
  //  m_log << MSG::DEBUG<< "No ROD Header containers found"<< endreq;
  //  return StatusCode::SUCCESS;
  //}

  //=============================================
  //   ROD Payload plots
  //=============================================

  std::vector<int> errors(NumberOfStatusBins);
  std::vector<int> errorsROB(19);
  std::vector<int> errorsUnpack(19);
  std::vector<int> crateErr(14);
  std::vector<int> noFragmentFlags(80, 1);
  std::vector<int> noPayloadFlags(56, 1);
  std::vector<const RodHeaderCollection*> cols;
  if (rodTES)     cols.push_back(rodTES);
  if (cpRoibTES)  cols.push_back(cpRoibTES);
  if (jepRoibTES) cols.push_back(jepRoibTES);
  std::vector<const RodHeaderCollection*>::const_iterator colIter =
                                                                 cols.begin();
  std::vector<const RodHeaderCollection*>::const_iterator colIterEnd =
                                                                   cols.end();
  for (; colIter != colIterEnd; ++colIter) {
    RodHeaderCollection::const_iterator iter    = (*colIter)->begin();
    RodHeaderCollection::const_iterator iterEnd = (*colIter)->end();
    for (; iter != iterEnd; ++iter) {
      const LVL1::RODHeader* header = *iter;
      const int crate = header->crate();
      const int slink = header->sLink();
      const int dataType = header->dataType();
      const int rod = crate + dataType*6;
      const int nData = header->payloadSize();
      const int pos = rod*4 + slink;
      // Skip obviously corrupt data
      if (crate > 13 || slink > 3 || nData < 0 ||
                        nData > 10000 || pos >= 80) continue;
      noFragmentFlags[pos] = 0;
      if (pos < 56 && nData > 0) noPayloadFlags[pos] = 0;
      m_sumPayloads1[pos] += nData;
      m_sumPayloads2[pos] += nData;
      // Status bits
      TH2F* hist = m_h_ROD_PP_stat;
      int val = pos;
      if (pos >= 72) {
        hist = m_h_ROD_RoI_stat;
	val = (pos-72)/2 + 8;
      } else if (pos >= 56) {
        hist = m_h_ROD_RoI_stat;
	val = (pos-56)/2;
      } else if (pos >= 48) {
        hist = m_h_ROD_CPJEP_stat;
	val = pos-48 + 8;
      } else if (pos >= 32) {
        hist = m_h_ROD_CPJEP_stat;
	val = (pos-32)/2;
      }
      // gLinkError is actually OR'ed with cmmParityError
      //if (header->gLinkError()) {
      if (header->cmmParityError()) {
        hist->Fill(GLink, val);
	errors[GLink] = 1;
	crateErr[crate] |= (1 << GLink);
      }
      //if (header->cmmParityError()) {
      //  hist->Fill(CMMParity, val);
      //  errors[CMMParity] = 1;
      //  crateErr[crate] |= (1 << CMMParity);
      //}
      if (header->lvdsLinkError()) {
        hist->Fill(LVDSLink, val);
	errors[LVDSLink] = 1;
	crateErr[crate] |= (1 << LVDSLink);
      }
      if (header->rodFifoOverflow()) {
        hist->Fill(FIFOOverflow, val);
	errors[FIFOOverflow] = 1;
	crateErr[crate] |= (1 << FIFOOverflow);
      }
      if (header->dataTransportError()) {
        hist->Fill(DataTransport, val);
	errors[DataTransport] = 1;
	crateErr[crate] |= (1 << DataTransport);
      }
      if (header->gLinkTimeout()) {
        hist->Fill(Timeout, val);
	errors[Timeout] = 1;
	crateErr[crate] |= (1 << Timeout);
      }
      if (header->bcnMismatch()) {
        hist->Fill(BCNMismatch, val);
	errors[BCNMismatch] = 1;
	crateErr[crate] |= (1 << BCNMismatch);
      }
      if (header->triggerTypeTimeout()) hist->Fill(TriggerType, val);
      if (pos >= 56 && header->limitedRoISet()) hist->Fill(LimitedRoI, val);
    }
  }

  // Update average payload plots

  ++m_events;
  const int events1 = 20;
  int events2 = m_events % events1;
  if (events2 == 0) events2 = events1;
  int events3 = events1;
  if (m_events <= events1) events3 = 0;
  const int events4 = events2 + events3;
  const double error1 = 1./sqrt(m_events);
  const double error2 = 1./sqrt(events4);
  for (int i = 0; i < 80; ++i) {
    const double average1 = m_sumPayloads1[i]/m_events;
    const double average2 = (m_sumPayloads2[i]+m_sumPayloads2[i+80])/events4;
    if (events2 == events1) {
      m_sumPayloads2[i+80] = m_sumPayloads2[i];
      m_sumPayloads2[i] = 0;
    }
    if (i >= 72) {
      if (i%2 == 0) {
        int bin = (i-72)/2 + 9;
        m_h_ROD_RoI->SetBinContent(bin, average1);
	m_h_ROD_RoI->SetBinError(bin, error1);
	if (online) {
          m_h_ROD_RoI->SetBinContent(13+bin, average2);
	  m_h_ROD_RoI->SetBinError(13+bin, error2);
        }
      }
    } else if (i >= 56) {
      if (i%2 == 0) {
        int bin = (i-56)/2 + 1;
        m_h_ROD_RoI->SetBinContent(bin, average1);
	m_h_ROD_RoI->SetBinError(bin, error1);
	if (online) {
          m_h_ROD_RoI->SetBinContent(13+bin, average2);
	  m_h_ROD_RoI->SetBinError(13+bin, error2);
        }
      }
    } else if (i >= 48) {
      int bin = i-48 + 1;
      m_h_ROD_JEP->SetBinContent(bin, average1);
      m_h_ROD_JEP->SetBinError(bin, error1);
      if (online) {
        m_h_ROD_JEP->SetBinContent(9+bin, average2);
        m_h_ROD_JEP->SetBinError(9+bin, error2);
      }
    } else if (i >= 32) {
      if (i%2 == 0) {
        int bin = (i-32)/2 + 1;
        m_h_ROD_CP->SetBinContent(bin, average1);
	m_h_ROD_CP->SetBinError(bin, error1);
	if (online) {
          m_h_ROD_CP->SetBinContent(9+bin, average2);
	  m_h_ROD_CP->SetBinError(9+bin, error2);
        }
      }
    } else {
      int bin = i + 1;
      m_h_ROD_PP->SetBinContent(bin, average1);
      m_h_ROD_PP->SetBinError(bin, error1);
      if (online) {
        m_h_ROD_PP->SetBinContent(33+bin, average2);
        m_h_ROD_PP->SetBinError(33+bin, error2);
      }
    }
  }

  // Update missing ROD fragments and payloads

  for (int pos = 0; pos < 80; ++pos) {
    if (noFragmentFlags[pos] || (pos < 56 && noPayloadFlags[pos])) {
      TH2F* hist = m_h_ROD_PP_stat;
      int val = pos;
      int crate = pos/4;
      if (crate > 13) crate -= 6;
      if (pos >= 72) {
        if (pos%2) continue;
	hist = m_h_ROD_RoI_stat;
	val = (pos-72)/2 + 8;
      } else if (pos >= 56) {
        if (pos%2) continue;
        hist = m_h_ROD_RoI_stat;
	val = (pos-56)/2;
      } else if (pos >= 48) {
        hist = m_h_ROD_CPJEP_stat;
        val = pos-48 + 8;
      } else if (pos >= 32) {
	if (pos%2) continue;
        hist = m_h_ROD_CPJEP_stat;
        val = (pos-32)/2;
      }
      if (noFragmentFlags[pos]) {
        hist->Fill(NoFragment, val);
	errors[NoFragment] = 1;
	crateErr[crate] |= (1 << NoFragment);
      } else {
        hist->Fill(NoPayload, val);
        errors[NoPayload] = 1;
	crateErr[crate] |= (1 << NoPayload);
      }
    }
  }

  // Update ROB Status and Unpacking Errors

  if (errVecTES && !errVecTES->empty()) {
    const int robMap[32] = {18,  0,  1,  2,  3,  4, 18, 18,
                            18, 18, 18, 18, 18, 18, 18, 18,
			    18,  5,  6,  7,  8,  9, 10, 11,
			    18, 12, 13, 14, 15, 16, 17, 18};
    ROBErrorCollection::const_iterator robIter  = errVecTES->begin();
    ROBErrorCollection::const_iterator robIterE = errVecTES->end();
    unsigned int numRobErr = *robIter;
    ++robIter;
    while (robIter != robIterE) {
      const unsigned int sourceId = *robIter;
      ++robIter;
      if (robIter != robIterE) {
        const int crate = sourceId & 0xf;
        const int slink = (sourceId >> 4) & 0x3;
        const int dataType = (sourceId >> 7) & 0x1;
        const int rob = crate + dataType*6;
        const int pos = rob*4 + slink;
	unsigned int err = *robIter;
	++robIter;
	if (err == 0) continue;
        TH2F* hist = (numRobErr) ? m_h_ROD_PP_rob : m_h_ROD_PP_unp;
        int val = pos;
        if (pos >= 72) {
          hist = (numRobErr) ? m_h_ROD_RoI_rob : m_h_ROD_RoI_unp;
	  val = (pos-72)/2 + 8;
        } else if (pos >= 56) {
          hist = (numRobErr) ? m_h_ROD_RoI_rob : m_h_ROD_RoI_unp;
	  val = (pos-56)/2;
        } else if (pos >= 48) {
          hist = (numRobErr) ? m_h_ROD_CPJEP_rob : m_h_ROD_CPJEP_unp;
	  val = pos-48 + 8;
        } else if (pos >= 32) {
          hist = (numRobErr) ? m_h_ROD_CPJEP_rob : m_h_ROD_CPJEP_unp;
	  val = (pos-32)/2;
        }
	if (numRobErr) {
	  for (int bit = 0; bit < 32; ++bit) {
	    const int robErr = (err >> bit) & 0x1;
	    if (robErr) {
	      hist->Fill(robMap[bit], val);
	      errorsROB[robMap[bit]] = 1;
	    }
	  }
	  crateErr[crate] |= (1 << ROBStatusError);
	  numRobErr--;
        } else {
	  if (err > 17) err = 18;
	  hist->Fill(err, val);
	  errorsUnpack[err] = 1;
	  crateErr[crate] |= (1 << UnpackingError);
        }
      }
    }
  }

  // Update summary plots

  for (int i = 0; i < NumberOfStatusBins; ++i) {
    if (errors[i]) m_h_ROD_summary->Fill(i);
  }
  for (int i = 0; i < 19; ++i) {
    if (errorsROB[i]) m_h_ROB_summary->Fill(i);
  }
  for (int i = 1; i < 19; ++i) {
    if (errorsUnpack[i]) m_h_Unp_summary->Fill(i);
  }

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

  // CPM Error data
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

  // ROD Error data

  for (int crate = 0; crate < ppmCrates+cpmCrates+jemCrates; ++crate) {
    int err = crateErr[crate];
    if (err == 0) continue;
    //if (err & 0x7f) m_h_global->Fill(RODStatus, crate);
    if (err & 0x3f) m_h_global->Fill(RODStatus, crate);
    if (((err >> NoFragment) & 0x1) || ((err >> NoPayload) & 0x1))
                    m_h_global->Fill(RODMissing, crate);
    if ((err >> ROBStatusError) & 0x1) m_h_global->Fill(ROBStatus, crate);
    if ((err >> UnpackingError) & 0x1) m_h_global->Fill(Unpacking, crate);
  }


  m_log << MSG::DEBUG << "Leaving fillHistograms" << endreq;

  return StatusCode::SUCCESS;

}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloRodMonTool::procHistograms(bool isEndOfEventsBlock,
                                  bool isEndOfLumiBlock, bool isEndOfRun)
/*---------------------------------------------------------*/
{
  m_log << MSG::DEBUG << "procHistograms entered" << endreq;

  if (isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun) {
  }

  return StatusCode::SUCCESS;
}

TH1F* TrigT1CaloRodMonTool::book1F(const std::string& name,
                                   const std::string& title,
                                   int nx, double xmin, double xmax)
{
  TH1F *hist = new TH1F(name.c_str(), title.c_str(), nx, xmin, xmax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    m_log << MSG::WARNING << "Could not register histogram : " 
	  << name << endreq;
  }
  hist->SetOption("hist");
  hist->SetStats(kFALSE);
  
  return hist;
}

TH2F* TrigT1CaloRodMonTool::book2F(const std::string& name,
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

void TrigT1CaloRodMonTool::setLabelsStatus(TH1* hist)
{
  hist->GetXaxis()->SetBinLabel(1+GLink,         "GLinkError");
  //hist->GetXaxis()->SetBinLabel(1+CMMParity,     "CMMParityError");
  hist->GetXaxis()->SetBinLabel(1+LVDSLink,      "LVDSLinkError");
  hist->GetXaxis()->SetBinLabel(1+FIFOOverflow,  "RODFIFOOverflow");
  hist->GetXaxis()->SetBinLabel(1+DataTransport, "DataTransportError");
  hist->GetXaxis()->SetBinLabel(1+Timeout,       "GLinkTimeout");
  hist->GetXaxis()->SetBinLabel(1+BCNMismatch,   "BCNMismatch");
  hist->GetXaxis()->SetBinLabel(1+TriggerType,   "TriggerTypeTimeout");
  hist->GetXaxis()->SetBinLabel(1+LimitedRoI,    "LimitedRoISet");
  hist->GetXaxis()->SetBinLabel(1+NoFragment,    "No ROD Fragment");
}

void TrigT1CaloRodMonTool::setLabelsCSL(TH1* hist, bool xAxis, int firstBin,
                                     int lastBin, int binIncr, int slinkIncr)
{
  int crate = 0;
  int slink = 0;
  for (int bin = firstBin; bin <= lastBin; bin += binIncr) {
    std::ostringstream cnum;
    cnum << crate << "/" << slink;
    if (xAxis) hist->GetXaxis()->SetBinLabel(bin, cnum.str().c_str());
    else       hist->GetYaxis()->SetBinLabel(bin, cnum.str().c_str());
    slink += slinkIncr;
    if (slink >= 4) {
      slink -= 4;
      ++crate;
    }
  }
}

void TrigT1CaloRodMonTool::setLabelsROBStatus(TH1* hist)
{
  hist->GetXaxis()->SetBinLabel(1,  "BCIDCheck");
  hist->GetXaxis()->SetBinLabel(2,  "EL1IDCheck");
  hist->GetXaxis()->SetBinLabel(3,  "Timeout");
  hist->GetXaxis()->SetBinLabel(4,  "DataError");
  hist->GetXaxis()->SetBinLabel(5,  "Overflow");
  hist->GetXaxis()->SetBinLabel(6,  "FragmentSize");
  hist->GetXaxis()->SetBinLabel(7,  "DataBlock");
  hist->GetXaxis()->SetBinLabel(8,  "ControlWord");
  hist->GetXaxis()->SetBinLabel(9,  "MissingBOF");
  hist->GetXaxis()->SetBinLabel(10, "MissingEOF");
  hist->GetXaxis()->SetBinLabel(11, "HeaderMarker");
  hist->GetXaxis()->SetBinLabel(12, "Format");
  hist->GetXaxis()->SetBinLabel(13, "Sequence");
  hist->GetXaxis()->SetBinLabel(14, "Transmission");
  hist->GetXaxis()->SetBinLabel(15, "Truncated");
  hist->GetXaxis()->SetBinLabel(16, "Short");
  hist->GetXaxis()->SetBinLabel(17, "Lost");
  hist->GetXaxis()->SetBinLabel(18, "Pending");
  hist->GetXaxis()->SetBinLabel(19, "Unknown");
}

void TrigT1CaloRodMonTool::setLabelsUnpacking(TH1* hist)
{
  hist->GetXaxis()->SetBinLabel(1,  "MissingHeader");
  hist->GetXaxis()->SetBinLabel(2,  "VersionFormat");
  hist->GetXaxis()->SetBinLabel(3,  "Channels");
  hist->GetXaxis()->SetBinLabel(4,  "MissingData");
  hist->GetXaxis()->SetBinLabel(5,  "CrateNumber");
  hist->GetXaxis()->SetBinLabel(6,  "ModuleNumber");
  hist->GetXaxis()->SetBinLabel(7,  "SliceNumber");
  hist->GetXaxis()->SetBinLabel(8,  "RoIWord");
  hist->GetXaxis()->SetBinLabel(9,  "RODnstatus");
  hist->GetXaxis()->SetBinLabel(10, "DuplicateROB");
  hist->GetXaxis()->SetBinLabel(11, "Version");
  hist->GetXaxis()->SetBinLabel(12, "Format");
  hist->GetXaxis()->SetBinLabel(13, "ComprVersion");
  hist->GetXaxis()->SetBinLabel(14, "ComprSlices");
  hist->GetXaxis()->SetBinLabel(15, "DataTruncated");
  hist->GetXaxis()->SetBinLabel(16, "SourceID");
  hist->GetXaxis()->SetBinLabel(17, "WordID");
  hist->GetXaxis()->SetBinLabel(18, "Unknown");
}
