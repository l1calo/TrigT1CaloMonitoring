// ********************************************************************
//
// NAME:        CMMMon.cxx
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of the JEP on CMM level
//
// ********************************************************************

#include <cmath>
#include <vector>

#include "LWHists/LWHist.h"
#include "LWHists/TH1F_LW.h"
#include "LWHists/TH2F_LW.h"
#include "LWHists/TH2I_LW.h"
#include "LWHists/TProfile2D_LW.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/StatusCode.h"
#include "SGTools/StlVectorClids.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigT1CaloMonitoring/CMMMon.h"
#include "TrigT1CaloMonitoring/TrigT1CaloMonErrorTool.h"
#include "TrigT1CaloMonitoring/TrigT1CaloLWHistogramTool.h"

#include "TrigT1CaloEvent/CMMJetHits.h"
#include "TrigT1CaloEvent/CMMEtSums.h"
#include "TrigT1CaloEvent/CMMRoI.h"
#include "TrigT1CaloUtils/QuadLinear.h"
#include "TrigT1CaloUtils/DataError.h"
#include "TrigT1CaloUtils/CrateEnergy.h"

#include "TrigT1Interfaces/TrigT1CaloDefs.h"


// *********************************************************************
// Public Methods
// *********************************************************************

/*---------------------------------------------------------*/
CMMMon::CMMMon( const std::string & type, const std::string & name,
		const IInterface* parent )
  : ManagedMonitorToolBase( type, name, parent ),
    m_errorTool("TrigT1CaloMonErrorTool"),
    m_histTool("TrigT1CaloLWHistogramTool")
/*---------------------------------------------------------*/
{
  // This is how you declare the parameters to Gaudi so that
  // they can be over-written via the job options file

  declareProperty( "CMMJetHitsLocation",
           m_CMMJetHitsLocation = LVL1::TrigT1CaloDefs::CMMJetHitsLocation);
  declareProperty( "CMMEtSumsLocation",
           m_CMMEtSumsLocation  = LVL1::TrigT1CaloDefs::CMMEtSumsLocation);  
  declareProperty( "CMMRoILocation",
           m_CMMRoILocation     = LVL1::TrigT1CaloDefs::CMMRoILocation);

  declareProperty( "PathInRootFile", m_PathInRootFile="L1Calo/JEM_CMM");
  declareProperty( "ErrorPathInRootFile",
                   m_ErrorPathInRootFile="L1Calo/JEM_CMM/Errors/Hardware");
}

/*---------------------------------------------------------*/
CMMMon::~CMMMon()
/*---------------------------------------------------------*/
{
}

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

/*---------------------------------------------------------*/
StatusCode CMMMon::initialize()
/*---------------------------------------------------------*/
{
  msg(MSG::INFO) << "Initializing " << name() << " - package version "
                 << PACKAGE_VERSION << endreq;

  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;

  sc = m_errorTool.retrieve();
  if( sc.isFailure() ) {
    msg(MSG::ERROR) << "Unable to locate Tool TrigT1CaloMonErrorTool"
                    << endreq;
    return sc;
  }

  sc = m_histTool.retrieve();
  if( sc.isFailure() ) {
    msg(MSG::ERROR) << "Unable to locate Tool TrigT1CaloLWHistogramTool"
                    << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode CMMMon::bookHistograms( bool isNewEventsBlock, 
				   bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  msg(MSG::DEBUG) << "in CMMMon::bookHistograms" << endreq;
  
  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
  
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }
  
  if ( isNewEventsBlock|| isNewLumiBlock) { }

  if( isNewRun ) {	

    MonGroup CMM_inputThresh( this, m_PathInRootFile+"/Input/Thresholds",
                                                                expert, run );
    MonGroup CMM_inputEnergy( this, m_PathInRootFile+"/Input/EnergySums",
                                                                expert, run );
    MonGroup CMM_jet( this, m_PathInRootFile+"/Output/Jet", expert, run );
    MonGroup CMM_energy( this, m_PathInRootFile+"/Output/Energy", expert, run );
    MonGroup CMM_RoI( this, m_PathInRootFile+"/Output/RoI", shift, run );
    MonGroup CMM_error( this, m_ErrorPathInRootFile, shift, run );
    MonGroup CMM_errorDetail( this, m_ErrorPathInRootFile, expert, run );
    MonGroup CMM_errorEvents( this, m_ErrorPathInRootFile, expert, run, "",
                                                                "eventSample" );

    m_histTool->setMonGroup(&CMM_inputThresh);

    m_h_CMMJetHits_JEM_MainHits = m_histTool->bookMainJetThresholds(
      "cmm_1d_thresh_MainHits",
      "Main Jet Multiplicity per Threshold  --  CMM input");
    m_h_CMMJetHits_JEM_FwdHitsRight = m_histTool->bookForwardJetThresholds(
      "cmm_1d_thresh_FwdHitsRight",
      "Forward Right Jet Multiplicity per Threshold  --  CMM input");
    m_h_CMMJetHits_JEM_FwdHitsLeft = m_histTool->bookBackwardJetThresholds(
      "cmm_1d_thresh_FwdHitsLeft",
      "Forward Left Jet Multiplicity per Threshold  --  CMM input");

    m_histTool->setMonGroup(&CMM_inputEnergy);

    m_h_CMMEtSums_JEM_Ex = m_histTool->bookJEMQuadLinear(
      "cmm_1d_energy_SubSumsEx", "CMM E_{x}^{JEM}  --  CMM input;Ex [GeV]");
    m_h_CMMEtSums_JEM_Ey = m_histTool->bookJEMQuadLinear(
      "cmm_1d_energy_SubSumsEy", "CMM E_{y}^{JEM}  --  CMM input;Ey [GeV]");
    m_h_CMMEtSums_JEM_Et = m_histTool->bookJEMQuadLinear(
      "cmm_1d_energy_SubSumsEt", "CMM E_{t}^{JEM}  --  CMM input;Et [GeV]");

    //-------------------------- CMM output to DAQ ---------------------------

    m_histTool->setMonGroup(&CMM_jet);

    m_h_CMMJetHits_MainJets = m_histTool->bookMainJetThresholds(
      "cmm_1d_thresh_TotalMainHits",
      "Main Jet Multiplicity per Threshold  --  CMM DAQ");
    m_h_CMMJetHits_FwdJetsRight = m_histTool->bookForwardJetThresholds(
      "cmm_1d_thresh_TotalFwdHitsRight",
      "Forward Right Jet Multiplicity per Threshold  --  CMM DAQ");
    m_h_CMMJetHits_FwdJetsLeft = m_histTool->bookBackwardJetThresholds(
      "cmm_1d_thresh_TotalFwdHitsLeft",
      "Forward Left Jet Multiplicity per Threshold  --  CMM DAQ");
    m_h_CMMJetHits_EtMap = m_histTool->bookJetEtThresholds(
      "cmm_1d_thresh_JetEtHits",
      "JetEt Multiplicity per Threshold  --  CMM DAQ");

    m_histTool->setMonGroup(&CMM_energy);

    m_h_CMMEtSums_MissingEtMap = m_histTool->bookMissingEtThresholds(
      "cmm_1d_energy_MissingEtHits",
      "MissingEt Multiplicity per Threshold  --  CMM DAQ");
    m_h_CMMEtSums_SumEtMap = m_histTool->bookSumEtThresholds(
      "cmm_1d_energy_SumEtHits",
      "SumEt Multiplicity per Threshold  --  CMM DAQ");
    m_h_CMMEtSums_MissingEtSigMap = m_histTool->bookMissingEtSigThresholds(
      "cmm_1d_energy_MissingEtSigHits",
      "MissingEtSig Multiplicity per Threshold  --  CMM DAQ");

    m_h_CMMEtSums_Ex = m_histTool->bookJEMQuadLinear("cmm_1d_energy_TotalEx",
      "E_{x}^{CMM}  --  CMM DAQ;Ex [GeV]", 8);
    m_h_CMMEtSums_Ey = m_histTool->bookJEMQuadLinear("cmm_1d_energy_TotalEy",
      "E_{y}^{CMM}  --  CMM DAQ;Ey [GeV]", 8);
    m_h_CMMEtSums_Et = m_histTool->bookJEMQuadLinear("cmm_1d_energy_TotalEt",
      "SumE_{t}^{CMM}  --  CMM DAQ;Et [GeV]", 8);
    m_h_CMMEtSums_Overflow = m_histTool->bookProfile2D("cmm_2d_energy_Overflow",
      "CMM Energy Overflow Rates", 3, 0., 3., 3, 0., 3.);
    LWHist::LWHistAxis* axis = m_h_CMMEtSums_Overflow->GetXaxis();
    axis->SetBinLabel(1, "Ex");
    axis->SetBinLabel(2, "Ey");
    axis->SetBinLabel(3, "Et");
    axis = m_h_CMMEtSums_Overflow->GetYaxis();
    axis->SetBinLabel(1, "Remote");
    axis->SetBinLabel(2, "Local");
    axis->SetBinLabel(3, "Total");

    //--------------------------- CMM output to RoI --------------------------

    m_histTool->setMonGroup(&CMM_RoI);

    m_h_CMMRoI_JetEtHits = m_histTool->bookJetEtThresholds(
      "cmm_1d_roi_JetEtHits",
      "JetEt Multiplicity per Threshold  --  CMM RoI");
    m_h_CMMRoI_MissingEtHits = m_histTool->bookMissingEtThresholds(
      "cmm_1d_roi_MissingEtHits",
      "MissingEt Multiplicity per Threshold  --  CMM RoI");
    m_h_CMMRoI_SumEtHits = m_histTool->bookSumEtThresholds(
      "cmm_1d_roi_SumEtHits",
      "SumEt Multiplicity per Threshold  --  CMM RoI");
    m_h_CMMRoI_MissingEtSigHits = m_histTool->bookMissingEtSigThresholds(
      "cmm_1d_roi_MissingEtSigHits",
      "MissingEtSig Multiplicity per Threshold  --  CMM RoI");

    m_h_CMMRoI_Ex = m_histTool->bookJEMQuadLinear("cmm_1d_roi_Ex",
      "E_{x}^{CMM}  --  CMM RoI;Ex [GeV]", 8);
    m_h_CMMRoI_Ey = m_histTool->bookJEMQuadLinear("cmm_1d_roi_Ey",
      "E_{y}^{CMM}  --  CMM RoI;Ey [GeV]", 8);
    m_h_CMMRoI_Et = m_histTool->bookJEMQuadLinear("cmm_1d_roi_Et",
      "SumE_{t}^{CMM}  --  CMM RoI;Et [GeV]", 8);

    //---------------------------- S-Link errors -----------------------------

    m_histTool->setMonGroup(&CMM_errorDetail);

    m_h_CMMJet_error = m_histTool->bookJEMSubStatusVsCrate(
      "cmm_2d_thresh_Status", "Errors from CMM Jet SubStatus Word");
    m_h_CMMEnergy_error = m_histTool->bookJEMSubStatusVsCrate(
      "cmm_2d_energy_Status", "Errors from CMM Energy SubStatus Word");
    m_h_CMMJet_parity = m_histTool->book2F("cmm_2d_thresh_Parity",
      "CMM Jet Parity Errors;Module or Remote;Crate", 18, 0., 18., 2, 0., 2.);
    m_histTool->numbers(m_h_CMMJet_parity, 0, 15);
    m_h_CMMJet_parity->GetXaxis()->SetBinLabel(17, "Main");
    m_h_CMMJet_parity->GetXaxis()->SetBinLabel(18, "Fwd");
    m_histTool->numbers(m_h_CMMJet_parity, 0, 1, 1, 0, false);
    m_h_CMMEnergy_parity = m_histTool->book2F("cmm_2d_energy_Parity",
      "CMM Energy Parity Errors;Module or Remote;Crate",
      19, 0., 19., 2, 0., 2.);
    m_histTool->numbers(m_h_CMMEnergy_parity, 0, 15);
    m_h_CMMEnergy_parity->GetXaxis()->SetBinLabel(17, "Ex");
    m_h_CMMEnergy_parity->GetXaxis()->SetBinLabel(18, "Ey");
    m_h_CMMEnergy_parity->GetXaxis()->SetBinLabel(19, "Et");
    m_histTool->numbers(m_h_CMMEnergy_parity, 0, 1, 1, 0, false);

    m_h_CMMRoI_error = m_histTool->book1F("cmm_1d_roi_Parity", 
      "CMM RoI Parity Errors",
      NumberOfRoIParityBins, 0., NumberOfRoIParityBins);
    axis = m_h_CMMRoI_error->GetXaxis();
    axis->SetBinLabel(1+ExParity,    "Ex, MissingEtSig");
    axis->SetBinLabel(1+EyParity,    "Ey, SumEt");
    axis->SetBinLabel(1+EtParity,    "Et, MissingEt");
    axis->SetBinLabel(1+JetEtParity, "JetEt");

    m_h_TriggeredSlice = m_histTool->book1F("cmm_1d_TriggeredSlices",
      "Comparison of CMM Jet and Energy triggered slice numbers;Difference",
      5, 0., 5.);
    m_histTool->numbers(m_h_TriggeredSlice, 0, 4);

    m_histTool->setMonGroup(&CMM_error);

    //Error Summary for all CMMs in system
    m_h_CMM_ErrorSummary = m_histTool->book1F("cmm_1d_ErrorSummary",
      "Error Summary of CMM Jet, Energy and RoI path",
      NumberOfSummaryBins, 0., NumberOfSummaryBins);

    m_histTool->setMonGroup(&CMM_errorEvents);

    m_h_CMM_Events = m_histTool->bookEventNumbers("cmm_2d_ErrorEventNumbers",
      "JEM-CMM Error Event Numbers",
      NumberOfSummaryBins, 0., NumberOfSummaryBins);

    axis = m_h_CMM_ErrorSummary->GetXaxis();
    for (int i = 0; i < 2; ++i) {
      axis->SetBinLabel(1+JetStatus,    "Jet status");
      axis->SetBinLabel(1+EnergyStatus, "Energy status");
      axis->SetBinLabel(1+JetParity,    "Jet parity");
      axis->SetBinLabel(1+EnergyParity, "Energy parity");
      axis->SetBinLabel(1+RoIParity,    "RoI parity");
      axis = m_h_CMM_Events->GetYaxis();
    }
  }
  
  return StatusCode::SUCCESS;
}


/*---------------------------------------------------------*/
StatusCode CMMMon::fillHistograms()
/*---------------------------------------------------------*/
{
  const bool debug = msgLvl(MSG::DEBUG);

  // Skip events believed to be corrupt

  if (m_errorTool->corrupt()) {
    if (debug) msg(MSG::DEBUG) << "Skipping corrupt event" << endreq;
    return StatusCode::SUCCESS;
  }

  // Error vector for global overview
  std::vector<int> overview(2);

  // triggered slice numbers
  int j_num_slice = -1;
  int e_num_slice = -1;

  using LVL1::DataError;

  // =========================================================================
  // ================= Container: CMM Jet Hits ===============================
  // =========================================================================

  // retrieve CMM Jet Hits from Storegate
  const CMMJetHitsCollection* CMMJetHits = 0;
  StatusCode sc = evtStore()->retrieve(CMMJetHits, m_CMMJetHitsLocation);
  if (sc == StatusCode::FAILURE) {
    msg(MSG::INFO) << "No CMM JetHits found in TES at "
                   << m_CMMJetHitsLocation << endreq;
    return StatusCode::SUCCESS;
  }  

  if (debug) {
    msg(MSG::DEBUG) << "--------------  CMM Jet Hits ---------------"<<endreq;  
  }

  CMMJetHitsCollection::const_iterator it_CMMJetHits ;
  
  // Step over all cells
  for (it_CMMJetHits = CMMJetHits->begin(); it_CMMJetHits != CMMJetHits->end();
                                                          ++it_CMMJetHits) {	  
    const int crate  = (*it_CMMJetHits)-> crate();
    const int dataID = (*it_CMMJetHits)-> dataID();
    const unsigned int jetHits = (*it_CMMJetHits)-> Hits();

    if (debug) {
      msg(MSG::DEBUG) << "CMMJetHits Crate: " << crate
                      << " dataID: "          << dataID
		      << "   Hits: ";
    }
      
    // -----------------------------------------------------------------------
    // --------- Histos with distribution of JEM Hit Multiplicities ----------
    // -----------------------------------------------------------------------
    
    //input data from JEMs have dataID 0..15
    if (dataID < 16) {

      const bool forward = (dataID%8 == 0 || dataID%8 == 7);
      const int nBits = (forward) ?  2 : 3;
      m_histTool->fillThresholds(m_h_CMMJetHits_JEM_MainHits, jetHits, 8,
                                                                       nBits);
      if (forward) {
        TH1F_LW* hist = (dataID%8 == 0) ? m_h_CMMJetHits_JEM_FwdHitsLeft
	                                : m_h_CMMJetHits_JEM_FwdHitsRight;
        m_histTool->fillThresholds(hist, (jetHits >> 16), 4, nBits);
      }

      if (debug) {
        const int nHits = (forward) ? 12 : 8;
	msg(MSG::DEBUG) << m_histTool->thresholdString(jetHits, nHits, nBits)
	                << endreq;
      }
      
    // -----------------------------------------------------------------------
    // --------- Histos with distribution of CMM Hit Multiplicities ----------
    // -----------------------------------------------------------------------

    } else {

      if (dataID == LVL1::CMMJetHits::TOTAL_MAIN) {
        m_histTool->fillThresholds(m_h_CMMJetHits_MainJets, jetHits, 8, 3);
      } else if (dataID == LVL1::CMMJetHits::TOTAL_FORWARD) {
        m_histTool->fillThresholds(m_h_CMMJetHits_FwdJetsLeft, jetHits, 4, 2);
	m_histTool->fillThresholds(m_h_CMMJetHits_FwdJetsRight, (jetHits >> 8),
	                                                                4, 2);
      } else if (dataID == LVL1::CMMJetHits::ET_MAP) {
        m_histTool->fillThresholds(m_h_CMMJetHits_EtMap, jetHits, 4, 1);
      }

      if (debug) {
        const int nHits = (dataID == LVL1::CMMJetHits::ET_MAP) ? 4 : 8;
        const int nBits = (dataID == LVL1::CMMJetHits::TOTAL_MAIN ||
	                   dataID == LVL1::CMMJetHits::LOCAL_MAIN ||
		           dataID == LVL1::CMMJetHits::REMOTE_MAIN) ? 3 :
                              ((dataID == LVL1::CMMJetHits::ET_MAP) ? 1 : 2);
	msg(MSG::DEBUG) << m_histTool->thresholdString(jetHits, nHits, nBits)
	                << endreq;
      }
    }

    // -----------------------------------------------------------------------
    // ----------------- Histos with SubStatus Word errors -------------------
    // -----------------------------------------------------------------------
       
    if (j_num_slice < 0) j_num_slice = (*it_CMMJetHits)->peak();

    const DataError err((*it_CMMJetHits)->Error());
	  
    //Error summary plots
    //substatus word
    const int status = (err.error() >> LVL1::DataError::GLinkParity) & 0xff;
    if (status) {
      for (int bit = 0; bit < 8; ++bit) {
        if ((status >> bit) & 0x1) m_h_CMMJet_error->Fill(bit, crate);
      }
      m_h_CMM_ErrorSummary->Fill(JetStatus);
      m_histTool->fillEventNumber(m_h_CMM_Events, JetStatus);
      overview[crate] |= (1 << JetStatus);
    }

    if (dataID < 16 || dataID == LVL1::CMMJetHits::REMOTE_MAIN
                    || dataID == LVL1::CMMJetHits::REMOTE_FORWARD) {
      // Parity
      if (err.get(DataError::Parity)) {
        const int xpos = (dataID < 16)
	                 ? dataID
	                 : (dataID == LVL1::CMMJetHits::REMOTE_MAIN) ? 16 : 17;
	m_h_CMMJet_parity->Fill(xpos, crate);
	m_h_CMM_ErrorSummary->Fill(JetParity);
	m_histTool->fillEventNumber(m_h_CMM_Events, JetParity);
	overview[crate] |= (1 << JetParity);
      }
    }

  }
  
  // =========================================================================
  // ================= Container: CMM Et Sums ================================
  // =========================================================================

  // retrieve CMM Et Sums from Storegate
  const CMMEtSumsCollection* CMMEtSums = 0;
  sc = evtStore()->retrieve(CMMEtSums, m_CMMEtSumsLocation);
  if (sc == StatusCode::FAILURE) {
    msg(MSG::INFO) << "No CMMEtSums found in TES at "
                   << m_CMMEtSumsLocation << endreq ;
    return StatusCode::SUCCESS;
  }
  
  msg(MSG::DEBUG) << "-------------- CMM Et Sums ---------------" << endreq;
  
  // Step over all cells 
  CMMEtSumsCollection::const_iterator it_CMMEtSums ;
  for (it_CMMEtSums = CMMEtSums->begin(); it_CMMEtSums != CMMEtSums->end();
                                                       ++it_CMMEtSums) {	  
    const int crate  = (*it_CMMEtSums)->crate();
    const int dataID = (*it_CMMEtSums)->dataID();
    const unsigned int rawEx = (*it_CMMEtSums)->Ex();
    const unsigned int rawEy = (*it_CMMEtSums)->Ey();
    const unsigned int rawEt = (*it_CMMEtSums)->Et();
    const int exError = (*it_CMMEtSums)->ExError();
    const int eyError = (*it_CMMEtSums)->EyError();
    const int etError = (*it_CMMEtSums)->EtError();

    // -----------------------------------------------------------------------
    // -------------- Histos with distribution of JEM Energies ---------------
    // -----------------------------------------------------------------------
    // JEM energy sums, dataID < 16
    if (dataID < 16) {
      // note: JEM energies are compressed -> use QuadLinear to expand!
      const int ex = LVL1::QuadLinear::Expand(rawEx);
      const int ey = LVL1::QuadLinear::Expand(rawEy);
      const int et = LVL1::QuadLinear::Expand(rawEt);
	  
      if (ex > 0) m_h_CMMEtSums_JEM_Ex->Fill(ex, 1.);
      if (ey > 0) m_h_CMMEtSums_JEM_Ey->Fill(ey, 1.);
      if (et > 0) m_h_CMMEtSums_JEM_Et->Fill(et, 1.);
    }
      
    // -----------------------------------------------------------------------
    // ---------- Histos with distribution of total Energy per system --------
    // -----------------------------------------------------------------------
    // total energy sums
    if (dataID == LVL1::CMMEtSums::TOTAL && crate == 1) {
      // Use CrateEnergy object to decode 15-bit twos-complement format
      LVL1::CrateEnergy cen(crate, rawEt, rawEx, rawEy,
                            etError&0x1, exError&0x1, eyError&0x1);
      const int ex = std::abs(cen.ex());
      const int ey = std::abs(cen.ey());
      const int et = rawEt;

      if (ex > 0 && !cen.exOverflow()) m_h_CMMEtSums_Ex->Fill(ex, 1.);
      if (ey > 0 && !cen.eyOverflow()) m_h_CMMEtSums_Ey->Fill(ey, 1.);
      if (et > 0 && !cen.etOverflow()) m_h_CMMEtSums_Et->Fill(et, 1.);

      if (debug) {
        msg(MSG::DEBUG) << "    Ex: " << ex << "; Ey: " << ey << "; Et " << et
                        << endreq;
        msg(MSG::DEBUG) << "raw Ex: " << rawEx << "; Ey: " << rawEy << "; Et "
                        << rawEt << endreq;
      }
    }
    //MissingEt/SumEt/MissingEtSig Hitmaps
    if ((dataID == LVL1::CMMEtSums::MISSING_ET_MAP ||
         dataID == LVL1::CMMEtSums::SUM_ET_MAP ||
	 dataID == LVL1::CMMEtSums::MISSING_ET_SIG_MAP) && crate == 1) {
      const int nHits = 8;
      TH1F_LW* hist = (dataID == LVL1::CMMEtSums::MISSING_ET_MAP)
                       ? m_h_CMMEtSums_MissingEtMap
		       : (dataID == LVL1::CMMEtSums::SUM_ET_MAP)
		       ? m_h_CMMEtSums_SumEtMap : m_h_CMMEtSums_MissingEtSigMap;
      m_histTool->fillThresholds(hist, rawEt, nHits, 1);

      if (debug) {
        if (dataID == LVL1::CMMEtSums::MISSING_ET_MAP) {
	  msg(MSG::DEBUG) << "MissingEt Hits: ";
        } else if (dataID == LVL1::CMMEtSums::SUM_ET_MAP) {
	  msg(MSG::DEBUG) << "SumEt Hits: ";
	} else msg(MSG::DEBUG) << "MissingEtSig Hits: ";
	msg(MSG::DEBUG) << m_histTool->thresholdString(rawEt, nHits) << endreq;
      }
    }
      
    if (e_num_slice < 0) {
      e_num_slice = (*it_CMMEtSums)-> peak();
      if (j_num_slice >= 0) {
         m_h_TriggeredSlice->Fill(std::abs(e_num_slice - j_num_slice));
      }
    }

    // -----------------------------------------------------------------------
    // --------------- Histos with SubStatus Word errors ---------------------
    // -----------------------------------------------------------------------
    //Error summary plots
    //substatus word
    const DataError exerr(exError);
    const DataError eyerr(eyError);
    const DataError eterr(etError);
    const int status = (eterr.error() >> LVL1::DataError::GLinkParity) & 0xff;
    if (status) {
      for (int bit = 0; bit < 8; ++bit) {
        if ((status >> bit) & 0x1) m_h_CMMEnergy_error->Fill(bit, crate);
      }
      m_h_CMM_ErrorSummary->Fill(EnergyStatus);
      m_histTool->fillEventNumber(m_h_CMM_Events, EnergyStatus);
      overview[crate] |= (1 << EnergyStatus);
    }

    if (dataID < 16 || dataID == LVL1::CMMEtSums::REMOTE) {
      // Parity
      if (eterr.get(DataError::Parity) || exerr.get(DataError::Parity) ||
          eyerr.get(DataError::Parity)) {
        if (dataID < 16) m_h_CMMEnergy_parity->Fill(dataID, crate);
	else {
	  if (exerr.get(DataError::Parity)) {
	    m_h_CMMEnergy_parity->Fill(16, crate);
	  }
	  if (eyerr.get(DataError::Parity)) {
	    m_h_CMMEnergy_parity->Fill(17, crate);
          }
	  if (eterr.get(DataError::Parity)) {
	    m_h_CMMEnergy_parity->Fill(18, crate);
	  }
	}
	m_h_CMM_ErrorSummary->Fill(EnergyParity);
	m_histTool->fillEventNumber(m_h_CMM_Events, EnergyParity);
	overview[crate] |= (1 << EnergyParity);
      }
    }

    // Overflow - not an error, plot rate
    if (crate == 1 && (dataID == LVL1::CMMEtSums::REMOTE ||
                       dataID == LVL1::CMMEtSums::LOCAL  ||
	               dataID == LVL1::CMMEtSums::TOTAL)) {
      const double ypos = (dataID == LVL1::CMMEtSums::REMOTE) ? 0.
                        : (dataID == LVL1::CMMEtSums::LOCAL)  ? 1. : 2.;
      m_h_CMMEtSums_Overflow->Fill(0., ypos, exerr.get(DataError::Overflow));
      m_h_CMMEtSums_Overflow->Fill(1., ypos, eyerr.get(DataError::Overflow));
      m_h_CMMEtSums_Overflow->Fill(2., ypos, eterr.get(DataError::Overflow));
    }

  }

  // =========================================================================
  // ================= Container: CMM RoI ====================================
  // =========================================================================
  
  // retrieve RoI information from Storegate
  const LVL1::CMMRoI* CR = 0;
  sc = evtStore()->retrieve (CR, m_CMMRoILocation);
  if (sc == StatusCode::FAILURE) {
    msg(MSG::INFO) << "No CMM RoI found in TES at " << m_CMMRoILocation
                   << endreq;
    return StatusCode::SUCCESS;    
  }

  if (debug) {
    msg(MSG::DEBUG) << "-------------- CMM RoI ---------------" << endreq;
  }

  // -------------------------------------------------------------------------
  // -------------- Histos filled with CMM RoI information -------------------
  // -------------------------------------------------------------------------

  const int rawEx = (CR)->ex();
  const int rawEy = (CR)->ey();
  const int et    = (CR)->et();
  const int exError = (CR)->exError();
  const int eyError = (CR)->eyError();
  const int etError = (CR)->etError();
  const int jetEtHits = (CR)->jetEtHits();
  const int sumEtHits = (CR)->sumEtHits();
  const int missingEtHits = (CR)->missingEtHits();
  const int missingEtSigHits = (CR)->missingEtSigHits();

  m_histTool->fillThresholds(m_h_CMMRoI_JetEtHits, jetEtHits, 4, 1);
  m_histTool->fillThresholds(m_h_CMMRoI_SumEtHits, sumEtHits, 8, 1);
  m_histTool->fillThresholds(m_h_CMMRoI_MissingEtHits, missingEtHits, 8, 1);
  m_histTool->fillThresholds(m_h_CMMRoI_MissingEtSigHits, missingEtSigHits,
                                                                      8, 1);

  if (debug) {
    msg(MSG::DEBUG) << "JetEtHits: "
                    << m_histTool->thresholdString(jetEtHits, 4)
		    << "; SumEtHits: "
		    << m_histTool->thresholdString(sumEtHits, 8)
		    << "; MissingEtHits: "
		    << m_histTool->thresholdString(missingEtHits, 8)
		    << "; MissingEtSigHits: "
		    << m_histTool->thresholdString(missingEtSigHits, 8)
		    << endreq;
  }

  // Use CrateEnergy object to decode 15-bit twos-complement format
  LVL1::CrateEnergy cen(1, et, rawEx, rawEy,
                        etError&0x1, exError&0x1, eyError&0x1);
  const int ex = std::abs(cen.ex());
  const int ey = std::abs(cen.ey());

  if (ex > 0 && !cen.exOverflow()) m_h_CMMRoI_Ex->Fill(ex,1);
  if (ey > 0 && !cen.eyOverflow()) m_h_CMMRoI_Ey->Fill(ey,1);
  if (et > 0 && !cen.etOverflow()) m_h_CMMRoI_Et->Fill(et,1);
 
  if (debug) {
    msg(MSG::DEBUG) << "    Ex: " << ex << "; Ey: " << ey << "; Et " << et
                    << endreq;
    msg(MSG::DEBUG) << "raw Ex: " << rawEx << "; Ey: "<< rawEy << "; Et " << et
                    << endreq;
    msg(MSG::DEBUG) << "CMM Slice numbers: Jet: " << j_num_slice
                    << " Energy: " << e_num_slice << endreq;
  }
  
  // errors
  const DataError exerr(exError);
  const DataError eyerr(eyError);
  const DataError eterr(etError);
  const DataError jetEterr((CR)-> jetEtError());

  // Parity (Ex,MissingEtSigMap)
  if (exerr.get(DataError::Parity)) m_h_CMMRoI_error->Fill(ExParity);
  // Parity (Ey,SumEtMap)
  if (eyerr.get(DataError::Parity)) m_h_CMMRoI_error->Fill(EyParity);
  // Parity (Et,MissingEtMap)
  if (eterr.get(DataError::Parity)) m_h_CMMRoI_error->Fill(EtParity);
  // Parity (JetEtMap)
  if (jetEterr.get(DataError::Parity)) m_h_CMMRoI_error->Fill(JetEtParity);
      
  //Error summary plots
  //parity
  if (exerr.get(DataError::Parity) || eyerr.get(DataError::Parity) ||
      eterr.get(DataError::Parity) || jetEterr.get(DataError::Parity)) {
    m_h_CMM_ErrorSummary->Fill(RoIParity);
    m_histTool->fillEventNumber(m_h_CMM_Events, RoIParity);
    overview[1] |= (1 << RoIParity);
           
  }

  // Write overview vector to StoreGate
  std::vector<int>* save = new std::vector<int>(overview);
  sc = evtStore()->record(save, "L1CaloJEMCMMErrorVector");
  if (sc != StatusCode::SUCCESS) {
    msg(MSG::ERROR) << "Error recording JEM CMM error vector in TES "
                    << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode CMMMon::procHistograms( bool isEndOfEventsBlock, 
				   bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
  msg(MSG::DEBUG) << "in procHistograms" << endreq ;

  if( isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun ) { }
	
  return StatusCode::SUCCESS;
}
