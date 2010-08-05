// ********************************************************************
//
// NAME:        JEMMon.cxx
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of the JEP on JEM level
//
// ********************************************************************

#include <cmath>
#include <set>
#include <sstream>

#include "LWHists/LWHist.h"
#include "LWHists/TH1F_LW.h"
#include "LWHists/TH2F_LW.h"
#include "LWHists/TH2I_LW.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/StatusCode.h"
#include "SGTools/StlVectorClids.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigT1CaloMonitoring/JEMMon.h"
#include "TrigT1CaloMonitoring/TrigT1CaloMonErrorTool.h"
#include "TrigT1CaloMonitoring/TrigT1CaloLWHistogramTool.h"

#include "TrigT1CaloEvent/JEMHits.h"
#include "TrigT1CaloEvent/JEMEtSums.h"
#include "TrigT1CaloEvent/JetElement.h"
#include "TrigT1CaloEvent/JEMRoI.h"
#include "TrigT1CaloUtils/QuadLinear.h"
#include "TrigT1CaloUtils/DataError.h"
#include "TrigT1CaloUtils/CoordToHardware.h"
#include "TrigT1Interfaces/Coordinate.h"
#include "TrigT1Interfaces/JEPRoIDecoder.h"
#include "TrigT1Interfaces/TrigT1CaloDefs.h"
#include "TrigT1Interfaces/Coordinate.h"
#include "TrigConfL1Data/L1DataDef.h"

/*---------------------------------------------------------*/
JEMMon::JEMMon( const std::string & type, const std::string & name,
		const IInterface* parent )
  : ManagedMonitorToolBase( type, name, parent ),
    m_errorTool("TrigT1CaloMonErrorTool"),
    m_histTool("TrigT1CaloLWHistogramTool")
/*---------------------------------------------------------*/
{
  // This is how you declare the parameters to Gaudi so that
  // they can be over-written via the job options file

  declareProperty( "JetElementLocation",
         m_JetElementLocation = LVL1::TrigT1CaloDefs::JetElementLocation); 
  declareProperty( "JEMHitsLocation",
         m_JEMHitsLocation    = LVL1::TrigT1CaloDefs::JEMHitsLocation) ;
  declareProperty( "JEMEtSumsLocation",
         m_JEMEtSumsLocation  = LVL1::TrigT1CaloDefs::JEMEtSumsLocation) ;
  declareProperty( "JEMRoILocation",
         m_JEMRoILocation     = LVL1::TrigT1CaloDefs::JEMRoILocation) ;

  declareProperty( "NumberOfSlices", m_SliceNo = 5);
  declareProperty( "MaxEnergyRange", m_MaxEnergyRange = 1024) ;

  declareProperty( "PathInRootFile", m_PathInRootFile = "L1Calo/JEM") ;
  declareProperty( "ErrorPathInRootFile",
                   m_ErrorPathInRootFile = "L1Calo/JEM/Errors/Hardware") ;

}


/*---------------------------------------------------------*/
JEMMon::~JEMMon()
/*---------------------------------------------------------*/
{
}

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

/*---------------------------------------------------------*/
StatusCode JEMMon::initialize()
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
StatusCode JEMMon::bookHistograms( bool isNewEventsBlock, 
				   bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  msg(MSG::DEBUG) << "in JEMMon::bookHistograms" << endreq;
 
  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  if ( isNewEventsBlock|| isNewLumiBlock) { }

  if( isNewRun ) {	

    MonGroup JetElements_expert(this, m_PathInRootFile+"/Input", expert, run);
    MonGroup JetElements_shift(this, m_PathInRootFile+"/Input", shift, run);
    MonGroup JEM_Thresholds(this, m_PathInRootFile+"/Output/Thresholds",
                                                                 expert, run);
    MonGroup JEM_EnergySums(this, m_PathInRootFile+"/Output/EnergySums",
                                                                 expert, run);
    MonGroup JEM_RoI(this, m_PathInRootFile+"/Output/RoI", shift, run);
    MonGroup JEM_Error(this, m_ErrorPathInRootFile, shift, run );
    MonGroup JEM_ErrorDetail(this, m_ErrorPathInRootFile, expert, run );
    MonGroup JEM_ErrorEvents(this, m_ErrorPathInRootFile, expert, run, "",
                                                              "eventSample" );

    //-------------------------- JetElements histos --------------------------

    m_histTool->setMonGroup(&JetElements_expert);

    m_h_je_emeta = m_histTool->bookJEMEta("jem_em_1d_jetEl_Eta",
      "em TowerSum distribution per #eta  --  JEM input");
    m_h_je_hadeta = m_histTool->bookJEMEta("jem_had_1d_jetEl_Eta",
      "had TowerSum distribution per #eta  --  JEM input");

    m_h_je_emphi = m_histTool->book1F("jem_em_1d_jetEl_Phi",
      "em TowerSum distribution per #phi  --  JEM input;phi", 32, 0., 2.*M_PI);
    m_h_je_hadphi = m_histTool->book1F("jem_had_1d_jetEl_Phi",
      "had TowerSum distribution per #phi  --  JEM input;phi", 32, 0., 2.*M_PI);

    int jebins = 256;
    if (m_MaxEnergyRange < jebins) jebins = m_MaxEnergyRange;
    m_h_je_emenergy = m_histTool->book1F("jem_em_1d_jetEl_Energy",
      "TowerSum EM energy distribution  --  JEM input;em energy [GeV]",
      jebins-1, 1, m_MaxEnergyRange);
    m_h_je_hadenergy = m_histTool->book1F("jem_had_1d_jetEl_Energy",
      "TowerSum HAD energy distribution  --  JEM input;had energy [GeV]",
      jebins-1, 1, m_MaxEnergyRange);

    // number of triggered slice
    m_h_je_triggeredSlice = m_histTool->book1F("JE_TriggeredSlice",
      "Number of the Triggered Slice for JE;#Slice", m_SliceNo, 0, m_SliceNo);
    m_histTool->numbers(m_h_je_triggeredSlice, 0, m_SliceNo-1);

    m_histTool->setMonGroup(&JetElements_shift);

    m_h_je_energy_emHitMap = m_histTool->bookJEMEtaVsPhi(
      "jem_em_2d_etaPhi_jetEl_HitMapWeighted",
      "#eta - #phi map of EM TowerSum weighted with energy  --  JEM input"); 
    m_h_je_energy_hadHitMap = m_histTool->bookJEMEtaVsPhi(
      "jem_had_2d_etaPhi_jetEl_HitMapWeighted",
      "#eta - #phi map of HAD TowerSum weighted with energy  --  JEM input"); 

    std::string name, title;
    std::stringstream buffer;
	  
    m_h_je_emHitMap.clear();
    m_h_je_hadHitMap.clear();
    for (int i = 0; i < m_SliceNo; i++) {
      buffer.str("");
      buffer << i;
	      
      name  = "jem_em_2d_etaPhi_jetEl_HitMapSlice" + buffer.str();
      title = "#eta - #phi map of EM TowerSum for Timeslice " + buffer.str()
                                                         +  "  --  JEM input";
      m_h_je_emHitMap.push_back(m_histTool->bookJEMEtaVsPhi(name,title));
	      
      name  = "jem_had_2d_etaPhi_jetEl_HitMapSlice" + buffer.str();
      title = "#eta - #phi map of HAD TowerSum for Timeslice " + buffer.str()
                                                          +  "  --  JEM input";
      m_h_je_hadHitMap.push_back(m_histTool->bookJEMEtaVsPhi(name,title));
    }
      
    //---------------------------- DAQ histos -----------------------------

    m_histTool->setMonGroup(&JEM_Thresholds);

    m_h_JEMHits_MainHits = m_histTool->bookMainJetThresholds(
      "jem_1d_thresh_MainHits",
      "Main Jet Hit Multiplicity per Threshold  --  JEM DAQ");
    m_h_JEMHits_FwdHitsRight = m_histTool->bookForwardJetThresholds(
      "jem_1d_thresh_FwdHitsRight",
      "Fwd Right Jet Hit Multiplicity per Threshold  --  JEM DAQ");
    m_h_JEMHits_FwdHitsLeft = m_histTool->bookBackwardJetThresholds(
      "jem_1d_thresh_FwdHitsLeft",
      "Fwd Left Jet Hit Multiplicity per Threshold  --  JEM DAQ");
      
    m_h_JEMDAQ_Hits_Map = m_histTool->bookJEMCrateModuleVsThresholds(
      "jem_2d_thresh_HitsPerJem", "HitMap of Hits per JEM");

    m_histTool->setMonGroup(&JEM_EnergySums);
      
    m_h_JEMEtSums_Ex = m_histTool->bookJEMQuadLinear("jem_1d_energy_SubSumsEx",
      "JEM E_{x}^{JEM}  --  JEM DAQ;Ex [GeV]");
    m_h_JEMEtSums_Ey = m_histTool->bookJEMQuadLinear("jem_1d_energy_SubSumsEy",
      "JEM E_{y}^{JEM}  --  JEM DAQ;Ey [GeV]");
    m_h_JEMEtSums_Et = m_histTool->bookJEMQuadLinear("jem_1d_energy_SubSumsEt",
      "JEM E_{t}^{JEM}  --  JEM DAQ;Et [GeV]");
      
    //---------------------------- RoI histos -----------------------------

    m_histTool->setMonGroup(&JEM_RoI);

    m_h_JEMRoI_MainHits = m_histTool->bookMainJetThresholds(
      "jem_1d_roi_MainHits",
      "Main Jet Hit Multiplicity per Threshold  --  JEM RoI");
    m_h_JEMRoI_FwdHitsRight = m_histTool->bookForwardJetThresholds(
      "jem_1d_roi_FwdHitsRight",
      "Forward Right Jet Hit Multiplicity per Threshold  --  JEM RoI");
    m_h_JEMRoI_FwdHitsLeft = m_histTool->bookBackwardJetThresholds(
      "jem_1d_roi_FwdHitsLeft",
      "Forward Left Jet Hit Multiplicity per Threshold  --  JEM RoI");
    m_h_JEMRoI_sat = m_histTool->bookJEMRoIEtaVsPhi(
      "jem_2d_etaPhi_roi_Saturation", "JEM RoI Saturation");

    //----------------------- HitThreshold per Eta-Phi -----------------------

    m_h_JEMRoI_MainThreshPerEtaPhi.clear();
    m_h_JEMRoI_FwdThreshPerEtaPhi.clear();
    std::vector<std::string> jetNames;
    std::vector<std::string> jfNames;
    std::vector<std::string> jbNames;
    m_histTool->thresholdNames(TrigConf::L1DataDef::jetType(), jetNames);
    m_histTool->thresholdNames(TrigConf::L1DataDef::jfType(),  jfNames);
    m_histTool->thresholdNames(TrigConf::L1DataDef::jbType(),  jbNames);
    for (int i = 0; i < 8; i++) {
      buffer.str("");
      buffer << i;
      name  = "jem_2d_etaPhi_roi_MainThresh" + buffer.str();
      title = "#eta - #phi Map of Main Hits passing Threshold "+ jetNames[i]
                                                               +"  --  JEM RoI";
      m_h_JEMRoI_MainThreshPerEtaPhi.push_back(
            m_histTool->bookJEMRoIEtaVsPhi(name.c_str(), title.c_str()));
      if (i >= 4) continue;
      name  = "jem_2d_etaPhi_roi_FwdThresh" + buffer.str();
      title = "#eta - #phi Map of Fwd Hits passing Threshold "+ jfNames[i];
      if (jfNames[i] != jbNames[i]) title += "/" + jbNames[i];
      title +="  --  JEM RoI";
      m_h_JEMRoI_FwdThreshPerEtaPhi.push_back(
            m_histTool->bookJEMRoIEtaVsPhi(name.c_str(), title.c_str()));
    }

    //--------------------------- Error Histos -------------------------------

    m_histTool->setMonGroup(&JEM_ErrorDetail);

    m_h_je_error = m_histTool->bookJEMSubStatusVsCrateModule("jem_2d_Status",
      "Error reports from JEM SubStatus Word");
    m_h_je_em_parity = m_histTool->bookJEMEtaVsPhi(
      "jem_em_2d_etaPhi_jetEl_Parity", "Jet Element EM Parity Errors");
    m_h_je_had_parity = m_histTool->bookJEMEtaVsPhi(
      "jem_had_2d_etaPhi_jetEl_Parity", "Jet Element Had Parity Errors");
    m_h_je_em_link = m_histTool->bookJEMEtaVsPhi(
      "jem_em_2d_etaPhi_jetEl_LinkDown", "Jet Element EM Link Down Errors");
    m_h_je_had_link = m_histTool->bookJEMEtaVsPhi(
      "jem_had_2d_etaPhi_jetEl_LinkDown", "Jet Element Had Link Down Errors");

    m_h_JEMRoI_error = m_histTool->bookJEMRoIEtaVsPhi(
      "jem_2d_etaPhi_roi_Parity", "JEM RoI Parity Errors");
	
    m_histTool->setMonGroup(&JEM_Error);

    m_h_JEM_ErrorSummary = m_histTool->book1F("jem_1d_ErrorSummary",
      "Summary of JEM Data Errors", NumberOfSummaryBins, 0.,
                                                          NumberOfSummaryBins);

    m_histTool->setMonGroup(&JEM_ErrorEvents);

    m_h_JEM_Events = m_histTool->bookEventNumbers("jem_2d_ErrorEventNumbers",
      "JEM Error Event Numbers", NumberOfSummaryBins, 0., NumberOfSummaryBins);

    LWHist::LWHistAxis* axis = m_h_JEM_ErrorSummary->GetXaxis();
    for (int i = 0; i < 2; ++i) {
      axis->SetBinLabel(1+EMParity,  "EM parity");
      axis->SetBinLabel(1+HadParity, "Had parity");
      axis->SetBinLabel(1+EMLink,    "EM link");
      axis->SetBinLabel(1+HadLink,   "Had link");
      axis->SetBinLabel(1+JEMStatus, "JEM status");
      axis->SetBinLabel(1+RoIParity, "RoI parity");
      axis = m_h_JEM_Events->GetYaxis();
    }
       
    m_histTool->unsetMonGroup();
  }
    
  return StatusCode::SUCCESS;
}



/*---------------------------------------------------------*/
StatusCode JEMMon::fillHistograms()
/*---------------------------------------------------------*/
{
  bool debug = msgLvl(MSG::DEBUG);

  // Skip events believed to be corrupt

  if (m_errorTool->corrupt()) {
    if (debug) msg(MSG::DEBUG) << "Skipping corrupt event" << endreq;
    return StatusCode::SUCCESS;
  }

  // Error vector for global overview
  std::vector<int> overview(2);

  using LVL1::DataError;

  // =========================================================================
  // ================= Container: JetElements ================================
  // =========================================================================

  // retrieve JetElements
  const JECollection* jetElements;
  StatusCode sc = evtStore()->retrieve(jetElements, m_JetElementLocation);

  if(sc == StatusCode::FAILURE) {
    msg(MSG::INFO) << "No JetElements found in TES at " << m_JetElementLocation
                   << endreq;
    return StatusCode::SUCCESS;
  }
         
  // Step over all cells 
  LVL1::CoordToHardware ToHW;
  JECollection::const_iterator it_je ;
  for (it_je = jetElements->begin(); it_je != jetElements->end(); ++it_je) {
    double eta = (*it_je)->eta();
    double phi = (*it_je)->phi();
    LVL1::Coordinate coord(phi, eta);
    int crate  = ToHW.jepCrate(coord);
    int module = ToHW.jepModule(coord);
    int cord   = ToHW.jepCoordinateWord(coord);
    int emEnergy  = (*it_je)->emEnergy();
    int hadEnergy = (*it_je)->hadEnergy();
	  
    if (debug) {
      msg(MSG::VERBOSE) << "JE has coords (eta,phi): " << eta << ", " << phi
                        << " and energies (Em,Had): " << emEnergy << ", "
  		        << hadEnergy << " HW Crate:" << crate
		        << " Module: " << module << " " << cord << endreq;
    }
         
    if (emEnergy > 0) { 
      m_h_je_emeta->Fill(eta, 1.);
      m_histTool->fillJEMPhi(m_h_je_emphi, eta, phi, 1.);
      m_h_je_emenergy->Fill(emEnergy, 1.);
      m_histTool->fillJEMEtaVsPhi(m_h_je_energy_emHitMap, eta, phi, emEnergy);
    }
    if (hadEnergy > 0) { 
      m_h_je_hadeta -> Fill(eta, 1.);
      m_histTool->fillJEMPhi(m_h_je_hadphi, eta, phi, 1.);
      m_h_je_hadenergy->Fill(hadEnergy, 1.);
      m_histTool->fillJEMEtaVsPhi(m_h_je_energy_hadHitMap, eta, phi, hadEnergy);
    }
      
    // number of triggered slice
    m_h_je_triggeredSlice->Fill((*it_je)->peak(),1);
      
    // ----------------- HitMaps per time slice ------------------------------
    const std::vector<int>& emEnergyVec((*it_je)->emEnergyVec());
    const std::vector<int>& hadEnergyVec((*it_je)->hadEnergyVec());
    int slicesEm  = emEnergyVec.size();
    int slicesHad = hadEnergyVec.size();
    for (int i = 0; i < m_SliceNo; i++) {
      if (i < slicesEm && emEnergyVec[i] > 0) {
        m_histTool->fillJEMEtaVsPhi(m_h_je_emHitMap[i], eta, phi , 1.);
      }
      if (i < slicesHad && hadEnergyVec[i] > 0) {
        m_histTool->fillJEMEtaVsPhi(m_h_je_hadHitMap[i], eta, phi, 1.);
      } 
    }

    // ----------------- Error Histos ----------------------------------------
    DataError err((*it_je)->emError());
    DataError haderr((*it_je)->hadError());

    int ypos = crate*16 + module;
    // EM Parity
    if (err.get(DataError::Parity)) {
      m_histTool->fillJEMEtaVsPhi(m_h_je_em_parity, eta, phi);
      m_histTool->fillEventNumber(m_h_JEM_Events, EMParity);
      m_h_JEM_ErrorSummary->Fill(EMParity);
      overview[crate] |= (1 << EMParity);
    }
    // HAD Parity
    if (haderr.get(DataError::Parity)) {
      m_histTool->fillJEMEtaVsPhi(m_h_je_had_parity, eta, phi);
      m_histTool->fillEventNumber(m_h_JEM_Events, HadParity);
      m_h_JEM_ErrorSummary->Fill(HadParity);
      overview[crate] |= (1 << HadParity);
    }
    // PPM Link down: em.
    if (err.get(DataError::LinkDown)) {
      m_histTool->fillJEMEtaVsPhi(m_h_je_em_link, eta, phi);
      m_histTool->fillEventNumber(m_h_JEM_Events, EMLink);
      m_h_JEM_ErrorSummary->Fill(EMLink);
      overview[crate] |= (1 << EMLink);
    }
    // PPM Link down: had.
    if (haderr.get(DataError::LinkDown)) {
      m_histTool->fillJEMEtaVsPhi(m_h_je_had_link, eta, phi);
      m_histTool->fillEventNumber(m_h_JEM_Events, HadLink);
      m_h_JEM_ErrorSummary->Fill(HadLink);
      overview[crate] |= (1 << HadLink);
    }
	  
    //Errors from substatus word from ROD: JEM
    const int status = (err.error() >> LVL1::DataError::GLinkParity) & 0xff;
    if (status) {
      for (int bit = 0; bit < 8; ++bit) {
        if ((status >> bit) & 0x1) m_h_je_error->Fill(bit, ypos);
      }
      m_histTool->fillEventNumber(m_h_JEM_Events, JEMStatus);
      m_h_JEM_ErrorSummary->Fill(JEMStatus);
      overview[crate] |= (1 << JEMStatus);
    }
  }

  // =========================================================================
  // ================= Container: JEM Hits ===================================
  // =========================================================================

  // retrieve JEMHits collection from storegate
  const JEMHitsCollection* JEMHits;
  sc = evtStore()->retrieve(JEMHits, m_JEMHitsLocation);
  if (sc == StatusCode::FAILURE) {
    msg(MSG::INFO) << "No JEMHits found in TES at " << m_JEMHitsLocation
                   << endreq ;
    return StatusCode::SUCCESS;
  }
  
  if (debug) {
    msg(MSG::DEBUG) << "-------------- JEM Hits ---------------" << endreq;
  }
  
  // Step over all cells and process
  JEMHitsCollection::const_iterator it_JEMHits ;
  for (it_JEMHits = JEMHits->begin(); it_JEMHits != JEMHits->end();
                                                              ++it_JEMHits ) {	  
    int crate  = (*it_JEMHits)->crate();
    int module = (*it_JEMHits)->module();
    int xpos   = crate*16 + module;
    bool forward = (*it_JEMHits)->forward();
    unsigned int jetHits = (*it_JEMHits)->JetHits();

    int nBits = (forward) ? 2 : 3;
    m_histTool->fillThresholds(m_h_JEMHits_MainHits, jetHits, 8, nBits);
    m_histTool->fillXVsThresholds(m_h_JEMDAQ_Hits_Map, xpos, jetHits, 8, nBits);
    if (forward) {
      int fwdHits = jetHits >> 16;
      int offset  = (module%8 == 0) ? 8 : 12;
      TH1F_LW* fwdHist = (module%8 == 0) ? m_h_JEMHits_FwdHitsLeft
                                         : m_h_JEMHits_FwdHitsRight;
      m_histTool->fillThresholds(fwdHist, fwdHits, 4, nBits);
      m_histTool->fillXVsThresholds(m_h_JEMDAQ_Hits_Map, xpos, fwdHits, 4,
                                                                nBits, offset);
    }

    if (debug) {
      msg(MSG::DEBUG) << "Crate: "<< crate << "  Module: " << module
 	  << "  JetHits: "
	  << m_histTool->thresholdString(jetHits, (forward) ? 12 : 8, nBits)
	  << endreq;
    }
  }   
  
  // =========================================================================
  // ================= Container: JEM Et Sums ================================
  // =========================================================================

  const JEMEtSumsCollection* JEMEtSums;
  sc = evtStore()->retrieve(JEMEtSums, m_JEMEtSumsLocation);
  if (sc == StatusCode::FAILURE) {
    msg(MSG::INFO) << "No JEMEtSums found in TES at " << m_JEMEtSumsLocation
                   << endreq ;
    return StatusCode::SUCCESS;
  }

  if (debug) {
    msg(MSG::DEBUG) << "-------------- JEM Et Sums ---------------" << endreq;
  }

  // Step over all cells
  JEMEtSumsCollection::const_iterator it_JEMEtSums ;
  LVL1::QuadLinear expand;

  for (it_JEMEtSums = JEMEtSums->begin(); it_JEMEtSums != JEMEtSums->end();
                                                             ++it_JEMEtSums) {	       
    // note: the energy values are compressed -> expand!
    int ex = expand.Expand((*it_JEMEtSums)->Ex());
    int ey = expand.Expand((*it_JEMEtSums)->Ey());
    int et = expand.Expand((*it_JEMEtSums)->Et());

    if (ex != 0) m_h_JEMEtSums_Ex->Fill(ex, 1.); 
    if (ey != 0) m_h_JEMEtSums_Ey->Fill(ey, 1.); 
    if (et != 0) m_h_JEMEtSums_Et->Fill(et, 1.); 

    if (debug) {
      msg(MSG::DEBUG) << " JEMEtSums Crate: " << (*it_JEMEtSums)->crate()
                      << "  Module: "         << (*it_JEMEtSums)->module()
                      << "   Ex: "            <<  ex
	              << "   Ey: "            <<  ey 
	              << "   Et: "            <<  et
		      << "   Et compressed: " << (*it_JEMEtSums)-> Et()
		      << endreq;
    }
  }

  // =========================================================================
  // ================= Container: JEM RoI ====================================
  // =========================================================================

  const JemRoiCollection* JEMRoIs = 0;
  sc = evtStore()->retrieve (JEMRoIs, m_JEMRoILocation);
  if (sc == StatusCode::FAILURE) {
    msg(MSG::INFO) << "No JEM RoIs found in TES at" << m_JEMRoILocation
                   << endreq;
    return StatusCode::SUCCESS;    
  }

  if (debug) {
    msg(MSG::DEBUG) << "-------------- JEM RoIs ---------------" << endreq;
  }

  // Step over all cells
  JemRoiCollection::const_iterator it_JEMRoIs ;

  for (it_JEMRoIs = JEMRoIs->begin(); it_JEMRoIs != JEMRoIs->end();
                                                          ++it_JEMRoIs) {	  
    int crate   = (*it_JEMRoIs)->crate();
    int module  = (*it_JEMRoIs)->jem();
    int forward = (*it_JEMRoIs)->forward();
    int roiHits = (*it_JEMRoIs)->hits();
    LVL1::JEPRoIDecoder decoder;
    LVL1::CoordinateRange coordRange =
                                  decoder.coordinate((*it_JEMRoIs)->roiWord());
    double eta = coordRange.eta();
    double phi = coordRange.phi();
      
    int nHits = 8;
    TH1F_LW* hist = m_h_JEMRoI_MainHits;
    if (forward) {
      nHits = 4;
      hist = (module%8 == 0) ? m_h_JEMRoI_FwdHitsLeft : m_h_JEMRoI_FwdHitsRight;
    }
    m_histTool->fillThresholds(hist, roiHits, nHits, 1);

    for (int thr = 0; thr < nHits; ++thr) {
      int hit = (roiHits >> thr) & 0x1;
      if (hit) {
        TH2F_LW* hist2 = (forward) ? m_h_JEMRoI_FwdThreshPerEtaPhi[thr]
	                           : m_h_JEMRoI_MainThreshPerEtaPhi[thr];
        m_histTool->fillJEMRoIEtaVsPhi(hist2, eta, phi, 1.);
      }
    }

    if (debug) {
      msg(MSG::DEBUG) << "JEMRoI Word: "
                      << MSG::hex << (*it_JEMRoIs)->roiWord() << MSG::dec
	              << "; Crate: "   << crate   << "; JEM: " << module
	              << "; forward: " << forward << "; Hits: "
	              << m_histTool->thresholdString(roiHits, nHits)
	              << endreq;
    }
      
    DataError err((*it_JEMRoIs)->error());

    if (err.get(DataError::Parity)) {
      m_histTool->fillJEMRoIEtaVsPhi(m_h_JEMRoI_error, eta, phi);
      m_histTool->fillEventNumber(m_h_JEM_Events, RoIParity);
      m_h_JEM_ErrorSummary->Fill(RoIParity);
      overview[crate] |= (1 << RoIParity);
    }
    // saturation
    if (err.get(DataError::Overflow)) {
      m_histTool->fillJEMRoIEtaVsPhi(m_h_JEMRoI_sat, eta, phi);
    }
  }

  // Write overview vector to StoreGate
  std::vector<int>* save = new std::vector<int>(overview);
  sc = evtStore()->record(save, "L1CaloJEMErrorVector");
  if (sc != StatusCode::SUCCESS) {
    msg(MSG::ERROR) << "Error recording JEM error vector in TES " << endreq;
    return sc;
  }

  if (debug) {
    msg(MSG::DEBUG) << "--------------------------------------" << endreq;
  }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode JEMMon::procHistograms( bool isEndOfEventsBlock, 
				   bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
  msg(MSG::DEBUG) << "in procHistograms" << endreq ;

  if( isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun ) { }
	
  return StatusCode::SUCCESS;
}
