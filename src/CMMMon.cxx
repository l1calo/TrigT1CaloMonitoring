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

#include <sstream>
#include <vector>
#include <set>

#include "GaudiKernel/IJobOptionsSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "StoreGate/StoreGateSvc.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "SGTools/StlVectorClids.h"

#include <TROOT.h>
#include <TColor.h> 
#include <TCanvas.h> 
#include "TH1.h"
#include "TStyle.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <iostream>
#include <cstdlib>

#include "TrigT1CaloMonitoring/CMMMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"
#include "TrigT1CaloMonitoring/TrigT1CaloMonErrorTool.h"

#include "TrigT1CaloEvent/CMMRoI.h"
#include "TrigT1CaloUtils/QuadLinear.h"
#include "TrigT1CaloUtils/DataError.h"
#include "TrigT1CaloUtils/CrateEnergy.h"

#include "TrigT1Interfaces/TrigT1CaloDefs.h"
#include "TrigT1Interfaces/Coordinate.h"

#include "AthenaMonitoring/AthenaMonManager.h"

namespace LVL1 {
  class CMMRoI;
}


// *********************************************************************
// Public Methods
// *********************************************************************

/*---------------------------------------------------------*/
CMMMon::CMMMon( const std::string & type, const std::string & name,
		const IInterface* parent )
  : ManagedMonitorToolBase( type, name, parent ),
    m_errorTool("TrigT1CaloMonErrorTool")
/*---------------------------------------------------------*/
{
  // This is how you declare the parameters to Gaudi so that
  // they can be over-written via the job options file

  declareProperty( "CMMJetHitsLocation", m_CMMJetHitsLocation =  LVL1::TrigT1CaloDefs::CMMJetHitsLocation) ;
  declareProperty( "CMMEtSumsLocation", m_CMMEtSumsLocation =  LVL1::TrigT1CaloDefs::CMMEtSumsLocation) ;  
  declareProperty( "CMMRoILocation", m_CMMRoILocation =  LVL1::TrigT1CaloDefs::CMMRoILocation) ;
  declareProperty( "MaxEnergyRange", m_MaxEnergyRange = 32768) ;
  declareProperty( "Offline", m_Offline = 1) ;


  declareProperty( "PathInRootFile", m_PathInRootFile="L1Calo/JEM_CMM") ;
  declareProperty( "ErrorPathInRootFile", m_ErrorPathInRootFile="L1Calo/JEM_CMM/Errors/Hardware") ;
  declareProperty( "TypeOfData", m_DataType="") ;
}


/*---------------------------------------------------------*/
CMMMon::~CMMMon()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode CMMMon::initialize()
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );

  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;

  sc = m_errorTool.retrieve();
  if( sc.isFailure() ) {
    mLog << MSG::ERROR << "Unable to locate Tool TrigT1CaloMonErrorTool"
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
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in CMMMon::bookHistograms" << endreq;
  
  /** get a handle of StoreGate for access to the Event Store */
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  
  if (sc.isFailure()) 
    {
      mLog << MSG::ERROR
	   << "Unable to retrieve pointer to StoreGateSvc"
	   << endreq;
      return sc;
    }

  /*
  ManagedMonitorToolBase::LevelOfDetail_t LevelOfDetail=shift;
  if (m_DataType=="Sim") LevelOfDetail = expert;

  MonGroup CMM_DAQ ( this, (m_PathInRootFile+"_DAQ").c_str(), expert, run );
  HistoBooker DAQ_Booker(&CMM_DAQ, &mLog, m_DataType);

  MonGroup CMM_input ( this, (m_PathInRootFile + "_input").c_str(), expert, run );
  HistoBooker input_Booker(&CMM_input, &mLog, m_DataType);
  
  MonGroup CMM_RoI ( this, (m_PathInRootFile + "_RoI").c_str(), LevelOfDetail, run );
  HistoBooker RoI_Booker(&CMM_RoI, &mLog, m_DataType);
  
  MonGroup CMM_transmission ( this, (m_ErrorPathInRootFile ).c_str(), shift, run );
  HistoBooker transmission_Booker(&CMM_transmission, &mLog, "");
  */
  MonGroup CMM_inputThresh ( this, m_PathInRootFile+"/Input/Thresholds", expert, run );
  HistoBooker inputThresh_Booker(&CMM_inputThresh, &mLog, "");

  MonGroup CMM_inputEnergy ( this, m_PathInRootFile+"/Input/EnergySums", expert, run );
  HistoBooker inputEnergy_Booker(&CMM_inputEnergy, &mLog, "");

  MonGroup CMM_jet ( this, m_PathInRootFile+"/Output/Jet", expert, run );
  HistoBooker jet_Booker(&CMM_jet, &mLog, "");

  MonGroup CMM_energy ( this, m_PathInRootFile+"/Output/Energy", expert, run );
  HistoBooker energy_Booker(&CMM_energy, &mLog, "");

  MonGroup CMM_RoI ( this, m_PathInRootFile+"/Output/RoI", shift, run );
  HistoBooker RoI_Booker(&CMM_RoI, &mLog, "");

  MonGroup CMM_transmission ( this, m_ErrorPathInRootFile, shift, run );
  HistoBooker transmission_Booker(&CMM_transmission, &mLog, "");
  
  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
  
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }
  
  if ( isNewEventsBlock|| isNewLumiBlock) { }

  if( isNewRun ) 
    {	
      m_NoEvents=0;
      //----------------------------------  CMM Input data from JEMs -----------------------------
      //m_h_CMMJetHits_JEM_MainHits=inputThresh_Booker.book1F("MainHits_CMM_input", "Main Jet Multiplicity per Threshold  --  CMM input", 8, -0.5, 7.5, "Threshold No.", "N");
      //m_h_CMMJetHits_JEM_FwdHitsRight=inputThresh_Booker.book1F("FwdHitsRight_CMM_input", "Forward Right Jet Multiplicity per Threshold  --  CMM input",4 , -0.5, 3.5, "Threshold No.", "N");
      //m_h_CMMJetHits_JEM_FwdHitsLeft=inputThresh_Booker.book1F("FwdHitsLeft_CMM_input", "Forward Left Jet Multiplicity per Threshold  --  CMM input", 4 , -0.5, 3.5,  "Threshold No.", "N");
      m_h_CMMJetHits_JEM_MainHits=inputThresh_Booker.book1F("cmm_1d_thresh_MainHits", "Main Jet Multiplicity per Threshold  --  CMM input", 8, -0.5, 7.5, "Threshold No.", "N");
      m_h_CMMJetHits_JEM_FwdHitsRight=inputThresh_Booker.book1F("cmm_1d_thresh_FwdHitsRight", "Forward Right Jet Multiplicity per Threshold  --  CMM input",4 , -0.5, 3.5, "Threshold No.", "N");
      m_h_CMMJetHits_JEM_FwdHitsLeft=inputThresh_Booker.book1F("cmm_1d_thresh_FwdHitsLeft", "Forward Left Jet Multiplicity per Threshold  --  CMM input", 4 , -0.5, 3.5,  "Threshold No.", "N");

      // Choose binning to match the encoding
      const int eRange = 256;  //range of encoded value
      const int dRange = 4096; //range of decoded value
      LVL1::QuadLinear expand;
      std::set<int> sorted;
      for (int i = 0; i < eRange; ++i) {
	int val = expand.Expand(i);
        if (val != 0) sorted.insert(val);
      }
      double binedges[eRange+2];
      int nbins = 0;
      std::set<int>::const_iterator iter  = sorted.begin();
      std::set<int>::const_iterator iterE = sorted.end();
      for (; iter != iterE; ++iter) {
        binedges[nbins] = *iter;
	mLog<<MSG::DEBUG<<"binedge["<<nbins<<"] = "<< *iter<<endreq;
	++nbins;
      }
      binedges[0] = 1;
      binedges[nbins] = dRange;
      mLog<<MSG::DEBUG<<"binedge["<<nbins<<"] = "<< dRange<<endreq;
      //m_h_CMMEtSums_JEM_Ex=inputEnergy_Booker.book1F("Ex_CMM_input", "CMM E_{x}^{JEM}  --  CMM input", nbins, 0, dRange, "Ex [GeV]", "N");
      //m_h_CMMEtSums_JEM_Ey=inputEnergy_Booker.book1F("Ey_CMM_input", "CMM E_{y}^{JEM}  --  CMM input", nbins, 0, dRange, "Ey [GeV]", "N");
      //m_h_CMMEtSums_JEM_Et=inputEnergy_Booker.book1F("Et_CMM_input", "CMM E_{t}^{JEM}  --  CMM input", nbins, 0, dRange, "Et [GeV]", "N");
      m_h_CMMEtSums_JEM_Ex=inputEnergy_Booker.book1F("cmm_1d_energy_SubSumsEx", "CMM E_{x}^{JEM}  --  CMM input", nbins, 1, dRange, "Ex [GeV]", "N");
      m_h_CMMEtSums_JEM_Ey=inputEnergy_Booker.book1F("cmm_1d_energy_SubSumsEy", "CMM E_{y}^{JEM}  --  CMM input", nbins, 1, dRange, "Ey [GeV]", "N");
      m_h_CMMEtSums_JEM_Et=inputEnergy_Booker.book1F("cmm_1d_energy_SubSumsEt", "CMM E_{t}^{JEM}  --  CMM input", nbins, 1, dRange, "Et [GeV]", "N");
      m_h_CMMEtSums_JEM_Ex->SetBins(nbins, binedges);
      m_h_CMMEtSums_JEM_Ey->SetBins(nbins, binedges);
      m_h_CMMEtSums_JEM_Et->SetBins(nbins, binedges);
      

      //---------------------------------- CMM output to DAQ -----------------------------
      //m_h_CMMJetHits_MainJets = jet_Booker.book1F("TotalMainHits_CMM_DAQ", "Main Jet Multiplicity per Threshold  --  CMM DAQ", 8, -0.5, 7.5, "Threshold No.", "N");
      //m_h_CMMJetHits_FwdJetsRight = jet_Booker.book1F("TotalFwdHitsRight_CMM_DAQ", "Forward Right Jet Multiplicity per Threshold  --  CMM DAQ", 4 , -0.5, 3.5, "Threshold No.", "N");
      //m_h_CMMJetHits_FwdJetsLeft = jet_Booker.book1F("TotalFwdHitsLeft_CMM_DAQ", "Forward Left Jet Multiplicity per Threshold  --  CMM DAQ", 4 , -0.5, 3.5,  "Threshold No.", "N");
      //m_h_CMMJetHits_EtMap = jet_Booker.book1F("JetEtHits_CMM_DAQ", "JetEt Multiplicity per Threshold  --  CMM DAQ", 4 ,-0.5, 3.5, "Threshold No.", "N");
      //m_h_CMMEtSums_MissingEtMap = energy_Booker.book1F("MissingEtHits_CMM_DAQ", "MissingEt Multiplicity per Threshold  --  CMM DAQ", 8, -0.5, 7.5, "Threshold No.", "N");
      //m_h_CMMEtSums_SumEtMap = energy_Booker.book1F("SumEtHits_CMM_DAQ", "SumEt Multiplicity per Threshold  --  CMM DAQ", 4, -0.5, 3.5, "Threshold No.", "N");
      m_h_CMMJetHits_MainJets = jet_Booker.book1F("cmm_1d_thresh_TotalMainHits", "Main Jet Multiplicity per Threshold  --  CMM DAQ", 8, -0.5, 7.5, "Threshold No.", "N");
      m_h_CMMJetHits_FwdJetsRight = jet_Booker.book1F("cmm_1d_thresh_TotalFwdHitsRight", "Forward Right Jet Multiplicity per Threshold  --  CMM DAQ", 4 , -0.5, 3.5, "Threshold No.", "N");
      m_h_CMMJetHits_FwdJetsLeft = jet_Booker.book1F("cmm_1d_thresh_TotalFwdHitsLeft", "Forward Left Jet Multiplicity per Threshold  --  CMM DAQ", 4 , -0.5, 3.5,  "Threshold No.", "N");
      m_h_CMMJetHits_EtMap = jet_Booker.book1F("cmm_1d_thresh_JetEtHits", "JetEt Multiplicity per Threshold  --  CMM DAQ", 4 ,-0.5, 3.5, "Threshold No.", "N");
      m_h_CMMEtSums_MissingEtMap = energy_Booker.book1F("cmm_1d_energy_MissingEtHits", "MissingEt Multiplicity per Threshold  --  CMM DAQ", 8, -0.5, 7.5, "Threshold No.", "N");
      m_h_CMMEtSums_SumEtMap = energy_Booker.book1F("cmm_1d_energy_SumEtHits", "SumEt Multiplicity per Threshold  --  CMM DAQ", 4, -0.5, 3.5, "Threshold No.", "N");

      for (int i = nbins; i >= 0; --i) binedges[i+1] = binedges[i]*8;
      nbins++;
      //m_h_CMMEtSums_Ex = energy_Booker.book1F("Ex_CMM_DAQ" "E_{x}^{CMM}  --  CMM DAQ", nbins, 0, binedges[nbins], "Ex [GeV]", "N");
      //m_h_CMMEtSums_Ey = energy_Booker.book1F("Ey_CMM_DAQ", "E_{y}^{CMM}  --  CMM DAQ", nbins, 0, binedges[nbins], "Ey [GeV]", "N");
      //m_h_CMMEtSums_Et = energy_Booker.book1F("Et_CMM_DAQ", "SumE_{t}^{CMM}  --  CMM DAQ", nbins, 0, binedges[nbins], "Et [GeV]", "N");
      m_h_CMMEtSums_Ex = energy_Booker.book1F("cmm_1d_energy_TotalEx", "E_{x}^{CMM}  --  CMM DAQ", nbins, 1, binedges[nbins], "Ex [GeV]", "N");
      m_h_CMMEtSums_Ey = energy_Booker.book1F("cmm_1d_energy_TotalEy", "E_{y}^{CMM}  --  CMM DAQ", nbins, 1, binedges[nbins], "Ey [GeV]", "N");
      m_h_CMMEtSums_Et = energy_Booker.book1F("cmm_1d_energy_TotalEt", "SumE_{t}^{CMM}  --  CMM DAQ", nbins, 1, binedges[nbins], "Et [GeV]", "N");
      m_h_CMMEtSums_Ex->SetBins(nbins, binedges);
      m_h_CMMEtSums_Ey->SetBins(nbins, binedges);
      m_h_CMMEtSums_Et->SetBins(nbins, binedges);


      //---------------------------------- CMM output to RoI -----------------------------
      //m_h_CMMRoI_JetEtHits =RoI_Booker.book1F("JetEtHits_CMM_RoI","JetEt Multiplicity per Threshold  --  CMM RoI", 4, -0.5,3.5,"Threshold No.","N");
      //m_h_CMMRoI_MissingEtHits =RoI_Booker.book1F("MissingEtHits_CMM_RoI","MissingEt Multiplicity per Threshold  --  CMM RoI", 8, -0.5,7.5,"Threshold No.","N");
      //m_h_CMMRoI_SumEtHits =RoI_Booker.book1F("SumEtHits_CMM_RoI","SumEt Multiplicity per Threshold  --  CMM RoI", 4, -0.5,3.5,"Threshold No.","N");

      //m_h_CMMRoI_Ex = RoI_Booker.book1F("Ex_CMM_RoI", "E_{x}^{CMM}  --  CMM RoI", nbins, 0, binedges[nbins], "Ex [GeV]", "N");
      //m_h_CMMRoI_Ey = RoI_Booker.book1F("Ey_CMM_RoI", "E_{y}^{CMM}  --  CMM RoI", nbins, 0, binedges[nbins], "Ey [GeV]", "N");
      //m_h_CMMRoI_Et = RoI_Booker.book1F("Et_CMM_RoI", "SumE_{t}^{CMM}  --  CMM RoI", nbins, 0, binedges[nbins], "Et [GeV]", "N");
      m_h_CMMRoI_JetEtHits =RoI_Booker.book1F("cmm_1d_roi_JetEtHits","JetEt Multiplicity per Threshold  --  CMM RoI", 4, -0.5,3.5,"Threshold No.","N");
      m_h_CMMRoI_MissingEtHits =RoI_Booker.book1F("cmm_1d_roi_MissingEtHits","MissingEt Multiplicity per Threshold  --  CMM RoI", 8, -0.5,7.5,"Threshold No.","N");
      m_h_CMMRoI_SumEtHits =RoI_Booker.book1F("cmm_1d_roi_SumEtHits","SumEt Multiplicity per Threshold  --  CMM RoI", 4, -0.5,3.5,"Threshold No.","N");

      m_h_CMMRoI_Ex = RoI_Booker.book1F("cmm_1d_roi_Ex", "E_{x}^{CMM}  --  CMM RoI", nbins, 1, binedges[nbins], "Ex [GeV]", "N");
      m_h_CMMRoI_Ey = RoI_Booker.book1F("cmm_1d_roi_Ey", "E_{y}^{CMM}  --  CMM RoI", nbins, 1, binedges[nbins], "Ey [GeV]", "N");
      m_h_CMMRoI_Et = RoI_Booker.book1F("cmm_1d_roi_Et", "SumE_{t}^{CMM}  --  CMM RoI", nbins, 1, binedges[nbins], "Et [GeV]", "N");
      m_h_CMMRoI_Ex->SetBins(nbins, binedges);
      m_h_CMMRoI_Ey->SetBins(nbins, binedges);
      m_h_CMMRoI_Et->SetBins(nbins, binedges);

       


      if (m_DataType=="BS")
	{
	  //---------------------------------- S-Link errors -----------------------------
	  //m_h_CMMJet_error=transmission_Booker.book2F("CMMJet_errors", "Errors from CMM Jet SubStatus Word",9,0.5,9.5,37,0.5,37.5,"","");
	  m_h_CMMJet_error=transmission_Booker.book2F("cmm_2d_thresh_Status", "Errors from CMM Jet SubStatus Word",9,0.5,9.5,37,0.5,37.5,"","");
	  m_h_CMMJet_error->SetStats(kFALSE);
	  
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(1, "Parity");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(3, "GLinkParity");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(4, "GLinkProtocol");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(5, "BCNMismatch");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(6, "FIFOOverflow");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(7, "ModuleError");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(8, "GLinkDown");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(9, "GLinkTimeout");
	  

	  //m_h_CMMEnergy_error=transmission_Booker.book2F("CMMEnergy_errors", "Errors from CMM Energy SubStatus Word",9,0.5,9.5,37,0.5,37.5,"","");
	  m_h_CMMEnergy_error=transmission_Booker.book2F("cmm_2d_energy_Status", "Errors from CMM Energy SubStatus Word",9,0.5,9.5,37,0.5,37.5,"","");
	  m_h_CMMEnergy_error->SetStats(kFALSE);
	 
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(1, "Parity");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(3, "GLinkParity");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(4, "GLinkProtocol");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(5, "BCNMismatch");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(6, "FIFOOverflow");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(7, "ModuleError");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(8, "GLinkDown");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(9, "GLinkTimeout");
	  

 	  std::string name;
	  std::stringstream buffer;
     
	  for (int i = 0; i < 16; i++)
	    {
	      buffer.str("");
	      buffer<<i;
	      name = "JEM " + buffer.str();
	      m_h_CMMJet_error->GetYaxis()->SetBinLabel((i+1), name.c_str());
	      m_h_CMMEnergy_error->GetYaxis()->SetBinLabel((i+1), name.c_str());
	      
	      buffer.str("");
	      buffer<<i;
	      name = "JEM " + buffer.str();
	      m_h_CMMJet_error->GetYaxis()->SetBinLabel((i+1+19), name.c_str());
	      m_h_CMMEnergy_error->GetYaxis()->SetBinLabel((i+1+19), name.c_str());
	    }
	  m_h_CMMJet_error->GetYaxis()->SetBinLabel(17, "CMM 1 ");
	  m_h_CMMJet_error->GetYaxis()->SetBinLabel(18, "Crate 0: ");
	  m_h_CMMJet_error->GetYaxis()->SetBinLabel(36, "CMM 1 ");
	  m_h_CMMJet_error->GetYaxis()->SetBinLabel(37, "Crate 1: ");

	  m_h_CMMEnergy_error->GetYaxis()->SetBinLabel(17, "CMM 0");
	  m_h_CMMEnergy_error->GetYaxis()->SetBinLabel(18, "Crate 0: ");
	  m_h_CMMEnergy_error->GetYaxis()->SetBinLabel(36, "CMM 0");
	  m_h_CMMEnergy_error->GetYaxis()->SetBinLabel(37, "Crate 1: ");
	

	  //m_h_CMMRoI_error=transmission_Booker.book1F("CMMRoI_errors", "CMM RoI Parity and Overflow",8,0.5,8.5,"");
	  m_h_CMMRoI_error=transmission_Booker.book1F("cmm_1d_roi_Parity", "CMM RoI Parity and Overflow",8,0.5,8.5,"");
	  m_h_CMMRoI_error->SetStats(kFALSE);

	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(1, "Parity (Ex)");
	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(2, "Parity (Ey, #SigmaEtMap)");
	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(3, "Parity (Et,Et_{Miss}Map)");
	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(4, "Parity (JetEtMap)");
	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(5, "Comp of #slice");
	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(6, "Overflow (Ex)");
	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(7, "Overflow (Ey)");
	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(8, "Overflow (Et)");
	  
	  //m_h_TriggeredSlice=transmission_Booker.book1F("TriggeredSlice", "Comparison of the triggered slice number",3,0.5,3.5,"Difference","N");
	  m_h_TriggeredSlice=transmission_Booker.book1F("cmm_1d_TriggeredSlices", "Comparison of CMM Jet and Energy triggered slice numbers",3,0.5,3.5,"Difference","N");
  
	 //Error Summary for all CMMs in system
	  //m_h_CMM_ErrorSummary = transmission_Booker.book1F("CMM_ErrorSummary", "Error Summary of CMM Jet, Energy and RoI path",
	  m_h_CMM_ErrorSummary = transmission_Booker.book1F("cmm_1d_ErrorSummary", "Error Summary of CMM Jet, Energy and RoI path",
	  3,0.5,3.5,"","Entries");	 
	  m_h_CMM_ErrorSummary->SetStats(kFALSE);
	  m_h_CMM_ErrorSummary->GetXaxis()->SetBinLabel(1,"CMM Status");
	  m_h_CMM_ErrorSummary->GetXaxis()->SetBinLabel(2,"Parity flags");
	  m_h_CMM_ErrorSummary->GetXaxis()->SetBinLabel(3,"Other");
	}
    }
  
  return StatusCode::SUCCESS;
}


/*---------------------------------------------------------*/
StatusCode CMMMon::fillHistograms()
  /*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );

  // Skip events believed to be corrupt

  if (m_errorTool->corrupt()) {
    mLog << MSG::DEBUG << "Skipping corrupt event" << endreq;
    return StatusCode::SUCCESS;
  }

  Helper Help;
  m_NoEvents++;

  // Error vector for global overview
  std::vector<int> overview(2);

  // triggered slice numbers
  int j_num_slice = -1;
  int e_num_slice = -1;

  using LVL1::DataError;

  // =============================================================================================
  // ================= Container: CMM Jet Hits ===================================================
  // =============================================================================================
  // retrieve CMM Jet Hits from Storegate
  const CMMJetHitsCollection* CMMJetHits;
  StatusCode sc = m_storeGate->retrieve(CMMJetHits, m_CMMJetHitsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No CMM JetHits found in TES at "  << m_CMMJetHitsLocation << endreq ;
      return StatusCode::SUCCESS;
    }  
  mLog<<MSG::DEBUG<<"--------------  "<< m_DataType<<" CMM Jet Hits ---------------"<<endreq;  
  CMMJetHitsCollection::const_iterator it_CMMJetHits ;
  
  // Step over all cells
  for( it_CMMJetHits  = CMMJetHits ->begin(); it_CMMJetHits < CMMJetHits -> end(); ++it_CMMJetHits )
    {	  
      //put CMM Jet Hit into string (for further processing)
      std::string CMMHit = Help.Binary((*it_CMMJetHits)->Hits(),24);
      
      mLog<<MSG::DEBUG<<"CMMJetHits Crate: " << (*it_CMMJetHits)-> crate()<< " dataID: "<< (*it_CMMJetHits)-> dataID() <<"   Hits: "<< (*it_CMMJetHits)-> Hits()
	  << " Hits(binary): " << CMMHit <<endreq;
      
      // ------------------------------------------------------------------------------------------
      // ----------------- Histos with distribution of JEM Hit Multiplicities ---------------------
      // ------------------------------------------------------------------------------------------
      //input data from JEMs have dataID 0..15
      if ((*it_CMMJetHits)-> dataID()<16) 
	{
	  //Fwd Hits left have dataID 0 or 8
	  if (((*it_CMMJetHits)-> dataID()==0)or((*it_CMMJetHits)-> dataID()==8) ) 
	    {
	      Help.FillHitsHisto(m_h_CMMJetHits_JEM_FwdHitsLeft, CMMHit, 0, 4, 8, 2, &mLog);
	      Help.FillHitsHisto(m_h_CMMJetHits_JEM_MainHits, CMMHit, 0, 8, 0, 2, &mLog);
	    }
	  else
	    {
	      //Fwd Hits right have dataID 7 or 15
	      if (((*it_CMMJetHits)-> dataID()==7)or((*it_CMMJetHits)-> dataID()==15) ) 
		{
		  Help.FillHitsHisto(m_h_CMMJetHits_JEM_FwdHitsRight, CMMHit, 0, 4, 8, 2, &mLog);
		  Help.FillHitsHisto(m_h_CMMJetHits_JEM_MainHits, CMMHit, 0, 8, 0, 2, &mLog);
		}
	      //Main Hits for all other modules
	      else 
		{
		  Help.FillHitsHisto(m_h_CMMJetHits_JEM_MainHits, CMMHit, 0, 8, 0, 3, &mLog);
		}
	    }
	}
      // ------------------------------------------------------------------------------------------
      // ----------------- Histos with SubStatus Word errors -----------------------------
      // ------------------------------------------------------------------------------------------
       
      //only for Bytestream data
      if (m_DataType=="BS")
	{
	  j_num_slice = (*it_CMMJetHits)-> peak();

	  DataError err((*it_CMMJetHits)-> Error());
	  //input data from JEMs have dataID 0..15   ---  fill only parity errors
	  
	  int crate = (*it_CMMJetHits)->crate();
	  int module = (*it_CMMJetHits)-> dataID();
	  
	  
	  //Error summary plots
	  //substatus word
	  if (err.get(DataError::GLinkParity) || err.get(DataError::GLinkProtocol) || err.get(DataError::BCNMismatch) ||
	      err.get(DataError::FIFOOverflow) || err.get(DataError::ModuleError) || err.get(DataError::GLinkDown) ||
	      err.get(DataError::GLinkTimeout))
	     {
	       m_h_CMM_ErrorSummary->Fill(1,1);
	       overview[crate] |= 1;
	     }
	  	  
	  int ypos = crate*19 + 1 + ((module<16) ? module : 16);
	  if (module<16 || module == LVL1::CMMJetHits::REMOTE_MAIN || module == LVL1::CMMJetHits::REMOTE_FORWARD)
	    {
	      // Parity
	      if (err.get(DataError::Parity))
	        {
	          m_h_CMMJet_error->Fill(1,ypos);
	          m_h_CMM_ErrorSummary->Fill(2);
	          overview[crate] |= (1 << 1);
	        }
	     
	    }
	  
	  // set L1CaloSubStatus for both Crate and System CMM
          ypos = crate*19 + 1 + 16;			      
	  if (err.get(DataError::GLinkParity))   m_h_CMMJet_error->Fill(3,ypos);
	  if (err.get(DataError::GLinkProtocol)) m_h_CMMJet_error->Fill(4,ypos);
	  if (err.get(DataError::BCNMismatch))   m_h_CMMJet_error->Fill(5,ypos);
	  if (err.get(DataError::FIFOOverflow))  m_h_CMMJet_error->Fill(6,ypos);
	  if (err.get(DataError::ModuleError))   m_h_CMMJet_error->Fill(7,ypos);
	  if (err.get(DataError::GLinkDown))     m_h_CMMJet_error->Fill(8,ypos);
	  if (err.get(DataError::GLinkTimeout))  m_h_CMMJet_error->Fill(9,ypos);
	      	      
	}
      
      // ------------------------------------------------------------------------------------------
      // ----------------- Histos with distribution of CMM Hit Multiplicities ---------------------
      // ------------------------------------------------------------------------------------------
      //main total jets 
      if ((*it_CMMJetHits)-> dataID() == LVL1::CMMJetHits::TOTAL_MAIN)  
	{
	  Help.FillHitsHisto(m_h_CMMJetHits_MainJets, CMMHit, 0, 8, 0, 3, &mLog);
	  mLog<<MSG::DEBUG<<"Total Jet Hits: " << CMMHit <<endreq;

	}
      
      //fwd total jets 
      if ((*it_CMMJetHits)-> dataID() == LVL1::CMMJetHits::TOTAL_FORWARD)  
	{
	  CMMHit= Help.Binary((*it_CMMJetHits)->Hits(),16); //total fwd jets only 16 bit long!
	  mLog<<MSG::DEBUG<<"Right|Left Total Jets Hits: " << CMMHit <<endreq;
	  
	  Help.FillHitsHisto(m_h_CMMJetHits_FwdJetsLeft, CMMHit, 0, 4, 0, 2, &mLog);
	  Help.FillHitsHisto(m_h_CMMJetHits_FwdJetsRight, CMMHit, 0, 4, 4, 2, &mLog);
	}
      
      //JetEtSum Hitmap
      if ((*it_CMMJetHits)-> dataID() == LVL1::CMMJetHits::ET_MAP)  
	{
	  CMMHit= Help.Binary((*it_CMMJetHits)->Hits(),4);
	  mLog<<MSG::DEBUG<<"JetEt Hits: " << CMMHit <<endreq;
	  
	  Help.FillHitsHisto(m_h_CMMJetHits_EtMap, CMMHit, 0, 4, 0, 1, &mLog);
	}
    }  
  
  // =============================================================================================
  // ================= Container: CMM Et Sums ====================================================
  // =============================================================================================
  LVL1::QuadLinear expand;
  int ex;
  int ey;
  int et;
  
  // retrieve CMM Et Sums from Storegate
  const CMMEtSumsCollection* CMMEtSums;
  sc = m_storeGate->retrieve(CMMEtSums, m_CMMEtSumsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No CMMEtSums found in TES at " << m_CMMEtSumsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  mLog<<MSG::DEBUG<<"--------------  "<< m_DataType<<" CMM Et Sums ---------------"<<endreq;
  
  // Step over all cells 
  CMMEtSumsCollection::const_iterator it_CMMEtSums ;
  for( it_CMMEtSums  = CMMEtSums ->begin(); it_CMMEtSums < CMMEtSums -> end(); ++it_CMMEtSums )
    {	  
     
      // ------------------------------------------------------------------------------------------
      // ----------------- Histos with distribution of JEM Energies -------------------------------
      // ------------------------------------------------------------------------------------------
      // JEM energy sums, dataID < 16
      if  ((*it_CMMEtSums)-> dataID()<16) 
	{
	  // note: JEM energies are compressed -> use QuadLinear to expand!
	  ex=expand.Expand( (*it_CMMEtSums)-> Ex());
	  ey=expand.Expand( (*it_CMMEtSums)-> Ey());
	  et=expand.Expand( (*it_CMMEtSums)-> Et() );
	  
	  if (ex>0) m_h_CMMEtSums_JEM_Ex -> Fill( ex, 1.);
	  if (ey>0) m_h_CMMEtSums_JEM_Ey -> Fill( ey, 1.);
	  if (et>0) m_h_CMMEtSums_JEM_Et -> Fill( et, 1.);
	}
      
      // -----------------------------------------------------------------------------------------
      // ----------------- Histos with distribution of total Energy per system -------------------
      // ------------------------------------------------------------------------------------------
      // total energy sums
      if  ((*it_CMMEtSums)-> dataID()==LVL1::CMMEtSums::TOTAL and (*it_CMMEtSums)->crate()==1) 
	{
	  // Use CrateEnergy object to decode 15-bit twos-complement format
	  LVL1::CrateEnergy cen((*it_CMMEtSums)->crate(), (*it_CMMEtSums)->Et(), (*it_CMMEtSums)->Ex(), (*it_CMMEtSums)->Ey(),
	                        ((*it_CMMEtSums)->EtError()&0x1), ((*it_CMMEtSums)->ExError()&0x1), ((*it_CMMEtSums)->EyError()&0x1));
	  int Ex = std::abs(cen.ex());
	  int Ey = std::abs(cen.ey());

	  if (Ex>0 && !cen.exOverflow()) m_h_CMMEtSums_Ex -> Fill( Ex, 1.);
	  if (Ey>0 && !cen.eyOverflow()) m_h_CMMEtSums_Ey -> Fill( Ey, 1.);
	  if ((*it_CMMEtSums)-> Et()>0 && !cen.etOverflow()) m_h_CMMEtSums_Et -> Fill( (*it_CMMEtSums)-> Et(), 1.);
	  mLog<<MSG::DEBUG<<"       Ex: "<<Ex<<"; Ey: "<<Ey<<"; Et "<<(*it_CMMEtSums)->Et()<<endreq;
	  mLog<<MSG::DEBUG<<"signed Ex: "<<(*it_CMMEtSums)->Ex()<<"; Ey: "<<(*it_CMMEtSums)->Ey()<<"; Et "<<(*it_CMMEtSums)->Et()<<endreq;

	}
      
      //MissingEt Hitmap
      if ((*it_CMMEtSums)-> dataID() == LVL1::CMMEtSums::MISSING_ET_MAP and (*it_CMMEtSums)->crate()==1)  
	{
	  std::string CMMHit= Help.Binary((*it_CMMEtSums)->Et(),8);
	  mLog<<MSG::DEBUG<<"MissingEt Hits: " << CMMHit <<endreq;
	  
	  Help.FillHitsHisto(m_h_CMMEtSums_MissingEtMap, CMMHit, 0, 8, 0, 1, &mLog);
	}
      
      //SumEt Hitmap
      if ((*it_CMMEtSums)-> dataID() == LVL1::CMMEtSums::SUM_ET_MAP and (*it_CMMEtSums)->crate()==1)  
	{
	  std::string CMMHit= Help.Binary((*it_CMMEtSums)->Et(),4);
	  mLog<<MSG::DEBUG<<"SumEt Hits: " << CMMHit <<endreq;
	  
	  Help.FillHitsHisto(m_h_CMMEtSums_SumEtMap, CMMHit, 0, 4, 0, 1, &mLog);
	}
      
      
      //only for Bytestream data
      if (m_DataType=="BS")
	{
	   e_num_slice = (*it_CMMEtSums)-> peak();
	   if (j_num_slice >= 0) m_h_TriggeredSlice->Fill(std::abs(e_num_slice - j_num_slice));

	  // ------------------------------------------------------------------------------------------
	  // ----------------- Histos with SubStatus Word errors -----------------------------
	  // ------------------------------------------------------------------------------------------
	  DataError eterr((*it_CMMEtSums)-> EtError());
	  
	  int crate = (*it_CMMEtSums)->crate();
	  int module = (*it_CMMEtSums)-> dataID();
	  
	  //Error summary plots
	  //substatus word
	  if (eterr.get(DataError::GLinkParity) || eterr.get(DataError::GLinkProtocol) || eterr.get(DataError::BCNMismatch) ||
	      eterr.get(DataError::FIFOOverflow) || eterr.get(DataError::ModuleError) || eterr.get(DataError::GLinkDown) ||
	      eterr.get(DataError::GLinkTimeout))
	     {
	       m_h_CMM_ErrorSummary->Fill(1,1);
	       overview[crate] |= 1;
	     }
	  	  
	  int ypos = crate*19 + 1 + ((module<16) ? module : 16);
	  if (module<16 || module == LVL1::CMMEtSums::REMOTE)
	    {
	      // Parity
	      DataError exerr((*it_CMMEtSums)-> ExError());
	      DataError eyerr((*it_CMMEtSums)-> EyError());
	      if (eterr.get(DataError::Parity) || exerr.get(DataError::Parity) || eyerr.get(DataError::Parity))
	        {
	          m_h_CMMEnergy_error->Fill(1,ypos);
	          m_h_CMM_ErrorSummary->Fill(2);
	          overview[crate] |= (1 << 1);
	        }
	     
	    }
	  
	  // set L1CaloSubStatus for both Crate and System CMM
          ypos = crate*19 + 1 + 16;			      
	  if (eterr.get(DataError::GLinkParity))   m_h_CMMEnergy_error->Fill(3,ypos);
	  if (eterr.get(DataError::GLinkProtocol)) m_h_CMMEnergy_error->Fill(4,ypos);
	  if (eterr.get(DataError::BCNMismatch))   m_h_CMMEnergy_error->Fill(5,ypos);
	  if (eterr.get(DataError::FIFOOverflow))  m_h_CMMEnergy_error->Fill(6,ypos);
	  if (eterr.get(DataError::ModuleError))   m_h_CMMEnergy_error->Fill(7,ypos);
	  if (eterr.get(DataError::GLinkDown))     m_h_CMMEnergy_error->Fill(8,ypos);
	  if (eterr.get(DataError::GLinkTimeout))  m_h_CMMEnergy_error->Fill(9,ypos);
	     
	}
    }

  // =============================================================================================
  // ================= Container: CMM RoI ========================================================
  // =============================================================================================
  
  // retrieve RoI information from Storegate
  const LVL1::CMMRoI* CR = 0;
  sc = m_storeGate->retrieve (CR, m_CMMRoILocation);
  if (sc==StatusCode::FAILURE)
    {
      mLog <<MSG::INFO<<"No CMM RoI found in TES at "<< m_CMMRoILocation<<endreq;
      return StatusCode::SUCCESS;    
    }

  mLog<<MSG::DEBUG<<"-------------- "<< m_DataType<<" CMM RoI ---------------"<<endreq;

  // ------------------------------------------------------------------------------------------
  // ----------------- Histos filled with CMM RoI information ---------------------------------
  // ------------------------------------------------------------------------------------------

  mLog<<MSG::DEBUG<<"JetEtHits: "<<Help.Binary((CR)->jetEtHits(),4)<<"; SumEtHits: "<<Help.Binary((CR)->sumEtHits(),4)<<"; MissingEtHits: "<<Help.Binary((CR)->missingEtHits(),8)<<endreq;

  // Jet Et Hits
  std::string CMMRoIHit = Help.Binary((CR)-> jetEtHits(),4);
  Help.FillHitsHisto(m_h_CMMRoI_JetEtHits, CMMRoIHit, 0, 4, 0, 1, &mLog);

  // Sum Et Hits
  CMMRoIHit = Help.Binary((CR)-> sumEtHits(),4);
  Help.FillHitsHisto(m_h_CMMRoI_SumEtHits, CMMRoIHit, 0, 4, 0, 1, &mLog);

  // Missing Et Hits
  CMMRoIHit = Help.Binary((CR)-> missingEtHits(),8);
  Help.FillHitsHisto(m_h_CMMRoI_MissingEtHits , CMMRoIHit, 0, 8, 0, 1, &mLog);
 

  // Use CrateEnergy object to decode 15-bit twos-complement format
  LVL1::CrateEnergy cen(1, (CR)->et(), (CR)->ex(), (CR)->ey(),
                        ((CR)->etError()&0x1), ((CR)->exError()&0x1), ((CR)->eyError()&0x1));
  int Ex = std::abs(cen.ex());
  int Ey = std::abs(cen.ey());
  
  mLog<<MSG::DEBUG<<"       Ex: "<<Ex<<"; Ey: "<<Ey<<"; Et "<<(CR)->et()<<endreq;
  mLog<<MSG::DEBUG<<"signed Ex: "<<(CR)-> ex()<<"; Ey: "<<(CR)-> ey()<<"; Et "<<(CR)-> et()<<endreq;

  if(Ex>0 && !cen.exOverflow()) m_h_CMMRoI_Ex -> Fill( Ex,1);
  if(Ey>0 && !cen.eyOverflow()) m_h_CMMRoI_Ey -> Fill( Ey,1);
  if((CR)->et()>0 && !cen.etOverflow()) m_h_CMMRoI_Et -> Fill( (CR)->et(),1);
  
 
  mLog<<MSG::DEBUG<<"CMM Slice numbers: "<<"Jet: "<<j_num_slice<<" Energy: "<<e_num_slice<<endreq;

  
  // errors
  if (m_DataType=="BS")
    {
      DataError exerr((CR)-> exError());
      DataError eyerr((CR)-> eyError());
      DataError eterr((CR)-> etError());
      DataError jetEterr((CR)-> jetEtError());

      // Parity (Ex)
      if (exerr.get(DataError::Parity)) m_h_CMMRoI_error->Fill(1);
      // Parity (Ey,SumEtMap)
      if (eyerr.get(DataError::Parity)) m_h_CMMRoI_error->Fill(2);
      // Parity (Et,MissingEtMap)
      if (eterr.get(DataError::Parity)) m_h_CMMRoI_error->Fill(3);
      // Parity (JetEtMap)
      if (jetEterr.get(DataError::Parity)) m_h_CMMRoI_error->Fill(4);
      
      //----------------Comparison on slice number-----------
      //if ((j_num_slice - e_num_slice)==1) m_h_CMMRoI_error->Fill(5,1);
      if (j_num_slice >= 0 && e_num_slice >= 0 && j_num_slice != e_num_slice) m_h_CMMRoI_error->Fill(5,1);
      //-----------------------------------------------------
      // Overflow (Ex)
      if (exerr.get(DataError::Overflow)) m_h_CMMRoI_error->Fill(6);
      // Overflow (Ey)
      if (eyerr.get(DataError::Overflow)) m_h_CMMRoI_error->Fill(7);
      // Overflow (Et)
      if (eterr.get(DataError::Overflow)) m_h_CMMRoI_error->Fill(8);
      
      //Error summary plots
     //parity
     if (exerr.get(DataError::Parity) || eyerr.get(DataError::Parity) || eterr.get(DataError::Parity) || jetEterr.get(DataError::Parity))
     {
       m_h_CMM_ErrorSummary->Fill(2,1);
       overview[1] |= 0x2;
       
     }
           
    }

  // Write overview vector to StoreGate
  std::vector<int>* save = new std::vector<int>(overview);
  sc = m_storeGate->record(save, "L1CaloJEMCMMErrorVector");
  if (sc != StatusCode::SUCCESS)
    {
      mLog << MSG::ERROR << "Error recording JEM CMM error vector in TES "
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
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in procHistograms" << endreq ;

  if( isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun ) 
    {  
    }
	
  /*
  if(m_Offline==1)
    {
      if (m_DataType=="BS")
	{
	  if( isEndOfRun ) { 
	    std::stringstream buffer;
	    buffer.str("");
	    buffer<<m_NoEvents;
	    std::string title;
	    
	    title = m_h_CMMJet_error-> GetTitle();
	    title=title + " | #events: " + buffer.str();
	    m_h_CMMJet_error->SetTitle(title.c_str());
	    
	    title = m_h_CMMEnergy_error-> GetTitle();
	    title=title + " | #events: " + buffer.str();
	    m_h_CMMEnergy_error->SetTitle(title.c_str());
	    
	    title = m_h_CMMRoI_error-> GetTitle();
	    title=title + " | #events: " + buffer.str();
	    m_h_CMMRoI_error->SetTitle(title.c_str());
	  }
	}
    }
  */
  return StatusCode::SUCCESS;
}
