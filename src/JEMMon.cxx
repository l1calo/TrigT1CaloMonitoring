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

#include <sstream>

#include "GaudiKernel/IJobOptionsSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "StoreGate/StoreGateSvc.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include <TROOT.h>
#include <TColor.h> 
#include <TCanvas.h> 
#include "TH1.h"
#include "TStyle.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <iostream>

#include "TrigT1CaloMonitoring/JEMMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

#include "TrigT1Calo/JEMRoI.h"
#include "TrigT1Calo/QuadLinear.h"
#include "TrigT1Calo/DataError.h"
#include "TrigT1Calo/CoordToHardware.h"
#include "TrigT1Interfaces/Coordinate.h"

#include "TrigT1Interfaces/JEPRoIDecoder.h"
#include "TrigT1Interfaces/TrigT1CaloDefs.h"
#include "TrigT1Interfaces/Coordinate.h"

#include "AthenaMonitoring/AthenaMonManager.h"

/*---------------------------------------------------------*/
JEMMon::JEMMon( const std::string & type, const std::string & name,
		const IInterface* parent )
  : ManagedMonitorToolBase( type, name, parent )
/*---------------------------------------------------------*/
{
  // This is how you declare the parameters to Gaudi so that
  // they can be over-written via the job options file

  declareProperty( "JetElementLocation", m_JetElementLocation = LVL1::TrigT1CaloDefs::JetElementLocation); 
  declareProperty( "JEMHitsLocation", m_JEMHitsLocation =  LVL1::TrigT1CaloDefs::JEMHitsLocation) ;
  declareProperty( "JEMEtSumsLocation", m_JEMEtSumsLocation=   LVL1::TrigT1CaloDefs::JEMEtSumsLocation) ;
  declareProperty( "JEMRoILocation", m_JEMRoILocation =  LVL1::TrigT1CaloDefs::JEMRoILocation) ;
  declareProperty( "NumberOfSlices", m_SliceNo = 5);
  declareProperty( "MaxEnergyRange", m_MaxEnergyRange = 50) ;
  declareProperty( "Offline", m_Offline = 1) ;

  declareProperty( "PathInRootFile", m_PathInRootFile="Stats/JEM") ;
  declareProperty( "ErrorPathInRootFile", m_ErrorPathInRootFile="Stats/L1Calo/Errors") ;
  declareProperty( "TypeOfData", m_DataType="") ;

}


/*---------------------------------------------------------*/
JEMMon::~JEMMon()
/*---------------------------------------------------------*/
{
}


/*---------------------------------------------------------*/
StatusCode JEMMon::bookHistograms( bool isNewEventsBlock, 
				   bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in JEMMon::bookHistograms" << endreq;

  /** get a handle of StoreGate for access to the Event Store */
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  if (sc.isFailure()) 
    {
      mLog << MSG::ERROR
	   << "Unable to retrieve pointer to StoreGateSvc"
	   << endreq;
      return sc;
    }
 
  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }
	
  ManagedMonitorToolBase::LevelOfDetail_t LevelOfDetail=shift;
  if (m_DataType=="Sim") LevelOfDetail = expert;

  MonGroup JetElements_expert (this,(m_PathInRootFile+"_input").c_str() ,expert, run);
  HistoBooker expert_Booker(&JetElements_expert, &mLog, m_DataType);

  MonGroup JetElements_shift (this,(m_PathInRootFile+"_input").c_str() ,LevelOfDetail, run);
  HistoBooker shift_Booker(&JetElements_shift, &mLog, m_DataType);
  
  MonGroup JEM_DAQ ( this, (m_PathInRootFile+"_DAQ").c_str(), expert, run );
  HistoBooker DAQ_Booker(&JEM_DAQ, &mLog, m_DataType);

  MonGroup JEM_RoI ( this, (m_PathInRootFile+"_RoI").c_str(), LevelOfDetail, run );
  HistoBooker RoI_Booker(&JEM_RoI, &mLog, m_DataType);

  MonGroup JEM_Error( this, (m_ErrorPathInRootFile).c_str(), shift, run );
  HistoBooker Error_Booker(&JEM_Error, &mLog, "");

  if ( isNewEventsBlock|| isNewLumiBlock) { }

  if( isNewRun ) 
    {	
      Helper Help;
      m_NoEvents=0;

      //---------------------------- JetElements histos -----------------------------

      m_h_je_emeta = expert_Booker.book1F("emeta_JEM_input", "em TowerSum distribution per #eta  --  JEM input" , 50, -5, 5, "#eta" , "N");
      m_h_je_emeta ->SetBins(32,Help.JEEtaBinning());
      m_h_je_emeta ->SetStats(kFALSE);
      m_h_je_hadeta = expert_Booker.book1F("hadeta_JEM_input", "had TowerSum distribution per #eta  --  JEM input" , 50, -5, 5, "#eta" , "N");
      m_h_je_hadeta ->SetBins(32,Help.JEEtaBinning());
      m_h_je_hadeta ->SetStats(kFALSE);

      m_h_je_emphi = expert_Booker.book1F("emphi_JEM_input", "em TowerSum distribution per #phi  --  JEM input", 32, 0, 6.4, "#phi" , "N");
      m_h_je_emphi->SetBins(32,Help.JEPhiBinning());
      m_h_je_hadphi = expert_Booker.book1F("hadphi_JEM_input", "had TowerSum distribution per #phi  --  JEM input", 32, 0, 6.4, "#phi" , "N");
      m_h_je_hadphi->SetBins(32,Help.JEPhiBinning());
      m_h_je_emphi ->SetStats(kFALSE);
      m_h_je_hadphi ->SetStats(kFALSE);

      m_h_je_emenergy = expert_Booker.book1F("EmEnergy_JEM_input", "TowerSum EM energy distribution  --  JEM input", m_MaxEnergyRange, 0, m_MaxEnergyRange, "em energy [GeV]" , "N");
      m_h_je_hadenergy  = expert_Booker.book1F("HadEnergy_JEM_input", "TowerSum HAD energy distribution  --  JEM input", m_MaxEnergyRange, 0, m_MaxEnergyRange, "had energy [GeV]" , "N");
      m_h_je_emenergy->SetStats(kFALSE);
      m_h_je_hadenergy->SetStats(kFALSE);

      m_h_je_energy_emHitMap = shift_Booker.book2F("JE_EM_HitMap_energy_JEM_input", "#eta - #phi map of EM TowerSum weighted with energy  --  JEM input", 50, -5, 5, 32, 0, 6.4 , "#eta", "#phi");
      m_h_je_energy_emHitMap->SetBins(32,Help.JEEtaBinning(),32,Help.JEPhiBinning());
      m_h_je_energy_emHitMap->SetStats(kFALSE);
     

      m_h_je_energy_hadHitMap = shift_Booker.book2F("JE_HAD_HitMap_energy_JEM_input", "#eta - #phi map of HAD TowerSum weighted with energy  --  JEM input", 50, -5, 5, 32, 0, 6.4, "#eta", "#phi");	  
      m_h_je_energy_hadHitMap->SetBins(32,Help.JEEtaBinning(),32,Help.JEPhiBinning());
      m_h_je_energy_hadHitMap->SetStats(kFALSE);
      

      
      if (m_DataType=="BS")
	{
	  std::string name,title;
	  std::stringstream buffer;
	  
	  for (int i = 0; i < m_SliceNo; i++)
	    {
	      buffer.str("");
	      buffer<<i;
	      
	      name = "TowerSum_EM_HitMap_" + buffer.str() + "_JEM_input";
	      title = "#eta - #phi map of EM TowerSum for Timeslice " + buffer.str() +  "  --  JEM input";
	      m_h_je_emHitMap[i]=shift_Booker.book2F(name,title,50, -5, 5, 32, 0, 6.4, "#eta", "#phi");	  
	      m_h_je_emHitMap[i]->SetBins(32,Help.JEEtaBinning(),32,Help.JEPhiBinning());
	      m_h_je_emHitMap[i]->SetStats(kFALSE);
	      

	      buffer.str("");
	      buffer<<i;
	      
	      name = "TowerSum_HAD_HitMap_" + buffer.str() + "_JEM_input";
	      title = "#eta - #phi map of HAD TowerSum for Timeslice " + buffer.str() +  "  --  JEM input";
	      m_h_je_hadHitMap[i]=shift_Booker.book2F(name,title,50, -5, 5, 32, 0, 6.4, "#eta", "#phi");	  
	      m_h_je_hadHitMap[i]->SetBins(32,Help.JEEtaBinning(),32,Help.JEPhiBinning());
	      m_h_je_hadHitMap[i]->SetStats(kFALSE);
	      
	    }

	  // ----------------------------------- Error Histos ------------------------------------------------------
	  m_h_je_error = Error_Booker.book2F("JEM_Error","Error reports from JEM SubStatus Word",11,0.5,11.5,35,0.5,35.5,"","");
	 
	  m_h_je_error->GetXaxis()->SetBinLabel(1, "EM Parity");
	  m_h_je_error->GetXaxis()->SetBinLabel(2, "HAD Parity");
	  m_h_je_error->GetXaxis()->SetBinLabel(3, "Link down (em)");
	  m_h_je_error->GetXaxis()->SetBinLabel(4, "Link down (had)");
	  m_h_je_error->GetXaxis()->SetBinLabel(5, "GLinkParity");
	  m_h_je_error->GetXaxis()->SetBinLabel(6, "GLinkProtocol");
	  m_h_je_error->GetXaxis()->SetBinLabel(7, "BCNMismatch");
	  m_h_je_error->GetXaxis()->SetBinLabel(8, "FIFOOverflow");
	  m_h_je_error->GetXaxis()->SetBinLabel(9, "ModuleError");
	  m_h_je_error->GetXaxis()->SetBinLabel(10, "GLinkDown");
	  m_h_je_error->GetXaxis()->SetBinLabel(11, "GLinkTimeout");
	  m_h_je_error->SetStats(kFALSE);
         
	  
	  for (int i = 0; i < 16; i++)
	    {
	      buffer.str("");
	      buffer<<i;
	      name = "JEM " + buffer.str();
	      m_h_je_error->GetYaxis()->SetBinLabel((i+1), name.c_str());
	      
	      buffer.str("");
	      buffer<<i;
	      name = "JEM " + buffer.str();
	      m_h_je_error->GetYaxis()->SetBinLabel((i+19), name.c_str());
	    }
	  m_h_je_error->GetYaxis()->SetBinLabel(17, "Crate 0: ");
	  m_h_je_error->GetYaxis()->SetBinLabel(35, "Crate 1: ");
	}


      // number of triggered slice: not sure if we want this still
      // m_h_je_triggeredSlice=expert_Booker.book1F("JE_TriggeredSlice","Number of the Triggered Slice for JE",7,-0.5,6.5,"#Slice");
      
      //---------------------------- DAQ histos -----------------------------
      m_h_JEMHits_MainHits = DAQ_Booker.book1F("MainHits_JEM_DAQ", "Main Jet Hit Multiplicity per Threshold  --  JEM DAQ", 8, -0.5, 7.5,  "Threshold No.", "N");
      m_h_JEMHits_FwdHitsRight = DAQ_Booker.book1F("FwdHitsRight_JEM_DAQ", "Fwd Right Jet Hit Multiplicity per Threshold  --  JEM DAQ", 4, -0.5, 3.5, "Threshold No.", "N");
      m_h_JEMHits_FwdHitsLeft = DAQ_Booker.book1F("FwdHitsLeft_JEM_DAQ", "Fwd Left Jet Hit Multiplicity per Threshold  --  JEM DAQ", 4, -0.5, 3.5, "Threshold No.", "N");
      
      m_h_JEMEtSums_Ex = DAQ_Booker.book1F("Ex_JEM_DAQ", "JEM E_{x}^{JEM}  --  JEM DAQ", m_MaxEnergyRange, 0,m_MaxEnergyRange, "Ex [GeV]", "N");
      m_h_JEMEtSums_Ey = DAQ_Booker.book1F("Ey_JEM_DAQ", "JEM E_{y}^{JEM}  --  JEM DAQ", m_MaxEnergyRange, 0,m_MaxEnergyRange, "Ey [GeV]", "N");
      m_h_JEMEtSums_Et = DAQ_Booker.book1F("Et_JEM_DAQ", "JEM E_{t}^{JEM}  --  JEM DAQ", m_MaxEnergyRange, 0,m_MaxEnergyRange, "Et [GeV]", "N");
      
      //---------------------------- RoI histos -----------------------------
      m_h_JEMRoI_MainHits = RoI_Booker.book1F("MainHits_JEM_RoI", "Main Jet Hit Multiplicity per Threshold  --  JEM RoI", 8, -0.5, 7.5,  "Threshold No.", "N");      
      m_h_JEMRoI_FwdHitsRight = RoI_Booker.book1F("FwdHitsRight_JEM_RoI", "Forward Right Jet Hit Multiplicity per Threshold  --  JEM RoI", 4, -0.5, 3.5, "Threshold No.", "N");
      m_h_JEMRoI_FwdHitsLeft = RoI_Booker.book1F("FwdHitsLeft_JEM_RoI", "Forward Left Jet Hit Multiplicity per Threshold  --  JEM RoI", 4, -0.5, 3.5, "Threshold No.", "N");
      
      m_h_JEMDAQ_Hits_Map= DAQ_Booker.book2F("Hits_Map_JEM","HitMap of Hits per JEM",12,-0.5,11.5,35,0.5,35.5,"","");
      m_h_JEMDAQ_Hits_Map->SetStats(kFALSE);
      
	  for (int i = 0; i < 16; i++)
	    {

	      std::string name,title,xbin_name;
	      std::stringstream buffer;
	      //order of the thresholds: main thresholds bins 1-8 + fwd thresholds bins 9-12
              buffer.str("");
	      buffer<<i;
	      name = "JEM " + buffer.str();
              if (i<12) {
	      xbin_name = "Thrh" + buffer.str();
	      m_h_JEMDAQ_Hits_Map->GetXaxis()->SetBinLabel((i+1),xbin_name.c_str());
              
	      }
               
	      m_h_JEMDAQ_Hits_Map->GetYaxis()->SetBinLabel((i+1), name.c_str());
	         
	      m_h_JEMDAQ_Hits_Map->GetYaxis()->SetBinLabel((i+1+18), name.c_str());

	    }
	m_h_JEMDAQ_Hits_Map->GetYaxis()->SetBinLabel(17, "Crate 0: ");
	m_h_JEMDAQ_Hits_Map->GetYaxis()->SetBinLabel(35, "Crate 1: ");

      //---------------------------- HitThreshold per Eta-Phi -----------------------------
      for (int i=0;i<8;i++)
	{
	  std::string name,title;
	  std::stringstream buffer;

	  buffer.str("");
	  buffer<<i;
	  name = "MainThresh" + buffer.str()+"_EtaPhi_JEM_RoI";
	  title="#eta - #phi Map of Main Hits passing Threshold "+ buffer.str()+"  --  JEM RoI";

	  m_h_JEMRoI_MainThreshPerEtaPhi[i] = RoI_Booker.book2F(name.c_str(), title.c_str(), 50, -5, 5, 32,0,6.4, "#eta", "#phi");
	  m_h_JEMRoI_MainThreshPerEtaPhi[i]->SetBins(32,Help.JEEtaBinning(),32,Help.JEPhiBinning());
	  m_h_JEMRoI_MainThreshPerEtaPhi[i]->SetStats(kFALSE);
	
	}
      for (int i=0;i<4;i++)
	{
	  std::string name,title;
	  std::stringstream buffer;

	  buffer.str("");
	  buffer<<i;
	  name = "FwdThresh" + buffer.str()+"_EtaPhi_JEM_RoI";
	  title="#eta - #phi Map of Fwd Hits passing Threshold "+ buffer.str()+"  --  JEM RoI";

	  m_h_JEMRoI_FwdThreshPerEtaPhi[i] = RoI_Booker.book2F(name.c_str(), title.c_str(), 50, -5, 5, 32,0,6.4, "#eta", "#phi");
	  m_h_JEMRoI_FwdThreshPerEtaPhi[i]->SetBins(32,Help.JEEtaBinning(),32,Help.JEPhiBinning());
	  m_h_JEMRoI_FwdThreshPerEtaPhi[i]->SetStats(kFALSE);
	}

      //-------------------------------- Error Histos ------------------------------------------------------
      if (m_DataType=="BS")
	{
	  std::string name,title;
	  std::stringstream buffer;

	  m_h_JEMRoI_error = Error_Booker.book2F("JEMROI_Error","JEMRoI Parity and Saturation",5,0.5,5.5,35,0.5,35.5,"","");
	  m_h_JEMRoI_error->GetXaxis()->SetBinLabel(1, "Parity (Main Jets)");
	  m_h_JEMRoI_error->GetXaxis()->SetBinLabel(2, "Parity (Fwd Jets)");
	  
	  m_h_JEMRoI_error->GetXaxis()->SetBinLabel(4, "Saturation (Main Jets)");
	  m_h_JEMRoI_error->GetXaxis()->SetBinLabel(5, "Saturation (Fwd Jets)");
	  m_h_JEMRoI_error->SetStats(kFALSE);	   
	  	  
	  for (int i = 0; i < 16; i++)
	    {
	      buffer.str("");
	      buffer<<i;
	      name = "JEM " + buffer.str();
	      m_h_JEMRoI_error->GetYaxis()->SetBinLabel((i+1), name.c_str());
	      
	      buffer.str("");
	      buffer<<i;
	      name = "JEM " + buffer.str();
	      m_h_JEMRoI_error->GetYaxis()->SetBinLabel((i+1+18), name.c_str());
	    }
	  m_h_JEMRoI_error->GetYaxis()->SetBinLabel(17, "Crate 0: ");
	  m_h_JEMRoI_error->GetYaxis()->SetBinLabel(35, "Crate 1: ");
	
	m_h_JEM_ErrorSummary = Error_Booker.book1F("JEM Error Summary","Summary of Data Errors", 3,-0.5,2.5,"","Entries");
	m_h_JEM_ErrorSummary ->GetXaxis()->SetBinLabel(1, "Jet element errors");
	m_h_JEM_ErrorSummary ->GetXaxis()->SetBinLabel(2, "Status error");
	m_h_JEM_ErrorSummary ->GetXaxis()->SetBinLabel(3, "Jet Hit errors");
	
	}
       
    }
    
  return StatusCode( StatusCode::SUCCESS );
}



/*---------------------------------------------------------*/
StatusCode JEMMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  Helper Help;
  m_NoEvents++;

  // =============================================================================================
  // ================= Container: JetElements ====================================================
  // =============================================================================================
  // retrieve JetElements
  const JECollection* jetElements;
  StatusCode sc = m_storeGate->retrieve(jetElements, m_JetElementLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JetElements found in TES at " << m_JetElementLocation << endreq ;
      return StatusCode::SUCCESS;
    }
         
  // Step over all cells 
  JECollection::const_iterator it_je ;
  for( it_je = jetElements ->begin(); it_je < jetElements->end(); ++it_je )
    {	  
    
          LVL1::Coordinate coord((*it_je)->phi(),(*it_je)->eta());
	  LVL1::CoordToHardware ToHW;
	  int crate = ToHW.jepCrate(coord);
	  int module=ToHW.jepModule(coord);
	  int cord=ToHW.jepCoordinateWord(coord);
	  
      mLog << MSG::VERBOSE<<m_DataType <<" JE has coords ("<<(*it_je)->phi()<<", "<<(*it_je)->eta()
	   << " and energies : "<<(*it_je)->emEnergy()<<", "<<(*it_je)->hadEnergy()<<" (Em,Had)"<<" HW Crate:"<<crate<<"Module: "<<module<<" "<<cord<<endreq;

         
      if ((*it_je)->emEnergy()!=0) 
	{ 
	  m_h_je_emeta -> Fill( (*it_je)-> eta(), 1.);
	  m_h_je_emphi->Fill( (*it_je)->phi() , 1.);
	}
      if ((*it_je)->hadEnergy()!=0) 
	{ 
	  m_h_je_hadeta -> Fill( (*it_je)-> eta(), 1.);
	  m_h_je_hadphi->Fill( (*it_je)->phi() , 1.);
	}
      
      m_h_je_emenergy->Fill( (*it_je)->emEnergy() , 1.);
      m_h_je_hadenergy->Fill( (*it_je)->hadEnergy() , 1.);
      
      m_h_je_energy_emHitMap->Fill( (*it_je)->eta(),(*it_je)->phi() , (*it_je)->emEnergy());
      m_h_je_energy_hadHitMap->Fill( (*it_je)->eta(),(*it_je)->phi() ,(*it_je)->hadEnergy() ); 
      
      // number of triggered slice
      //m_h_je_triggeredSlice->Fill((*it_je)->peak(),1);
      
      // ------------------------------------------------------------------------------------------
      // ----------------- Histos filled only for BS data -----------------------------------------
      // ------------------------------------------------------------------------------------------
      if (m_DataType=="BS")
	{
	  // ----------------- HitMaps per time slice -----------------------------------------
	  for (int i = 0; i < m_SliceNo; i++)
	    {
	      if (i < static_cast<int> (((*it_je)->emEnergyVec()).size()))
		{
		  if ((*it_je)->emEnergyVec()[i] != 0) m_h_je_emHitMap[i]-> Fill( (*it_je)->eta(),(*it_je)->phi() ,1);
		  if ((*it_je)->hadEnergyVec()[i] != 0) m_h_je_hadHitMap[i]-> Fill( (*it_je)->eta(),(*it_je)->phi() ,1);
		} 
	    }

	  // ----------------- Error Histos -----------------------------------------
	  LVL1::DataError err((*it_je)->emError());
	  LVL1::DataError haderr((*it_je)->hadError());
	  LVL1::CoordToHardware ToHW;
	  LVL1::Coordinate coord((*it_je)->phi(),(*it_je)->eta());

	  int crate = ToHW.jepCrate(coord);
	  int module=ToHW.jepModule(coord);

	  // EM Parity
	  m_h_je_error->Fill(1,(crate*18 + module+1),err.get(1));
	  // HAD Parity
	  m_h_je_error->Fill(2,(crate*18 + module+1),haderr.get(1));
	  // PPM Link down: em.
	  m_h_je_error->Fill(3,(crate*18 + module+1),err.get(2));
	  // PPM Link down: had.
	  m_h_je_error->Fill(4,(crate*18 + module+1),haderr.get(2));
	  
	  mLog << MSG::VERBOSE<<"link down bit (em) "<< err.get(2)<<endreq;
	  mLog << MSG::VERBOSE<<"link down bit (had) "<< haderr.get(2)<<endreq;
	  mLog << MSG::VERBOSE<<"Parity em "<< err.get(1)<<endreq;
	  mLog << MSG::VERBOSE<<"Parity had "<< haderr.get(1)<<endreq;

	  // GLinkParity
	  m_h_je_error->Fill(5,(crate*18 + module+1),err.get(16));
          // GLinkProtocol
	  m_h_je_error->Fill(6,(crate*18 + module+1),err.get(17));
	  // BCNMismatch
	  m_h_je_error->Fill(7,(crate*18 + module+1),err.get(18));
	  // FIFOOverflow
	  m_h_je_error->Fill(8,(crate*18 + module+1),err.get(19));
	  // Module Error
	  m_h_je_error->Fill(9,(crate*18 + module+1),err.get(20));
	  // GLinkDown
	  m_h_je_error->Fill(10,(crate*18 + module+1),err.get(22));
	  // GLinkTimeout
	  m_h_je_error->Fill(11,(crate*18 + module+1),err.get(23));
	  
	  
	  
	  //Filling the Error Summary histogram
	 //Jet element errors 
	  m_h_JEM_ErrorSummary->Fill(1,err.get(1));
	  m_h_JEM_ErrorSummary->Fill(1,haderr.get(1));
	  m_h_JEM_ErrorSummary->Fill(1,err.get(2));
	  m_h_JEM_ErrorSummary->Fill(1,haderr.get(2));
	 //Errors from substatus word from ROD: JEM
	  if (err.get(16)!=0 or err.get(17)!=0 or err.get(18)!=0 or err.get(19)!=0 or err.get(20)!=0 or err.get(22)!=0 or err.get(23)!=0 ) {
	  m_h_JEM_ErrorSummary->Fill(2, 1);
	  }
	 
	  
	}   
    }
  

  // =============================================================================================
  // ================= Container: JEM Hits =======================================================
  // =============================================================================================
  // retrieve JEMHits collection from storegate
  const JEMHitsCollection* JEMHits;
  sc = m_storeGate->retrieve(JEMHits, m_JEMHitsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JEMHits found in TES at " << m_JEMHitsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  mLog<<MSG::DEBUG<<"-------------- "<< m_DataType<<" JEM Hits ---------------"<<endreq;
  
  // Step over all cells and process
  JEMHitsCollection::const_iterator it_JEMHits ;
  for( it_JEMHits  = JEMHits ->begin(); it_JEMHits < JEMHits -> end(); ++it_JEMHits )
    {	  
      // the binary hit information is represented by an integer number, 
      //that is converted to a string in order
      // to get the real binary information
      // later the multiplicities of the several thresholds are retrieved from this string
      std::string JEMHit = Help.Binary((*it_JEMHits)-> JetHits(),24);

      mLog<<MSG::DEBUG<<"Crate: "<< (*it_JEMHits)->crate()<<"  Module: "<<(*it_JEMHits)->module()
	  << "  JetHits: " <<JEMHit<<   endreq;
      
      //Main Jets
      if  ((*it_JEMHits)->forward()==0) 
	{
	  Help.FillHitsHisto(m_h_JEMHits_MainHits,JEMHit, 0, 8, 0, 3, &mLog);
	}
      
      //fwd jets a bit complicated!
      // fwd and main hits are contained in the same hitword
      if  ((*it_JEMHits)->forward()==1) 	{
	// JEMs No 0 and 8 are processing forward left hits,
	// JEMs No 7 and 15 forward right hits
	//left fwd hits
	if (((*it_JEMHits)-> module()==0) or((*it_JEMHits)-> module()==8) )
	  {
	    Help.FillHitsHisto(m_h_JEMHits_FwdHitsLeft, JEMHit, 0, 4, 8, 2, &mLog);
	  }
	//right fwd hits
	if (((*it_JEMHits)-> module()==7) or((*it_JEMHits)-> module()==15) )
	  {
	    Help.FillHitsHisto(m_h_JEMHits_FwdHitsRight, JEMHit, 0, 4, 8, 2, &mLog);
	  }
	// main hits
	Help.FillHitsHisto(m_h_JEMHits_MainHits , JEMHit, 0, 8, 0, 2, &mLog);
      }
   
      //m_h_JEMDAQ_Hit_Map:

      if ((*it_JEMHits)-> JetHits() !=0) {

	if ((*it_JEMHits)->forward()==0) {
	  if ((*it_JEMHits)->crate() == 0) {

	    for(int i=0; i<8;i++)
	    {
	      if ((Help.Multiplicity(JEMHit,3*i,1))!=0){
		m_h_JEMDAQ_Hits_Map->Fill(i,(*it_JEMHits)->module()+1,1);
	      }
	    }
	  }
	   if ((*it_JEMHits)->crate() == 1) {
	    
	     for(int i=0; i<8;i++) {
	       if ((Help.Multiplicity(JEMHit,3*i,1))!=0){
	       m_h_JEMDAQ_Hits_Map->Fill(i,(*it_JEMHits)->module()+19,1);
	       }
	     }
	   }
	}


        if ((*it_JEMHits)->forward()==1) {
	  if ((*it_JEMHits)->crate() == 0) {
	    for(int i=0; i<8;i++)
	      {
		if ((Help.Multiplicity(JEMHit,2*i,1))!=0){
		  m_h_JEMDAQ_Hits_Map->Fill(i,(*it_JEMHits)->module()+1,1);
		}
	      }
	    	  
	    for(int i=8; i<12;i++)
	      {
		if ((Help.Multiplicity(JEMHit,2*i,1))!=0){
		  m_h_JEMDAQ_Hits_Map->Fill(i,(*it_JEMHits)->module()+1,1);
		}
	      }
	  }

	  if ((*it_JEMHits)->crate() == 1) {
	    for(int i=0; i<8;i++)
	      {
		if ((Help.Multiplicity(JEMHit,2*i,1))!=0){
		  m_h_JEMDAQ_Hits_Map->Fill(i,(*it_JEMHits)->module()+19,1);
		}
	      }
	    	  
	    for(int i=8; i<12;i++)
	      {
		if ((Help.Multiplicity(JEMHit,2*i,1))!=0){
		  m_h_JEMDAQ_Hits_Map->Fill(i,(*it_JEMHits)->module()+19,1);
		}
	      }
	  }
	}
      }
    }
	     

 
  
  // =============================================================================================
  // ================= Container: JEM Et Sums ====================================================
  // =============================================================================================
  const JEMEtSumsCollection* JEMEtSums;
  sc = m_storeGate->retrieve(JEMEtSums, m_JEMEtSumsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JEMEtSums found in TES at "<< m_JEMEtSumsLocation << endreq ;
      return StatusCode::SUCCESS;
    }

  mLog<<MSG::DEBUG<<"-------------- "<< m_DataType<<" JEM Et Sums ---------------"<<endreq;

  // Step over all cells
  JEMEtSumsCollection::const_iterator it_JEMEtSums ;
  LVL1::QuadLinear expand;
  //  expand = new QuadLinear();
  int ex;
  int ey;
  int et;

  for( it_JEMEtSums  = JEMEtSums ->begin(); it_JEMEtSums < JEMEtSums -> end(); ++it_JEMEtSums )
    {	       
      // note: the energy values are compressed -> expand!
      ex=expand.Expand( (*it_JEMEtSums)-> Ex());
      ey=expand.Expand( (*it_JEMEtSums)-> Ey());
      et=expand.Expand( (*it_JEMEtSums)-> Et() );

      m_h_JEMEtSums_Ex -> Fill(ex, 1.); 
      m_h_JEMEtSums_Ey -> Fill(ey, 1.); 
      m_h_JEMEtSums_Et -> Fill(et, 1.); 
      mLog <<MSG::DEBUG<< " JEMEtSums Crate: "<<(*it_JEMEtSums)->crate()<<"  Module: "<<(*it_JEMEtSums)->module()
	   <<"   Ex: "<<  ex
	   <<"   Ey: "<<  ey 
	   <<"   Et: "<<  et  << "    Et compressed: "<< (*it_JEMEtSums)-> Et() <<endreq;
    }   


  // ==============================================================================================
  // ================= Container: JEM RoI =========================================================
  // ==============================================================================================
  const JemRoiCollection* JEMRoIs = 0;
  sc = m_storeGate->retrieve (JEMRoIs, m_JEMRoILocation);
  if (sc==StatusCode::FAILURE)
    {
      mLog <<MSG::INFO<<"No JEM RoIs found in TES at"<< m_JEMRoILocation<<endreq;
      return StatusCode::SUCCESS;    
    }

  mLog<<MSG::DEBUG<<"-------------- "<< m_DataType<<" JEM RoIs ---------------"<<endreq;

  // Step over all cells
  JemRoiCollection::const_iterator it_JEMRoIs ;

  for( it_JEMRoIs  = JEMRoIs ->begin(); it_JEMRoIs < JEMRoIs -> end(); ++it_JEMRoIs )
    {	  
      std::string JEMRoIHits ;
      //Main Jets
      if  ((*it_JEMRoIs)->forward()==0) 
	{
	  JEMRoIHits = Help.Binary((*it_JEMRoIs)-> hits(),8);
	  Help.FillHitsHisto(m_h_JEMRoI_MainHits , JEMRoIHits, 0, 8, 0, 1, &mLog);

	  // RoI HitMaps per threshold
	  LVL1::JEPRoIDecoder decoder;
	  LVL1::CoordinateRange coordRange = decoder.coordinate((*it_JEMRoIs)->roiWord());
	  double eta = coordRange.eta();
	  double phi = coordRange.phi();
	  
	  for (int i=0; i<8;i++)
	    {
	      if ((Help.Multiplicity(JEMRoIHits,i,1))!=0)
		{
		  m_h_JEMRoI_MainThreshPerEtaPhi[i]->Fill(eta,phi,1);
		}
	    }
	}      
      else 
	{
	  JEMRoIHits = Help.Binary((*it_JEMRoIs)-> hits(),4);
	  
	  // JEMs No 0 and 8 are processing forward left hits,
	  // JEMs No 7 and 15 forward right hits
	  //left fwd hits
	  if (((*it_JEMRoIs)-> jem()==0) or((*it_JEMRoIs)-> jem()==8) )
	    {
	      Help.FillHitsHisto(m_h_JEMRoI_FwdHitsLeft, JEMRoIHits, 0, 4, 0, 1, &mLog);
	    }
	  //right fwd hits
	  if (((*it_JEMRoIs)-> jem()==7) or((*it_JEMRoIs)-> jem()==15) )
	    {
	      Help.FillHitsHisto(m_h_JEMRoI_FwdHitsRight, JEMRoIHits, 0, 4, 0, 1, &mLog);
	    }

	  // RoI HitMaps per threshold
	  LVL1::JEPRoIDecoder decoder;
	  LVL1::CoordinateRange coordRange = decoder.coordinate((*it_JEMRoIs)->roiWord());
	  double eta = coordRange.eta();
	  double phi = coordRange.phi();
	  
	  for (int i=0; i<4;i++)
	    {
	      if ((Help.Multiplicity(JEMRoIHits,i,1))!=0)
		{
		  m_h_JEMRoI_FwdThreshPerEtaPhi[i]->Fill(eta,phi,1);
		}
	    }
	}
      
      mLog<<MSG::DEBUG<<"JEMRoI Word: "<<Help.Binary((*it_JEMRoIs)->roiWord(),32)<<endreq;
      mLog<<MSG::DEBUG<<"Crate: "<<(*it_JEMRoIs)->crate()<<"; JEM: "<<(*it_JEMRoIs)->jem()
	  <<"; forward: "<<(*it_JEMRoIs)->forward() <<"; Hits: "<<JEMRoIHits<<endreq;
      
      
      if (m_DataType=="BS")
	{
	  LVL1::DataError err((*it_JEMRoIs)->error());

	  int crate=(*it_JEMRoIs)->crate();
	  int module = (*it_JEMRoIs)->jem();

	  // Parity (Main Jets)
	  if ((*it_JEMRoIs)->forward()==0) m_h_JEMRoI_error->Fill(1,(crate*18 +module +1 ),err.get(1));
	  // Parity (Fwd Jets)
	  if ((*it_JEMRoIs)->forward()==1) m_h_JEMRoI_error->Fill(2,(crate*18 +module +1 ),err.get(1));
	  
	// Saturation (Main Jets)
	  if ((*it_JEMRoIs)->forward()==0) m_h_JEMRoI_error->Fill(4,(crate*18 +module +1 ),err.get(0));
	// Saturation (Fwd Jets)
	  if ((*it_JEMRoIs)->forward()==1) m_h_JEMRoI_error->Fill(5,(crate*18 +module +1 ),err.get(0));
	  
	 
	
	//Filling the Error Summary histogram
	 //Jet errors 
	  m_h_JEM_ErrorSummary->Fill(3,err.get(1));
	
	}     	
 

      
    }

  mLog<<MSG::DEBUG<<"--------------------------------------"<<endreq;

   return StatusCode( StatusCode::SUCCESS );
}

/*---------------------------------------------------------*/
StatusCode JEMMon::procHistograms( bool isEndOfEventsBlock, 
				   bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in procHistograms" << endreq ;

  if( isEndOfEventsBlock || isEndOfLumiBlock ) 
    {  
    }
	
  if(m_Offline==1)
    {
      if (m_DataType=="BS")
	{
	  if( isEndOfRun ) { 
	    std::stringstream buffer;
	    buffer.str("");
	    buffer<<m_NoEvents;
	    std::string title;
	    
	    title = m_h_JEMRoI_error-> GetTitle();
	    title=title + " | #events: " + buffer.str();
	    m_h_JEMRoI_error->SetTitle(title.c_str());
	  }
	}
    }
  return StatusCode( StatusCode::SUCCESS );
}
