// ********************************************************************
//
// NAME:     TrigT1CaloMonTool.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
//
// ********************************************************************

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/ISvcLocator.h"

#include "TString.h"

#include "StoreGate/StoreGateSvc.h"
#include "SGTools/StlVectorClids.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "TrigT1CaloMonitoring/PPrSpareMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"
#include "TrigT1CaloMonitoring/TrigT1CaloMonErrorTool.h"

#include "TrigT1CaloEvent/TriggerTower_ClassDEF.h"
#include "TrigT1CaloEvent/TriggerTowerCollection.h"
#include "TrigT1CaloUtils/DataError.h"


/*---------------------------------------------------------*/
PPrSpareMon::PPrSpareMon(const std::string & type, const std::string & name,
					 const IInterface* parent)
  : ManagedMonitorToolBase ( type, name, parent ),
    m_errorTool("TrigT1CaloMonErrorTool")
/*---------------------------------------------------------*/
{
  declareProperty("BS_TriggerTowerContainer",  m_TriggerTowerContainerName = "TriggerTowersSpare");
  declareProperty("ADCHitMap_Thresh",  m_TT_ADC_HitMap_Thresh = 40);

  declareProperty("PathInRootFile", m_PathInRootFile="L1Calo/PPM/SpareChannels") ;
  declareProperty("ErrorPathInRootFile", m_ErrorPathInRootFile="L1Calo/PPM/SpareChannels/Errors") ;
  declareProperty("OnlineTest", m_onlineTest = false,
                  "Test online code when running offline");

  // Maximum possible number of ADC slices
  m_SliceNo=15;
}

/*---------------------------------------------------------*/
PPrSpareMon::~PPrSpareMon()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode PPrSpareMon::initialize()
/*---------------------------------------------------------*/
{
  MsgStream log( msgSvc(), name() );

  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;

  sc = m_errorTool.retrieve();
  if( sc.isFailure() ) {
    log << MSG::ERROR << "Unable to locate Tool TrigT1CaloMonErrorTool"
                      << endreq;
    return sc;
  }
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode PPrSpareMon::bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  MsgStream log( msgSvc(), name() );
  log << MSG::DEBUG << "in PPrSpareMon::bookHistograms" << endreq;

  /** get a handle of StoreGate for access to the Event Store */
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  if (sc.isFailure()) 
    {
      log << MSG::ERROR
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


  MonGroup TT_ADC( this, m_PathInRootFile+"/ADC", shift, run );
  HistoBooker ADC_Booker(&TT_ADC, &log, "");

  MonGroup TT_Error( this, m_ErrorPathInRootFile, shift, run );
  HistoBooker Error_Booker(&TT_Error, &log, "");

  MonGroup TT_ErrorDetail( this, m_ErrorPathInRootFile+"/Detail", expert, run );
  HistoBooker ErrorDetail_Booker(&TT_ErrorDetail, &log, "");

  if ( isNewEventsBlock|| isNewLumiBlock) { }

  if( isNewRun )

    //if( isNewEventsBlock || isNewLumiBlock ) 
    {	
      Helper Help;

      std::string name,title;
      std::stringstream buffer;
	
	  

      //---------------------------- ADC Hitmaps for Triggered Timeslice -----------------------------
    

      buffer.str("");
      buffer<<m_TT_ADC_HitMap_Thresh;

      title="Spare Channels Hit Map of FADC > "+ buffer.str() + " for Triggered Timeslice";
      m_h_TT_HitMap_ADC=ADC_Booker.book2F("ppmspare_2d_tt_adc_HitMap",title,64,32.0,96.0,64,0.,64.,"crate/module","submodule/channel");
      setHitmapLabels(m_h_TT_HitMap_ADC);
      title="Spare Channels Profile Map of FADC for Triggered Timeslice";
      m_p_TT_HitMap_ADC=ADC_Booker.bookProfile2D("ppmspare_2d_tt_adc_ProfileMap",title,64,32.0,96.0,64,0.,64.,"crate/module","submodule/channel");
      setHitmapLabels(m_p_TT_HitMap_ADC);


      //-------------------------Summary of Errors-----------------------------------------------

      m_h_TT_Error=Error_Booker.book1F("ppmspare_1d_ErrorSummary","Spare Channels Summary of Errors",7,0.5,7.5,""); //without Ppm fw errors

      m_h_TT_Error->GetXaxis()->SetBinLabel(1, "GLinkParity");
      m_h_TT_Error->GetXaxis()->SetBinLabel(2, "GLinkProtocol");
      m_h_TT_Error->GetXaxis()->SetBinLabel(3, "FIFOOverflow");
      m_h_TT_Error->GetXaxis()->SetBinLabel(4, "ModuleError");
      m_h_TT_Error->GetXaxis()->SetBinLabel(5, "GLinkDown");
      m_h_TT_Error->GetXaxis()->SetBinLabel(6, "GLinkTimeout");
      m_h_TT_Error->GetXaxis()->SetBinLabel(7, "BCNMismatch");

      //---------------------------- SubStatus Word errors -----------------------------
      // divided in: crate, ROD status and PPm fw errors
      
       //L1Calo Substatus word
      m_h_TT_error_Crate_25=Error_Booker.book2F("ppmspare_2d_Status25","Spare Channels: Errors from TT SubStatus Word (crates 2-5)",7,0.5,7.5,71,0.5,71.5,"","");
            
      m_h_TT_error_Crate_25->GetXaxis()->SetBinLabel(1, "GLinkParity");
      m_h_TT_error_Crate_25->GetXaxis()->SetBinLabel(2, "GLinkProtocol");
      m_h_TT_error_Crate_25->GetXaxis()->SetBinLabel(3, "FIFOOverflow");
      m_h_TT_error_Crate_25->GetXaxis()->SetBinLabel(4, "ModuleError");
      m_h_TT_error_Crate_25->GetXaxis()->SetBinLabel(5, "GLinkDown");
      m_h_TT_error_Crate_25->GetXaxis()->SetBinLabel(6, "GLinkTimeout");
      m_h_TT_error_Crate_25->GetXaxis()->SetBinLabel(7, "BCNMismatch");
      m_h_TT_error_Crate_25->SetStats(kFALSE);
      

      //error bit field from ASIC data
      m_h_fwPpmError_Crate_25=Error_Booker.book2F("ppmspare_2d_ErrorField25","Spare Channels: Errors from ASIC error field (crates 2-5)",8,0.5,8.5,71,0.5,71.5,"","");

      m_h_fwPpmError_Crate_25->GetXaxis()->SetBinLabel(1, "ChannelDisabled");
      m_h_fwPpmError_Crate_25->GetXaxis()->SetBinLabel(2, "MCMAbsent");
      m_h_fwPpmError_Crate_25->GetXaxis()->SetBinLabel(3, "Timeout");
      m_h_fwPpmError_Crate_25->GetXaxis()->SetBinLabel(4, "ASICFull");
      m_h_fwPpmError_Crate_25->GetXaxis()->SetBinLabel(5, "EventMismatch");
      m_h_fwPpmError_Crate_25->GetXaxis()->SetBinLabel(6, "BunchMismatch");
      m_h_fwPpmError_Crate_25->GetXaxis()->SetBinLabel(7, "FIFOCorrupt");
      m_h_fwPpmError_Crate_25->GetXaxis()->SetBinLabel(8, "PinParity");
      m_h_fwPpmError_Crate_25->SetStats(kFALSE);

    
      for (int i=1; i<17; i+=2)
	{
	  buffer.str("");
	  buffer<<i-1;
	  
	  name = "PPM " + buffer.str();
	  m_h_TT_error_Crate_25->GetYaxis()->SetBinLabel(i,    name.c_str());
	  m_h_TT_error_Crate_25->GetYaxis()->SetBinLabel(i+18, name.c_str());
	  m_h_TT_error_Crate_25->GetYaxis()->SetBinLabel(i+36, name.c_str());
	  m_h_TT_error_Crate_25->GetYaxis()->SetBinLabel(i+54, name.c_str());

	  m_h_fwPpmError_Crate_25->GetYaxis()->SetBinLabel(i,    name.c_str());
	  m_h_fwPpmError_Crate_25->GetYaxis()->SetBinLabel(i+18, name.c_str());
	  m_h_fwPpmError_Crate_25->GetYaxis()->SetBinLabel(i+36, name.c_str());
	  m_h_fwPpmError_Crate_25->GetYaxis()->SetBinLabel(i+54, name.c_str());

	}

	  m_h_TT_error_Crate_25->GetYaxis()->SetBinLabel(17,    "Crate 2");
	  m_h_TT_error_Crate_25->GetYaxis()->SetBinLabel(17+18, "Crate 3");
	  m_h_TT_error_Crate_25->GetYaxis()->SetBinLabel(17+36, "Crate 4");
	  m_h_TT_error_Crate_25->GetYaxis()->SetBinLabel(17+54, "Crate 5");

	  m_h_fwPpmError_Crate_25->GetYaxis()->SetBinLabel(17,    "Crate 2");
	  m_h_fwPpmError_Crate_25->GetYaxis()->SetBinLabel(17+18, "Crate 3");
	  m_h_fwPpmError_Crate_25->GetYaxis()->SetBinLabel(17+36, "Crate 4");
	  m_h_fwPpmError_Crate_25->GetYaxis()->SetBinLabel(17+54, "Crate 5");

      m_h_ErrorDetails.clear();
      std::vector<std::string> errNames;
      errNames.push_back("Channel0Disabled");
      errNames.push_back("Channel1Disabled");
      errNames.push_back("Channel2Disabled");
      errNames.push_back("Channel3Disabled");
      errNames.push_back("MCMAbsent");
      errNames.push_back("");
      errNames.push_back("Timeout");
      errNames.push_back("ASICFull");
      errNames.push_back("EventMismatch");
      errNames.push_back("BunchMismatch");
      errNames.push_back("FIFOCorrupt");
      errNames.push_back("PinParity");
      for (int error = 0; error < 12; error+=2) 
        {
          for (int crate = 2; crate < 6; crate+=2)
            {
	      buffer.str("");
	      buffer<<crate;
	      std::string name = "ppmspare_2d_"+errNames[error]+errNames[error+1]+"Crate"+buffer.str();
	      std::string title = "ASIC Errors "+errNames[error]+" "+errNames[error+1]+" for Crates "+buffer.str();
	      buffer.str("");
	      buffer<<(crate+1);
	      name += buffer.str();
	      title += "-"+buffer.str();
	      TH2F* hist = 0;
	      if (error != 4) hist = ErrorDetail_Booker.book2F(name,title,32,0,32,32,0,32,"MCM","Crate/Module");
	      else            hist = ErrorDetail_Booker.book2F(name,title,16,0,16,32,0,32,"MCM","Crate/Module");
	      m_h_ErrorDetails.push_back(hist);
	      hist->SetStats(kFALSE);
	      for (int mcm = 0; mcm < 16; mcm+=2)
	        {
		  if (mcm == 0)
		    {
		      hist->GetXaxis()->SetBinLabel(1, errNames[error].c_str());
		      if (error != 4) hist->GetXaxis()->SetBinLabel(17, errNames[error+1].c_str());
                    }
                  else
		    {
		      buffer.str("");
		      buffer<<mcm;
		      hist->GetXaxis()->SetBinLabel(1+mcm, buffer.str().c_str());
		      if (error != 4) hist->GetXaxis()->SetBinLabel(17+mcm, buffer.str().c_str());
                    }
                }
              for (int cr = 0; cr < 2; ++cr)
	        {
		  for (int module = 0; module < 16; module+=2)
		    {
		      buffer.str("");
		      buffer<<(cr+crate)<<"/"<<module;
		      hist->GetYaxis()->SetBinLabel(1+cr*16+module, buffer.str().c_str());
                    }
                }
            }
        }

	  
      //---------------------------- number of triggered slice -----------------------------
      m_h_TT_triggeredSlice=ADC_Booker.book1F("ppmspare_1d_tt_adc_TriggeredSlice","Spare Channels Number of the Triggered Slice",m_SliceNo,-0.5,m_SliceNo-0.5,"#Slice");
      
	}	

  if ( isNewLumiBlock )
    {
    }
    
  return StatusCode::SUCCESS;
	}

/*---------------------------------------------------------*/
StatusCode PPrSpareMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());

  log << MSG::DEBUG << "in fillHistograms()" << endreq;

  // Skip events believed to be corrupt

  if (m_errorTool->corrupt()) {
    log << MSG::DEBUG << "Skipping corrupt event" << endreq;
    return StatusCode::SUCCESS;
  }

  // Error vector for global overview
  std::vector<int> overview(8);

  //Retrieve TriggerTowers from SG
  StatusCode sc;
  const TriggerTowerCollection* TriggerTowerTES = 0; 
  if (m_storeGate->contains<TriggerTowerCollection>(m_TriggerTowerContainerName)) {
    sc = m_storeGate->retrieve(TriggerTowerTES, m_TriggerTowerContainerName); 
  } else sc = StatusCode::FAILURE;
  if( sc.isFailure() ) 
    {
      log << MSG::DEBUG << "No TriggerTower found in TES at "<< m_TriggerTowerContainerName<< endreq ;
      return StatusCode::SUCCESS;
    }

    
  // =============================================================================================
  // ================= Container: TriggerTower ===================================================
  // =============================================================================================

  TriggerTowerCollection::const_iterator TriggerTowerIterator    = TriggerTowerTES->begin(); 
  TriggerTowerCollection::const_iterator TriggerTowerIteratorEnd = TriggerTowerTES->end(); 
 
  for (; TriggerTowerIterator != TriggerTowerIteratorEnd; ++TriggerTowerIterator) 
    {
    
	  

     //---------------------------- ADC HitMaps -----------------------------


      double crateModule      = (*TriggerTowerIterator)->eta();
      double submoduleChannel = (*TriggerTowerIterator)->phi();
      const int adc = (*TriggerTowerIterator)->emADC()[(*TriggerTowerIterator)->emADCPeak()];
      if (adc > m_TT_ADC_HitMap_Thresh) m_h_TT_HitMap_ADC->Fill(crateModule, submoduleChannel, 1);
      m_p_TT_HitMap_ADC->Fill(crateModule, submoduleChannel, adc);
    

      //---------------------------- SubStatus Word errors -----------------------------

      using LVL1::DataError;
      DataError error((*TriggerTowerIterator)-> emError());

      int icm       = crateModule;
      int isc       = submoduleChannel;
      int crate     = icm/16;
      int module    = icm%16;
      int submodule = isc/4;
      int channel   = isc%4;
   
      //Summary

      if (error.get(DataError::GLinkParity))   m_h_TT_Error->Fill(1);
      if (error.get(DataError::GLinkProtocol)) m_h_TT_Error->Fill(2);
      if (error.get(DataError::FIFOOverflow))  m_h_TT_Error->Fill(3);
      if (error.get(DataError::ModuleError))   m_h_TT_Error->Fill(4);
      if (error.get(DataError::GLinkDown))     m_h_TT_Error->Fill(5);
      if (error.get(DataError::GLinkTimeout))  m_h_TT_Error->Fill(6);
      if (error.get(DataError::BCNMismatch))   m_h_TT_Error->Fill(7);
      

	//---------------- per crate and module --------------------  m_h_TT_error_Crate_25
      if (crate > 1 && crate < 6)
        {
          int ypos = (module+1)+((crate-2)*18);
          if (error.get(DataError::ChannelDisabled)) m_h_fwPpmError_Crate_25->Fill(1,ypos);
          if (error.get(DataError::MCMAbsent))       m_h_fwPpmError_Crate_25->Fill(2,ypos);
          if (error.get(DataError::Timeout))         m_h_fwPpmError_Crate_25->Fill(3,ypos);
          if (error.get(DataError::ASICFull))        m_h_fwPpmError_Crate_25->Fill(4,ypos);
          if (error.get(DataError::EventMismatch))   m_h_fwPpmError_Crate_25->Fill(5,ypos);
          if (error.get(DataError::BunchMismatch))   m_h_fwPpmError_Crate_25->Fill(6,ypos);
          if (error.get(DataError::FIFOCorrupt))     m_h_fwPpmError_Crate_25->Fill(7,ypos);
          if (error.get(DataError::PinParity))       m_h_fwPpmError_Crate_25->Fill(8,ypos);
      		  
          if (error.get(DataError::GLinkParity))     m_h_TT_error_Crate_25->Fill(1,ypos);
          if (error.get(DataError::GLinkProtocol))   m_h_TT_error_Crate_25->Fill(2,ypos);
          if (error.get(DataError::FIFOOverflow))    m_h_TT_error_Crate_25->Fill(3,ypos);
          if (error.get(DataError::ModuleError))     m_h_TT_error_Crate_25->Fill(4,ypos);
          if (error.get(DataError::GLinkDown))       m_h_TT_error_Crate_25->Fill(5,ypos);
          if (error.get(DataError::GLinkTimeout))    m_h_TT_error_Crate_25->Fill(6,ypos);
          if (error.get(DataError::BCNMismatch))     m_h_TT_error_Crate_25->Fill(7,ypos);

          // Detailed plots by MCM
          ypos = (crate%2)*16+module;
          if (error.get(DataError::ChannelDisabled)) m_h_ErrorDetails[(channel/2)*4+(crate-2)/2]->Fill((channel%2)*16+submodule, ypos);
          if (error.get(DataError::MCMAbsent))       m_h_ErrorDetails[4+(crate-2)/2]->Fill(submodule, ypos);
          if (error.get(DataError::Timeout))         m_h_ErrorDetails[6+(crate-2)/2]->Fill(submodule, ypos);
          if (error.get(DataError::ASICFull))        m_h_ErrorDetails[6+(crate-2)/2]->Fill(16+submodule, ypos);
          if (error.get(DataError::EventMismatch))   m_h_ErrorDetails[8+(crate-2)/2]->Fill(submodule, ypos);
          if (error.get(DataError::BunchMismatch))   m_h_ErrorDetails[8+(crate-2)/2]->Fill(16+submodule, ypos);
          if (error.get(DataError::FIFOCorrupt))     m_h_ErrorDetails[10+(crate-2)/2]->Fill(submodule, ypos);
          if (error.get(DataError::PinParity))       m_h_ErrorDetails[10+(crate-2)/2]->Fill(16+submodule, ypos);
        }

      if (error.get(DataError::ChannelDisabled) || error.get(DataError::MCMAbsent)) overview[crate] |= 1;

      if (error.get(DataError::Timeout)       || error.get(DataError::ASICFull)      ||
          error.get(DataError::EventMismatch) || error.get(DataError::BunchMismatch) ||
          error.get(DataError::FIFOCorrupt)   || error.get(DataError::PinParity)) overview[crate] |= (1 << 1);

      if (error.get(DataError::GLinkParity)  || error.get(DataError::GLinkProtocol) ||
          error.get(DataError::FIFOOverflow) || error.get(DataError::ModuleError)   ||
          error.get(DataError::GLinkDown)    || error.get(DataError::GLinkTimeout)  ||
          error.get(DataError::BCNMismatch)) overview[crate] |= (1 << 2);
     
      // number of triggered slice
      m_h_TT_triggeredSlice->Fill((*TriggerTowerIterator)->emADCPeak(),1);

	}	     
     
  // Write overview vector to StoreGate
  std::vector<int>* save = new std::vector<int>(overview);
  sc = m_storeGate->record(save, "L1CaloPPMSpareErrorVector");
  if (sc != StatusCode::SUCCESS)
    {
      log << MSG::ERROR << "Error recording PPMSpare error vector in TES "
          << endreq;
      return sc;
    }

  
  return StatusCode( StatusCode::SUCCESS );
}


   
/*---------------------------------------------------------*/
StatusCode PPrSpareMon::procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in procHistograms" << endreq ;

  if( isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun ) 
    {  
    }
	
  return StatusCode( StatusCode::SUCCESS );
}


/*---------------------------------------------------------*/
void PPrSpareMon::setHitmapLabels(TH2* hist) {
/*---------------------------------------------------------*/

  for (int crate = 2; crate < 6; ++crate) {
    for (int module = 0; module < 16; module+=4) {
      std::ostringstream cnum;
      cnum << crate << "/" << module;
      hist->GetXaxis()->SetBinLabel((crate-2)*16 + module + 1, cnum.str().c_str());
    }
  }
  for (int submodule = 0; submodule < 16; ++submodule) {
    std::ostringstream cnum;
    cnum << submodule << "/0";
    hist->GetYaxis()->SetBinLabel(submodule*4 + 1, cnum.str().c_str());
  }
  hist->SetStats(kFALSE);
}
