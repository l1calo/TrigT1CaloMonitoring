// ********************************************************************
//
// NAME:     PPrPerfMon.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
// 
// Checks:
//    
// ================= "PPr: ADC -> BCID, LUT" ==============================================
// Check	Functionality	Status	
//            TriggerTowerTool(PPrDAQ.ADC) = PPrDAQ.(BCDI,LUT) (not for all Timeslices exact!)     
//
// ********************************************************************

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ITHistSvc.h"

#include <TH1D.h>
#include <TH2D.h>

#include "TString.h"

#include "StoreGate/StoreGateSvc.h"

#include "TrigT1CaloMonitoring/PPrPerfMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "TrigT1Calo/TriggerTowerCollection.h"
#include "TrigT1Calo/TrigT1CaloDict.h"
#include "TrigT1Calo/TriggerTower_ClassDEF.h"

#include "GaudiKernel/IService.h"
#include "GaudiKernel/IToolSvc.h"


/*---------------------------------------------------------*/
PPrPerfMon::PPrPerfMon(const std::string & type, const std::string & name,
					 const IInterface* parent)
  : ManagedMonitorToolBase ( type, name, parent )
/*---------------------------------------------------------*/
{
 declareProperty( "PathInRootFile", m_PathInRootFile="Stats/CMM") ;
  declareProperty( "DataType", m_DataType="") ;
}

/*---------------------------------------------------------*/
PPrPerfMon::~PPrPerfMon()
/*---------------------------------------------------------*/
{
}


/*---------------------------------------------------------*/
StatusCode PPrPerfMon::bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  MsgStream log( msgSvc(), name() );
  log << MSG::DEBUG << "in PPrPerfMon::bookHistograms" << endreq;

  /** get a handle of StoreGate for access to the Event Store */
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  if (sc.isFailure()) 
    {
      log << MSG::ERROR << "Unable to retrieve pointer to StoreGateSvc" << endreq;
      return sc;
    }
  
  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  if ( isNewEventsBlock|| isNewLumiBlock) { }

  if(isNewRun  ) 
    {
    }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode PPrPerfMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  
  log << MSG::DEBUG << "in fillHistograms()" << endreq;


  // =============================================================================================
  // ================= "PPr: ADC -> BCID, LUT" ==============================================
  // =============================================================================================

  // ---------------------------------------------------------------------------------------------
  // Check	Functionality	Status	
  //            TriggerTowerTool(PPrDAQ.ADC) = PPrDAQ.(BCDI,LUT) (not for all Timeslices exact!)
  // ---------------------------------------------------------------------------------------------





  return StatusCode( StatusCode::SUCCESS );

}
/*---------------------------------------------------------*/
StatusCode PPrPerfMon::procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
  if( isEndOfEventsBlock || isEndOfLumiBlock ) 
    {  
    }
	
  if( isEndOfRun ) 
    {   
    }
  
  return StatusCode( StatusCode::SUCCESS );
}

