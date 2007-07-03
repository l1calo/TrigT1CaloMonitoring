// ********************************************************************
//
// NAME:     TrigT1CaloMonTool.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Ethan-Etienne Woehrling (eew@hep.ph.bham.ac.uk)
//           Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
//
// ********************************************************************

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ITHistSvc.h"

#include <TH1D.h>
#include <TH2D.h>

#include "TString.h"

#include "StoreGate/StoreGateSvc.h"

#include "CaloEvent/CaloTower.h"
#include "CaloEvent/CaloTowerContainer.h"

#include "CaloEvent/CaloCell.h"
#include "CaloEvent/CaloCellContainer.h"

#include "TrigT1CaloMonitoring/TriggerTowerMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "TrigT1Calo/TriggerTowerCollection.h"
#include "TrigT1Calo/TrigT1CaloDict.h"
#include "TrigT1Calo/TriggerTower_ClassDEF.h"

//#include "TrigT1CaloCalibTools/L1CaloTTIdTools.h"


/*---------------------------------------------------------*/
TriggerTowerMon::TriggerTowerMon(const std::string & type, const std::string & name,
					 const IInterface* parent)
  : ManagedMonitorToolBase ( type, name, parent )
/*---------------------------------------------------------*/
{
  declareProperty("BS_TriggerTowerContainer",  m_TriggerTowerContainerName = "LVL1TriggerTowers");
  
  declareProperty( "PathInRootFile", m_PathInRootFile="Stats/CMM") ;
  declareProperty( "DataType", m_DataType="") ;
}

/*---------------------------------------------------------*/
TriggerTowerMon::~TriggerTowerMon()
/*---------------------------------------------------------*/
{
}


/*---------------------------------------------------------*/
StatusCode TriggerTowerMon::bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  MsgStream log( msgSvc(), name() );
  log << MSG::DEBUG << "in TriggerTowerMon::bookHistograms" << endreq;

  /** get a handle of StoreGate for access to the Event Store */
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  if (sc.isFailure()) 
    {
      log << MSG::ERROR
	   << "Unable to retrieve pointer to StoreGateSvc"
	   << endreq;
      return sc;
    }

  // Get a pointer to DetectorStore services
  sc = service("DetectorStore", m_detStore);
  if (sc.isFailure()) 
    {
      log << MSG::ERROR << "Cannot access DetectorStore" << endreq;
      return StatusCode::FAILURE;;
    }

  // Retrieve the CaloIdManager from the detector store
  sc = m_detStore->retrieve(m_caloMgr);
  if (sc.isFailure()) 
    {
      log << MSG::ERROR << "Unable to retrieve CaloIdManager from DetectorStore" << endreq;
      return StatusCode::FAILURE;
    }
  
  // Use the CaloIdManager to get a pointer to an instance of the CaloLVL1_ID helper
  m_lvl1Helper = m_caloMgr->getLVL1_ID();
  if (!m_lvl1Helper) 
    {
      log << MSG::ERROR << "Could not access CaloLVL1_ID helper" << endreq;
      return StatusCode::FAILURE;
    }


  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  MonGroup TT_EmADCPeak( this, (m_PathInRootFile+"/EmADCPeak").c_str(), expert, eventsBlock );
  HistoBooker* EmADCPeak_Booker = new HistoBooker(&TT_EmADCPeak, &log, m_DataType);

  MonGroup TT_HadADCPeak( this, (m_PathInRootFile+"/HadADCPeak").c_str(), expert, eventsBlock );
  HistoBooker* HadADCPeak_Booker = new HistoBooker(&TT_HadADCPeak, &log, m_DataType);

  MonGroup TT_HitMaps( this, (m_PathInRootFile+"/HitMaps").c_str(), expert, eventsBlock );
  HistoBooker* HitMaps_Booker = new HistoBooker(&TT_HitMaps, &log, m_DataType);

  MonGroup TT_EmLUTPeak( this, (m_PathInRootFile+"/EmLUTPeak").c_str(), expert, eventsBlock );
  HistoBooker* EmLUTPeak_Booker = new HistoBooker(&TT_EmLUTPeak, &log, m_DataType);

  MonGroup TT_HadLUTPeak( this, (m_PathInRootFile+"/HadLUTPeak").c_str(), expert, eventsBlock );
  HistoBooker* HadLUTPeak_Booker = new HistoBooker(&TT_HadLUTPeak, &log, m_DataType);

  if( isNewEventsBlock || isNewLumiBlock ) 
    {	

      Helper* Help = new Helper();


   //Energy distribution per Channel
      std::vector<Identifier>::const_iterator tower_it = m_lvl1Helper->tower_begin();
      std::string name,title;
      std::stringstream buffer, etabuffer,phibuffer;
      int emTT=0;
      int hadTT=0;

      for(;tower_it!=m_lvl1Helper->tower_end();++tower_it) 
	{
	  Identifier towerId = (*tower_it);
	  buffer.str("");
	  etabuffer.str("");
	  phibuffer.str("");
	  buffer<<towerId;
	  etabuffer<<IDeta(towerId);
	  phibuffer<<IDphi(towerId);

	  if (m_lvl1Helper->sampling(towerId)==0) //EM TT
	    {
	      name = "emADCTT_" + buffer.str();
	      title = "emADC TT: eta | phi = " + etabuffer.str() + " | " + phibuffer.str();
	      m_h_TT_EmADCPeak[towerId]=EmADCPeak_Booker->book1F(name,title,50,-0.5,49.5,"Et [GeV]");

	      name = "emLUTTT_" + buffer.str();
	      title = "emLUT TT: eta | phi = " + etabuffer.str() + " | " + phibuffer.str();
	      m_h_TT_EmLUTPeak[towerId]=EmLUTPeak_Booker->book1F(name,title,255,-0.5,254.5,"Et [GeV]");
	      emTT+=1;
	    }

	  if (m_lvl1Helper->sampling(towerId)==1) //Had TT
	    {
	      name = "hadADCTT_" + buffer.str();
	      title = "hadADC TT: eta | phi = " + etabuffer.str() + " | " + phibuffer.str();
	      m_h_TT_HadADCPeak[towerId]=HadADCPeak_Booker->book1F(name,title,50,-0.5,49.5,"Et [GeV]");

	      name = "hadLUTTT_" + buffer.str();
	      title = "hadLUT TT: eta | phi = " + etabuffer.str() + " | " + phibuffer.str();
	      m_h_TT_HadLUTPeak[towerId]=HadLUTPeak_Booker->book1F(name,title,255,-0.5,254.5,"Et [GeV]");
	      hadTT+=1;
	    }

	  //TT Conversion Check:
	  double etaCheck = IDeta(towerId);
	  double phiCheck = IDphi(towerId);

	  int detside= pos_neg_z(etaCheck);
	  int detregion = region(etaCheck);
	  int eta=etaIndex(etaCheck);
	  int phi=phiIndex(etaCheck, phiCheck);

	  Identifier towerIdCheck;

	  towerIdCheck = m_lvl1Helper->tower_id(detside, m_lvl1Helper->sampling(towerId), detregion,eta,phi );  

	  if (towerIdCheck=!towerId)


	    {
	      log << MSG::ERROR << "TT conversion failed for TTId = "<< towerId<<endreq;
	    }
	}

      log << MSG::DEBUG << "em TT " <<emTT<< endreq;
      log << MSG::DEBUG << "had TT " <<hadTT<< endreq;

   m_h_TT_EmHitMap_30GeV=HitMaps_Booker->book2F("TTEm30GeV","TT Em E>30GeV",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_TT_EmHitMap_30GeV->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning()); 
  
   m_h_TT_EmHitMap_10GeV=HitMaps_Booker->book2F("TTEm10GeV","TT Em E>10GeV",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_TT_EmHitMap_10GeV->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning()); 
  
   m_h_TT_EmHitMap_5GeV=HitMaps_Booker->book2F("TTEm5GeV","TT Em E>5GeV",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_TT_EmHitMap_5GeV->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning()); 
  
   m_h_TT_EmHitMap_2GeV=HitMaps_Booker->book2F("TTEm2GeV","TT Em E>2GeV",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_TT_EmHitMap_2GeV->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());  
 
   m_h_TT_EmHitMap_1GeV=HitMaps_Booker->book2F("TTEm1GeV","TT Em E>1GeV",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_TT_EmHitMap_1GeV->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning()); 
  

   m_h_TT_HadHitMap_30GeV=HitMaps_Booker->book2F("TTHad30GeV","TT Had E>30GeV",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_TT_HadHitMap_30GeV->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());   

   m_h_TT_HadHitMap_10GeV=HitMaps_Booker->book2F("TTHad10GeV","TT Had E>10GeV",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_TT_HadHitMap_10GeV->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());   

   m_h_TT_HadHitMap_5GeV=HitMaps_Booker->book2F("TTHad5GeV","TT Had E>5GeV",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_TT_HadHitMap_5GeV->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());   

   m_h_TT_HadHitMap_2GeV=HitMaps_Booker->book2F("TTHad2GeV","TT Had E>2GeV",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_TT_HadHitMap_2GeV->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());   

   m_h_TT_HadHitMap_1GeV=HitMaps_Booker->book2F("TTHad1GeV","TT Had E>1GeV",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_TT_HadHitMap_1GeV->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());   
    }

  if( isNewRun ) { }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TriggerTowerMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  
  log << MSG::DEBUG << "in fillHistograms()" << endreq;

  //Retrieve TriggerTowers from SG
  const TriggerTowerCollection* TriggerTowerTES = 0; 
  StatusCode sc = m_storeGate->retrieve(TriggerTowerTES, m_TriggerTowerContainerName); 
  if( (sc==StatusCode::FAILURE) ) 
    {
      log << MSG::INFO
	   << "No TriggerTower found in TES at "
	   << m_TriggerTowerContainerName
	   << endreq ;
      return StatusCode::SUCCESS;
    }


  // =========================================================================================
  // ================================== TriggerTower Plots ===================================
  // =========================================================================================

  // TTs are only filled/stored if they contain energy, 
  // so looping over all TTs will only loop over those with energy

  TriggerTowerCollection::const_iterator TriggerTowerIterator    = TriggerTowerTES->begin(); 
  TriggerTowerCollection::const_iterator TriggerTowerIteratorEnd = TriggerTowerTES->end(); 

  for (; TriggerTowerIterator != TriggerTowerIteratorEnd; ++TriggerTowerIterator) 
    {
      Identifier EmTowerId,HadTowerId;

      int detside= pos_neg_z((*TriggerTowerIterator)->eta());
      int detregion = region((*TriggerTowerIterator)->eta());
      int eta=etaIndex((*TriggerTowerIterator)->eta());
      int phi=phiIndex((*TriggerTowerIterator)->eta(), (*TriggerTowerIterator)->phi());

      EmTowerId = m_lvl1Helper->tower_id(detside, 0, detregion,eta,phi );  //em
      int EmEnergy = (*TriggerTowerIterator)->emADC()[(*TriggerTowerIterator)->emPeak()];
      if (EmEnergy > 49)
	{
	  m_h_TT_EmADCPeak[EmTowerId]->Fill(49,1);
	} 
      else
	{
	  m_h_TT_EmADCPeak[EmTowerId]->Fill(EmEnergy,1);
	}  
   
      EmEnergy = (*TriggerTowerIterator)->emLUT()[(*TriggerTowerIterator)->emPeak()];
      if (EmEnergy > 49)
	{
	  m_h_TT_EmLUTPeak[EmTowerId]->Fill(49,1);
	} 
      else
	{
	  m_h_TT_EmLUTPeak[EmTowerId]->Fill(EmEnergy,1);
	}     
      if (EmEnergy>30)
	{
	  m_h_TT_EmHitMap_30GeV->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (EmEnergy>10)
	{
	  m_h_TT_EmHitMap_10GeV->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (EmEnergy>5)
	{
	  m_h_TT_EmHitMap_5GeV->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (EmEnergy>2)
	{
	  m_h_TT_EmHitMap_2GeV->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (EmEnergy>1)
	{
	  m_h_TT_EmHitMap_1GeV->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}

      
      HadTowerId = m_lvl1Helper->tower_id(detside, 1, detregion,eta,phi );  //had
      int HadEnergy = (*TriggerTowerIterator)->hadADC()[(*TriggerTowerIterator)->hadPeak()];
      if (HadEnergy > 49)
	{
	  m_h_TT_HadADCPeak[HadTowerId]->Fill(49,1);
	} 
      else
	{
	  m_h_TT_HadADCPeak[HadTowerId]->Fill(HadEnergy,1);
	}  
   
      HadEnergy = (*TriggerTowerIterator)->hadLUT()[(*TriggerTowerIterator)->hadPeak()];
      if (HadEnergy > 49)
	{
	  m_h_TT_HadLUTPeak[HadTowerId]->Fill(49,1);
	} 
      else
	{
	  m_h_TT_HadLUTPeak[HadTowerId]->Fill(HadEnergy,1);
	}     
      if (HadEnergy>30)
	{
	  m_h_TT_HadHitMap_30GeV->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (HadEnergy>10)
	{
	  m_h_TT_HadHitMap_10GeV->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (HadEnergy>5)
	{
	  m_h_TT_HadHitMap_5GeV->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (HadEnergy>2)
	{
	  m_h_TT_HadHitMap_2GeV->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (HadEnergy>1)
	{
	  m_h_TT_HadHitMap_1GeV->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}

    //TT Conversion Check:
    double EmEtaCheck = IDeta(EmTowerId);
    double EmPhiCheck = IDphi(EmTowerId);
    double HadEtaCheck = IDeta(HadTowerId);
    double HadPhiCheck = IDphi(HadTowerId);

    if ((EmEtaCheck=!(*TriggerTowerIterator)->eta()) or (HadEtaCheck=!(*TriggerTowerIterator)->eta()))
      {
	log << MSG::ERROR << "TT Id conversion failed for eta = "<<(*TriggerTowerIterator)->eta() <<endreq;
      }
    if ((EmPhiCheck=!(*TriggerTowerIterator)->phi()) or (HadPhiCheck=!(*TriggerTowerIterator)->phi()))
      {
	log << MSG::ERROR << "TT Id conversion failed for phi = "<<(*TriggerTowerIterator)->phi() <<endreq;
      }

    }

  
  return StatusCode( StatusCode::SUCCESS );

}
/*---------------------------------------------------------*/
StatusCode TriggerTowerMon::procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun )
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

#define BASEDETA 0.1
#define BASEDPHI 0.098175
#define FCALDETA 0.425

#define ETAMAXREGION0 2.5
#define ETAMAXREGION1 3.1
#define ETAMAXREGION2 3.2
#define ETAMAXREGION3 4.9

#define ETAMIN -4.9
#define PHIMIN 0.

/*---------------------------------------------------------*/
int TriggerTowerMon::pos_neg_z(double eta) const 
/*---------------------------------------------------------*/
{
  int detside = 1;
  if (eta<0) detside=-1;

  return detside;

}

/*---------------------------------------------------------*/
int TriggerTowerMon::region(double eta) const 
/*---------------------------------------------------------*/
{
        if(fabs(eta)<ETAMAXREGION0) {
                return 0;
                
        } else if(fabs(eta)<ETAMAXREGION1) {
                return 1;

        } else if (fabs(eta)<ETAMAXREGION2) {
                return 2;

        } else {
                return 3;                
        }
}

/*---------------------------------------------------------*/
double TriggerTowerMon::etaWidth(double eta) const 
/*---------------------------------------------------------*/
{
        if(fabs(eta)<ETAMAXREGION0) {
                return BASEDETA;
                
        } else if(fabs(eta)<ETAMAXREGION1) {
                return BASEDETA*2.;

        } else if (fabs(eta)<ETAMAXREGION2) {
                return BASEDETA;

        } else {
                return FCALDETA;                
        }
}

/*---------------------------------------------------------*/
double TriggerTowerMon::phiWidth(double eta) const 
/*---------------------------------------------------------*/
{
        if(fabs(eta)<ETAMAXREGION0) {
                return BASEDPHI;
                
        } else if(fabs(eta)<ETAMAXREGION1) {
                return BASEDPHI*2.;

        } else if (fabs(eta)<ETAMAXREGION2) {
                return BASEDPHI*2.;

        } else {
                return BASEDPHI*4.;     
        }
}

/*---------------------------------------------------------*/
int TriggerTowerMon::etaIndex(double eta) const 
/*---------------------------------------------------------*/
{
        int etacenter;       
        double deta = etaWidth(eta);

        if(fabs(eta)<ETAMAXREGION0) {
                etacenter = floor(fabs(eta)/deta);
                
        } else if(fabs(eta)<ETAMAXREGION1) {
	  //int sign = (int) (eta/fabs(eta));
                etacenter = (floor((fabs(eta)-ETAMAXREGION0)/deta));

        } else if (fabs(eta)<ETAMAXREGION2) {
	  //int sign = (int) (eta/fabs(eta));
                etacenter = (floor((fabs(eta)-ETAMAXREGION1)/deta));

        } else {
	  //int sign = (int) (eta/fabs(eta));
                etacenter = (floor((fabs(eta)-ETAMAXREGION2)/deta));
        }
        return etacenter;
}

/*---------------------------------------------------------*/
int TriggerTowerMon::phiIndex(double eta, double phi) const 
/*---------------------------------------------------------*/
{
        double dphi = phiWidth(eta);
        int phicenter = floor(phi/dphi);
        return phicenter;
}

/*---------------------------------------------------------*/
double TriggerTowerMon::IDeta(const Identifier& id) 
/*---------------------------------------------------------*/
{
  int region = m_lvl1Helper->region(id);
  int ieta = m_lvl1Helper->eta(id);
  int sign = m_lvl1Helper->pos_neg_z(id);

  double gran[4] = {0.1, 0.2, 0.1, 0.425};
  double offset[4] = {0., 2.5, 3.1, 3.2};
  double eta;

  if (region>=0 && region<=3) {
    eta = sign* ( ( (ieta+0.5) * gran[region] ) + offset[region] );
  }
  else {
    eta = 0.;
  }
  return eta;
}

/*---------------------------------------------------------*/
double TriggerTowerMon::IDphi(const Identifier& id) 
/*---------------------------------------------------------*/
{
  Identifier regId = m_lvl1Helper->region_id(id);

  double phiMax = m_lvl1Helper->phi_max(regId);
  int iphi = m_lvl1Helper->phi(id);
#ifndef M_PI
  double M_PI = acos (-1.0);
#endif

  double phi = (iphi+0.5)*2.*M_PI/(phiMax+1.);
  return phi;
}
