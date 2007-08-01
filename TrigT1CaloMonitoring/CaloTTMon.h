// ********************************************************************
//
// NAME:     CaloTTMon.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//	     
//
// ********************************************************************
#ifndef CALOTTMON_H
#define CALOTTMON_H

#include <map>
#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Units/SystemOfUnits.h"


#include "TH1.h"
#include "TH2.h"

class CaloTTMon: public ManagedMonitorToolBase
{

 public:
  
  CaloTTMon(const std::string & type, const std::string & name,
		  const IInterface* parent);
    

  virtual ~CaloTTMon();

  virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );

private:
  std::string m_CaloTTContainerName;
  int m_CaloTT_HitMap_Thresh0;
  int m_CaloTT_HitMap_Thresh1;
  int m_CaloTT_HitMap_Thresh2;

  std::string m_DataType;
  std::string m_PathInRootFile;
      
protected:
   /** a handle on Store Gate for access to the Event Store */
   StoreGateSvc* m_storeGate;

   // StoreGate service
   StoreGateSvc* m_detStore;

   //CaloTower 
   TH1F* m_h_CaloTower_phi;
   TH1F* m_h_CaloTower_eta;
   TH2F* m_h_CaloTower_HitMap;
   TH1F* m_h_CaloTower_et;

   //CaloCell 
   TH1F* m_h_CaloCell_phi;
   TH1F* m_h_CaloCell_eta;
   TH2F* m_h_CaloCell_HitMap;
   TH1F* m_h_CaloCell_et;

  //LUT Hitmaps per threshold
   TH2F* m_h_CaloTT_HitMap_emLUT_Thresh0;
   TH2F* m_h_CaloTT_HitMap_emLUT_Thresh1;
   TH2F* m_h_CaloTT_HitMap_emLUT_Thresh2;

   TH2F* m_h_CaloTT_HitMap_hadLUT_Thresh0;
   TH2F* m_h_CaloTT_HitMap_hadLUT_Thresh1;
   TH2F* m_h_CaloTT_HitMap_hadLUT_Thresh2;

   //distribution of LUT peak per detector region
   TH1F* m_h_CaloTT_emLUT;
   TH1F* m_h_CaloTT_emLUT_eta;
   TH1F* m_h_CaloTT_emLUT_phi;
   /*TH1F* m_h_CaloTT_emLUT_barrel;
   TH1F* m_h_CaloTT_emLUT_EC;
   TH1F* m_h_CaloTT_emLUT_Quadrant[NoQuadrant];
   TH1F* m_h_CaloTT_emLUT_DetSide[Side];*/

   TH1F* m_h_CaloTT_hadLUT; 
   TH1F* m_h_CaloTT_hadLUT_eta;
   TH1F* m_h_CaloTT_hadLUT_phi;
   /*TH1F* m_h_CaloTT_hadLUT_barrel;
   TH1F* m_h_CaloTT_hadLUT_EC;
   TH1F* m_h_CaloTT_hadLUT_Quadrant[NoQuadrant];
   TH1F* m_h_CaloTT_hadLUT_DetSide[Side];*/
          
};

#endif
