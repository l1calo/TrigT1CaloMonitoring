// ********************************************************************
//
// NAME:     TriggerTowerMon.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Ethan-Etienne Woehrling (eew@hep.ph.bham.ac.uk)
//           Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//	     
//
// ********************************************************************
#ifndef TRIGGERTOWERMON_H
#define TRIGGERTOWERMON_H

#include <map>
#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CaloIdentifier/CaloIdManager.h"
#include "CaloIdentifier/CaloLVL1_ID.h"
#include "Identifier/Identifier.h"
#include "TrigT1CaloCalibTools/L1CaloTTIdTools.h"


#include "TH1.h"
#include "TH2.h"
//class StoreGateSvc;

class TriggerTowerMon: public ManagedMonitorToolBase
{

 public:
  
  TriggerTowerMon(const std::string & type, const std::string & name,
		  const IInterface* parent);
    

  virtual ~TriggerTowerMon();

  virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );

private:
  std::string m_TriggerTowerContainerName;
  int m_TT_HitMap_Thresh0;
  int m_TT_HitMap_Thresh1;
  int m_TT_HitMap_Thresh2;
  int m_TT_ADC_HitMap_Thresh;
  int m_SliceNo;
  bool m_TT_DistPerChannel;
  bool m_TT_DistPerChannelAndTimeSlice;

  std::string m_DataType;
  std::string m_PathInRootFile;
      
protected:
   /** a handle on Store Gate for access to the Event Store */
   StoreGateSvc* m_storeGate;

   // StoreGate service
   StoreGateSvc* m_detStore;;
   // Calorimeter Id manager
   const CaloIdManager* m_caloMgr;
   // CaloLVL1_ID Id helper
   const CaloLVL1_ID* m_lvl1Helper;
   const L1CaloTTIdTools* m_l1CaloTTIdTools;

   // histos per channel
  std::map <Identifier, TH1F*>  m_h_TT_EmADCPeak;
  std::map <Identifier, TH1F*>  m_h_TT_HadADCPeak;
  std::map <Identifier, TH1F*>  m_h_TT_EmLUTPeak;
  std::map <Identifier, TH1F*>  m_h_TT_HadLUTPeak;

  //ADC Hitmaps per TimeSlice
  std::map <int,TH2F*> m_h_TT_HitMap_emADC;
  std::map <int,TH2F*> m_h_TT_HitMap_hadADC;

  //LUT Hitmaps per threshold
   TH2F* m_h_TT_HitMap_emLUT_Thresh0;
   TH2F* m_h_TT_HitMap_emLUT_Thresh1;
   TH2F* m_h_TT_HitMap_emLUT_Thresh2;

   TH2F* m_h_TT_HitMap_hadLUT_Thresh0;
   TH2F* m_h_TT_HitMap_hadLUT_Thresh1;
   TH2F* m_h_TT_HitMap_hadLUT_Thresh2;

   //distribution of LUT peak per detector region
   TH1F* m_h_TT_emLUT;
   TH1F* m_h_TT_emLUT_eta;
   TH1F* m_h_TT_emLUT_phi;
   /*TH1F* m_h_TT_emLUT_barrel;
   TH1F* m_h_TT_emLUT_EC;
   TH1F* m_h_TT_emLUT_Quadrant[NoQuadrant];
   TH1F* m_h_TT_emLUT_DetSide[Side];*/

   TH1F* m_h_TT_hadLUT; 
   TH1F* m_h_TT_hadLUT_eta;
   TH1F* m_h_TT_hadLUT_phi;
   /*TH1F* m_h_TT_hadLUT_barrel;
   TH1F* m_h_TT_hadLUT_EC;
   TH1F* m_h_TT_hadLUT_Quadrant[NoQuadrant];
   TH1F* m_h_TT_hadLUT_DetSide[Side];*/

   // number of triggered slice
   TH1F* m_h_TT_triggeredSlice_em;
   TH1F* m_h_TT_triggeredSlice_had;
     
};

#endif
