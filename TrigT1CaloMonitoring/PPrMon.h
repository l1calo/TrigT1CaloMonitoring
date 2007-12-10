// ********************************************************************
//
// NAME:     PPrMon.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Ethan-Etienne Woehrling (eew@hep.ph.bham.ac.uk)
//           Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//	     
//
// ********************************************************************
#ifndef PPRMON_H
#define PPRMON_H

#include <map>
#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CaloIdentifier/CaloIdManager.h"
#include "CaloIdentifier/CaloLVL1_ID.h"
#include "Identifier/Identifier.h"
#include "TrigT1CaloCalibTools/L1CaloTTIdTools.h"
#include "CaloTriggerTool/CaloTriggerTowerService.h"

#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
//class StoreGateSvc;

class PPrMon: public ManagedMonitorToolBase
{

 public:
  
  PPrMon(const std::string & type, const std::string & name,
		  const IInterface* parent);
    

  virtual ~PPrMon();

  virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );

private:

  double recTime(const std::vector<int>& vFAdc);
  int FADCSum(const std::vector<int>& vFAdc) ;

  std::string m_TriggerTowerContainerName;
  int m_TT_HitMap_Thresh0;
  int m_TT_HitMap_Thresh1;
  int m_TT_HitMap_Thresh2;
  int m_TT_ADC_HitMap_Thresh;
  int m_SliceNo;
  bool m_TT_DistPerChannel;
  bool m_TT_DistPerChannelAndTimeSlice;
  int m_MaxEnergyRange;
  bool m_Offline;
  int m_TT_ADC_Pedestal;
  bool m_TT_ADCTimingPerChannel;
  int m_HADFADCCut;
  int m_EMFADCCut;

  std::string m_DataType;
  std::string m_PathInRootFile;
  std::string m_ErrorPathInRootFile;
  std::string m_EventPathInRootFile;
     
protected:
   /** a handle on Store Gate for access to the Event Store */
   StoreGateSvc* m_storeGate;

   // StoreGate service
   StoreGateSvc* m_detStore;
   // Calorimeter Id manager
   const CaloIdManager* m_caloMgr;
   // CaloLVL1_ID Id helper
   const CaloLVL1_ID* m_lvl1Helper;
   const L1CaloTTIdTools* m_l1CaloTTIdTools;
   
   CaloTriggerTowerService* m_ttSvc;
   // TTOnlineID Id helper
   const TTOnlineID* m_l1ttonlineHelper;
   

   // histos per channel
  std::map <Identifier, TH1F*>  m_h_TT_EmADCPeak;
  std::map <Identifier, TH1F*>  m_h_TT_HadADCPeak;
  std::map <Identifier, TH1F*>  m_h_TT_EmLUTPeak;
  std::map <Identifier, TH1F*>  m_h_TT_HadLUTPeak;
  std::map <int,TProfile*> m_h_TT_HitMap_emADCChannel_timing;
  std::map <int,TProfile*> m_h_TT_HitMap_hadADCChannel_timing;

  //ADC Hitmaps per TimeSlice
  std::map <int,TH2F*> m_h_TT_HitMap_emADC;
  std::map <int,TH2F*> m_h_TT_HitMap_hadADC;

  //timing HitMaps
  TProfile2D* m_h_TT_ADC_emTiming_signal;
  TProfile2D* m_h_TT_ADC_hadTiming_signal;

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

   // error
   int m_NoEvents;
   TH1F* m_h_TT_emerror;
   TH2F* m_h_TT_error_Crate_03;
   TH2F* m_h_TT_error_Crate_47;
   TH1F* m_h_TT_haderror;
   TH2F* m_h_TT_em_GLinkDown;
   TH2F* m_h_TT_em_GLinkTimeout;
   TH2F* m_h_TT_had_GLinkDown;
   TH2F* m_h_TT_had_GLinkTimeout;

   // number of triggered slice
   TH1F* m_h_TT_triggeredSlice_em;
   TH1F* m_h_TT_triggeredSlice_had;
   
   TH1F* m_h_NumberEvents;  
};

#endif
