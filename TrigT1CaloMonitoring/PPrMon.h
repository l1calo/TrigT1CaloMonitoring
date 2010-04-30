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
#include "GaudiKernel/ToolHandle.h"
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

class TrigT1CaloMonErrorTool;
class TrigT1CaloHistogramTool;

class PPrMon: public ManagedMonitorToolBase
{

 public:
  
  PPrMon(const std::string & type, const std::string & name,
		  const IInterface* parent);
    

  virtual ~PPrMon();

  virtual StatusCode initialize();
  virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );

private:

  enum CaloPartitions { LArFCAL1C, LArEMECC, LArOverlapC, LArEMBC, LArEMBA, LArOverlapA, LArEMECA, LArFCAL1A,
                        LArFCAL23C, LArHECC, TileEBC, TileLBC, TileLBA, TileEBA, LArHECA, LArFCAL23A, MaxPartitions };

  double recTime(const std::vector<int>& vFAdc);
  int FADCSum(const std::vector<int>& vFAdc) ;
  int partition(int layer, double eta);
  std::string partitionName(int part);

  std::string m_TriggerTowerContainerName;
  int m_TT_HitMap_Thresh0;
  int m_TT_HitMap_Thresh1;
  int m_TT_HitMap_Thresh2;
  int m_TT_HitMap_ThreshMax;
  int m_TT_HitMap_LumiBlocks;
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
  bool m_onlineTest;
  



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

   // Tool to retrieve bytestream errors
   ToolHandle<TrigT1CaloMonErrorTool> m_errorTool;
   ToolHandle<TrigT1CaloHistogramTool> m_histTool;
   

   // histos per channel
  std::map <Identifier, TH1F*>  m_h_TT_EmADCPeak;
  std::map <Identifier, TH1F*>  m_h_TT_HadADCPeak;
  std::map <Identifier, TH1F*>  m_h_TT_EmLUTPeak;
  std::map <Identifier, TH1F*>  m_h_TT_HadLUTPeak;
  std::map <int,TProfile*> m_h_TT_HitMap_emADCChannel_timing;
  std::map <int,TProfile*> m_h_TT_HitMap_hadADCChannel_timing;

  //ADC Hitmaps per TimeSlice
  
  /* reducing time slices, only time slice 00100
  std::map <int,TH2F*> m_h_TT_HitMap_emADC;
  std::map <int,TH2F*> m_h_TT_HitMap_hadADC;
  */

  TH2F* m_h_TT_HitMap_emADC_00100;
  TH2F* m_h_TT_HitMap_hadADC_00100;
  TH1F* m_h_dist_had_max;
  TH1F* m_h_dist_em_max;

  TProfile2D* m_p_TT_HitMap_emADC_00100;
  TProfile2D* m_p_TT_HitMap_hadADC_00100;

  //timing HitMaps
  TProfile2D* m_h_TT_ADC_emTiming_signal;
  TProfile2D* m_h_TT_ADC_hadTiming_signal;
  
  std::vector<TProfile*> m_h_TT_SignalProfile;

  //LUT Hitmaps per threshold
  std::vector<TH2F*> m_h_TT_HitMap_emLUT_Thresh;
  std::vector<TH2F*> m_h_TT_HitMap_hadLUT_Thresh;
  TProfile2D* m_p_TT_HitMap_emLUT_etAv;
  TProfile2D* m_p_TT_HitMap_hadLUT_etAv;

   //distribution of LUT peak per detector region
   TH1F* m_h_TT_emLUT;
   TH1F* m_h_TT_emLUT_eta;
   TH1F* m_h_TT_emLUT_phi;
   

   TH1F* m_h_TT_hadLUT; 
   TH1F* m_h_TT_hadLUT_eta;
   TH1F* m_h_TT_hadLUT_phi;

   TH1F* m_h_TT_BCLUT;
   TH2F* m_h_TT_BCID;
   

   // error
   int m_NoEvents;
   TH1F* m_h_TT_Error; // Error summary
   TH2F* m_h_TT_error_Crate_03; // just ROD sub-status word
   TH2F* m_h_TT_error_Crate_47; //just ROD sub-status word
   TH2F* m_h_BCNmis_Crate_03; //expert: BCN mismatch
   TH2F* m_h_BCNmis_Crate_47; //expert: BCN mismatch
   TH2F* m_h_fwPpmError_Crate_03; //not implemented yet
   TH2F* m_h_fwPpmError_Crate_47; //      "  
   std::vector<TH2F*> m_h_ErrorDetails; //ASIC errors by MCM


   // number of triggered slice
   TH1F* m_h_TT_triggeredSlice_em;
   TH1F* m_h_TT_triggeredSlice_had;
   
 
   TH1F* m_h_NumberEvents;  
};

#endif
