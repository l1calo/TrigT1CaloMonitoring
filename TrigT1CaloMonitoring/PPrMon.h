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

#include <string>
#include <vector>

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "GaudiKernel/ToolHandle.h"

class TH1F;
class TH2F;
class TH2I;
class TProfile;
class TProfile2D;

class StatusCode;

class TrigT1CaloMonErrorTool;
class TrigT1CaloHistogramTool;

namespace LVL1 {
  class IL1TriggerTowerTool;
}

class PPrMon: public ManagedMonitorToolBase
{

 public:
  
  PPrMon(const std::string & type, const std::string & name,
		                   const IInterface* parent);

  virtual ~PPrMon();

  virtual StatusCode initialize();
  virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock,
                                                            bool isNewRun );
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms( bool isEndOfEventsBlock,
                                     bool isEndOfLumiBlock, bool isEndOfRun );

private:

  enum CaloPartitions { LArFCAL1C, LArEMECC, LArOverlapC, LArEMBC, LArEMBA,
                        LArOverlapA, LArEMECA, LArFCAL1A, LArFCAL23C, LArHECC,
			TileEBC, TileLBC, TileLBA, TileEBA, LArHECA, LArFCAL23A,
			MaxPartitions };

  double recTime(const std::vector<int>& vFAdc, int cut);
  int partition(int layer, double eta);
  std::string partitionName(int part);

  std::string m_TriggerTowerContainerName;
  int m_TT_HitMap_ThreshMax;
  int m_TT_HitMap_LumiBlocks;
  int m_TT_ADC_HitMap_Thresh;
  int m_SliceNo;
  int m_MaxEnergyRange;
  int m_TT_ADC_Pedestal;
  int m_HADFADCCut;
  int m_EMFADCCut;
  bool m_onlineTest;

  std::string m_PathInRootFile;
  std::string m_ErrorPathInRootFile;
  std::string m_EventPathInRootFile;
     
  // Tool to retrieve bytestream errors
  ToolHandle<TrigT1CaloMonErrorTool>    m_errorTool;
  ToolHandle<TrigT1CaloHistogramTool>   m_histTool;
  ToolHandle<LVL1::IL1TriggerTowerTool> m_ttTool; 

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
  TH2F* m_h_fwPpmError_Crate_03; //not implemented yet
  TH2F* m_h_fwPpmError_Crate_47; //      "  
  std::vector<TH2F*> m_h_ErrorDetails; //ASIC errors by MCM
  TH2I* m_h_TT_EventNumbers;
  TH2I* m_h_TT_ASICEventNumbers;

  // number of triggered slice
  TH1F* m_h_TT_triggeredSlice_em;
  TH1F* m_h_TT_triggeredSlice_had;
   
  TH1F* m_h_NumberEvents;  
};

#endif
