// ********************************************************************
//
// NAME:     PPrSpareMon.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Peter Faulkner
//	     
//
// ********************************************************************
#ifndef PPRSPAREMON_H
#define PPRSPAREMON_H

#include <string>
#include <vector>

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "GaudiKernel/ToolHandle.h"

class TH1F;
class TH2F;
class TProfile2D;
class StatusCode;

class TrigT1CaloMonErrorTool;
class TrigT1CaloHistogramTool;

class PPrSpareMon: public ManagedMonitorToolBase
{

 public:
  
  PPrSpareMon(const std::string & type, const std::string & name,
 		                        const IInterface* parent);
    
  virtual ~PPrSpareMon();

  virtual StatusCode initialize();
  virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock,
                                                            bool isNewRun );
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms( bool isEndOfEventsBlock,
                                     bool isEndOfLumiBlock, bool isEndOfRun );

private:

  std::string m_TriggerTowerContainerName;
  int m_TT_ADC_HitMap_Thresh;
  int m_SliceNo;
  bool m_onlineTest;

  std::string m_PathInRootFile;
  std::string m_ErrorPathInRootFile;
  std::string m_EventPathInRootFile;
     
  // Tool to retrieve bytestream errors
  ToolHandle<TrigT1CaloMonErrorTool> m_errorTool;
  ToolHandle<TrigT1CaloHistogramTool> m_histTool;

  //ADC Hitmaps per TimeSlice
  TH2F* m_h_TT_HitMap_ADC;
  TProfile2D* m_p_TT_HitMap_ADC;

  // error
  TH1F* m_h_TT_Error; // Error summary
  TH2F* m_h_TT_error_Crate_25; // just ROD sub-status word
  TH2F* m_h_fwPpmError_Crate_25;
  std::vector<TH2F*> m_h_ErrorDetails; //ASIC errors by MCM
  TH2I* m_h_TT_EventNumbers;
  TH2I* m_h_TT_ASICEventNumbers;

  // number of triggered slice
  TH1F* m_h_TT_triggeredSlice;
   
};

#endif
