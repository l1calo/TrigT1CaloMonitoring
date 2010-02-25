// ********************************************************************
//
// NAME:     PPrSpareMon.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Ethan-Etienne Woehrling (eew@hep.ph.bham.ac.uk)
//           Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//	     
//
// ********************************************************************
#ifndef PPRSPAREMON_H
#define PPRSPAREMON_H

#include <map>
#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/ToolHandle.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
class StoreGateSvc;

class TrigT1CaloMonErrorTool;

class PPrSpareMon: public ManagedMonitorToolBase
{

 public:
  
  PPrSpareMon(const std::string & type, const std::string & name,
		  const IInterface* parent);
    

  virtual ~PPrSpareMon();

  virtual StatusCode initialize();
  virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );

private:

  void setHitmapLabels(TH2* hist);

  std::string m_TriggerTowerContainerName;
  int m_TT_ADC_HitMap_Thresh;
  int m_SliceNo;
  bool m_onlineTest;
  


  std::string m_PathInRootFile;
  std::string m_ErrorPathInRootFile;
  std::string m_EventPathInRootFile;
     
protected:
   /** a handle on Store Gate for access to the Event Store */
   StoreGateSvc* m_storeGate;

   // Tool to retrieve bytestream errors
   ToolHandle<TrigT1CaloMonErrorTool> m_errorTool;
   

  //ADC Hitmaps per TimeSlice
  
  TH2F* m_h_TT_HitMap_ADC;

  TProfile2D* m_p_TT_HitMap_ADC;

   // error
   TH1F* m_h_TT_Error; // Error summary
   TH2F* m_h_TT_error_Crate_25; // just ROD sub-status word
   TH2F* m_h_fwPpmError_Crate_25;
   std::vector<TH2F*> m_h_ErrorDetails; //ASIC errors by MCM


   // number of triggered slice
   TH1F* m_h_TT_triggeredSlice;
   
};

#endif
