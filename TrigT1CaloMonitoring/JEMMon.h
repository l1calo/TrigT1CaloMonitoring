// ********************************************************************
//
// NAME:        JEMMon.h
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of the JEP on JEM level
//
// ********************************************************************

#ifndef JEMMon_H
#define JEMMon_H

#include <string>
#include <vector>

#include "GaudiKernel/ToolHandle.h"

#include "DataModel/DataVector.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"

class TH1F_LW;
class TH2F_LW;
class TH2I_LW;

class StatusCode;

class TrigT1CaloMonErrorTool;
class TrigT1CaloLWHistogramTool;

namespace LVL1 {
  class JEMHits;
  class JEMEtSums;
  class JetElement;
  class JEMRoI;
}


class JEMMon : public ManagedMonitorToolBase
{
public:

   JEMMon( const std::string & type, const std::string & name,
   	                             const IInterface* parent ); 

   virtual ~JEMMon();

   virtual StatusCode initialize();
   virtual StatusCode bookHistograms( bool isNewEventsBlock,
	                              bool isNewLumiBlock, bool isNewRun );
   virtual StatusCode fillHistograms();
   virtual StatusCode procHistograms( bool isEndOfEventsBlock,
	                              bool isEndOfLumiBlock, bool isEndOfRun );

private:

   typedef DataVector<LVL1::JetElement> JECollection;
   typedef DataVector<LVL1::JEMHits> JEMHitsCollection;
   typedef DataVector<LVL1::JEMEtSums> JEMEtSumsCollection;
   typedef DataVector<LVL1::JEMRoI> JemRoiCollection;

   // Tool to retrieve bytestream errors
   ToolHandle<TrigT1CaloMonErrorTool>    m_errorTool;
   ToolHandle<TrigT1CaloLWHistogramTool> m_histTool;

   /** location of data */
   std::string m_JetElementLocation;
   std::string m_JEMHitsLocation;
   std::string m_JEMEtSumsLocation;   
   std::string m_JEMRoILocation;

   int m_SliceNo;
   int m_MaxEnergyRange;

   std::string m_PathInRootFile;   
   std::string m_ErrorPathInRootFile;

   /** Histos */
   TH1F_LW* m_h_je_emeta;
   TH1F_LW* m_h_je_hadeta;
   TH1F_LW* m_h_je_emphi;
   TH1F_LW* m_h_je_hadphi;
   TH1F_LW* m_h_je_emenergy;
   TH1F_LW* m_h_je_hadenergy; 

   // HitMaps
   TH2F_LW* m_h_je_energy_emHitMap;
   TH2F_LW* m_h_je_energy_hadHitMap;
   std::vector<TH2F_LW*> m_h_je_emHitMap;
   std::vector<TH2F_LW*> m_h_je_hadHitMap;

   // error maps
   int m_NoEvents;
   TH2F_LW* m_h_je_error;
   
   // number of triggered slice
   TH1F_LW* m_h_je_triggeredSlice;

   // JEM Hits
   TH1F_LW* m_h_JEMHits_MainHits;
   TH1F_LW* m_h_JEMHits_FwdHitsRight;
   TH1F_LW* m_h_JEMHits_FwdHitsLeft;
   TH2F_LW* m_h_JEMDAQ_Hits_Map;
   // JEM Et Sums
   TH1F_LW*  m_h_JEMEtSums_Ex;
   TH1F_LW*  m_h_JEMEtSums_Ey;
   TH1F_LW*  m_h_JEMEtSums_Et; 

   // JEM RoI
   TH1F_LW* m_h_JEMRoI_MainHits;
   TH1F_LW* m_h_JEMRoI_FwdHitsRight;
   TH1F_LW* m_h_JEMRoI_FwdHitsLeft;
  
   std::vector<TH2F_LW*> m_h_JEMRoI_MainThreshPerEtaPhi;
   std::vector<TH2F_LW*> m_h_JEMRoI_FwdThreshPerEtaPhi;
 
   // errors and saturation
   TH2F_LW* m_h_JEMRoI_error;
   
   // Error Summary
   TH1F_LW* m_h_JEM_ErrorSummary;
   TH2I_LW* m_h_JEM_Events;
};


#endif
