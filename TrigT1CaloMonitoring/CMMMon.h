// ********************************************************************
//
// NAME:        CMMMon.h
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of the JEP on CMM level
//
// ********************************************************************

#ifndef CMMMon_H
#define CMMMon_H

#include <string>

#include "GaudiKernel/ToolHandle.h"

#include "DataModel/DataVector.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"

class TH1F;
class TH2F;
class TH2I;

class StatusCode;

namespace LVL1 {
  class CMMJetHits;
  class CMMEtSums;
}

class TrigT1CaloMonErrorTool;
class TrigT1CaloHistogramTool;

class CMMMon : public ManagedMonitorToolBase
{
public:

   typedef DataVector<LVL1::CMMJetHits> CMMJetHitsCollection;
   typedef DataVector<LVL1::CMMEtSums>  CMMEtSumsCollection;
  
   CMMMon( const std::string & type, const std::string & name,
	                             const IInterface* parent ); 

   virtual ~CMMMon();

   virtual StatusCode initialize();
   virtual StatusCode bookHistograms( bool isNewEventsBlock,
                                      bool isNewLumiBlock, bool isNewRun );
   virtual StatusCode fillHistograms();
   virtual StatusCode procHistograms( bool isEndOfEventsBlock,
                                      bool isEndOfLumiBlock, bool isEndOfRun );

private:

   // Tool to retrieve bytestream errors
   ToolHandle<TrigT1CaloMonErrorTool> m_errorTool;
   ToolHandle<TrigT1CaloHistogramTool> m_histTool;

   /** location of data */
   std::string m_CMMJetHitsLocation;
   std::string m_CMMEtSumsLocation;
   std::string m_CMMRoILocation;

   std::string m_PathInRootFile;   
   std::string m_ErrorPathInRootFile;

   /** Histos */   
   // CMM Jet Hits
   TH1F* m_h_CMMJetHits_MainJets;
   TH1F* m_h_CMMJetHits_FwdJetsRight;
   TH1F* m_h_CMMJetHits_FwdJetsLeft;
   TH1F* m_h_CMMJetHits_EtMap;
   // JEM Hits
   TH1F* m_h_CMMJetHits_JEM_MainHits;
   TH1F* m_h_CMMJetHits_JEM_FwdHitsRight;
   TH1F* m_h_CMMJetHits_JEM_FwdHitsLeft;
   TH1F* m_h_CMMJetHits_JEM_Crate0ParityError;
   TH1F* m_h_CMMJetHits_JEM_Crate1ParityError;

   // CMM Et Sums
   TH1F* m_h_CMMEtSums_Ex;
   TH1F* m_h_CMMEtSums_Ey;
   TH1F* m_h_CMMEtSums_Et;
   TH1F* m_h_CMMEtSums_MissingEtMap;
   TH1F* m_h_CMMEtSums_SumEtMap;
   // JEM Et Sums
   TH1F*  m_h_CMMEtSums_JEM_Ex;
   TH1F*  m_h_CMMEtSums_JEM_Ey;
   TH1F*  m_h_CMMEtSums_JEM_Et; 

   // CMM RoI
   TH1F* m_h_CMMRoI_JetEtHits;
   TH1F* m_h_CMMRoI_SumEtHits;
   TH1F* m_h_CMMRoI_MissingEtHits;

   TH1F* m_h_CMMRoI_Ex;
   TH1F* m_h_CMMRoI_Ey;
   TH1F* m_h_CMMRoI_Et;

   //errors
   TH2F* m_h_CMMJet_error;
   TH2F* m_h_CMMEnergy_error;
   TH1F* m_h_CMMRoI_error;
   TH1F* m_h_CMM_ErrorSummary;
   TH1F* m_h_TriggeredSlice;	  
   TH2I* m_h_CMM_Events;
};


#endif
