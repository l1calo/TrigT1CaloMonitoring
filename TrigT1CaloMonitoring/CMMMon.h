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

#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "TH1.h"
#include "TH2.h"

#include "DataModel/DataVector.h"
#include "TrigT1Calo/CMMJetHits.h"
#include "TrigT1Calo/CMMEtSums.h"


#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"


class CMMMon : public ManagedMonitorToolBase
{
public:
        typedef DataVector<LVL1::CMMJetHits> CMMJetHitsCollection;
        typedef DataVector<LVL1::CMMEtSums> CMMEtSumsCollection;


	CMMMon( const std::string & type, const std::string & name,
	                 const IInterface* parent ); 

	virtual ~CMMMon();

	virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
	virtual StatusCode fillHistograms();
	virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );

protected:

   /** a handle on Store Gate for access to the Event Store */
   StoreGateSvc* m_storeGate;

   /** location of data */
   std::string m_CMMJetHitsLocation;
   std::string m_CMMEtSumsLocation;
   std::string m_CMMRoILocation;

   std::string m_DataType;   
   std::string m_PathInRootFile;   

  /** Histos */   
   // CMM Jet Hits
   TH1* m_h_CMMJetHits_MainJets;
   TH1* m_h_CMMJetHits_FwdJetsRight;
   TH1* m_h_CMMJetHits_FwdJetsLeft;
   TH1* m_h_CMMJetHits_EtMap;

   // CMM Et Sums
   TH1* m_h_CMMEtSums_Ex;
   TH1* m_h_CMMEtSums_Ey;
   TH1* m_h_CMMEtSums_Et;
   TH1* m_h_CMMEtSums_MissingEtMap;
   TH1* m_h_CMMEtSums_SumEtMap;

   // CMM RoI
   TH1* m_h_CMMRoI_JetEtHits;

};


#endif
