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

#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "TH1.h"
#include "TH2.h"

#include "DataModel/DataVector.h"
#include "TrigT1Calo/JEMHits.h"
#include "TrigT1Calo/JEMEtSums.h"

#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"

namespace LVL1 {
  class JEMRoI;
}


class JEMMon : public ManagedMonitorToolBase
{
public:
        typedef DataVector<LVL1::JEMHits> JEMHitsCollection;
	typedef DataVector<LVL1::JEMEtSums> JEMEtSumsCollection;
	typedef DataVector<LVL1::JEMRoI> JemRoiCollection;
       

	JEMMon( const std::string & type, const std::string & name,
	                 const IInterface* parent ); 

	virtual ~JEMMon();

	virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
	virtual StatusCode fillHistograms();
	virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );

protected:

   /** a handle on Store Gate for access to the Event Store */
   StoreGateSvc* m_storeGate;

   /** location of data */
   std::string m_JEMHitsLocation;
   std::string m_JEMEtSumsLocation;   
   std::string m_JEMRoILocation;

   std::string m_DataType;   
   std::string m_PathInRootFile;   
   std::string m_ErrorPathInRootFile;

   /** Histos */

   // JEM Hits
   TH1F* m_h_JEMHits_MainHits;
   TH1F* m_h_JEMHits_FwdHitsRight;
   TH1F* m_h_JEMHits_FwdHitsLeft;

   // JEM Et Sums
   TH1F*  m_h_JEMEtSums_Ex;
   TH1F*  m_h_JEMEtSums_Ey;
   TH1F*  m_h_JEMEtSums_Et; 

   // JEM RoI
   TH1F* m_h_JEMRoI_MainHits;
   TH1F* m_h_JEMRoI_FwdHitsRight;
   TH1F* m_h_JEMRoI_FwdHitsLeft;

   std::map <int, TH2F*> m_h_JEMRoI_MainThreshPerEtaPhi;
   std::map <int, TH2F*> m_h_JEMRoI_FwdThreshPerEtaPhi;
 
   // errors and saturation
   TH2F* m_h_JEMRoI_error;
};


#endif
