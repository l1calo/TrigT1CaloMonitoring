#ifndef TRIGT1CALOBYTESTREAM_TRIGT1CALOMONERRORTOOL_H
#define TRIGT1CALOBYTESTREAM_TRIGT1CALOMONERRORTOOL_H

#include <string>
#include <vector>

#include "AthenaBaseComps/AthAlgTool.h"

class IInterface;
class InterfaceID;
class StatusCode;

/** Tool to retrieve ROB status and ROD unpacking errors from StoreGate.
 *
 *  Forces retrieve of all containers to ensure error vector is complete.
 *
 *  @author Peter Faulkner
 */

class TrigT1CaloMonErrorTool : public AthAlgTool {

 public:
   TrigT1CaloMonErrorTool(const std::string& type, const std::string& name,
                             const IInterface* parent);
   virtual ~TrigT1CaloMonErrorTool();

   /// AlgTool InterfaceID
   static const InterfaceID& interfaceID();

   virtual StatusCode initialize();
   virtual StatusCode finalize();

   /// Retrieve error vector
   StatusCode retrieve(const std::vector<unsigned int>*& errColl);
   /// Return true if current event has any corruption errors
   bool corrupt();

 private:

   /// Trigger Tower container StoreGate key
   std::string m_triggerTowerLocation;
   /// CPM core tower container StoreGate key
   std::string m_cpmTowerLocation;
   /// CPM overlap tower container StoreGate key
   std::string m_cpmTowerLocationOverlap;
   /// CPM hits container StoreGate key
   std::string m_cpmHitsLocation;
   /// CMM-CP hits container StoreGate key
   std::string m_cmmCpHitsLocation;
   /// CPM RoI container StoreGate key
   std::string m_cpmRoiLocation;
   /// Core Jet Element container StoreGate key
   std::string m_jetElementLocation;
   /// Overlap Jet Element container StoreGate key
   std::string m_jetElementLocationOverlap;
   /// JEM hits container StoreGate key
   std::string m_jemHitsLocation;
   /// CMM-Jet hits container StoreGate key
   std::string m_cmmJetHitsLocation;
   /// JEM RoI container StoreGate key
   std::string m_jemRoiLocation;
   /// CMM RoI container StoreGate key
   std::string m_cmmRoiLocation;
   /// JEM Et sums container StoreGate key
   std::string m_jemEtSumsLocation;
   /// CMM Et sums container StoreGate key
   std::string m_cmmEtSumsLocation;
   /// ROD header container StoreGate key
   std::string m_rodHeaderLocation;
   /// CP RoIB ROD header container StoreGate key
   std::string m_cpRoibRodHeaderLocation;
   /// JEP RoIB ROD header container StoreGate key
   std::string m_jepRoibRodHeaderLocation;
   /// ROB and Unpacking Error vector StoreGate key
   std::string m_robErrorVectorLocation;
   /// Flag corrupt events
   bool m_flagCorruptEvents;

};

#endif
