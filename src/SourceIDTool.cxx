#include <sstream>
#include <algorithm>
#include <math.h>
#include <functional>
#include <iostream>

#include <numeric>
#include <utility>
#include "TrigT1CaloMonitoring/SourceIDTool.h"

#include "GaudiKernel/IInterface.h"

#include "TrigT1Calo/DataError.h"
#include "TrigT1Calo/TriggerTower.h"
#include "TrigT1Calo/TriggerTowerKey.h"

#include "TrigT1CaloByteStream/ChannelCoordinate.h"
#include "TrigT1CaloByteStream/L1CaloRodStatus.h"
#include "TrigT1CaloByteStream/L1CaloSubBlock.h"
#include "TrigT1CaloByteStream/L1CaloUserHeader.h"
#include "TrigT1CaloByteStream/ModifySlices.h"
#include "TrigT1CaloByteStream/PpmByteStreamTool.h"
#include "TrigT1CaloByteStream/PpmCrateMappings.h"
#include "TrigT1CaloByteStream/PpmSubBlock.h"


// Constructor
/*---------------------------------------------------------*/
SourceIDTool::SourceIDTool() 
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
SourceIDTool::~SourceIDTool()
/*---------------------------------------------------------*/
{
}


/*---------------------------------------------------------*/
void SourceIDTool::SourceIDs() 
/*---------------------------------------------------------*/
{
  std::cout<<std::endl;
  std::cout<<std::endl;
  //---------------------------------------------------------
  //------------------ PPr ----------------------------------
  //---------------------------------------------------------
  m_subDetector = eformat::TDAQ_CALO_PREPROC;
  m_srcIdMap    = new LVL1BS::L1CaloSrcIdMap();
  m_ppmMaps     = new LVL1BS::PpmCrateMappings();
  m_crates      = m_ppmMaps->crates();
  m_modules     = m_ppmMaps->modules();
  m_channels    = m_ppmMaps->channels();

  const int maxlinks = m_srcIdMap->maxSlinks();
  std::cout <<"Preprozessor ROB IDs:"<<std::endl;
  for (int crate = 0; crate < m_crates; ++crate) {
    for (int slink = 0; slink < maxlinks; ++slink) {
      const int daqOrRoi = 0;
      const uint32_t rodId = m_srcIdMap->getRodID(crate, slink, daqOrRoi,  m_subDetector);
      std::cout <<std::hex<<"0x"<<rodId<<std::endl;
    }
  } 

  std::cout<<std::endl;
  //---------------------------------------------------------
  //------------------ CP DAQ -------------------------------
  //---------------------------------------------------------
  int m_crateOffsetHw  = 8;

  m_subDetector = eformat::TDAQ_CALO_CLUSTER_PROC_DAQ;
  m_cpmMaps     = new LVL1BS::CpmCrateMappings();
  m_crates      = m_cpmMaps->crates();
  m_modules     = m_cpmMaps->modules();
  m_channels    = m_cpmMaps->channels();

  std::cout <<"CP DAQ ROB IDs:"<<std::endl;
  int maxCrates = m_crates + m_crateOffsetHw;
  int maxSlinks = m_srcIdMap->maxSlinks();
  for (int hwCrate = m_crateOffsetHw; hwCrate < maxCrates; ++hwCrate) {
    for (int slink = 0; slink < maxSlinks; ++slink) {
      const int daqOrRoi = 0;
      const uint32_t rodId = m_srcIdMap->getRodID(hwCrate, slink, daqOrRoi,
                                                             m_subDetector);
      std::cout <<std::hex<<"0x"<<rodId<<std::endl;
   }
  } 
  
  std::cout<<std::endl;
  //---------------------------------------------------------
  //------------------ CP RoI -------------------------------
  //---------------------------------------------------------
  
  m_subDetector = eformat::TDAQ_CALO_CLUSTER_PROC_ROI;
  m_crates      = LVL1BS::CpmCrateMappings::crates();
  m_modules     = LVL1BS::CpmCrateMappings::modules();

  std::cout <<"CP RoI ROB IDs:"<<std::endl;
  maxCrates = m_crates + m_crateOffsetHw;
  maxSlinks = m_srcIdMap->maxSlinks();
  for (int hwCrate = m_crateOffsetHw; hwCrate < maxCrates; ++hwCrate) {
    for (int slink = 0; slink < maxSlinks; ++slink) {
      const int daqOrRoi = 1;
      const uint32_t rodId = m_srcIdMap->getRodID(hwCrate, slink, daqOrRoi,
                                                             m_subDetector);
      std::cout <<std::hex<<"0x"<<rodId<<std::endl;
   }
  }
  
  std::cout<<std::endl;
  //---------------------------------------------------------
  //------------------ JEP DAQ ------------------------------
  //---------------------------------------------------------
  m_crateOffsetHw  = 12;
  
  m_subDetector = eformat::TDAQ_CALO_JET_PROC_DAQ;
  m_jemMaps     = new LVL1BS::JemCrateMappings();
  m_crates      = m_jemMaps->crates();
  m_modules     = m_jemMaps->modules();
  m_channels    = m_jemMaps->channels();
  
  std::cout <<"JEP DAQ ROB IDs:"<<std::endl;
  maxCrates = m_crates + m_crateOffsetHw;
  maxSlinks = m_srcIdMap->maxSlinks();
  for (int hwCrate = m_crateOffsetHw; hwCrate < maxCrates; ++hwCrate) {
    for (int slink = 0; slink < maxSlinks; ++slink) {
      const int daqOrRoi = 0;
      const uint32_t rodId = m_srcIdMap->getRodID(hwCrate, slink, daqOrRoi,
                                                             m_subDetector);
      std::cout <<std::hex<<"0x"<<rodId<<std::endl;
   }
  }
  
  std::cout<<std::endl;
 //---------------------------------------------------------
  //------------------ JEP RoI ------------------------------
  //---------------------------------------------------------
  m_subDetector = eformat::TDAQ_CALO_JET_PROC_ROI;
  m_crates      = LVL1BS::JemCrateMappings::crates();
  m_modules     = LVL1BS::JemCrateMappings::modules();

  std::cout <<"JEP RoI ROB IDs:"<<std::endl;
  maxCrates = m_crates + m_crateOffsetHw;
  maxSlinks = m_srcIdMap->maxSlinks();
  for (int hwCrate = m_crateOffsetHw; hwCrate < maxCrates; ++hwCrate) {
    for (int slink = 0; slink < maxSlinks; ++slink) {
      const int daqOrRoi = 1;
      const uint32_t rodId = m_srcIdMap->getRodID(hwCrate, slink, daqOrRoi,
						  m_subDetector);
 
      std::cout <<std::hex<<"0x"<<rodId<<std::endl;
    }
  } 
  std::cout<<std::endl;
  std::cout<<std::endl;
  

}


