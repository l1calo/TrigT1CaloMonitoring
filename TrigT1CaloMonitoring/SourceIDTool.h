// ********************************************************************
//
// NAME:        SourceIDTool.h
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
//
// ********************************************************************

#ifndef SourceIDTool_H
#define SourceIDTool_H

#include "GaudiKernel/StatusCode.h"
#include <stdint.h>

#include <map>
#include <string>
#include <vector>

#include "DataModel/DataVector.h"
#include "eformat/SourceIdentifier.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"

#include "TrigT1CaloByteStream/PpmCrateMappings.h"
#include "TrigT1CaloByteStream/JemCrateMappings.h"
#include "TrigT1CaloByteStream/CpmCrateMappings.h"
#include "TrigT1CaloByteStream/L1CaloSrcIdMap.h"


class SourceIDTool 
{
 public:
  
  SourceIDTool(); 
  virtual ~SourceIDTool();
  void SourceIDs();
    
 private:
  /// Number of channels per module (may not all be used)
  int m_channels;
  /// Number of crates
  int m_crates;
  /// Number of modules per crate (may not all exist)
  int m_modules;
  /// Number of slinks per crate when writing out bytestream
  int m_slinks;
  
  /// Sub-detector type
  eformat::SubDetector m_subDetector;
  /// Source ID converter
  LVL1BS::L1CaloSrcIdMap* m_srcIdMap;
  /// crate mappings
  LVL1BS::PpmCrateMappings* m_ppmMaps;
  LVL1BS::CpmCrateMappings* m_cpmMaps;
  LVL1BS::JemCrateMappings* m_jemMaps;
  
};
 

#endif
