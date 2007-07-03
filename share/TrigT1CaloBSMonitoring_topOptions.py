
####################################################################################################################
##################################### Bytestream Optionen ##########################################################
####################################################################################################################

#theApp.TopAlg  = [ ]
 
include( "TriggerRelease/TrigReadBS_topOptions.py" )

include ( "TrigT1CaloByteStream/ReadPpmBS_jobOptions.py" )
ByteStreamAddressProviderSvc = Service( "ByteStreamAddressProviderSvc" )
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::TriggerTower>/BS_TriggerTowers" ]

#include ( "TrigT1CaloByteStream/TestReadJepBS_jobOptions.py" )
include ( "TrigT1CaloByteStream/ReadJepBS_jobOptions.py" )
ByteStreamAddressProviderSvc = Service( "ByteStreamAddressProviderSvc" )
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::JetElement>/BS_JetElements" ]
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::JEMHits>/BS_JEMHits" ]
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::JEMEtSums>/BS_JEMEtSums" ]
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::CMMJetHits>/BS_CMMJetHits" ]
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::CMMEtSums>/BS_CMMEtSums" ]

include ( "TrigT1CaloByteStream/ReadJepRoiBS_jobOptions.py" )
ByteStreamAddressProviderSvc = Service( "ByteStreamAddressProviderSvc" )
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::JEMRoI>/BS_JEMRoIs" ]
ByteStreamAddressProviderSvc.TypeNames += [ "LVL1::CMMRoI/BS_CMMRoIs" ]

include ( "TrigT1CaloByteStream/ReadCpBS_jobOptions.py" )
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::CPMTower>/BS_CPMTowers" ]
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::CPMHits>/BS_CPMHits" ]
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::CMMCPHits>/BS_CMMCPHits" ]

include ( "TrigT1CaloByteStream/ReadCpmRoiBS_jobOptions.py" )
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::CPMRoI>/BS_CPMRoIs" ]

####################################################################################################################
########################################## L1Calo Simulation #######################################################
####################################################################################################################

#--------------------------------------------------------------
# TrigT1Calo Algorithms Private Options

#--------------------------------------------------------------

#load relevant libraries
theApp.Dlls += [ "TrigT1Calo" ]

#services
theApp.CreateSvc += [ "LVL1::TrigT1CaloInfoLoader/TrigT1CaloInfoLoader" ]
TrigT1CaloInfoLoader = Service( "TrigT1CaloInfoLoader" )
TrigT1CaloInfoLoader.TrigT1CaloInfoName = "TrigT1CaloInfo"

# TrigT1Calo : create trigger towers and run trigger sim.
theApp.TopAlg += ["LVL1::CPMTowerMaker/CPMTowerMaker" ]
theApp.TopAlg += ["LVL1::JetElementMaker/JetElementMaker"]       
theApp.TopAlg += ["LVL1::EmTauTrigger/EmTauTrigger" ]
theApp.TopAlg += ["LVL1::EnergyTrigger/EnergyTrigger" ]          
theApp.TopAlg += ["LVL1::JetTrigger/JetTrigger" ]         
theApp.TopAlg += ["LVL1::ROD/ROD" ]          
theApp.TopAlg += ["LVL1::JEPCMMMaker/JEPCMMMaker" ]
theApp.TopAlg += ["LVL1::CPCMMMaker/CPCMMMaker" ]

#*************************************************************
# Set input and output locations  
#*************************************************************
# CPMTowerMaker
CPMTowerMaker = Algorithm( "CPMTowerMaker" )
CPMTowerMaker.TriggerTowerLocation = "BS_TriggerTowers"
CPMTowerMaker.CPMTowerLocation = "Sim_CPMTowers"

# JetElementMaker
JetElementMaker = Algorithm( "JetElementMaker" )
JetElementMaker.TriggerTowerLocation ="BS_TriggerTowers"  
JetElementMaker.JetElementLocation ="Sim_JetElements"

# EmTauTrigger
EmTauTrigger = Algorithm( "EmTauTrigger" )
EmTauTrigger.TriggerTowerLocation = "BS_TriggerTowers"
EmTauTrigger.CPMHitsLocation = "Sim_CPMHits"

# JetTrigger
JetTrigger = Algorithm( "JetTrigger" )
JetTrigger.JetElementLocation="BS_JetElements"
JetTrigger.JEMHitsLocation = "Sim_JEMHits"

# EnergyTrigger
EnergyTrigger = Algorithm( "EnergyTrigger" )
EnergyTrigger.JetElementLocation="BS_JetElements"
EnergyTrigger.JEMEtSumsLocation="Sim_JEMEtSums"

# JEP CMMs
JEPCMMMaker = Algorithm( "JEPCMMMaker" )
JEPCMMMaker.JetElementLocation="BS_JetElements"
JEPCMMMaker.JEMHitsLocation = "BS_JEMHits"
JEPCMMMaker.JEMEtSumsLocation = "BS_JEMEtSums"
JEPCMMMaker.CMMJetHitsLocation = "Sim_CMMJetHits"
JEPCMMMaker.CMMEtSumsLocation = "Sim_CMMEtSums"
JEPCMMMaker.JEMRoILocation = "Sim_JEMRoIs"
JEPCMMMaker.CMMRoILocation = "Sim_CMMRoIs"
#JEPCMMMaker.OutputLevel = VERBOSE

# CP CMMs
CPCMMMaker = Algorithm( "CPCMMMaker" )
CPCMMMaker.CPMTowerLocation = "BS_CPMTowers"
CPCMMMaker.CPMHitsLocation = "BS_CPMHits"
CPCMMMaker.CMMCPHitsLocation = "Sim_CMMCPHits"
CPCMMMaker.CPMRoILocation = "Sim_CPMRoIs"

# ROD
ROD = Algorithm( "ROD" )
#ROD.JetEtRoILocation = "Sim_JetEtROIs"


include ("TrigT1CaloMonitoring/TrigT1CaloConfiguration.py")

####################################################################################################################
##################################### Monitoring ###################################################################
####################################################################################################################


# *************************
# Application configuration
# *************************

#theApp.setup( MONTECARLO )
theApp.Dlls   += [ "AthenaMonitoring" ]
theApp.Dlls   += [ "TrigT1CaloMonitoring" ]
theApp.TopAlg += [ "AthenaMonManager/PrimaryManager" ]

# ************************
# Monitoring configuration
# ************************

## Setup the output file(s):
THistSvc = Algorithm( "THistSvc" )

## The string "TestMon" in the argument below is the 'FileKey'
## used by Athena to access the output file internally
THistSvc.Output = ["CaloBSMon DATAFILE='CaloBSMon.root' OPT='RECREATE'"]

## AthenaMonManager is the Algorithm that manages many classes inheriting
## from ManagedMonitorToolBase
monMan = Algorithm( "PrimaryManager" )

## Add all the ManagedMonitorToolBase objects
monMan.AthenaMonTools += [ "TriggerTowerMon/L1TriggerTowerMonTool" ]


monMan.AthenaMonTools += [ "JetElementMon/BS_L1JetElementMonTool" ]
monMan.AthenaMonTools += [ "JetElementMon/Sim_L1JetElementMonTool" ]

monMan.AthenaMonTools += [ "JEMMon/BS_L1JEMMonTool" ]
monMan.AthenaMonTools += [ "JEMMon/Sim_L1JEMMonTool" ]

monMan.AthenaMonTools += [ "CMMMon/BS_L1CMMMonTool" ]
monMan.AthenaMonTools += [ "CMMMon/Sim_L1CMMMonTool" ]

monMan.AthenaMonTools += [ "SimBSMon/SimBSMonTool" ]

#monMan.AthenaMonTools += [ "TrigT1CaloBSMonTool/L1CaloBSTool" ]
monMan.AthenaMonTools += [ "TrigT1CaloCpmMonTool/L1BSCPMMonTool" ]
monMan.AthenaMonTools += [ "TrigT1CaloCpmMonTool/L1SimCPMMonTool" ]

## get a handle on the ToolSvc
#from AthenaCommon.AppMgr import ToolSvc as toolSvc
ToolSvc = Algorithm( "ToolSvc" )

####################### TriggerTower ################################
#ToolSvc.L1TriggerTowerMonTool.DataType = "M3"  
ToolSvc.L1TriggerTowerMonTool.BS_TriggerTowerContainer = "BS_TriggerTowers"
ToolSvc.L1TriggerTowerMonTool.PathInRootFile = "Stats/TriggerTower"
#ToolSvc.L1TriggerTowerMonTool.OutputLevel = DEBUG

####################### JetElements ################################
ToolSvc.BS_L1JetElementMonTool.DataType = "BS"  #BS or Sim data?
ToolSvc.BS_L1JetElementMonTool.JetElementLocation = "BS_JetElements"
ToolSvc.BS_L1JetElementMonTool.PathInRootFile = "Stats/JEM/BS/input"
#ToolSvc.BS_L1JetElementMonTool.OutputLevel = DEBUG

ToolSvc.Sim_L1JetElementMonTool.DataType = "Sim"  #BS or Sim data?
ToolSvc.Sim_L1JetElementMonTool.JetElementLocation = "Sim_JetElements"
ToolSvc.Sim_L1JetElementMonTool.PathInRootFile = "Stats/JEM/Sim/input"
#ToolSvc.Sim_L1JetElementMonTool.OutputLevel = DEBUG

####################### JEMs ################################
ToolSvc.BS_L1JEMMonTool.DataType = "BS"  #BS or Sim data?
ToolSvc.BS_L1JEMMonTool.JEMHitsLocation = "BS_JEMHits"
ToolSvc.BS_L1JEMMonTool.JEMEtSumsLocation = "BS_JEMEtSums"
ToolSvc.BS_L1JEMMonTool.JEMRoILocation = "BS_JEMRoIs"
ToolSvc.BS_L1JEMMonTool.PathInRootFile = "Stats/JEM/BS"
#ToolSvc.BS_L1JEMMonTool.OutputLevel = DEBUG

ToolSvc.Sim_L1JEMMonTool.DataType = "Sim"  #BS or Sim data?
ToolSvc.Sim_L1JEMMonTool.JEMHitsLocation = "Sim_JEMHits"
ToolSvc.Sim_L1JEMMonTool.JEMEtSumsLocation = "Sim_JEMEtSums"
ToolSvc.Sim_L1JEMMonTool.JEMRoILocation = "Sim_JEMRoIs"
ToolSvc.Sim_L1JEMMonTool.PathInRootFile = "Stats/JEM/Sim"
#ToolSvc.Sim_L1JEMMonTool.OutputLevel = DEBUG

####################### CMMs ################################
ToolSvc.BS_L1CMMMonTool.DataType = "BS"  #BS or Sim data?
ToolSvc.BS_L1CMMMonTool.CMMJetHitsLocation = "BS_CMMJetHits"
ToolSvc.BS_L1CMMMonTool.CMMEtSumsLocation = "BS_CMMEtSums"
ToolSvc.BS_L1CMMMonTool.CMMRoILocation = "BS_CMMRoIs"
ToolSvc.BS_L1CMMMonTool.JEMHitsLocation = "BS_JEMHits"
ToolSvc.BS_L1CMMMonTool.JEMEtSumsLocation = "BS_JEMEtSums"
ToolSvc.BS_L1CMMMonTool.PathInRootFile = "Stats/CMM/BS"
#ToolSvc.BS_L1CMMMonTool.OutputLevel = DEBUG

ToolSvc.Sim_L1CMMMonTool.DataType = "Sim"  #BS or Sim data?
ToolSvc.Sim_L1CMMMonTool.CMMJetHitsLocation = "Sim_CMMJetHits"
ToolSvc.Sim_L1CMMMonTool.CMMEtSumsLocation = "Sim_CMMEtSums"
ToolSvc.Sim_L1CMMMonTool.CMMRoILocation = "Sim_CMMRoIs"
ToolSvc.Sim_L1CMMMonTool.JEMHitsLocation = "Sim_JEMHits"
ToolSvc.Sim_L1CMMMonTool.JEMEtSumsLocation = "Sim_JEMEtSums"
ToolSvc.Sim_L1CMMMonTool.PathInRootFile = "Stats/CMM/Sim"
#ToolSvc.Sim_L1CMMMonTool.OutputLevel = DEBUG

####################### Calorimeter ################################
#ToolSvc.L1CaloBSTool.DataType = "BS"  #BS or Sim data?
#ToolSvc.L1CaloBSTool.BS_TriggerTowerContainer = "BS_TriggerTowers"
#ToolSvc.L1CaloBSTool.BS_JetElementContainer = "BS_JetElements"
#ToolSvc.L1CaloBSTool.PathInRootFile = "Stats/Calo"
#ToolSvc.L1CaloBSTool.OutputLevel = DEBUG

####################### CPMs ################################
ToolSvc.L1BSCPMMonTool.HistogramPrefix = "BS"
ToolSvc.L1BSCPMMonTool.TriggerTowerLocation = "BS_TriggerTowers"
ToolSvc.L1BSCPMMonTool.CPMTowerLocation = "BS_CPMTowers"
ToolSvc.L1BSCPMMonTool.CPMHitsLocation = "BS_CPMHits"
ToolSvc.L1BSCPMMonTool.CMMCPHitsLocation = "BS_CMMCPHits"
ToolSvc.L1BSCPMMonTool.CPMRoILocation = "BS_CPMRoIs"
#ToolSvc.L1BSCPMMonTool.OutputLevel = DEBUG

ToolSvc.L1SimCPMMonTool.HistogramPrefix = "Sim"
ToolSvc.L1SimCPMMonTool.TriggerTowerLocation = "BS_TriggerTowers"
ToolSvc.L1SimCPMMonTool.CPMTowerLocation = "Sim_CPMTowers"
ToolSvc.L1SimCPMMonTool.CPMHitsLocation = "Sim_CPMHits"
ToolSvc.L1SimCPMMonTool.CMMCPHitsLocation = "Sim_CMMCPHits"
ToolSvc.L1SimCPMMonTool.CPMRoILocation = "Sim_CPMRoIs"
#ToolSvc.L1SimCPMMonTool.OutputLevel = DEBUG

####################### Simulation ################################
ToolSvc.SimBSMonTool.BS_JEMEtSumsLocation = "BS_JEMEtSums"
ToolSvc.SimBSMonTool.Sim_JEMEtSumsLocation = "Sim_JEMEtSums"
ToolSvc.SimBSMonTool.BS_JEMHitsLocation = "BS_JEMHits"
ToolSvc.SimBSMonTool.Sim_JEMHitsLocation = "Sim_JEMHits"
ToolSvc.SimBSMonTool.BS_CMMEtSumsLocation = "BS_CMMEtSums"
ToolSvc.SimBSMonTool.Sim_CMMEtSumsLocation = "Sim_CMMEtSums"
ToolSvc.SimBSMonTool.BS_CMMJetHitsLocation = "BS_CMMJetHits"
ToolSvc.SimBSMonTool.Sim_CMMJetHitsLocation = "Sim_CMMJetHits"
ToolSvc.SimBSMonTool.PathInRootFile = "Stats/JEM/SimBS"
#ToolSvc.SimBSMonTool.OutputLevel = DEBUG


## FileKey must match that given to THistSvc
monMan.FileKey = "CaloBSMon"

## Set global monitoring parameters: see the AthenaMonManager class
## in the Control/AthenaMonitoring package
monMan.ManualDataTypeSetup = True
monMan.DataType            = "collisions"
monMan.Environment         = "tier0"
monMan.ManualRunLBSetup    = True
monMan.Run                 = 1
monMan.LumiBlock           = 1

####################################################################################################################
###################################### Ende ########################################################################
####################################################################################################################
