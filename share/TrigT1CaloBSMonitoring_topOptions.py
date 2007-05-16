
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

include ( "TrigT1CaloByteStream/ReadJepBS_jobOptions.py" )
ByteStreamAddressProviderSvc = Service( "ByteStreamAddressProviderSvc" )
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::JEMRoI>/BS_JEMRoIs" ]
ByteStreamAddressProviderSvc.TypeNames += [ "LVL1::CMMRoI/BS_CMMRoIs" ]

#include ( "TrigT1CaloByteStream/ReadCpBS_jobOptions.py" )

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
theApp.TopAlg += ["LVL1::JetElementMaker/JetElementMaker"]       
theApp.TopAlg += ["LVL1::EnergyTrigger/EnergyTrigger" ]          
theApp.TopAlg += ["LVL1::JetTrigger/JetTrigger" ]         
theApp.TopAlg += ["LVL1::ROD/ROD" ]          
theApp.TopAlg += ["LVL1::JEPCMMMaker/JEPCMMMaker" ]

#*************************************************************
# Set input and output locations  
#*************************************************************
# JetElementMaker
JetElementMaker = Algorithm( "JetElementMaker" )
JetElementMaker.TriggerTowerLocation ="BS_TriggerTowers"  
JetElementMaker.JetElementLocation ="Sim_JetElements"

# JetTrigger
JetTrigger = Algorithm( "JetTrigger" )
JetTrigger.JetElementLocation="Sim_JetElements"
#JetTrigger.JetROIOutputLocation = "Sim_JetROIs"
JetTrigger.JEMHitsLocation = "Sim_JEMHits"

# EnergyTrigger
EnergyTrigger = Algorithm( "EnergyTrigger" )
EnergyTrigger.JetElementLocation="Sim_JetElements"
EnergyTrigger.JEMEtSumsLocation="Sim_JEMEtSums"

# JEP CMMs
JEPCMMMaker = Algorithm( "JEPCMMMaker" )
JEPCMMMaker.JetElementLocation="Sim_JetElements"
JEPCMMMaker.JEMHitsLocation = "Sim_JEMHits"
JEPCMMMaker.JEMEtSumsLocation = "Sim_JEMEtSums"
JEPCMMMaker.CMMJetHitsLocation = "Sim_CMMJetHits"
JEPCMMMaker.CMMEtSumsLocation = "Sim_CMMEtSums"
#JEPCMMMaker.EnergyRoILocation = "Sim_EnergyROIs"
#JEPCMMMaker.JetEtRoILocation = "Sim_JetEtROIs"
#JEPCMMMaker.JetRoILocation = "Sim_JetROIs"
JEPCMMMaker.JEMRoILocation = "Sim_JEMRoIs"
JEPCMMMaker.CMMRoILocation = "Sim_CMMRoIs"
JEPCMMMaker.OutputLevel = VERBOSE


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
#monMan.AthenaMonTools += [ "JetElementMon/L1JetElementMonTool" ]
#monMan.AthenaMonTools += [ "JEMMon/L1JEMMonTool" ]
#monMan.AthenaMonTools += [ "CMMMon/L1CMMMonTool" ]
#monMan.AthenaMonTools += [ "TrigT1CaloBSMonTool/L1CaloBSTool" ]
monMan.AthenaMonTools += [ "JEMRoIMon/L1BSJEMRoIMonTool" ]
monMan.AthenaMonTools += [ "JEMRoIMon/L1SimJEMRoIMonTool" ]
#monMan.AthenaMonTools += [ "CMMRoIMon/L1BSCMMRoIMonTool" ]
monMan.AthenaMonTools += [ "CMMRoIMon/L1SimCMMRoIMonTool" ]

## get a handle on the ToolSvc
#from AthenaCommon.AppMgr import ToolSvc as toolSvc
ToolSvc = Algorithm( "ToolSvc" )

#ToolSvc.L1JetElementMonTool.BS_JetElementLocation = "BS_JetElements"
#ToolSvc.L1JetElementMonTool.Sim_JetElementLocation = "Sim_JetElements"

#ToolSvc.L1JEMMonTool.BS_JEMHitsLocation = "BS_JEMHits"
#ToolSvc.L1JEMMonTool.BS_JEMEtSumsLocation = "BS_JEMEtSums"
#ToolSvc.L1JEMMonTool.Sim_JEMHitsLocation = "Sim_JEMHits"
#ToolSvc.L1JEMMonTool.Sim_JEMEtSumsLocation = "Sim_JEMEtSums"
#ToolSvc.L1JEMMonTool.OutputLevel = DEBUG

#ToolSvc.L1CMMMonTool.BS_CMMJetHitsLocation = "BS_CMMJetHits"
#ToolSvc.L1CMMMonTool.BS_CMMEtSumsLocation = "BS_CMMEtSums"
#ToolSvc.L1CMMMonTool.Sim_CMMJetHitsLocation = "Sim_CMMJetHits"
#ToolSvc.L1CMMMonTool.Sim_CMMEtSumsLocation = "Sim_CMMEtSums"
#ToolSvc.L1CMMMonTool.OutputLevel = DEBUG

#ToolSvc.L1CaloBSTool.BS_TriggerTowerContainer = "BS_TriggerTowers"
#ToolSvc.L1CaloBSTool.BS_JetElementContainer = "BS_JetElements"
#ToolSvc.L1CaloBSTool.OutputLevel = DEBUG

#JEM RoI
ToolSvc.L1BSJEMRoIMonTool.DataType = "BS"  #BS or Sim data?
ToolSvc.L1BSJEMRoIMonTool.JEMRoILocation = "BS_JEMRoIs"
ToolSvc.L1BSJEMRoIMonTool.PathInRootFile = "Stats/BS_JEMRoI"
#ToolSvc.L1BSJEMRoIMonTool.OutputLevel = DEBUG

ToolSvc.L1SimJEMRoIMonTool.DataType = "Sim"  #BS or Sim data?
ToolSvc.L1SimJEMRoIMonTool.JEMRoILocation = "Sim_JEMRoIs"
ToolSvc.L1SimJEMRoIMonTool.PathInRootFile = "Stats/Sim_JEMRoI"
#ToolSvc.L1SimJEMRoIMonTool.OutputLevel = DEBUG

#CMM RoI
#ToolSvc.L1BSCMMRoIMonTool.DataType = "BS"  #BS or Sim data?
#ToolSvc.L1BSCMMRoIMonTool.CMMRoILocation = "BS_CMMRoIs"
#ToolSvc.L1BSCMMRoIMonTool.PathInRootFile = "Stats/BS_CMMRoI"
#ToolSvc.L1BSCMMRoIMonTool.OutputLevel = DEBUG

ToolSvc.L1SimCMMRoIMonTool.DataType = "Sim"  #BS or Sim data?
ToolSvc.L1SimCMMRoIMonTool.CMMRoILocation = "Sim_CMMRoIs"
ToolSvc.L1SimCMMRoIMonTool.PathInRootFile = "Stats/Sim_CMMRoI"
ToolSvc.L1SimCMMRoIMonTool.OutputLevel = DEBUG

#ToolSvc.TrigT1JetMonTool.OutputLevel = DEBUG

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
