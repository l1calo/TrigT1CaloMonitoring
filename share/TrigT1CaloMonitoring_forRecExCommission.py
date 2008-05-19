Offline= not ATLASCosmicFlags.doOnline
CompareWithSimulation=True

#MaxEnergyRange is set individually

#================================= TriggerMenu ===================================
if not ATLASCosmicFlags.doCTPMon:
    log.info("will setup LVL1ConfigSvc and add instance to ServiceMgr")
    from TrigConfigSvc.TrigConfigSvcConfig import LVL1ConfigSvc
    LVL1ConfigSvc = LVL1ConfigSvc('LVL1ConfigSvc')
    LVL1ConfigSvc.ConfigSource = "XML"
    #LVL1ConfigSvc.XMLFile = "/afs/cern.ch/user/n/neusiedl/myLVL1config.xml"
    LVL1ConfigSvc.XMLFile = "L1MenuM5.xml"
    LVL1ConfigSvc.CreateLegacyObjects = True
    LVL1ConfigSvc.DumpTTVmap = False
    # LVL1ConfigSvc.OutputLevel = VERBOSE
    svcMgr += LVL1ConfigSvc
    theApp.CreateSvc += [ "TrigConf::LVL1ConfigSvc/LVL1ConfigSvc" ]

#================================= L1Calo-Simulation =============================
from TrigT1CaloSim.TrigT1CaloSimConf import LVL1__JetElementMaker
from AthenaCommon.AlgSequence import AlgSequence
myjob = AlgSequence()
myjob += LVL1__JetElementMaker( 'JetElementMaker' )
#myjob.JetElementMaker.TriggerTowerLocation ="TriggerTowers"
myjob.JetElementMaker.JetElementLocation ="Sim_JetElements"
#job.JetElementMaker.OutputLevel = VERBOSE

if CompareWithSimulation:
    include ("TrigT1CaloMonitoring/TrigT1CaloMonitoring_L1CaloSimulation.py")

#================================= Monitoring configuration ======================
topSequence += AthenaMonManager( "L1CaloMonManager" )
L1CaloMan = topSequence.L1CaloMonManager

## get a handle on the ToolSvc
from AthenaCommon.AppMgr import ToolSvc as toolSvc

#=================================================================================
#================================= PPr ===========================================
#=================================================================================
include("CaloConditions/CaloConditions_jobOptions.py")
include("LArDetDescr/LArDetDescr_joboptions.py")

#--------------------------------- PPM -------------------------------------------
from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import PPrMon
L1PPrMonTool = PPrMon(
    name = "L1PPrMonTool",
    BS_TriggerTowerContainer = "TriggerTowers",
    #DistPerChannel = True,
    DistPerChannel = False,
    DistPerChannelAndTimeSlice = False,
    LUTHitMap_Thresh0 = 1,
    LUTHitMap_Thresh1 = 3,
    LUTHitMap_Thresh2 = 7,
    ADCHitMap_Thresh = 50,
    MaxEnergyRange = 256,
    ADCTimingPerChannel = False,
    EMFADCCut = 40,
    HADFADCCut = 40,
    ADCPedestal = 35,
    PathInRootFile = "L1Calo/1_PPr",
    ErrorPathInRootFile = "L1Calo/01_Errors_PPr",
    EventPathInRootFile = "L1Calo",
    Offline = Offline,
    #OutputLevel = VERBOSE,
    #OutputLevel = INFO,
    )
ToolSvc += L1PPrMonTool
L1CaloMan.AthenaMonTools += [ L1PPrMonTool ]

#---------------------------- Performance Checks -----------------------------------

#=================================================================================
#=================================== JEP =========================================
#=================================================================================

#------------------------------------ JEM ----------------------------------------
from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import JEMMon
BS_L1JEMMonTool = JEMMon(
    name = "BS_L1JEMMonTool",
    TypeOfData = "BS",  #BS or Sim data?
    JetElementLocation = "JetElements",
    JEMHitsLocation = "JEMHits",
    JEMEtSumsLocation = "JEMEtSums",
    JEMRoILocation = "JEMRoIs",
    MaxEnergyRange = 1024,
    PathInRootFile = "L1Calo/2_JEP_JEM",
    ErrorPathInRootFile = "L1Calo/02_Errors_JEM",
    Offline = Offline,
    #OutputLevel = VERBOSE,
    )
ToolSvc += BS_L1JEMMonTool
L1CaloMan.AthenaMonTools += [ BS_L1JEMMonTool ]

#if CompareWithSimulation:
#    L1CaloMan.AthenaMonTools += [ "JEMMon/Sim_L1JEMMonTool" ]
#    toolSvc.Sim_L1JEMMonTool.DataType = "Sim"  #BS or Sim data?
#    toolSvc.Sim_L1JEMMonTool.JetElementLocation = "Sim_JetElements"
#    toolSvc.Sim_L1JEMMonTool.JEMHitsLocation = "Sim_JEMHits"
#    toolSvc.Sim_L1JEMMonTool.JEMEtSumsLocation = "Sim_JEMEtSums"
#    toolSvc.Sim_L1JEMMonTool.JEMRoILocation = "Sim_JEMRoIs"
#    toolSvc.Sim_L1JEMMonTool.MaxEnergyRange = MaxEnergyRange
#    toolSvc.Sim_L1JEMMonTool.PathInRootFile = "L1Calo/Sim/2_JEP_JEM"
#toolSvc.Sim_L1JEMMonTool.OutputLevel = DEBUG

#----------------------------------- CMM ------------------------------------------
from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import CMMMon
BS_L1CMMMonTool = CMMMon (
    name = "BS_L1CMMMonTool",
    TypeOfData = "BS",  #BS or Sim data?
    CMMJetHitsLocation = "CMMJetHits",
    CMMEtSumsLocation = "CMMEtSums",
    CMMRoILocation = "CMMRoIs",
    MaxEnergyRange = 32767,
    PathInRootFile = "L1Calo/3_JEP_CMM",
    ErrorPathInRootFile = "L1Calo/03_Errors_CMM",
    Offline = Offline,
    #OutputLevel = VERBOSE,
    )
ToolSvc += BS_L1CMMMonTool
L1CaloMan.AthenaMonTools += [ BS_L1CMMMonTool ]

#if CompareWithSimulation:
#    L1CaloMan.AthenaMonTools += [ "CMMMon/Sim_L1CMMMonTool" ]
#    toolSvc.Sim_L1CMMMonTool.DataType = "Sim"  #BS or Sim data?
#    toolSvc.Sim_L1CMMMonTool.CMMJetHitsLocation = "Sim_CMMJetHits"
#    toolSvc.Sim_L1CMMMonTool.CMMEtSumsLocation = "Sim_CMMEtSums"
#   toolSvc.Sim_L1CMMMonTool.CMMRoILocation = "Sim_CMMRoIs"
#   toolSvc.Sim_L1CMMMonTool.MaxEnergyRange = MaxEnergyRange
#   toolSvc.Sim_L1CMMMonTool.PathInRootFile = "L1Calo/Sim/3_JEP_CMM"
#   toolSvc.Sim_L1CMMMonTool.OutputLevel = DEBUG

#--------------------- Transmission and Performance ------------------------------
from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import JEPTransPerfMon
JEPTransPerfMonTool = JEPTransPerfMon (
    name = "JEPTransPerfMonTool",
    BS_JetElementLocation = "JetElements",
    #BS_TriggerTowerLocation = "TriggerTowers",
    BS_TriggerTowerLocation = "Sim_JetElements",
    NoLUTSlices=1,

    BS_JEMHitsLocation = "JEMHits",
    Sim_JEMHitsLocation = "Sim_JEMHits",
    BS_JEMEtSumsLocation = "JEMEtSums",
    Sim_JEMEtSumsLocation = "Sim_JEMEtSums",
    BS_JEMRoILocation = "JEMRoIs",
    Sim_JEMRoILocation = "Sim_JEMRoIs",

    BS_CMMJetHitsLocation = "CMMJetHits",
    Sim_CMMJetHitsLocation = "Sim_CMMJetHits",
    BS_CMMEtSumsLocation = "CMMEtSums",
    Sim_CMMEtSumsLocation = "Sim_CMMEtSums",
    BS_CMMRoILocation = "CMMRoIs",
    Sim_CMMRoILocation = "Sim_CMMRoIs",

    PathInRootFile = "L1Calo/3_JEP_TransmissionAndPerformance",
    Offline = Offline,
    CompareWithSimulation = CompareWithSimulation,

    #OutputLevel = VERBOSE,
    OutputLevel = WARNING,
    )
ToolSvc += JEPTransPerfMonTool
L1CaloMan.AthenaMonTools += [ JEPTransPerfMonTool ]

#=================================================================================
#===================================== CP ========================================
#=================================================================================
from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import TrigT1CaloCpmMonTool
L1BSCPMMonTool = TrigT1CaloCpmMonTool (
    name = "L1BSCPMMonTool",
    TriggerTowerLocation = "TriggerTowers",
    CPMTowerLocation = "CPMTowers",
    CPMHitsLocation = "CPMHits",
    CMMCPHitsLocation = "CMMCPHits",
    CPMRoILocation = "CPMRoIs",
    RootDirectory = "L1Calo",
    MaxEnergyRange = 50,
    #SingleDirectory = False,
    Offline = Offline,
    #OutputLevel = DEBUG,
    )
ToolSvc += L1BSCPMMonTool
L1CaloMan.AthenaMonTools += [ L1BSCPMMonTool ]

from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import CPMSimBSMon
CPMSimBSMonTool = CPMSimBSMon("CPMSimBSMonTool",
                  CompareWithSimulation = CompareWithSimulation)
ToolSvc += CPMSimBSMonTool
L1CaloMan.AthenaMonTools += [ CPMSimBSMonTool ]
#ToolSvc.CPMSimBSMonTool.OutputLevel = DEBUG

#=================================================================================
# FileKey must match that given to THistSvc
L1CaloMan.FileKey = "GLOBAL"
