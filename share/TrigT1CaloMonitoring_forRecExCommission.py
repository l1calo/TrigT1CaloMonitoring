Offline= not athenaCommonFlags.isOnline
CompareWithSimulation=False
if (globalflags.DataSource() == "data") :
    CompareWithSimulation=True
    
#MaxEnergyRange is set individually

#================================= TriggerMenu ===================================
from AthenaCommon.AppMgr import ServiceMgr as svcMgr
#if hasattr(svcMgr,"LVL1ConfigSvc"):
#    log.info("ServiceMgr already includes LVL1ConfigSvc")
#else:
#    log.warning( "Did not find a configured LVL1 configuration service!" )
#    log.info("will setup LVL1ConfigSvc and add instance to ServiceMgr")
#    from TrigConfigSvc.TrigConfigSvcConfig import LVL1ConfigSvc
#    LVL1ConfigSvc = LVL1ConfigSvc('LVL1ConfigSvc')
#    svcMgr += LVL1ConfigSvc
#    svcMgr.LVL1ConfigSvc.ConfigSource = "XML"
#    from TriggerJobOpts.TriggerFlags import TriggerFlags as tf  
#    #tf.inputLVL1configFile = "LVL1config_SingleBeam_v1_7-bit_trigger_types_20080905.xml"
#    tf.inputLVL1configFile = "LVL1config_SingleBeam_v1_7-bit_trigger_types.xml"
#    svcMgr.LVL1ConfigSvc.XMLFile = tf.inputLVL1configFile()
#    svcMgr.LVL1ConfigSvc.CreateLegacyObjects = True
#    svcMgr.LVL1ConfigSvc.DumpTTVmap = False
#    # LVL1ConfigSvc.OutputLevel = VERBOSE
#    theApp.CreateSvc += [ "TrigConf::LVL1ConfigSvc/LVL1ConfigSvc" ]

if globalflags.DataSource() == "data":
  #================================= L1Calo-Simulation =============================
  from TrigT1CaloSim.TrigT1CaloSimConf import LVL1__JetElementMaker
  from AthenaCommon.AlgSequence import AlgSequence
  myjob = AlgSequence()
  myjob += LVL1__JetElementMaker( 'JetElementMaker_Mon' )
  myjob.JetElementMaker_Mon.JetElementLocation ="Sim_JetElements"


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

if globalflags.DataSource() == "data":
  #--------------------- Transmission and Performance ------------------------------
  from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import JEPTransPerfMon
  JEPTransPerfMonTool = JEPTransPerfMon (
      name = "JEPTransPerfMonTool",

      JEPHitsTool = "LVL1::L1JEPHitsTools/L1JEPHitsTools_Mon",
      JetTool = "LVL1::L1JetTools/L1JetTools_Mon",
      JEPEtSumsTool = "LVL1::L1JEPEtSumsTools/L1JEPEtSumsTools_Mon",

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
      OutputLevel = INFO,
      )
  ToolSvc += JEPTransPerfMonTool
  L1CaloMan.AthenaMonTools += [ JEPTransPerfMonTool ]

  from TrigT1CaloTools.TrigT1CaloToolsConf import LVL1__L1JEPHitsTools
  L1JEPHitsTools = LVL1__L1JEPHitsTools("L1JEPHitsTools_Mon")
  L1JEPHitsTools.LVL1ConfigSvc="TrigConf::TrigConfigSvc/TrigConfigSvc"
  ToolSvc += L1JEPHitsTools
  from TrigT1CaloTools.TrigT1CaloToolsConf import LVL1__L1JetTools
  L1JetTools = LVL1__L1JetTools("L1JetTools_Mon")
  L1JetTools.LVL1ConfigSvc="TrigConf::TrigConfigSvc/TrigConfigSvc"
  ToolSvc += L1JetTools
  from TrigT1CaloTools.TrigT1CaloToolsConf import LVL1__L1EtTools
  L1EtTools = LVL1__L1EtTools("L1EtTools_Mon")
  L1EtTools.LVL1ConfigSvc="TrigConf::TrigConfigSvc/TrigConfigSvc"
  ToolSvc += L1EtTools
  from TrigT1CaloTools.TrigT1CaloToolsConf import LVL1__L1JEPEtSumsTools
  L1JEPEtSumsTools = LVL1__L1JEPEtSumsTools("L1JEPEtSumsTools_Mon",
                                            EtTool = "LVL1::L1EtTools/L1EtTools_Mon")
  L1JEPEtSumsTools.LVL1ConfigSvc="TrigConf::TrigConfigSvc/TrigConfigSvc"
  ToolSvc += L1JEPEtSumsTools

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
    MaxEnergyRange = 256,
    #OutputLevel = DEBUG,
    )
ToolSvc += L1BSCPMMonTool
L1CaloMan.AthenaMonTools += [ L1BSCPMMonTool ]

from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import CPMSimBSMon
CPMSimBSMonTool = CPMSimBSMon("CPMSimBSMonTool",
                              EmTauTool = "LVL1::L1EmTauTools/L1EmTauTools_Mon",
                              CompareWithSimulation = CompareWithSimulation)
ToolSvc += CPMSimBSMonTool
L1CaloMan.AthenaMonTools += [ CPMSimBSMonTool ]
#ToolSvc.CPMSimBSMonTool.IgnoreTowersEM  = [ 1890,              #LUT readout
#              4082, 4083, 4146, 4147, 4210, 4211, 4274, 4275 ] #LVDS channel
#ToolSvc.CPMSimBSMonTool.IgnoreTowersHad = [ 3473, 3643, 4824 ] #LUT readout
#ToolSvc.CPMSimBSMonTool.OutputLevel = DEBUG

from TrigT1CaloTools.TrigT1CaloToolsConf import LVL1__L1EmTauTools
L1EmTauTools = LVL1__L1EmTauTools("L1EmTauTools_Mon")
L1EmTauTools.LVL1ConfigSvc="TrigConf::TrigConfigSvc/TrigConfigSvc"
ToolSvc += L1EmTauTools

#=================================================================================
#===================================== ROD =======================================
#=================================================================================
from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import TrigT1CaloRodMonTool
L1BSRODMonTool = TrigT1CaloRodMonTool (
    name = "L1BSRODMonTool",
    #OutputLevel = DEBUG,
    )
ToolSvc += L1BSRODMonTool
L1CaloMan.AthenaMonTools += [ L1BSRODMonTool ]

#=================================================================================
# FileKey must match that given to THistSvc
if not 'DQMonFlags' in dir():
   print "TrigT1CaloMonitoring_forRecExCommission.py: DQMonFlags not yet imported - I import them now"
   from AthenaMonitoring.DQMonFlags import DQMonFlags
L1CaloMan.FileKey             = DQMonFlags.monManFileKey()
L1CaloMan.Environment         = DQMonFlags.monManEnvironment()
L1CaloMan.ManualDataTypeSetup = DQMonFlags.monManManualDataTypeSetup()
L1CaloMan.DataType            = DQMonFlags.monManDataType()
