####################################################################################
#
#             Lvl 1 Calo Trigger Monitoring Job Options File
#
####################################################################################

#-----------------------------------
# Force a given Detector description
#-----------------------------------
DetDescrVersion="ATLAS-DC3-02"

#-----------------------------------
# Pool Converters
#-----------------------------------
include( "ParticleEventAthenaPool/AOD_PoolCnv_jobOptions.py") # may not need this
include( "TrkEventAthenaPool/TrkEventAthenaPool_joboptions.py" ) # or this
include( "TrigT1EventAthenaPool/TrigT1EventAthenaPool_joboptions.py" )
include( "CaloIdCnv/CaloIdCnv_joboptions.py" )
include( "LArAthenaPool/LArRecAthenaPool_ReadTest_jobOptions.py")

#-----------------------------------
#  Set Input File
#-----------------------------------
EventSelector.InputCollections = [ 
"PythiaZee.ESD.root"
]

#-----------------------------------
#  DLLs
#-----------------------------------

theApp.Dlls   += [ "TruthParticleAlgs" ]
theApp.Dlls   += [ "AnalysisTools" ]
theApp.Dlls   += [ "UserAnalysis" ]
theApp.Dlls   += [ "TileRecAlgs" ]
theApp.Dlls   += [ "LArAthenaPoolPoolCnv" ]
theApp.Dlls   += [ "TileEventAthenaPoolPoolCnv" ]
theApp.Dlls   += [ "LArClusterRec" ]

#
# Set Algorithm
#

theApp.Dlls += [ "AthenaMonitoring", "TrigT1CaloMonitoring" ]
theApp.Dlls += [ "TrigT1CaloMonitoring" ]
theApp.TopAlg += [ "AthenaMon/TrigT1CaloMon"]
TrigT1MonTool = Algorithm( "TrigT1CaloMon" )
TrigT1MonTool.AthenaMonTools = ["TrigT1CaloMonTool/TrigT1CaloMon"]

#
#  Output stuff
#

# THistService for native root in Athena
theApp.Dlls += [ "GaudiSvc" ]
THistSvc = Algorithm( "THistSvc" )
## Note that the THistSvc doesn't do a "recreate", so ...
RootHistOutputFileName="MonitorTrigT1.root"
if os.path.exists(RootHistOutputFileName):
   os.remove(RootHistOutputFileName)
THistSvc.Output = ["MonitorOutPut DATAFILE='"+RootHistOutputFileName+"' OPT='NEW'"]

# Number of Events to process
theApp.EvtMax = 200


# Set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )
MessageSvc = Service( "MessageSvc" )
MessageSvc.OutputLevel = INFO
