###############################################################
#
# Set Flags for Running the Trigger in readBS Mode
# includes Trigger_topOptions.py to do the real work.
# Reads from file RawEvent.re (file must already exits
# run athena with TrigWriteBS_Flags.py to produce it)
# 
# use:
# athena.py -bs TrigReadBS_Flags.py
# athena.py -bs -c onlyCalo=True -bs TrigReadBS_Flags.py
# athena.py -bs -c onlyID=True TrigReadBS_Flags.py
# athena.py -bs -c onlyMuon=True TrigReadBS_Flags.py
#==============================================================

# Set detector description
DetDescrVersion = "Rome-Initial"

#BSRDOInput=["RawEvent.re"]
BSRDOInput=["daq_Athena_DC2_0005144_file01.data"]

# Set flags that determine the input data source
# Set only one flag True:
readZEBRA=False
readPool=False 
readBS=True     
# These flags can be set in addition:
writeBS=False     
transientBS=False

include ("TriggerRelease/TriggerFlags.py")
TriggerFlags.fakeLVL1=False
TriggerFlags.doLVL1=False
TriggerFlags.doNtuple=True
TriggerFlags.useOfflineSpacePoints=False
useROOTNtuple=False # False : hbook ntuple

# Set Truth flags
doTruth=False    # needs to be true if creating Fake RoI 
 
#-------------------------------------------------------------
# End of setting flags
#-------------------------------------------------------------

#include( "TriggerRelease/Trigger_topOptions.py" )
include( "TriggerRelease/TrigReadBS_topOptions.py" )


# Specify the Converters
ByteStreamCnvSvc = Service( "ByteStreamCnvSvc" )
ByteStreamCnvSvc.InitCnvs += [ "DataVector<LVL1::TriggerTower>" ]
ByteStreamCnvSvc.InitCnvs += [ "DataVector<LVL1::JetElement>" ]
ByteStreamCnvSvc.InitCnvs += [ "DataVector<LVL1::JEMHits>" ]
ByteStreamCnvSvc.InitCnvs += [ "DataVector<LVL1::JEMEtSums>" ]
# DLLs
theApp.Dlls += [ "TrigT1CaloByteStream" ]
# Algs



#theApp.TopAlg = [ "PpmTester" ]
#Tester = Algorithm( "PpmTester" )

#
# Set Algorithm
#

theApp.Dlls += [ "AthenaMonitoring", "TrigT1CaloMonitoring" ]
theApp.Dlls += [ "TrigT1CaloMonitoring" ]
theApp.TopAlg += [ "AthenaMon/TrigT1CaloMon"]
TrigT1MonTool = Algorithm( "TrigT1CaloMon" )
TrigT1MonTool.AthenaMonTools = ["TrigT1CaloBSMonTool/TrigT1CaloMon"]

ByteStreamAddressProviderSvc = Service( "ByteStreamAddressProviderSvc" )
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::TriggerTower>/LVL1TriggerTowers" ]
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::JetElement>/LVL1JetElements" ]
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::JEMHits>/JEMHits" ]
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::JEMEtSums>/JEMEtSums" ]
    

ByteStreamAddressProviderSvc = Service( "ByteStreamAddressProviderSvc" )
ByteStreamAddressProviderSvc.TypeNames += [ "DataVector<LVL1::TriggerTower>/LVL1TriggerTowers" ]

#
#  Output stuff
#

# THistService for native root in Athena
theApp.Dlls += [ "GaudiSvc" ]
THistSvc = Algorithm( "THistSvc" )
## Note that the THistSvc doesn't do a "recreate", so ...
RootHistOutputFileName="MonitorTrigT1BS.root"
if os.path.exists(RootHistOutputFileName):
   os.remove(RootHistOutputFileName)
THistSvc.Output = ["MonitorOutPut DATAFILE='"+RootHistOutputFileName+"' OPT='NEW'"]
    
MessageSvc.OutputLevel = DEBUG
MessageSvc.Format = "% F%18W%S%7W%R%T %0W%M"
#MessageSvc.OutputLevel = INFO
#NovaCnvSvc.OutputLevel = INFO

# Number of events to be processed 
theApp.EvtMax = 10


