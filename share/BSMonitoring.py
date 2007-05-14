####################################################################################################################
##################################### Input / Output  ##############################################################
####################################################################################################################

# file produced by re-running LVL1 using offline release 12.
BSRDOInput=["daq_Athena_DC2_0005012_file01.data"]     #jetjet
#BSRDOInput=["daq_Athena_DC2_0007003_file01.data"]      #single e
DetDescrVersion = "ATLAS-DC3-02"

MessageSvc.OutputLevel = INFO
#NovaCnvSvc.OutputLevel = INFO
#IOVDbSvc.OutputLevel = INFO
# Number of events to be processed 
theApp.EvtMax = 5
MessageSvc.debugLimit = 1000000

####################################################################################################################
##################################### Konfiguration Simulation #####################################################
####################################################################################################################

DoCaloMon=False
DoBSMon=True 

include ("TriggerMonitoring/Simulation_topOptions.py")

####################################################################################################################
##################################### Monitoring ###################################################################
####################################################################################################################

#unpack BS, let the TrigT1Calo-Simulation run on the BS, do the Monitoring
include ("TriggerMonitoring/TrigT1CaloBSMonitoring_topOptions.py")
