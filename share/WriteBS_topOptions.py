####################################################################################################################
##################################### Bytestream Optionen ##########################################################
####################################################################################################################


include( "TriggerRelease/TrigWriteBS_topOptions.py" )

StreamBS = Algorithm( "StreamBS" )
StreamBS.ItemList    = [ ]
 
include ( "TrigT1CaloByteStream/WritePpmBS_jobOptions.py" )
include ( "TrigT1CaloByteStream/WriteJepBS_jobOptions.py" )
include ( "TrigT1CaloByteStream/WriteCpBS_jobOptions.py" )

include ( "TrigT1CaloByteStream/WriteJepRoiBS_jobOptions.py" )
include ( "TrigT1CaloByteStream/WriteCpmRoiBS_jobOptions.py" )

# Algs
theApp.Dlls += [ "TrigT1CaloByteStream" ]

