#!/bin/bash

min=19151031
#min=19155022
max=19155022

for runId in `cat goodRunNumber_3GeV.txt`
do
    if [   $((runId))  -le  $((max)) ]
    then
        if [   $((runId))  -ge  $((min)) ]
        then
            echo "**********************************  PROCESSING SUBMISSION FOR RUN ${runId}  **********************************"
            rm -rf production/$runId
            mkdir production/$runId

            #star-submit-template -template submitByRun.xml -entities RUN=$runId
            star-submit-template -template submitByRun_perFile.xml -entities RUN=$runId,listOfFiles=/star/u/ypwang/disk01/analysisFlow/polarization/FXT3_85_EPD/filelist_full/${runId}.txt
            echo "**********************************  ENDING SUBMISSION FOR RUN ${runId}  **********************************"
        fi
    fi
done

