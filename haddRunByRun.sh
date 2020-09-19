#!/bin/bash

min=$1
#max=19151084
max=$2

for runId in `cat goodRunNumber_3GeV.txt`
do
    if [   $((runId))  -le  $((max)) ]
    then
        if [   $((runId))  -ge  $((min)) ]
        then
            mkdir -p production_perRun/${runId}
            hadd production_perRun/${runId}/${runId}.root production/${runId}/*.root
            fsize=$(ls -l production_perRun/${runId}/${runId}.root | awk '{print $5}')
            if [   $((fsize))  -le  1000 ]
            then
                echo " file size of ${runId}.root: $fsize "
                rm production_perRun/${runId}/${runId}.root
            fi
        fi
    fi
done
