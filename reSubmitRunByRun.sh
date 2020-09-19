#!/bin/bash

starver SL19e

rm filelist_reSubmit/*

min=19151031
#min=19155022
max=19155022

dir=/star/u/ypwang/disk01/analysisFlow/polarization/FXT3_85_EPD/polarization_ld

for runId in `cat ../goodRunNumber_3GeV.txt`
do
    if [   $((runId))  -le  $((max)) ]
    then
        if [   $((runId))  -ge  $((min)) ]
        then
            for a in $(cat /star/u/ypwang/disk01/analysisFlow/polarization/FXT3_85_EPD/filelist_full/${runId}.txt)
            #for a in $(cat /star/u/ypwang/disk01/analysisFlow/polarization/FXT3_85_EPD/filelist_full_LocalNfs/${runId}.txt)
            do
                picoTree=`basename $a`
                lambdaTree=`echo ${runId}_$picoTree | sed 's/.picoDst.root/_polarization.root/'`

                if [ -e $dir/production/${runId}/${lambdaTree} ]; then
                    fsize=$(ls -l $dir/production/${runId}/${lambdaTree} | awk '{print $5}')
                    if [   $((fsize))  -le  270000 ]
                    then
                        echo $a >> $dir/filelist_reSubmit/${runId}.list
                        rm $dir/production/${runId}/${lambdaTree}
                    fi
                else
                    echo $a >> $dir/filelist_reSubmit/${runId}.list
                fi
            done
            
            if [ -f $dir/filelist_reSubmit/${runId}.list ]; then
                nJobs=`sed -n '$=' $dir/filelist_reSubmit/${runId}.list`
                echo "RESUBMIT $nJobs JOBS FOR ${runId}!"
                star-submit-template -template submitByRun_perFile.xml -entities RUN=$runId,listOfFiles=$dir/filelist_reSubmit/$runId.list
            fi
        fi
    fi
done

condor_q
