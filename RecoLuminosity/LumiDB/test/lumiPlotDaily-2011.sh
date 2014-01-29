#!/bin/sh
currentdir=`pwd`
sarch="slc5_ia32_gcc434"
workdir="/afs/cern.ch/user/l/lumipro/scratch0/mydev/CMSSW_3_11_0"
authdir="/afs/cern.ch/user/l/lumipro"
overviewdir=$currentdir
operationdir=$currentdir
physicsdir=$currentdir
publicresultdir=$currentdir
logpath=$currentdir
withcorrection="--with-correction"

logname="lumiPlotDaily-2011.log"
logfilename="$logpath/$logname"
dbConnectionString="oracle://cms_orcoff_prod/cms_lumi_prod"
#physicsselectionFile="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/StreamExpress/goodrunlist_json.txt"
beamenergy="3500"
beamfluctuation="0.2"
beamstatus="stable"

source /afs/cern.ch/cms/cmsset_default.sh;
cd $workdir
export SCRAM_ARCH="$sarch";
eval `scramv1 runtime -sh`
touch $logfilename
date >> $logfilename
lumiPlotFiller.py -c $dbConnectionString -P $authdir -o $overviewdir -beamstatus $beamstatus -beamenergy $beamenergy -beamfluctuation $beamfluctuation --withTextOutput $withcorrection total2011vstime >> $logfilename 
sleep 1
lumiPlotFiller.py -c $dbConnectionString -P $authdir -o $publicresultdir -beamstatus $beamstatus -beamenergy $beamenergy -beamfluctuation $beamfluctuation --withTextOutput  $withcorrection total2011vstime >> $logfilename 
sleep 1
lumiPlotFiller.py -c $dbConnectionString -P $authdir -o $overviewdir -beamstatus $beamstatus -beamenergy $beamenergy -beamfluctuation $beamfluctuation --withTextOutput  $withcorrection perday2011 >> $logfilename 
sleep 1;
lumiPlotFiller.py -c $dbConnectionString -P $authdir -o $publicresultdir -beamstatus $beamstatus -beamenergy $beamenergy -beamfluctuation $beamfluctuation --withTextOutput  $withcorrection perday2011 >> $logfilename 
sleep 1;
lumiPlotFiller.py -c $dbConnectionString -P $authdir -o $operationdir --withTextOutput  $withcorrection instpeak2011vstime >> $logfilename 
sleep 1
lumiPlotFiller.py -c $dbConnectionString -P $authdir -o $publicresultdir --withTextOutput  $withcorrection instpeak2011vstime >> $logfilename 
sleep 1
lumiPlotFiller.py -c $dbConnectionString -P $authdir -o $operationdir -beamstatus $beamstatus  -beamenergy $beamenergy -beamfluctuation $beamfluctuation --withTextOutput  $withcorrection total2011vsfill >> $logfilename 
sleep 1
lumiPlotFiller.py -c $dbConnectionString -P $authdir -o $operationdir -beamstatus $beamstatus  -beamenergy $beamenergy -beamfluctuation $beamfluctuation --withTextOutput  $withcorrection total2011vsrun >> $logfilename 
#sleep 1
#lumiPlotFiller.py -c $dbConnectionString -P $authdir -o $physicsdir -i $physicsselectionFile --withTextOutput physicsvstime >> $logfilename
#sleep 1
#lumiPlotFiller.py -c $dbConnectionString -P $authdir -o $physicsdir -i $physicsselectionFile --withTextOutput physicsperday >> $logfilename
cd $currentdir
