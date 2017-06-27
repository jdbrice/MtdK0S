#!/bin/bash

echo "job_$1.root"
echo "filelist=$2"

./pico.app config/PicoDstSkimmer.xml --PDS.output.TFile:url=/global/project/projectdirs/star/pwg/starlfs/jdb/MtdK0S/job_$1.root --PDS.input.dst:url=$2


