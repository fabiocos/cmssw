#!/bin/sh -e

# Pass in a message and the command's status
die() {
    if (( $2 == 0 )); then
            echo "$1: expected a nonzero status, got 0" >&2
            #exit 1
    fi

    echo "$1: status $2" >&2
    #exit 0
}

TEST_DIR=$CMSSW_BASE/src/DataFormats/TrackerRecHit2D/test
DATA_DIR=$CMSSW_BASE/src/DataFormats/TrackerRecHit2D/data/

(cmsRun ${TEST_DIR}/readMTDTrackingRecHits_cfg.py filePrepend=${DATA_DIR}) ||
    die 'Failure using readMTDTrackingRecHits_cfg.py' "$?"
