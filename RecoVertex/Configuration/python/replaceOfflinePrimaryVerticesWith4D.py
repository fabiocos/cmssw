import FWCore.ParameterSet.Config as cms
#
# makes OfflinePrimaryVertices equivalent to OfflinePrimaryVerticesWithBS
# by changing the input vertices collection of the sorted PV
# see file https://github.com/cms-sw/cmssw/blob/master/RecoVertex/Configuration/python/RecoVertex_cff.py
#
def swapOfflinePrimaryVerticesToUse4D(process):
    if hasattr(process,'offlinePrimaryVertices'):
        process.offlinePrimaryVertices = process.offlinePrimaryVertices4D.clone()
    if hasattr(process,'offlinePrimaryVerticesWithBS'):
        process.offlinePrimaryVerticesWithBS = process.offlinePrimaryVertices4DWithBS.clone()

    return process
