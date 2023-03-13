import os
import re
import math
import HistogramTools as HT
from itertools import chain
from functools import partial

from bamboo.plots import Plot, SummedPlot
from bamboo import treefunctions as op
from bamboo import treedecorators as td
from bamboo import scalefactors
from bamboo.analysisutils import makePileupWeight
from bamboo.analysisutils import configureRochesterCorrection, configureJets, configureType1MET
from bamboo.scalefactors import get_correction

import utils


def getYearFromEra(era):
    """ Go from '2017'/'2018UL'/'2016ULpreVFP' to '17'/'18'/'16' """
    return re.search(r"20([0-9]+).*", era).group(1)


#### Configure NanoAOD decoration

def getNanoAODDescription(era, isMC, doRocCor=True):
    groups = ["PV_", "Flag_", "HLT_", "MET_"]
    collections = ["nFatJet","nJet", "nMuon", "nElectron","nTau","nTrigObj"]
    varReaders = []

    if doRocCor:
        varReaders.append(td.nanoRochesterCalc)

    if isMC:
        mcGroups = ["Pileup_", "GenMET_", "Generator_"]
        mcCollections = ["nGenDressedLepton", "nGenJet", "nGenPart","nGenJetAK8"]
        # NOTE: not adding the Jet module to data because the JECs in nanoAODv9 are up-to-date... but this could change!
        # NOTE: not including MET since we're not using it... but this could change
        # NOTE: no corrections to softdrop FatJet=("pt", "mass", "msoftdrop")
        varReaders.append(td.CalcCollectionsGroups(Jet=("pt", "mass"), MET=("pt", "phi"),FatJet=("pt", "mass")))
        return td.NanoAODDescription(groups=groups + mcGroups, collections=collections + mcCollections, systVariations=varReaders)
    else:
        # NOTE: JECs in nanoAODv9 are up-to-date..         
        varReaders.append(td.CalcCollectionsGroups(Jet=("pt", "mass"), MET=("pt", "phi"), FatJet=("pt", "mass")))
        return td.NanoAODDescription(groups=groups, collections=collections, systVariations=varReaders)


JECTagDatabase = {
    "2016ULpreVFP":{
        "MC":"Summer19UL16APV_V7_MC",
        "B_preVFP":"Summer19UL16APV_RunBCD_V7_DATA",
        "C_preVFP":"Summer19UL16APV_RunBCD_V7_DATA",
        "D_preVFP":"Summer19UL16APV_RunBCD_V7_DATA",
        "E_preVFP":"Summer19UL16APV_RunEF_V7_DATA",
        "F_preVFP":"Summer19UL16APV_RunEF_V7_DATA"
     },
    "2016ULpostVFP":{
        "MC":"Summer19UL16_V7_MC",
        "F_postVFP":"Summer19UL16_RunFGH_V7_DATA",
        "G_postVFP":"Summer19UL16_RunFGH_V7_DATA",
        "H_postVFP":"Summer19UL16_RunFGH_V7_DATA"
     },
    "2017UL":{
        "MC":"Summer19UL17_V5_MC",
        "B":"Summer19UL17_RunB_V5_DATA",
        "C":"Summer19UL17_RunC_V5_DATA",
        "D":"Summer19UL17_RunD_V5_DATA",
        "E":"Summer19UL17_RunE_V5_DATA",
        "F":"Summer19UL17_RunF_V5_DATA"
     },
    "2018UL":{
        "MC":"Summer19UL18_V5_MC",
        "A":"Summer19UL18_RunA_V5_DATA",
        "B":"Summer19UL18_RunB_V5_DATA",
        "C":"Summer19UL18_RunC_V5_DATA",
        "D":"Summer19UL18_RunD_V5_DATA"
     },
}
JERTagDatabase = {
    "2016ULpreVFP": "Summer20UL16APV_JRV3_MC",
    "2016ULpostVFP": "Summer20UL16_JRV3_MC",
    "2017UL": "Summer19UL17_JRV2_MC",
    "2018UL": "Summer19UL18_JRV2_MC"
}

def configureJetMETCorrections(tree, era, isNotWorker, isMC, backend, sampleName, splitJER=True, jesScheme="Merged", configMET=False, configFatJet=False):
    if isMC:
        # NOTE: not including variations that would only affect MET since we're not using it... but this could change
        exclJetSysts = []
        if splitJER:
            exclJetSysts += [f"jer{x}" for x in range(2, 6) ]
        if jesScheme == "Merged":
            sources = "Merged"
            exclJetSysts += ["jesTotal"]
            exclJetSysts += [ f"jesHF_{e}" for e in ["2016", "2017", "2018"] ]
            exclJetSysts += [ f"jesEC2_{e}" for e in ["2016", "2017", "2018"] ]
            exclJetSysts.append("jesHF")
        elif jesScheme == "All":
            # here we specify explicitly the sources we use, to avoid having the jetMet 
            # calculator pre-compute all the sources we won't use use in the end
            sources = ["AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias",
                       "Fragmentation", "SinglePionECAL", "SinglePionHCAL",
                       "FlavorQCD",
                       "FlavorPureGluon", "FlavorPureQuark", "FlavorPureCharm", "FlavorPureBottom",
                       "TimePtEta", "RelativeJEREC1", "RelativePtBB", "RelativePtEC1",
                       "RelativeBal", "RelativeSample", "RelativeFSR", "RelativeStatFSR", "RelativeStatEC",
                       "PileUpDataMC", "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1"]
        else:
            sources = ["Total"]
        exclJetSysts += list(chain.from_iterable([ (f"{j}up", f"{j}down") for j in exclJetSysts ]))

        configureJets(tree._Jet, "AK4PFchs", jec=JECTagDatabase[era]["MC"], smear=JERTagDatabase[era],
                      jecLevels=[], # NOTE: not re-applying the JEC, only computing uncertainties!
                      jesUncertaintySources=sources, regroupTag="V2", enableSystematics=lambda v: v not in exclJetSysts,
                      splitJER=splitJER, addHEM2018Issue=(era == "2018UL"), mayWriteCache=isNotWorker,
                      isMC=isMC, backend=backend, uName=sampleName)

        if configMET:
            configureType1MET(tree._MET, jec=JECTagDatabase[era]["MC"], smear=JERTagDatabase[era],
                              jesUncertaintySources="Merged", regroupTag="V2", splitJER=True, addHEM2018Issue=(era == "2018UL"), mayWriteCache=isNotWorker,
                              isMC=isMC, backend=backend, uName=sampleName)

        if configFatJet:
            configureJets(tree._FatJet, "AK8PFPuppi", jec=JECTagDatabase[era], smear=JERTagDatabase[era],
                            jecLevels=[], # NOTE: not re-applying the JEC, only computing uncertainties!
                            genMatchDR=0.4,
                            jesUncertaintySources=sources, regroupTag="V2", enableSystematics=lambda v: v not in exclJetSysts,
                            splitJER=splitJER, addHEM2018Issue=(era == "2018UL"), mayWriteCache=isNotWorker,
                            isMC=isMC, backend=backend, uName=sampleName)

    else:
        sources = None
        runEra = utils.getRunEra(sampleName)
        """"      
        configureJets(tree._Jet, "AK4PFchs", jec=JECTagDatabase[era][runEra],
#                      jecLevels=[], # NOTE: not re-applying the JEC, only computing uncertainties!,
                      # NOTE: if left out the recommendations are used: L1FastJet, L2Relative, L3Absolute, and also L2L3Residual for data
                      mayWriteCache=isNotWorker, 
                      isMC=isMC, backend=backend, uName=sampleName)
        
        if configMET:
            configureType1MET(tree._MET, jec=JECTagDatabase[era][runEra],
                                mayWriteCache=isNotWorker,
                                isMC=isMC, backend=backend, uName=sampleName)
        
        if configFatJet:
            configureJets(tree._FatJet, "AK8PFPuppi", jec=JECTagDatabase[era], smear=JERTagDatabase[era],
                            jecLevels=[], # NOTE: not re-applying the JEC, only computing uncertainties!
                            mayWriteCache=isNotWorker,
                            isMC=isMC, backend=backend, uName=sampleName)
        """        
#### Reco-level object definitions


def flagDef(flags, era, isMC):
    # from https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    cuts = [
        flags.goodVertices,
        flags.globalSuperTightHalo2016Filter,
        flags.HBHENoiseFilter,
        flags.HBHENoiseIsoFilter,
        flags.EcalDeadCellTriggerPrimitiveFilter,
        flags.BadPFMuonDzFilter,
        flags.BadPFMuonFilter,
        flags.eeBadScFilter
    ]
    if era == '2017UL' or era == '2018UL':
        cuts.append(flags.ecalBadCalibFilter)
    return cuts

def addMuonRocCor(be, origMuons, era, sample, isMC):
    if era == "2016ULpreVFP":
        era = "2016aUL"
    if era == "2016ULpostVFP":
        era = "2016bUL"
    rochesterFile = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "Rochester", f"RoccoR{era}.txt")
    configureRochesterCorrection(origMuons, rochesterFile, isMC=isMC, backend=be)

def muonTriggerDef(HLT, sample, era, isMC):
    cuts = []
    if era == '2016ULpreVFP' or era == '2016ULpostVFP':
        cuts.append(op.OR(HLT.IsoMu24, HLT.IsoTkMu24))
    if era == '2017UL':
        cuts.append(HLT.IsoMu27)
    if era == '2018UL':
        cuts.append(HLT.IsoMu24)
    if not isMC:
        cuts.append(op.c_bool("SingleMuon" in sample))
    return cuts

def eleTriggerDef(TrigObj, HLT, ele, sample, era, isMC):
    cuts = []
    if era == '2016ULpreVFP' or era == '2016ULpostVFP':
        cuts.append(HLT.Ele27_WPTight_Gsf)
    if era == '2017UL':
        cuts.append(op.OR(
            op.AND(HLT.Ele32_WPTight_Gsf_L1DoubleEG, op.rng_any(TrigObj, lambda obj: op.AND(op.deltaR(obj.p4, ele.p4) < 0.1, obj.filterBits & 1024))),
            HLT.Ele28_eta2p1_WPTight_Gsf_HT150))
    if era == '2018UL':
        cuts.append(op.OR(HLT.Ele32_WPTight_Gsf, HLT.Ele28_eta2p1_WPTight_Gsf_HT150))
    if not isMC:
        if era == '2016ULpreVFP' or era == '2016ULpostVFP' or era == '2017UL':
            cuts.append(op.c_bool("SingleElectron" in sample))
        if era == '2018UL':
            cuts.append(op.c_bool("EGamma" in sample)) # only called EGamma in 2018
    return cuts

def muonDef(era):
    def muonDefImpl(mu):
        if era == '2016ULpreVFP' or era == '2016ULpostVFP' or era == '2018UL':
            ptCut = 26.
        elif era == '2017UL':
            ptCut = 29.
        return op.AND(
            mu.pt > ptCut,
            op.abs(mu.eta) < 2.4,
            # tight ID
            mu.tightId,
            # tight deltabeta ISO R=0.4
            mu.pfRelIso04_all < 0.15,
        )
    return muonDefImpl

def vetoMuonDef(mu):
    return op.AND(
        mu.pt > 15., op.abs(mu.eta) < 2.4,
        # loose ID
        mu.looseId,
        # loose deltabeta ISO R=0.4
        mu.pfRelIso04_all < 0.25,
    )

def eleDef(era):
    def eleDefImpl(ele):
        absEtaSC = op.abs(ele.eta + ele.deltaEtaSC)
        ele_pt = ele.pt # / ele.eCorr # uncomment to use uncalibrated electron pt
        if era == '2016ULpreVFP' or era == '2016ULpostVFP':
            eraCut = op.AND(ele_pt > 29, op.abs(ele.eta) < 2.5)
        if era == '2017UL' or era == '2018UL':
            eraCut = op.OR(
                op.AND(ele_pt > 34., op.abs(ele.eta) < 2.5),
                op.AND(ele_pt > 30., op.abs(ele.eta) < 2.1),
            )
        return op.AND(
            eraCut,
            op.OR(absEtaSC < 1.4442, absEtaSC > 1.566),
            # "average" d0 and dz cuts, to be tuned?
            op.OR(
                # barrel
                op.AND(absEtaSC <= 1.479, op.abs(ele.dxy) < 0.05, op.abs(ele.dz) < 0.1),
                # endcap
                op.AND(absEtaSC > 1.479, op.abs(ele.dxy) < 0.1, op.abs(ele.dz) < 0.2),
            ),
            # tight cut-based ID
            ele.cutBased == 4,
        )
    return eleDefImpl

def vetoEleDef(ele):
    ele_pt = ele.pt # / ele.eCorr # uncomment to use uncalibrated electron pt
    return op.AND(
        ele_pt > 15., op.abs(ele.eta) < 2.5,
        # veto cut-based ID
        ele.cutBased >= 1,
    )

def jetDef(jet):
    return op.AND(
        jet.pt > 30., op.abs(jet.eta) < 2.4,
        # tight lepton veto jet ID
        jet.jetId & 4,
        # loose puID for jets with pt < 50
        op.OR(jet.pt > 50, jet.puId > 0)
    )

def fatjetDef(fatjet):
    return op.AND(
        fatjet.pt > 200., op.abs(fatjet.eta) < 2.4,
        fatjet.jetId & 2,
    )

# Clean jets from leptons -> since we have veto'd extra loose leptons we
# don't have to use vetoElectrons/Muons at this point
# However we NEED rng_any since we don't know how many electrons/muons we have (could be 0 or 1)
# Also make sure the jets are sorted by Pt (not guaranteed since JER is applied)
def cleanJets(jets, muons, electrons, sort=True):
    jets = op.select(jets, lambda jet: op.AND(
            op.NOT(op.rng_any(electrons, lambda ele: op.deltaR(jet.p4, ele.p4) < 0.4)),
            op.NOT(op.rng_any(muons, lambda mu: op.deltaR(jet.p4, mu.p4) < 0.4))
        ))
    if sort:
        jets = op.sort(jets, lambda j: -j.pt)
    return jets

def cleanFatJets(fatjets, muons, electrons, sort=True):
    fatjets = op.select(fatjets, lambda fatjet: op.AND(
            op.NOT(op.rng_any(electrons, lambda ele: op.deltaR(fatjet.p4, ele.p4) < 0.8)),
            op.NOT(op.rng_any(muons, lambda mu: op.deltaR(fatjet.p4, mu.p4) < 0.8))
        ))
    if sort:
        fatjets = op.sort(fatjets, lambda fj: -fj.pt)
    return fatjets

bTagWorkingPoints = {
    "2016ULpreVFP": {
        "btagDeepFlavB": { # DeepJet
            "L": 0.0508,
            "M": 0.2598,
            "T": 0.6502
        },
        "btagDeepB": { # DeepCSV
            "L": 0.2027,
            "M": 0.6001,
            "T": 0.8819
        },
    },
    "2016ULpostVFP": {
        "btagDeepFlavB": {# DeepJet
            "L": 0.0480,
            "M": 0.2489,
            "T": 0.6377
        },
        "btagDeepB": { # DeepCSV
            "L": 0.1918,
            "M": 0.5847,
            "T": 0.8767
        },
    },
    "2017UL": {
        "btagDeepFlavB": {# DeepJet
            "L": 0.0532,
            "M": 0.3040,
            "T": 0.7476
        },
        "btagDeepB": { # DeepCSV
            "L": 0.1355,
            "M": 0.4506,
            "T": 0.7738
        },
    },
    "2018UL": {
        "btagDeepFlavB": {# DeepJet
            "L": 0.0490,
            "M": 0.2783,
            "T": 0.7100
        },
        "btagDeepB": { # DeepCSV
            "L": 0.1208,
            "M": 0.4168,
            "T": 0.7665
        },
    },
}

def bTagDef(jets, era, wp="M", tagger="btagDeepFlavB"):
    return op.select(jets, lambda jet: getattr(jet, tagger) >= bTagWorkingPoints[era][tagger][wp])

def lightTagDef(jets, era, wp="M", tagger="btagDeepFlavB"):
    return op.select(jets, lambda jet: getattr(jet, tagger) < bTagWorkingPoints[era][tagger][wp])

class corrMET(object):
    def __init__(self, rawMET, pv, sample, era, isMC): #From https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection_withUL17andUL18andUL16.h (preliminary)
        if isMC:
            if era == "2016ULpreVFP":
                xcorr = (0.188743, -0.136539)
                ycorr = (-0.0127927, -0.117747)
            if era == "2016ULpostVFP":
                xcorr = (0.153497, 0.231751)
                ycorr = (-0.00731978, -0.243323)
            if era == "2017UL":
                xcorr = (0.300155, -1.90608)
                ycorr = (-0.300213, 2.02232)
            if era == "2018UL":
                xcorr = (-0.183518, -0.546754)
                ycorr = (-0.192263, 0.42121)
        else:
            runEra = utils.getRunEra(sample)
            if era == "2016ULpreVFP" or era == "2016ULpostVFP":
                if runEra == "B_preVFP":
                    xcorr = (0.0214894, 0.188255)
                    ycorr = (-0.0876624, -0.812885)
                if runEra == "C_preVFP":
                    xcorr = (0.032209, -0.067288)
                    ycorr = (-0.113917, -0.743906)
                if runEra == "D_preVFP":
                    xcorr = (0.0293663, -0.21106)
                    ycorr = (-0.11331, -0.815787)
                if runEra == "E_preVFP":
                    xcorr = (0.0132046, -0.20073)
                    ycorr = (-0.134809, -0.679068)
                if runEra == "F_preVFP":
                    xcorr = (0.0543566, -0.816597)
                    ycorr = (-0.114225, -1.17266)
                if runEra == "F_postVFP":
                    xcorr = (-0.134616, 0.89965)
                    ycorr = (-0.134616, 1.0385)
                if runEra == "G_postVFP":
                    xcorr = (-0.121809, 0.584893)
                    ycorr = (-0.0558974, -0.891234)
                if runEra == "H_postVFP":
                    xcorr = (-0.0868828, 0.703489)
                    ycorr = (-0.0888774, -0.902632)
            if era == "2017UL":
                if runEra == "B":
                    xcorr = (0.211161, -0.419333)
                    ycorr = (-0.251789, 1.28089)
                if runEra == "C":
                    xcorr = (0.185184, 0.164009)
                    ycorr = (-0.200941, 0.56853)
                if runEra == "D":
                    xcorr = (0.201606, -0.426502)
                    ycorr = (-0.188208, 0.58313)
                if runEra == "E":
                    xcorr = (0.162472, -0.176329)
                    ycorr = (-0.138076, 0.250239)
                if runEra == "F":
                    xcorr = (0.210639, -0.72934)
                    ycorr = (-0.198626, -1.028)
            if era == "2018UL":
                if runEra == "A":
                    xcorr = (-0.263733, 1.91115)
                    ycorr = (-0.0431304, 0.112043)
                if runEra == "B":
                    xcorr = (-0.400466, 3.05914)
                    ycorr = (-0.146125, 0.533233)
                if runEra == "C":
                    xcorr = (-0.430911, 1.42865)
                    ycorr = (-0.0620083, 1.46021)
                if runEra == "D":
                    xcorr = (-0.457327, 1.56856)
                    ycorr = (-0.0684071, 0.928372)
        METxcorr = xcorr[0] * pv.npvs + xcorr[1]
        METycorr = ycorr[0] * pv.npvs + ycorr[1]
        corrMETx = rawMET.pt * op.cos(rawMET.phi) + METxcorr
        corrMETy = rawMET.pt * op.sin(rawMET.phi) + METycorr
        self.pt = op.sqrt(corrMETx**2 + corrMETy**2)
        atan = op.atan(corrMETy / corrMETx)
        self.phi = op.multiSwitch(
                (corrMETx > 0, atan),
                (corrMETy > 0, atan + math.pi),
                atan - math.pi
            )

# scale factors

def pogEraFormat(era):
    return era.replace("UL", "") + "_UL"

def localizePOGSF(era, POG, fileName):
    return os.path.join("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration", "POG", POG, pogEraFormat(era), fileName)

def makePUWeight(tree, era, selection):
    goldenJSON = f"Collisions{getYearFromEra(era)}_UltraLegacy_goldenJSON"
    puTuple = (localizePOGSF(era, "LUM", "puWeights.json.gz"), goldenJSON)
    return makePileupWeight(puTuple, tree.Pileup.nTrueInt, systName="pileup", sel=selection)

# maps name of systematic to name of correction inside of jsons
leptonSFLib = {
    "electron_ID": "UL-Electron-ID-SF",
    "electron_reco": "UL-Electron-ID-SF",
    "electron_trigger": "EleTriggerSF",
    "muon_reco": "NUM_TrackerMuons_DEN_genTracks",
    "muon_ID": "NUM_TightID_DEN_TrackerMuons",
    "muon_iso": "NUM_TightRelIso_DEN_TightIDandIPCut",
    "muon_trigger": {
        "2016ULpreVFP": "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
        "2016ULpostVFP": "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
        "2017UL": "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight",
        "2018UL": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
    },
}

def getLeptonSF(era, systName):
    if systName == "muon_trigger":
        corrName = leptonSFLib[systName][era]
    else:
        corrName = leptonSFLib[systName]

    if "muon" in systName:
        path = localizePOGSF(era, "MUO", "muon_Z.json.gz")
    elif systName == "electron_trigger":
        path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..", "..", "scale-factors", "eleTrigSFs", era + "_EleTriggerSF_NanoAODv9_v0.json") 
    elif "electron" in systName:
        path = localizePOGSF(era, "EGM", "electron.json.gz")

    return path, corrName

def getScaleFactor(era, noSel, systName, defineOnFirstUse=True):
    fileName, correction = getLeptonSF(era, systName)

    if "muon" in systName:
        etaParam = "abseta"
        etaExpr = lambda mu: op.abs(mu.eta)
    elif "electron" in systName:
        etaParam = "eta"
        etaExpr = lambda el: el.eta + el.deltaEtaSC
    else:
        raise ValueError("Only muon or electron SFs are handled here!")

    if "muon" in systName:
        return get_correction(fileName, correction, params={"pt": lambda mu: mu.pt, etaParam: etaExpr, "year": pogEraFormat(era)},
                              systParam="ValType", systNomName="sf",
                              systVariations={f"{systName}up": "systup", f"{systName}down": "systdown"},
                              defineOnFirstUse=defineOnFirstUse, sel=noSel)
    elif systName == "electron_trigger":
        return get_correction(fileName, correction, params={"pt": lambda el: el.pt, etaParam: etaExpr},
                              systParam="sf", systNomName="central",
                              systVariations=("up", "down"), systName=systName,
                              defineOnFirstUse=defineOnFirstUse, sel=noSel)
    else:
        wp = "Tight" if systName == "electron_ID" else "RecoAbove20"
        return get_correction(fileName, correction, params={"pt": lambda el: el.pt, etaParam: etaExpr, "year": era.replace("UL", ""), "WorkingPoint": wp},
                              systParam="ValType", systNomName="sf",
                              systVariations={f"{systName}up": "sfup", f"{systName}down": "sfdown"},
                              defineOnFirstUse=defineOnFirstUse, sel=noSel)

def getL1PrefiringSystematic(tree):
    return op.systematic(tree.L1PreFiringWeight_Nom, name="L1prefire", up=tree.L1PreFiringWeight_Up, down=tree.L1PreFiringWeight_Dn)

# Return exclusive muon selection (without and with applied trigger), with trigger, ID, and iso scale factors (for MC)
def buildMuonSelections(tree, noSel, muons, vetoMuons, electrons, vetoElectrons, sample, era, isMC):
    scaleFactors = []
    if isMC:
        muonRecoSF = getScaleFactor(era, noSel, systName="muon_reco")
        muonIDSF = getScaleFactor(era, noSel, systName="muon_ID")
        muonIsoSF = getScaleFactor(era, noSel, systName="muon_iso")
        scaleFactors = [ muonRecoSF(muons[0]), muonIDSF(muons[0]), muonIsoSF(muons[0]) ]

    oneMuSel = noSel.refine("muon",
                    cut=op.AND(
                        op.rng_len(muons) == 1,
                        op.rng_len(vetoMuons) == 1,
                        op.rng_len(vetoElectrons) == 0
                    ),
                    weight=scaleFactors
                )
    triggerSFWeights = []
    if isMC:
        muonTriggerSF = getScaleFactor(era, oneMuSel, systName="muon_trigger")
        triggerSFWeights.append(muonTriggerSF(muons[0]))
        triggerSFWeights.append(getL1PrefiringSystematic(tree))
    oneMuTriggerSel = oneMuSel.refine("muonTrigger",
                                    cut=muonTriggerDef(tree.HLT, sample, era, isMC),
                                    weight=triggerSFWeights)

    return oneMuSel, oneMuTriggerSel

# Return exclusive electron selection (without and with applied trigger), with trigger, ID/iso, and reco scale factors (for MC)
def buildElectronSelections(tree, noSel, muons, vetoMuons, electrons, vetoElectrons, sample, era, isMC):
    scaleFactors = []
    if isMC:
        eleRecoSF = getScaleFactor(era, noSel, systName="electron_reco")
        eleIDSF = getScaleFactor(era, noSel, systName="electron_ID")
        scaleFactors = [ eleRecoSF(electrons[0]), eleIDSF(electrons[0]) ]

    oneEleSel = noSel.refine("electron",
                    cut=op.AND(
                        op.rng_len(vetoMuons) == 0,
                        op.rng_len(vetoElectrons) == 1,
                        op.rng_len(electrons) == 1
                    ),
                    weight=scaleFactors
                )
    triggerSFWeights = []
    if isMC:
        eleTriggerSF = getScaleFactor(era, oneEleSel, systName="electron_trigger")
        triggerSFWeights.append(eleTriggerSF(electrons[0]))
        triggerSFWeights.append(getL1PrefiringSystematic(tree))
        if era == "2017UL":
            # HLT Z_vtx correction
            triggerSFWeights.append(0.991)
    oneEleTriggerSel = oneEleSel.refine("electronTrigger",
                                    cut=eleTriggerDef(tree.TrigObj, tree.HLT, electrons[0], sample, era, isMC),
                                    weight=triggerSFWeights)

    return oneEleSel, oneEleTriggerSel

def get_bTagSF_fixWP(wp, flav, era, sel, use_nominal_jet_pt=False, heavy_method="comb",
                    syst_prefix="btagSF_deepjet_fixWP_", decorr_wps=False, decorr_eras=True, full_scheme=False, full_scheme_mapping=None, defineOnFirstUse=True,
                    prepostVFPfix=False):
    params = {
        "pt": lambda j: op.forSystematicVariation(j.pt, "nominal") if use_nominal_jet_pt else j.pt,
        "abseta": lambda j: op.abs(j.eta), "working_point": wp, "flavor": flav
    }
    systName = syst_prefix + ("light" if flav == 0 else "heavy")
    wpDecorrTag = f"_{wp}" if decorr_wps else ""
    systVariations = {}
    for d in ("up", "down"):
        if not decorr_eras and not full_scheme:
            systVariations[f"{systName}{d}"] = d
        elif decorr_eras and (not full_scheme or flav == 0):
            systVariations[f"{systName}{d}"] = f"{d}_correlated"
            systVariations[f"{systName}{wpDecorrTag}_{era}{d}"] = f"{d}_uncorrelated"
        elif full_scheme and flav > 0:
            systVariations[f"{syst_prefix}statistic{wpDecorrTag}_{era}{d}"] = f"{d}_statistic"
            for var,varBTV in full_scheme_mapping.items():
                if varBTV is None:
                    systVariations[f"{syst_prefix}{var}{d}"] = f"{d}_{var}"
                else:
                    systVariations[f"{var}{d}"] = f"{d}_{varBTV}"

    method = "incl" if flav == 0 else heavy_method

    # FIXME temporary for 2016: use preVFP mistag SFs in postVFP era
    if prepostVFPfix and era == "2016ULpostVFP" and flav == 0:
        era = "2016ULpreVFP"
    return get_correction(localizePOGSF(era, "BTV", "btagging.json.gz"), f"deepJet_{method}", params=params,
                          systParam="systematic", systNomName="central",
                          systVariations=systVariations, sel=sel, defineOnFirstUse=defineOnFirstUse)

def makeBtagSFWPs(era, sel, wps, use_nominal_jet_pt=False, defineOnFirstUse=True, **kwargs):
    """wps from looser to tighter"""
    get_bTagSF = partial(get_bTagSF_fixWP, era=era, sel=sel, use_nominal_jet_pt=use_nominal_jet_pt, defineOnFirstUse=defineOnFirstUse, **kwargs)

    # Functions selecting the right SF depending on jet flavour
    # The flavour can't be passed as parameter to correctionlib because the
    # uncertainties depend on it.
    bTagSF = {}
    for wp in wps:
        bTagSF[wp] = lambda j, _wp=wp: op.multiSwitch(
            (j.hadronFlavour == 5, get_bTagSF(_wp, 5)(j)),
            (j.hadronFlavour == 4, get_bTagSF(_wp, 4)(j)),
            get_bTagSF(_wp, 0)(j))

    # functions defined in bamboo
    wFail = op.extMethod("scalefactorWeightForFailingObject", returnType="double")

    bTagEff_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scale-factors", "btagEff", f"btagEff_deepJet_{era.replace('UL', '')}.json")
    def get_bTagEff(wp):
        params = {
            "pt": lambda j: op.forSystematicVariation(j.pt, "nominal") if use_nominal_jet_pt else j.pt, "eta": lambda j: j.eta,
            "jetFlavor": lambda j: j.hadronFlavour, "workingPoint": wp
        }
        return get_correction(bTagEff_file, "btagEff", params=params, sel=sel,
                             defineOnFirstUse=defineOnFirstUse)
    bTagEff = { wp: get_bTagEff(wp) for wp in wps }

    def jet_SF(j):
        factors = []
        # [L,M,T] -> (if discr >= T, then: )
        tightest = wps[-1]
        factors.append( (j.btagDeepFlavB >= bTagWorkingPoints[era]["btagDeepFlavB"][tightest], bTagSF[tightest](j)) )
        # [L,M,T] -> (elif discr >= M, then: ), (elif discr >= L, then: )
        for i in range(len(wps)-1, 0, -1):
            tighter = wps[i]
            looser = wps[i-1]
            factors.append( (j.btagDeepFlavB >= bTagWorkingPoints[era]["btagDeepFlavB"][looser], wFail(bTagSF[tighter](j), bTagEff[tighter](j), bTagSF[looser](j), bTagEff[looser](j))) )
        # [L,M,T] -> (else: )
        loosest = wps[0]
        factors.append( wFail(bTagSF[loosest](j), bTagEff[loosest](j)) )
        return op.multiSwitch(*factors)

    # method 1a: product over jets, factors depend on discriminator value
    return jet_SF

def get_bTagSF_itFit(flav, era, sel, use_nominal_jet_pt=False, syst_prefix="btagSF_deepjet_shape_", decorr_eras=True, jet_mapping=None):
    systListUnCorr = []
    if flav == 4:
        systListCorr = ["cferr1", "cferr2"]
    else:
        systListCorr = ["hf", "lf"]
        statSyst = ["hfstats1", "hfstats2", "lfstats1", "lfstats2"]
        if decorr_eras:
            systListUnCorr = statSyst
        else:
            systListCorr += statSyst

    systVariations = {}
    for d in ("up", "down"):
        for var in systListCorr:
            systVariations[f"{syst_prefix}{var}{d}"] = f"{d}_{var}"
        for var in systListUnCorr:
            systVariations[f"{syst_prefix}{var}_{era}{d}"] = f"{d}_{var}"
        if jet_mapping is not None:
            for var,varBTV in jet_mapping.items():
                if varBTV is None:
                    systVariations[f"{syst_prefix}{var}{d}"] = f"{d}_{var}"
                else:
                    systVariations[f"{var}{d}"] = f"{d}_{varBTV}"

    params = {
        "pt": lambda j: op.forSystematicVariation(j.pt, "nominal") if use_nominal_jet_pt else j.pt, "abseta": lambda j: op.abs(j.eta),
        "discriminant": lambda j: j.btagDeepFlavB, "flavor": flav
    }
    return get_correction(localizePOGSF(era, "BTV", "btagging.json.gz"), "deepJet_shape", params=params,
                          systParam="systematic", systNomName="central",
                          systVariations=systVariations, sel=sel)

def makeBtagSFItFit(era, sel, **kwargs):
    # The flavour can't be passed as parameter to correctionlib because the
    # uncertainties depend on it.
    get_bTagSF = partial(get_bTagSF_itFit, era=era, sel=sel, **kwargs)
    bTagSF = lambda j: op.multiSwitch(
        (j.hadronFlavour == 5, get_bTagSF(5)(j)),
        (j.hadronFlavour == 4, get_bTagSF(4)(j)),
        get_bTagSF(0)(j))
    return bTagSF
