import os
import re
import json
import yaml
import tarfile, tempfile
from copy import deepcopy

import logging
logger = logging.getLogger("utils")

from bamboo.plots import Plot, SummedPlot, CutFlowReport
from bamboo import treefunctions as op
from bamboo import scalefactors
from bamboo.root import ROOT
from bamboo.workflow import SampleTask
from bamboo import treeproxies as tp

import HistogramTools as HT


# some hacks to make sure we get "raw strings" dumped in the YAML file
def str_is_latex(data):
    return ("$" in data or "\\" in data or "{" in data)
def yaml_latex_representer(dumper, data):
    return dumper.represent_scalar(u"tag:yaml.org,2002:str", data, style="'" if str_is_latex(data) else None)
yaml.add_representer(str, yaml_latex_representer)


def getOpts(uname, **kwargs):
    label = None
    opts = {}
    if "1ele" in uname:
        label = "1 electron"
    elif "2ele" in uname:
        label = "2 electrons"
    elif "1mu" in uname:
        label = "1 muon"
    elif "2mu" in uname:
        label = "2 muons"
    elif "1lep" in uname:
        label = "1 lepton (e/#mu)"
    elif "2lep" in uname:
        label = "2 leptons (ee/#mu#mu)"
    if "gen" in uname:
        label = "Gen: " + label
    if "0j" in uname:
        label += ", #geq 0 jets"
    elif "1j" in uname:
        label += ", #geq 1 jets"
    elif "2j" in uname:
        label += ", #geq 2 jets"
    elif "3j" in uname:
        label += ", #geq 3 jets"
    elif "4j" in uname:
        label += ", #geq 4 jets"
    elif "5j" in uname:
        label += ", #geq 5 jets"
    elif "6j" in uname:
        label += ", #geq 6 jets"
    elif "7j" in uname:
        label += ", #geq 7 jets"
    if "0fj" in uname:
        label += ", #geq 0 fatjets"
    elif "1fj" in uname:
        label += ", #geq 1 fatjets"
    elif "2fj" in uname:
        label += ", #geq 2 fatjets"
    if "1b" in uname:
        label += ", #geq 1 b"
    elif "2b" in uname:
        label += ", #geq 2 b"
    elif "3b" in uname:
        label += ", #geq 3 b"
    elif "4b" in uname:
        label += ", #geq 4 b"
    if "3l" in uname:
        label += ", #geq 3 l"
    if label:
        opts = {
            "labels": [{"text": label, "position": [0.205, 0.912]}]
        }
    opts.update(kwargs)
    return opts

def getCounter(i):
    if i <= 0:
        return str(i)
    if i == 1:
        return "1st"
    if i == 2:
        return "2nd"
    if i == 3:
        return "3rd"
    if i >= 4:
        return "{}th".format(i)

def getRunEra(sample):
    """Return run era (A/B/...) for data sample"""
    result = re.search(r'Run201.([A-Z]?)', sample)
    if result is None:
        raise RuntimeError("Could not find run era from sample {}".format(sample))
    if "preVFP" in sample:
        suffix = "_preVFP"
    elif "postVFP" in sample:
        suffix = "_postVFP"
    else:
        suffix = ""
    runEra = result.group(1) + suffix
    return runEra

def getAnyDataSample(samples, mcChecker, era=None):
    """Return any data from the sample list.
    If an era is specified, find a dataset from that era.
    mcChecker a function: mcChecker(sample) is True if sample is MC - typically pass analysisModule.isMC"""
    for smpNm,smp in samples.items():
        if not mcChecker(smpNm):
            if (era is None) or (era and smp["era"] == era):
                return smpNm
def getAnyMCSample(samples, mcChecker, era=None, noSystSample=True, signalSample=True):
    """Return any MC from the sample list.
    If an era is specified, find a sample from that era.
    mcChecker a function: mcChecker(sample) is True if sample is MC - typically pass analysisModule.isMC
    If noSystSample is True, systematic variation samples are not considered
    If signalSample is True, the returned sample must be a signal"""
    for smpNm,smp in samples.items():
        if mcChecker(smpNm):
            if (era is None) or (era and smp["era"] == era):
                if noSystSample and "syst" in smp:
                    continue
                if signalSample and not smp.get("is_signal", False):
                    continue
                if smp.get("gen_only", False):
                    continue
                return smpNm

def makeMergedPlots(categDef, newCat, name, binning, var=None, **kwargs):
    """ Make a series of plots which will be merged.
    - categDef can either be e.g.:
        - `[("mu", muSelection), ("el", elSelection), ...]`, in which case the same variable `var` is used for all sub-categories
        - `[("mu", muSelection, muVar), ("el", elSelection, elVar), ...]`, for cases where the variable is different for each sub-category
    - `newCat`: name of the merged category
    - `name`: name of the merged plot (-> full plot name is newCat_name)
    - `var`: variable to plot (if it is the same for all categories)
    - `binning`: binning to be used
    Any further named args will be forwarded to the plot constructor.
    The variables can also be iterables for multi-dimensional plots.
    """
    plotsToAdd = []
    plotopts = kwargs.pop("plotopts", {})

    for cat in categDef:
        if len(cat) == 2:
            (catName, catSel), catVar = cat, var
        elif len(cat) == 3:
            catName, catSet, catVar = cat
        else:
            raise Exception(f"{cat} should have 2 or 3 entries")
        thisName = f"{catName}_{name}"
        if not hasattr(catVar, "__len__") or isinstance(catVar, tp.VectorProxy):
            plotType = Plot.make1D
        elif len(catVar) == 2:
            plotType = Plot.make2D
        elif len(catVar) == 3:
            plotType = Plot.make3D
        newPlotopts = getOpts(catName, **plotopts)
        plotsToAdd.append(plotType(thisName, catVar, catSel, binning, plotopts=newPlotopts, **kwargs))

    newPlotopts = getOpts(newCat, **plotopts)

    return plotsToAdd + [SummedPlot(f"{newCat}_{name}", plotsToAdd, plotopts=newPlotopts, **kwargs)]


#### Common tasks (systematics, sample splittings...)
def getTopPtWeight(tree):
    def top_pt_weight(pt):
        return 0.103 * op.exp(-0.0118 * pt) - 0.000134 * pt + 0.973

    lastCopy = op.select(tree.GenPart, lambda p: (p.statusFlags >> 13) & 1)
    tops = op.select(lastCopy, lambda p: p.pdgId == 6)
    antitops = op.select(lastCopy, lambda p: p.pdgId == -6)
    weight = op.switch(op.AND(op.rng_len(tops) >= 1, op.rng_len(antitops) >= 1),
                       op.sqrt(top_pt_weight(tops[0].pt) * top_pt_weight(antitops[0].pt)),
                       1.)

    return weight


def applyTopPtReweighting(tree, sel):
    logger.info("Applying Top Pt reweighting everywhere, with a 'noTopPt' variation without it")

    return sel.refine("topPt", weight=op.systematic(getTopPtWeight(tree), noTopPt=op.c_float(1.)))


import numpy as np
@ROOT.Numba.Declare(['RVec<float>', 'float', 'int'], 'float')
def computeHessianPDFUncertainty(weights, SF=1., max_index=0):
    if len(weights) < 2:
        return 0.
    weights = np.asarray(weights)
    if (max_index != 0) and (max_index != -1):
        weights = weights[:max_index+1]
    return np.sqrt(np.sum((SF * weights[1:] - weights[0])**2))


NO_MESCALE_WGTS = ["ttbb_4fs_openloops_sherpa"]
NO_PS_WGTS = ["tt_herwig7", "ttbb_4fs_openloops_sherpa"]


def addTheorySystematics(plotter, sample, sampleCfg, tree, noSel, qcdVar=True, PSISR=True, PSFSR=True, pdf_mode="simple", topPt=True):
    signal_tag = sampleCfg.get("signal_tag")

    SF = 1.  # FxFx has a strange bug where ME and PDF (Hessian) weights are divided by 2
    if signal_tag == "ttjets_amcatnloFXFX":
        SF = 2.
    SF = op.c_float(SF)

    ##### QCD muR/muF ######
    if qcdVar:
        if hasattr(tree, "LHEScaleWeight") and signal_tag not in NO_MESCALE_WGTS:
            if plotter.qcdScaleVarMode == "separate":
                # for muF and muR separately
                logger.info("Adding separate muF/muR systematics")
                noSel = noSel.refine("qcdMuF", weight=op.systematic(op.c_float(1.), name="qcdMuF", up=SF * tree.LHEScaleWeight[5], down=SF * tree.LHEScaleWeight[3]))
                noSel = noSel.refine("qcdMuR", weight=op.systematic(op.c_float(1.), name="qcdMuR", up=SF * tree.LHEScaleWeight[7], down=SF * tree.LHEScaleWeight[1]))
                noSel = noSel.refine("qcdMuRF", weight=op.systematic(op.c_float(1.), name="qcdMuRF", up=SF * tree.LHEScaleWeight[8], down=SF * tree.LHEScaleWeight[0]))

            elif plotter.qcdScaleVarMode == "combined":
                # for taking envelope of 7-point variation
                qcdScaleVariations = plotter.getQcdScaleVariations(sample, sampleCfg)
                if qcdScaleVariations:
                    logger.info("Adding 7-point muF/muR systematics")
                    qcdScaleVariations = { varNm: SF * tree.LHEScaleWeight[varIdx] for (varIdx, varNm) in qcdScaleVariations }
                    qcdScaleSyst = op.systematic(op.c_float(1.), name="qcdScale", **qcdScaleVariations)
                    noSel = noSel.refine("qcdScale", weight=qcdScaleSyst)
        else:
            logger.warning("LHEScaleWeight not present in tree, muF/muR systematics will not be added")

    isTopNano = "topNanoAOD" in sampleCfg.get("db", "")

    ###### Parton Shower #####
    if hasattr(tree, "PSWeight") and signal_tag not in NO_PS_WGTS:
        if PSISR:
            logger.info("Adding PS ISR systematics")
            if isTopNano:
                psISRSyst = op.systematic(op.c_float(1.), name="psISR", up=tree.PSWeight[27], down=tree.PSWeight[26])
            else:
                psISRSyst = op.systematic(op.c_float(1.), name="psISR", up=tree.PSWeight[0], down=tree.PSWeight[2])
            noSel = noSel.refine("psISR", weight=psISRSyst)

        if PSFSR:
            logger.info("Adding PS FSR systematics")
            if isTopNano:
                psFSRSyst = op.systematic(op.c_float(1.), name="psFSR", up=tree.PSWeight[5], down=tree.PSWeight[4])
            else:
                psFSRSyst = op.systematic(op.c_float(1.), name="psFSR", up=tree.PSWeight[1], down=tree.PSWeight[3])
            noSel = noSel.refine("psFSR", weight=psFSRSyst)
    else:
        logger.warning("PSWeight not present in tree, PS ISR and PS FSR systematics will not be added")

    ###### PDF ######
    pdf_mc = sampleCfg.get("pdf_mc", False)
    nPdfVars = (0, 101) if pdf_mc else (1, 103)
    if hasattr(tree, "LHEPdfWeight"):
        if pdf_mode == "full" and sampleCfg.get("pdf_full", False):
            logger.info("Adding full PDF systematics")
            pdfVars = { f"pdf{i}": SF * tree.LHEPdfWeight[i] for i in range(*nPdfVars) }
        else:
            logger.info("Adding simplified PDF systematics")
            if pdf_mc:
                pdfSigma = op.rng_stddev(tree.LHEPdfWeight)
            else:
                sigmaCalc = op.extMethod("Numba::computeHessianPDFUncertainty", returnType="float")
                max_index = op.c_int(100)
                pdfSigma = sigmaCalc(tree.LHEPdfWeight, SF, max_index)
            pdfVars = { "pdfup": tree.LHEPdfWeight[0] + pdfSigma, "pdfdown": tree.LHEPdfWeight[0] - pdfSigma }
            if not pdf_mc: #also include the alphaS variations separately
                alphaSVars = { "pdfAlphaSup": SF * tree.LHEPdfWeight[101], "pdfAlphaSdown": SF * tree.LHEPdfWeight[102] } # last two weights should be alpha_S variations
                pdfVars.update( alphaSVars )
        if pdf_mc:
            logger.info("This sample has MC PDF variations")
        else:
            logger.info("This sample has Hessian PDF variations")
        noSel = noSel.refine("PDF", weight=op.systematic(op.c_float(1.), **pdfVars))
    else:
        logger.warning("LHEPdfWeight not present in tree, PDF systematics will not be added")


    ###### Top Pt ######
    if topPt and "subprocess" in sampleCfg:
        logger.info("Adding top pt reweighting as additional variation")
        noSel = noSel.refine("topPt", weight=op.systematic(op.c_float(1.), topPt=getTopPtWeight(tree)))

    return noSel


def ttJetFlavCuts(subProc, tree):
    if subProc == "ttbb":
        return (tree.genTtbarId % 100) >= 53
    if subProc == "ttb":
        return op.in_range(50, tree.genTtbarId % 100, 53)
    if subProc == "ttcc":
        return op.in_range(40, tree.genTtbarId % 100, 46)
    if subProc == "ttjj":
        return (tree.genTtbarId % 100) < 41
    if subProc == "ttB":
        return (tree.genTtbarId % 100) >= 51

def splitTTjetFlavours(cfg, tree, noSel):
    subProc = cfg["subprocess"]
    return noSel.refine(subProc, cut=ttJetFlavCuts(subProc, tree))

def removeSignalOverlapPlotIt(mainSignalTag, samplesCfg):
    """Keep a single signal prediction before running plotIt"""
    if not any(smpCfg.get("signal_tag", None) == mainSignalTag for smpCfg in samplesCfg.values()):
        logger.warning(f"Not removing signal overlap in plotIt because main requested signal {mainSignalTag} is not in sample list!")
        return
    toRemove = []
    for smpName,smpCfg in samplesCfg.items():
        if "signal_tag" in smpCfg and smpCfg["signal_tag"] != mainSignalTag:
            toRemove.append(smpName)
    if toRemove:
        logger.info(f"Will remove the following from plotIt list because of signal overlap: {toRemove}")
    for nm in toRemove:
        samplesCfg.pop(nm)

def getSumw(resultsFile, smpCfg, readCounters=None):
    if "generated-events" in smpCfg:
        if isinstance(smpCfg["generated-events"], str):
            genEvts = readCounters(resultsFile)[smpCfg["generated-events"]]
        else:
            genEvts = smpCfg["generated-events"]
    else:
        genEvts = None
    return genEvts

def normalizeAndSumSamples(eras, samples, inDir, outPath, readCounters=lambda f: -1.):
    """
    Produce file containing the sum of all the histograms over the processes, 
    after normalizing the processes by their cross section, sum of weights and luminosity.
    Note: The systematics are handled but are expected to be SAME for all processes and eras.
    A separate output file is produced for each era (`outPath_era.root`), as well as a total one (`outPath_run2.root`).
    """
    for era in eras:
        lumi = eras[era]["luminosity"]
        mergedHists = {}
        for proc,cfg in samples.items():
            if cfg["era"] != era: continue
            if "syst" in cfg: continue
            xs = cfg["cross-section"]
            tf = HT.openFileAndGet(os.path.join(inDir, proc + ".root"))
            sumWgt = getSumw(tf, cfg, readCounters)
            keyList = tf.GetListOfKeys()
            for key in keyList:
                hist = key.ReadObj()
                if not hist.InheritsFrom("TH1"): continue
                hist.Scale(lumi * xs / sumWgt)
                name = hist.GetName()
                if name not in mergedHists:
                    mergedHists[name] = hist.Clone()
                    mergedHists[name].SetDirectory(0)
                else:
                    mergedHists[name].Add(hist)
            tf.Close()
        mergedFile = HT.openFileAndGet(outPath + "_" + era + ".root", "recreate")
        for hist in mergedHists.values():
            hist.Write()
        mergedFile.Close()
    os.system("hadd -f " + outPath + "_run2.root " + " ".join([ f"{outPath}_{era}.root" for era in eras ]))

def produceMEScaleEnvelopes(plots, scaleVariations, path):
    if not scaleVariations:
        return

    logger.info(f"Producing QCD scale envelopes for file {path}")

    tf = HT.openFileAndGet(path, "update")
    listOfKeys = [ k.GetName() for k in tf.GetListOfKeys() ]

    for plot in plots:
        if isinstance(plot, CutFlowReport):
            continue
        # Compute envelope histograms for QCD scale variations
        nominal = tf.Get(plot.name)
        variations = []
        for var in scaleVariations:
            varName = "{}__{}".format(plot.name, var)
            if varName in listOfKeys:
                variations.append(tf.Get(varName))
        if len(variations) != len(scaleVariations):
            logger.warning("Did not find {} variations for plot {} in file {}".format(len(scaleVariations), plot.name, path))
            continue
        if not variations:
            continue
        up,down = HT.getEnvelopeHistograms(nominal, variations)
        up.Write(f"{plot.name}__qcdScaleup", ROOT.TObject.kOverwrite)
        down.Write(f"{plot.name}__qcdScaledown", ROOT.TObject.kOverwrite)

    tf.Close()

def producePDFEnvelopes(plots, task, resultsdir):
    import numpy as np

    sample = task.name
    smpCfg = task.config
    path = os.path.join(resultsdir, task.outputFile)

    logger.info(f"Producing PDF uncertainty envelopes for sample {sample}")

    def sigmaFromReplicasMC(replicas):
        return np.std(replicas, axis=0, ddof=1)

    def sigmaFromReplicasHessian(residuals):
        sq = residuals**2
        return np.sqrt(sq.sum(axis=0))

    def buildTHFromNP(nom_th, values, name):
        new_th = nom_th.Clone(name)
        assert(nom_th.GetNcells() == len(values))
        for i in range(len(values)):
            new_th.SetBinContent(i, values[i])
        return new_th

    tf = HT.openFileAndGet(path, "update")
    listOfKeys = [ k.GetName() for k in tf.GetListOfKeys() ]
    nVar = 0

    for plot in plots:
        if isinstance(plot, CutFlowReport):
            continue
        # Compute envelope histograms for PDF variations
        nominal = tf.Get(plot.name)
        def isPDFVar(name):
            return name.startswith(f"{plot.name}__pdf") and not (
                name.endswith("up") or name.endswith("down"))
        variations = [ tf.Get(varName) for varName in listOfKeys if isPDFVar(varName) ]

        if not variations:
            logger.warning(f"Did not find PDF variations for plot {plot.name} in file {path}")
            continue
        nVar = len(variations)

        replica_values = np.vstack([ np.array(h) for h in variations ])
        nom_values = np.array(nominal)
        if smpCfg.get("pdf_mc", False):
            # PDF MC set
            sigma = sigmaFromReplicasMC(replica_values)
        else:
            # Hessian MC set
            sigma = sigmaFromReplicasHessian(replica_values - nom_values)
        up = buildTHFromNP(nominal, nom_values + sigma, f"{plot.name}__pdfup")
        down = buildTHFromNP(nominal, np.clip(nom_values - sigma, 0., None), f"{plot.name}__pdfdown")
        up.Write("", ROOT.TObject.kOverwrite)
        down.Write("", ROOT.TObject.kOverwrite)

    logger.info(f"Found {nVar} PDF variations for sample {sample}")

    tf.Close()

def postProcSystSamples(taskList, samples, resultsdir, cutFlowPrefix="yields_", readCounters=lambda f: -1.):
    """
    Rename the systematics samples so they can be interpreted by plotIt
    Also fix the normalization of their histograms
    Finally, remove them from the config sample list for plotIt
    """
    for task in taskList:
        if not "syst" in task.config: continue

        nominalSample = task.config["syst"][1]
        systVar = task.config["syst"][0]

        # always remove from list propagated to plotIt since otherwise plotIt counts it as another process
        samples.pop(task.name)

        if nominalSample not in samples:
            logger.warning(f"For syst {systVar} could not find nominal sample {nominalSample}, skipping renormalization step")
            continue

        tfSyst = HT.openFileAndGet(os.path.join(resultsdir, task.outputFile))
        tfNom = HT.openFileAndGet(os.path.join(resultsdir, nominalSample + ".root"), "update")
        wgtRatio = getSumw(tfNom, samples[nominalSample], readCounters) / getSumw(tfSyst, task.config, readCounters)

        for k in tfSyst.GetListOfKeys():
            # skip cutflowreports for now, as systematics are not entirely supported yet
            if k.GetName().startswith(cutFlowPrefix): continue
            hist = k.ReadObj()
            if not hist.InheritsFrom("TH1"): continue
            name = hist.GetName() + "__" + systVar
            hist.SetName(name)
            hist.Scale(wgtRatio)
            hist.Write("", ROOT.TObject.kOverwrite)
        tfNom.cd()
        tfNom.Close()
        tfSyst.Close()

def fillSampleTemplate(template, selEras=None, split_processes=True):
    outTemplate = {}

    # Expand eras
    for name,sample in template.items():
        if "dbs" in sample:
            for era,das in sample["dbs"].items():
                era = str(era)
                if selEras is not None and era not in selEras:
                    continue
                thisSample = deepcopy(sample)
                if "syst" in thisSample:
                    syst,nom = thisSample["syst"]
                    newName = f"{nom}__{era}__{syst}"
                    thisSample["syst"][1] = f"{nom}__{era}"
                else:
                    newName = f"{name}__{era}"
                thisSample["db"] = das
                thisSample["era"] = era
                thisSample.pop("dbs")
                outTemplate[newName] = thisSample
        else:
            # Keep same format of __era also for data
            era = sample["era"]
            newName = f"{name}__{era}"
            outTemplate[newName] = sample

    if not split_processes:
        return outTemplate

    # Expand split processes
    for name, sample in list(outTemplate.items()):
        if "subprocesses" in sample:
            for subProc in sample["subprocesses"]:
                newProc = deepcopy(sample)
                newProc.pop("subprocesses")
                newProc["subprocess"] = subProc
                if subProc in newProc.pop("signal_subprocesses", []):
                    newProc["is_signal"] = True
                if not newProc.get("is_signal", False):
                    newProc.pop("signal_tag", None)
                newProc["group"] = subProc
                newProc["stat-group"] = subProc
                # be careful with the systematics samples (they have two '__' separators)
                newName = name.split("__")[0] + "_" + subProc + "__" + "__".join(name.split("__")[1:])
                if "syst" in newProc:
                    newProc["syst"][1] = newProc["syst"][1].split("__")[0] + "_" + subProc + "__" + newProc["syst"][1].split("__")[1]
                outTemplate[newName] = newProc
            del outTemplate[name]

    return outTemplate
