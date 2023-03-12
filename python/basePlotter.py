import os
import json
import yaml
import glob
import shutil
from itertools import chain

import logging
logger = logging.getLogger("base plotter")

from bamboo.analysismodules import NanoAODModule, SkimmerModule, HistogramsModule

from bamboo.analysisutils import parseAnalysisConfig

import definitions as defs
import utils

class baseModule(NanoAODModule):
    """ Base class for all plotters, prepares everything needed for both reco-level and gen-level analysis """
    def __init__(self, args):
        super().__init__(args)

        self.doSysts = self.args.systematic
        # FIXME find a better way of doing that... add an argument?
        # either "separate" (muR/muF variations)
        # or "combined" (7-point envelope)
        self.qcdScaleVarMode = "separate"

    def getQcdScaleVariations(self, sample, sampleCfg):
        """Just to know what histogram envelopes to take during postprocessing:
            only for MC, if systematics are enabled, if the sample has no buggy weights,
            and if we take the QCD scale envelope!"""
        if self.doSysts and self.isMC(sample) and self.qcdScaleVarMode == "combined":
            return [ (i, f"qcdScalevar{i}") for i in [0, 1, 3, 5, 7, 8] ]
        return []

    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument("--test", action="store_true", help="Test max. 1 MC and 1 data file (to test before launching a large job submission)")
        parser.add_argument("-s", "--systematic", action="store_true", help="Produce systematic variations")
        parser.add_argument("--syst-samples", action="store_true", help="Process systematic variation samples")
        parser.add_argument("-r", "--reduce-split", type=int, default=0, help="Reduce number of jobs by factor X")
        parser.add_argument("--samples", nargs='*', required=True, help="Sample template YML file")
        parser.add_argument("--roc-corr", action="store_true", help="Enable muon Rochester correction")
        parser.add_argument("--pdf-mode", choices=["simple", "full"], default="simple", help="Full PDF uncertainties (100 histogram variations) or simple (event-based envelope) (only if systematics enabled) (default: %(default)s)")
        parser.add_argument("--top-pt", action="store_true", help="Apply top pt reweighting (and add a 'noTopPt' variation where it is not applied)")

    def _loadSampleLists(self, analysisCfg):
        # fill sample template using JSON files
        if "samples" not in analysisCfg:
            eras = self.args.eras[1]
            samples = {}
            # make sure we use absolute paths as this argument will be used by the worker jobs
            self.args.samples = [ os.path.abspath(p) for p in self.args.samples ]
            for tmpPath in self.args.samples:
                with open(tmpPath) as f_:
                    template = yaml.load(f_, Loader=yaml.SafeLoader)
                    samples.update(utils.fillSampleTemplate(template, eras))
            analysisCfg["samples"] = samples

    def customizeAnalysisCfg(self, analysisCfg):
        self._loadSampleLists(analysisCfg)
        samples = analysisCfg["samples"]

        # reduce job splitting
        if self.args.reduce_split:
            for smp in samples.values():
                smp["split"] *= self.args.reduce_split

        # if we're not doing systematics, remove the systematics samples from the list
        if not (self.doSysts and self.args.syst_samples):
            for smp in list(samples.keys()):
                if "syst" in samples[smp]:
                    samples.pop(smp)

        if not self.args.distributed or self.args.distributed != "worker":
            if self.args.test:
                # only keep 1 MC (if possible, a signal) and 1 data file, for testing the plotter
                chosenEra = self.args.eras[1][0] if self.args.eras[1] else None
                foundMC = utils.getAnyMCSample(samples, self.isMC, era=chosenEra, signalSample=True)
                if foundMC is None:
                    logger.info("No signal sample found for testing, falling back on background sample")
                    foundMC = utils.getAnyMCSample(samples, self.isMC, era=chosenEra, signalSample=False)
                if foundMC is None:
                    logger.warning("No MC sample found for testing!")
                else:
                    logger.info(f"Found MC sample for testing: {foundMC}")
                    chosenEra = samples[foundMC]["era"]
                foundData = utils.getAnyDataSample(samples, self.isMC, era=chosenEra)
                if foundData is None:
                    logger.warning("No data sample found for testing!")
                else:
                    logger.info(f"Found data sample for testing: {foundData}")
                    if chosenEra is None:
                        chosenEra = samples[foundData]["era"]
                for smpNm in list(samples.keys()):
                    if smpNm != foundMC and smpNm != foundData:
                        samples.pop(smpNm)
                # only keep 1 file per sample
                self.args.maxFiles = 1
                # adjust the eras in the analysis config
                for era in list(analysisCfg["eras"].keys()):
                    if era != chosenEra:
                        analysisCfg["eras"].pop(era)
                logger.info(f"Testing mode: only using one file; only running on era {chosenEra}; using data: {foundData}; using MC: {foundMC}")

            # back up analysis config - not really needed since git config is stored, but JIC
            newCfg = os.path.join(self.args.output, "full_analysis.yml")
            os.makedirs(self.args.output, exist_ok=True)
            yaml.add_representer(str, utils.yaml_latex_representer)
            with open(newCfg, "w") as f_:
                yaml.dump(analysisCfg, f_, default_flow_style=False, indent=4)

    def prepareTree(self, tree, sample=None, sampleCfg=None, splitTT=True):
        era = sampleCfg["era"]
        tree,noSel,be,lumiArgs = super().prepareTree(tree, sample=sample, sampleCfg=sampleCfg, description=defs.getNanoAODDescription(era, self.isMC(sample), doRocCor=self.args.roc_corr), backend="lazy")

        be.addDependency(includePath=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "include"))
        be.addDependency(headers="HistogramEvaluator.h")

        if self.isMC(sample):
            # if it's a systematics sample, turn off other systematics
            sample_doSysts = self.doSysts
            if "syst" in sampleCfg:
                sample_doSysts = False

            noSel = noSel.refine("mcWeight", weight=tree.genWeight, autoSyst=sample_doSysts)

            if sample_doSysts:
                noSel = utils.addTheorySystematics(self, sample, sampleCfg, tree, noSel,
                                                   pdf_mode=self.args.pdf_mode,
                                                   # only add additional top pt reweighting variation
                                                   # if we're not applying it everywhere
                                                   topPt=not self.args.top_pt)

            if "subprocess" in sampleCfg and self.args.top_pt:
                # apply top pt reweighting for every variation
                noSel = utils.applyTopPtReweighting(tree, noSel)

            if "subprocess" in sampleCfg and splitTT:
                logger.info(f"Adding ttbar category cuts for {sampleCfg['subprocess']}")
                noSel = utils.splitTTjetFlavours(sampleCfg, tree, noSel)

        else: ## DATA
            pass

        return tree,noSel,be,lumiArgs

    def defineObjects(self, tree, noSel, sample=None, sampleCfg=None):
        pass

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None,
                    makeBU=True, removeBatch=True):
        if makeBU:
            # create backup directory with all merged results
            # if something screws up in the postprocessing, just copy it back as `results`
            resultsdir_bu = os.path.join(workdir, "results_backup")
            if not os.path.isdir(resultsdir_bu):
                shutil.copytree(resultsdir, resultsdir_bu)

        if removeBatch:
            # remove batch output files (not needed anymore since they've been copied/merged into the results directory)
            batchOut = os.path.join(workdir, "batch", "output")
            if os.path.isdir(batchOut):
                shutil.rmtree(batchOut)


class basePlotter(baseModule, HistogramsModule):
    """ Base class for all plotters, prepares everything needed for both reco-level and gen-level analysis """
    def __init__(self, args):
        super().__init__(args)

        self.plotDefaults = {
            "y-axis"           : "Events",
            "log-y"            : "both",
            "y-axis-show-zero" : True,
            "save-extensions"  : ["pdf"],
            "show-ratio"       : True,
            "sort-by-yields"   : True,
        }

        # used for yields tables/cutflowreports
        self.yieldsPrefix = "yields"

    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument("--main-signal", default="powheg_4FS", help="Main signal sample, used for control plots (default: %(default)s)")

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None, 
                    makeBU=True, removeBatch=True, createEnvelope=True, moveSystHists=True, removeSignalOverlap=True, 
                    sampleHint=None, fileHint=None):
        super().postProcess(taskList, config, workdir, resultsdir, makeBU, removeBatch)

        # always use a MC file to get plot list, and a signal sample if possible
        chosenEra = self.args.eras[1][0] if self.args.eras[1] else None
        if sampleHint is None:
            sampleHint = utils.getAnyMCSample(config["samples"], self.isMC, era=chosenEra)
            if sampleHint is None:
                sampleHint = utils.getAnyMCSample(config["samples"], self.isMC, era=chosenEra, signalSample=False)
            fileHint = None
        if sampleHint is not None and fileHint is None:
            fileHint = self.sampleFilesResolver(sampleHint, config["samples"][sampleHint])[0]

        self.plotList = self.getPlotList(resultsdir=resultsdir, config=config, sampleHint=sampleHint, fileHint=fileHint)

        if createEnvelope and self.doSysts:
            for task in taskList:
                if self.isMC(task.name) and "syst" not in task.config:
                    # create QCD scale envelopes - if needed (only if 7-point variations are done)
                    qcdScaleVariations = self.getQcdScaleVariations(task.name, task.config)
                    utils.produceMEScaleEnvelopes(self.plotList, qcdScaleVariations, os.path.join(resultsdir, task.outputFile))
                    if self.args.pdf_mode == "full" and task.config.get("pdf_full", False):
                        utils.producePDFEnvelopes(self.plotList, task, resultsdir)

        if moveSystHists:
            # renormalize, copy and rename histograms from systematic variation files into nominal files
            # also remove the systematic samples from the plotIt list

            if self.doSysts:
                utils.postProcSystSamples(taskList, config["samples"], resultsdir, cutFlowPrefix=self.yieldsPrefix, readCounters=self.readCounters)

        if removeSignalOverlap:
            # remove overlap between alternative signal contributions for plotIt
            utils.removeSignalOverlapPlotIt(self.args.main_signal, config["samples"])

    def runPlotIt(self, taskList, config=None, workdir=None, resultsdir=None):
        # run plotIt as defined in HistogramsModule
        HistogramsModule.postProcess(self, taskList, config, workdir, resultsdir)

class baseSkimmer(baseModule, SkimmerModule):
    pass
