#!/usr/bin/env python

import os
import sys
import gc
import subprocess
import multiprocessing
import argparse

from multiprocessing.dummy import Pool
from datetime import datetime
import signal

from . import DBManager as DB
from .RNASeq import RNASeqManager
from .Bowtie import Bowtie
from .RPKM import RPKMs
from .MasterRPKM import MasterRPKM
from .tsv_splitter import Splitter
from . import fs_autocomplete


class Pipeline():
    '''The RNASeq Pipeline was created to automate the analysis & sequencing of high-volume RNASeq data (164 libraries as of January 31st, 2014)
    the Pipeline utilizes sqlite3 to log run information & the relevant paths/variables necessary for each run.
    The final product of running the pipeline is a collection of individual RPKM/Hit files for each library and an aggregate RPKM file and an aggregate Hit file
    representing the values calculated by the Pipeline.

    1. All raw data is run against the database of existing RPKMS. All missing libraries will be aligned to the specified Glymafile.
    2. The Bowties are aligned using values entered by the user in the prompt.
    3. The Hit values of each model are counted in order to calculate the RPKM values.
    4. The RPKM/Hit values are written into their own respective files as well as an aggregate ".tsv" file.
    (5.) The User may choose to then apply their tsv_splitter function to allow for the import restraints of microsoft excel.

    ALL FUNCTIONS MAY BE RUN INDIVIDUALLY
    as the amount of data and the time necessary to process each library has drastically increased, errors with storage space & memory have become increasingly
    problematic. For this reason, all functions of the pipeline can be run individually to avoid wasting time on redundant functions.
    '''

    def __init__(self, *args, **kwargs):
        self.config = DB.DBM.config
        self.pool = Pool(processes=3)
        self.manager = dict()
        self.RNA = RNASeqManager()

    def runPipeline(self):
        startTime = datetime.now()
        self.queueBowties()
        self.queueRPKMs()
        self.Aggregate()
        duration = datetime.now() - startTime
        print "Pipeline took:"
        print(duration)

    def queueBowties(self):
        self.bowtiePool = Pool(processes=3)
        self.manager = self.getBowties()
        results = []
        startTime = datetime.now()
        for lib, bwt in self.manager.iteritems():
            print lib
            result = self.bowtiePool.apply_async(bwt.runBowtie)
            results.append(result)
        self.bowtiePool.close()
        self.bowtiePool.join()
        return

    def getBowties(self):
        cmds = dict()
        for lib in self.RNA.MissingLibraries:
            library = self.RNA.RNASeqDir[lib]
            query_path = os.path.join(library.raw_folder, library.raw_file)
            output_path = os.path.join(
                self.config["Bowtie Output"],
                library.raw_file +
                '.albwt')
            library.bowtie_file = output_path
            bowtie = Bowtie(
                mismatches=self.config["-v"],
                report_all_alignments=self.config["-a"],
                suppress_alignments_above=self.config["-m"],
                query=query_path,
                reference=self.config["Reference Path"],
                output=output_path,
            )
            cmds[lib] = bowtie
        return cmds

    def getRPKMs(self):
        cmds = dict()
        for lib in self.RNA.MissingLibraries:
            library = self.RNA.RNASeqDir[lib]
            rpkm = RPKMs(bowtiefile=library.bowtie_file,
                         glymafile=self.config["GlymaFile"],
                         outputfile=os.path.join(self.config["RPKM Output"],
                                                 library.bowtie_file.rpartition('/')[-1].rstrip('.albwt') + '.rpkm'),
                         libraryname=library.name,
                         bowtie_model_column=self.config["Bowtie Column"],
                         glyma_model_column=self.config["Glyma Column"],
                         )
            cmds[lib] = rpkm
        return cmds

    def queueRPKMs(self):
        for lib, rpkm in self.getRPKMs().iteritems():
            try:
                rpkm.runRPKM()
            except:
                pdb.set_trace()
        return

    def RPKMbyDirectory(self, directory=None):
        '''finds all bowtie files in a directory and turns them into rpkm files'''
        if directory is None:
            directory = fs_autocomplete.get_input(
                "Please enter the bowtie directory: ")
        for key in os.listdir(directory):
            if key.rpartition('.')[-1] == 'albwt':
                print self.config["RPKM Output"]
                if self.config["RPKM Output"] is None or self.config == '':
                    print "no RPKM directory"
                    break
                rpkm = RPKMs(
                    bowtiefile=os.path.join(
                        directory,
                        key),
                    glymafile=self.config["GlymaFile"],
                    outputfile=os.path.join(
                        self.config["RPKM Output"],
                        key.rpartition('.')[0] +
                        '.rpkm'),
                    libraryname=self.cmdlibname(key),
                    bowtie_model_column=self.config["Bowtie Column"],
                    glyma_model_column=self.config["Glyma Column"],
                )
                print rpkm.outputfile
                rpkm.runRPKM()
        self.pool.close()
        self.pool.join()
        return

    def alarmHandler(self, selfsignum, frame):
        raise AlarmException

    def cmdlibname(self, entry, timeout=10):
        '''allows user to enter their own library name for automated rpkm search'''
        signal.signal(signal.SIGALRM, self.alarmHandler)
        signal.alarm(timeout)
        try:
            libname = raw_input(
                "Please enter the library name for %s: " %
                entry)
            signal.alarm(0)
            return libname
        except AlarmException:
            libname = entry.partition('_')[0]
        signal.signal(signal.SIGALRM, signal.SIG_IGN)
        return libname

    def Aggregate(self):
        master = MasterRPKM(
            annotationpath=self.config["Annotation Path"],
            rpkm_path=self.config["RPKM Path"],
        )
        master._run_Aggregate()
        return


class AlarmException(Exception):
    pass


class PipeArgs:

    '''manages command line input to determine what functions to run'''

    def __init__(self):
        self.parser = self.parser()

    def parser(self):
        parser = argparse.ArgumentParser()
        functions = parser.add_mutually_exclusive_group()
        functions.add_argument(
            "-pipeline",
            "-pipeline",
            action="store_true",
            help="Run the full RNASeq Pipeline from beginning to end. ")
        functions.add_argument(
            "-RPKM",
            "-rpkm",
            action="store_true",
            help="Only run an individual RPKM (add options -directory or -single to choose between multiple or single runs")
        functions.add_argument(
            "-aggregate",
            "-agg",
            action="store_true",
            help="Aggregates the RPKM & Hit files based on the default RPKM directory")
        functions.add_argument(
            "-split",
            "-Split",
            action="store_true",
            help="splits a provided aggregate file at the default library")
        functions.add_argument(
            "-update",
            "-updateDB",
            action="store_true",
            help="Update the default options for the pipeline")
        rpkmargs = parser.add_mutually_exclusive_group()
        rpkmargs.add_argument(
            "-dir",
            "-directory",
            action="store_true",
            help="create rpkm files for all the BOWTIE(albwt) files in a provided folder")
        rpkmargs.add_argument(
            "-single",
            action="store_true",
            help="create rpkm file for a single BOWTIE(albwt) file")
        parser.add_argument(
            '-view',
            '-viewoptions',
            action="store_true",
            help="print the default options for the pipeline stored in the SQLdatabase")
        if len(sys.argv) == 1:
            parser.print_help()
            sys.exit(1)
        return parser.parse_args()

    def run(self):
        p = self.parser
        if p.view:
            DB.DBM.printOptions()
        if p.pipeline:
            Pipeline().runPipeline()
        if p.RPKM:
            if p.dir:
                Pipeline().RPKMbyDirectory()
            elif p.single:
                rpkm = RPKMs()
                rpkm.cmdRPKM()
        if p.aggregate:
            Pipeline().Aggregate()
        if p.split:
            resp = fs_autocomplete.get_input(
                "Please enter the file you wish to split: ")
            Splitter(entryfile=resp)
        if p.update:
            DB.DBM.getOptions()
        gc.collect()
        return


def main():
    PipeArgs().run()

if __name__ == '__main__':
    main()
