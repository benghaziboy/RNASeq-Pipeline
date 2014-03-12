import os
import sys
import subprocess
import sqlite3 as db
import random

from collections import defaultdict
from datetime import datetime

from . import fs_autocomplete


'''

v1 needs:
    - sqlite3
    - ftp program
    - tar extraction program
    - manage multiple jobs at one time


'''


class Bowtie:

    def __init__(
            self,
            mismatches=3,
            report_all_alignments=True,
            suppress_alignments_above=25,
            query=None,
            reference=None,
            output=None,
            bowtieFolder=None,
            indexFolder=None,
            **kwargs):

        self.query = query
        self.reference = reference
        self.output = output
        self.mismatches = mismatches
        self.report_all_alignments = report_all_alignments
        self.suppress_alignments_above = suppress_alignments_above
        self.bowtieFolder = bowtieFolder
        self.defaultIndex = indexFolder
        self.duration = None

    def updateKwargs(self, kwargs):
        '''set kwargs as attributes of the Bowtie object'''
        for key, value in kwargs.iteritems():
            setattr(self, key, value)
        return

    def gui(self):
        '''user prompt'''
        self.query = fs_autocomplete.get_input(
            "please enter in a FASTA file: ")
        self.reference = fs_autocomplete.get_input(
            "Please enter in a reference file: ")
        self.output = fs_autocomplete.get_input(
            "Please enter in the output file: ")
        self.runBowtie()

    def runBowtie(self):
        '''call Bowtie function
            report - report all sequenced alignemnts
            fastAQ - determines input filetype
            checkReference & checkOutput confirm inputs are valid    
        '''
        startTime = datetime.now()
        app = os.path.join(self.bowtieFolder, 'bowtie')
        report = ''
        fastAQ = self.checkQuery(self.query)
        self.checkReference(self.reference)
        self.checkOutput(self.output)
        if self.report_all_alignments:
            report = '-a'
        cmd = [app,
               fastAQ,
               '-p',
               '8',
               '-v',
               str(self.mismatches),
               report,
               '-m',
               str(self.suppress_alignments_above),
               self.reference,
               self.query,
               self.output]
        run = subprocess.call(cmd, stdout=subprocess.PIPE)
        self.duration = datetime.now() - startTime()
        return

    def checkQuery(self, Query):
        '''reads query file to determine if the file is FASTA (-f) or FASTQ (-q)'''
        with open(Query) as fh:
            fchar = fh.readline().strip()[0]
            if fchar == '>':
                return '-f'
            elif fchar == '@':
                return '-q'
            else:
                while True:
                    resp = fs_autocomplete.get_input(
                        "The provided query file is neither FASTA or FASTQ. Please enter a new one: ")
                    self.query = resp
                    self.checkQuery(self.query)

    def checkOutput(self, Output):
        '''recursively steps through the output path, creating any missing directories'''
        while True:
            if os.path.isdir(os.path.dirname(Output)):
                return
            else:
                try:
                    os.mkdir(os.path.dirname(Output))
                except:
                    self.checkOutput(os.path.dirname(Output))

    def checkReference(self, Reference):
        '''check if the provided reference already exists or needs to be created.
        If reference is not present, the "bowtie-build" application is run on the raw file'''
        ReferenceSp = Reference.split('/')
        if self.folderReferenceCheck(self.defaultIndex, ReferenceSp[-1]):
            return
        elif self.folderReferenceCheck(os.path.dirname(Reference), ReferenceSp[-1]):
            return
        else:
            self.bowtieBuildInit(Reference)
            return

    def folderReferenceCheck(self, Folder, ReferenceName):
        '''checks provided index folder for reference file'''
        indexList = os.listdir(Folder)
        if '.'.join([ReferenceName, '1', 'ebwt']) in indexList:
            self.reference = os.path.join(self.defaultIndex, ReferenceName)
            return True
        else:
            return False

    def bowtieBuildInit(self, Reference):
        '''asks user if they want to build a bowtie index from the provided file and where they want it to be stored'''
        if os.path.isfile(Reference):
            while True:
                buildQ = raw_input(
                    "A pre-existing bowtie index for the provided file could not be found. Do you want to build one? [y/n]: ")
                buildQ = buildQ.strip().lower()
                if buildQ == "y":
                    folderQ = fs_autocomplete.get_input(
                        "What folder do you want to store the reference file in? (Leave blank for default [%s]): " %
                        self.defaultIndex)
                    if folderQ == "":
                        folderQ = self.defaultIndex
                    nameQ = raw_input(
                        "What would you like to call the index file?: ")
                    self.runBowtieBuild(
                        Reference,
                        os.path.join(
                            folderQ,
                            nameQ))
                    return
                elif buildQ == "n":
                    print "GoodBye."
                    sys.exit[0]

    def runBowtieBuild(self, Reference, Output):
        '''create Bowtie Reference'''
        app = os.path.join(self.bowtieFolder, 'bowtie-build')
        subprocess.call([app, Reference, Output])
        self.reference = Output
        return
