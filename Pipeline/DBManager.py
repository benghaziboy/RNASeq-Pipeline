import os
import sys
import sqlite3 as db

from collections import defaultdict

from . import PathCheck
from . import fs_autocomplete


class DBM:

    def __init__(self):
        self.pipeline_options = {
            "-v": "mismatches",
            "Bowtie Path": "Folder containing bowtie application",
            "Reference Path": "path to standardized index for pipeline",
            "-m": "suppress all alignments above (int) value",
            "-a": "report all valid alignments per read",
            "-p": "number of cores used by bowtie",
            "Bowtie Output": "Folder for the bowtie output",
            "RPKM Output": "Folder for the rpkm output",
            "Aggregate Output": "Folder for the aggregate file output",
            "Annotation Path": "Location of annotation files",
            "Split": "(int) value of which library to split the RPKM/Hit aggregate files",
            "Raw Path": "Folder containing raw RNASeq files",
            "RPKM Path": "Folder containing pre-existing RPKM files",
            "GlymaFile": "File used for RPKM calculation",
            "Glyma Column": "Column in GlymaFile that contains the sequence",
            "Bowtie Column": "Column in Bowtie File that contains the sequence",
        }
        self.initializeSqlite()

    def initializeSqlite(self):
        self.conn = db.connect(
            os.path.join(
                os.path.split(__file__)[0],
                'rnaseq.sqlite'))
        self.curs = self.conn.cursor()
        self.initializeDatabase()

    def checkFolderOptions(self):
        for option in ["Bowtie Output", "RPKM Output", "Aggregate Output", "RPKM Path"]:
            if option == "Aggregate Output":
                [PathCheck.check_folder(os.path.join(self.config[option], s))
                 for s in ["Hits/", "RPKMs/"]]
#                 [PathCheck.check_folder(os.path.join(self.config[option], s) for s in ["Hits", "RPKMs"])]
            else:
                PathCheck.check_folder(str(self.config[option]))
        return

    def initializeDatabase(self):
        try:
            self.curs.execute("select * from pipelineconfig")
        except Exception as e:
            self.curs.execute(
                "create table pipelineconfig (optionname VARCHAR(1000), optionvalue VARCHAR(1000));")
        config = defaultdict(lambda: "")
        self.curs.execute("select * from pipelineconfig")
        for row in self.curs:
            config[row[0]] = row[1]
        self.config = config

    def updateOption(self, option):
        self.getOption(option)
        self.curs.execute("delete from pipelineconfig")
        self.conn.commit()
        self.curs.execute(
            "insert into pipelineconfig values ('%s', '%s');" %
            (option, self.config[option]))
        self.conn.commit()

    def getOptions(self):
        [self.getOption(option) for option in self.pipeline_options.keys()]
        self.curs.execute("delete from pipelineconfig")
        self.conn.commit()
        for option in self.config:
            self.curs.execute(
                "insert into pipelineconfig values ('%s', '%s');" %
                (option, self.config[option]))
        self.conn.commit()
        self.checkFolderOptions()

    def getOption(self, option):
        resp = fs_autocomplete.get_input(
            "%s [%s]: " %
            (option, self.config[option]))
        if resp:
            self.config[option] = resp
        while not self.checkOption(option):
            print "Please try again."
            self.config[option] = fs_autocomplete.get_input(
                "%s: [] " %
                (option, ))
        return

    def checkOption(self, option):
        if self.config[option] is None or self.config[option] == '':
            return False
        if option.lower() in ["Bowtie Path", "Reference Path"]:
            return os.path.isdir(
                self.config[option]) or os.path.isfile(
                self.config[option])
        if option.lower() in ["Raw Path", "Annotation Path"]:
            return os.path.isdir(self.config[option])
        if option == "-v":
            try:
                return 0 <= int(self.config[option]) <= 3
            except:
                return False
        if option == "-m" or option == "-p" or option == "Split":
            try:
                return isinstance(int(self.config[option]), int)
            except:
                return False
        if option == "-a":
            if self.config[option].lower() in ['true', 't', 'y', 'yes']:
                self.config[option] = True
                return True
            else:
                return False
        else:
            return True

    def printOptions(self):
        for key, value in self.config.iteritems():
            print "[%s]: %s" % (key, value)
        return

DBM = DBM()
config = DBM.config
