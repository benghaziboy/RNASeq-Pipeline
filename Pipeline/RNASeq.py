import os
import sys

from collections import defaultdict
from . import DBManager as DB


class RNASeq(object):
    '''RNASeq object to record all applied variables and output values'''


    RNASeqCount = 0

    def __init__(self, *args, **kwargs):
        self.config = DB.DBM.config
        self.attrs = {
            "name": None,
            "id": None,
            "raw_file": None,
            "raw_folder": self.config["Raw Path"],
            "rpkm_file": None,
            "rpkm_folder": self.config["RPKM Path"],
            "bowtie_file": None,
        }
        self.manageKwargs(kwargs)
        self.initAttrs()
        self.raw_file = self.get_file_by_id(self.raw_folder)
        self.rpkm_file = self.get_file_by_id(self.rpkm_folder)
        RNASeq.RNASeqCount += 1

    def initAttrs(self):
        for key, value in self.attrs.iteritems():
            if hasattr(self, key):
                continue
            else:
                setattr(self, key, value)
        return

    def manageKwargs(self, kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        return

    def get_file_by_id(self, folder):
        file_list = os.listdir(folder)
        related_files = [
            x for x in file_list if self.id.strip() == x.split("_")[0].strip()]
        if len(related_files) > 1:
            print related_files
            print "Too many files associated with %s in %s. Exiting" % (self.id, folder)
            sys.exit(0)
        elif len(related_files) < 1:
            print "No File found for %s in %s." % (self.id, folder)
            return None
        else:
            return related_files[0]

    def __str__(self):
        return self.id


class RNASeqManager:

    def __init__(self, *args, **kwargs):
        self.attrs = {
            "raw_folder": None,
            "rpkm_folder": None,
        }
        self.manageKwargs(kwargs)
        self.initAttrs()
        self.MissingLibraries = []
        self.RNASeqDir = self.get_complete_RNASeq()

    def initAttrs(self):
        for key, value in self.attrs.iteritems():
            if hasattr(self, key):
                continue
            else:
                setattr(self, key, value)
        return

    def manageKwargs(self, kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        return

    def get_complete_RNASeq(self):
        RNASeqDict = dict()
        for key in os.listdir(self.raw_folder):
            id = key.split('_')[0]
            RNASeqDict[id] = RNASeq(id=id)
            if not RNASeqDict[id].rpkm_file:
                self.MissingLibraries.append(id)
        return RNASeqDict
