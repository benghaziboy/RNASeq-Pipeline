import os
import sys

from collections import defaultdict

from .RNASeq import RNASeq


class RNASeqManager:

    def __init__(self, *args, **kwargs):
        '''local settings for RNASeq pipeline in UIUC-Crop Sciences lab'''
        self.attrs = {
            "raw_folder": "/home/vodkinlab/RNASeq/raw_sequences/",
            "rpkm_folder": "/home/vodkinlab/RNASeq/bowtie_v3m25a_nobeststrata_rpkms/",
        }
        self.manageKwargs(kwargs)
        self.initAttrs()
        self.RNASeqDir = defaultdict(lambda: None)
        self.get_complete_RNASeq()

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
        for key in os.listdir(self.raw_folder):
            id = key.split('_')[0]
            self.RNASeqDir[id] = RNASeq(id=id)
        return
