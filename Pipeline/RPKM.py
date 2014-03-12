import os
import sys

from collections import defaultdict

from . import fs_autocomplete
from . import PathCheck

import pdb


'''calculates Hits & Rpkms for RNASeq Libraries

    function:
        m - total # of mapped reads for every sequence
        h - total # of reads for the given sequence/model
        l - length of given sequence

        (h/(l/1000))/(m/1000000)

        '''


class RPKMs(object):
    '''Reads Per Kilobase of Exon Model per Million MappedReads'''

    def __init__(self, *args, **kwargs):
        self.attrs = {
            "bowtiefile": None,
            "glymafile": None,
            "outputfile": None,
            "bowtie_model_column": 2,
            "glyma_model_column": 2,
            "libraryname": None}
        self.defaultInitDict()
        self.kwargs = kwargs
        self.manageKwargs()

    def checklibname(self):
        if self.libraryname is None:
            try:
                self.libraryname = self.bowtiefile.partition('_')[0]
            except:
                self.libraryname = "RNASeq"

    def manageKwargs(self):
        for key, value in self.kwargs.items():
            if key in ["bowtie_model_column", "glyma_model_column"]:
                setattr(self, key, value)
            else:
                setattr(self, key, str(value))
        return

    def defaultInitDict(self):
        '''sets default for various attributes'''
        for key, value in self.attrs.iteritems():
            setattr(self, key, value)
        return

    def getHits(self):
        ''' The hit for each glyma model is calculated by simple reading the alignment output and incrementing the number of successful alignments made by each model.'''
        hitDict = defaultdict(lambda: 0)
        with open(self.bowtiefile) as fh:
            for line in fh:
                try:
                    model = line.strip().split('\t')[int(self.bowtie_model_column)]
                    hitDict[model] += 1
                except:
                    print fh
                    break
        return hitDict

    def getLengths(self):
        ''' The lengths of each sequence is necessary for calculating the RPKM value of each model and thus the file is read
        for each sequence and recorded.'''
        lengthDict = defaultdict(lambda: 0)
        with open(self.glymafile) as fh:
            for line in fh:
                sp = line.strip().split('\t')
                model = sp[0]
                try:
                    lengthDict[model] = len(sp[1])
                except:
                    raise "Cannot get lengths"
        return lengthDict

    def getRPKMs(self, hitDict, lengthDict):
        rpkmDict = defaultdict(lambda: 0)
        mappedReads = sum(hitDict.values())
        mappedReadsPerM = mappedReads / 1000000.0
        for model in hitDict:
            try:
                rpkmDict[model] = hitDict[model] / \
                    (lengthDict[model] / 1000.0) / mappedReadsPerM
            except:
                print "can't create rpkm dict"
                raise
        return rpkmDict

    def writeRPKMs(self, lengthDict, hitDict, rpkmDict):
        print "writing"
        self.checklibname()
        if self.outputfile:
            PathCheck.check_folder(self.outputfile)
        with open(self.outputfile, 'w') as fh:
            fh.write('\t'.join(["Model Name",
                                self.libraryname + " rpkm ",
                                self.libraryname + " hits " + "\n"]))
            for model in lengthDict:
                fh.write(
                    "\t".join([model, str(rpkmDict[model]), str(hitDict[model])]))
                fh.write(os.linesep)
            print self.outputfile
        fh.close()
        return

    def runRPKM(self):
        hitDict = self.getHits()
        lengthDict = self.getLengths()
        rpkmDict = self.getRPKMs(hitDict, lengthDict)
        self.writeRPKMs(lengthDict, hitDict, rpkmDict)
        return

    def cmdRPKM(self):
        '''prompt for user to run individual rpkm'''
        for key, value in self.attrs.iteritems():
            resp = fs_autocomplete.get_input("%s [%s]: " % (key, value))
            if resp.strip() != '':
                setattr(self, key, resp)
        self.runRPKM()
        return


RPKM = RPKMs()


def main():
    RPKMs().cmdRPKM()

if __name__ == '__main__':
    main()
