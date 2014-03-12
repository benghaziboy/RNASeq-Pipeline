'''takes all RNASeq rpkm files, aggregates, and annotates them'''
import os
import sys

from collections import defaultdict
from datetime import datetime
import gc

from . import DBManager as DB
from .DefaultList import DefaultList
from . import PathCheck

import pdb


class MasterRPKM(object):
    '''reads the rpkm files for every library, aggregates, and annotates them'''

    def __init__(self, *args, **kwargs):
        self.config = DB.DBM.config
        self.attrs = {
            'masterfile': str(
                os.path.join(
                    self.config["Aggregate Output"],
                    'RNASeq-%s' %
                    datetime.now().strftime("%H%M-%m%d%y"))),
            'annotationpath': self.config["Annotation Path"],
            'rpkm_path': self.config["RPKM Path"],
            'newest_library': latest_lib,
            'masterHit': str(
                os.path.join(
                    os.path.join(
                        self.config["Aggregate Output"],
                        "Hits"),
                    'RNASeq-%s' %
                    datetime.now().strftime("%H%M-%m%d%y"))),
            'masterRPKM': str(
                os.path.join(
                    os.path.join(
                        self.config["Aggregate Output"],
                        "RPKMs"),
                    'RNASeq-%s' %
                    datetime.now().strftime("%H%M-%m%d%y"))),
        }
        self.kwargs = kwargs
        self.defaultInit()
        self.manageKwargs()
        self.missinglibraries = list()
        self.sortedlibraries = list()

    def defaultInit(self):
        '''sets default for various attributes from the given dictionary'''
        for key, value in self.attrs.iteritems():
            setattr(self, key, value)
        return

    def manageKwargs(self):
        for key, value in self.kwargs.items():
            if key in ["bowtie_model_column", "glyma_model_column"]:
                setattr(self, key, value)
            else:
                setattr(self, key, str(value))
        return

    def _get_RNASeq_Dict(self):
        """"returns a dictionary with the keys being RNASeq Library Numbers and the values being a tuple of the id # and filename
            e.g., RNASeqDict['R01'] = (01, 'R01_RNASeq1_W43')"""
        RNASeqDict = defaultdict(lambda: (' ', ' '))
        for key in os.listdir(self.rpkm_path):
            sp = key.split('_')
            id_number = int(sp[0].lstrip('R'))
            if id_number > self.newest_library:
                self.newest_library = id_number
            RNASeqDict[sp[0]] = ("_".join(sp[:3]), key)
            if id_number > self.newest_library:
                self.newest_library = id_number
        return RNASeqDict

    def _get_annotation_Dicts(self):
        '''creates dictionaries for each annotation file except for the phytozome file "Gmax_109_annotation_info.txt" which contains a unique annotation format'''
        annotation_ref = {
            "nr": "nr_annotation.tsv",
            "swiss": "Swiss_annotation.tsv",
            "trembl": "trEMBL_annotation.tsv",
            "keyword": "Keywords.tsv",
            "pfam": "PFAMAnnotation.tsv",
            "phyto": "v6/Gmax_109_annotation_info.txt"
        }
        annotation_files = dict(zip(annotation_ref.keys(), [
                                os.path.join(self.annotationpath, f) for f in annotation_ref.values()]))
        annotation_dict = defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(
                    lambda: DefaultList(
                        default="None"))))
        annotation_dict = self._get_phytozome_annotations(
            annotation_files["phyto"],
            annotation_dict)
        for annogroup, filename in annotation_files.iteritems():
            if annogroup != "phyto":
                with open(filename) as fh:
                    fh.readline()
                    lmodel = ''
                    for line in fh:
                        annotation = "None"
                        sp = line.strip().partition("\t")
                        model, content = sp[0].strip(), sp[2].strip()
                        if model != lmodel:
                            i = 1
                            lmodel = model
                        if content != "" or content != " ":
                            annotation = content
                        annotation_dict[model][i][annogroup].append(annotation)
                        i += 1
            else:
                continue
        for model in annotation_dict:
            for it in annotation_dict[model]:
                if 'phyto' not in annotation_dict[model][it]:
                    annotation_dict[model][it]['phyto'] = ["None"] * 8
        return annotation_dict

    def _get_phytozome_annotations(self, filename, annotations_dict):
        with open(filename) as fh:
            fh.readline()
            anno = ["None"] * 8
            for line in fh:
                sp = line.strip().split('\t')
                model = sp[0]
                for i, x in enumerate(sp[1:]):
                    if x.strip() != '':
                        anno[i] = x
                annostring = '\t'.join(anno)
                for i in range(1, 11):
                    annotations_dict[model][i]["phyto"] = anno
        return annotations_dict

    def _get_model_lengths(self):
        cds_file = os.path.join(self.annotationpath, "cds_Length.tsv")
        cdna_file = os.path.join(self.annotationpath, "cDNA_Length.tsv")
        cds_dict = defaultdict(lambda: 0)
        cdna_dict = defaultdict(lambda: 0)
        with open(cds_file) as fh:
            fh.readline()
            for l in fh:
                sp = l.strip().split("\t")
                cds_dict[sp[0]] = sp[1]
        with open(cdna_file) as fh:
            fh.readline()
            for l in fh:
                sp = l.strip().split("\t")
                cdna_dict[sp[0]] = sp[1]
        return cds_dict, cdna_dict

    def _get_rpkm_dict(self):
        RNASeq_libs = self._get_RNASeq_Dict()
        rpkm_dict = dict(zip([r[0] for r in RNASeq_libs.values()], [
                         defaultdict(lambda:0) for r in RNASeq_libs.values()]))
        hit_dict = dict(zip([r[0] for r in RNASeq_libs.values()], [
                        defaultdict(lambda:0) for r in RNASeq_libs.values()]))
        LibNum = 1
        print self.newest_library
        while LibNum <= self.newest_library:
            try:
                RNASeq_id = "R%s" % str(LibNum).zfill(2)
                RNASeq_name = RNASeq_libs[RNASeq_id][0]
                RNASeq_filepath = os.path.join(
                    self.rpkm_path,
                    RNASeq_libs[RNASeq_id][1])
                with open(os.path.join(RNASeq_filepath)) as fh:
                    fh.readline()
                    for line in fh:
                        sp = line.strip().split("\t")
                        model, rpkm, hit = sp[0], sp[1], sp[2]
                        rpkm_dict[RNASeq_name][model] = rpkm
                        hit_dict[RNASeq_name][model] = hit
                self.sortedlibraries.append(RNASeq_name)
            except:
                self.missinglibraries.append("R%s" % str(LibNum).zfill(2))
            LibNum += 1
        return rpkm_dict, hit_dict

    def _run_Aggregate(self):
        '''create aggregate files depending on if writeHits & writeRPKM are True'''
        annotation_dict = self._get_annotation_Dicts()
        cds_dict, cdna_dict = self._get_model_lengths()
        rpkmdicts, hitdicts = self._get_rpkm_dict()
        pre_header = ["Model", "cds Length", "cDNA Length", "Hit Number"]
        sorted_annotation = [
            "keyword",
            "nr",
            "swiss",
            "trembl",
            "phyto",
            "pfam"]
        app_header = [
            "Keywords",
            "nr Annotation",
            "Swiss annotation",
            "trEMBL annotation",
            "PFAM annotation",
            "Panther annotation",
            "KOG annotation",
            "KEGG ec",
            "KEGG Orthology",
            "Best Arabidopsis Hit name",
            "Best Arabidopsis Hit symbol",
            "Best Arabidopsis Hit defline",
            "PFAM Annotation"]
        with open(self.masterHit + ".Hits.R%s.tsv" % str(self.newest_library).zfill(2), 'w') as Hitfh:
            with open(self.masterRPKM + ".RPKM.R%s.tsv" % str(self.newest_library).zfill(2), 'w') as RPKMfh:
                Hitfh.write("\t".join(pre_header +
                                      [l +
                                       " Hits" for l in self.sortedlibraries] +
                                      app_header) +
                            "\n")
                RPKMfh.write("\t".join(
                    pre_header + [l + " RPKMs" for l in self.sortedlibraries] + app_header) + "\n")
                for model in rpkmdicts[rpkmdicts.keys()[0]].keys():
                    for i in range(1, 11):
                        Hitfh.write("\t".join([model, cds_dict[model], cdna_dict[model], str(i)] +
                                              [hitdicts[k][model] for k in self.sortedlibraries] +
                                              ['\t'.join(annotation_dict[model][i][k]) for k in sorted_annotation]) +
                                    "\n")
                        RPKMfh.write("\t".join([model, cds_dict[model], cdna_dict[model], str(i)] +
                                               [rpkmdicts[k][model] for k in self.sortedlibraries] +
                                               ['\t'.join(annotation_dict[model][i][k]) for k in sorted_annotation]) +
                                     "\n")
        gc.collect()

        return
