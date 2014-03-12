import os
import sys


class Splitter(object):
    '''Once aggregate files became too cumbersome, it became necessary to split the output using varied split points
    receives a tsv argument and library # to determine at which point the file should be divided'''

    def __init__(
            self,
            entryfile=None,
            cutofflibrary=100,
            nfiles=2,
            *args,
            **kwargs):
        self.entryfile = entryfile
        self.cutofflibrary = cutofflibrary
        self.nfiles = nfiles
        self.split_list = list()
        self.leftend, self.splitpoint, self.rightstart = self._get_header_ranges(
        )
        self._get_new_files()
        self._write_files()

    def _write_files(self):
        with open(self.entryfile) as fh:
            with open(self.split_list[0], 'w') as fw1:
                with open(self.split_list[1], 'w') as fw2:
                    for line in fh:
                        sp = line.strip().split('\t')
                        fw1.write(
                            '\t'.join(sp[:self.splitpoint] + sp[self.rightstart:] + ["\n"]))
                        fw2.write(
                            '\t'.join(sp[:self.leftend - 1] + sp[self.splitpoint:] + ["\n"]))
        return

    def _get_new_files(self):
        i = 1
        while i <= self.nfiles:
            fp = self.entryfile.rpartition('.')
            newname = "%s.%s.%s" % (fp[0], str(i), fp[2])
            self.split_list.append(newname)
            i += 1
        return

    def _get_list(self):
        return self.split_list

    def _get_header_ranges(self):
        leftend, splitpoint, rightstart = 1, 0, 0
        with open(self.entryfile) as fh:
            past_libraries = False
            headers = fh.readline().split('\t')
            for key in headers:
                rightstart += 1
                if "R%s" % self.cutofflibrary in key:
                    splitpoint = headers.index(key)

                if self._is_library(key):
                    past_libraries = True
                    continue
                else:
                    if past_libraries:
                        return leftend, splitpoint, rightstart
                    else:
                        leftend += 1

    def _is_library(self, header):
        '''checks if provided header is a library'''
        id = header.split('_')
        if len(id) < 2:
            return False
        elif id[0][0] != 'R':
            return False
        else:
            try:
                id_num = int(id[0][1])
                return True
            except:
                return False
