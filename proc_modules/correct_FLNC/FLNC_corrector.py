#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time : 2020/5/6 17:49
# @Author : Yuancong Wang
# @Site : 
# @File : FLNC_corrector.py
# @Software: PyCharm
"""
Description: This routine using high-throughput short-read sequencing data to correct the long-read sequencing data.
Output: The corrected long-read sequencing data (corrected FLNC CCSs)
Input: path to the RNA-seq data file; path to the FLNC CCSs file; path to the prooveread software
Other notes: using proovread
"""

import os
import subprocess
from funcs.funcs import get_abspath_from_input_path


class FlncCorrect:

    def __init__(self):
        self._flnc_fpath = None
        self._flnc_fa_fpath = None
        self._out_dpath_uppername = None
        self._short_read_fpaths = None
        self._proovread_fpath = None
        self._movie_name = None
        self._out_dpath = None

    def fit(self, flnc_fpath, short_read_fpaths, proovread_fpath):

        for spath in short_read_fpaths:
            assert os.path.exists(spath), \
                'The short-read file is not existed, please make sure that you give the right path of each file.'
        self._short_read_fpaths = short_read_fpaths
        assert os.path.exists(flnc_fpath), \
            'The FLNC CCS file is not existed, please make sure that you give the right path of the file.'
        self._flnc_fpath = flnc_fpath
        self._movie_name = os.path.basename(flnc_fpath).split('.')[0]
        if not proovread_fpath:
            self._proovread_fpath = 'proovread'
        else:
            self._proovread_fpath = get_abspath_from_input_path(proovread_fpath)
        self._out_dpath_uppername = os.path.dirname(self._flnc_fpath)
        return self

    def bam2fa(self):
        """
        convert the .bam to .fasta format
        :return: path of the fasta file
        """
        out_file_name = self._movie_name + '.flnc.fasta'
        self._flnc_fa_fpath = os.path.join(self._out_dpath_uppername, out_file_name)
        if not os.path.exists(self._flnc_fa_fpath):
            cmd_bam2fa = '''
            samtools view {flnc_bam_path} | awk \
                    '{{OFS="\\t"; print ">"$1"\\n"$10}}' - > {flnc_fa_path}
            '''.format(flnc_bam_path=self._flnc_fpath, flnc_fa_path=self._flnc_fa_fpath)
            if subprocess.check_call(cmd_bam2fa, shell=True) != 0:
                raise SystemCommandError
        return self

    def correct_flnc(self):
        """
        using short-read sequencing data (RNA-seq) data to correct long-read sequencing data
        :return: the directory path of the corrected sequences file
        """
        prooread_out_dpath_basename = self._movie_name + '_pread'
        self._out_dpath = os.path.join(self._out_dpath_uppername, prooread_out_dpath_basename)
        cmd_correct_seq = '{path_of_proovread} -l {long_read_seq_fpath} -s {short_read_seq_fpaths} \
            -p {pre} > proovread.log'.format(path_of_proovread=self._proovread_fpath,
                                             long_read_seq_fpath=self._flnc_fa_fpath,
                                             short_read_seq_fpaths=' '.join(self._short_read_fpaths),
                                             pre=self._out_dpath)
        if subprocess.check_call(cmd_correct_seq, shell=True) != 0:
            raise SystemCommandError
        return self._out_dpath
