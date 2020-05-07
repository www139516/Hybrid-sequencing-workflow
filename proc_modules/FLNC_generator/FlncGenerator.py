#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time : 2020/5/6 10:04
# @Author : Yuancong Wang
# @Site : 
# @File : FlncGenerator.py
# @Software: PyCharm
"""
Description: This class processes the subreads from long-read sequencing to generate full length non-chimeric CCSs
  (FLNC CCSs).
Output: FLNC CCSs file
Input: Subreads.bam file
Other notes: Using IsoSeq3 to process subreads
"""

import os
import subprocess
from funcs.funcs import get_abspath_from_input_path


class FlncGenerator:

    def __init__(self):
        self._subread_fpath = None
        self._primer_fpath = None
        self._out_dpath = None
        self._movie_name = None
        self._out_ccs_fpath = None
        self._out_flccs_fpath = None
        self._in_flnc_fpath = None
        self._out_flnc_fpath = None

    def fit(self, subread_fpath, primer_fpath, out_dpath=''):
        """
        set the values for this class
        :param subread_fpath: the path of the subreads file
        :param primer_fpath: the path of the file containing the primer sequences used in the PacBio-sequencing
        :param out_dpath: the out put directory path
        :return: self
        """
        wkdir = os.getcwd()
        abs_subread_fpath = get_abspath_from_input_path(subread_fpath)
        abs_primer_fpath = get_abspath_from_input_path(primer_fpath)
        if out_dpath:
            abs_out_dpath = get_abspath_from_input_path(out_dpath)
        else:
            abs_out_dpath = os.path.join(wkdir, 'out')

        assert os.path.exists(abs_subread_fpath), \
            'The subread file is not existed, please make sure that you have the subread file.'
        assert os.path.exists(abs_primer_fpath), \
            'The primer file is not existed, please make sure that you have the primer file.'

        if not os.path.exists(abs_out_dpath):
            os.mkdir(abs_out_dpath)

        self._subread_fpath = abs_subread_fpath
        self._primer_fpath = abs_primer_fpath
        self._out_dpath = abs_out_dpath
        self._movie_name = os.path.basename(self._subread_fpath).split(".")[0]

        return self

    def sub2ccs(self):
        """
        generate ccs based on subreads.bam
        :return: None
        """
        out_file_name = self._movie_name + '.ccs.bam'
        self._out_ccs_fpath = os.path.join(self._out_dpath, out_file_name)

        cmd_subreads2ccs = 'ccs {input_subreads} {output} '.format(input_subreads=self._subread_fpath,
                                                                   output=self._out_ccs_fpath) + '--noPolish --minPasses 1'
        if subprocess.check_call(cmd_subreads2ccs, shell=True) != 0:
            raise SystemCommandError

    def ccs2flccs(self):
        """
        generate full length ccs based on ccs
        :return: None
        """
        out_file_name = self._movie_name + '.flccs.bam'
        self._out_flccs_fpath = os.path.join(self._out_dpath, out_file_name)
        cmd_ccs2flccs = 'lima {input_ccs} {barcode} {output_flccs} '.format(input_ccs=self._out_ccs_fpath,
                                                                            barcode=self._primer_fpath,
                                                                            output_flccs=self._out_flccs_fpath) + \
                        '--isoseq --no-pbi --peek-guess'
        if subprocess.check_call(cmd_ccs2flccs, shell=True) != 0:
            raise SystemCommandError

    def flccs2flnc(self):
        """
        generate FLNC CCSs based on FLCCS
        :return: the path of the FLNC CCS file
        """
        in_file_name = self._movie_name + '.flccs.primer_5p--primer_3p.bam'
        self._in_flnc_fpath = os.path.join(self._out_dpath, in_file_name)
        out_file_name = self._movie_name + '.flnc.bam'
        self._out_flnc_fpath = os.path.join(self._out_dpath, out_file_name)
        cmd_flccs2flnc = "isoseq3 refine {input_flccs} {barcode} {output_flnc} ".format(input_flccs=self._in_flnc_fpath,
                                                                                        barcode=self._primer_fpath,
                                                                                        output_flnc=self._out_flnc_fpath) + \
                         "--require-polya"
        if subprocess.check_call(cmd_flccs2flnc, shell=True) != 0:
            raise SystemCommandError

        return self._out_flnc_fpath



