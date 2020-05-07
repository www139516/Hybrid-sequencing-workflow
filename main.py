#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time : 2020/5/6 11:07
# @Author : Yuancong Wang
# @Site : 
# @File : main.py.py
# @Software: PyCharm
# E-mail: wangyuancong@163.com
'''
Description: The one-step pipeline to get the high-conficence full length cDNA sequences
Output: The validated FLNC CCSs
Input: the "subreads.bam", the "primer.fq", the output dir (the current dir if not give)
Other notes: None
'''
import click
import subprocess
from proc_modules.FLNC_generator.FlncGenerator import FlncGenerator
from proc_modules.correct_FLNC.FLNC_corrector import FlncCorrect


@click.command()
@click.option("--subreads",
              help="The path to the subreads.bam, the subreads.pbi and subreads.xml should be in the same directory")
@click.option("--output", help="The directory of the output files")
@click.option("--primer", help="The sequence of primer used in the PacBio sequencing")
@click.option('--proov', help='The path of the proovread exec file, usually under the bin/. Do not need give \
                              if the proovread is already added to the $PATH', default='')
@click.option('--sread', help='The path of the short-read sequencing file (RNA-seq), seperated by one space', nargs=-1)
def main(subreads, output, primer, proov, sread):

    flnc_gen = FlncGenerator()
    flnc_gen.fit(subreads, primer, output)
    flnc_gen.sub2ccs()
    flnc_gen.ccs2flccs()
    flnc_fpath = flnc_gen.flccs2flnc()

    flnc_corrector = FlncCorrect()
    flnc_corrector = flnc_corrector.fit(flnc_fpath, sread, proov)
    pread_out_dpath = flnc_corrector.correct_flnc()


if __name__ == '__main__':
    main()
