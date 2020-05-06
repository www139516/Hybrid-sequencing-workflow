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
from proc_modules.FLNC_generator.FlncGenerator import FlncGenerator


@click.command()
@click.option("--subreads",
              help="The path to the subreads.bam, the subreads.pbi and subreads.xml should be in the same directory")
@click.option("--output", help="The directory of the output files")
@click.option("--primer", help="The sequence of primer used in the PacBio sequencing")
def main(subreads, output, primer):

    flnc_gen = FlncGenerator()
    flnc_gen.fit(subreads, primer, output)
    flnc_gen.sub2ccs()
    flnc_gen.ccs2flccs()
    flnc_fpath = flnc_gen.flccs2flnc()


if __name__ == '__main__':
    main()
