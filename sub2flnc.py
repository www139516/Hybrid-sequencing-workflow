# coding: UTF-8
import __future__
import glob
import os
import subprocess
import shutil
import re
import click



def sub_flnc (subreads, output, primer_seq):
	curr_dir = os.getcwd()
	movie_name = os.path.basename(subreads).split(".")[0]
	if "/" in subreads:
		subreads = subreads
	else:
		subreads = curr_dir + "/" + subreads

	if "/" in output:
		output = output
	else:
		output = curr_dir + "/" + output

	if "/" in primer_seq:
		primer_seq = primer_seq
	else:
		primer_seq = curr_dir + "/" + primer_seq

	input_subreads = subreads
	output_ccs = output + "/" + movie_name + ".ccs.bam"
	
	cmd_subreads2ccs = 'ccs {input_subreads} {output} '.format(input_subreads=input_subreads, output=output_ccs) + '--noPolish --minPasses 1'
	if subprocess.check_call(cmd_subreads2ccs, shell=True) != 0:
		raise SystemCommandError

	output_flccs = output + "/" + movie_name + ".flccs.bam"
	cmd_ccs2flccs = 'lima {input_ccs} {barcode} {output_flccs} '.format(input_ccs=output_ccs, barcode=primer_seq, output_flccs=output_flccs) + '--isoseq --no-pbi --peek-guess'
	if subprocess.check_call(cmd_ccs2flccs, shell=True) != 0:
		raise SystemCommandError

	input_flnc = output + "/" + movie_name + ".flccs.5p_primer_5p--3p_primer_3p.bam"
	output_flnc = output + "/" + movie_name + ".flnc.bam"
	cmd_flccs2flnc = "isoseq3 refine {input_flccs} {barcode} {output_flnc} ".format(input_flccs=input_flnc, barcode=primer_seq, output_flnc=output_flnc) + "--require-polya"
	if subprocess.check_call(cmd_flccs2flnc, shell=True) != 0:
		raise SystemCommandError

@click.command()
@click.option("--subreads", help="The path to the subreads.bam, the subreads.pbi and subreads.xml should be in the same directory")
@click.option("--output", help="The directory of the output files")
@click.option("--primer", help="The sequence of primer used in the PacBio sequencing")
def main (subreads, output, primer):
	sub_flnc (subreads, output, primer)


if __name__ == '__main__':
	main()
