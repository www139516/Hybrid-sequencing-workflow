import glob
import os
import subprocess
import shutil



def accur(Read):
	if 'NM:i:' in Read:
		Len=len(Read.split('\t')[9])
		MisMatch=int(re.match(".+NM:i:(\d+)", Read).group(1))
		Accur=MisMatch/Len
	return Accur

def sam_filter_sort(ori_sam, processed_sam):
	filtering_sam=open('filtering_sam.sam','w')
	for line in open(ori_sam, 'r'):
		if line[0]=='@':
			filtering_sam.write('{0}'.format(line))
		else:
			if int(line.split('\t')[4]) >= 30 and accur(line) < 0.15:
				filtering_sam.write('{0}'.format(line))
	filtering_sam.close()
	cmd='sort -k 3,3 -k 4,4n {0} > {1}'.format('filtering_sam.sam', processed_sam)
	if subprocess.check_call(cmd, shell=True) != 0:
		raise SystemCommandError
	

def main(args):
	sorted_sam = sam_filter_sort(args.input, args.output)



if __name__ == "__main__":
	from argparse import ArgumentParser
	parser = ArgumentParser(description='collapsing the error-corrected FLNC CCSs into splicing isoforms and validated them by checking their chained splicing junctions with high confidence short reads')
	parser.add_argument('-input', required=True, help='Unsorted sam file')
	parser.add_argument('-output', required=True, help='The output file of sorted sam file')
	args = parser.parse_args()
	main(args)
