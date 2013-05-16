#!/usr/bin/env python 

import sys 
import json 

def main(file_in, file_out, biom):
	print 'Converting ' + file_in + ' to ' + file_out
	fin_handle = open(file_in, 'U')
	fout_handle = open(file_out, 'w')
	fb_handle = open(biom, 'U')

	o = json.loads(fb_handle.read())
	names = [u['id'] for u in o['columns']]
	labels = fin_handle.read().split()

	fout_handle.write('#SampleID\tClass\n')
	for n,l in map(None, names, labels):
		if float(l) == 1.0:
			fout_handle.write(n+'\tC1'+'\n')
		else:
			fout_handle.write(n+'\tC2'+'\n')
	
	fb_handle.close()
	fin_handle.close()
	fout_handle.close()
	return None


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3])