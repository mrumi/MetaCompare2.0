from sys import argv
from sys import exit
import os

from annotation import generate_annotation, process_annotation
from calculation import calculate_score
from visualization import visual_file_creation

def getopts(argv):
	opts = {}  # Empty dictionary to store key-value pairs.
	while argv:  # While there are arguments left to parse...
		if argv[0][0] == '-' and argv[0][1] != 'h':  # Found a "-name value" pair.
			opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
		argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
	return opts

if __name__ == '__main__':
	myargs = getopts(argv)
	if '-h' in myargs or len(myargs) == 0:  # Example usage.
		print('\nUsage: ./metacmp.py -c filename1.fa [-t 64 -b 1] \n')
		print('\t-c: Specify FASTA file containing assembled contigs')
		print('\t-t: Specify the number of threads will be used in executing blast (default: 64).')
		print('\t-b: Specify the pipeline to execute [0: both (default), 1: ecological risk score, 2: human health risk score ].')
		print('\t-o: output file path')
		print()
		exit()


	if '-t' in myargs:
		nthread = myargs['-t']
	else:
		nthread = '64'

	if not '-o' in myargs:
		myargs['-o'] = ''
		
	if not '-b' in myargs:
		myargs['-b'] = '0'
		
	pipeline = int(myargs['-b'])
	mge_len_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metacmpDB/mobileAll_len.txt")
	sample_name = os.path.splitext(os.path.basename(myargs['-c']))[0]
	out_file = os.path.join(myargs['-o'], sample_name + "_out.txt")	
		
	annotated_data = generate_annotation(myargs['-c'], myargs['-o'], nthread)	
	
	if pipeline == 1:		
		data_to_be_processed = [annotated_data[0], annotated_data[1], annotated_data[2]]
		pathogens = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metacmpDB/pathogen_list.txt")
		filtered_data = process_annotation(data_to_be_processed, mge_len_file, pathogens)
		result = calculate_score(myargs['-c'], filtered_data, pipeline)
	elif pipeline == 2:		
		data_to_be_processed = [annotated_data[3], annotated_data[1], annotated_data[2]]
		pathogens = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metacmpDB/eskape.txt")
		filtered_data = process_annotation(data_to_be_processed, mge_len_file, pathogens)
		result = calculate_score(myargs['-c'], filtered_data, pipeline)
	else:
		data_to_be_processed = [annotated_data[0], annotated_data[1], annotated_data[2]]
		pathogens = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metacmpDB/pathogen_list.txt")
		filtered_data_e = process_annotation(data_to_be_processed, mge_len_file, pathogens)
		result_e = calculate_score(myargs['-c'], filtered_data_e, 1)
		
		data_to_be_processed = [annotated_data[3], annotated_data[1], annotated_data[2]]
		pathogens = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metacmpDB/eskape.txt")
		filtered_data_h = process_annotation(data_to_be_processed, mge_len_file, pathogens)
		result_h = calculate_score(myargs['-c'], filtered_data_h, 2)
		
		result = result_e.append(result_h)
	
	result.to_csv(out_file, header=True, index=None, sep = "\t")
			
	
	
