#!/usr/bin/python

import os
from Bio import SeqIO

def get_best_hit(data):
	# get the rows with highest bit-score for each id
	idx = data.groupby(['id'])['bit'].transform(max) == data['bit']
	best_hit_id = data[idx]

	# get the rows with lowest e-value for each id
	idx = best_hit_id.groupby(['id'])['eval'].transform(min) == best_hit_id['eval']
	best_hit_eval = best_hit_id[idx]

	# get the first result after scrutinizing by bit-score and e-value
	best_hit = best_hit_eval.drop_duplicates(['id'])
				
	# add scaffold id
	scaffold_id = best_hit['id'].apply(lambda x: "_".join(x.split("_")[:-1]))
	best_hit['scaffold_id'] = scaffold_id
	
	return best_hit
	
def visual_file_creation(all_data, out_path, prodigal_output, pipeline):
	
	def generate_visualization(arg_data, mge_data, pat_data, names):
		if not (os.path.exists(os.path.join(out_path, names[0])) or os.path.exists(os.path.join(out_path, names[1])) or os.path.exists(os.path.join(out_path, names[2]))):
			try:
				print("Creating Visualization files for web interface")         

				# Get best hits for ARGs
				arg_best_hit = get_best_hit(arg_data)
				
				# Get Best hits for MGEs				
				mge_best_hit = get_best_hit(mge_data)
				
				# Get best hits for PATHOGENs
				# get the rows with highest bit-score for each id
				idx = pat_data.groupby(['id'])['bit'].transform(max) == pat_data['bit']
				pat_best_hit = pat_data[idx]

				# get the first result after scrutinizing by bit-score and e-value
				pat_best_hit = pat_best_hit.drop_duplicates(['id'])
				pat_best_hit = pat_best_hit.rename(columns = {'id':'scaffold_id'})
				
				#############################################
				# Changes using prodigal output

				# 1. Fix qStart and qEnd in gene hits
				# 2. Add prodigal gene start and end
				# 3. Add prodigal gene length
				#############################################
				scaffold_indexed = SeqIO.to_dict(SeqIO.parse(prodigal_output, "fasta"))
				arg_scaffolds = [scaffold_indexed.get(key) for key in arg_best_hit.id.values]
				mge_scaffolds = [scaffold_indexed.get(key) for key in mge_best_hit.id.values]				

				for item in arg_scaffolds:
					start = int(item.description.split('#')[1])
					end = item.description.split('#')[2]
					length = int(len(item))

					try:
						row_idx = arg_best_hit.loc[arg_best_hit['id'] == item.id].index[0].item()
						arg_best_hit.loc[row_idx, 'qStart'] = arg_best_hit.loc[row_idx, 'qStart'].item() + start - 1
						arg_best_hit.loc[row_idx, 'qEnd'] = arg_best_hit.loc[row_idx, 'qEnd'].item() + start - 1
						arg_best_hit.loc[row_idx, 'pg_length'] = length
						arg_best_hit.loc[row_idx, 'pg_start'] = start
						arg_best_hit.loc[row_idx, 'pg_end'] = end
																			
					except:
						None
				
				for item in mge_scaffolds:
					start = int(item.description.split('#')[1])
					end = item.description.split('#')[2]
					length = int(len(item))
					try:
						row_idx = mge_best_hit.loc[mge_best_hit['id'] == item.id].index[0].item()
						mge_best_hit.loc[row_idx, 'qStart'] = mge_best_hit.loc[row_idx, 'qStart'].item() + start - 1
						mge_best_hit.loc[row_idx, 'qEnd'] = mge_best_hit.loc[row_idx, 'qEnd'].item() + start - 1
						mge_best_hit.loc[row_idx, 'pg_length'] = length
						mge_best_hit.loc[row_idx, 'pg_start'] = start
						mge_best_hit.loc[row_idx, 'pg_end'] = end
																			
					except:
						None
				
				#########################################
				# End of prodigal change
				#########################################
				# Export dataframe results to csv file
				arg_best_hit.to_csv(os.path.join(out_path, names[0]))
				mge_best_hit.to_csv(os.path.join(out_path, names[1]))
				pat_best_hit.to_csv(os.path.join(out_path, names[2]))
				print("Visualization files created successfully!")
			except:
				print("Visualization files could not be created!")

		else:
			print("Visualization files already exists.")
	if pipeline == 1:		
		names = ["arg_result.csv", "mge_result.csv", "pat_result.csv"]
	elif pipeline == 2:
		names = ["arg_hh_result.csv", "mge_hh_result.csv", "pat_hh_result.csv"]
	generate_visualization(all_data[0], all_data[1], all_data[2], names)
			
