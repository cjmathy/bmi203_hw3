import io
import methods

import numpy as np
import os
import itertools
import random
import time

#Read in sequences to align, storing them as a dictionary
#(e.g.: sequences['sequences/prot-0014.fa'] = 'DACEQAAIQCVESACESLCTEGEDRTGCYMYIYSNCPPYV')
sequences = io.read_sequences(os.getcwd())

#Read in positive and negative pairs as dictionaries.
pospairs = io.read_pairs(os.getcwd()+'/Pospairs.txt')
negpairs = io.read_pairs(os.getcwd()+'/Negpairs.txt')

#Iterate through all score matrices
scoremat_files = ['BLOSUM50','BLOSUM62','PAM100','PAM250','MATIO']
scoremat_file = scoremat_files[4]
scoremat, aa_dict, aa_list = io.read_scoring_matrix(scoremat_file)

#Define d = penalty to open a gap, e = penalty to extend a gap. Iterate through various d,e values.
penalties = d,e = 20,3

# Align all pairs of positive sequences
# with open(output, 'a') as f: f.write("~~Positive Pairs\n")
pos_scores = []	
pos_alignments = []
for pair in pospairs:
	alignment, score = methods.align(pair,sequences,scoremat,penalties,aa_dict)
	pos_scores.append(score)
	pos_alignments.append(alignment)
	# io.write_alignment(pair,alignment,score,output,pospairs,negpairs)

# Align all pairs of negative sequences
# with open(output, 'a') as f: f.write("~~Negative Pairs\n")
neg_scores = []	
neg_alignments = []
for pair in negpairs:
	alignment, score = methods.align(pair,sequences,scoremat,penalties,aa_dict)
	neg_scores.append(score)
	neg_alignments.append(alignment)
	# io.write_alignment(pair,alignment,score,output,pospairs,negpairs)

# threshold, fp_rate = methods.get_fp_rate(pos_scores,neg_scores)

output = io.prepare_output('MATIO_OPT')

best_obj_fn = 0
best_matrix = np.copy(scoremat)

kmax = 50000

for k in range(1,kmax):
	if k%10 == 0:
		print "k=%i" %k
	pos_scores = methods.calc_score(pos_alignments,scoremat,penalties,aa_dict)
	neg_scores = methods.calc_score(neg_alignments,scoremat,penalties,aa_dict)
	obj_fn = 0
	for fpr in [0,0.1,0.2,0.3]:
		obj_fn += methods.get_tp_rate(pos_scores,neg_scores,fpr)
	delta = obj_fn - best_obj_fn
	if delta >= 0:
		best_obj_fn = obj_fn
		print best_obj_fn
		best_matrix = np.copy(scoremat)
		curr = np.copy(scoremat)
	elif np.exp(delta/(float(k)/kmax)) > random.random():
		curr = np.copy(scoremat)
	else: scoremat = np.copy(curr)
	
	i = random.randint(0,scoremat.shape[0]-1)
	j = random.randint(0,scoremat.shape[1]-1)
	scoremat[i,j] += 3*np.random.normal()
	scoremat[j,i] = scoremat[i,j]


print best_obj_fn
# with open(output, 'a') as f: f.write("\nthreshold = %f, false positive rate = %f" %(threshold,fp_rate))

with open(output, 'w') as f:
	f.write(' ')
	for aa in aa_list:
		f.write("  %s" %aa)
	for i in range(best_matrix.shape[0]):
		f.write('\n')
		for j in range(best_matrix.shape[1]):
			f.write(' %f' %best_matrix[i][j])
	
io.make_roc_curve(pos_scores,neg_scores,'MATIO_OPT')
# io.make_roc_curve(pos_scores,neg_scores,scoremat_file)


