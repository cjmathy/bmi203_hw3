import numpy as np
import random

def align(pair,sequences,scoring_matrix,penalties,aa_dict):
	seq1, seq2 = sequences[pair[0]], sequences[pair[1]]

	#Create matrix F (containing maximum scores up to cell F(i,j)) and tracking matrix P (containing the tracking flag describing whether the step into the cell was a match or a gap, or if the cell is the end of a local alignment.
	F = np.zeros((len(seq1),len(seq2)))
	P = np.zeros((len(seq1),len(seq2)))

	F, P = fill_scoring_matrix(F,P,seq1,seq2,scoring_matrix,aa_dict,penalties)

	return traceback(F,P,seq1,seq2)

def fill_scoring_matrix(F,P,seq1,seq2,scoring_matrix,aa_dict,penalties):
	d,e = penalties

	it = np.nditer(F, flags=['multi_index'])
	while not it.finished:
		i,j = it.multi_index
		aa1_index = aa_dict[seq1[i]]
		aa2_index = aa_dict[seq2[j]]

		#Set penalties: opening a gap has penalty d, extending a gap has penalty e, and we do not allow a gap to be opened introduced in one sequence after a gap has occured in the other (because this would effectively allow for jumps from local alignment to local alignment).
		
		#If the cell directly above was arrived at by introducing a gap in seq1, then do not allow a gap in seq2
		if P[i-1,j] == 1:
			gap_in_seq1,gap_in_seq2 = e,float("inf")
		
		#If the cell directly to the left was arrived at by introducing a gap in seq2, then do not allow a gap in seq1
		elif P[i,j-1] == -1:		
			gap_in_seq1,gap_in_seq2 = float("inf"),e
		
		else:
			gap_in_seq1,gap_in_seq2 = d,d
		possible_scores	= [0,
							F[i-1,j-1]+scoring_matrix[aa1_index][aa2_index],
						  	F[i,j-1]-gap_in_seq1,
						  	F[i-1,j]-gap_in_seq2
						  ]

		#Find maximum of possible_scores to decide on move
		F[i,j] = max(possible_scores)
		if 	 F[i,j] == possible_scores[1]:	P[i,j] = 0	#If match, set tracking flag to 0
		elif F[i,j] == possible_scores[2]:	P[i,j] = 1	#If gap in seq1, set tracking flag to 1
		else							 :	P[i,j] = -1	#If gap in seq2 (or if max was 0, and alignment ended), set tracking flag to -1
		
		it.iternext()

	F[:,0] = 0
	F[0,:] = 0

	return F,P

def traceback(F,P,seq1,seq2):
	
	#Find max value
	max_index = F.argmax()
	i,j = np.unravel_index(max_index, F.shape)
	score = F[i,j]

	alignment = ('','')

	while F[i,j] > 0 and i is not 0 and j is not 0:

		if P[i,j] == 0: 	#There must have been a match of the amino acid at seq1[i] and seq2[j]
			alignment = seq1[i]+alignment[0],seq2[j]+alignment[1]
			i -= 1
			j -= 1
		elif P[i,j] == 1: 	#There must have been a gap in seq1
			alignment = "_"+alignment[0],seq2[j]+alignment[1]
			j -= 1
		elif P[i,j] == -1: #There must have been a gap in seq2
			alignment = seq1[i]+alignment[0],"_"+alignment[1]
			i -= 1

	return alignment, score

def get_fp_rate(pos_scores,neg_scores):
	#We want to set our threshold as the value of the element in the pos_scores list that is greater than or equal to exactly 30% of the scores in the list. Thus, we sort the pos_scores, and take the value of the 15th element out of the 50 in the list to be the threshold.
	pos_scores.sort()
	threshold = pos_scores[14]
	fp = 0
	for score in neg_scores:
		if score > threshold:
			fp+=1 
	fp_rate = fp/50.0
	return threshold, fp_rate

def get_tp_rate(pos_scores,neg_scores,fpr):
	
	neg_scores.sort()
	index = int(50*(1-fpr))-1
	threshold = neg_scores[index]
	tp = 0
	for score in pos_scores:
		if score > threshold:
			tp+=1
	tp_rate = tp/50.0

	return tp_rate

def calc_score(alignments,scoremat,penalties,aa_dict):
	d,e = penalties
	scores = []
	for alignment in alignments:
		score = 0
		for index in range(len(alignment[0])):
			score_flag = False
			aa1 = alignment[0][index]
			aa2 = alignment[1][index]
			if aa1 == '_':
				if ext_flag == True:
					score -= e
					score_flag = True
				else:
					score -= d
					ext_flag = True
					score_flag = True
			else: i = aa_dict[aa1]
			
			if aa2 == '_':
				if ext_flag == True:
					score -= e
					score_flag = True
				else:
					score -= d
					ext_flag = True
					score_flag = True
			else: j = aa_dict[aa2]

			if (score_flag == False):
				score += scoremat[i,j]
				ext_flag = False
		scores.append(score)
	return scores










