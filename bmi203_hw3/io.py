import glob
import numpy as np
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

def read_scoring_matrix(filename):
    with open(filename,'r') as f:
        line = f.readline()
        while(line.startswith("#")): line = f.readline()

        #This will be the row containing the amino acid labels. Make a dictionary mapping the amino acid letter to the index of the scoring matrix
        aa_labels = list(line.split())
        aa_dict = {}
        for aa in aa_labels: aa_dict[aa] = aa_labels.index(aa)

        #Convert every row from a string into a list of integers, then add it to a numpy array to create the scoring matrix
        matrix = np.zeros((len(aa_labels),len(aa_labels)))
        for i in range(len(aa_labels)): matrix[i,:] = list(map(int, f.readline().split()))
    return matrix, aa_dict, aa_labels

def read_sequences(cwd):
    filenames = glob.glob(cwd + '/sequences/*.fa')
    sequences = {}
    for filename in filenames:
        name = filename[-22:]   #e.g.: name = 'sequences/prot-0014.fa'
        with open(filename, 'r') as f:
            seq = ''
            f.readline() #skips header row
            for line in f.readlines(): seq += line.strip()
            seq = seq.upper()
        sequences[name] = seq
    return sequences

def read_pairs(filename):
    pairs = []
    with open(filename, 'r') as f:
        for line in f: pairs.append(line.split())
    return pairs

def prepare_output(filename):
    with open(filename, 'w') as f: f.write("")
    return filename

def write_alignment(pair,alignment,score,filename,pospairs,negpairs):
    with open(filename, 'a') as f:
        f.write("#New Alignment\n")
        f.write("Sequence 1: %s\n" %pair[0])
        f.write("Sequence 2: %s\n" %pair[1])
        f.write("%s\n" %alignment[0])
        f.write("%s\n" %alignment[1])
        f.write("Score = %f\n" %score)
    return

def make_roc_curve(pos_scores,neg_scores,matrix):
    #Create a label array, y, with 1's for positives and 0's for negatives, and a scoring array with the scores for each labeled pair.
    y = np.array([1]*len(pos_scores)+[0]*len(neg_scores))
    scores = np.array(pos_scores + neg_scores)

    #Using scikit-learn's roc_curve and auc, calculate the parameters necessary for an ROC curve, and use them to compute the area under the curve
    fpr,tpr,thresholds = roc_curve(y,scores)
    roc_auc = auc(fpr, tpr)

    #plot and save an ROC curve.
    plt.figure()
    plt.title('Receiver Operating Characteristic, %s' %matrix)
    plt.plot(fpr,tpr, 'b', label='AUC = %0.2f'% roc_auc)
    plt.legend(loc='lower right')
    plt.plot([0,1],[0,1],'r--')
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.savefig('roc_%s.png' %matrix)

    return
