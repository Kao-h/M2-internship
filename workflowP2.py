import pandas as pd
import glob
import re
import numpy as np
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import os
from collections import Counter



RFam_ID_chainDF = pd.read_csv("/home/kaouther/RNA_Proj/coconetFinal/RFam_ID_chain.csv")
rfams = list(RFam_ID_chainDF['Rfam'])
pdbs = list(RFam_ID_chainDF['PDB'])
# create a dict key ids value chain in MSA
rfam_pdb_dict = {k: v for k, v in zip(rfams, pdbs)}

# reverse keys and values in dico

CocoRfams =  {v.lower(): k for k, v in rfam_pdb_dict.items()}


def make_vector(l1, l2) :
    """
    :param: l1, l2 2 lists
    :return: list of the fraction of total pairs of each possible di-nucleotides pair.

    This function converts 2 lists in a dict-nucleotides frequecies list
    """
    pair_dict = Counter({
        ('G', 'A'): 1,
        ('G', 'G'): 1,
        ('G', 'C'): 1,
        ('G', 'U'): 1,
        ('G', '-'): 1,

        ('C', 'C'): 1,
        ('C', 'U'): 1,
        ('C', 'G'): 1,
        ('C', 'A'): 1,
        ('C', '-'): 1,

        ('U', 'A'): 1,
        ('U', 'G'): 1,
        ('U', '-'): 1,
        ('U', 'U'): 1,
        ('U', 'C'): 1,

        ('A', '-'): 1,
        ('A', 'G'): 1,
        ('A', 'U'): 1,
        ('A', 'C'): 1,
        ('A', 'A'): 1,

        ('-', '-'): 1,
        ('-', 'C'): 1,
        ('-', 'U'): 1,
        ('-', 'G'): 1,
        ('-', 'A'): 1})

    d = Counter(list(zip(l1, l2)))
    a = d + pair_dict

    for key, value in a.items():
        pair_dict[key] = round(((value - 1) / (sum(a.values()) - 25)), 4)

    return list(pair_dict.values())


def read_file(filename):
    with open(filename, 'r') as rf:
        return [l.strip() for l in rf.readlines()]

def split(sequence):
    """
    :param: sequence
    :return: list of caracters
    this function splits sequence as chars
    """
    return [char for char in sequence]


def makeSeqList(filename):
    """
    :param: filename input MSA file in fasta format
    :return: a list of sequences like ['seq1','seq2', etc]

    this function takes a msa in fasta format and makes a list of sequences like ['seq1','seq2', etc]
    """
    FASTAFile = read_file(filename)
    FASTADict = {}
    FASTALabel = ''
    for line in FASTAFile:
        if '>' in line:
            FASTALabel = line
            FASTADict[FASTALabel] = ''
        else:
            FASTADict[FASTALabel] += line.replace('M', '-').replace('Y', '-').replace('N', '-').replace('R',
                                                                                                        '-').replace(
                'W', '-').replace('H', '-').replace('K', '-').replace('S', '-').replace('n', '-').replace('B',
                                                                                                          '-').replace(
                'D', '-').replace('V', '-')
    SeqListStr = []
    for i in FASTADict:
        SeqListStr.append(FASTADict[i])
    SeqList = []
    for sequence in SeqListStr:
        SeqList.append(split(sequence))

    return SeqList
# SeqList = makeSeqList(MSAfilePath)

def make_array(SeqList):
    '''
    takes a list of sequences and makes them into an array that represents a msa like [
                                                                                        [a,u,c,g],
                                                                                        [a,u,g,c],
                                                                                        [a,g,c,u],
                                                                                        ]
    '''

    rows = len(SeqList)
    ResiduesPerSequence = len(SeqList[0])
    SequenceMatrix = np.zeros([rows, ResiduesPerSequence], dtype=list)
    cc = int(ResiduesPerSequence * (ResiduesPerSequence))
    ResiduePairMatrix = np.zeros([cc, rows], dtype=list)
    MatrixCount = 0
    for i in SeqList:
        SequenceMatrix[MatrixCount] = i
        MatrixCount += 1

    return SequenceMatrix


def make_dinucleotide_freq_array(array):
    """
    :param: array: a representation of an MSA as a liste of list of nucleotides
    :return: an array with the di-nucleotide frequencies

    this function makes an lxlx25 array where each position i,j in the lxl array represents
    position i compared to position j in the rna sequence.
    Then at each position i,j it puts the dinucleotide frequency vector
    """
    num_col = int(np.shape(array)[1])
    cov_matrix = np.zeros([num_col, num_col, 25])
    for i in range(0, num_col):
        for j in range(0, num_col):
            if j > i:
                cov_matrix[i, j, :] = make_vector(array[:, i], array[:, j])
                cov_matrix[j, i, :] = make_vector(array[:, i], array[:, j])

    return cov_matrix


def make_dinucleotide_freq_df(path):
    """
    :param: path: folderpath to the MSAs
    :return: dataframe representing the lxlx25 dinucleotide frequencies matrix
    takes all the lxlx25 dinucleotide frequency arrays and puts them in a dataframe
    """
    li = []
    path = path
    all_files = list(glob.glob(path + "/RF*"))
    for filename in all_files:
        try:
            step = make_array(makeSeqList(filename))
            a = make_dinucleotide_freq_array(make_array(makeSeqList(filename)))
            # print(a)
            col_names = ['positions', 'covariance', 'file']
            df = pd.DataFrame(columns=col_names)

            for x in range(len(a[0])):
                for y in range(len(a[0])):
                    d = {'positions': str(x + 1) + ':' + str(y + 1), 'covariance': np.array(a[x, y, :]),
                             'file': rfam_pdb_dict[filename.split('/')[6].split('.')[0]]}

                    df = df.append(d, True)
            li.append(df)
            print(filename.split('/')[6].split('.')[0])

        except:
            print(filename, 'not working')
    frame = pd.concat(li, axis=0, ignore_index=True)
    return frame


df = make_dinucleotide_freq_df('/home/kaouther/RNA_Proj/coconetFinal/MSA')
#write df to excel
df.to_csv('/home/kaouther/RNA_Proj/coconetFinal/MSA/dinucleotidefreqDF.csv', index=False)
#floyds df
df = pd.read_csv("/home/kaouther/RNA_Proj/coconetFinal/labeleddinucfreq_march25.csv")
df = df.set_index('index')
df['NtCs'] = 0

df = df.drop('bond_type', axis=1)
# change the index of the df as a combo of the position and the pdb ID
df['index'] = df['positions']+df['file'].upper()
df = df.set_index('index')
# add a column for NtCs
df['NtCs'] = 0

#merge light DF containing data from DNATCO about NtCs

files = glob.glob(os.path.join("/home/kaouther/RNA_Proj/coconetFinal/lightDF/light*"))
# joining files with concat and read_csv
mergedDFs = pd.concat(map(pd.read_csv, files), ignore_index=True)
mergedDFs.to_csv("/home/kaouther/RNA_Proj/coconetFinal/lightDF/mergedDFs.csv", index=False)

#NtCs df
ef = pd.read_csv("/home/kaouther/RNA_Proj/coconetFinal/lightDF/mergedDFs.csv")

# create some format index as df
indexlist= []
for i in range(len(ef)):
    index = ef['SeqPos-SeqPos'][i].split('-')[0]+':'+ef['SeqPos-SeqPos'][i].split('-')[1]+ ef['PDB_ID'][i].upper()
    indexlist.append(index)
ef['index'] = indexlist
ef = ef.set_index('index')

#hot encoding of the NtCs for better readability
l = list(set(ef['NtCs']))
l2 = [item for item in range(1, len(l)+1)]
encoding_dict = {k: v for k, v in zip(l, l2)}
ef.replace({'NtCs': encoding_dict}, inplace = True)

#save encoding dict in a file
with open('/home/kaouther/RNA_Proj/coconetFinal/NtCs_encoding_dict.csv', 'w') as f:
    for key in encoding_dict.keys():
        f.write("%s,%s\n"%(key,encoding_dict[key]))

#add encoded column to df
for index, row in ef.iterrows():
    if index in df.index:

        df.loc[index,'NtCs'] = int(row[6])

df = df.drop(df[df.index =='16:446HAG'].index)

# export to csv
df.to_csv("/home/kaouther/RNA_Proj/coconetFinal/LabeledNtCsDataset5April.csv")


#counter for all possible di-nucleotide pairs
pair_dict = Counter({
        ('G', 'A'): 1,
        ('G', 'G'): 1,
        ('G', 'C'): 1,
        ('G', 'U'): 1,
        ('G', '-'): 1,

        ('C', 'C'): 1,
        ('C', 'U'): 1,
        ('C', 'G'): 1,
        ('C', 'A'): 1,
        ('C', '-'): 1,

        ('U', 'A'): 1,
        ('U', 'G'): 1,
        ('U', '-'): 1,
        ('U', 'U'): 1,
        ('U', 'C'): 1,

        ('A', '-'): 1,
        ('A', 'G'): 1,
        ('A', 'U'): 1,
        ('A', 'C'): 1,
        ('A', 'A'): 1,

        ('-', '-'): 1,
        ('-', 'C'): 1,
        ('-', 'U'): 1,
        ('-', 'G'): 1,
        ('-', 'A'): 1})

#create labels for covariance vector
listLabels = list(pair_dict.keys())
listLabelsDF = []

for label in listLabels:
    colLabel= ''.join(label)
    listLabelsDF.append(colLabel)

covarianceColumn = list(df['covariance'])

covarianceDF = pd.DataFrame(columns=['covariance','list'])
covarianceDF['covariance'] =  covarianceColumn
for i in range(len(covarianceDF)):
    covarianceDF.iloc[i]['list'] = covarianceDF.iloc[i]['covariance'].split()

def vectorToFloatlist(vector):
    """
    :param vector: input a vector of strings
    :return: a list of floats representing the covariance
    this function converts a list of string to a list of floats
    """

    for i in range(len(vector)):
        if "[" in vector[i]:
            new = vector[i].replace("[", "")
            vector[i] = new
        if "]" in vector[i]:
            new = vector[i].replace("]", "")
            vector[i] = new
        if "\n" in vector[i]:
            new = vector[i].replace("\n", "")
    if '' in vector:
        vector.remove('')

    for i in range(len(vector)):
        floate = float(vector[i])
        vector[i] = floate

    return vector

vectorListcolumn = list(covarianceDF['list'])
newVectorList = []
for vector in vectorListcolumn:

    newVector = vectorToFloatlist(vector)
    newVectorList.append(newVector)

#replace covariance with list of corresponding floats
df['covariance'] = newVectorList
# dict keys co variance floats values list of di-nuc freq
dico = {}
for col in listLabelsDF:
    dico[col] = []
for i in range(len(df)):
    labeldict =  dict(zip(listLabelsDF, df.iloc[i]["covariance"]))
    for col in labeldict.keys():
        dico[col].append(labeldict[col])
# add 25 new columns to df
for col in listLabelsDF:
    df[col] = list(dico[col])
df = df.drop('covariance',axis=1)
df.to_csv("/home/kaouther/RNA_Proj/coconetFinal/TESTLabeled_dataset.csv", index=False)