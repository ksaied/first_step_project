# -*- coding: utf-8 -*-
"""
Created on Sat Dec  2 11:21:26 2017

@author: Karim
"""
###############################################################################
#Processing data
###############################################################################

import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np

font = {'family' : 'monospace',
        'weight' : 'normal',
        'size'   : 22}

plt.rc('font', **font)  # pass in the font dict as kwargs
hfont = {'fontname':'Helvetica'}

#The current tables contain all from human that have secondary structure
best_profiles = pd.read_csv('/Users/karimsaied/Desktop/Biology_v2/master/first_step_project/data/csv_tables_human_v2/best_profiles_human_v2.csv')
genes = pd.read_csv('/Users/karimsaied/Desktop/Biology_v2/master/first_step_project/data/csv_tables_human_v2/genes_human_v2.csv')
motifs = pd.read_csv('/Users/karimsaied/Desktop/Biology_v2/master/first_step_project/data/csv_tables_human_v2/motifs_human_v2.csv')
ensembl = pd.read_csv('/Users/karimsaied/Desktop/Biology_v2/master/first_step_project/data/csv_tables_human_v2/ensembl_human_v2.csv')
conserv = pd.read_csv('/Users/karimsaied/Desktop/Biology_v2/master/first_step_project/data/csv_tables_human_v2/gene_conserv_human_v2.csv')
site_conserv = pd.read_csv('/Users/karimsaied/Desktop/Biology_v2/master/first_step_project/data/csv_tables_human_v2/site_conserv_human_v2.csv')

#Useless columns are removed
best_profiles = best_profiles.drop([
        'ID',
        'BEST_PROFILE', 
        'ESTIMATED_W1_NULL_MODEL', 
        'ESTIMATED_W2_NULL_MODEL', 
        'ESTIMATED_LOG_LIKELIHOOD_NULL_MODEL',
        'S_COEV',
        'D_COEV',
        'R1_COEV',
        'R2_COEV',
        'ESTIMATED_LOG_LIKELIHOOD_COEV'
        ], axis=1)

ensembl = ensembl.drop([
        'ALIGNEMENT',
        'TREE'
        ], axis=1)

#Creation of a list of values for each data frame -> matrices
genes_val = genes.values.tolist()
best_profiles_val = best_profiles.values.tolist()
ensembl_val = ensembl.values.tolist()
motifs_val = motifs.values.tolist()
conserv_val = conserv.values.tolist()
site_conserv_val = site_conserv.values.tolist()

#In genes_val, replacement of the first element of each row by the id corresponding to ID_ENS in best_profiles.
for i in genes_val:
    for j in ensembl_val:
        if i[0] == j[1]:
            i[0] = j[0]

#In motifs_val, replacement of the first element of each row by the id corresponding to ID_ENS in best_profiles.
for i in motifs_val:
    for j in ensembl_val:
        if i[0] == j[1]:
            i[0] = j[0]

#In conserv_val, replacement of the first element of each row by the id corresponding to ID_ENS in best_profiles.
for i in conserv_val:
    for j in ensembl_val:
        if i[0] == j[1]:
            i[0] = j[0]

#In site_conserv_val, replacement of the first element of each row by the id corresponding to ID_ENS in best_profiles.
for i in site_conserv_val:
    for j in ensembl_val:
        if i[0] == j[1]:
            i[0] = j[0]
            
#In proc_val, replacement of the first element of each row by the id corresponding to ID_ENS in best_profiles.
for i in proc_val:
    for j in ensembl_val:
        if i[0] == j[1]:
            i[0] = j[0]
            
#Selection of genes that have a secondary structure --> equal to genes_val
'''motifs_yes = []
for line in genes_val:
    if line[3] == 'Yes':
        motifs_yes.append(line)
'''
#Converting float in integer from best_profiles_val
for row in best_profiles_val:
    row[0], row[1], row[2] = int(row[0]), int(row[1]), int(row[2])

#Positions in best_profiles (best_profiles_val) have to be converted in term of amino acids positions.
def codon(nuc_pos): #nucleotide position as argument
    codon_pos = math.ceil(nuc_pos/3) # It returns the smallest integer greater than or equal to x
    return codon_pos

for i in best_profiles_val:
    i[1], i[2] = codon(i[1]), codon(i[2]) # Due to the conversion from nucleotide to amino acids, two pairs can be present more than once.

###############################################################################
#First plot (histogram)
###############################################################################
#get_array(2603, motifs_val)    

def get_array(ens_id, motifs_val):
    aux = []
    array = []
    for row in motifs_val:
        if ens_id == row[0]:
            aux.append(row)
        
    array = ['' for row in range(10000)]   
    
    for row in aux:
        if row[2] == 'HELIX':
            for j in range(row[3], row[4]+1):
                array[j-1] = 'H'
        elif row[2] == 'STRAND':
            for j in range(row[3], row[4]+1):
                array[j-1] = 'S'
        elif row[2] == 'TURN':
            for j in range(row[3], row[4]+1):
                array[j-1] = 'T'
                
    return array

##########################################################
    
def get_motifs(pos1, pos2, struct):
    motif_1 = struct[pos1-1] 
    motif_2 = struct[pos2-1] 
   
    sol = ''
    
    if motif_1 == 'H' and motif_2 == 'H':
        sol = 'H-H'
    elif motif_1 == 'H' and motif_2 == 'S':
        sol = 'H-S'
    elif motif_1 == 'S' and motif_2 == 'H':
        sol = 'H-S'
    elif motif_1 == 'T' and motif_2 == 'T':
        sol = 'T-T'
    elif motif_1 == 'T' and motif_2 == 'H':
        sol = 'H-T'
    elif motif_1 == 'H' and motif_2 == 'T':
        sol = 'H-T'
    elif motif_1 == 'S' and motif_2 == 'S':
        sol = 'S-S'
    elif motif_1 == 'S' and motif_2 == 'T':
        sol = 'S-T'
    elif motif_1 == 'T' and motif_2 == 'S':
        sol = 'S-T'
    elif motif_1 == 'T' and motif_2 == 'T':
        sol = 'T-T'
    '''else: 
        if motif_1 == 'H' and motif_2 == '':
            sol = 'H-X'
        elif motif_1 == 'S' and motif_2 == '':
            sol = 'S-X'
        elif motif_1 == 'T' and motif_2 == '':
            sol = 'T-X'
        elif motif_1 == '' and motif_2 == 'H':
            sol = 'X-H'
        elif motif_1 == '' and motif_2 == 'S':
            sol = 'X-S'
        elif motif_1 == '' and motif_2 == 'T':
            sol = 'X-T'
        elif motif_1 == '' and motif_2 == '':
            sol = 'X-X'
    '''
    return sol

##########################################################

list_motifs = []

for i in best_profiles_val:
    try:
        struct = get_array(i[0], motifs_val)
    except: 
        print('Error in get_array(i[0], motifs_val)')
        print (i)
    try:
        list_motifs.append(get_motifs(i[1], i[2], struct))
    except:
        print('Error in get_motifs(i[1], i[2], struct)')
        print(i)

motifs_patterns = []
for i in list_motifs:
    if i == 'H-H':
        motifs_patterns.append(i)
    elif i == 'H-S':
        motifs_patterns.append(i)
    elif i == 'H-T':
        motifs_patterns.append(i)
    elif i == 'S-S':
        motifs_patterns.append(i)
    elif i == 'S-T':
        motifs_patterns.append(i)
    elif i == 'T-T':
        motifs_patterns.append(i)

struct_patterns = pd.DataFrame(motifs_patterns)
struct_patterns.to_csv('motifs_patterns.csv', sep=',')
motifs_patterns.count('S-S') # To count the occurrences of 'H-H'

plt.hist(motifs_patterns, bins = 20, color='#808080')
plt.xlabel('Secondary structure patterns', **hfont)
plt.ylabel('Coevolving positions pairs', **hfont)
plt.title('Distribution of positions pairs in the secondary structures', **hfont)
plt.show()

###############################################################################
#Boxplot - DAIC against motifs pattern
###############################################################################
# Instead, show DAIC of positions that are in motifs H-H, H-S, etc.
'''struct_yes = []
struct_no = []
motifs_ens_id = [i[0] for i in motifs_yes] # Storage of the ens_ids that have a secondary structure
for i in best_profiles_val:
    if i[0] in motifs_ens_id:
        struct_yes.append(i[3])
    else:
        struct_no.append(i[3])

data_motifs = [struct_yes, struct_no]
'''

DAIC_HH = []
DAIC_HS = []
DAIC_HT = []
DAIC_HX = []
DAIC_SS = []
DAIC_ST = []
DAIC_SX = []
DAIC_TT = []
DAIC_TX = []
DAIC_XX = []
pre_data = []
for i in best_profiles_val:
    buffer = ''
    struct = get_array(i[0], motifs_val)
    which_motifs = get_motifs(i[1], i[2], struct)
    pre_data.append([i[0], i[1], i[2], which_motifs, i[3]])

for i in pre_data:
    if i[3] == 'H-H':
        DAIC_HH.append(i[4])
    if i[3] == 'H-S':
        DAIC_HS.append(i[4])
    if i[3] == 'H-T':
        DAIC_HT.append(i[4])
    if i[3] == 'H-X' or i[3] == 'X-H':
        DAIC_HX.append(i[4])
    if i[3] == 'S-S':
        DAIC_SS.append(i[4])
    if i[3] == 'S-T':
        DAIC_ST.append(i[4])
    if i[3] == 'S-X' or i[3] == 'X-S':
        DAIC_SX.append(i[4])
    if i[3] == 'T-T':
        DAIC_TT.append(i[4])
    if i[3] == 'T-X' or i[3] == 'X-T':
        DAIC_TX.append(i[4])
    if i[3] == 'X-X':
        DAIC_XX.append(i[4])
        
data_motifs = [
        DAIC_HH,
        DAIC_HS,
        DAIC_HT,
        DAIC_HX,
        DAIC_SS,
        DAIC_ST,
        DAIC_SX,
        DAIC_TT,
        DAIC_TX,
        DAIC_XX
        ]

# Create a figure instance
fig = plt.figure(1)
fig.suptitle('Distribution of DAIC in different secondary structure patterns', fontsize=22, fontweight='bold', **hfont)

# Create an axes instance
ax = fig.add_subplot(111)
# "111" means "1x1 grid, first subplot" and "234" means "2x3 grid, 4th subplot".
# ax.set_aspect(0.05) # or some other float

# Create the boxplot
ax.boxplot(data_motifs, widths=0.6)

# Add x labels
ax.set_xticklabels(['H-H', 'H-S', 'H-T', 'H-X', 'S-S', 'S-T', 'S-X', 'T-T', 'T-X', 'X-X'], **hfont)
ax.set_ylabel('DAIC', **hfont)
ax.set_ylim([0,150])

###############################################################################
#DAIC mean per motifs
###############################################################################

DAIC_HH = []
DAIC_HS = []
DAIC_HT = []
DAIC_HX = []
DAIC_SS = []
DAIC_ST = []
DAIC_SX = []
DAIC_TT = []
DAIC_TX = []
DAIC_XX = []
pre_data = []
for i in best_profiles_val:
    buffer = ''
    struct = get_array(i[0], motifs_val)
    which_motifs = get_motifs(i[1], i[2], struct)
    pre_data.append([i[0], i[1], i[2], which_motifs, i[3]])

for i in pre_data:
    if i[3] == 'H-H':
        DAIC_HH.append(i[4])
    if i[3] == 'H-S':
        DAIC_HS.append(i[4])
    if i[3] == 'H-T':
        DAIC_HT.append(i[4])
    if i[3] == 'H-X' or i[3] == 'X-H':
        DAIC_HX.append(i[4])
    if i[3] == 'S-S':
        DAIC_SS.append(i[4])
    if i[3] == 'S-T':
        DAIC_ST.append(i[4])
    if i[3] == 'S-X' or i[3] == 'X-S':
        DAIC_SX.append(i[4])
    if i[3] == 'T-T':
        DAIC_TT.append(i[4])
    if i[3] == 'T-X' or i[3] == 'X-T':
        DAIC_TX.append(i[4])
    if i[3] == 'X-X':
        DAIC_XX.append(i[4])

# Appending the DAIC means for each motifs pattern
DAIC_means_per_motif = []
DAIC_means_per_motif.append(sum(DAIC_HH)/len(DAIC_HH)) 
DAIC_means_per_motif.append(sum(DAIC_HS)/len(DAIC_HS))
DAIC_means_per_motif.append(sum(DAIC_HT)/len(DAIC_HT))
DAIC_means_per_motif.append(sum(DAIC_HX)/len(DAIC_HX))
DAIC_means_per_motif.append(sum(DAIC_SS)/len(DAIC_SS))
DAIC_means_per_motif.append(sum(DAIC_ST)/len(DAIC_ST))
DAIC_means_per_motif.append(sum(DAIC_SX)/len(DAIC_SX))
DAIC_means_per_motif.append(sum(DAIC_TT)/len(DAIC_TT))
DAIC_means_per_motif.append(sum(DAIC_TX)/len(DAIC_TX))
DAIC_means_per_motif.append(sum(DAIC_XX)/len(DAIC_XX))

labels = ['H-H', 'H-S', 'H-T', 'H-X', 'S-S', 'S-T', 'S-X', 'T-T', 'T-X', 'X-X']
indices2 = np.arange(10)
plt.bar(indices2, DAIC_means_per_motif, color='#808080')
plt.xticks(indices2, labels)
plt.ylabel('DAIC mean', **hfont)
plt.title('DAIC mean per motif', **hfont)
plt.ylim([0, 35])
plt.show()

daic = [i[4] for i in pre_data]
import scipy.stats as stats
norm2 = stats.boxcox(daic)
stats.normaltest(norm2)

###############################################################################
#DAIC against conservation
###############################################################################

daic = []
list_conserv = []

for i in best_profiles_val:
    for j in conserv_val:
        if i[0] == j[0]:
            if i[3] < 150: # Threshold is set because some DAIC value are huge and it flattens the graph
                daic.append(i[3])
                list_conserv.append(j[1])


plt.scatter(list_conserv, daic)

###############################################################################
#Listing all biological processes
###############################################################################
'''
distinct_processes = []

for i in proc_val:
    if i[2] not in distinct_processes:
        distinct_processes.append(i[2])
'''
###############################################################################
#Displaying the genes that have DAIC > 100
###############################################################################
'''
high_daic = []
for i in best_profiles_val:
    if i[3] > 100:
        high_daic.append(i)

df = pd.DataFrame(high_daic)
df.columns =  ['ensembl_id', 'pos1', 'pos2', 'DAIC'] 

print(df.groupby('ensembl_id').size())
'''
###############################################################################
#Conservation of genes according to their amount of related pairs. See in R
###############################################################################
#Pairs with DAIC > 20

pairs_conserv = []
for i in conserv_val:
    count = 0
    for j in best_profiles_val:
        if i[0] == j[0]:
            count += 1
    pairs_conserv.append([i[0], count, i[1]])

coev_conserv = pd.DataFrame(pairs_conserv)
coev_conserv.columns = ['gene id', 'pairs', 'conservation']
coev_conserv.to_csv('coev_conserv_v2.csv', sep = ',')

###############################################################################
#Building data for machine learning
############################################################################### 

ml_data = []

for i in best_profiles_val:
    buffer = ''
    site_buffer = ''
    struct = get_array(i[0], motifs_val)
    which_motifs = get_motifs(i[1], i[2], struct)
    for j in conserv_val:
        if i[0] == j[0]:
            buffer = j[1]
    '''        
    for k in site_conserv_val:
        if i[1] and i[2] == k[1] and k[2]:
            site_buffer = k[3]
    '''
    ml_data.append([i[0], i[1], i[2], which_motifs, i[3], buffer])

ml_data_v2 = []   
for i in ml_data: #Removing the entries with NA values.
    if i[3] != '':
        ml_data_v2.append(i)




for i in ml_data_v2:
    for j in site_conserv_val:
        if i[0] == j[0]:
            if (i[1] and i[2]) == (j[1] and j[2]):
                print(str(j[3]))
                ml_data_v2.extend(str(j[3]))
                break




data = pd.DataFrame(ml_data_v2)
data.columns =  ['ensembl_id', 'pos1', 'pos2', 'motifs', 'DAIC', 'gene_conserv'] 

data = data.assign(site_conserv = site_conservation)

data.to_csv('ml_data_v2.csv', sep = ',')



