# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 22:33:17 2017

@author: Karim
"""
import mysql.connector
import urllib
import os
import math

def codon(nuc_pos): #nucleotide position as argument
    codon_pos = math.ceil(nuc_pos/3) # It returns the smallest integer greater than or equal to x
    return codon_pos

###############################################################################

def get_header(file_data): #change from file to file_data
    list_headers = []

    for i in file_data: #change from file to file_data
        if i.find('ENS') == 0:
            a = i.split('\n')
            list_headers.append(a[0])

    for header in list_headers:
        req = urllib.request.Request('http://www.uniprot.org/uniprot/?query=%s&sort=score' % header, headers={'User-Agent': 'Mozilla/5.0'})
        url = urllib.request.urlopen(req)
        url_read = url.readlines()
        raw_webpage = b''.join(url_read).decode('utf-8')
        webpage = raw_webpage.split('\n')
        url.close()
        for line in webpage:
            if 'Sorry, no results found for your search term' in line:
                list_headers.remove(header) #headers relative to all species

    return list_headers

###############################################################################
'''
def get_sequence(msa):
    ens_id = []
    sequences = []
    length = []
    msa = msa.split('>')
    for i in msa2:
        if i.find('ENSP0') == 0 or i.find('ENSMUS') == 0:
            ens_id.append(''.join(i.split('\n')[0]))
            sequences.append(''.join(i.split('\n')[1:-2])) # -2 to avoid taking the double \n at the end of the sequence.
            sequences = [i.replace('-', '') for i in sequences]

    for i in sequences:
        length.append(int(len(i)/3))

    return ens_id, sequences, length
'''
###############################################################################
'''
def seq_length_comparision(msa):
    ens, sequences, nuc_length = get_sequence(msa2)
    list_entries = []
    for header in ens:
        req = urllib.request.Request('http://www.uniprot.org/uniprot/?query=%s&sort=score' % header, headers={'User-Agent': 'Mozilla/5.0'})
        url = urllib.request.urlopen(req)
        url_read = url.readlines()
        raw_webpage = b''.join(url_read).decode('utf-8')
        webpage = raw_webpage.split('\n')
        url.close()
        for line in webpage:
            if line.find('</script>') == 0:
                x = line.split(' ')
                list_entries.append(''.join(x[1][4:-1]))

    aa_length = []
    for entry in list_entries:
        req = urllib.request.Request('http://www.uniprot.org/uniprot/%s.txt' % entry, headers={'User-Agent': 'Mozilla/5.0'})
        url = urllib.request.urlopen(req)
        url_read = url.readlines()
        raw_txt_file = b''.join(url_read).decode('utf-8')
        txt_file = raw_txt_file.split('\n')
        url.close()
        for line in txt_file:
            if line.find('ID') == 0:
                aa_length.append(int(line.split()[3]))

    length = ''
    prob = ''
    for i, j in zip(nuc_length, aa_length):
        if i == j:
            length = 'OK'
'''
###############################################################################

def get_msa_features(msa):
    seq_align = []
    species = []
    tmp = msa.split('>')
    for i in tmp:
        if 'ENS' in i:
            a = i.split('\n')
            species.append(a[0])
            seq_align.append(''.join(a[1:]))

    return len(species), len(seq_align[0]), species, seq_align

###############################################################################

def get_conservation(msa):
    species_num, site_num, species, sequence_align = get_msa_features(msa)
    siteIdentity = 0 # It gives the numbers of fully conserved sites.
    for col in range(site_num):
        nuc = [0,0,0,0]
        for row in range(species_num):
            if sequence_align[row][col] == 'A':
                nuc[0] += 1
            elif sequence_align[row][col] == 'C':
                nuc[1] += 1
            elif sequence_align[row][col] == 'G':
                nuc[2] += 1
            elif sequence_align[row][col] == 'T':
                nuc[3] += 1

        if max([x/int(species_num) for x in nuc]) == 1:
            siteIdentity += 1

    return siteIdentity/int(site_num) # % of identity in a sequence alignment

###############################################################################

def get_conservation_pairs(msa, pos1, pos2):
    species_num, site_num, species, sequence_align = get_msa_features(msa)
    nucl1 = [0,0,0,0]
    nucl2 = [0,0,0,0]

    for row in range(species_num):
        if sequence_align[row][pos1] == 'A':
            nucl1[0] += 1
        elif sequence_align[row][pos1] == 'C':
            nucl1[1] += 1
        elif sequence_align[row][pos1] == 'G':
            nucl1[2] += 1
        elif sequence_align[row][pos1] == 'T':
            nucl1[3] += 1

        if sequence_align[row][pos2] == 'A':
            nucl2[0] += 1
        elif sequence_align[row][pos2] == 'C':
            nucl2[1] += 1
        elif sequence_align[row][pos2] == 'G':
            nucl2[2] += 1
        elif sequence_align[row][pos2] == 'T':
            nucl2[3] += 1

    buffer1 = max(nucl1)/species_num
    buffer2 = max(nucl2)/species_num
    conserv_pairs = (buffer1 + buffer2)/2

    return conserv_pairs # % of identity for 2 pairs of sites in a MSA

###############################################################################

def get_uniprot_entries(ens_headers):
    list_uniprot_entries = []
    for header in ens_headers:
        req = urllib.request.Request('http://www.uniprot.org/uniprot/?query=%s&sort=score' % header, headers={'User-Agent': 'Mozilla/5.0'})
        url = urllib.request.urlopen(req)
        url_read = url.readlines()
        raw_webpage = b''.join(url_read).decode('utf-8')
        webpage = raw_webpage.split('\n')
        url.close()
        for line in webpage:
            if line.find('</script>') == 0 and '_HUMAN' in line:
                x = line.split(' ')
                list_uniprot_entries.append(''.join(x[1][4:-1])) #Entry relative to human only

    return list_uniprot_entries

###############################################################################

def get_entry_names(uniprot_entries):
    list_entry_names = []
    for i in uniprot_entries:
        req = urllib.request.Request('http://www.uniprot.org/uniprot/%s.txt' % i, headers={'User-Agent': 'Mozilla/5.0'})
        url = urllib.request.urlopen(req)
        url_read = url.readlines()
        raw_txt_file = b''.join(url_read).decode('utf-8')
        txt_file = raw_txt_file.split('\n')
        url.close()
        for line in txt_file:
            if line.find('ID') == 0:
                a = line.split()
                list_entry_names.append(a[1])

    return list_entry_names

###############################################################################

def get_features(uniprot_entries):
    mol_f = []
    bio_p = []
    motifs = []
    for i in uniprot_entries:
        req = urllib.request.Request('http://www.uniprot.org/uniprot/%s.txt' % i, headers={'User-Agent': 'Mozilla/5.0'})
        url = urllib.request.urlopen(req)
        url_read = url.readlines()
        raw_txt_file = b''.join(url_read).decode('utf-8')
        txt_file = raw_txt_file.split('\n')
        url.close()
        function = False
        process = False
        structure = False
        for line in txt_file:
            if line.find('DR') == 0 and 'GO;' in line:
                if 'F:' in line:
                    function = True
                elif 'P:' in line:
                    process = True
            elif line.find('FT') == 0:
                if 'HELIX' in line or 'STRAND' in line or 'TURN' in line:
                    structure = True

        if function:
            mol_f.append('Yes')
        else:
            mol_f.append('No')
        if process:
            bio_p.append('Yes')
        else:
            bio_p.append('No')
        if structure:
            motifs.append('Yes')
        else:
            motifs.append('No')

    return mol_f, bio_p, motifs

###############################################################################
# Creation of a list of all possible function and processes

'''def get_all_functions(uniprot_entries):
    all_molecular_func = []
    all_biological_proc = []
    for i in uniprot_entries:
        req = urllib.request.Request('http://www.uniprot.org/uniprot/%s.txt' % i, headers={'User-Agent': 'Mozilla/5.0'})
        url = urllib.request.urlopen(req)
        url_read = url.readlines()
        raw_txt_file = b''.join(url_read).decode('utf-8')
        txt_file = raw_txt_file.split('\n')
        url.close()
        for line in txt_file:
            if line.find('DR') == 0 and 'GO;' in line:
                if 'F:' in line:
                    #a = line.split(';')
                    all_molecular_func.append(line.split(';')[2][3:])
                elif 'P:' in line:
                    #a = line.split(';')
                    all_biological_proc.append(line.split(';')[2][3:])

    return all_molecular_func, all_biological_proc'''

###############################################################################

def get_functions(uniprot_entries):
    molecular_func = []
    uni_entries = []
    for i in uniprot_entries:
        req = urllib.request.Request('http://www.uniprot.org/uniprot/%s.txt' % i, headers={'User-Agent': 'Mozilla/5.0'})
        url = urllib.request.urlopen(req)
        url_read = url.readlines()
        raw_txt_file = b''.join(url_read).decode('utf-8')
        txt_file = raw_txt_file.split('\n')
        url.close()
        for line in txt_file:
            if line.find('DR') == 0 and 'GO;' in line:
                if 'F:' in line:
                    molecular_func.append(line.split(';')[2][3:])
                    uni_entries.append(i)

    return uni_entries, molecular_func

###############################################################################

def get_processes(uniprot_entries):
    biological_proc = []
    uni_entries = []
    for i in uniprot_entries:
        req = urllib.request.Request('http://www.uniprot.org/uniprot/%s.txt' % i, headers={'User-Agent': 'Mozilla/5.0'})
        url = urllib.request.urlopen(req)
        url_read = url.readlines()
        raw_txt_file = b''.join(url_read).decode('utf-8')
        txt_file = raw_txt_file.split('\n')
        url.close()
        for line in txt_file:
            if line.find('DR') == 0 and 'GO;' in line:
                if 'P:' in line:
                    biological_proc.append(line.split(';')[2][3:])
                    uni_entries.append(i)

    return uni_entries, biological_proc

###############################################################################

def get_motifs(uniprot_entries):
    struct = []
    aa_start = []
    aa_stop = []
    uni_entries = []
    for i in uniprot_entries:
        req = urllib.request.Request('http://www.uniprot.org/uniprot/%s.txt' % i, headers={'User-Agent': 'Mozilla/5.0'})
        url = urllib.request.urlopen(req)
        url_read = url.readlines()
        raw_txt_file = b''.join(url_read).decode('utf-8')
        txt_file = raw_txt_file.split('\n')
        url.close()
        for line in txt_file:
            if line.find('FT') == 0:
                if 'HELIX' in line:
                    a = line.split()
                    struct.append(a[1]) #Variable matifs already assigned!
                    aa_start.append(a[2])
                    aa_stop.append(a[3])
                    uni_entries.append(i)
                elif 'STRAND' in line:
                    a = line.split()
                    struct.append(a[1]) #Variable matifs already assigned!
                    aa_start.append(a[2])
                    aa_stop.append(a[3])
                    uni_entries.append(i)
                elif 'TURN' in line:
                    a = line.split()
                    struct.append(a[1]) #Variable matifs already assigned!
                    aa_start.append(a[2])
                    aa_stop.append(a[3])
                    uni_entries.append(i)

    return uni_entries, struct, aa_start, aa_stop

###############################################################################

def get_ens_names():
    import mysql.connector
    cnx = mysql.connector.connect(host = '***', user = '***', password = '***', database = 'coevolution')
    cur = cnx.cursor()
    cur.execute('''SELECT NAME FROM ensembl WHERE ID IN (
            SELECT DISTINCT ID_ENS FROM best_profiles)''')

    ens_names = cur.fetchall()

    cnx.close()
    cur.close()

    list_ens_names = [i[0] + '.nt' for i in ens_names]

    return list_ens_names

###############################################################################
#update_tables('/Users/karimsaied/Desktop/Biology_v2/master/first_step_project/data/full_selectome_v06_Euteleostomi-nt_unmasked')

def update_tables(directory):
    cnx = mysql.connector.connect(host = '***', user = '***', password = '***', database = 'coevolution')
    cur = cnx.cursor()
    ens_names = get_ens_names()
    disable_files = []
    count_files = 0
    for file in os.listdir(directory):
        ensembl_id = ''
        count_files = count_files + 1
        if file.endswith('.nt.fas') and file[0:-4] in ens_names:
            a = file.split('.')
            ensembl_id = '%s.%s.%s' % (a[0], a[1], a[2]) #Store the Ensembl id (part of the file name) in a list
            filepath = os.path.join(directory, file)
            msa_ens = open(filepath, 'r') #filepath
            msa = msa_ens.read()
            msa_ens.close()
            file_data = msa.split('>')

            try:
                ens_headers = get_header(file_data)
                conservation = get_conservation(msa)
                pairs, site_conservation = get_pairs(ensembl_id, msa)
                uniprot_entries = get_uniprot_entries(ens_headers)
                #all_functions, all_processes = get_all_functions(uniprot_entries)
                entries_func, molecular_function = get_functions(uniprot_entries)
                entries_proc, biological_processes = get_processes(uniprot_entries)
                entry_names = get_entry_names(uniprot_entries)
                mol_function, bio_process, struct = get_features(uniprot_entries)
                entries, motifs, start, stop = get_motifs(uniprot_entries)

                for entry in range(len(entries)):
                    update_motif = 'INSERT INTO Motifs (Ensembl_id, UniProt_id, Sec_struct_pattern, Aa_start, Aa_stop) VALUES (%s, %s, %s, %s, %s)'
                    cur.execute(update_motif, (ensembl_id, entries[entry], motifs[entry], start[entry], stop[entry]))
                    cnx.commit()

                for entry in range(len(entries_func)):
                    update_functions = 'INSERT INTO Functions (Ensembl_id, UniProt_id, Molecular_function) VALUES (%s, %s, %s)'
                    cur.execute(update_functions, (ensembl_id, entries_func[entry], molecular_function[entry]))
                    cnx.commit()

                for entry in range(len(entries_proc)):
                    update_processes = 'INSERT INTO Processes (Ensembl_id, UniProt_id, Biological_process) VALUES (%s, %s, %s)'
                    cur.execute(update_processes, (ensembl_id, entries_proc[entry], biological_processes[entry]))
                    cnx.commit()

                for i in range(len(uniprot_entries)):
                    update_table = 'INSERT INTO Genes (Ensembl_id, UniProt_id, Gene_name, Secondary_structure, Biological_process, Molecular_function) VALUES (%s, %s, %s, %s, %s, %s)'
                    cur.execute(update_table, (ensembl_id, uniprot_entries[i], entry_names[i], struct[i], bio_process[i], mol_function[i]))
                    cnx.commit()

                for i in range(len(site_conservation)):
                    position1 = ''
                    position2 = ''
                    position1 = pairs[i][0]
                    position2 = pairs[i][1]
                    update_site_conserv = 'INSERT INTO Sites_conserv (Ensembl_id, Position1, Position2, Sites_conservation) VALUES (%s, %s, %s, %s)'
                    cur.execute(update_site_conserv, (ensembl_id, position1, position2, site_conservation[i]))
                    cnx.commit()

                update_conservation = 'INSERT INTO Conservation (Ensembl_id, Conservation) VALUES (%s, %s)'
                cur.execute(update_conservation, (ensembl_id, conservation))
                cnx.commit()

            except:
                print("Something was wrong")
                disable_files.append(file)
            print(count_files)
    print("**********Files to process again**********")
    print(disable_files)

    cur.close()
    cnx.close()

###############################################################################
'''
import panda as pd
import os

def get_msa(directory, best_profiles_path): # The arguments are paths.
    ens_names = get_ens_names()
    ens_id = get_ens_id_conserv()
    ensembl_names = []
    best_profiles_nuc = pd.read_csv(best_profiles_path)
    best_profiles_nuc = best_profiles_nuc.drop([
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

    best_profiles_val = best_profiles_nuc.values.tolist()

    for row in best_profiles_val:
        row[0], row[1], row[2] = int(row[0]), int(row[1]), int(row[2])


    #list_msa = []
    for file in os.listdir(directory):
        msa = []
        if file.endswith('.nt.fas') and file[0:-4] in ens_names:
            ensembl_names.append(file[0:-4])
            a = file.split('.')
            filepath = os.path.join(directory, file)
            msa_ens = open('ENSGT00390000001046.Euteleostomi.004.nt.fas', 'r') #filepath
            raw_msa = msa_ens.read()
            msa_ens.close()
            file_data = raw_msa.split('>')
            for i in file_data:
                msa.append(''.join(i.split('\n')[1:-2])) # -2 to avoid taking the double \n at the end of the sequence.

            msa = msa[1:]
            #list_msa.append(msa)
            break
'''
############################################

def get_pairs(name_file, msa): # msa should correspond to raw_msa (see main function)
    import mysql.connector
    pairs_pos = []
    a = "SELECT POSITION1, POSITION2 FROM best_profiles WHERE ID_ENS IN (SELECT ID FROM ensembl WHERE NAME = '%s')" % name_file
    cnx = mysql.connector.connect(host = '***', user = '***', password = '***', database = 'coevolution')
    cur = cnx.cursor()
    cur.execute(a)

    pairs_pos = cur.fetchall()

    cnx.close()
    cur.close()

    #return pairs_pos

    sites_conserv = []
    pairs = []
    for i in pairs_pos:
        sites_conserv.append(get_conservation_pairs(msa, i[0], i[1]))
        pairs.append([i[0], i[1]])

    return pairs, sites_conserv
