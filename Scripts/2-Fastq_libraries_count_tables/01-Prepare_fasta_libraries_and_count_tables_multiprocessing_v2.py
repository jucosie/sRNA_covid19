#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
- CREATE A UNIQUE SRNA IDENTIFIERS DATABASE.

Read the library fasta files and check for unique sequences removing sRNAs with 
indeterminacies (N) and keeping sRNAs within the length range of 12 to 34
 ntds. 
Then, build a big database assigning a unique identifier to each unique sRNA 
sequence. Example: Homo sapiens will be HSAXXXXXXXX.

*If the database already exists, it can be used as a parameter.

--------- FILTER 1

- FILTER THE LIBRARIES AND RENAME THE SRNAS SEQUENCES WITH THE UNIQUE SRNA 
IDENTIFIERS DATABASE.

First, read the library fasta files again, remove sRNAs with indeterminacies (N) 
and keeps sRNAs within the length range of 12 to 34 ntds. Then, rename the sRNA 
sequences across all the libraries using the unique sRNA identifiers database 
and save the new files.

*If you put the database in the command line as a parameter and a sequence is not
in the database, the pipeline will assign it a identifier and it will be put in 
the database

--------- FILTER 2 (OPTIONAL)

- FILTER THE LIBRARIES BY RFAM.

Align the libraries on the RFAM database and remove the sequences which align.

---------

- CREATE THE ABSOLUTE COUNT TABLE OF EACH LIBRARY.

Read the libraries filtered and create a absolute count table (csv format). The 
table has two columns: seq and counts.


- CREATE THE RPM COUNTS TABLE OF EACH LIBRARY.

Read the absolute count tables and create a rpm count table (csv format). The 
table has two columns: seq and RPMs. The RPMs = absolute count * 1000000 / size 
library.

- JOIN ALL THE RPMs TABLES AND ALL THE ABSOLUTE TABLES IN TWO TABLES.

Moreover, the pipeline creates a table with the mean between replicates.

---------

Created on Tue Feb 16 10:00:16 2021

Last Modified:
    - Thu Apr 01 11:05:27 2021 --> RFAM and table union.

@author: pasviber - Pascual Villalba Bermell
"""

## IMPORT MODULES

import os
from Bio import SeqIO
import numpy as np
import pandas as pd
import argparse
import multiprocessing
import itertools
from time import time
import sys


## FUNCTIONS

def ReadSummarySamples (path_summary):
    """
    This function reads a tsv table with the information of 
    each library (fasta filename, time, group, replica). 

    Parameters
    ----------
    path_summary : String
        path of the summary samples table.

    Returns
    -------
    summary : Dictionary
        Dictionary that stores the information of each library to construct 
        subsequent names.

    """
    
    f_in = open(path_summary, 'r')
    summary = {}
    for line in f_in:
        line = line.strip().split("\t")
        key = line[0]
        value = line[1:]
        summary[key] = value
    
    return summary
    

def AllUniqueSequencesListPROCESSING (libraries_list, path_read, num_process):
    """
    This function calls in parallel the function  AllUniqueSequencesList 
    function to determine which are the unique sequences in a library list.
    """
    
    number_files = len(libraries_list)
    distribution = number_files//num_process
    
    processes = []
    manager = multiprocessing.Manager()
    global_unique_seqs = manager.list()
        
    pi = 0
    for i in range(num_process):
        
        if i == num_process-1:
            pf = number_files
        else:
            pf = pi + distribution
        
        processes.append(multiprocessing.Process(target = AllUniqueSequencesList , args=(libraries_list[pi:pf], path_read, global_unique_seqs) )) 
        processes[i].start()
        print('Process', i, 'launched.')
        
        pi = pf
        
    for p in processes:
        p.join()
    
    unique_sequences = list(set(list(itertools.chain.from_iterable(global_unique_seqs))))
    unique_sequences.sort()
    
    print('\nNumber of accumulated unique sRNAs:', len(unique_sequences))
    print('\nDONE.')
    
    return unique_sequences


def AllUniqueSequencesList (libraries_list, path_read, global_unique_seqs):
    """
    This function goes through all the fasta files in a list and saving 
    in a global list all the unique sequences present in these files.

    Parameters
    ----------
    libraries_list : List
        List containing the fasta files from which we have to extract the unique sequences.
        As it is parallelized, each process will call this function with a different list of 
        fasta files.
    path_read : String
        Absolute path to the folder containing the fasta libraries of each sample. 
    global_unique_seqs : List of lists
        List that will contain different lists where the results of each launched process are stored.
        It is a shared list between processes.

    Returns
    -------
    global_unique_seqs : List of lists
        Explained in parameters.
    """
    
    ## List which will save all the unique sRNA sequences of the fasta files list.
    unique_seqs = []
    ## Iterate in fasta file list.
    for name in libraries_list:
        fasta = SeqIO.parse(path_read + name, 'fastq')
        ## List which will save all the unique sRNA sequences of the one fasta file.
        filtered_seqs = []
        ## Filter the fasta file by indeterminacies and length range.
        for sequence in fasta:
            seq = ''.join(sequence.seq)
            if 'N' in seq or 'n' in seq: # Remove the indeterminacies (N).
                continue
            elif len(seq) < 12 or len(seq) > 34: # Check the length of the sRNA.
                continue
            ## Save the seq filtered.
            filtered_seqs.append(seq)
        ## Add the seqs filtered to the main list and remove the sRNA duplicated 
        ## among files.
        unique_seqs = list(set(unique_seqs + filtered_seqs))
    
    ## The global list is shared among processes, Then we must save the main list,
    ## obtained in the function, in the global list.
    global_unique_seqs.append(unique_seqs)
    
    return global_unique_seqs


def CreationUniqueSRNAsDatabase (unique_sequences, prefix, max_digits):
    """
    This function is in charge of creating a database (dictionary) of all the 
    unique sequences present among all the libraries. The key will be the 
    sequence and the value will be the code generated for that sequence which 
    will consist of the prefix chosen for the species and normally 8 digits 
    starting from 1 to the number of sequences in the database.

    Parameters
    ----------
    unique_sequences : List
        List of the unique sRNA sequences of all the fasta libraries obtained 
        in the previous step and sorted in alphabetical order.
    prefix : String
        Prefix provided to create the unique identifier. 
    max_digits : Integer
        The number of digits that the unique identifier will contain (default 8).

    Returns
    -------
    database : Dictionary
        Dictionary built as a database of unique sRNA sequences.
    """
    
    print('Building the database...')
    ## Dictionary where we save the database of unique sRNAs
    database = {}
    ## Iterate in the unique sequences list length to provide the number which 
    ## will be assigned to the sRNA, localizated in that position minus 1 in the 
    ## list of unique sRNA sequences, together to a specific number of previous 
    ## zeros (the maximum number of digits is by default 8).
    for index in range(1, len(unique_sequences) + 1):
        n_digits = len(str(index))
        n_zeros = max_digits - n_digits
        value = str(prefix + '0' * n_zeros + str(index)) # e.g. >CUCUSA000000001
        key = unique_sequences[index - 1]
        database[key] = value
        
    print('\nDONE.')
        
    return database


def WriteDatabase (database, path, specie):
    """
    This function is in charge of creating a fasta file and a tsv file from the 
    created database of unique sequences.

    Parameters
    ----------
    database : Dictionary
        Dictionary with the database of unique sRNA sequences. The 
        key is the sequence and the value is the unique identifier created.
    path : String
        Absolute path of the directory where the fasta file with the database 
        will be saved. 
    specie : String
        Species we are working with. Example: homo sapiens.

    Returns
    -------
    name_1 : String
        Absolute path of the fasta file name where the database (FASTA) has 
        been saved.
    """
    
    print('Writing the database in a fasta file and a tsv file...')
    ## Create the name path of the database and open the file to write.
    name = path + 'all_sRNA_' + specie.replace(" ", "_") + '.tsv'
    tsv_database = open(name, 'w')
    ## Write the database fasta file and tsv file.
    for key in list(database.keys()): 
        tsv_database.write(database[key] + '\t' + key + '\n') # e.g. >CUCUSA000000001
    ## Close the fasta file.
    tsv_database.close()
    
    print('\nDONE.')
    print('\nName file (TSV):', 'all_sRNA_' + specie.replace(" ", "_") + '.tsv')
    print('Path:', path)
    
    return name


def ReadDatabase (path_database):
    """
    This function reads the database of unique identifiers entered through 
    the command line.

    Parameters
    ----------
    path_database : string
        Path of the unique identifier database.

    Returns
    -------
    db_dict : Dictionary
        Dictionary of the database with sequences as keys and identifiers 
        as values.
    """
    
    db = open(path_database, 'r')
    db_dict = {}
    
    for line in db:
        line = line.strip().split("\t")
        key = line[1]
        value = line[0]
        db_dict[key] = value
    
    return db_dict


def FilterFastaLibrariesANDRenameSRNAsPROCESSING (libraries_list, path_read, path_write, database, path_database, summary, specie, prefix, max_digits, num_process):
    """
    This function calls in parallel the function FilterFastaLibrariesANDRenameSRNAs 
    to filter the libraries and rename the sequences with the new unique identifiers.
    """
    
    number_files = len(libraries_list)
    distribution = number_files//num_process
    processes = []
        
    pi = 0
    for i in range(num_process):
        
        if i == num_process-1:
            pf = number_files
        else:
            pf = pi + distribution
        
        processes.append(multiprocessing.Process(target = FilterFastaLibrariesANDRenameSRNAs , args=(libraries_list[pi:pf], path_read, path_write, database, path_database, summary, specie, prefix, max_digits) )) 
        processes[i].start()
        print('Process', i, 'launched.')
        
        pi = pf
        
    for p in processes:
        p.join()
    
    print('\nDONE')
    

def FilterFastaLibrariesANDRenameSRNAs (libraries_list, path_read, path_write, database, path_database, summary, specie, prefix, max_digits):
    """
    This function filters the libraries by indeterminacy and sequence length range
    (12-34 ntds) (criteria also followed for the creation of the database), renaming 
    the filtered sequences with the new identifier using the database and to save 
    the new libraries in a new fasta file and in a new folder.

    Parameters
    ----------
    libraries_list : List
        List containing the fasta files in which we will have to filter and rename their sRNAs.  
        As it is parallelized, each process will call this function with a different list of fasta files.
    path_read : String
        Absolute path of the directory where the fasta files of the unfiltered libraries are located.
    path_write : String
        Absolute path of the directory where the new fasta files of the libraries will be saved.
    database : Dictionary
        Dictionary of the database.
    path_database : String
        Absolute path to the fasta file containing the database.
    specie : String
        Species we are working with.
    prefix : String
        Species prefix used to create the unique identifier of each sequence.
    max_digits : Integer
        Maximum number of digits that a sequence identifier can contain.

    Returns
    -------
    None.
    """
    
    ## Iterate the library fasta files list.
    for name in libraries_list:
        
        sample = "-".join(["_".join(summary[name][:-1]), summary[name][-1]])
        
        fasta_final = open(path_write + 'filtered_cleaned_' + specie.replace(" ", "_") + '_' + sample + '.fa', 'w')
        input_seq_iterator = SeqIO.parse(path_read + name, "fastq")
        
        ## Iterate the seqs of each library fasta file.
        for record in input_seq_iterator:
            seq = ''.join(record.seq)
            ## Filter by indeterminacies (N)
            if 'N' in seq or 'n' in seq:
                continue
            ## Filter by range of sRNA sequence.
            elif len(seq) < 12 or len(seq) > 34:
                continue
            ## If the seq passes the filter.
            else:
                ## Try to rename the sRNA sequence using the database.
                ## Example of identifier: CUCUSA00000001
                try:
                    new_id = str('>' + database[seq]) 
                ## If the seq is not in the database (This never happens if we
                ## have programmed the code correctly).
                except:
                    print(seq)
                    print("sRNA not found database. Appending to database...")
                    sys.exit()
                    last_id = list(database.values())[-1]
                    new_last_number = int(last_id[len(prefix):]) + 1 # Get rid of string CUCUSA from the id, keep only integer.
                    srna_file = open(path_database, 'a')
                    n_digits = len(str(new_last_number))
                    n_zeros = max_digits - n_digits
                    id_database = prefix + '0' * n_zeros + str(new_last_number)
                    srna_file.write(id_database + '\t' + seq + '\n')
                    database[seq] = id_database
                    srna_file.close()
                    new_id = '>' + id_database
            ## Write the new_id of the sRNAs of each library and the seq.
            fasta_final.write(new_id + '\n') # e.g. >CUCUSA000000001
            fasta_final.write(seq + '\n')
        ## Close the file.
        fasta_final.close()


def FilterRFAM (library_list, rfam_database, path_write):
    """
    This function is in charge of filtering the libraries by eliminating those 
    sequences that map against the RFAM database. 

    Parameters
    ----------
    library_list : list
        List of previously filtered fasta libraries with new identifiers.
    rfam_database : string
        Path to the RFAM database generated outside this script by another script.
        another script.
    path_write : string
        Path to the RFAM directory where all RFAM filtered libraries will be stored.

    Returns
    -------
    None.

    """
    path_align = path_write + "alignments/"
    os.system("mkdir -p " + path_align)
    os.chdir(path_align)
    path_index = path_align + "Index/"
    os.system("mkdir -p " + path_index)
    print("\nPercentage of sequences aligned on RFAM database:")
    for library in library_list:
        
        name = library.split("_")[-2] + "_" + library.split("_")[-1].split(".")[-2]
        
        os.system("rm *.fasta")
        os.system("rm *.sam")
        os.system("bowtie-build " + rfam_database + " ./Index/rfam")
        os.system("bowtie ./Index/rfam --best -v 0 -k 1 --norc -f " + library + " -S " + name + ".sam --al aligned_" + name + ".fasta --un ./../rfam_filtered_" + name + ".fasta")
        
        align = open(path_align + "aligned_" + name + ".fasta", "r")
        unalign = open(path_write + "rfam_filtered_" + name + ".fasta", "r")
        n1 = len([line for line in align if line.strip()[0] == ">"])
        n2 = len([line for line in unalign if line.strip()[0] == ">"])
        align.close()
        unalign.close()
        
        total = n1 + n2
        percentage_aligned = (n1/total)*100
        
        print("--" + name + ": " + str(percentage_aligned) + "%")
    
    os.system("rm *.fasta")
    os.system("rm *.sam")
    os.chdir(path_write)
    os.system("rm -R alignments")


def CreateDirectory01 (libraries_list, path, specie):
    """
    This function is in charge of creating the directory where the absolute 
    count tables of each library will be stored.

    Parameters
    ----------
    libraries_list : List
        List of the new built fasta files of each library.
    path : String
        Absolute path of the directory where the absolute count 
        tables will be stored.
    specie : String
        Species.

    Returns
    -------
    None.

    """
    
    try:
        os.mkdir(path)
    except FileExistsError:
        print('The folder', path.split("/")[-2], 'already exists.')
        pass
    
    for name in libraries_list:
        
        time = name.split("_")[-2]
        group = name.split("_")[-1].split("-")[-2]
        
        try:
            os.mkdir(path + time)
        except FileExistsError:
            print('The folder', time, 'already exists.')
            pass
        
        try:
            os.mkdir(path + time + "/" + group)
        except FileExistsError:
            print('The folder', group, 'already exists.')
            pass
    
    print('\nDONE.')


def AllAbsoluteCountsPROCESSING (libraries_list, path_read, path_write, specie, num_process):
    """
    This function is in charge of calling in parallel to the AllAbsoluteCounts 
    function to calculate the table of absolute counts of each library.
    """
    number_files = len(libraries_list)
    distribution = number_files//num_process
    
    processes = []
        
    pi = 0
    for i in range(num_process):
        
        if i == num_process-1:
            pf = number_files
        else:
            pf = pi + distribution
        
        processes.append(multiprocessing.Process(target = AllAbsoluteCounts , args=(libraries_list[pi:pf], path_read, path_write, specie) )) 
        processes[i].start()
        print('Process', i, 'launched.')
        
        pi = pf
        
    for p in processes:
        p.join()
    
    print('\nDONE')
    
    
def AllAbsoluteCounts (libraries_list, path_read, path_write, specie):
    """
    This function computes the table of absolute counts for each library 
    from the libraries filtered by indeterminacy and sequence length range.

    Parameters
    ----------
    libraries_list : List
        List that will contain the new libraries' fasta files. As it is parallelized, 
        each process will call this function with a different list of fasta files.
    path_read : String
        Absolute path to the directory where the new library fasta files are located.
    path_write : String
        Absolute path of the directory built in the previous step to store the absolute 
        count tables.
    specie : String
        Species.

    Returns
    -------
    None.
    
    """
    ## Iterate in the library fasta file list.
    for name in libraries_list:
        ## Create the path and the name where the count table will be saved.
        time = name.split("_")[-2]
        group = name.split("_")[-1].split("-")[-2]
        replica = name.split("_")[-1].split(".")[-2].split("-")[-1]
        new_name = group + "-" + str(replica) + ".abs_counts.csv"
        route_file_in = path_read + name
        route_file_out = path_write + time + "/" + group + "/" + new_name
        ## Create the count table using bash code. It's faster.
        command = "echo \"seq,counts\" > " + route_file_out
        os.system(command)
        command = "grep -v '>' " + route_file_in + " | sort | uniq -c | sort -nr | awk 'BEGIN{FS=\" \"; OFS=\",\"} {print $2, $1}' >> " + route_file_out
        os.system(command)


def CreateDirectory02 (libraries_list, path):
    """
    This function builds the directory where the RPM tables of each 
    library will be stored.

    Parameters
    ----------
    libraries_list : List
        List of the absolute paths of the csv files of the absolute count 
        tables.
    path : String
        Absolute path of the directory where the RPM tables will be stored.

    Returns
    -------
    None.

    """
    try:
        os.mkdir(path)
    except FileExistsError:
        print('The folder', path.split("/")[-2], 'already exists.')
        pass
    
    for route_file in libraries_list:
        
        time = route_file.split("/")[-3]
        group = route_file.split("/")[-2]
        
        try:
            os.mkdir(path + time)
        except FileExistsError:
            print('The folder', time, 'already exists.')
            pass
        
        try:
            os.mkdir(path + time + "/" + group)
        except FileExistsError:
            print('The folder', group, 'already exists.')
            pass
    
    print('DONE.')


def AllRPMCountsPROCESSING (libraries_list, path_write, num_process):
    """
    This function calls in parallel the AllRPMCounts function to calculate 
    the RPM table for each library.
    """
    number_files = len(libraries_list)
    distribution = number_files//num_process
    
    processes = []
        
    pi = 0
    for i in range(num_process):
        
        if i == num_process-1:
            pf = number_files
        else:
            pf = pi + distribution
        
        processes.append(multiprocessing.Process(target = AllRPMCounts , args=(libraries_list, path_write) )) 
        processes[i].start()
        print('Process', i, 'launched.')
        
        pi = pf
        
    for p in processes:
        p.join()
    
    print('\nDONE.')
    
    
def AllRPMCounts (libraries_list, path_write):
    """
    This function calculates the RPM table for each library from the absolute count 
    tables constructed in the previous step.

    Parameters
    ----------
    libraries_list : List
        List containing the absolute paths of the csv files of the absolute count 
        tables. As it is parallelized, each process will call this function with a 
        different list of csv files. 
    path_write : String
        Absolute path of the directory where the csv files of the RPM tables will be saved.

    Returns
    -------
    None.

    """
    ## Iterate in the route files list of absolute counts tables.
    for route_file_in in libraries_list:
        ## Create the path name of the rpm table.
        time = route_file_in.split("/")[-3]
        type_ = route_file_in.split("/")[-2]
        sample = route_file_in.split("/")[-1].split(".")[-3]
        new_name = sample + ".rpm_counts.csv"
        route_file_out = path_write + time + "/" + type_ + "/" + new_name
        ## Create the RPMs table.
        count_tab = pd.read_csv(route_file_in, sep = ",")
        abs_ = np.array(list(count_tab['counts']), dtype = float)
        rpms = (abs_*1000000)/sum(abs_)
        count_tab.insert(value = rpms, column='RPMs', loc=2)
        count_tab = count_tab.iloc[:,[0,2]]
        count_tab.columns = ['seq', 'RPMs']
        count_tab.to_csv(route_file_out, index = None)


def FusionTables (list_abs, list_RPMs, path_write, mode):
    """
    This function merges all absolute count tables into a single table and merges 
    all RPM tables into another single table. In addition, it also creates a table 
    for each case by calculating the average of the replicates.

    Parameters
    ----------
    list_abs : List
        List of absolute paths of the absolute counts table for each library.
    list_RPMs : List
        List of absolute paths from the RPMs table of each library.
    path_write : String
        Path of the directory where all the created tables will be stored.
    mode : string
        It is the union mode used for the union of the replicas of the same sample.
        sample. It can be inner (used for the PCA) or outer (used for the Differential 
        Expression Analysis).The union between different samples (groups) will always 
        be outer.

    Returns
    -------
    None.

    """
    ## Create a dict with all the absolute count tables opened. And create a time
    ## list (T1, T2...) and a group list (mock, HSVd).
    dic_abs = {}
    times = []
    groups = []
    for route in list_abs:
        list_route = route.split("/")
        time = list_route[-3]
        if time not in times:
            times.append(time)
        sample = list_route[-1].split(".")[0]
        group = sample.split("-")[0]
        if group not in groups:
            groups.append(group)
        name = time + "_" + sample
        table = pd.read_csv(route, sep = ",")
        dic_abs[name] = table
    ## Create a dict with all the RPMs tables opened.
    dic_RPMs = {}
    for route in list_RPMs:
        list_route = route.split("/")
        time = list_route[-3]
        sample = list_route[-1].split(".")[0]
        name = time + "_" + sample
        table = pd.read_csv(route, sep = ",")
        dic_RPMs[name] = table
    
    keys = sorted(list(dic_abs.keys()))
    times = sorted(times)
    groups = sorted(groups)
    interm_abs = {}
    interm_RPMs = {}
    samples = []
    
    ## Iterate in the times list and groups list.
    for t in times:
        for g in groups:
            i = 0
            ## Select the replicates of a specific time and group.
            replicas = [key for key in keys if t in key and g in key]
            if replicas != []:
                ## Obtain the name sample. Remove the number of the replicate.
                sample = replicas[0].split("-")[0]
                samples.append(sample)
                ## Join all the replicates in one table.
                for replica in replicas:
                    if i == 0:
                        fusion_abs = dic_abs[replica]
                        fusion_abs = fusion_abs.rename(columns = {'counts':replica})
                        fusion_RPMs = dic_RPMs[replica]
                        fusion_RPMs = fusion_RPMs.rename(columns = {'RPMs':replica})
                    else:
                        table_abs = dic_abs[replica]
                        table_abs = table_abs.rename(columns = {'counts':replica})
                        fusion_abs = fusion_abs.merge(table_abs, how=mode, left_on='seq', right_on='seq')
                        table_RPMs = dic_RPMs[replica]
                        table_RPMs = table_RPMs.rename(columns = {'RPMs':replica})
                        fusion_RPMs = fusion_RPMs.merge(table_RPMs, how=mode, left_on='seq', right_on='seq')

                    i += 1
                ## Save in a dict the replicate tables joined.
                interm_abs[sample] = fusion_abs
                interm_RPMs[sample] = fusion_RPMs
    
    ## Join all the samples in one only table. And calculate the mean to create
    ## another list.
    for x in range(len(samples)):
        ncols = len(interm_abs[samples[x]].columns)
        if x == 0:
            Global_abs = interm_abs[samples[x]]
            Global_mean_abs = Global_abs.iloc[:,[i for i in range(1, ncols)]].mean(axis=1)
            Global_mean_abs = pd.DataFrame(data=Global_mean_abs)
            Global_mean_abs.insert(value = Global_abs['seq'], column='seq', loc=0)
            Global_mean_abs.columns = ['seq', samples[x]]
            
            Global_RPMs = interm_RPMs[samples[x]]
            Global_mean_RPMs = Global_RPMs.iloc[:,[i for i in range(1, ncols)]].mean(axis=1)
            Global_mean_RPMs = pd.DataFrame(data=Global_mean_RPMs)
            Global_mean_RPMs.insert(value = Global_RPMs['seq'], column='seq', loc=0)
            Global_mean_RPMs.columns = ['seq', samples[x]]
            
        else:
            table_abs = interm_abs[samples[x]]
            table_mean_abs = table_abs.iloc[:,[i for i in range(1, ncols)]].mean(axis=1)
            table_mean_abs = pd.DataFrame(data=table_mean_abs)
            table_mean_abs.insert(value = table_abs['seq'], column='seq', loc=0)
            table_mean_abs.columns = ['seq', samples[x]]
            
            table_RPMs = interm_RPMs[samples[x]]
            table_mean_RPMs = table_RPMs.iloc[:,[i for i in range(1, ncols)]].mean(axis=1)
            table_mean_RPMs = pd.DataFrame(data=table_mean_RPMs)
            table_mean_RPMs.insert(value = table_RPMs['seq'], column='seq', loc=0)
            table_mean_RPMs.columns = ['seq', samples[x]]
            
            Global_abs = Global_abs.merge(table_abs, how='outer', left_on='seq', right_on='seq')
            Global_RPMs = Global_RPMs.merge(table_RPMs, how='outer', left_on='seq', right_on='seq')
            Global_mean_abs = Global_mean_abs.merge(table_mean_abs, how='outer', left_on='seq', right_on='seq')
            Global_mean_RPMs = Global_mean_RPMs.merge(table_mean_RPMs, how='outer', left_on='seq', right_on='seq')
    
    ## Save the tables in the directory.
    Global_abs.to_csv(path_write + "Fusion_abs-" + mode + ".csv", index = None)
    Global_RPMs.to_csv(path_write + "Fusion_RPMs-" + mode + ".csv", index = None)
    Global_mean_abs.to_csv(path_write + "Fusion_abs-MEAN-" + mode + ".csv", index = None)
    Global_mean_RPMs.to_csv(path_write + "Fusion_RPMs-MEAN-" + mode + ".csv", index = None)
    


## MAIN PROGRAM
        
def main():
    """
    Programa principal
    """
    
    parser = argparse.ArgumentParser(prog='FILTAB V2', 
                                     description='''This parallelized program creates a database \
                                     of unique sRNA sequence specific of specie, filter \
                                     and rename the sRNA sequences of each library, filter by \
                                     RFAM (optional), create the absolute and RPMs table of each library, \
                                     and join the all the tables in one only table.''',  
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-l", "--path-libraries", type=str, nargs=1,  
                        help="Absolute path of libraries folder.")
    parser.add_argument("-r", "--path-results", type=str, nargs=1, 
                        help="Absolute path of results folder.")
    parser.add_argument("-f", "--folder-initial-libraries", type=str, nargs=1, 
                        help="Folder that contains the libraries (fasta files) only \
                        filtered by quality and adapters")
    parser.add_argument("-s", "--specie", type=str, nargs=2, 
                        help="Name of the specie. Example: cucumis sativus.")
    parser.add_argument("-p", "--prefix-specie", type=str, nargs=1, 
                        help="Prefix that makes up the sRNA unique identifier. \
                        Example: CUCUSA.")
    parser.add_argument("-m", "--max-digits", type=int, nargs='?', const=8, default=8,
                        help="Maximum number of digits that make up the sRNA \
                        unique identifier (default: %(default)s)")
    parser.add_argument("-t", "--threads", type=int, nargs='?', const=1, default=1,
                        help="Number of threads (default: %(default)s)")
    parser.add_argument("-d", "--path-database", type=str, nargs=1,  
                        help="Absolute path of sRNAs database.")
    parser.add_argument("-x", "--path-summary", type=str, nargs=1,  
                        help="Absolute path of summary samples table.")
    parser.add_argument("-z", "--path-rfam", type=str, nargs=1,  
                        help="Absolute path of RFAM database.")
    
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    #args = parser.parse_args(["-l", "/mnt/f/", "-r", "/mnt/f/", "-f", "hola", "-s", "cucumis", "sativus", "-p", "CUCUSA", "-m", "8", "-x", "/mnt/h/"])
    
    try:
        path_database = args.path_database[0]
        D = 1
    except:
        print("\nIMPORTANT: No database has been inserted.")
        D = 0
        
    try:
        path_rfam = args.path_rfam[0]
        R = 1
    except:
        print("\nIMPORTANT: No RFAM database has been inserted.")
        R = 0
    
    try:
        path_lib = args.path_libraries[0]
        path_res = args.path_results[0]
        folder_read_1 = args.folder_initial_libraries[0]
        specie = " ".join(args.specie)
        prefix_specie = args.prefix_specie[0]
        max_digits = args.max_digits
        num_process = args.threads
        path_summary = args.path_summary[0]
    except:
        print("ERROR: You have inserted a wrong parameter or you are missing a parameter.")
        parser.print_help()
        sys.exit()
    
    
    ##########################################################################
    ## Paths.
    path_read_1 = path_lib + folder_read_1 + "/"
    folder_write_1 = "04-" + specie.replace(" ", "_") + "_libraries_clean"
    path_write_1 = path_lib + folder_write_1 + "/"

    ## List of fasta library files.
    fasta_files_list_1 = os.listdir(path_read_1)
    for file in fasta_files_list_1:
        if file.endswith(".fq"):
            continue
        else:
            fasta_files_list_1.remove(file)
    
    ## Summary samples dictionary.
    summary = ReadSummarySamples (path_summary)
    
    ##########################################################################
    ## Filter and rename the sequences with unique identifiers.
    if D == 0:
        
        ######################################################################
        print('\nLOOKING FOR THE UNIQUE SRNAS SEQUENCES ACROSS ALL THE LIBRARIES.\n')
        print('Looking for the unique sequences...')
        
        t_inicial = time()
        unique_sequences = AllUniqueSequencesListPROCESSING (fasta_files_list_1, path_read_1, num_process)
        t_final = time()
        t_ejecucion = t_final - t_inicial
        print('Time:', t_ejecucion)
        
        ######################################################################
        print('\nCREATING A DATABASE OF UNIQUE SRNAS WITH A UNIQUE IDENTIFIERS.\n')
        database = CreationUniqueSRNAsDatabase (unique_sequences, prefix_specie, max_digits)
        path_database = WriteDatabase (database, path_lib, specie)
        
    elif D == 1:
        print('\nREADING THE DATABASE OF UNIQUE SRNAS WITH A UNIQUE IDENTIFIERS.\n')
        database = ReadDatabase (path_database)
    
    ##########################################################################
    print('\nFILTERING THE LIBRARY FASTA FILES AND RENAME THE SRNAS WITH THE DATABASE.\n')
    
    try:
        os.mkdir(path_write_1)
    except FileExistsError:
        print('The folder', path_write_1.split("/")[-2], 'already exists.')
        pass
    
    t_inicial = time()
    FilterFastaLibrariesANDRenameSRNAsPROCESSING (fasta_files_list_1, path_read_1, path_write_1, database, path_database, summary, specie, prefix_specie, max_digits, num_process)
    t_final = time()
    t_ejecucion = t_final - t_inicial
    print('Time:', t_ejecucion)
    
    ##########################################################################
    ## Filter by RFAM (optional) and calculate the absolute and RPMs counts.
    if R == 0:
        
        ######################################################################
        ## Paths.
        path_read_abs = path_write_1
        folder_write_abs = "01-library_absolute_counts"
        path_write_abs = path_res + folder_write_abs + "/"
    
        ## List of fasta library files
        fasta_files_list_abs = os.listdir(path_read_abs)
        for file in fasta_files_list_abs:
            if file.endswith(".fa"):
                continue
            else:
                fasta_files_list_abs.remove(file)
        
        ######################################################################
        print('\nCREATING THE ABSOLUTE COUNT TABLE OF EACH LIBRARY.\n')
        
        print('Building the directory...\n')
        CreateDirectory01 (fasta_files_list_abs, path_write_abs, specie)
        
        print('\nCalculating the absolute counts...')
        t_inicial = time()
        AllAbsoluteCountsPROCESSING (fasta_files_list_abs, path_read_abs, path_write_abs, specie, num_process)
        t_final = time()
        t_ejecucion = t_final - t_inicial
        print('Time:', t_ejecucion)
        
        ######################################################################
        ## Paths.
        path_read_RPMs = path_write_abs
        folder_write_RPMs = "02-library_RPM_counts"
        path_write_RPMs = path_res + folder_write_RPMs + "/"
    
        ## List of csv library files
        csv_files_list_RPMs = []
        for root, dirs, files in os.walk(path_read_RPMs):
            for file in files:
                if file.endswith(".csv"):
                     csv_files_list_RPMs.append(os.path.join(root, file))
        
        ######################################################################
        print('\nCREATING THE RPM TABLE OF EACH LIBRARY.\n')
        
        print('Building the directory...\n')
        CreateDirectory02 (csv_files_list_RPMs, path_write_RPMs)
        
        print('\nCalculating the RPMs counts...')
        t_inicial = time()
        AllRPMCountsPROCESSING (csv_files_list_RPMs, path_write_RPMs, num_process)
        t_final = time()
        t_ejecucion = t_final - t_inicial
        print('Time:', t_ejecucion)
        
        ######################################################################
        ## Join the tables all True or not.
        ## Paths.
        path_read_join_abs = path_write_abs
        path_read_join_RPMs = path_write_RPMs
        folder_write_join = "03-Fusion_tables"
        path_write_join = path_res + folder_write_join + "/"
        ## List of csv library files
        csv_files_list_abs = []
        for root, dirs, files in os.walk(path_read_join_abs):
            for file in files:
                if file.endswith(".csv"):
                     csv_files_list_abs.append(os.path.join(root, file))
        csv_files_list_RPMs = []
        for root, dirs, files in os.walk(path_read_join_RPMs):
            for file in files:
                if file.endswith(".csv"):
                     csv_files_list_RPMs.append(os.path.join(root, file))
        
        ######################################################################
        print('\nJOINING THE RPMs AND ABSOLUTE TABLES IN TWO TABLES.\n')
        
        print('Building the directory...\n')
        try:
            os.mkdir(path_write_join)
        except FileExistsError:
            print('The folder', path_write_join.split("/")[-2], 'already exists.')
            pass
        
        print('\nJoining the tables mode = inner between replicates...')
        FusionTables (csv_files_list_abs, csv_files_list_RPMs, path_write_join, "inner")
        print('\nJoining the tables mode = outer between replicates...')
        FusionTables (csv_files_list_abs, csv_files_list_RPMs, path_write_join, "outer")
        
    elif R == 1:
        
        ######################################################################
        ## Paths.
        path_read_2 = path_write_1
        folder_write_2 = "05-RFAM_filtered"
        path_write_2 = path_lib + folder_write_2 + "/"
    
        ## List of fasta library files
        fasta_files_list_2 = []
        for root, dirs, files in os.walk(path_read_2):
            for file in files:
                if file.endswith(".fa"):
                     fasta_files_list_2.append(os.path.join(root, file))
        
        ######################################################################
        print('\nFILTER BY RFAM DATABASE.\n')
        print('Filtering...\n')
        FilterRFAM (fasta_files_list_2, path_rfam, path_write_2)
        
        ######################################################################
        ## Paths.
        path_read_abs = path_write_2
        folder_write_abs = "01-library_absolute_counts"
        path_write_abs = path_res + folder_write_abs + "/"
    
        ## List of fasta library files
        fasta_files_list_abs = os.listdir(path_read_abs)
        for file in fasta_files_list_abs:
            if file.endswith(".fasta"):
                continue
            else:
                fasta_files_list_abs.remove(file)
        
        ######################################################################
        print('\nCREATING THE ABSOLUTE COUNT TABLE OF EACH LIBRARY.\n')
        
        print('Building the directory...\n')
        CreateDirectory01 (fasta_files_list_abs, path_write_abs, specie)
        
        print('\nCalculating the absolute counts...')
        t_inicial = time()
        AllAbsoluteCountsPROCESSING (fasta_files_list_abs, path_read_abs, path_write_abs, specie, num_process)
        t_final = time()
        t_ejecucion = t_final - t_inicial
        print('Time:', t_ejecucion)
        
        ######################################################################
        ## Paths.
        path_read_RPMs = path_write_abs
        folder_write_RPMs = "02-library_RPM_counts"
        path_write_RPMs = path_res + folder_write_RPMs + "/"
    
        ## List of csv library files
        csv_files_list_RPMs = []
        for root, dirs, files in os.walk(path_read_RPMs):
            for file in files:
                if file.endswith(".csv"):
                     csv_files_list_RPMs.append(os.path.join(root, file))
        
        ######################################################################
        print('\nCREATING THE RPM TABLE OF EACH LIBRARY.\n')
        
        print('Building the directory...\n')
        CreateDirectory02 (csv_files_list_RPMs, path_write_RPMs)
        
        print('\nCalculating the RPMs counts...')
        t_inicial = time()
        AllRPMCountsPROCESSING (csv_files_list_RPMs, path_write_RPMs, num_process)
        t_final = time()
        t_ejecucion = t_final - t_inicial
        print('Time:', t_ejecucion)
        
        ######################################################################
        ## Join the tables all True or not.
        ## Paths.
        path_read_join_abs = path_write_abs
        path_read_join_RPMs = path_write_RPMs
        folder_write_join = "03-Fusion_tables"
        path_write_join = path_res + folder_write_join + "/"
        ## List of csv library files
        csv_files_list_abs = []
        for root, dirs, files in os.walk(path_read_join_abs):
            for file in files:
                if file.endswith(".csv"):
                     csv_files_list_abs.append(os.path.join(root, file))
        csv_files_list_RPMs = []
        for root, dirs, files in os.walk(path_read_join_RPMs):
            for file in files:
                if file.endswith(".csv"):
                     csv_files_list_RPMs.append(os.path.join(root, file))
        
        ######################################################################
        print('\nJOINING THE RPMs AND ASBOLUTE TABLES IN TWO TABLES.\n')
        
        print('Building the directory...\n')
        try:
            os.mkdir(path_write_join)
        except FileExistsError:
            print('The folder', path_write_join.split("/")[-2], 'already exists.')
            pass
        
        print('\nJoining the tables mode = inner between replicates...')
        FusionTables (csv_files_list_abs, csv_files_list_RPMs, path_write_join, "inner")
        print('\nJoining the tables mode = outer between replicates...')
        FusionTables (csv_files_list_abs, csv_files_list_RPMs, path_write_join, "outer")
        
    
## CALL THE MAIN PROGRAM.
if __name__ == '__main__':
    """
    Call to the main program.
    """
    main()
