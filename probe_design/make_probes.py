#!/usr/bin/env python
# coding: utf-8


# DNA_probe functiuons 

import os
import sys
import numpy as np
import pandas as pd
import gzip
from scipy.optimize import curve_fit
import seaborn as sns
import math as m
import scipy.stats as stats
import pickle
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from shutil import copyfile
import matplotlib.colors as mcolors
import random
import copy
import matplotlib.patches as mpatches
from collections import Counter
import matplotlib.pyplot as plt
from fuzzywuzzy import fuzz
import textdistance
import time
from scipy.optimize import curve_fit
import statistics as stat
import plotly.graph_objects as go
import plotly.offline as pyo
import itertools
import csv
import gc
import shelve
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing

def find_nth_n_str(string,
                   leter,
                   n):
    m=0
    for c, char in enumerate(string):
        if m<n:
            if char == leter:
                m+=1
            else:
                pass
        else:
            return(c-1)
        
def reverse_complement(dna):
    complement = {'A': 'T','a': 'T', 'C': 'G','c': 'G', 'G':'C','g': 'C','T': 'A','t': 'A','N':'N','n':'N'}
    return ''.join([complement[base] for base in dna[::-1]])

def save_pickle(obj,
                file_path,
                file_name):
    """ saves objects as pickle files """
    path = os.path.normpath(os.path.join(file_path,file_name+'.pk'))
    pkl = open(path,'wb')
    pickle.dump(obj,pkl,protocol=4)
    pkl.close()
    
def load_pickle(file_path):
    """ loads objects from pickle files """
    if os.path.exists(file_path) != True:
        obj = None
    else:
        pkl = open(file_path,'rb')
        obj = pickle.load(pkl)
        pkl.close()
    return(obj)

def Average(lst): 
    return sum(lst) / len(lst)

# =============================================================================

# =============================================================================

# =============================================================================

def hamming(seq1,seq2):
    distance = 0
    for base1, base2 in zip(seq1,seq2):
        if base1 != base2:
            distance += 1
    return(distance)

def Levenshtein_ratio(seq1,seq2):
    ratio = fuzz.ratio(seq1,seq2)
    return(ratio) 

def getKmer(seq, k):
    if k>len(seq):
         return([])
    kmer = []
    for i in range(len(seq)-k):
        kmer.append(seq[i:i+k])
    return(kmer)

def calq_plot_seq_diversity(UMIs,plot = True):
    number_uniqu = []
    read_num = 0
    read_numList = []
    reads_uni = {}
    uniq_read_count = 0
    t0 = time.time()
    totalreads = len(readDF['all_reads2'].tolist())
    random.seed('DAF')
    reads2 =  random.sample(UMIs,k = len(UMIs))

    for read in reads2:
        if time.time()>t0:
            print('fraction done = '+str(read_num/totalreads))
            t0+=30
        read_num += 1
        if reads_uni.get(read) == None:
            reads_uni[read] = 1
            uniq_read_count+=1
        else:
            continue
        read_numList.append(read_num)
        number_uniqu.append(uniq_read_count)
    
    if plot == True:
        fig, ax = plt.subplots(figsize=(7,7))
        ax = plt.scatter(read_numList,number_uniqu,s = 15,color = 'black')
        fig.set_facecolor('white')
        plt.xlabel('read #',fontsize=16)
        plt.ylabel('unique read',fontsize=16)
        #plt.xscale('log')
        #plt.yscale('log')
        plt.show()
    return(read_numList,number_uniqu)

def Average(lst): 
    return sum(lst) / len(lst)

# =============================================================================

def getter_setter_gen(name, type_):
    '''Used to check type of class objects. '''
    def getter(self):
        return getattr(self, "__" + name)
    def setter(self, value):
        if not isinstance(value, type_):
            raise TypeError(f"{name} attribute must be set to an instance of {type_}")
        setattr(self, "__" + name, value)
    return property(getter, setter)

def auto_attr_check(cls):
    '''Class decorator to set type of class objects. '''
    new_dct = {}
    for key, value in cls.__dict__.items():
        if isinstance(value, type):
            value = getter_setter_gen(key, value)
        new_dct[key] = value
    # Creates a new class, using the modified dictionary as the class dict:
    return type(cls)(cls.__name__, cls.__bases__, new_dct)

# =============================================================================

class gene(object):
    """A class that contains data on a single gene
    
    Properties:
    ID = str ID used for genes
    TS_ID = list all transcript IDs for gene 
    chr = str Chromosome in form 'chr#'
    strand = str 1 or -1
    start = int start base on geneome
    stop = int stop base on geneome
    seq = str sequence of gene
    full_name = str 
    symbol = str
    biotype = str
    exp = float expresion of gene 
    
    Methods:
    info() returns above information about gene
    
    """
    ID = str
    TS_ID = list
    chr = str
    strand = str
    start = int
    stop = int
    seq = str 
    full_name = str 
    symbol = str
    biotype = str
    exp = float
    def __init__(self):
        self.ID = ''
        self.TS_ID = []
        self.chr = ''
        self.strand = ''
        self.start = int(-1)
        self.stop = int(-1)
        self.seq = ''
        self.full_name = '' # geneID: [chromosome, strand, start, stop] 
        self.symbol = ''
        self.biotype = ''
        self.exp = 0.0
    def info(self):
        print('ID - '+str(self.ID))
        print('TS_ID - '+str(self.TS_ID))
        print('chromosome(chr) - '+str(self.chr))
        print('strand - '+str(self.strand))
        print('start - '+str(self.start))
        print('stop - '+str(self.stop))
        print('full_name - '+str(self.full_name))
        print('symbol - '+str(self.symbol))
        print('biotype - '+str(self.biotype))
        print('expression(exp) - '+str(self.exp))

# =============================================================================

class transcript(object):
    """A class that contains data on a single transcript
    
    Properties:
    seq = str spliced sequence of transcript
    start = int start base on gene
    stop = int start base on gene
    ID = str ID used for transcript
    len = int length of spliced transcript in bases 
    gene_ID = str ID used for gene that transcript belongs to
    exp = float expresion of transcript 
    in_exon_index_start = list ordered list of exon start position in order of rank in transcript
    in_exon_index_stop = list ordered list of exon stop position in order of rank in transcript
    ex_exon_index_start = list ordered list of exon start position in order of rank in geneome
    ex_exon_index_stop = list ordered list of exon stop position in order of rank in geneome
    constitutive_exon = list ordered list of exons that exist in all transcripts of gene
    biotype = str biological classification of transcript type 
    
    Methods:
    info() returns above information about transcript
    
    """ 
    seq = str
    start = int
    stop = int
    ID = str
    len = int
    gene_ID = str
    exp = float
    in_exon_index_start = list 
    in_exon_index_stop = list 
    ex_exon_index_start = list
    ex_exon_index_stop = list
    constitutive_exon = list 
    biotype = str
    
    def __init__(self):
        self.seq = ''
        self.start = int(-1)
        self.stop = int(-1)
        self.ID = ''
        self.len = len(self.seq)
        self.gene_ID = ''
        self.exp = 0.0
        self.in_exon_index_start = [] 
        self.in_exon_index_stop = [] 
        self.ex_exon_index_start = []
        self.ex_exon_index_stop = []
        self.constitutive_exon = [] 
        self.biotype = ''
    def info(self):
        print('ID - '+str(self.ID))
        print('gene_ID - '+str(self.gene_ID))
        print('start - '+str(self.start))
        print('stop - '+str(self.stop))
        print('len - '+str(self.len))
        print('exp - '+str(self.exp))
        print('in_exon_index_start - '+str(self.in_exon_index_start))
        print('in_exon_index_stop - '+str(self.in_exon_index_stop))
        print('ex_exon_index_start - '+str(self.ex_exon_index_start))
        print('ex_exon_index_stop - '+str(self.ex_exon_index_stop))
        print('constitutive_exon - '+str(self.constitutive_exon))
        print('biotype - '+str(self.biotype))

# =============================================================================
        
class Transcriptome(object):
    """ Transcriptome - class to store objects of class transcript and gene. 
    Has methods to calculate thermodynamic properties and homology of all transcripts enclosed
    
    Properties:
    TS = dict transcript ID : transcript() 
    G = dict  gene ID : gene()
    IDs = dict gene IDs and transcript IDs : []
    Gsym2ID = dict gene symbol : gene ID
    TDP_db_path = str file path to the database of thermodynamic properties
    OT_db_path = str file path to the database of off target calulations 
    cell_type = str cell type used for expresion values 
    
    Methods:
    get_TDP(file_path, replace = False, verbose = False)
        file_path: str path of file location to store thermodynamic properties
        replace: bool recalulate thermodynamic properties and replace if present in file_path
        verbose: bool print progress of calulation every 20s 
    get_OT_index(file_path, OT_dict_file_name = 'OT_dict', OT_length = 17, 
                 replace = False, verbose = False, frag = 10)
        file_path: str path of file location to store off target calulations
        OT_dict_file_name: str name of file to make and store off target calulations
        OT_length: str length of off target subsections to calulate off target score for
        replace: bool recalulate off target score and replace if present in file_path
        verbose: bool print progress of calulation every 20s 
        frag: int number of fragments to break up transcriptome for calulations 
        threads: int number of cores to use if -1 use all cores
    subset(subset_list = None ,index = [0,0])
        subset_list: list of transcript IDs, gene IDs, biotypes, and/or 
        gene symbols to subset on
        index: list of two numbers of indicies of transcripts to subset 
        returns subseted transcriptome
    append(TS = None,Gene = None,Transcriptome2 = None)
        TS: list of transcript objects to append to transcriptome 
        Gene: list of gene objects to append to transcriptome
        Transcriptome2: Transcriptome object to append to transcriptome
        returns appended transcriptome
    G_biotypes()
        returns dict of biotype: count 
    TS_biotypes()
        returns dict of biotype: count 
    genes()
        returns a list of all gene IDs in in transcriptome 
    transcripts()
        returns a list of all transcript IDs in in transcriptome
    n_genes()
        returns number of genes in transcriptome
    n_transcripts():
        returns number of transcripts in transcriptome
    info()
        returns cell_type, number of transcripts, number of genes, 
        thermodynamic properties database path, and off target properties path
        
    """
    TS = dict 
    G = dict 
    IDs = dict 
    Gsym2ID = dict
    TDP_db_path = str
    OT_db_path = str
    cell_type = str
    
    def __init__(self): 
        self.TS = dict() 
        self.G = dict() 
        self.IDs = dict() 
        self.Gsym2ID = dict()
        self.TDP_db_path = ''
        self.OT_db_path = ''
        self.cell_type = ''

    
    def get_TDP(self, file_path, replace = False, verbose = False): # see get_TDP_CG function 
        if self.TDP_db_path == '':
            print(isinstance(self,Transcriptome))
            DB_path = get_TDP_CG(self,file_path,verbose=verbose)
            self.TDP_db_path = DB_path
        elif self.TDP_db_path != '' and replace == False:
            print('TDP_db_path already exist, change_file path or set replace = True')
        elif self.TDP_db_path != '' and replace == True:
            DB_path = get_TDP_CG(self,file_path,verbose=verbose)
            self.TDP_db_path = DB_path
    def get_OT_index(self, file_path, OT_dict_file_name = 'OT_dict', OT_length = 17, replace = False, verbose = False, frag = 10, threads = -1):
        if self.OT_db_path == '':
            OT_dict = Make_OT_dict_Full(self, file_path = file_path, OT_dict_file_name = OT_dict_file_name ,OT_length = OT_length, verbose = verbose, frag = frag, threads = threads) # see Make_OT_dict_Full function 
            OT_db_path = get_OT_dict_Iso(self, OT_dict, file_path, verbose = verbose, frag = frag,threads = threads) # see  get_OT_dict_Iso function 
            self.OT_db_path = OT_db_path
        elif self.OT_db_path != '' and replace == False:
            print('OT_db_path already exist, change file_path or set replace = True')
        elif self.OT_db_path != '' and replace == True:
            OT_dict = Make_OT_dict_Full(self, file_path = file_path, OT_dict_file_name = OT_dict_file_name ,OT_length = OT_length, verbose = verbose, frag = frag, threads = threads)
            OT_db_path = get_OT_dict_Iso(self, OT_dict, file_path, verbose = verbose, frag = frag,threads = threads)
            self.OT_db_path = OT_db_path
    def subset(self, subset_list = None ,index = [0,0]): # see Subset_Transcriptom function 
        subset = Subset_Transcriptome(self,subset_list = subset_list,index = index)
        print('database at TDP_db_path and OT_db_path are not subseted')
        return(subset)
    def append(self,TS = None,Gene = None,Transcriptome2 = None):
        appended = append_transcriptome(self,TS = TS,Gene = Gene,Transcriptome2 = Transcriptome2)
        print('database at TDP_db_path and OT_db_path are not appended. Must recalulate to update')
        return(appended)
    def G_biotypes(self):
        BT = dict()
        for G in self.G:
            G=self.G[G]
            bioty = G.biotype
            if BT.get(bioty) == None:
                BT[bioty] = 1
            else:
                BT[bioty] += 1
        return(BT)
    def TS_biotypes(self):
        BT = dict()
        for TS in self.TS:
            TS=self.TS[TS]
            bioty = TS.biotype
            if BT.get(bioty) == None:
                BT[bioty] = 1
            else:
                BT[bioty] += 1
        return(BT)
    def genes(self):
        gids = list(self.G.keys())   
        return(gids)
    def transcripts(self):
        tsid = list(self.TS.keys())
        return(tsid)
    def n_genes(self):
        glen = len(self.G)
        return(glen)
    def n_transcripts(self):
        tslen = len(self.TS)
        return(tslen)
    def info(self):
        print('cell_type - '+ self.cell_type)
        print('Number of transcripts - '+ str(self.n_transcripts()))
        print('Number of genes - '+ str(self.n_genes()))
        print('Thermodynamic properties(TDP_db_path) @ '+ self.TDP_db_path)
        print('Off target properties(OT_db_path) @ '+ self.OT_db_path)

# =============================================================================

#cell_expr_iso_dirt: {ensembl_transcriptID:expression(TPM)}

#cdna_fasta: from ensembl website - Homo_sapiens.GRCh38.cdna.all.fa

#geneome_fasta: from ensembl biomart in FASTA header form (unspliced gene) sequence (sep = |):  
  #geneID(no version number)
  #gene name
  #discription
  #chr
  #start(bp)
  #stop(bp)
  #gene type

#exon_CSV: from ensembl biomart in CSV form (sep = ,): 
 #Gene stable ID (no version number)
 #Transcript stable ID (no version number)
 #Strand,Transcript start (bp)
 #Transcript end (bp)
 #Exon region start (bp)
 #Exon region end (bp)

def Build_Transcriptome(cell_expr_iso_dirt, 
                        cdna_fasta, # in format:  Gene stable ID,Transcript stable ID,Gene type,Constitutive exon,Exon rank in transcript,Exon region start (bp),Exon region end (bp)
                        geneome_fasta,
                        exon_CSV,
                        cell_type = '',
                        expression_cutoff = 0.0001):
    
    """ extract and build a transcriptome from ensemble files 
    
    Build_Transcriptome(cell_expr_iso_dirt, cdna_fasta, geneome_fasta, exon_CSV,
                        cell_type, expression_cutoff = 0.0001)
        cell_expr_iso_dirt: dict of transcript ID: expression
        cdna_fasta: str file path 
        geneome_fasta: str file path # geneID, symbol, gene name, chromosome, start, stop, type?
        exon_CSV: str file path #Gene stable ID,Transcript stable ID,Strand,Transcript start (bp),Transcript end (bp),Exon region start (bp),Exon region end (bp)
        cell_type: str cell type
        expression_cutoff: cut off below which transcripts will be excluded
        
        returns transcriptome object
    
    """
    #cdna_fasta: >ENSG00000000419|DPM1|dolichyl-phosphate mannosyltransferase subunit 1, catalytic [Source:HGNC Symbol;Acc:HGNC:3005]|20|50934867|50958555|protein_coding
    # geneID, symbol, gene name, chromosome, start, stop, type?
    
    #exon_CSV: Gene stable ID,Transcript stable ID,Strand,Transcript start (bp),Transcript end (bp),Exon region start (bp),Exon region end (bp)
    
    
    if expression_cutoff < 0:
       print('expression_cutoff must be larger than or = to 0')
       return(print('expression_cutoff must be larger than or = to 0'))
   
    t = Transcriptome()
    t.cell_type = cell_type
    #'>ENST00000551848.1 cdna chromosome:GRCh38:12:32106835:32108796:1 gene:ENSG00000151746.13 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:BICD1 description:BICD cargo adaptor 1 [Source:HGNC Symbol;Acc:HGNC:1049]
    #>ENSMUST00000178537.2 cdna chromosome:GRCm39:6:41510135:41510146:1 gene:ENSMUSG00000095668.2 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:Trbd1 description:T cell receptor beta, D region 1 [Source:MGI Symbol;Acc:MGI:4439571]

    for ts in cell_expr_iso_dirt:
        t.TS[ts] = transcript()
    print('reading transcript data')
    RNA_file=open(cdna_fasta,'r')
    line = RNA_file.readlines() #all of RNA_file, read in FASTA from ensemble
    # go line by line of FASTA, extract header info or add onto sequence 
    fasta_TS = dict()
    fasta_G = dict()
    for item in line: # line is a list of lines 
        if item[0] == '>':# for each new FASTA entry
            trascript_ID=item[1:find_nth_n_str(item,'.',1)]
            
            fasta_TS[trascript_ID] = 0
            
            fasta_G[item[find_nth_n_str(item,':',6)+1:find_nth_n_str(item,'.',2)]] = 0
   
            exp = cell_expr_iso_dirt.get(trascript_ID,0)
            if exp < expression_cutoff:
                if t.TS.get(trascript_ID) == None:
                    ks = False
                    continue
                else:
                    del t.TS[trascript_ID]
                    ks = False
                    continue
            else:
                ks = True
                t.TS[trascript_ID].ID = trascript_ID
                gene_ID = item[find_nth_n_str(item,':',6)+1:find_nth_n_str(item,'.',2)]
                t.G[gene_ID] = gene()
                t.G[gene_ID].ID = gene_ID
                t.TS[trascript_ID].gene_ID = gene_ID
                t.TS[trascript_ID].exp = cell_expr_iso_dirt[trascript_ID]
                t.TS[trascript_ID].biotype  = item[find_nth_n_str(item,':',8)+1:find_nth_n_str(item,' ',6)]
                this_strand = item[find_nth_n_str(item,':',5)+1:find_nth_n_str(item,' ',3)]
                t.G[gene_ID].strand = this_strand 
                if this_strand == '1':
                    t.TS[trascript_ID].start = int(item[find_nth_n_str(item,':',3)+1:find_nth_n_str(item,':',4)])
                    t.TS[trascript_ID].stop = int(item[find_nth_n_str(item,':',4)+1:find_nth_n_str(item,':',5)])
                elif this_strand == '-1':
                    t.TS[trascript_ID].stop = int(item[find_nth_n_str(item,':',3)+1:find_nth_n_str(item,':',4)])
                    t.TS[trascript_ID].start = int(item[find_nth_n_str(item,':',4)+1:find_nth_n_str(item,':',5)])

                ch = item[find_nth_n_str(item,':',2)+1:find_nth_n_str(item,':',3)]
                t.G[gene_ID].chr = ch

        else:
            if ks == True:
                item=item.replace('\n', '')
                t.TS[trascript_ID].seq += item
            else:
                continue

    RNA_file.close()

    del_list = []
    for ts in t.TS:
        length = len(t.TS[ts].seq)
        if length == 0:
            del_list.append(ts)
    for ts in del_list:
        del t.TS[ts]
        
    print('reading exon data')
    exon_f=open(exon_CSV,'r')
    line = exon_f.readlines()
    print('building exons')
    #Gene stable ID,Transcript stable ID,Strand,Transcript start (bp),Transcript end (bp),Exon region start (bp),Exon region end (bp)

    exon_TS = dict()
    exon_G = dict()
    for l in line: # line is a list of lines
        item = l.replace('\n', '')
        gene_id = item[0:find_nth_n_str(item,',',1)]
        transcriptID = item[find_nth_n_str(item,',',1)+1:find_nth_n_str(item,',',2)]
        
        exon_TS[transcriptID]=0
        exon_G[gene_id]=0
        
        if t.G.get(gene_id) == None:
            continue
        else:
            transcriptID = item[find_nth_n_str(item,',',1)+1:find_nth_n_str(item,',',2)]
            if t.TS.get(transcriptID) == None:
                continue
            else:
                strand = item[find_nth_n_str(item,',',2)+1:find_nth_n_str(item,',',3)]
                startG = int(item[find_nth_n_str(item,',',5)+1:find_nth_n_str(item,',',6)])
                stopG = int(item[find_nth_n_str(item,',',6)+1:])
                if strand == '1':
                    t.TS[transcriptID].ex_exon_index_start.append(startG)
                    t.TS[transcriptID].ex_exon_index_stop.append(stopG)
                elif strand == '-1':
                    t.TS[transcriptID].ex_exon_index_start.append(stopG)
                    t.TS[transcriptID].ex_exon_index_stop.append(startG)
                
    exon_f.close()

    del_list = []
    
    for ts in t.TS:
        length1 = len(t.TS[ts].ex_exon_index_start)
        length2 = len(t.TS[ts].ex_exon_index_stop)
        if length1 == 0 or length2 == 0:
            del_list.append(ts)
        else:
            strand = t.G[t.TS[ts].gene_ID].strand
            if strand == '1':
                t.TS[ts].ex_exon_index_start.sort()
                t.TS[ts].ex_exon_index_stop.sort()
                t.TS[ts].in_exon_index_start = [t.TS[ts].ex_exon_index_start[0] - t.TS[ts].start]
                t.TS[ts].in_exon_index_stop = []
                for i, exonstart in enumerate(t.TS[ts].ex_exon_index_start):
                    exon_len = abs(t.TS[ts].ex_exon_index_stop[i] - exonstart)+ 1 # add 1 to keep indexing correct
                    t.TS[ts].in_exon_index_start.append(t.TS[ts].in_exon_index_start[i]+exon_len)
                    t.TS[ts].in_exon_index_stop.append(t.TS[ts].in_exon_index_start[i+1])
                    
                t.TS[ts].in_exon_index_start = t.TS[ts].in_exon_index_start[:-1]
                t.TS[ts].in_exon_index_start.sort()
                t.TS[ts].in_exon_index_stop.sort()
                
            elif strand =='-1':
                t.TS[ts].ex_exon_index_start.sort(reverse = True)
                t.TS[ts].ex_exon_index_stop.sort(reverse = True)
                t.TS[ts].in_exon_index_start = [t.TS[ts].start - t.TS[ts].ex_exon_index_start[0]]
                t.TS[ts].in_exon_index_stop = []
                for i, exonstart in enumerate(t.TS[ts].ex_exon_index_start):
                    exon_len =  abs(exonstart - t.TS[ts].ex_exon_index_stop[i])+1 # add 1 to keep indexing correct
                    t.TS[ts].in_exon_index_start.append(t.TS[ts].in_exon_index_start[i]+exon_len)
                    t.TS[ts].in_exon_index_stop.append(t.TS[ts].in_exon_index_start[i+1])
                    
                t.TS[ts].in_exon_index_start = t.TS[ts].in_exon_index_start[:-1]
                t.TS[ts].in_exon_index_start.sort()
                t.TS[ts].in_exon_index_stop.sort()

    for ts in del_list:
        del t.TS[ts]
      #ENSG00000000003|TSPAN6|tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:11858]|X|100627108|100639991|protein_coding
    geneome_file = open(geneome_fasta,'r')

      # go line by line of FASTA, extract header info or add onto sequence 
    ks = False
    print('reading genome file')
    lines = geneome_file.readlines()
    print('building geneome')
    geneome = dict()
    for item in lines: # line is a list of lines
          item = item.replace('\n', '')
          if item == '':
              continue
          elif item[0] == '>':# for each new FASTA entry
              G_IDt=item[1:find_nth_n_str(item,'|',1)]
              geneome[G_IDt] = 0
              if t.G.get(G_IDt) == None: 
                  ks = False
                  continue
              else:
                  ks = True
                  G_ID = G_IDt
                  t.G[G_ID].symbol = item[find_nth_n_str(item,'|',1)+1:find_nth_n_str(item,'|',2)]
                  t.G[G_ID].full_name = item[find_nth_n_str(item,'|',2)+1:find_nth_n_str(item,'|',3)]
                  t.G[G_ID].start = int(item[find_nth_n_str(item,'|',4)+1:find_nth_n_str(item,'|',5)])
                  t.G[G_ID].stop = int(item[find_nth_n_str(item,'|',5)+1:find_nth_n_str(item,'|',6)])
                  t.G[G_ID].biotype = item[find_nth_n_str(item,'|',6)+1:].replace('\n', '')
          else:
              if ks == True:
                  t.G[G_ID].seq += item
              else:
                  continue
    geneome_file.close()


    print('checking data')
    del_TS_list = []
    del_G_list = []
    for ts in t.TS:
          GID =  t.TS[ts].gene_ID
          t.TS[ts].len = len(t.TS[ts].seq)
          if t.G.get(GID) == None:
              del_TS_list.append(ts)
          else:
              t.G[GID].TS_ID.append(ts)
    for G in t.G:
          if len(t.G[G].TS_ID) == 0:
              del_G_list.append(G)
          else:
              t.Gsym2ID[t.G[G].symbol] = G
              for ts in t.G[G].TS_ID:
                  t.G[G].exp += t.TS[ts].exp
    for ts in del_TS_list:
          del t.TS[ts]
    for g in del_G_list:
          del t.G[g]
    
    commonTS = dict()
    noncommonTS = dict()
    
    for ts in exon_TS:
        if cell_expr_iso_dirt.get(ts) != None and fasta_TS.get(ts) != None:
            commonTS[ts] = 0
        else:
            noncommonTS[ts] = 0

    commonG = dict()
    noncommonG = dict()
    
    for g in geneome:
        if fasta_G.get(g) != None and exon_G.get(g)!= None:
            commonG[g] = 0
        else:
            noncommonG[g]=0

    #print('cell_expr_iso_dirt = '+str(len(cell_expr_iso_dirt)))
    #print('exon_TS = '+str(len(exon_TS)))
    #print('fasta_TS = '+str(len(fasta_TS)))
    
    #print('geneome = '+str(len(geneome)))
    #print('fasta_G = '+str(len(fasta_G)))
    #print('exon_G = '+str(len(exon_G)))
    
    #print('common TS = '+str(len(commonTS)))
    #print('uncommon TS = '+str(len(noncommonTS)))
    #print('common G = '+str(len(commonG)))
    #print('uncommon G = '+str(len(noncommonG)))
    #gc.collect()

    print(t.info())
    return(t)
    
# =============================================================================

def get_size(obj, seen=None):
    """Recursively finds size of objects
    
    get_size(obj, seen=None)
        obj: object to calulate size of 
        seen: ?
        returns size of object in bytes
    
    """
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size

# =============================================================================

def save_pickle(obj,
                file_path,
                file_name):
    """ saves objects as pickle files, protocol=4
    
    save_pickle(obj, file_path,file_name)
        obj: object to be pickled 
        file_path: str path to save location 
        file_name: str name of file to be saved
        returns None if file_path does not exist
    
    """
     ################### check inputs #########################  
    file_path2check = os.path.normpath(file_path)
    if type(file_path) != str:
        raise TypeError('file_path must be a str not a '+str(type(file_path)))
   
    if os.path.exists(file_path2check) == False:
        raise Exception(file_path2check+' does not exist')
    
    if type(file_name) != str:
        raise TypeError('file_name must be a str not a '+str(type(file_name)))
     #########################################################  
    
    path = os.path.normpath(os.path.join(file_path,file_name+'.pk'))
    pkl = open(path,'wb')
    pickle.dump(obj,pkl,protocol=4)
    pkl.close()
    
# =============================================================================    

def load_pickle(file_path):
    """ loads objects from pickle files
    
    load_pickle(file_path)
        file_path: path of pickled file to load
        retuns object or None if file_path does not exist
        
    """
     ################### check inputs #########################
    if type(file_path) != str:
        raise TypeError('file_path must be a str not a '+str(type(file_path)))
   
    file_path = os.path.normpath(file_path)
    if os.path.exists(file_path) == False:
        raise Exception(file_path+' does not exist')
    ##########################################################
        
    else:
        pkl = open(file_path,'rb')
        obj = pickle.load(pkl)
        pkl.close()
    return(obj)
  
# =============================================================================
 
def save_ISO_or_TDP_dict(A_dict,
                         file_path):
    """ recursively saves items in a dictionary as pickle file
    
    save_ISO_or_TDP_dict(A_dict, file_path)
        A_dict: dict in form key1:dict(key2:value)
        file_path: str file path to save separated dict in pickle form 

    """
     ################### check inputs #########################
    if type(file_path) != str:
        raise TypeError('file_path must be a str not a '+str(type(file_path)))
    
    file_path = os.path.normpath(file_path)
    if os.path.exists(file_path) == False:
        raise Exception(file_path+' does not exist')
    
    if type(A_dict) != dict:
        raise TypeError('A_dict must be a dict not a '+str(type(A_dict)))
    ##########################################################
        
    for G_ID in A_dict:
        iso_info = A_dict[G_ID]
        save_pickle(iso_info , file_path , G_ID)

# =============================================================================

def get_TDP_CG(transcriptome,
               file_path,
               verbose = False):
    """ calulates local entropy, enthalpy, and CG content of all transcripts in a transcriptome 
    
    get_TDP_CG(transcriptome, file_path, verbose = False)
        transcriptome: transcriptome object to calulate thermodynamic properties for
        file_path: str file path to save thermodynamic properties
        verbose: bool print progress of calulation every 20s 
        returns path to thermodynamic properties database
    
    """
    ################### check inputs #########################
    if type(file_path) != str:
        raise TypeError('file_path must be a str not a '+str(type(file_path)))
    
    file_path = os.path.normpath(file_path)
    if os.path.exists(file_path) == False:
        raise Exception(file_path+' does not exist')
    
    #if isinstance(transcriptome,Transcriptome) != True:
    #    raise TypeError('transcriptome must be a transcriptome not a '+str(type(transcriptome)))
   
    if type(verbose) != bool:
        raise TypeError('verbose must be a True or False not a '+str(type(verbose)))
    ##########################################################
        
    # make local names for Transcriptome data
    RNA_list = []
    trascript_ID = []
    gene_ID = []
    TS_length = []
    
    for Ts in transcriptome.TS:
        Ts = transcriptome.TS[Ts]
        RNA_list.append(Ts.seq)
        trascript_ID.append(Ts.ID)
        gene_ID.append(Ts.gene_ID)
        TS_length.append(len(Ts.seq))
    
    for G in transcriptome.G:
        G = transcriptome.G[G]
        RNA_list.append(G.seq)
        trascript_ID.append('NSG'+G.ID)
        gene_ID.append(G.ID)
        TS_length.append(len(G.seq))

    RNA_list_len = len(trascript_ID)    
    TDP_CG=dict() # make dictionary of gene IDs
    # dictionary of thermodynamic properties from SantaLucia J Jr, Hicks D (2004) The 
    # thermodynamics of DNA structural motifs. Annu Rev Biophys Biomol
    # Struct 33(1):415–440
    TDP_d = {'AA':[-7.6,-21.3], 
         'TT':[-7.6,-21.3],
         'AT':[-7.2,-20.4],
         'TA':[-7.2,-21.3],
         'CA':[-8.5,-22.7],
         'TG':[-8.5,-22.7],
         'AC':[-8.4,-22.4],
         'GT':[-8.4,-22.4],
         'CT':[-7.8,-21.0],
         'AG':[-7.8,-21.0],
         'GA':[-8.2,-22.2],
         'TC':[-8.2,-22.2],
         'CG':[-10.6,-27.2],
         'GC':[-9.8,-24.4],
         'GG':[-8.0,-19.9],
         'CC':[-8.0,-19.9],}
    # from merFISH MATlab code https://github.com/ZhuangLab/MERFISH_analysis/blob/master/probe_construction/TRDesigner.m
    # H = [-7.6 -8.4 -7.8 -7.2 -8.5 -8.0 -10.6 -7.8 ...   % AA/TT GT/CA CT/GA AT/TA CA/GT GG/CC CG/GC CT/GA
    #         -8.2 -9.8 -8.0 -8.4 -7.2 -8.2 -8.5 -7.6];       % GA/CT GC/CG GG/CC GT/CA TA/AT GA/CT CA/GT AA/TT
    #     S = [-21.3 -22.4 -21.0 -20.4 -22.7 -19.9 -27.2 -21.0 ...    % AA/TT GT/CA CT/GA AT/TA CA/GT GG/CC CG/GC CT/GA
    #         -22.2 -24.4 -19.9 -22.4 -21.3 -22.2 -22.7 -21.3];       % GA/CT GC/CG GG/CC GT/CA TA/AT GA/CT CA/GT AA/TT
    
    c=0 # counter
    t1=0
    cg_dict = {'C':1,'G':1}
    # loop over transcripts 
    for RNA_index, RNA in enumerate(RNA_list):
        G_ID = gene_ID[RNA_index] # get IDs and data for each transcript
        TS_ID = trascript_ID[RNA_index]
        length = TS_length[RNA_index]
        # initiate list of enthalpy, entrop, and CG  
        TDP_H=[]
        TDP_S=[]
        CG=[]
        
        # verbose print progress 
        t = time.time()
        if verbose == True and t > t1:
            print('thermodynamic properties calculated for ' + str(c)+'/'+str(RNA_list_len) + ' sequences')
            t1 = t + 20
            
        # if gene ID not in dictionary, add it        
        c+=1
        if TDP_CG.get(G_ID) == None: 
            TDP_CG[G_ID]=dict()
        else:
            pass

        # for each base, add to CG list, H, and S 
        for b, base in enumerate(RNA):
            CG.append(cg_dict.get(base,0)) # add 1 to list if C or G
            if b == length-1: 
                pass
            else:
                part = RNA[b:b+2]
                TDP_H.append(TDP_d.get(part,[0,0])[0])
                TDP_S.append(TDP_d.get(part,[0,0])[1])
        TDPCG=[TDP_H,TDP_S,CG]
        TDP_CG[G_ID][TS_ID] = TDPCG #TDP_CG=[TDP_H,TDP_S,CG]
    print('saveing thermodynamic properties')
    #save it, make db file if not already their 
    db_path = os.path.join(file_path,'db')
    db_path = os.path.normpath(db_path)
    if os.path.exists(db_path) == False:
        os.makedirs(db_path)
    number = 0
    while True:
        TDP_CG_path = os.path.normpath(os.path.join(db_path, 'TDP_CG_' +str(number)))
        if os.path.exists(TDP_CG_path) == True:
            number +=1
        else:
            break
    TDP_CG_path = os.path.normpath(os.path.join(db_path, 'TDP_CG_' +str(number)))
    os.makedirs(TDP_CG_path)
    save_ISO_or_TDP_dict(TDP_CG,TDP_CG_path)
    transcriptome.TDP_db_path = TDP_CG_path
    del TDP_CG
    gc.collect()
    return(TDP_CG_path)

# =============================================================================

def DNA2int(DNA, N = 'A'):
    ''' Converts DNA into intiger form.
    
        DNA: str DNA with bases A,T,G,C,a,t,g,c 
        N: 'random', None, or bases A,T,G,C,a,t,g,c 
        returns an intager repersentation of DNA, if DNA has base not in above 
        list, return 0
        
    '''
    ################### check inputs #########################
    if type(DNA) != str:
        raise TypeError('DNA must be a str not a '+str(type(DNA)))
    if type(N) != type(None) and N != 'random' and type(N) != str :
        raise TypeError('N must be a bool, \'random\', or a single charature in \'ACTG\' not '+str(N))
    if type(N) == str and N!= 'random' and len(N) != 1:
        raise Exception('if N given as str, length must be 1')
    if type(N) == str and N!= 'random' and N not in 'ACTGactg':
        raise Exception('if N, must be a DNA base in \'ACTGactg\'')
        
    ##########################################################
        
    d2b = {'A':'00','C':'01','T':'10','G':'11','a':'00','c':'01','t':'10','g':'11'}
    bs = '1'
    # binary 000000 -> 0 not AAAA so all codes have to start with a 1 so no ambugiuity
    for base in DNA:
        if d2b.get(base) == None:
            if N != None and base in 'Nn':
                if N == 'random':
                    base = random.choice(['A','T','C','G'])
                elif type(N) == str:
                    base = N
            else:
                raise Exception('base <'+base+'> not in libuary (A,T,G,C,a,t,g,c)') 
        bs+=d2b[base]
    intager = int(bs, 2)
    return(intager)

# =============================================================================

def int2DNA(intager):
    ''' Converts intager into DNA form.
    
        DNA: int DNA 
        returns an DNA str with bases A,T,G,C or None if intager is 0
        
    '''
    ################### check inputs ######################### 
    if type(intager) != int:
        raise TypeError('intager must be a int not a '+str(type(intager)))
    ##########################################################
        
    else:
        b2d = {'00':'A','01':'C','10':'T','11':'G'}
        # 3 because the bin finction gives '0b1(the code for seq)' the leading 1 is because leading 0 are ambiguous 
        bi = bin(intager)[3:]
        seq = ''
        for sec in range(int(len(bi)/2)):
            loc = sec*2
            seq += b2d[bi[loc:loc+2]]
        return(seq)
    
# =============================================================================

def Make_OT_dict_Full(transcriptome,
                      file_path,
                      OT_dict_file_name = 'OT_dict',
                      OT_length = 17,
                      verbose = False,
                      frag = 10,
                      threads = -1):
    """ builds a dict of off target sequences with associated penalty for transcriptome 
    
    Make_OT_dict_Full(transcriptome, file_path, OT_dict_file_name = 'OT_dict',
                      OT_length = 17, verbose = False, frag = 10)
        transcriptome: transcriptome objuct to calulate off target values for
        file_path: str path of file location to store off target calulations
        OT_dict_file_name: str name of file to make and store off target calulations
        OT_length: str length of off target subsections to calulate off target score for
        verbose: bool print progress of calulation every 20s 
        frag: int number of fragments to break up transcriptome for calulations 
        threads: number of threads to use in calulation, default is all cores if -1
        returns file path as str for OT_dict_files

    """
    ################### check inputs #########################
    if type(frag) != int:
        raise TypeError('frag must be of type int not '+ str(type(frag)))    
    
    # check value of OT_length
    if type(OT_length) != int:
        raise TypeError('OT_length must be a int not a '+str(type(OT_length)))
    
    if OT_length <= 0:
        raise Exception('OT_length must be an intiger grater than 0')
   
    # check that OT_dict_file_name is a string 
    if type(OT_dict_file_name) != str:
        raise TypeError('OT_dict_file_name must be of type str not '+ str(type(OT_dict_file_name)))
   
    # check that OT_dict_file_name exist
    if os.path.exists(file_path) == False:
        raise Exception('file_path '+file_path+' does not exist')
    
    # check that Transcriptome is a Transcriptome object 
    #if isinstance(transcriptome,Transcriptome) != True:
    #    raise TypeError('transcriptome must be of type Transcriptome not '+ str(type(transcriptome)))
    ##########################################################
        
    else:
        RNA_list = [] # get trascripts seqs
        exp = []  # get expression of trascripts
        RNA_len = []
        extend = OT_length -1
    
        for Ts in transcriptome.TS:
            Ts = transcriptome.TS[Ts]
            RNA_list.append(Ts.seq)
            exp.append(Ts.exp)
            RNA_len.append(len(Ts.seq))
            gene = transcriptome.G[Ts.gene_ID]
            strand = gene.strand
            
            # get extended intron sequences from gene seq
            exon_start = Ts.ex_exon_index_start
            exon_stop = Ts.ex_exon_index_stop
            extended_intorn_start_stop = [] # list of local index on gene sequence [(start1, stop1),(start2, stop2),...]
            gene_len = len(gene.seq)
            
            if strand == '1' or strand != '-1':
              gene_start = gene.start 
            else:
              gene_start = gene.stop
            
            # get location(index) of extended introns 
            for i, intron_start in enumerate(exon_stop[:-1]):
              inton_stop = exon_start[i+1]
              extended_intron_start = max(0,abs(intron_start-gene_start)-extend) # index on local gene seq
              extended_intron_stop = min(gene_len,abs(inton_stop-gene_start)+extend) # index on local gene seq
              extended_intorn_start_stop.append((extended_intron_start,extended_intron_stop))
             
             # get seq of extended introns and add to list
            for indexs in extended_intorn_start_stop:
              start =  indexs[0]
              stop = indexs[1]
              extended_intron = gene.seq[start:stop]
              RNA_list.append(extended_intron)
              exp.append(Ts.exp)
              RNA_len.append(len(extended_intron))
        RNAl = len(RNA_list)
        
        # make a db file to store all the calulated values
        db_path = os.path.normpath(os.path.join(file_path,'db'))
        if os.path.exists(db_path) == False:
            os.makedirs(db_path)
        number = 0
        while True:
            OT_dict_path = os.path.normpath(os.path.join(db_path, OT_dict_file_name +'_' + str(number)))
            if os.path.exists(OT_dict_path) == True:
                number +=1
            else:
                break
        OT_dict_path = os.path.normpath(os.path.join(db_path, OT_dict_file_name +'_' + str(number)))
        os.makedirs(OT_dict_path)

        c =0
        t0 = time.time()
        t1=0
        print(len(RNA_list))
        print(time.time())

        OT_dict = dict()
        OT_dict_len =0
        dict_num = 0
        
        #frag data
        chunk_size = round(len(RNA_list)/frag)
        index_start = []
        for x in range(frag):
            index_start.append(x*chunk_size)
        index_stop = index_start[1:]
        index_stop.append(len(RNA_list))
        
        chunk_list_RNA = []
        chunk_list_exp = []
        RNA_len_list = []
        dict_list_num = []
        OT_dict_path_list = []
        for s, start in enumerate(index_start):
            chunk_list_RNA.append(RNA_list[start:index_stop[s]]) 
            chunk_list_exp.append(exp[start:index_stop[s]]) 
            RNA_len_list.append(RNA_len[start:index_stop[s]])
            dict_list_num.append(s)
            OT_dict_path_list.append(OT_dict_path)
        print(dict_list_num)
        
        # function to calulate off target on fraged data 
        def make_OT(RNA_list,exp,RNA_len,dict_list_num,OT_dict_path):
           print(dict_list_num)
           OT_dict = dict()
           t1 = time.time()
           len_RNA = len(RNA_list)
           for RNA_index, RNA in enumerate(RNA_list):
               t = time.time()
               #if t>t1:
               #    print(str(RNA_index)+'/'+str(len_RNA)+' - '+str(dict_list_num))
               #    t1 += 30 
               exp_ = exp[RNA_index]
               RNA_len_ = RNA_len[RNA_index]
               s = RNA_len_-OT_length + 1
               for b in range(s):
                   subseq = DNA2int(RNA[b:b+OT_length], N = 'random')
                   if OT_dict.get(subseq)==None:
                       OT_dict[subseq] = exp_
                   else:
                       OT_dict[subseq] += exp_
           print('saveing_'+str(dict_list_num))
           save_pickle(OT_dict,OT_dict_path,'OT_dict_'+str(dict_list_num))
           
        # set up threads/cores 
        if threads == -1: 
            threads = multiprocessing.cpu_count()
        pool = ThreadPool(threads)
        
        # run calulations on mulitiple cores 
        results = pool.starmap(make_OT, zip(chunk_list_RNA,chunk_list_exp,RNA_len_list,dict_list_num,OT_dict_path_list),1)#
        pool.close() 
        pool.join()
    
    OT_dict['OT_length'] = OT_length
    
    print('saveing OT_dict')
    save_pickle(OT_dict, OT_dict_path , 'OT_length' ) # save as pickle file
    print()
    print(OT_dict_path)
    gc.collect()
    return(OT_dict_path)

# =============================================================================
        
def get_OT_dict_Iso(transcriptome,
                    OT_dict_dict,
                    file_path,
                    verbose = False,
                    frag = 10,
                    threads = -1):
    """ builds a dict of off tager sequences with associated penalty for isoforms 
    
    get_OT_dict_Iso(transcriptome, OT_dict_dict, file_path, verbose = False,
                    frag = 10)
        transcriptome: transcriptome object to calulate off target values for
        OT_dict_dict: str file path to calulated off target values from Make_OT_dict_Full
        file_path: str path of file location to store off target calulations
        verbose: bool print progress of calulation every 20s 
        frag: int number of fragments to break up transcriptome for calulations
        returns str file path as str for ISO_dict_path_files
    
    """
    ################### check inputs #########################  
    if type(frag) != int:
        raise TypeError('frag must be of type int not '+ str(type(frag)))  
    
    # if isinstance(transcriptome,Transcriptome) != True:
    #     raise TypeError('transcriptome must be of form Transcriptome not '+ str(type(transcriptome)))

    if type(OT_dict_dict) != str:
        raise TypeError('OT_dict_dict must be of form str not '+ str(type(OT_dict_dict)))

    if type(file_path) != str:
        raise TypeError('file_path must be of form str not '+ str(type(file_path)))

    if os.path.exists(os.path.normpath(OT_dict_dict)) == False:
        raise Exception('OT_dict_dict '+os.path.normpath(OT_dict_dict)+' does not exist')

    if os.path.exists(os.path.normpath(file_path)) == False:
        raise Exception('file_path '+os.path.normpath(file_path)+' does not exist')
    ##########################################################

    # get list of file names of offtarget dicts 
    OT_dict_dict = os.path.normpath(OT_dict_dict)
    db_list = [f for f in os.listdir(OT_dict_dict) if not f.startswith('.') and not f.startswith('OT_length')]
    db_list1 = []
    for x in range(frag):
        db_list1.append([db_list[x]])
    db_list = db_list1
    OT_length_db = load_pickle(os.path.normpath(os.path.join(OT_dict_dict,'OT_length.pk')))
    OT_length = OT_length_db['OT_length']
    del OT_length_db

    RNA_list = [] 
    trascript_ID = []
    gene_ID = []
    exp = []
    
    extend=OT_length
    for G in transcriptome.G:
        G = transcriptome.G[G]
        for Ts in G.TS_ID:    
            Ts = transcriptome.TS[Ts]
            RNA_list.append(Ts.seq)
            trascript_ID.append(Ts.ID)
            gene_ID.append(Ts.gene_ID)
            exp.append(Ts.exp)
            gene = transcriptome.G[Ts.gene_ID]
            exon_start = Ts.ex_exon_index_start
            exon_stop = Ts.ex_exon_index_stop
            extended_intorn_start_stop = [] # list of local index on gene sequence [(start1, stop1),(start2, stop2),...]
            gene_len = len(gene.seq)
            strand = gene.strand
            if strand == '1' or strand != '-1':
                gene_start = gene.start 
            else:
                gene_start = gene.stop
                
            # get location(index) of extended introns
            for i, intron_start in enumerate(exon_stop[:-1]):
                inton_stop = exon_start[i+1]
                extended_intron_start = max(0,abs(intron_start-gene_start)-extend) # index on local gene seq
                extended_intron_stop = min(gene_len,abs(inton_stop-gene_start)+extend) # index on local gene seq
                extended_intorn_start_stop.append((extended_intron_start,extended_intron_stop))
                
            # get extended intron sequences from gene seq
            for indexs in extended_intorn_start_stop:
                start =  indexs[0]
                stop = indexs[1]
                extended_intron = gene.seq[start:stop]
                RNA_list.append(extended_intron)
                trascript_ID.append(Ts.ID)
                gene_ID.append(gene.ID)
                exp.append(Ts.exp)
    
    if verbose == True:            
        print('building all')
        
    RNA_len = len(RNA_list)
    #save it, make db file if not already their 
    db_path = os.path.normpath(os.path.join(file_path,'db'))
    if os.path.exists(db_path) == False:
        os.makedirs(db_path)
    number = 0
    while True:
        ISO_dict_path = os.path.normpath(os.path.join(db_path,'OT_dict_iso' + str(number)))
        if os.path.exists(ISO_dict_path) == True:
            number +=1
        else:
            break
    
    ISO_dict_path = os.path.normpath(os.path.join(db_path,'OT_dict_iso' + str(number)))
    os.makedirs(ISO_dict_path)
    
    def build_all(RNA_list,exp,gene_ID_,OT_length,ISO_dict_path):
        count = 0
        ISO_dict=dict()# made dictanary for all gene IDs
        ISO_dict['all']= dict()
        iso_local = ISO_dict
        RNA_len_last_index = RNA_len-1
        for RNA_index, RNA in enumerate(RNA_list): # loop over all transcripts
            count += 1
            exp_of_TS = exp[RNA_index] # get expression data 
            G_ID = gene_ID_ #[RNA_index]# get gene ID data
            s = len(RNA)-OT_length+1 # number of seqs with OT_length in transcript
            for b in range(s): #loop over trascript for all s 
                seq = RNA[b:b+(OT_length)]
                seq = DNA2int(seq, N = 'A')
                if iso_local['all'].get(seq) == None:
                    iso_local['all'][seq] = exp_of_TS
                else:
                    iso_local['all'][seq] += exp_of_TS
                
            # if RNA_index == RNA_len_last_index or G_ID != gene_ID[RNA_index+1]:
            #     save_pickle(iso_local,ISO_dict_path,G_ID)
            #     ISO_dict=dict()# made dictanary for all gene IDs
            #     ISO_dict['all']= dict()
            #     iso_local = ISO_dict
        save_pickle(iso_local,ISO_dict_path,gene_ID_)
        return(count)
    
    # frag data
    #RNA_list,exp,gene_ID,OT_length,ISO_dict_path
    #RNA_list_ = []
    #exp_ = []
    #gene_ID_ = []
    #OT_length_ = []
    #ISO_dict_path_ = []
    
    #last_GID = gene_ID[0]
    #start_GID_I = 0
    #gene_ID_.append(last_GID)
    
    
    #for G , G_ID in enumerate(gene_ID):
        #if G_ID != last_GID and G != len(gene_ID)-1:
        #    RNA_list_.append(RNA_list[start_GID_I:G])
        #    exp_.append(exp[start_GID_I:G])
        #    gene_ID_.append(G_ID)
        #    OT_length_.append(OT_length)
        #    ISO_dict_path_.append(ISO_dict_path)
        #    last_GID = G_ID
        #    start_GID_I = G
        #elif G == len(gene_ID)-1:
            #print('####################################')
        #    RNA_list_.append(RNA_list[start_GID_I:])
        #    exp_.append(exp[start_GID_I:])
        #    gene_ID_.append(G_ID)
        #    OT_length_.append(OT_length)
        #    ISO_dict_path_.append(ISO_dict_path)
        #else:
            #continue
    RNA_list_ = [] #list of list 
    exp_ = []
    gene_ID_ = []
    OT_length_ = []
    ISO_dict_path_ = []
    
    # build DF 
    frag_DF = pd.DataFrame(data = {'gene_ID':gene_ID,'RNA_list':RNA_list,'exp':exp},
                           columns=['gene_ID','RNA_list','exp'])
    gene_IDUnique = list(set(gene_ID))
    
    Uc = 0
    for G_ID in gene_IDUnique:
        Uc += 1
        if Uc%100 ==0:
            print(Uc)
            
        RNA_list_.append(frag_DF[frag_DF['gene_ID'] == G_ID]['RNA_list'].tolist())
        exp_.append(frag_DF[frag_DF['gene_ID'] == G_ID]['exp'].tolist())
        gene_ID_.append(G_ID)
        OT_length_.append(OT_length)
        ISO_dict_path_.append(ISO_dict_path)
        
    #print(RNA_list_[0:5])
    #print(exp_[0:5])
    #print(len(gene_ID_))
    #print(gene_ID_)
    #RNA_list,exp,gene_ID,OT_length,ISO_dict_path            
                
    # set up threads/cores 
    if threads == -1: 
        threads = multiprocessing.cpu_count()
    else:
        threads = min(multiprocessing.cpu_count(),threads)
    
    if threads > 1:
        pool = ThreadPool(threads)
        # run calulations on mulitiple cores 
        results = pool.starmap(build_all, zip(RNA_list_,exp_,gene_ID_,OT_length_,ISO_dict_path_),1)#
        pool.close() 
        pool.join()
        #print(gene_ID)
    else:
        #print('###################################')
        #print(gene_ID_)
        ziped = zip(RNA_list_,exp_,gene_ID_,OT_length_,ISO_dict_path_)
        
        for A,B,C,D,E in ziped:
            BA = build_all(A,B,C,D,E)

    if verbose == True:
        print('building iso and gene spec')           
    # calulate iso index and gene index for each OT seq in each trascript
    RNA_len_last_index = len(RNA_list)-1
    for set_ in db_list:
        db_s = []
        for db in set_:
            path = os.path.normpath(os.path.join(OT_dict_dict,db))
            OT_dict = load_pickle(path)
            db_s.append(OT_dict)
            if verbose == True:
                print(path+' loaded')
        iso_local = load_pickle(os.path.normpath(os.path.join(ISO_dict_path,(gene_ID[0]+'.pk'))))
        #print(gene_ID[0])

        for RNA_index, RNA in enumerate(RNA_list):
                exp_of_TS = exp[RNA_index]
                TS_ID = trascript_ID[RNA_index]# get transcript ID dataTS_ID = trascript_ID[RNA_index]# get transcript ID data
                #print(TS_ID)
                G_ID = gene_ID[RNA_index]
                s = len(RNA)-OT_length+1 # number of seqs with OT_length in transcript
                if iso_local.get(TS_ID) == None:
                    iso_local[TS_ID]=dict() # make dict for each transcript
                for b in range(s): #loop over trascript for all s 
                    seq = RNA[b:b+OT_length] # sliding window 
                    seq = DNA2int(seq, N = 'A')
                    iso_pen = iso_local['all'][seq] # get its iso penelty 
                    
                    if iso_local[TS_ID].get(seq) == None:
                        TST_pen = 0
                        for db in db_s:    
                            TST_pen += db.get(seq,0)
                        iso_index = exp_of_TS/iso_pen
                        gene_index = TST_pen
                        iso_local[TS_ID][seq] = [iso_index,gene_index] # add to dictinary
                    else:
                        TST_pen = 0
                        for db in db_s:    
                            TST_pen += db.get(seq,0)
                            
                        iso_local[TS_ID][seq][1] += TST_pen
                        
                if RNA_index == RNA_len_last_index or  G_ID != gene_ID[RNA_index+1]:  
                    save_pickle(iso_local,ISO_dict_path,G_ID)
                    if RNA_index != RNA_len_last_index:
                        iso_local = load_pickle(os.path.normpath(os.path.join(ISO_dict_path,(gene_ID[RNA_index+1]+'.pk'))))
                        #print(gene_ID[RNA_index+1])

    if verbose == True:
        print('updating gene index')
    def update_gene_index(ISO_dict_path,G):    
    #for G in transcriptome.G:
        iso_local = load_pickle(os.path.normpath(os.path.join(ISO_dict_path,G+'.pk')))
        for TS in iso_local:
            if TS =='all':
                continue
            else:
                for seq in iso_local[TS]:
                    iso_pen = iso_local['all'][seq]
                    TST_pen = iso_local[TS][seq][1]
                    if TST_pen == 0:
                        gene_index = 1
                    else:
                        gene_index = iso_pen/TST_pen
                    iso_local[TS][seq][1] = gene_index
        save_pickle(iso_local,ISO_dict_path,G)
    
    pool = ThreadPool(threads)
    # run calulations on mulitiple cores 
    G_ = list(transcriptome.G.keys())
    ISO_dict_path_ = [ISO_dict_path]*len(G_)
    results = pool.starmap(update_gene_index, zip(ISO_dict_path_, G_),1)#
    pool.close() 
    pool.join()
    
    OT = {'OT_length': OT_length}
    save_pickle(OT,ISO_dict_path,'OT_length')

    transcriptome.OT_db_path = ISO_dict_path
    gc.collect()

    return(ISO_dict_path)

# =============================================================================

def List_intersection(lst1, lst2):
    # Use of hybrid method
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3

# =============================================================================

def Subset_Transcriptome(Old_Transcriptome, 
                         subset_list = None, # takes gene id. TS id, biotype, gene symbol
                         index = [0,0],
                         verbose = False):
    """ subsets transcriptome object 
    
    Subset_Transcriptome(Old_Transcriptome, subset_list = None, index = [0,0])
        Old_Transcriptome: Transcriptome object to be subsetted
        subset_list: list of transcript IDs, gene IDs, biotypes, and/or 
        gene symbols to subset on
        index: list of two numbers of indicies of transcripts to subset 
        returns transcriptome object subseted 
        
        prints things in subset_list not found in Old_Transcriptome
    
    """
    ################### check inputs #########################
    # if isinstance(Old_Transcriptome,Transcriptome) != True:
    #     raise TypeError('Old_Transcriptome must be of form Transcriptome not '+ str(type(Old_Transcriptome)))
    
    if type(subset_list) != list and type(subset_list) != type(None):
        raise TypeError('subset_list must be of form str or None not '+ str(type(subset_list)))
    
    if type(index) != list:
        raise TypeError('index must be of form list in form [start index, stop index] not '+ str(index))
    ##########################################################
        
    t = Transcriptome() # new trascriptome
    
    Ts_list_ind = []
    gene_list_ind = []
    could_not_find = []
    gene_sym_dict = dict()
    TS_biotype_dict = dict()
    if subset_list == None:
        pass
    else:
        for gene in Old_Transcriptome.G: # get symbol dict
            gene_sym_dict[Old_Transcriptome.G[gene].symbol] = gene

        for TS in Old_Transcriptome.TS: # make biotype dict
            if TS_biotype_dict.get(Old_Transcriptome.TS[TS].biotype) == None:
                TS_biotype_dict[Old_Transcriptome.TS[TS].biotype] = []
            TS_biotype_dict[Old_Transcriptome.TS[TS].biotype].append(TS)

        for item in subset_list:
            if Old_Transcriptome.G.get(item) != None: # if gene IDs given 
                gene_list_ind.append(item)
                ind = Old_Transcriptome.G[item].TS_ID
                for n in ind:
                    Ts_list_ind.append(n)
                    
            elif gene_sym_dict.get(item) != None: # if gene symbol given 
                gid = gene_sym_dict[item]
                gene_list_ind.append(gid)
                ind = Old_Transcriptome.G[gid].TS_ID
                for n in ind:
                    Ts_list_ind.append(n)
            
            elif Old_Transcriptome.TS.get(item) != None: # if TS id given
                gene_list_ind.append(Old_Transcriptome.TS[item].gene_ID)
                Ts_list_ind.append(item)
                
            elif TS_biotype_dict.get(item) != None: # if biotype given 
                TS_ids = TS_biotype_dict[item]
                for n in TS_ids:
                    Ts_list_ind.append(n)
                    gene_list_ind.append(Old_Transcriptome.TS[n].gene_ID)
            else:
                could_not_find.append(item)

    start = index[0]
    end = index[1]
    gid_list  = list(Old_Transcriptome.G.keys()) # indexing on gene list 
    for ind in range(start, end):
        gene_list_ind.append(gid_list[ind])
        for TS in Old_Transcriptome.G[gid_list[ind]].TS_ID:
            Ts_list_ind.append(TS)
    
    # make sure all transcripts have their included genes
    for tsID in Ts_list_ind:
        gene_list_ind.append(Old_Transcriptome.TS[tsID].gene_ID)
    
    # get only unique genes and TSids 
    Ts_list_ind = list(set(Ts_list_ind))
    gene_list_ind = list(set(gene_list_ind))
    
    
    # add gene() to new transcriptome, update TS_ID in each gene
    for gID in gene_list_ind:
        tempGene = Old_Transcriptome.G[gID]
        #clear TS_ID list 
        tempGene.TS_ID = []
        t.G[gID] = tempGene
        t.IDs[gID] = []
        
    #add transcript() to new transcriptome
    for tsID in Ts_list_ind:
        # add TS ids to TS_ID for gene
        Geen_IDforTS = Old_Transcriptome.TS[tsID].gene_ID
        t.G[Geen_IDforTS].TS_ID.append(tsID)

        t.TS[tsID] = Old_Transcriptome.TS[tsID]
        t.IDs[tsID] = []
    
    
    #inharit DB paths and cell type
    t.TDP_db_path = Old_Transcriptome.TDP_db_path
    t.OT_db_path = Old_Transcriptome.OT_db_path
    t.cell_type = Old_Transcriptome.cell_type
    
    # build IDs and Gsym2ID dicts 
    for GID in gene_list_ind:
        t.IDs[GID] = []
        t.Gsym2ID[GID] = t.G[GID].symbol
    for TID in Ts_list_ind:
        t.IDs[TID] = []
    
    if verbose == True:
        if len(could_not_find)>0:
            print('could not find')
            for n in could_not_find:
                print(n)
    else:
        pass
    #t = copy.deepcopy(t)
    gc.collect()
    return(t)

# =============================================================================

def append_transcriptome(Transcriptome1,
                         TS = None, # list of transcript objects named with TS ids
                         Gene = None, # list of gene objects names with GIds
                         Transcriptome2 = None):
    """ adds to transcriptome object 
    
    append(TS = None,Gene = None,Transcriptome2 = None)
        TS: list of transcript objects to append to transcriptome 
        Gene: list of gene objects to append to transcriptome
        Transcriptome2: Transcriptome object to append to transcriptome
        returns appended transcriptome
        
    """
    ################### check inputs #########################
    if isinstance(Transcriptome1,Transcriptome) != True:
        raise TypeError('Transcriptome1 must be of form Transcriptome not '+ str(type(Transcriptome1)))
    
    if isinstance(Transcriptome2,Transcriptome) != True and type(Transcriptome2) != type(None):
        raise TypeError('Transcriptome2 must be of form Transcriptome or None not '+ str(type(Transcriptome2)))

    if type(TS) != list and type(TS) != type(None):
        raise TypeError('TS must be of form list or None not '+ str(type(TS)))

    if type(Gene) != list and type(Gene) != type(None):
        raise TypeError('Gene must be of form list or None not '+ str(type(Gene)))
    ##########################################################

    #T3 = copy.deepcopy(Transcriptome1)
    T3 = Transcriptome1
    #Transcriptome2 = copy.deepcopy(Transcriptome2)
    
    if Transcriptome2 == None:
        pass
    else:
        for Ts in Transcriptome2.TS:
            if T3.TS.get(Ts) != None:
                print( Ts + ' already exist in transcriptome, now  ' + Ts + '.1')
                Transcriptome2.TS[Ts].ID += '.1'
            T3.TS[Transcriptome2.TS[Ts].ID] = Transcriptome2.TS[Ts]
        for gene in Transcriptome2.G:
            if T3.G.get(gene) == None:
                T3.G[gene] = Transcriptome2.G[gene]
    TS_checked = []
    gene_dict = dict()
    if type(Gene) == list:
        for gene in Gene:
            gene_dict[gene.ID] = gene
   
    # check all TS for corasponding Gene 
    if type(TS) == list: 
        for Ts in TS:
            if gene_dict.get(Ts.gene_ID) == None and T3.G.get(Ts.gene_ID) == None:
                print('no corresponding gene object for transcript ' + Ts + ' must include gene '+ Ts +' to append ' + Ts)
            else:
                TS_checked.append(Ts)

        for Ts in TS_checked:
            if T3.TS.get(Ts.ID) != None:
                print( Ts.ID + ' already exist in transcriptome, now ' + Ts.ID + '.1')
                Ts.ID = Ts.ID +'.1'
            T3.TS[Ts.ID] = Ts
            if T3.G.get(Ts.gene_ID) == None:
                T3.G[Ts.gene_ID] = gene_dict[Ts.gene_ID]
    else:
        pass
    for Ts in T3.TS:
        T3.G[T3.TS[Ts].gene_ID].TS_ID.append(Ts)
        T3.IDs[T3.TS[Ts].ID]=[]
    for gene in T3.G:
        T3.G[gene].TS_ID = list(set(T3.G[gene].TS_ID))
        T3.IDs[T3.G[gene].ID]=[]
    gc.collect()

    return(T3)
        
# =============================================================================

def find_nth_n_str(string, leter, n):
    
    """ finds index of the nth letter in a string 
    
    find_nth_n_str(string, leter,n)
        string: str to search 
        leter: str charater to find 
        n: int how many letters to skip-1
        
    """
    ################### check inputs #########################
    if type(string) != str:
        raise TypeError('string must be of type str not '+ str(type(string)))
    if type(leter) != str:
        raise TypeError('leter must be of type str not '+ str(type(leter)))
    if type(n) != int:
        raise TypeError('n must be of type str not '+ str(type(n)))
    ##########################################################
        
    m=0
    for c, char in enumerate(string):
        if m<n:
            if char == leter:
                m+=1
            else:
                pass
        else:
            return(c-1)

# =============================================================================

def isDNA(string):
    '''Function that checks if a string is a DNA string
    
    isDNA(string)
    string: str to check if DNA or has letters in [A,C,T,G,a,c,t,g]
    
    '''
    ################### check inputs #########################
    if type(string) != str:
        return(False)
    ##########################################################
        
    condition = all(i in 'ACTGactg' for i in string)
    return(condition)

# =============================================================================

def random_DNA(DNA_len):
    DNA = ''
    Bases = "ACTG"
    for x in range(DNA_len):
        DNA += random.sample(Bases,1)[0]
    return(DNA)

# =============================================================================

def get_split_probes(transcriptome,
               probe_length,
               probe_db_out_path,
               file_name = 'probe_set_db_',
               overhang_5 = None,
               overhang_3 = None,
               ban_seqs_dict = None,
               iso_form = False,
               verbose = False, 
               exon_only = True,
               ban_list = None,
               threads = -1): 
    
    """ builds a database of probes for a transcriptome object 
    
        get_split_probes(Transcriptome, tile_gene, probe_length, probe_db_out_path, overhang_5 = None,
               overhang_3 = None, ban_seqs_dict = None, iso_form = False, verbose = False, 
               exon_only = True,threads = -1)
            transcriptome: transcriptome object to get probes for
            probe_length: int length of probes to be constructed and split, must be even int
            probe_db_out_path: str file path to stor probes 
            overhang_5: str to append to the begininig of probes 
            overhang_3: str to append to the end of probes 
            ban_seqs_dict: dict of band sequences all of same length
            iso_form: bool if true isoform porbes calulated 
            verbose: bool print progress of calulation every 20s 
            exon_only: bool if true only claulate probes for exons
            threads: int number of cores to use, if -1 use all cores
            ban_list: list of DNA sequences strings <= probe_length
            return str file path to database of probes
            
    """
    # timer 
    #tstart = time.time()
    
    ################### check inputs #########################
    # if isinstance(transcriptome,Transcriptome) != True:
    #     raise TypeError('transcriptome must be of type Transcriptome not '+ str(type(transcriptome)))
    if type(file_name) != str:
        raise TypeError('file_name must be of type str not '+ str(type(file_name)))
    
    if type(probe_length) != int:
        raise TypeError('probe_length must be of type int not '+ str(type(probe_length)))
    
    if type(probe_db_out_path) != str:
        raise TypeError('probe_db_out_path must be of type str not '+ str(type(probe_db_out_path)))
    
    if os.path.exists(os.path.normpath(probe_db_out_path)) == False:
        raise Exception('probe_db_out_path '+os.path.normpath(probe_db_out_path)+' does not exist')
        
    if type(overhang_5) != str and type(overhang_5) != type(None):
        raise TypeError('overhang_5 must be of type str or None not '+ str(type(overhang_5)))
    
    if isDNA(overhang_5) == False and type(overhang_5) != type(None):
        raise Exception('overhang_5 must be a DNA seq with bases A,C,T,G,a,c,t, or g')
    
    if type(overhang_3) != str and type(overhang_3) != type(None):
        raise TypeError('overhang_3 must be of type str or None not '+ str(type(overhang_3)))
    
    if isDNA(overhang_3) == False and type(overhang_3) != type(None):
        raise Exception('overhang_3 must be a DNA seq with bases A,C,T,G,a,c,t, or g')
    
    if type(ban_seqs_dict) != dict and type(ban_seqs_dict) != type(None):
        raise TypeError('ban_seqs_dict must be of type dict or None not '+ str(type(ban_seqs_dict)))
        for ban in ban_seqs_dict:
            if type(ban) != str:
                raise TypeError('ban_seqs_dict keys must be of type str not '+ str(type(ban)))
            if isDNA(ban) == False:
                raise Exception('ban_seqs_dict keys must be DNA with bases in \'ACTGactg\ not '+ban)
    
    if type(iso_form) != bool:
        raise TypeError('iso_form must be of type bool not '+ str(type(iso_form)))
    
    if type(verbose) != bool:
        raise TypeError('verbose must be of type bool not '+ str(type(verbose)))
        
    if type(exon_only) != bool:
        raise TypeError('exon_only must be of type bool not '+ str(type(exon_only)))
    
    if transcriptome.OT_db_path == '':
        raise Exception('ISO_dict_db_path must be defined in transcriptome')
        
    if os.path.exists(os.path.normpath(transcriptome.OT_db_path)) == False:
        raise Exception('OT_db_path in transcriptome '+os.path.normpath(transcriptome.OT_db_path)+' does not exist')
    
    if transcriptome.TDP_db_path == '':
        raise Exception('TDP_db_path must be defined in transcriptome')
        
    if os.path.exists(os.path.normpath(transcriptome.TDP_db_path)) == False:
        raise Exception('TDP_db_path in transcriptome '+os.path.normpath(transcriptome.TDP_db_path)+' does not exist')
    
    if type(ban_list) != list and type(ban_list) != type(None):
        raise TypeError('ban_list must be of type list or None not '+ str(type(ban_list)))
    
    if type(ban_list) != type(None):
        for ban in ban_list:
            if type(ban) != str:
                raise TypeError('ban_list items must be of type str not '+ str(type(ban)))
            if isDNA(ban) == False:
                raise Exception('ban_list items must be DNA with bases in \'ACTGactg\' not '+ban)
            if len(ban) > probe_length:
                raise Exception('ban_list items must be <= probe_length not length '+str(len(ban)))
        
    if probe_length%2 != 0:
        raise Exception('split probes can only be made from probes of even length')

    ##########################################################
        
    overhang_5 = overhang_5.upper()
    overhang_3 = overhang_3.upper()
    ISO_dict_db_path = transcriptome.OT_db_path
    
    # get off target length 
    OT_length = load_pickle(os.path.normpath(os.path.join(ISO_dict_db_path,'OT_length.pk')))
    OT_length = OT_length['OT_length']
    if probe_length/2 < OT_length:
        raise Exception('probe_length must be larger than OT_length')
    
    # set up ban sequences length 
    if ban_seqs_dict == None:
        ban_seqs_per_probe = 0
        ban_seqs_len = probe_length
        oh5_len = 0
        oh3_len = 0
    else:
        ban_seqs_len = len(list(ban_seqs_dict.keys())[0])
        #overhang lengths for ban seq testing
        if type(overhang_5) == str:
          oh5_len = len(overhang_5)
        else:
          oh5_len = 0
        if type(overhang_3) == str:
          oh3_len = len(overhang_3)
        else:
          oh3_len = 0
        ban_seqs_per_probe = probe_length+min(oh5_len,ban_seqs_len)+min(oh3_len,ban_seqs_len)-ban_seqs_len

    if True:
        RNA_list_len= transcriptome.n_transcripts()
        TDP_CG_db_path = transcriptome.TDP_db_path

        # make db directory with unique name in probe_db_out_path
        number = 0
        while True:
            db_path = os.path.normpath(os.path.join(probe_db_out_path,file_name + str(number)))
            if os.path.exists(db_path) == True:
                number +=1
            else:
                break
        os.makedirs(db_path)
        
        # this is used in calulating the offtarget indexes 
        OT_seqs_per_probe = int((probe_length/2)-OT_length)
        OT_seqs_per_center = int(OT_length-1)
        center_start_index = OT_seqs_per_probe+1        
        # this is weird becaue the ot and ban tables are build of coading strands
        oh3 = reverse_complement(overhang_5[:-min(oh5_len,ban_seqs_len)])
        oh5 = reverse_complement(overhang_3[0:min(oh3_len,ban_seqs_len)])
        
        overhang_3_full = overhang_3
        overhang_5_full = overhang_5
        
        # set up empty list ot check if ban_list is None 
        if ban_list == None:
            ban_list = []
        
        # make probe set for each trascript, function for multi thread  
        def GetProbes(GID,
                      probe_length,
                      db_path,
                      OT_length,
                      oh3,
                      oh5,
                      ban_seqs_dict,
                      ban_seqs_per_probe,
                      ban_seqs_len,
                      TDP_CG_db_path,
                      iso_form,
                      verbose,
                      exon_only,
                      ban_list,
                      overhang_5_full,
                      overhang_3_full):
            
            probe_count = 0
            G = GID
            G = transcriptome.G[G]
            # load off target calulations and thrmo proporties 
            ISO_dict = load_pickle(os.path.normpath(os.path.join(ISO_dict_db_path,G.ID+'.pk')))
            thurmo = load_pickle(os.path.normpath(os.path.join(TDP_CG_db_path,G.ID+'.pk')))
            ALL = [[],[],[],[],[],transcriptome.G[G.ID].symbol,transcriptome.G[G.ID].chr,[],[]]
            gene_path = os.path.normpath(os.path.join(db_path,G.ID))
            if os.path.exists(gene_path) == False:
                os.makedirs(gene_path)
            else:
                if verbose == True:
                    print(gene_path+' already exist')
            
            for TS in G.TS_ID:
                RNA = transcriptome.TS[TS]
                G_ID = RNA.gene_ID
                TS_ID = RNA.ID

                #t = time.time()
                #if verbose == True and t > t1:
                #    print(str(c)+'/'+str(RNA_list_len) + ' probe sets made')
                #    t1 = t + 20

                #check gene iso_dict and TDP exist
                if ISO_dict == None or thurmo == None:
                    if verbose == True:
                        print('no TDP or specificity info for '+ TS)
                    continue
                if thurmo.get(TS_ID) == None:
                    if verbose == True:
                        print('no TDP or specificity info for '+ TS)
                    continue
                if ISO_dict.get(TS_ID) == None:
                    if verbose == True:
                        print('no TDP or specificity info for '+ TS)
                    continue

                probe_data = [[],[],[],[],[],transcriptome.G[G_ID].symbol,transcriptome.G[G_ID].chr,[],[]]
                 # add list [[seq],[[H,S]],[CG],[Iso_index],[gene_index],gene symbol, [chr], [start], [stop posiotional_index_on ch],]     
               
                s = len(RNA.seq)
                p = s-probe_length+1 # number of probe doubles
                split_probe_size = int(probe_length/2)
                for b in range(p):
                    #TDPCG=[TDP_H,TDP_S,CG] TDP_CG[G_ID][TS_ID] = TDPCG #TDP_CG=[TDP_H,TDP_S,CG]
                    start = b
                    end = b + probe_length

                    probe = RNA.seq[start:end] # this is the RC of final probe 

                    left_probe = probe[split_probe_size:] # because RC of probe 
                    
                    right_probe = probe[:split_probe_size]
                    
                    probe_with_full_oh = overhang_5_full+reverse_complement(probe)+overhang_3_full
                    left_probe_with_full_oh = overhang_5_full+reverse_complement(left_probe)
                    right_probe_with_full_oh = reverse_complement(right_probe)+overhang_3_full

                    #check band seqs with overhang
                    ban = False
                    for st in range(ban_seqs_per_probe):
                        probe_with_oh = oh5+probe+oh3 # oh are right with above coad 
                        to_check = probe_with_oh[st:(st+ban_seqs_len)]
                        if ban_seqs_dict.get(to_check) == None:
                            ban = False
                        else:
                            ban = True
                        if ban == True:
                          # if verbose == True:
                          #     print('ban found')
                          continue
                    if ban == True:
                        continue
                    
                    #check band list with full overhangs
                    ban = False
                    for bans in ban_list:
                        if bans in probe_with_full_oh:
                            ban = True
                        if ban == True:
                            # if verbose == True:
                            #     print('ban found')
                            continue
                    if ban == True:
                        continue
                    
                    # calulate global position
                    n_exon = len(RNA.in_exon_index_start) - 1 # number of exons 
                    if start >= RNA.in_exon_index_start[-1]: # check if in last exon 
                        s_exon = n_exon
                        e_exon = n_exon
                    # find which exon probe is in
                    else:     
                        s_exon = 'x'
                        e_exon = 'x'
                        for e, start_on_TS in enumerate(RNA.in_exon_index_start):
                            # check if in same exon 
                            if start >= start_on_TS and end < RNA.in_exon_index_start[e+1]: 
                                s_exon = e
                                e_exon = e
                                break
                            # check if in span 2 exons
                            elif start >= start_on_TS and start < RNA.in_exon_index_start[e+1]:
                                s_exon = e 
                                e_exon = e+1
                                break

                    strand = transcriptome.G[RNA.gene_ID].strand
                    if strand == 1:
                         strand = '1'
                    if strand == -1:
                         strand = '-1'
                    if strand not in ['-1','1']:
                        strand = '1'
                    if strand == '1':
                        if s_exon != e_exon: # if start and end in diffrent exon
                            start_global = [(start - RNA.in_exon_index_start[s_exon]) + RNA.ex_exon_index_start[s_exon], RNA.ex_exon_index_start[s_exon+1]]
                            end_global = [RNA.ex_exon_index_stop[s_exon],(end - RNA.in_exon_index_start[e_exon]) + RNA.ex_exon_index_start[e_exon]]
                        else: # if in same exons 
                            start_global = [(start - RNA.in_exon_index_start[s_exon]) + RNA.ex_exon_index_start[s_exon]] 
                            end_global = [(end - RNA.in_exon_index_start[e_exon]) + RNA.ex_exon_index_start[e_exon]] 
                    elif strand =='-1':
                        if s_exon != e_exon: # if start and end in diffrent exon
                            start_global = [RNA.ex_exon_index_stop[s_exon] + (RNA.in_exon_index_stop[s_exon] - start), RNA.ex_exon_index_start[s_exon+1]]
                            end_global = [RNA.ex_exon_index_stop[s_exon], RNA.ex_exon_index_start[e_exon] - (end - RNA.in_exon_index_start[e_exon]) ]
                        else: # if in same exons 
                            start_global = [RNA.ex_exon_index_start[s_exon]  - (start - RNA.in_exon_index_start[s_exon]) ] 
                            end_global = [RNA.ex_exon_index_start[e_exon] - (end - RNA.in_exon_index_start[e_exon])] 

                    #calulate Tm and CG for both probes left and right
                    Hl = sum(thurmo[TS_ID][0][start:end+1][split_probe_size:])
                    Sl = sum(thurmo[TS_ID][1][start:end+1][split_probe_size:])
                    CGl = sum(thurmo[TS_ID][2][start:end+1][split_probe_size:])/(split_probe_size)
                    
                    Hr = sum(thurmo[TS_ID][0][start:end+1][:split_probe_size])
                    Sr = sum(thurmo[TS_ID][1][start:end+1][:split_probe_size])
                    CGr = sum(thurmo[TS_ID][2][start:end+1][:split_probe_size])/(split_probe_size)
                    
                    # add ending comp for AT left
                    if start == 0:
                        fiveprime = 0
                    elif RNA.seq[start-1] in ('A','T'):
                        fiveprime = 1
                    else:
                        fiveprime = 0
                    if RNA.seq[start+split_probe_size+1] in ('A','T'):
                        threeprime = 1
                    else:
                        threeprime = 0
                    Hl = Hl + 0.2 + (2.2*fiveprime) + (2.2*threeprime)
                    Sl = Sl -5.7 + (6.9*fiveprime) + (6.9*threeprime)
                    HSl=[Hl,Sl]
                    
                    # add ending comp for AT right
                    if RNA.seq[start+split_probe_size-1] in ('A','T'):
                        fiveprime = 1
                    else:
                        fiveprime = 0
                    if end == s:
                        threeprime = 0
                    elif end+1 >= s:
                        threeprime = 0
                    elif RNA.seq[end+1] in ('A','T'):
                        threeprime = 1
                    else:
                        threeprime = 0
                    Hr = Hr + 0.2 + (2.2*fiveprime) + (2.2*threeprime)
                    Sr = Sr -5.7 + (6.9*fiveprime) + (6.9*threeprime)
                    HSr=[Hr,Sr]
                    
                    #get Iso_index and gene_index
                    if iso_form == True:
                        #left 
                        probe_iso_indexl = 0
                        probe_gene_indexl = 0
                        for sta in range(OT_seqs_per_probe):# [iso_pen,iso_index,gene_index]
                            seql = left_probe[sta:(sta+OT_length)]
                            seql = DNA2int(seql)
                            indexl = ISO_dict[TS_ID].get(seql,[1,1])
                            probe_iso_indexl += indexl[0]
                            probe_gene_indexl += indexl[1]
                            # normilize sums of probe_iso_index and  probe_gene_index
                        probe_iso_indexl = probe_iso_indexl/OT_seqs_per_probe
                        probe_gene_indexl = probe_gene_indexl/OT_seqs_per_probe
                        
                        #right 
                        probe_iso_indexr = 0
                        probe_gene_indexr = 0
                        for sta in range(OT_seqs_per_probe):# [iso_pen,iso_index,gene_index]
                            seqr = right_probe[sta:(sta+OT_length)]
                            seqr = DNA2int(seqr)
                            indexr = ISO_dict[TS_ID].get(seqr,[1,1])
                            probe_iso_indexr += indexr[0]
                            probe_gene_indexr += indexr[1]
                            # normilize sums of probe_iso_index and  probe_gene_index
                        probe_iso_indexr = probe_iso_indexr/OT_seqs_per_probe
                        probe_gene_indexr = probe_gene_indexr/OT_seqs_per_probe
                            # add to tanscript specfice probe list 
                            #[[seq],[Tm],[CG],[Iso_index],[gene_index],[posiotional_index]
                        
                        # center 
                        probe_iso_indexc = 0
                        probe_gene_indexc = 0
                        center = probe[center_start_index:int(center_start_index+(2*OT_length))]
                        for sta in range(OT_seqs_per_center):# [iso_pen,iso_index,gene_index]
                            seqc = center[sta:(sta+OT_length)]
                            seqc = DNA2int(seqc)
                            indexc = ISO_dict[TS_ID].get(seqc,[1,1])
                            probe_iso_indexc += indexc[0]
                            probe_gene_indexc += indexc[1]
                            # normilize sums of probe_iso_index and  probe_gene_index
                        probe_iso_indexc = probe_iso_indexc/OT_seqs_per_center
                        probe_gene_indexc = probe_gene_indexc/OT_seqs_per_center
                            # add to tanscript specfice probe list 
                            #[[seq],[Tm],[CG],[Iso_index],[gene_index],[posiotional_index]
                        
                        probe_data[0].append([DNA2int(reverse_complement(left_probe)),DNA2int(reverse_complement(right_probe))])
                        probe_count+=1
                        probe_data[1].append([HSl,HSr])
                        probe_data[2].append([CGl,CGr])
                        probe_data[3].append([probe_iso_indexl,probe_iso_indexr,probe_iso_indexc])
                        probe_data[4].append([probe_gene_indexl,probe_gene_indexr,probe_gene_indexc])
                        probe_data[7].append(start_global)
                        probe_data[8].append(end_global)
                        
                    else:
                        #left
                        probe_gene_indexl = 0
                        for sta in range(OT_seqs_per_probe):# [iso_pen,iso_index,gene_index]
                            seql = left_probe[sta:(sta+OT_length)]
                            seql = DNA2int(seql)
                            indexl = ISO_dict[TS_ID].get(seql,[1,1])
                            probe_gene_indexl += indexl[1]
                            # normilize sums of probe_gene_index
                        probe_gene_indexl = probe_gene_indexl/OT_seqs_per_probe
                        
                        #right 
                        probe_gene_indexr = 0
                        for sta in range(OT_seqs_per_probe):# [iso_pen,iso_index,gene_index]
                            seqr = right_probe[sta:(sta+OT_length)]
                            seqr = DNA2int(seqr)
                            indexr = ISO_dict[TS_ID].get(seqr,[1,1])
                            probe_gene_indexr += indexr[1]
                            # normilize sums of probe_gene_index
                        probe_gene_indexr = probe_gene_indexr/OT_seqs_per_probe
                        
                        # center 
                        probe_gene_indexc = 0
                        center = probe[center_start_index:int(center_start_index+(2*OT_length))]
                        for sta in range(OT_seqs_per_center):# [iso_pen,iso_index,gene_index]
                            seqc = center[sta:(sta+OT_length)]
                            seqc = DNA2int(seqc)
                            indexc = ISO_dict[TS_ID].get(seqc,[1,1])
                            probe_gene_indexc += indexc[1]
                            # normilize sums of probe_iso_index and  probe_gene_index
                        probe_gene_indexc = probe_gene_indexc/OT_seqs_per_center

                        
                        # add to gene probe list 
                        probe_data[0].append([DNA2int(reverse_complement(left_probe)),DNA2int(reverse_complement(right_probe))])
                        probe_count+=1
                        probe_data[1].append([HSl,HSr])
                        probe_data[2].append([CGl,CGr])
                        probe_data[3].append([0,0,0])
                        probe_data[4].append([probe_gene_indexl,probe_gene_indexr,probe_gene_indexc])
                        probe_data[7].append(start_global)
                        probe_data[8].append(end_global)
                
                if exon_only == True:
                    #if iso_form == True:
                    save_pickle(probe_data,gene_path,TS_ID)
                    all_obj = ALL
                    all_obj[0] += probe_data[0]
                    all_obj[1] += probe_data[1]
                    all_obj[2] += probe_data[2]
                    all_obj[3] += probe_data[3]
                    all_obj[4] += probe_data[4]
                    all_obj[7] += probe_data[7]
                    all_obj[8] += probe_data[8]
                    #save_pickle(ALL, gene_path,'all')
                    continue
                        #
                #save_pickle(probe_data,gene_path,TS_ID)
                #update all 
                        # all_obj = ALL
                        # all_obj[0] += probe_data[0]
                        # all_obj[1] += probe_data[1]
                        # all_obj[2] += probe_data[2]
                        # all_obj[3] += probe_data[3]
                        # all_obj[4] += probe_data[4]
                        # all_obj[7] += probe_data[7]
                        # all_obj[8] += probe_data[8]
                        # save_pickle(ALL, gene_path,'all')
                        

                #add probes for non-spliced genes
                
                strand = transcriptome.G[RNA.gene_ID].strand
                gene = transcriptome.G[G_ID]
                if strand == '1'or strand != '-1':
                  unspliced_seq = gene.seq[RNA.start-gene.start:RNA.stop-gene.start]
                else:
                  unspliced_seq = gene.seq[abs(RNA.start-gene.stop):abs(RNA.stop-gene.stop)]
                s = len(unspliced_seq)
                p = s-probe_length+1 # number of probes
                for b in range(p):
                    start = b
                    end = b + probe_length
                    probe = unspliced_seq[start:end]
                    
                    left_probe = probe[split_probe_size:] # because RC of probe 
                    right_probe = probe[:split_probe_size]
                    
                    probe_with_full_oh = overhang_5_full+reverse_complement(probe)+overhang_3_full
                    left_probe_with_full_oh = overhang_5_full+reverse_complement(left_probe)
                    right_probe_with_full_oh = reverse_complement(right_probe)+overhang_3_full


                    #check band seqs with overhang
                    ban = False
                    for st in range(ban_seqs_per_probe):
                        probe_with_oh = oh5+probe+oh3 # oh are right with above coad 
                        to_check = probe_with_oh[st:(st+ban_seqs_len)]
                        if ban_seqs_dict.get(to_check) == None:
                            ban = False
                        else:
                            ban = True
                        if ban == True:
                          # if verbose == True:
                          #     print('ban found')
                          continue
                    if ban == True:
                        continue
                    
                    #check band list with full overhangs
                    ban = False
                    for bans in ban_list:
                        if bans in probe_with_full_oh:
                            ban = True
                        if ban == True:
                            # if verbose == True:
                            #     print('ban found')
                            continue
                    if ban == True:
                        continue


                    # calulate global position
                    if strand == '1'or strand != '-1':
                      global_start = RNA.start + b
                      global_stop = global_start + probe_length
                    else:
                      global_start = RNA.start - b
                      global_stop = global_start - probe_length

                    #calulate Tm and CG
                    if strand == '1'or strand != '-1':
                      g_start = gene.start
                    else:
                      g_start = gene.stop
                    TDP_start = abs(global_start-g_start)
                    TDP_stop = TDP_start + probe_length
                    TDP_ID = 'NSG'+ G_ID
                                        
                    #calulate Tm and CG for both probes left and right
                    Hl = sum(thurmo[TDP_ID][0][TDP_start:TDP_stop+1][split_probe_size:])
                    Sl = sum(thurmo[TDP_ID][1][TDP_start:TDP_stop+1][split_probe_size:])
                    CGl = sum(thurmo[TDP_ID][2][TDP_start:TDP_stop+1][split_probe_size:])/(split_probe_size)
                    
                    Hr = sum(thurmo[TDP_ID][0][TDP_start:TDP_stop+1][:split_probe_size])
                    Sr = sum(thurmo[TDP_ID][1][TDP_start:TDP_stop+1][:split_probe_size])
                    CGr = sum(thurmo[TDP_ID][2][TDP_start:TDP_stop+1][:split_probe_size])/(split_probe_size)
                    

                    # add ending comp for AT left
                    if start == 0:
                        fiveprime = 0
                    elif unspliced_seq[start-1] in ('A','T'):
                        fiveprime = 1
                    else:
                        fiveprime = 0
                    if unspliced_seq[start+split_probe_size+1] in ('A','T'):
                        threeprime = 1
                    else:
                        threeprime = 0
                    Hl = Hl + 0.2 + (2.2*fiveprime) + (2.2*threeprime)
                    Sl = Sl -5.7 + (6.9*fiveprime) + (6.9*threeprime)
                    HSl=[Hl,Sl]
                    
                    # add ending comp for AT right
                    if unspliced_seq[start+split_probe_size-1] in ('A','T'):
                        fiveprime = 1
                    else:
                        fiveprime = 0
                    if end == s:
                        threeprime = 0
                    elif end+1 >= s:
                        threeprime = 0
                    elif unspliced_seq[end+1] in ('A','T'):
                        threeprime = 1
                    else:
                        threeprime = 0
                    Hr = Hr + 0.2 + (2.2*fiveprime) + (2.2*threeprime)
                    Sr = Sr -5.7 + (6.9*fiveprime) + (6.9*threeprime)
                    HSr=[Hr,Sr]
                    
                    if iso_form == True:
                        #left 
                        probe_iso_indexl = 0
                        probe_gene_indexl = 0
                        for sta in range(OT_seqs_per_probe):# [iso_pen,iso_index,gene_index]
                            seql = left_probe[sta:(sta+OT_length)]
                            seql = DNA2int(seql)
                            indexl = ISO_dict[TS_ID].get(seql,[1,1])
                            probe_iso_indexl += indexl[0]
                            probe_gene_indexl += indexl[1]
                            # normilize sums of probe_iso_index and  probe_gene_index
                        probe_iso_indexl = probe_iso_indexl/OT_seqs_per_probe
                        probe_gene_indexl = probe_gene_indexl/OT_seqs_per_probe
                        
                        #right 
                        probe_iso_indexr = 0
                        probe_gene_indexr = 0
                        for sta in range(OT_seqs_per_probe):# [iso_pen,iso_index,gene_index]
                            seqr = right_probe[sta:(sta+OT_length)]
                            seqr = DNA2int(seqr)
                            indexr = ISO_dict[TS_ID].get(seqr,[1,1])
                            probe_iso_indexr += indexr[0]
                            probe_gene_indexr += indexr[1]
                            # normilize sums of probe_iso_index and  probe_gene_index
                        probe_iso_indexr = probe_iso_indexr/OT_seqs_per_probe
                        probe_gene_indexr = probe_gene_indexr/OT_seqs_per_probe
                            # add to tanscript specfice probe list 
                            #[[seq],[Tm],[CG],[Iso_index],[gene_index],[posiotional_index]
                        
                        # center 
                        probe_iso_indexc = 0
                        probe_gene_indexc = 0
                        center = probe[center_start_index:int(center_start_index+(2*OT_length))]
                        for sta in range(OT_seqs_per_center):# [iso_pen,iso_index,gene_index]
                            seqc = center[sta:(sta+OT_length)]
                            seqc = DNA2int(seqc)
                            indexc = ISO_dict[TS_ID].get(seqc,[1,1])
                            probe_iso_indexc += indexc[0]
                            probe_gene_indexc += indexc[1]
                            # normilize sums of probe_iso_index and  probe_gene_index
                        probe_iso_indexc = probe_iso_indexc/OT_seqs_per_center
                        probe_gene_indexc = probe_gene_indexc/OT_seqs_per_center
                        
                        probe_data[0].append([DNA2int(reverse_complement(left_probe)),DNA2int(reverse_complement(right_probe))])
                        probe_count+=1
                        probe_data[1].append([HSl,HSr])
                        probe_data[2].append([CGl,CGr])
                        probe_data[3].append([probe_iso_indexl,probe_iso_indexr, probe_iso_indexc])
                        probe_data[4].append([probe_gene_indexl,probe_gene_indexr,probe_gene_indexc])
                        probe_data[7].append([global_start])
                        probe_data[8].append([global_stop])
                        
                    else:
                        #left
                        probe_gene_indexl = 0
                        for sta in range(OT_seqs_per_probe):# [iso_pen,iso_index,gene_index]
                            seql = left_probe[sta:(sta+OT_length)]
                            seql = DNA2int(seql)
                            indexl = ISO_dict[TS_ID].get(seql,[1,1])
                            probe_gene_indexl += indexl[1]
                            # normilize sums of probe_gene_index
                        probe_gene_indexl = probe_gene_indexl/OT_seqs_per_probe
                        
                        # right
                        probe_gene_indexr = 0
                        for sta in range(OT_seqs_per_probe):# [iso_pen,iso_index,gene_index]
                            seqr = right_probe[sta:(sta+OT_length)]
                            seqr = DNA2int(seqr)
                            indexr = ISO_dict[TS_ID].get(seqr,[1,1])
                            probe_gene_indexr += indexr[1]
                            # normilize sums of probe_gene_index
                        probe_gene_indexr = probe_gene_indexr/OT_seqs_per_probe
                        
                        # center 
                        probe_gene_indexc = 0
                        center = probe[center_start_index:int(center_start_index+(2*OT_length))]
                        for sta in range(OT_seqs_per_center):# [iso_pen,iso_index,gene_index]
                            seqc = center[sta:(sta+OT_length)]
                            seqc = DNA2int(seqc)
                            indexc = ISO_dict[TS_ID].get(seqc,[1,1])
                            probe_gene_indexc += indexc[1]
                            # normilize sums of probe_iso_index and  probe_gene_index
                        probe_gene_indexc = probe_gene_indexc/OT_seqs_per_center
                        
                        # add to gene probe list 
                        probe_data[0].append([DNA2int(reverse_complement(left_probe)),DNA2int(reverse_complement(right_probe))])
                        probe_count+=1
                        probe_data[1].append([HSl,HSr])
                        probe_data[2].append([CGl,CGr])
                        probe_data[3].append([0,0,0])
                        probe_data[4].append([probe_gene_indexl,probe_gene_indexr,probe_gene_indexc])
                        probe_data[7].append([global_start])
                        probe_data[8].append([global_stop])

                if iso_form == True:
                    save_pickle(probe_data,gene_path,TS_ID)
        
                #update all 
                all_obj = ALL
                all_obj[0] += probe_data[0]
                all_obj[1] += probe_data[1]
                all_obj[2] += probe_data[2]
                all_obj[3] += probe_data[3]
                all_obj[4] += probe_data[4]
                all_obj[7] += probe_data[7]
                all_obj[8] += probe_data[8]
  
            save_pickle(ALL, gene_path,'all')
            return(probe_count)

        # frag data for multicore 
        GID_ = list(transcriptome.G.keys())
        genes_n = len(GID_)
        probe_length_ = [probe_length] * genes_n
        db_path_ = [db_path] * genes_n
        OT_length_ = [OT_length] * genes_n
        oh3_ = [oh3] * genes_n
        oh5_ = [oh5] * genes_n
        ban_seqs_dict_ = [ban_seqs_dict] * genes_n
        ban_seqs_per_probe_ = [ban_seqs_per_probe] * genes_n
        ban_seqs_len_ = [ban_seqs_len] * genes_n
        TDP_CG_db_path_ = [TDP_CG_db_path] * genes_n
        iso_form_ = [iso_form] * genes_n
        verbose_ = [verbose] * genes_n
        exon_only_ = [exon_only] * genes_n
        ban_list_ = [ban_list] * genes_n
        overhang_5_full_ = [overhang_5_full] * genes_n
        overhang_3_full_ = [overhang_3_full] * genes_n
        ziped = zip(GID_,
                      probe_length_,
                      db_path_,
                      OT_length_,
                      oh3_,
                      oh5_,
                      ban_seqs_dict_,
                      ban_seqs_per_probe_,
                      ban_seqs_len_,
                      TDP_CG_db_path_,
                      iso_form_,
                      verbose_,
                      exon_only_,
                      ban_list_,
                      overhang_5_full_,
                      overhang_3_full_)
        
        # set up threads/cores 
        if threads == -1: 
            threads = multiprocessing.cpu_count()
        
        if threads != 1:
            pool = ThreadPool(threads)
            # run calulations on mulitiple cores 
            results = pool.starmap(GetProbes,zip(GID_,
                      probe_length_,
                      db_path_,
                      OT_length_,
                      oh3_,
                      oh5_,
                      ban_seqs_dict_,
                      ban_seqs_per_probe_,
                      ban_seqs_len_,
                      TDP_CG_db_path_,
                      iso_form_,
                      verbose_,
                      exon_only_,
                      ban_list_,
                      overhang_5_full_,
                      overhang_3_full_),1)#
            pool.close() 
            pool.join()
        else:
            results = []
            for a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,p1 in ziped:
                x = GetProbes(a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,p1)
                results.append(x)
        
        probe_count = sum(results)
        
        # calulate time to compleate making probes 
        #tend = time.time()
        #time_total = tend - tstart
        #if verbose == True:
        #    print(str(round((time_total/60),2))+' min to make ' + str(probe_count)+ ' probes')
            
        # save a metadata file 
        ban_seqs_used = False
        ban_seqs_path = 'None'
        if ban_seqs_dict != None:
            ban_seqs_used = True
            save_pickle(ban_seqs_dict, db_path, 'ban_seqs_dict')
            ban_seqs_path = os.path.normpath(os.path.join(db_path,'ban_seqs_dict.pk'))
        else:
            ban_seqs_path = ''
            
        ban_list_used = False
        if ban_list != []:
            ban_list_used = True
            save_pickle(ban_list, db_path, 'ban_seqs_list')
            ban_seqs_list_path = os.path.normpath(os.path.join(db_path,'ban_seqs_list.pk'))
        else: 
            ban_seqs_list_path = ''
        
        meta_path = os.path.normpath(os.path.join(db_path,'probe_metadata.txt'))
        meta_file = open(meta_path, "w") 
        to_write = ['transcriptome cell type: '+ transcriptome.cell_type,
                    'thermodynamic properties: '+transcriptome.TDP_db_path,
                    'off target properties: '+transcriptome.OT_db_path,
                    'probe_length: '+str(probe_length),
                    'probe count: '+str(probe_count),
                    'probe_db_out_path: '+db_path,
                    '5 prime overhang: '+str(overhang_5),
                    '3 prime overhang: '+str(overhang_3),
                    'isoform included: '+str(iso_form),
                    'ban_seqs_dict used: '+str(ban_seqs_used),
                    'ban_seqs_dict path: '+str(ban_seqs_path),
                    'ban_list used: '+str(ban_list_used),
                    'ban_list path: '+str(ban_seqs_list_path),
                    'split probes: True',
                    #'time to make probes (min): '+str(round((time_total/60),2)),
                    'cores used: '+str(threads)]
    
        for line in to_write:
          meta_file.write(line)
          meta_file.write("\n")
        meta_file.close() 

    if verbose == True:
        print('finished making probes')
    gc.collect()    
    return(db_path)

# =============================================================================
    
def get_probes(transcriptome,
               probe_length,
               probe_db_out_path,
               overhang_5 = None,
               overhang_3 = None,
               ban_seqs_dict = None,
               iso_form = False,
               verbose = False, 
               exon_only = True,
               ban_list = None,
               threads = -1): 
    
    """ builds a database of probes for a transcriptome object 
    
        get_probes(Transcriptome, probe_length, probe_db_out_path, overhang_5 = None,
               overhang_3 = None, ban_seqs_dict = None, iso_form = False, verbose = False, 
               exon_only = True,threads = -1)
            transcriptome: transcriptome object to get probes for
            probe_length: int length of probes to be constructed 
            probe_db_out_path: str file path to stor probes 
            overhang_5: str to append to the begininig of probes 
            overhang_3: str to append to the end of probes 
            ban_seqs_dict: dict of band sequences all of same length
            iso_form: bool if true isoform porbes calulated 
            verbose: bool print progress of calulation every 20s 
            exon_only: bool if true only claulate probes for exons
            threads: int number of cores to use, if -1 use all cores
            ban_list: list of DNA sequences strings <= probe_length
            return str file path to database of probes
            
    """
    # timer 
    tstart = time.time()
    
    ################### check inputs #########################
    # if isinstance(transcriptome,Transcriptome) != True:
    #     raise TypeError('transcriptome must be of type Transcriptome not '+ str(type(transcriptome)))
    
    if type(probe_length) != int:
        raise TypeError('probe_length must be of type int not '+ str(type(probe_length)))
    
    if type(probe_db_out_path) != str:
        raise TypeError('probe_db_out_path must be of type str not '+ str(type(probe_db_out_path)))
    
    if os.path.exists(os.path.normpath(probe_db_out_path)) == False:
        raise Exception('probe_db_out_path '+os.path.normpath(probe_db_out_path)+' does not exist')
        
    if type(overhang_5) != str and type(overhang_5) != type(None):
        raise TypeError('overhang_5 must be of type str or None not '+ str(type(overhang_5)))
    
    if isDNA(overhang_5) == False and type(overhang_5) != type(None):
        raise Exception('overhang_5 must be a DNA seq with bases A,C,T,G,a,c,t, or g')
    
    if type(overhang_3) != str and type(overhang_3) != type(None):
        raise TypeError('overhang_3 must be of type str or None not '+ str(type(overhang_3)))
    
    if isDNA(overhang_3) == False and type(overhang_3) != type(None):
        raise Exception('overhang_3 must be a DNA seq with bases A,C,T,G,a,c,t, or g')
    
    if type(ban_seqs_dict) != dict and type(ban_seqs_dict) != type(None):
        raise TypeError('ban_seqs_dict must be of type dict or None not '+ str(type(ban_seqs_dict)))
        for ban in ban_seqs_dict:
            if type(ban) != str:
                raise TypeError('ban_seqs_dict keys must be of type str not '+ str(type(ban)))
            if isDNA(ban) == False:
                raise Exception('ban_seqs_dict keys must be DNA with bases in \'ACTGactg\ not '+ban)
    
    if type(iso_form) != bool:
        raise TypeError('iso_form must be of type bool not '+ str(type(iso_form)))
    
    if type(verbose) != bool:
        raise TypeError('verbose must be of type bool not '+ str(type(verbose)))
        
    if type(exon_only) != bool:
        raise TypeError('exon_only must be of type bool not '+ str(type(exon_only)))
    
    if transcriptome.OT_db_path == '':
        raise Exception('ISO_dict_db_path must be defined in transcriptome')
        
    if os.path.exists(os.path.normpath(transcriptome.OT_db_path)) == False:
        raise Exception('OT_db_path in transcriptome '+os.path.normpath(transcriptome.OT_db_path)+' does not exist')
    
    if transcriptome.TDP_db_path == '':
        raise Exception('TDP_db_path must be defined in transcriptome')
        
    if os.path.exists(os.path.normpath(transcriptome.TDP_db_path)) == False:
        raise Exception('TDP_db_path in transcriptome '+os.path.normpath(transcriptome.TDP_db_path)+' does not exist')
    
    if type(ban_list) != list and type(ban_list) != type(None):
        raise TypeError('ban_list must be of type list or None not '+ str(type(ban_list)))
    
    if type(ban_list) != type(None):
        for ban in ban_list:
            if type(ban) != str:
                raise TypeError('ban_list items must be of type str not '+ str(type(ban)))
            if isDNA(ban) == False:
                raise Exception('ban_list items must be DNA with bases in \'ACTGactg\ not '+ban)
            if len(ban) > probe_length:
                raise Exception('ban_list items must be <= probe_length not length '+str(len(ban)))    

    ##########################################################
    if overhang_5 == None:
        overhang_5 = ''
    if overhang_3 == None:
        overhang_3 = ''
    overhang_5 = overhang_5.upper()
    overhang_3 = overhang_3.upper()
    ISO_dict_db_path = transcriptome.OT_db_path
    
    # get off target length 
    OT_length = load_pickle(os.path.normpath(os.path.join(ISO_dict_db_path,'OT_length.pk')))
    OT_length = OT_length['OT_length']
    if probe_length < OT_length:
        raise Exception('probe_length must be larger than OT_length')
    
    # set up ban sequences length 
    if ban_seqs_dict == None:
        ban_seqs_per_probe = 0
        ban_seqs_len = probe_length
        oh5_len = 0
        oh3_len = 0
    else:
        ban_seqs_len = len([ban_seqs_dict.keys()][0])
        #overhang lengths for ban seq testing
        if type(overhang_5) == str:
          oh5_len = len(overhang_5)
        else:
          oh5_len = 0
        if type(overhang_3) == str:
          oh3_len = len(overhang_3)
        else:
          oh3_len = 0
        ban_seqs_per_probe = probe_length+min(oh5_len,ban_seqs_len)+min(oh3_len,ban_seqs_len)-ban_seqs_len

    if True:
        RNA_list_len= transcriptome.n_transcripts()
        TDP_CG_db_path = transcriptome.TDP_db_path

        # make db directory with unique name in probe_db_out_path
        number = 0
        while True:
            db_path = os.path.normpath(os.path.join(probe_db_out_path,'probe_set_db_' + str(number)))
            if os.path.exists(db_path) == True:
                number +=1
            else:
                break
        os.makedirs(db_path)
        OT_seqs_per_probe = probe_length-OT_length
                
        # this is weird becaue the ot and ban tables are build of coading strands
        oh3 = reverse_complement(overhang_5[:-min(oh5_len,ban_seqs_len)])
        oh5 = reverse_complement(overhang_3[0:min(oh3_len,ban_seqs_len)])
        
        overhang_3_full = overhang_3
        overhang_5_full = overhang_5
        
        # set up empty list ot check if ban_list is None 
        if ban_list == None:
            ban_list = []
        
        # make probe set for each trascript, function for multi thread  
        def GetProbes(GID,
                      probe_length,
                      db_path,
                      OT_length,
                      oh3,
                      oh5,
                      ban_seqs_dict,
                      ban_seqs_per_probe,
                      ban_seqs_len,
                      TDP_CG_db_path,
                      iso_form,
                      verbose,
                      exon_only,
                      ban_list,
                      overhang_5_full,
                      overhang_3_full):
            
            probe_count = 0
            G = GID
            G = transcriptome.G[G]
            # load off target calulations and thrmo proporties 
            ISO_dict = load_pickle(os.path.normpath(os.path.join(ISO_dict_db_path,G.ID+'.pk')))
            thurmo = load_pickle(os.path.normpath(os.path.join(TDP_CG_db_path,G.ID+'.pk')))
            ALL = [[],[],[],[],[],transcriptome.G[G.ID].symbol,transcriptome.G[G.ID].chr,[],[]]
            gene_path = os.path.normpath(os.path.join(db_path,G.ID))
            if os.path.exists(gene_path) == False:
                os.makedirs(gene_path)
            else:
                if verbose == True:
                    print(gene_path+' already exist')
            
            for TS in G.TS_ID:
                RNA = transcriptome.TS[TS]
                G_ID = RNA.gene_ID
                TS_ID = RNA.ID

                #t = time.time()
                #if verbose == True and t > t1:
                #    print(str(c)+'/'+str(RNA_list_len) + ' probe sets made')
                #    t1 = t + 20

                #check gene iso_dict and TDP exist
                if ISO_dict == None or thurmo == None:
                    if verbose == True:
                        print('no TDP or specificity info for '+ TS)
                    continue
                if thurmo.get(TS_ID) == None:
                    if verbose == True:
                        print('no TDP or specificity info for '+ TS)
                    continue
                if ISO_dict.get(TS_ID) == None:
                    if verbose == True:
                        print('no TDP or specificity info for '+ TS)
                    continue

                probe_data = [[],[],[],[],[],transcriptome.G[G_ID].symbol,transcriptome.G[G_ID].chr,[],[]]
                 # add list [[seq],[[H,S]],[CG],[Iso_index],[gene_index],gene symbol, [chr], [start], [stop posiotional_index_on ch],]     
                s = len(RNA.seq)
                p = s-probe_length+1 # number of probes
                for b in range(p):
                    #TDPCG=[TDP_H,TDP_S,CG] TDP_CG[G_ID][TS_ID] = TDPCG #TDP_CG=[TDP_H,TDP_S,CG]
                    start = b
                    end = b + probe_length
                    probe = RNA.seq[start:end] # this is the RC of final probe
                    probe_with_full_oh = overhang_5_full+reverse_complement(probe)+overhang_3_full

                    #check band seqs with overhang
                    ban = False
                    for st in range(ban_seqs_per_probe):
                        probe_with_oh = oh5+probe+oh3 # oh are right with above coad 
                        to_check = probe_with_oh[st:(st+ban_seqs_len)]
                        if ban_seqs_dict.get(to_check) == None:
                            ban = False
                        else:
                            ban = True
                        if ban == True:
                          # if verbose == True:
                          #     print('ban found')
                          continue
                    if ban == True:
                            continue
                        
                        
                    #check band list with full overhangs
                    ban = False
                    for bans in ban_list:
                        if bans in probe_with_full_oh:
                            ban = True
                        if ban == True:
                            # if verbose == True:
                            #     print('ban found')
                            continue
                    if ban == True:
                            continue    
                    
                    # calulate global position
                    n_exon = len(RNA.in_exon_index_start) - 1 # number of exons 
                    if start >= RNA.in_exon_index_start[-1]: # check if in last exon 
                        s_exon = n_exon
                        e_exon = n_exon
                    # find which exon probe is in
                    else:     
                        s_exon = 'x'
                        e_exon = 'x'
                        for e, start_on_TS in enumerate(RNA.in_exon_index_start):
                            # check if in same exon 
                            if start >= start_on_TS and end < RNA.in_exon_index_start[e+1]: 
                                s_exon = e
                                e_exon = e
                                break
                            # check if in span 2 exons
                            elif start >= start_on_TS and start < RNA.in_exon_index_start[e+1]:
                                s_exon = e 
                                e_exon = e+1
                                break

                    strand = transcriptome.G[RNA.gene_ID].strand
                    if strand == '1':
                        if s_exon != e_exon: # if start and end in diffrent exon
                            start_global = [(start - RNA.in_exon_index_start[s_exon]) + RNA.ex_exon_index_start[s_exon], RNA.ex_exon_index_start[s_exon+1]]
                            end_global = [RNA.ex_exon_index_stop[s_exon],(end - RNA.in_exon_index_start[e_exon]) + RNA.ex_exon_index_start[e_exon]]
                        else: # if in same exons 
                            start_global = [(start - RNA.in_exon_index_start[s_exon]) + RNA.ex_exon_index_start[s_exon]] 
                            end_global = [(end - RNA.in_exon_index_start[e_exon]) + RNA.ex_exon_index_start[e_exon]] 
                    elif strand =='-1':
                        if s_exon != e_exon: # if start and end in diffrent exon
                            start_global = [RNA.ex_exon_index_stop[s_exon] + (RNA.in_exon_index_stop[s_exon] - start), RNA.ex_exon_index_start[s_exon+1]]
                            end_global = [RNA.ex_exon_index_stop[s_exon], RNA.ex_exon_index_start[e_exon] - (end - RNA.in_exon_index_start[e_exon]) ]
                        else: # if in same exons 
                            start_global = [RNA.ex_exon_index_start[s_exon]  - (start - RNA.in_exon_index_start[s_exon]) ] 
                            end_global = [RNA.ex_exon_index_start[e_exon] - (end - RNA.in_exon_index_start[e_exon])] 

                    #calulate Tm and CG
                    H = sum(thurmo[TS_ID][0][start:end+1])
                    S = sum(thurmo[TS_ID][1][start:end+1])
                    CG = sum(thurmo[TS_ID][2][start:end+1])/(probe_length)
                    # add ending comp for AT
                    if start == 0 :
                        fiveprime = 0
                    elif RNA.seq[start-1] in ('A','T'):
                        fiveprime = 1
                    else:
                        fiveprime = 0
                    if end == s:
                        threeprime = 0
                    elif end+1 >= s:
                        threeprime = 0
                    elif RNA.seq[end+1] in ('A','T'):
                        threeprime = 1
                    else:
                        threeprime = 0
                    H = H + 0.2 + (2.2*fiveprime) + (2.2*threeprime)
                    S = S -5.7 + (6.9*fiveprime) + (6.9*threeprime)
                    HS=[H,S]
                    #get Iso_index and gene_index
                    if iso_form == True:
                        probe_iso_index = 0
                        probe_gene_index = 0
                        for sta in range(OT_seqs_per_probe):# [iso_pen,iso_index,gene_index]
                            seq = probe[sta:(sta+OT_length)]
                            seq = DNA2int(seq)
                            index = ISO_dict[TS_ID].get(seq,[1,1])
                            probe_iso_index += index[0]
                            probe_gene_index += index[1]
                            # normilize sums of probe_iso_index and  probe_gene_index
                        probe_iso_index = probe_iso_index/OT_seqs_per_probe
                        probe_gene_index = probe_gene_index/OT_seqs_per_probe
                            # add to tanscript specfice probe list 
                            #[[seq],[Tm],[CG],[Iso_index],[gene_index],[posiotional_index]
                        probe_data[0].append(DNA2int(reverse_complement(probe)))
                        probe_count+=1
                        probe_data[1].append(HS)
                        probe_data[2].append(CG)
                        probe_data[3].append(probe_iso_index)
                        probe_data[4].append(probe_gene_index)
                        probe_data[7].append(start_global)# these are not [] un like below 
                        probe_data[8].append(end_global)
                        
                    else:
                        probe_gene_index = 0
                        for sta in range(OT_seqs_per_probe):# [iso_pen,iso_index,gene_index]
                            seq = probe[sta:(sta+OT_length)]
                            seq = DNA2int(seq)
                            index = ISO_dict[TS_ID].get(seq,[1,1])
                            probe_gene_index += index[1]
                            # normilize sums of probe_gene_index
                        probe_gene_index = probe_gene_index/OT_seqs_per_probe
                        # add to gene probe list 
                        probe_data[0].append(DNA2int(reverse_complement(probe)))
                        probe_count+=1
                        probe_data[1].append(HS)
                        probe_data[2].append(CG)
                        probe_data[3].append(0)
                        probe_data[4].append(probe_gene_index)
                        probe_data[7].append(start_global)# these are not [] un like below 
                        probe_data[8].append(end_global)

                if exon_only == True:
                    if iso_form == True:
                        save_pickle(probe_data,gene_path,TS_ID)
                #update all 
                    all_obj = ALL
                    all_obj[0] += probe_data[0]
                    all_obj[1] += probe_data[1]
                    all_obj[2] += probe_data[2]
                    all_obj[3] += probe_data[3]
                    all_obj[4] += probe_data[4]
                    all_obj[7] += probe_data[7]
                    all_obj[8] += probe_data[8]
                    save_pickle(ALL, gene_path,'all')
                    continue
                
                #add probes for non-spliced genes
                gene = transcriptome.G[G_ID]
                if strand == '1'or strand != '-1':
                  unspliced_seq = gene.seq[RNA.start-gene.start:RNA.stop-gene.start]
                else:
                  unspliced_seq = gene.seq[abs(RNA.start-gene.stop):abs(RNA.stop-gene.stop)]
                s = len(unspliced_seq)
                p = s-probe_length+1 # number of probes
                for b in range(p):
                    start = b
                    end = b + probe_length
                    probe = unspliced_seq[start:end]
                    probe_with_full_oh = overhang_5_full+reverse_complement(probe)+overhang_3_full

                    #check band seqs with overhang
                    ban = False
                    for st in range(ban_seqs_per_probe):
                        probe_with_oh = oh5+probe+oh3
                        to_check = probe_with_oh[st:(st+ban_seqs_len)]
                        if ban_seqs_dict.get(to_check) == None:
                            ban = False
                        else:
                            ban = True
                        if ban == True:
                          continue
                    if ban == True:
                            continue
                    
                    #check band list with full overhangs
                    ban = False
                    for bans in ban_list:
                        if bans in probe_with_full_oh:
                            ban = True
                        if ban == True:
                            # if verbose == True:
                            #     print('ban found')
                            continue
                    if ban == True:
                            continue

                    # calulate global position
                    if strand == '1'or strand != '-1':
                      global_start = RNA.start + b
                      global_stop = global_start + probe_length
                    else:
                      global_start = RNA.start - b
                      global_stop = global_start - probe_length

                    #calulate Tm and CG
                    if strand == '1'or strand != '-1':
                      g_start = gene.start
                    else:
                      g_start = gene.stop
                    TDP_start = abs(global_start-g_start)
                    TDP_stop = TDP_start + probe_length
                    TDP_ID = 'NSG'+ G_ID

                    H = sum(thurmo[TDP_ID][0][TDP_start:TDP_stop+1])
                    S = sum(thurmo[TDP_ID][1][TDP_start:TDP_stop+1])
                    CG = sum(thurmo[TDP_ID][2][TDP_start:TDP_stop+1])/(probe_length)
                    # add ending comp for AT
                    if start == 0 :
                        fiveprime = 0
                    elif unspliced_seq[start-1] in ('A','T'):
                        fiveprime = 1
                    else:
                        fiveprime = 0
                    if end == s:
                        threeprime = 0
                    elif end+1 >= s:
                        threeprime = 0
                    elif unspliced_seq[end+1] in ('A','T'):
                        threeprime = 1
                    else:
                        threeprime = 0
                    H = H + 0.2 + (2.2*fiveprime) + (2.2*threeprime)
                    S = S -5.7 + (6.9*fiveprime) + (6.9*threeprime)
                    HS=[H,S]
                    #get Iso_index and gene_index
                    if iso_form == True:
                        probe_iso_index = 0
                        probe_gene_index = 0
                        for sta in range(OT_seqs_per_probe):# [iso_pen,iso_index,gene_index]
                            seq = probe[sta:(sta+OT_length)]
                            seq = DNA2int(seq)
                            index = ISO_dict[TS_ID].get(seq,[1,1])
                            probe_iso_index += index[0]
                            probe_gene_index += index[1]
                            # normilize sums of probe_iso_index and  probe_gene_index
                        probe_iso_index = probe_iso_index/OT_seqs_per_probe
                        probe_gene_index = probe_gene_index/OT_seqs_per_probe
                            # add to tanscript specfice probe list 
                            #[[seq],[Tm],[CG],[Iso_index],[gene_index],[posiotional_index]
                        probe_data[0].append(DNA2int(reverse_complement(probe)))
                        probe_count+=1
                        probe_data[1].append(HS)
                        probe_data[2].append(CG)
                        probe_data[3].append(probe_iso_index)
                        probe_data[4].append(probe_gene_index)
                        probe_data[7].append([global_start])
                        probe_data[8].append([global_stop])
                    else:
                        probe_gene_index = 0
                        for sta in range(OT_seqs_per_probe):# [iso_pen,iso_index,gene_index]
                            seq = probe[sta:(sta+OT_length)]
                            seq = DNA2int(seq)
                            index = ISO_dict[TS_ID][seq]
                            probe_gene_index += index[1]
                            # normilize sums of probe_gene_index
                        probe_gene_index = probe_gene_index/OT_seqs_per_probe
                        # add to gene probe list 
                        probe_data[0].append(DNA2int(reverse_complement(probe)))
                        probe_count+=1
                        probe_data[1].append(HS)
                        probe_data[2].append(CG)
                        probe_data[3].append(0)
                        probe_data[4].append(probe_gene_index)
                        probe_data[7].append([global_start])
                        probe_data[8].append([global_stop])

                if iso_form == True:
                    save_pickle(probe_data,gene_path,TS_ID)
        
                #update all 
                all_obj = ALL
                all_obj[0] += probe_data[0]
                all_obj[1] += probe_data[1]
                all_obj[2] += probe_data[2]
                all_obj[3] += probe_data[3]
                all_obj[4] += probe_data[4]
                all_obj[7] += probe_data[7]
                all_obj[8] += probe_data[8]
  
            save_pickle(ALL, gene_path,'all')
            return(probe_count)

        # frag data for multicore 
        GID_ = list(transcriptome.G.keys())
        genes_n = len(GID_)
        probe_length_ = [probe_length] * genes_n
        db_path_ = [db_path] * genes_n
        OT_length_ = [OT_length] * genes_n
        oh3_ = [oh3] * genes_n
        oh5_ = [oh5] * genes_n
        ban_seqs_dict_ = [ban_seqs_dict] * genes_n
        ban_seqs_per_probe_ = [ban_seqs_per_probe] * genes_n
        ban_seqs_len_ = [ban_seqs_len] * genes_n
        TDP_CG_db_path_ = [TDP_CG_db_path] * genes_n
        iso_form_ = [iso_form] * genes_n
        verbose_ = [verbose] * genes_n
        exon_only_ = [exon_only] * genes_n
        ban_list_ = [ban_list] * genes_n
        overhang_5_full_ = [overhang_5_full] * genes_n
        overhang_3_full_ = [overhang_3_full] * genes_n
        
        # set up threads/cores 
        if threads == -1: 
            threads = multiprocessing.cpu_count()
        pool = ThreadPool(threads)
        
        # run calulations on mulitiple cores 
        results = pool.starmap(GetProbes,zip(GID_,
                      probe_length_,
                      db_path_,
                      OT_length_,
                      oh3_,
                      oh5_,
                      ban_seqs_dict_,
                      ban_seqs_per_probe_,
                      ban_seqs_len_,
                      TDP_CG_db_path_,
                      iso_form_,
                      verbose_,
                      exon_only_,
                      ban_list_,
                      overhang_5_full_,
                      overhang_3_full_),1)#
        pool.close() 
        pool.join()
        
        probe_count = sum(results)
        
        # calulate time to compleate making probes 
        tend = time.time()
        time_total = tend - tstart
        if verbose == True:
            print(str(round((time_total/60),2))+' min to make ' + str(probe_count)+ ' probes')
            
        # save a metadata file 
        ban_seqs_used = False
        ban_seqs_path = 'None'
        if ban_seqs_dict != None:
            ban_seqs_used = True
            save_pickle(ban_seqs_dict, db_path, 'ban_seqs_dict')
            ban_seqs_path = os.path.normpath(os.path.join(db_path,'ban_seqs_dict.pk'))
        else:
            ban_seqs_path = ''
            
        ban_list_used = False
        if ban_list != []:
            ban_list_used = True
            save_pickle(ban_list, db_path, 'ban_seqs_list')
            ban_seqs_list_path = os.path.normpath(os.path.join(db_path,'ban_seqs_list.pk'))
        else: 
            ban_seqs_list_path = ''
        
        meta_path = os.path.normpath(os.path.join(db_path,'probe_metadata.txt'))
        meta_file = open(meta_path, "w") 
        to_write = ['transcriptome cell type: '+ transcriptome.cell_type,
                    'thermodynamic properties: '+transcriptome.TDP_db_path,
                    'off target properties: '+transcriptome.OT_db_path,
                    'probe_length: '+str(probe_length),
                    'probe count: '+str(probe_count),
                    'probe_db_out_path: '+db_path,
                    '5 prime overhang: '+str(overhang_5),
                    '3 prime overhang: '+str(overhang_3),
                    'isoform included: '+str(iso_form),
                    'ban_seqs_dict used: '+str(ban_seqs_used),
                    'ban_seqs_dict path: '+str(ban_seqs_path),
                    'ban_list used: '+str(ban_list_used),
                    'ban_list path: '+str(ban_seqs_list_path),
                    'split probes: False',
                    'time to make probes (min): '+str(round((time_total/60),2)),
                    'cores used: '+str(threads)]
    
        for line in to_write:
          meta_file.write(line)
          meta_file.write("\n")
        meta_file.close() 

    if verbose == True:
        print('finished making probes')
    gc.collect()    
    return(db_path)

# =============================================================================

def probe_stats(probe_sets,file_path):
    Tm = []
    probe_iso_index = []
    probe_gene_index = []
    if type(probe_sets) == dict: 
        all_gene_ID = list(probe_sets.keys())
        type_ = dict
    elif os.path.exists(probe_sets) == False:
        print('probe_sets is not a dictinary or the file path does not exist')
        return(None)
    elif os.path.exists(probe_sets) == True:
        all_gene_ID = [f for f in os.listdir(probe_sets) if not f.startswith('.')]
        all_gene_ID.remove('probe_metadata.txt')
        type_ = 'db' 
    gene_ID=all_gene_ID
    for gene in gene_ID:
       if type_ == dict:
            d = probe_sets[gene]['all'] #[[seq],[Tm],[CG],[Iso_index],[gene_index],[posiotional_index]]
       elif type_ == 'db':
            d = load_pickle(os.path.normpath(os.path.join(probe_sets,gene,'all')))
       Tm.extend(d[1])
       probe_iso_index.extend(d[3])
       probe_gene_index.extend(d[4])
    save_pickle(Tm,file_path,'Tm.pk')
    save_pickle(probe_iso_index,file_path,'iso_index.pk')
    save_pickle(probe_gene_index ,file_path,'gene_index.pk')

# =============================================================================

def reverse_complement(dna,N = 'A'):
    complement = {'N':N,'n':N,'A':'T','a':'T','C':'G','c':'G','G':'C','g': 'C','T': 'A','t': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

# =============================================================================

def pick_best_probes(probe_sets,
                    output_path,
                    tm_low,
                    tm_high,
                    formamide_fraction,
                    monovalentSalt,
                    probeConc,
                    medadata_out_path = None,
                    output_file_name = 'filtered_probes', 
                    CG_low=0,
                    CG_high=1,
                    iso_low=0,
                    iso_high=1.1,
                    min_probe_spaceing = 0,
                    max_probes=float('inf'),
                    cross_hyb = 15,
                    splice = None, 
                    gene_ID = 'all',
                    Chem_corrected_Tm = False,
                    iso = False,
                    verbose = False,
                    ligase = 'splintR',
                    cross = dict(),
                    exon_intron_only = False,
                    threads = 4):
    
    """ filters and orders probes, returns dict of panda dataframs 
    
    pick_best_probes(probe_sets, medadata_out_path, tm_low, tm_high, molar_formamide,
                    monovalentSalt, probeConc, CG_low=0, CG_high=1, iso_low=0,
                    iso_high=1.1, min_probe_spaceing = 0, max_probes=30, cross_hyb = 15,
                    splice = None, gene_ID = 'all', Chem_corrected_Tm = False,
                    iso = False, verbose = False, cross = dict())
        probe_sets: str file path or dict of probes 
        medadata_out_path: str file path to store metadata file
        output_path: str file path to store dataframs of probes 
        tm_low: int or float lower bound for metling temp of probe 
        tm_high: int or float upper bound for metling temp of probe 
        formamide_fraction: int or float concentration of formamide in fraction 
        monovalentSalt: int or float mononvalent salt concentration in molar 0.79 = 2XSSC
        probeConc: int or float concentration of probes in molar
        CG_low: int or float 0<= CG_low <=1 lower bound for CG content of probe 
        CG_high: int or float 0<= CG_high <=1 upper bound for CG content of probe 
        iso_low: int or float 0<= iso_low <=1 lower bound for isoform index of probe
        iso_high: int or float 0<= iso_high <=1 upper bound for isoform index of probe
        min_probe_spaceing: int min number of bases between probes
        max_probes: int or float max number of probes to keep per gene or transcript 
        cross_hyb: int length of substring to check for cross hybridization  
        splice: ('spliced' or 'unspliced' or 'both' , int) probe kind, input mininmum overlap on splice site (no longer suported. only exactly on splice site),
            'both' indicates probe pairs that lay on exon/inton('unspliced') boundry and exon/exon('spliced') boundry
        gene_ID: 'all' or list of gene IDs to make probes for 
        Chem_corrected_Tm: bool calulate a chemicaly corrected Tm based on molar_formamide and monovalentSalt
        iso: bool prioritize isofrom index in picking probes 
        verbose: bool print progress of calulation every 20s 
        cross: dict probe:0 all probes to be checked against and added to as probes are made 
        exon_intron_only: only keep probes taht are not splicing probe/not across exons 
    
    """
    ################### check inputs #########################
    if type(probe_sets) != dict and type(probe_sets) != str:
        raise TypeError('probe_sets must be of type dict or str not '+ str(type(probe_sets)))
    
    if type(probe_sets) != dict and os.path.exists(os.path.normpath(probe_sets)) == False:
        raise Exception('probe_sets '+os.path.normpath(probe_sets)+' does not exist')
    
    if type(output_path) != str:
        raise TypeError('output_path must be of type str not '+ str(type(output_path)))
        
    if os.path.exists(os.path.normpath(output_path)) == False:
        raise Exception('output_path '+os.path.normpath(output_path)+' does not exist')
    
    if medadata_out_path == None or medadata_out_path == '':
        medadata_out_path = output_path
    
    if type(medadata_out_path) != str:
        raise TypeError('medadata_out_path must be of type str not '+ str(type(medadata_out_path)))
    
    if os.path.exists(os.path.normpath(medadata_out_path)) == False:
        raise Exception('medadata_out_path '+os.path.normpath(medadata_out_path)+' does not exist')
    
    if type(tm_low) != int and type(tm_low) != float:
        raise TypeError('tm_low must be of type int or float not '+ str(type(tm_low)))
    
    if type(tm_high) != int and type(tm_high) != float:
        raise TypeError('tm_high must be of type int or float not '+ str(type(tm_high)))
        
    if type(formamide_fraction) != int and type(formamide_fraction) != float:
        raise TypeError('molar_formamide must be of type int or float not '+ str(type(formamide_fraction)))

    if type(monovalentSalt) != int and type(monovalentSalt) != float:
        raise TypeError('monovalentSalt must be of type int or float not '+ str(type(monovalentSalt)))
    
    if type(probeConc) != int and type(probeConc) != float:
        raise TypeError('probeConc must be of type int or float not '+ str(type(probeConc)))
    
    if type(CG_low) != int and type(CG_low) != float:
        raise TypeError('CG_low must be of type int or float not '+ str(type(CG_low)))
    
    if CG_low < 0 or CG_low > 1:
        raise Exception('CG_low must be between 0 and 1 not'+ str(CG_low))
    
    if type(CG_high) != int and type(CG_high) != float:
        raise TypeError('CG_high must be of type int or float not '+ str(type(CG_high)))
    
    
    if CG_high < 0 or CG_high > 1:
        raise Exception('CG_high must be between 0 and 1 not'+ str(CG_high))
    
    #coment out for now to disregard
    #if type(iso_low) != int and type(iso_low) != float:
    #    raise TypeError('iso_low must be of type int or float not '+ str(type(iso_low)))
    #
    #if iso_low < 0 or iso_low > 1:
    #    raise Exception('iso_low must be between 0 and 1 not'+ str(iso_low))
    #
    #if type(iso_high) != int and type(iso_high) != float:
    #    raise TypeError('iso_high must be of type int or float not '+ str(type(iso_high)))
    #
    #if iso_high < 0 or iso_high > 1:
    #    raise Exception('iso_high must be between 0 and 1 not'+ str(iso_high))

    if type(min_probe_spaceing) != int:
        raise TypeError('min_probe_spaceing must be of type int not '+ str(type(min_probe_spaceing)))
    
    if type(max_probes) != int and type(max_probes) != float:
        raise TypeError('max_probes must be of type int not '+ str(type(max_probes)))

    if type(cross_hyb) != int:
        raise TypeError('cross_hyb must be of type int not '+ str(type(cross_hyb)))
        
    if type(splice) != type(None):
        if type(splice) != tuple and type(splice) != type(None) or len(splice) != 2 or splice[0] not in ['spliced','unspliced','both'] or type(splice[1]) != int :
            raise Exception('splice must be of type tuple in the form (\'spliced\' or \'unspliced\', int) '+ str(type(splice)))
    
    if type(gene_ID) != list and gene_ID != 'all':
        raise TypeError('gene_ID must be \'all\' or of type list not ' + str(type(gene_ID)))
      
    if type(Chem_corrected_Tm) != bool:
        raise TypeError('Chem_corrected_Tm must be of type bool not '+ str(type(Chem_corrected_Tm)))
    
    if type(iso) != bool:
        raise TypeError('iso must be of type bool not '+ str(type(iso)))
        
    if type(verbose) != bool:
        raise TypeError('verbose must be of type bool not '+ str(type(verbose)))
    
    if type(cross) != dict:
        raise TypeError('cross must be of type dict not '+ str(type(cross)))
        
    if type(exon_intron_only) != bool:
        raise TypeError('exon_intron_only must be of type bool not '+ str(type(exon_intron_only)))
        
    ##########################################################
    # conver formimide ration to molar 
    molar_formamide = formamide_fraction*(1133/45.041)
    
    # make output directory 
    clear = False
    out_path = os.path.normpath(os.path.join(output_path,output_file_name))
    if os.path.exists(out_path) == False:
            os.makedirs(out_path)
            clear = True
    count = 0    
    while True:    
        if clear==True:
            break
        if os.path.exists(out_path) == False:
            os.makedirs(out_path)
            break
        else:
            out_path = os.path.normpath(os.path.join(output_path,output_file_name+'_'+str(count)))
            count+=1
    
    if medadata_out_path ==  output_path:
        medadata_out_path = out_path
    
    # check if ther is a max probe directory 
    max_probes_dict = False
    if type(max_probes) == dict:
        max_probes_dict = True
    tstart = time.time()
   
    # right in middle of high and low Tm 
    ideal_t = ((tm_high - tm_low)/2) + tm_low
    
    #best_probes = dict()
    
    if type(probe_sets) == dict: 
        all_gene_ID = list(probe_sets.keys())
        type_ = dict

    elif os.path.exists(os.path.normpath(probe_sets)) == True:
        all_gene_ID = [f for f in os.listdir(probe_sets) if not f.startswith('.')]
        all_gene_ID.remove('probe_metadata.txt')
        if 'ban_seqs_dict.pk' in all_gene_ID:
            all_gene_ID.remove('ban_seqs_dict.pk')
        if 'ban_seqs_list.pk' in all_gene_ID:
            all_gene_ID.remove('ban_seqs_list.pk')
        type_ = 'db' 
        
    if gene_ID == 'all':
        gene_ID=all_gene_ID
    
    #get over hangs from meta, there is a better way to do this
    meta_path = os.path.normpath(os.path.join(probe_sets,'probe_metadata.txt'))
    meta_file = open(meta_path, "r")
    meta_file_data = meta_file.readlines()
    meta_file_data_list = []
    for line in meta_file_data:
      meta_file_data_list.append(line.replace('\n',''))
    meta_file.close()
    five = meta_file_data_list[6]
    three = meta_file_data_list[7]
    five_OH = five[find_nth_n_str(meta_file_data_list[5],':',1)+2:] 
    three_OH = three[find_nth_n_str(meta_file_data_list[6],':',1)+2:]
    
    #check if split probes
    split_text = meta_file_data_list[13][-5:]
    if 'T' in split_text:
        split = True
    elif 'F' in split_text:
        split = False
    else:
        split = False
        
    OHs =[]
    if five_OH == 'None' or five_OH == ' ':
        five_OH_add = ''
    else:
        if len(five_OH) >= cross_hyb:
          five_OH_add = five_OH[-(cross_hyb-1):]
          OHs.append(five_OH)
        else:
          five_OH_add = five_OH
    
    if three_OH == 'None' or three_OH == ' ':
        three_OH_add = ''
    else: 
        if len(three_OH) >= cross_hyb:
          three_OH_add = three_OH[0:(cross_hyb-1)]
          OHs.append(three_OH)
        else: 
          three_OH_add = three_OH
        
        
     # reverse complements of probes chuncks of length cross_hyb value
    # add overhangs to cross
    if len(OHs)>0:
      for seq in OHs:
          for b in range(len(seq)-cross_hyb):
              chunk = seq[b:b+cross_hyb]
              #print(chunk)
              cross[reverse_complement(chunk)] = 0

    num_gene = len(gene_ID)
    count = 0
    t1=0
    
    #for gene in gene_ID:
    def filter_probes(gene,
                     output_path,
                     tm_low,
                     tm_high,
                     formamide_fraction,
                     monovalentSalt,
                     probeConc,
                     output_file_name , 
                     CG_low,
                     CG_high,
                     iso_low,
                     iso_high,
                     min_probe_spaceing,
                     max_probes,
                     cross_hyb,
                     splice, 
                     gene_ID,
                     Chem_corrected_Tm,
                     iso,
                     verbose,
                     ligase,
                     exon_intron_only,
                     cross):
    
       #cross = {}
       #count+=1
       if type_ == dict:
            d = probe_sets[gene]['all'] #[[seq],[HS],[CG],[Iso_index],[gene_index],gene symbol, [chr], [start], [stop posiotional_index_on ch],]
       elif type_ == 'db':
            d = load_pickle(os.path.normpath(os.path.join(probe_sets,gene,'all.pk')))
            print(gene+' loaded')
       if d[0]==[]:
           print('d[0]==[]')
           return('')
      
       if split == True:
           probe_length = len(int2DNA(d[0][0][0]))*2
       else:
           probe_length = len(int2DNA(d[0][0]))
       cross_hyb_per_probe = probe_length - cross_hyb +1
       #calulate Tm from H and S 
       HS_list = d[1]
       Tml_list = []
       Tmr_list = []
       for i,HS in enumerate(HS_list):
         if split == True:
             Hl=HS[0][0]
             Sl=HS[0][1]
             Hr=HS[1][0]
             Sr=HS[1][1]
             
             CGl = d[2][i][0]
             CGr = d[2][i][1]
             
             # add salt correction
             Sl = Sl + 0.368*(probe_length-1)*m.log(monovalentSalt); #might need to correct for salt and 
             Sr = Sr + 0.368*(probe_length-1)*m.log(monovalentSalt); #might need to correct for salt and fomimide
             
             Tml = ((Hl*1000) / (Sl + (1.9872 * m.log(probeConc)))) - 273.15 # this is mostly correct 
             Tmr = ((Hr*1000) / (Sr + (1.9872 * m.log(probeConc)))) - 273.15 # this is mostly correct 
             
             if Chem_corrected_Tm == True:
                 Tml += (0.453 * CGl - 2.88) * molar_formamide
                 Tmr += (0.453 * CGr - 2.88) * molar_formamide
             Tml_list.append(Tml)
             Tmr_list.append(Tmr)
             
         else:
             H=HS[0]
             S=HS[1]
             CG = d[2][i]
             # add salt correction
             S = S + 0.368*(probe_length-1)*m.log(monovalentSalt); #might need to correct for salt and fomimide
             Tm = ((H*1000) / (S + (1.9872 * m.log(probeConc)))) - 273.15 # this is mostly correct 
             if Chem_corrected_Tm == True:
                 Tm += (0.453 * CG - 2.88) * molar_formamide
             Tml_list.append(Tm)
             Tmr_list.append(Tm)
             
       cgl = []
       cgr = []
       if split == True:
           for cg in d[2]:
               cgl.append(cg[0])
               cgr.append(cg[1])
       else:
           cgl = d[2]
           cgr = d[2]
           
       Isol = []
       Isor = []
       Isoc = []
       GeneSpesl = []
       GeneSpesr = []
       GeneSpesc = []
       for iso_,GenS in zip(d[3],d[4]):
           if split == True:
               Isol.append(iso_[0])
               Isor.append(iso_[1])
               Isoc.append(iso_[2])
               GeneSpesl.append(GenS[0])
               GeneSpesr.append(GenS[1])
               GeneSpesc.append(GenS[2])
           else:
               Isol.append(iso_)
               Isor.append(iso_)
               Isoc.append(iso_)
               GeneSpesl.append(GenS)
               GeneSpesr.append(GenS)
               GeneSpesc.append(GenS)
       
       DNA = []
       fullDNA = []
       CGinRight = []
       ########################### add another filter for left probes G at end!!!#######################
       # this paper says it has lower efficentcy (Efficient DNA ligation in DNA–RNA hybrid helices by Chlorella virus DNA ligase 
        #Gregory J. S. Lohman, Yinhua Zhang, Alexander M. Zhelkovsky, Eric J. Cantor, Thomas C. Evans, Jr)

### Also from same paper a dT/(phos)dA ligation junction had best reaction rate #######################
       for Seq in d[0]:
           if split == True:
               leftSeq = int2DNA(Seq[0])
               rihgtSeq = int2DNA(Seq[1])
               fullDNA.append(leftSeq+rihgtSeq)
               DNA.append([leftSeq,rihgtSeq])
               if 'C' in rihgtSeq[:2]:
                   CGinRight.append('T')
               elif 'G' in rihgtSeq[:2] or 'G' in leftSeq[-1:]:
                   CGinRight.append('T')
               else:
                   CGinRight.append('F')
                   
               ##################################################
   
    # The enzyme tolerates all base pair combinations at the ligation junction,
    # but is partially inhibited by dC/G and dG/C base pairs at the donor (phosphorylated)
    # side ligation junction, particularly when the +2 base was also a C/G base pair.
     #              not C/G--->d?d?<---not C/G
     # 5' dNdNdNdNdNdNdN    (P)dNdNdNdNdNdNdN 3'
     # 3' rNrNrNrNrNrNrN-------rNrNrNrNrNrNrN 5'
               
           else: 
               DNA.append(int2DNA(Seq))
               fullDNA.append(int2DNA(Seq))
               CGinRight.append('F')
       
       if ligase == 'splitnR' and split == True:
           ligaseCG = True
                
       # create list for first stop position and second stop position 
       start1 = []
       stop2 = []
    
       for st,sto in zip(d[7],d[8]):
         if type(st) == list:
             start1.append(st[0])
             if len(sto) > 1:
                stop2.append(sto[1])
             else:
                stop2.append(sto[0])
            
         elif type(st) == int:
             start1.append(st)
             stop2.append(sto)
            
       #make data frame of all probes 
       probes_per_set = len(d[0])
       n_gene_name = [gene] * probes_per_set
       d1 = {'gene_ID': n_gene_name ,'gene_sym':[d[5]] * probes_per_set,
             'seq': DNA,'Tm_l':Tml_list,'Tm_r':Tmr_list,'CGl':cgl,'CGr':cgr,
             'isoform_specificity_index_l':Isol,'isoform_specificity_index_r':Isor,
             'isoform_specificity_index_c':Isoc,'gene_specificity_index_l':GeneSpesl,
             'gene_specificity_index_r':GeneSpesr,'gene_specificity_index_c':GeneSpesc,
             'chr':d[6],'start':d[7],'stop':d[8],'start1':start1,'stop2':stop2,'seqfull': fullDNA, 'GCinRight':CGinRight}

       df = pd.DataFrame(d1, columns=['gene_ID','gene_sym', 'seq', 'Tm_l',
                                      'Tm_r','CGl','CGr','isoform_specificity_index_l',
                                      'isoform_specificity_index_r',
                                      'isoform_specificity_index_c',
                                      'gene_specificity_index_l',
                                      'gene_specificity_index_r',
                                      'gene_specificity_index_c',
                                      'chr','start','stop','start1','stop2','seqfull','GCinRight'])
       # imply strand
       if start1[0]<stop2[0]:
           strand = '1'
       else:
           strand = '-1'
        
       # spice site specficity prioriety
       if splice != None and splice[0] in ['spliced','unspliced','both']:
           spliced_indexes = []
           unspliced_starts = []
           unspliced_stops2 = []
           unspliced_indexes = []
           start_idexes = df['start'].tolist()
           stop_indexes = df['stop'].tolist()
           half_probe_len = probe_length/2
          
           # probes on bounderys will have 2 genomic starts, enforce probes exactly on bondry of spliec sites, check both probes
           for i, start in enumerate(start_idexes):
                  if strand == '1':
                      if len(start) == 2  and abs(start[0]-stop_indexes[i][0]) == half_probe_len and abs(start[1]-stop_indexes[i][1])+1 == half_probe_len:
                          spliced_indexes.append(i+1) # forward strand
                          unspliced_starts.append(start[0]+1)# forward strand
                          unspliced_stops2.append(stop_indexes[i][1]+1) # keep sedond stop index # forward strand
                  if strand == '-1':
                      if len(start) == 2  and abs(start[0]-stop_indexes[i][0]) == half_probe_len and abs(start[1]-stop_indexes[i][1]) == half_probe_len:
                          spliced_indexes.append(i) # 
                          unspliced_starts.append(start[0]-1)# the added base makes unsplied probes correct
                          unspliced_stops2.append(stop_indexes[i][1]) # keep sedond stop index # 
                    
                    
           if splice[0] in ['spliced','both'] :
               spliced_df = df.iloc[spliced_indexes]
                  
           if splice[0] in ['unspliced','both']:
               start1_ = np.array(start1)
               stop2_ = np.array(stop2)
           
               for start in unspliced_starts:
                   indices = np.where(start1_ == start)[0].tolist()
                   for num in indices: # filter for not exon/exon probes 
                       if num not in spliced_indexes:
                           unspliced_indexes.append(num)
                   
               for stop in unspliced_stops2:
                   indices = np.where(stop2_ == stop)[0].tolist()
                   for num in indices: # filter for not exon/exon probes 
                       if num not in spliced_indexes:
                           unspliced_indexes.append(num) 
                       
               unspliced_df = df.iloc[unspliced_indexes]
           
           if splice[0] == 'spliced':
               df = spliced_df
           elif splice[0] == 'unspliced':
               df = unspliced_df
           elif splice[0] == 'both':
             df = spliced_df.append(unspliced_df,ignore_index=True)
 
           df.drop('start1', axis = 1,inplace=True)
           df.drop('stop2', axis = 1,inplace=True)
       
       # filter for probe taht dont cross splicing boundry 
       if exon_intron_only == True:
           non_spaning_probe_indexes = []
           start_idexes = df['start'].tolist()
           for i, start in enumerate(start_idexes):
               if len(start) == 1:
                   non_spaning_probe_indexes.append(i)
           non_spaning_df = df.iloc[non_spaning_probe_indexes]
           df = non_spaning_df

       # filter for Tm and GC content and GCinRight
    
       if ligase == 'splintR':
           df = df[df.GCinRight=='F'] # keep only probes with no CG in right probes, or all if not split 
        
       if split == True:
           df = df[df.Tm_l>=tm_low]
           df = df[df.Tm_l<=tm_high]
           df = df[df.CGl>=CG_low]
           df = df[df.CGl<=CG_high]
           
           df = df[df.Tm_r>=tm_low]
           df = df[df.Tm_r<=tm_high]
           df = df[df.CGr>=CG_low]
           df = df[df.CGr<=CG_high]
       else:
           df = df[df.Tm_l>=tm_low]
           df = df[df.Tm_l<=tm_high]
           df = df[df.CGl>=CG_low]
           df = df[df.CGl<=CG_high]
       # add cloumb for distance to ideal Tm
       d_ideal_t =[]
       for tml,tmr in zip(list(df['Tm_l']),list(df['Tm_r'])):
           if split == True:
               tm = (tml+tmr)/2
           else:
               tm = tml
           d_ideal_t.append(-abs(ideal_t - tm))
           
       df['d_ideal_t'] = d_ideal_t
       # check for empty dataframs 
       if df.empty:
           
           d1 = {'gene_ID': gene ,'gene_sym':d[5], 'seq': 'no_probes_for_this_gene',
             'Tm_l':None,'Tm_r':None,'CGl':None,'CGr':None,
             'isoform_specificity_index_l':None,'isoform_specificity_index_r':None,
             'isoform_specificity_index_c':None,'gene_specificity_index_l':None,
             'gene_specificity_index_r':None,'gene_specificity_index_c':None,
             'chr':None,'start':None,'stop':None,'seqfull':None}

           df1 = pd.DataFrame(d1, columns=['gene_ID','gene_sym', 'seq', 'Tm_l','Tm_r','CGl','CGr',
                                      'isoform_specificity_index_l',
                                      'isoform_specificity_index_r',
                                      'isoform_specificity_index_c',
                                      'gene_specificity_index_l',
                                      'gene_specificity_index_r',
                                      'gene_specificity_index_c',
                                      'chr','start','stop','seqfull'],index=[0])
           
           # d1 = {'gene_ID':gene ,'gene_sym':d[5], 'seq': 'no_probes_for_this_gene',
           #       'Tm':None,'CG':None,'isoform_specificity_index':None,
           #       'gene_specificity_index':None,'chr':None,'start':None,'stop':None}

           # df1 = pd.DataFrame(d1, columns=['gene_ID','gene_sym', 'seq', 'Tm','CG',
           #                                 'isoform_specificity_index','gene_specificity_index',
           #                                 'chr','start','stop'],index=[0])
           
           #best_probes[gene] = df1
           save_pickle(df1,out_path,gene)
           
           if verbose == True:    
               print('no probes kept for gene '+gene+ ' ('+d[5]+') because CG and/or Tm parameters too stringent')
           return('')
       
       # order df by transcriptome specficity and Tm and isotype specficity. drop all duplicates if you want isoform specfic 
       if iso == True:
           df.sort_values(by=['gene_specificity_index_c','isoform_specificity_index_c','d_ideal_t'], inplace=True, ascending=False) # order on iso specficity, should get exon boundry probes 
           #df.drop_duplicates(subset='seqfull', keep=False ,inplace=True) # drop duplicats 
           df.drop_duplicates(subset='seqfull', keep='first' ,inplace=True) # drop duplicats 
            
       else:
           df['isoform_specificity_index_c'] = [1-x for x in df['isoform_specificity_index_c']] # invert isofrom index to make hi number less unique to a gene 
            # this is to select for more general probes that bind many isofroms 
           df.sort_values(by=['gene_specificity_index_c','isoform_specificity_index_c','d_ideal_t'], inplace=True, ascending=False)
           df.drop_duplicates(subset='seqfull', keep='first',inplace=True) # drop keep first duplicats 
            # sqitcht the index back to hi = more isofrom specfic 
           df['isoform_specificity_index_c'] = [1-x for x in df['isoform_specificity_index_c']]
           
       df.reset_index(drop=True,inplace=True)
       # get non overlapping, check for cross hibridization in libuary
       first_index = df['start'].tolist()
       sequences = df['seqfull'].tolist()

       if len(first_index) == 0:
           #continue
           #print('len(first_index) == 0')
           return('')
        
        # find first non-crosshyb probes starting from begings of ordered list
       for i, seq in enumerate(sequences):
         keep = True
         for b in range(cross_hyb_per_probe):
           check = seq[b:b+cross_hyb]
           if cross.get(check) != None:
             keep = False
             break
         if keep == True:
           first_good_index = i
           break
            
            
       #add first non-crosshyb reverse compliment chunckes to cross
       first_good_seq = sequences[first_good_index]
       for b in range(cross_hyb_per_probe):
         chunk = first_good_seq[b:b+cross_hyb]
         cross[reverse_complement(chunk)] = 0

       #get first start and stop positions 
       first_index_s = first_index[first_good_index]
       first_index_e = (df['stop'].tolist())[first_good_index]

       if first_index_s[0] > first_index_e[0]:
           strand = '-1'
       else:
           strand = '1'
       NOLi=[[first_index_s,first_index_e]]
       indexes =[first_good_index]
       stopL = df['stop'].values.tolist()
       
       if max_probes_dict == True:
           max_probes_this_gene = max_probes.get(gene, float('inf'))
       else:
           max_probes_this_gene = max_probes
       for i, start in enumerate(df['start'].values.tolist()):
           if len(indexes)>= max_probes_this_gene:
                break
           add = [0]
           #check for cross hybridizatio
           seq = sequences[i]

           keep = True
           for b in range(cross_hyb_per_probe):
             check = seq[b:b+cross_hyb]
             if cross.get(check) != None:
               keep = False
               break
           if keep == False:
              continue

           #check for non overlaping
           for good in NOLi:
                if strand =='1':
                    s = good[0]
                    stop = good[1]
                    for n, st in enumerate(s):
                        for sta in start:
                            if sta >= (st-probe_length-min_probe_spaceing) and sta <= (stop[n]+min_probe_spaceing):
                                add.append(1)
                                break 
                            
                elif strand == '-1':
                    s = good[1]
                    stop = good[0]
                    for n, st in enumerate(s):
                        for sta in stopL[i]:
                            if sta >= (st-probe_length-min_probe_spaceing) and sta <= (stop[n]+min_probe_spaceing):
                                add.append(1)
                                break
           add = int(sum(add))                 
           if add == 0:
               NOLi.append([start,stopL[i]])
               indexes.append(i)
               #add to chunks to cross 
               for b in range(cross_hyb_per_probe):
                  chunk = seq[b:b+cross_hyb]
                  cross[reverse_complement(chunk)] = 0
                    
       # keep overlaping probes to maintian both exon/exon and intron/exon probes              
       if splice != None and splice[0] == 'both':
           df = df#print('saveing all probes')
        
       #elif exon_intron_only == True:
       #    df = df

       else:
           df = df.iloc[indexes]
       
       #best_probes[gene] = df
       save_pickle(df,out_path,gene)
       print(gene+' saved') 
        
       del df
       #gc.collect()
       #t = time.time()
       #if verbose == True and t > t1:
       #    print(str(count)+'/'+str(num_gene) + ' probe sets filtered')
       #    t1 = t + 20
    # multip thread 
    # set up threads/cores 
    if threads == -1: 
        threads = multiprocessing.cpu_count()
    pool = ThreadPool(threads)
        
    # run calulations on mulitiple cores 
    g_len = len(gene_ID)
    results = pool.starmap(filter_probes, zip(gene_ID,
                                              [output_path]*g_len,
                                              [tm_low]*g_len,
                                              [tm_high]*g_len,
                                              [formamide_fraction]*g_len,
                                              [monovalentSalt]*g_len,
                                              [probeConc]*g_len,
                                              [output_file_name]*g_len, 
                                              [CG_low]*g_len,
                                              [CG_high]*g_len,
                                              [iso_low]*g_len,
                                              [iso_high]*g_len,
                                              [min_probe_spaceing]*g_len,
                                              [max_probes]*g_len,
                                              [cross_hyb]*g_len,
                                              [splice]*g_len, 
                                              [gene_ID]*g_len,
                                              [Chem_corrected_Tm]*g_len,
                                              [iso]*g_len,
                                              [verbose]*g_len,
                                              [ligase]*g_len,
                                              [exon_intron_only]*g_len,
                                              [cross]*g_len),1)#
    
    pool.close() 
    pool.join()   
    
    if splice != None:
      splice_bool = 'True, '+splice[0]+' transcripts, minimum splice cover: '+str(splice[1]) 
    else:
      splice_bool = 'False'
    
    num = 0
    while True:
        extended_meta_path = os.path.normpath(os.path.join(medadata_out_path,'selected_probe_metadata_'+str(num)+'.txt'))
        if os.path.exists(extended_meta_path) == True:
            num +=1
        else:
            break
    #extended_meta_path = os.path.normpath(os.path.join(medadata_out_path,'selected_probe_metadata_'+str(num)+'.txt'))
    meta_file = open(extended_meta_path, "w")
    tend = time.time()
    ttotal = round((tend-tstart)/60,2)
    to_write = meta_file_data_list + ['\n'] + ['__Probe selection parameters__',
                'output probes path: '+str(out_path),                             
                'number of genes: '+str(num_gene),
                'melting temperature range: '+str(tm_low)+' < Tm < '+str(tm_high),
                'CG content range: ' + str(CG_low) + ' < CG < ' +str(CG_high),
                'formamide concentration (M): '+str(molar_formamide),
                'mononvalent salt concentration (M): '+str(monovalentSalt),
                'probe concentration (M): '+str(probeConc),
                'chemically correct melting temperature: '+str(Chem_corrected_Tm), 
                'minimum probe spacing: '+str(min_probe_spaceing),
                'maximum probe per gene: '+str(max_probes),
                'cross hybridization length: '+str(cross_hyb),
                'isoform specificity index priority: '+str(iso),
                'isoform specificity index (ISI) range: '+str(iso_low)+' < ISI < '+str(iso_high),
                'splice site probes: '+splice_bool,
                'split probes: '+ str(split),
                'exon/intron only probes: '+str(exon_intron_only),
                '\n',
                'time to select probes (min): '+str(ttotal)]

    meta_file.write("\n")
    for line in to_write:
          meta_file.write(line)
          meta_file.write("\n")
    meta_file.close()
    save_pickle(cross,out_path,'cross')

    if verbose == True:
        print('finished picking best probes')
    gc.collect()
    return(out_path)

# =============================================================================

def probes_to_bed(probe_sets, #db path or dict
                  Transcriptome,
                  file_path,
                  probe_bed_file_name,
                  score = 'gene_index', #[HS],[CG],[Iso_index],[gene_index],
                  geneome_browser  = None, # ('IGV' or 'UCSD')
                  gene_ID = all,
                  verbose = False):
    
    """ filters and orders probes, returns dict of panda dataframs """
    best_probes = dict()
    
    if type(probe_sets) == dict: 
        all_gene_ID = list(probe_sets.keys())
        type_ = dict
    elif os.path.exists(probe_sets) == False:
        print('probe_sets is not a dictinary or the file path does not exist')
        return(None)
    elif os.path.exists(probe_sets) == True:
        all_gene_ID = [f for f in os.listdir(probe_sets) if not f.startswith('.')]
        type_ = 'db' 
    if gene_ID == all:
        gene_ID=all_gene_ID
    num_gene = len(gene_ID)
    count = 0
    t1=0
    for gene in gene_ID:
       count+=1
       if type_ == dict:
            d = probe_sets[gene]['all'] #[[seq],[Tm],[CG],[Iso_index],[gene_index],gene symbol, [chr], [start], [stop posiotional_index_on ch],]
       elif type_ == 'db':
            d = load_pickle(os.path.normpath(os.path.join(probe_sets,gene,'all.pk')))    
       if d[0]==[]:
           continue
       #make data frame 
       probes_per_set = len(d[0])
       n_gene_name = [gene] * probes_per_set
       d1 = {'gene_ID': n_gene_name ,'gene_sym':[d[5]] * probes_per_set, 
             'seq': d[0],'Tm':d[1],'CG':d[2],'isoform_specificity_index':d[3],
             'gene_specificity_index':d[4],'chr':('chr'+d[6]),'start':d[7],
             'stop':d[8]}
       df = pd.DataFrame(d1, columns=['gene_ID','gene_sym', 'seq', 'Tm','CG',
                                      'isoform_specificity_index','gene_specificity_index',
                                      'chr','start','stop'])

       # check for empty dataframs 
       if df.empty:
           d1 = {'gene_ID':gene ,'gene_sym':d[5], 'seq': 'no_probes_for_this_gene',
                 'Tm':None,'CG':None,'isoform_specificity_index':None,
                 'gene_specificity_index':None,'chr':None,'start':None,'stop':None}
           df1 = pd.DataFrame(d1, columns=['gene_ID','gene_sym', 'seq', 'Tm','CG',
                                           'isoform_specificity_index','gene_specificity_index',
                                           'chr','start','stop'],index=[0])
           best_probes[gene] = df1
           if verbose == True:    
               print('no probes kept for gene '+gene+ ' ('+d[5]+') because CG and/or Tm parameters too stringent')
           continue
       best_probes[gene] = df
       t = time.time()
       if verbose == True and t > t1:
           print(str(count)+'/'+str(num_gene) + ' probe sets filtered')
           t1 = t + 20

    ids = list(best_probes.keys())
    chr = []
    start = []
    stop = []
    name =[]
    score_ = []
    strand = []
    thickStart = []
    thickEnd = []
    itemRgb = []
    blockCount = []
    blockSizes = []
    blockStarts = []
    start_end = []
    for gID in ids:
        df = best_probes[gID]
        gene_sym = df['gene_sym'].values.tolist()[0]
        probe_len = len(df.iloc[0]['seq'])
        chromosome = df['chr'].values.tolist()
        chr += chromosome
        score_ += df[score].values.tolist()
        gene_specificity_index_local = df['gene_specificity_index'].values.tolist()
        isoform_specificity_index_local = df['isoform_specificity_index'].values.tolist()
        seq = df['seq'].values.tolist()
        Tm = df['Tm'].values.tolist()
        e = df['stop'].values.tolist()
        s = df['start'].values.tolist()

        if Transcriptome.TS.get(gID) == None:
            #this_chr = Transcriptome.G[gID].chr
            this_strand = Transcriptome.G[gID].strand
        else:
            #this_chr = Transcriptome.G[Transcriptome.TS[gID].gene_ID].chr
            this_strand = Transcriptome.G[Transcriptome.TS[gID].gene_ID].strand
        #chr += [this_chr] * rep
        
        if this_strand == '1' or this_strand !='-1':
          for n,st in enumerate(s):
              start.append(st[0])
              start_end_start = st[0]
              if len(st)>1:
                block1 = abs(st[0]- e[n][0])
                if block1 == probe_len:
                  blockCount.append(1)
                  blockStarts.append(str(0))
                  blockSizes.append(probe_len)
                  stop.append(e[n][0])
                  start_end_stop = e[n][0]
                else:
                  blockStarts.append(str(0) + ',' + str(abs(st[0]-st[1])) + ',')
                  blockCount.append(2)
                  blockSizes.append(str(block1)+','+str(probe_len - block1) + ',')
                  stop.append(e[n][-1])
                  start_end_stop = e[n][-1]
              else:
                blockCount.append(1)
                blockStarts.append(str(0))
                blockSizes.append(probe_len)
                stop.append(e[n][-1])
                start_end_stop = e[n][-1]
              start_end.append(chromosome[0] + str(start_end_start)+':'+str(start_end_stop))

        elif this_strand == '-1':
          for n,st in enumerate(s):
              start.append(e[n][-1])
              start_end_start = e[n][-1]
              if len(st)>1:
                block1 = abs(e[n][-1]- st[1])
                if block1 == probe_len:
                  blockCount.append(1)
                  blockStarts.append(str(0))
                  blockSizes.append(probe_len)
                  stop.append(st[1])
                  start_end_stop = st[1]
                else:
                  blockStarts.append(str(0) + ',' + str(abs(e[n][-1]-e[n][0])) + ',')
                  blockCount.append(2)
                  blockSizes.append(str(block1)+','+str(probe_len - block1) + ',')
                  stop.append(st[0])
                  start_end_stop = st[0]
              else:
                blockCount.append(1)
                blockStarts.append(str(0))
                blockSizes.append(probe_len)
                stop.append(st[0])
                start_end_stop = st[0] 
              start_end.append(chromosome[0] + str(start_end_start)+':'+str(start_end_stop))
            
        for index in seq:
          if geneome_browser == 'IGV':
              # name columbn in gff3 style https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
              name_full = 'Name='+gene_sym+';'
              name_full += 'ID='+gID+';'
              name_full += 'Gene_specficity='+str(gene_specificity_index_local[n])+';'
              name_full += 'Isoform_specficity='+str(isoform_specificity_index_local[n])+';'
              name_full += 'Melt_temp(C)='+str(Tm[n])+';'
              name_full += 'Sequence='+seq[n]+';'
          elif geneome_browser  == 'UCSD' or geneome_browser == None:
            name_full  = gene_sym
          name.append(name_full)
          
        rep = len(s)
        if this_strand == '1' or this_strand !='-1':
            st = '+'
        elif this_strand == '-1':
            st = '-'
        else:
            st = '.'
        strand += [st] * rep
        thickStart += ['0'] * rep
        thickEnd += ['0'] * rep
        itemRgb += ['0'] * rep
    data = {'chr':chr ,'start':start,'stop':stop, 'gene_sym':name,'score':score_,
            'strand':strand, 'thickStart':thickStart,'thickEnd':thickEnd,'itemRgb':itemRgb,'blockCount':blockCount,
            'blockSizes':blockSizes, 'blockStarts': blockStarts,'start_end':start_end}
    cols = ['chr','start','stop', 'gene_sym','score','strand', 
            'thickStart','thickEnd','itemRgb','blockCount', 'blockSizes', 'blockStarts','start_end']    
    df = pd.DataFrame(data,columns=cols)
    df.drop_duplicates(subset=['start_end'],keep="first", inplace=True)
    df.drop('start_end', axis = 1,inplace=True)

    df_list = df.values.tolist()
    count = 0
    while True:
        path = os.path.normpath(os.path.join(file_path,probe_bed_file_name+'_'+str(count)+'.bed'))
        if os.path.exists(path) == True:
            count +=1
        else:
            break
    bed_path = os.path.normpath(os.path.join(file_path,probe_bed_file_name+'_'+str(count)+'.bed'))
    with open(bed_path, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        if geneome_browser  == 'IGV': 
          tsv_writer.writerow(['#gffTags'])
        for row in df_list:
            tsv_writer.writerow(row)
    print('probe bed file at '+ bed_path)

# =============================================================================

def write2csv(probe_sets, file_path,csv_name,overhang5=None,overhang3=None):
    best_probes = probe_sets
    if os.path.exists(os.path.normpath(best_probes)) == True:
        all_gene_ID = [f for f in os.listdir(probe_sets) if not f.startswith('.')]
        all_gene_ID = [f for f in all_gene_ID if not f.endswith('.txt')]
        if 'cross.pk' in all_gene_ID:
            all_gene_ID.remove('cross.pk')
    ids = all_gene_ID

    data = {'gene_ID':[0], 'gene_sym':[0], 'seq':[0], 'Tm_l':[0], 'Tm_r':[0], 'CGl':[0], 'CGr':[0],
       'isoform_specificity_index_l':[0], 'isoform_specificity_index_r':[0],
       'isoform_specificity_index_c':[0], 'gene_specificity_index_l':[0],
       'gene_specificity_index_r':[0], 'gene_specificity_index_c':[0], 'chr':[0], 'start':[0],
       'stop':[0], 'seqfull':[0], 'GCinRight':[0], 'd_ideal_t':[0]}
            
    cols = ['gene_ID', 'gene_sym', 'seq', 'Tm_l', 'Tm_r', 'CGl', 'CGr',
       'isoform_specificity_index_l', 'isoform_specificity_index_r',
       'isoform_specificity_index_c', 'gene_specificity_index_l',
       'gene_specificity_index_r', 'gene_specificity_index_c', 'chr', 'start',
       'stop', 'seqfull', 'GCinRight', 'd_ideal_t']
    df = pd.DataFrame(data,columns=cols)
    for gID in ids:
        df_ = load_pickle(os.path.normpath(os.path.join(probe_sets,gID)))
        #gID = gID.replace('.pk','')
        if df_.seq.tolist()[0] == 'no_probes_for_this_gene':
            continue
        df = df.append(df_)
    
    df.reset_index(drop=True,inplace=True)
    df = df.drop(df.index[0])
    seq_count = len(df['seq'].tolist()[0])
    if seq_count == 2:
        split = True
    else:
        split = False
    
    if overhang5 != None or overhang3 != None:
      if overhang5 == None:
        overhang5 = ''
      if overhang3 == None:
        overhang3 = ''
      if split == True:
          df['seq_l+OHs'] = df['seq'].apply(lambda x: overhang5 + x[0])
          df['seq_r+OHs'] = df['seq'].apply(lambda x: x[1] + overhang3)
      if split == False:
          df['seq_l+OHs'] = df['seq'].apply(lambda x: overhang5 + x)
          df['seq_r+OHs'] = df['seq'].apply(lambda x: x + overhang3)
      if split == True:
          df['seq'] = df['seq'].apply(lambda x: x[0]+' ,'+x[1])
    
    df.rename(columns = {'seq': "targeting_seqs"}, inplace=True)
    df.drop('d_ideal_t',axis=1,inplace = True )
    df.drop('GCinRight',axis=1,inplace = True )
    s = df['start'].tolist()
    s_=[]
    for start in s:
        if type(start) == float:
            continue
        
        if len(start) == 2:
            s_.append(str(start[0])+' ,'+str(start[1]))
        else:
            s_.append(str(start[0]))
    
    e = df['stop'].tolist()

    e_=[]
    for end in e:
        if type(end) == float:
            continue
        if len(end) == 2:
            e_.append(str(end[0])+' ,'+str(end[1]))
        else:
            e_.append(str(end[0]))
            
    
    df['start'] = s_
    df['stop'] = e_
    
    
    count = 0
    while True:
        path = os.path.normpath(os.path.join(file_path,csv_name+'_'+str(count)+'.csv'))
        if os.path.exists(path) == True:
            count +=1
        else:
            break
    path = os.path.normpath(os.path.join(file_path,csv_name+'_'+str(count)+'.csv'))
    df.to_csv(path, index = False)
    print('csv at '+ path)
    return(path)

# =============================================================================

def write2bed(probe_sets,
              Transcriptome,
              file_path,
              probe_bed_file_name,
              transcriptom_bed_file_name,
              geneome_browser  = None): # ('IGV' or 'UCSD')
    #first make bed file for the probes
    best_probes = probe_sets
    
    if os.path.exists(os.path.normpath(best_probes)) == True:
        all_gene_ID = [f for f in os.listdir(probe_sets) if not f.startswith('.')]
        all_gene_ID = [f for f in all_gene_ID if not f.endswith('.txt')]
        if 'cross.pk' in all_gene_ID:
            all_gene_ID.remove('cross.pk')

    ids = all_gene_ID    
    chr_ = []
    start = []
    stop = []
    name =[]
    gene_specificity_index = []
    strand = []
    thickStart = []
    thickEnd = []
    itemRgb = []
    blockCount = []
    blockSizes = []
    blockStarts = []
    start_end = []
    for gID in ids:
        df = load_pickle(os.path.normpath(os.path.join(best_probes,gID)))
        # check for empyth DF
        if df.empty == True: 
            continue
        if df.seq.tolist()[0] == 'no_probes_for_this_gene':
            continue
        gID = gID.replace('.pk','')
        
        gene_sym = df['gene_sym'].values.tolist()[0]

        probe_len = len(df.iloc[0]['seqfull'])
        chromosome = df['chr'].values.tolist()[0]
        
        gene_specificity_index += df['gene_specificity_index_c'].values.tolist()
        gene_specificity_index_local = df['gene_specificity_index_c'].values.tolist()
        isoform_specificity_index_local = df['isoform_specificity_index_c'].values.tolist()
        seq = df['seq'].values.tolist()
        Tml = df['Tm_l'].values.tolist()
        Tmr = df['Tm_r'].values.tolist()
        Tm = []
        for Tml,Tmr in zip(Tml,Tmr):
            Tm.append(str(round(Tml,2))+','+str(round(Tmr,2)))
        #Tm = df['Tm'].values.tolist()
        e = df['stop'].values.tolist()
        s = df['start'].values.tolist()
        if s == [] or s == None:
            continue
        if Transcriptome.TS.get(gID) == None:
            #this_chr = Transcriptome.G[gID].chr
            this_strand = Transcriptome.G[gID].strand
        else:
            #this_chr = Transcriptome.G[Transcriptome.TS[gID].gene_ID].chr
            this_strand = Transcriptome.G[Transcriptome.TS[gID].gene_ID].strand
        #chr += [this_chr] * rep
        
        if this_strand == '1' or this_strand !='-1':
          for n,st in enumerate(s):
              start.append(st[0])
              start_end_start = st[0]
              if len(st)>1:
                block1 = abs(st[0]- e[n][0])
                if block1 == probe_len:
                  blockCount.append(1)
                  blockStarts.append(str(0))
                  blockSizes.append(probe_len)
                  stop.append(e[n][0])
                  start_end_stop = e[n][0]
                else:
                  blockStarts.append(str(0) + ',' + str(abs(st[0]-st[1])) + ',')
                  blockCount.append(2)
                  blockSizes.append(str(block1)+','+str(probe_len - block1) + ',')
                  stop.append(e[n][-1])
                  start_end_stop = e[n][-1]
              else:
                blockCount.append(1)
                blockStarts.append(str(0))
                blockSizes.append(probe_len)
                stop.append(e[n][-1])
                start_end_stop = e[n][-1]
              start_end.append(chromosome + str(start_end_start)+':'+str(start_end_stop))

        elif this_strand == '-1':
          for n,st in enumerate(s):
              start.append(e[n][-1])
              start_end_start = e[n][-1]
              if len(st)>1:
                block1 = abs(e[n][-1]- st[1])
                if block1 == probe_len:
                  blockCount.append(1)
                  blockStarts.append(str(0))
                  blockSizes.append(probe_len)
                  stop.append(st[1])
                  start_end_stop = st[1]
                else:
                  blockStarts.append(str(0) + ',' + str(abs(e[n][-1]-e[n][0])) + ',')
                  blockCount.append(2)
                  blockSizes.append(str(block1)+','+str(probe_len - block1) + ',')
                  stop.append(st[0])
                  start_end_stop = st[0]
              else:
                blockCount.append(1)
                blockStarts.append(str(0))
                blockSizes.append(probe_len)
                stop.append(st[0])
                start_end_stop = st[0] 
              start_end.append(chromosome[0] + str(start_end_start)+':'+str(start_end_stop))
            
        for n, index in enumerate(seq):
          if geneome_browser == 'IGV':
              # name columbn in gff3 style https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
              name_full = 'Name='+gene_sym+';'
              name_full += 'ID='+gID+';'
              name_full += 'Gene_specficity='+str(gene_specificity_index_local[n])+';'
              name_full += 'Isoform_specficity='+str(isoform_specificity_index_local[n])+';'
              name_full += 'Melt_temp(C)='+str(Tm[n])+';'
              name_full += 'Sequence='+str(seq[n][0]+'---'+seq[n][1])+';'
          elif geneome_browser  == 'UCSD' or geneome_browser == None:
            name_full  = gene_sym
          name.append(name_full)
          
        rep = len(s)
        if this_strand == '1' or this_strand !='-1':
            st = '+'
        elif this_strand == '-1':
            st = '-'
        else:
            st = '.'
        strand += [st] * rep
        thickStart += ['0'] * rep
        thickEnd += ['0'] * rep
        itemRgb += ['0'] * rep
        chr = ['chr'+chromosome] * rep
        chr_ += chr
    #forP = [chr_,start,stop,name,gene_specificity_index,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts,start_end]
    #for item in forP:
      #print(len(item))    
    data = {'chr':chr_ ,'start':start,'stop':stop, 'gene_sym':name,'gene_specificity_index':gene_specificity_index,
            'strand':strand, 'thickStart':thickStart,'thickEnd':thickEnd,'itemRgb':itemRgb,'blockCount':blockCount,
            'blockSizes':blockSizes, 'blockStarts': blockStarts,'start_end':start_end}
    cols = ['chr','start','stop', 'gene_sym','gene_specificity_index','strand', 'thickStart','thickEnd','itemRgb','blockCount', 'blockSizes', 'blockStarts','start_end']    
    df = pd.DataFrame(data,columns=cols)
    df.drop_duplicates(subset=['start_end'],keep="first", inplace=True)
    df.drop('start_end', axis = 1,inplace=True)

    df_list = df.values.tolist()
    count = 0
    while True:
        path = os.path.normpath(os.path.join(file_path,probe_bed_file_name+'_'+str(count)+'.bed'))
        if os.path.exists(path) == True:
            count +=1
        else:
            break
    bed_path = os.path.normpath(os.path.join(file_path,probe_bed_file_name+'_'+str(count)+'.bed'))

    with open(bed_path, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        if geneome_browser  == 'IGV': 
          tsv_writer.writerow(['#gffTags'])
        for row in df_list:
            tsv_writer.writerow(row)

    #make a bed file to build a transcriptome for IGV
    if transcriptom_bed_file_name != None:
      chrom = []
      chromStart = []
      chromEnd = []
      name = []
      score = []
      strand = []
      thickStart = []
      thickEnd = []
      itemRgb = []
      blockCount = []
      blockSizes = []
      blockStarts = []
    
      for ts in Transcriptome.TS:
        T = Transcriptome
        gID = T.TS[ts].gene_ID
        chrom.append('chr'+str(T.G[gID].chr))
        score.append(str(T.TS[ts].exp))
        this_strand = T.G[gID].strand
        if this_strand == '1':
            st = '+'
            startL = T.TS[ts].start
            stopL = T.TS[ts].stop
            starts = ','.join(str(e-startL) for e in T.TS[ts].ex_exon_index_start)
            Tstart = str(T.TS[ts].ex_exon_index_start[0])
            Tend = str(T.TS[ts].ex_exon_index_stop[-1])
            
        elif this_strand == '-1':
            st = '-'
            startL = T.TS[ts].stop
            stopL = T.TS[ts].start
            L_ex_exon_index_stop = T.TS[ts].ex_exon_index_stop
            starts = [str(e-startL) for e in L_ex_exon_index_stop]
            starts.reverse()
            starts = ','.join(starts)
            Tstart = str(T.TS[ts].ex_exon_index_stop[-1])
            Tend = str(T.TS[ts].ex_exon_index_start[0])
            
        else:
            st = '.'
            startL = T.TS[ts].start
            stopL = T.TS[ts].stop
            starts = ','.join(str(e-startL) for e in T.TS[ts].ex_exon_index_start)
            Tstart = str(T.TS[ts].ex_exon_index_start[0])
            Tend = str(T.TS[ts].ex_exon_index_stop[-1])

        blockStarts.append(starts)        
        thickStart.append(Tstart)
        thickEnd.append(Tend)   
        chromStart.append(str(startL))
        chromEnd.append(str(stopL))
        strand.append(st)
        # name columbn in gff3 style https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        name_full = 'Name='+T.G[gID].symbol+';'
        name_full += 'ID='+ts+';'
        name_full += 'Gene_ID='+gID+';'
        name_full += 'Expression='+str(T.TS[ts].exp)+';'
        name_full += 'Biotype='+T.TS[ts].biotype+';'
        name.append(name_full)
        
        itemRgb.append('.')
        blockCount.append(str(len(T.TS[ts].in_exon_index_start)))
        
        sizes = []
        
        for n, start_exon in enumerate(T.TS[ts].ex_exon_index_start):
            stop_exon = T.TS[ts].ex_exon_index_stop[n]
            sizes.append(str(abs(start_exon-stop_exon)))
        sizes.reverse()
        sizes = ','.join(sizes)
        blockSizes.append(sizes)

      #  #gffTags line to the beginning of a .bed file,
      data = {'chr':chrom ,'start':chromStart,'stop':chromEnd, 'gene_sym':name,'gene_specificity_index':score,'strand':strand, 'thickStart':thickStart,'thickEnd':thickEnd,'itemRgb':itemRgb,'blockCount':blockCount, 'blockSizes':blockSizes, 'blockStarts': blockStarts}
      cols = ['chr','start','stop', 'gene_sym','gene_specificity_index','strand', 'thickStart','thickEnd','itemRgb','blockCount', 'blockSizes', 'blockStarts']    
    
      df = pd.DataFrame(data,columns=cols)
      df_list = df.values.tolist()
      count = 0
      while True:
        path = os.path.normpath(os.path.join(file_path,transcriptom_bed_file_name+'_'+str(count)+'.bed'))
        if os.path.exists(path) == True:
            count +=1
        else:
            break
      bed_path2 = os.path.normpath(os.path.join(file_path,transcriptom_bed_file_name+'_'+str(count)+'.bed'))
      with open(bed_path2, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['#gffTags'])
        for row in df_list:
            tsv_writer.writerow(row)
            
      print('transcriptom bed file at '+ bed_path2)
    
def get_exon_sizes(transcriptome):
    exon_sizes = []
    for ts in transcriptome.TS:
        if len(transcriptome.TS[ts].in_exon_index_start)<2:
            continue
        for start,stop in zip([transcriptome.TS[ts].in_exon_index_start[0]],[transcriptome.TS[ts].in_exon_index_stop[0]]):
            exon_sizes.append(stop-start)
    return(exon_sizes)

def get_intron_sizes(transcriptome):
    intron_sizes = []
    for ts in transcriptome.TS:
        if len(transcriptome.TS[ts].in_exon_index_start)<2:
            continue
        #fram shift list 
        for intron_start,intron_stop in zip(transcriptome.TS[ts].ex_exon_index_stop[:-1],transcriptome.TS[ts].ex_exon_index_start[1:]):
             intron_sizes.append(abs(intron_stop - intron_start))
    return(intron_sizes)

def get_transcript_sizes(transcriptome):
    Ts_size = []
    for ts in transcriptome.TS:
        Ts_size.append(len(transcriptome.TS[ts].seq))
    return(Ts_size)
        
def get_biotypes(transcriptome):
    bt = {}
    for ts in transcriptome.TS:
        ts_bt = transcriptome.TS[ts].biotype
        if bt.get(ts_bt) == None:
            bt[ts_bt] = []
        bt[ts_bt].append(ts)
    return(bt)

def get_junctions(TS_object, len_on_each_side):
    junctions = []
    if len(TS_object.in_exon_index_start)<2:
        return(junctions)
    for junc in TS_object.in_exon_index_start[1:]:
        junctions.append(TS_object.seq[junc-len_on_each_side:junc+len_on_each_side])
    return(junctions)
    
def get_splice_sites(transcriptome, len_on_each_side = 10):
    
    splice_sites = dict()
    
    for ts in transcriptome.TS:
        G = transcriptome.TS[ts].gene_ID
        junctions = get_junctions(transcriptome.TS[ts],len_on_each_side)
        for jun in junctions:
            jun = DNA2int(jun, N = 'A')
            if splice_sites.get(jun) != None: 
                splice_sites[jun].append(G)
            else: 
                splice_sites[jun] = [G]
    return(splice_sites)

def get_gene_size(transcriptome):
    size = []
    for g in transcriptome.G:
        size.append(len(transcriptome.G[g].seq)/1000)
    return(size)

def get_num_exons(transcriptome):
    introns = []
    for ts in transcriptome.TS:
        introns.append(len(transcriptome.TS[ts].in_exon_index_start))
    return(introns)      


# first we need to build the Transcriptome object 

# you will need 3 files to start 
#1. cdna_fasta
#2. geneome_fasta
#3. Human_exons (csv)

#cdna_fasta: from ensembl website - Homo_sapiens.GRCh38.cdna.all.fa
#http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

#Human_geneome.txt
#geneome_fasta: from ensembl biomart in FASTA header form (unspliced gene) sequence (sep = |):  
  #geneID(no version number)
  #gene name
  #discription
  #chr
  #start(bp)
  #stop(bp)
  #gene type

#Human_exons.txt
#exon_CSV: from ensembl biomart in CSV form (sep = ,): 
 #Gene stable ID (no version number)
 #Transcript stable ID (no version number)
 #Strand
 #Transcript start (bp)
 #Transcript end (bp)
 #Exon region start (bp)
 #Exon region end (bp)

#cdna_fasta example: >ENSG00000000419|DPM1|dolichyl-phosphate mannosyltransferase subunit 1, catalytic [Source:HGNC Symbol;Acc:HGNC:3005]|20|50934867|50958555|protein_coding, geneID, symbol, gene name, chromosome, start, stop, type
#exon_CSV exampl: Gene stable ID,Transcript stable ID,Strand,Transcript start (bp),Transcript end (bp),Exon region start (bp),Exon region end (bp)

out_put = os.path.normpath('/Human_geneome/')

geneome_fasta = os.path.normpath('/Human_geneome//Human_geneome.txt')
cdna_fasta = os.path.normpath('/Human_geneome/Homo_sapiens.GRCh38.cdna.all.fa')
exon_CSV = os.path.normpath('/Human_geneome/Human_exons.txt')

# read exon file 
exon = pd.read_csv(exon_CSV)
# get IDs
ts_uni = list(set(exon['Transcript stable ID'].tolist()))

# set expression to 1 for all genes 
Human_exp = {}
for ts in ts_uni:
    Human_exp[ts] = 1

# build the Transcriptome object 
HumanALL = Build_Transcriptome(
    cell_expr_iso_dirt = Human_exp,
    cdna_fasta = cdna_fasta,
    geneome_fasta = geneome_fasta,
    exon_CSV = exon_CSV,
    cell_type = 'Human_all',
    expression_cutoff=0.00000001)

save_pickle(HumanALL, out_put, 'HumanALL_transcriptome')

HumanALL = load_pickle(os.path.normpath('/Human_geneome//HumanALL_transcriptome.pk'))

# calulate all the thermodynamic propoties of seqs
HumanALL.get_TDP(file_path = out_put, replace = False, verbose = True) 

# calulate times a subsequence is in trascriptome, this will take a while 
Human_OT_dict_Full = Make_OT_dict_Full(HumanALL,
                      file_path = out_put,
                      OT_dict_file_name = 'OT_dict',
                      OT_length = 15,
                      verbose = False,
                      frag = 20, # number of fragments to split up look up table
                      threads = 4)

# after this is calulated we can calualte tables for our specfic genes 
two_genes_TSome = HumanALL.subset(subset_list = ['ENST00000374980','ENST00000372134'])

# now we calulate specifisity metrics for these 
ot_iso_PBMC =  get_OT_dict_Iso(transcriptome = two_genes_TSome,
                                 OT_dict_dict = r'/Human_geneome/db/OT_dict_0/',
                                 file_path = r'/two_genes_probes/',
                                 verbose = True,
                                 frag = 20, # number of fragments from above
                                 threads = 4)


# set path in trascriptome object 
ot_iso = r'/two_genes_probes/db/OT_dict_iso0'
two_genes_TSome.OT_db_path = ot_iso 

# make probes for genes in two_genes
probes = get_split_probes(transcriptome = two_genes_TSome ,
               probe_length =60,
               probe_db_out_path = r'/two_genes_probes/',
               file_name = 'probe_set_',
               overhang_5 = '',
               overhang_3 = '',
               ban_seqs_dict = None,
               iso_form = True,
               verbose = True, 
               exon_only = True,
               ban_list = ['AAAA','CCCC','GGGG','TTTT'],
               threads = 16)

# pick best probes
pick_best_probes(probes,
                    output_path=r'/two_genes_probes/',
                    tm_low = 0,
                    tm_high = 200,
                    formamide_fraction = 0.3,
                    monovalentSalt = 0.79,
                    probeConc = 0.000000001,
                    medadata_out_path = r'/two_genes_probes/',
                    output_file_name = 'two_genes_probes', 
                    CG_low=0.3,
                    CG_high=0.7,
                    iso_low=0,
                    iso_high=1.1,
                    min_probe_spaceing = 0,
                    max_probes= 2, 
                    cross_hyb = 30,
                    splice = None, 
                    gene_ID = 'all',
                    Chem_corrected_Tm = False,
                    iso = False,
                    verbose = False,
                    ligase = 'splintR', # not splintR # dont inforce the probe to exclue C,G at ligation point
                    cross = dict(),
                    exon_intron_only = True,
                    threads = 16)



write2csv(r'/two_genes_probes/two_genes_probes', 
         r'/two_genes_probes/two_genes_probes',
          r'two_genes_probes',
          overhang5=None,overhang3=None)

