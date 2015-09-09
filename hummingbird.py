#! /project/home/sravishankar9/local/bin/python3.4
""" KXpress is a k-mer based rna expression quantification tool developed by
Shashidhar Ravishankar, at the Vannberg Lab, Georgia Institute of Technology.
    There are two modules to KXpress, namely the index and the express module.
The index module, takes as input the reference cDNA fasta file or probes list,
and creates a hash based index, for quick k-mer look up. Indexing needs to be
done once for a given k-mer size.
    The express module, takes as input the index file, and a sequence file in
fastq or fastqgz format and using KAnalyze, performs a quick k-mer count and 
calculates gene expression in term of k-mers per kilobase per million k-mers 
mapped."""
import os
import sys
import csv
import time
import math
import glob
import gzip
import shutil
import pickle
import profile
import argparse
import subprocess
import statistics as st
import logging
import tempfile
from operator import itemgetter
from Bio import SeqIO
import multiprocessing
from collections import Counter
from itertools import repeat
from natsort import natsorted
import numpy as np
import collections
from scipy import stats
global kanalyze
kanalyze = os.path.abspath(os.path.dirname(__file__))+'/kanalyze-0.9.7/'

def outlier_detection(transcript_index,kmer_count,transcript_order,names):
    logging.info('Detecting outliers')
    zscore_file = csv.writer(open('zscore_uniq.csv','w'),delimiter='\t')
    dist_file = csv.writer(open(names+'_dist_corrected.csv','w'),delimiter='\t')
    for genes in transcript_index:
        genemed = np.median(transcript_index[genes])
        mad = np.median([abs(vals - genemed) for vals in transcript_index[genes]])
        val_string = ','.join([str(vals) for vals in transcript_index[genes]])
        dist_file.writerow([transcript_order[genes],str(mad),val_string])
        zscores = [0.6745*(vals - genemed)/mad for vals in transcript_index[genes]]
        val_string = ','.join([str(vals) for vals in zscores])
        for i in range(len(zscores)):
            if zscores[i] >= 3.5:
                transcript_index[genes][i] = genemed
        val_string = ','.join([str(vals) for vals in transcript_index[genes]])
    return(transcript_index,kmer_count)

def call_expression(kc_chunk,index,transcript_index,kmer_med,kmer_mad):
    logging.info('Calling expression')
    kmer_count = 0
    for lines in open(index):
        line = lines.strip().split('\t')
        if int(line[0]) in kc_chunk: # and int(line[1]) == 1:
            kmer_count += kc_chunk[int(line[0])]
            for genes in line[2].split(','):
                gene, pos = genes.split(':')
                transcript_index[int(gene)][int(pos)-1] = kc_chunk[int(line[0])]    #/int(line[1])
    return (kmer_count)

def express(kc_chunk,index,transcript_order,transcript_length,klen,seq_len, kmer_count):
    for lines in open(index):
        line = lines.strip().split('\t')
        if int(line[1]) <= 1:
            for genes in line[2].split(','):
                try:
                    if int(line[0]) in kc_chunk:
                        transcript_index[int(genes)] += kc_chunk[int(line[0])]
                        kmer_set.add(int(line[0]))
                    else:
                        continue
                except KeyError:
                    if int(line[0]) in kc_chunk:
                        transcript_index[int(genes)] = kc_chunk[int(line[0])]
                        kmer_set.add(int(line[0]))
                    else:
                        continue
        else:
            if int(line[0]) in kc_chunk:
                res_trans = [int(gene) for gene in line[2].split(',')]
                transcript_rescue[int(line[0])] = [res_trans,kc_chunk[int(line[0])]]
            else:
                continue
    return (kmer_count)

def get_next(iterator):
    try:
        return(next(iterator).strip().split('\t'))
    except StopIteration:
        return(None)

def rescue(transcript_index,transcript_rescue):
     logging.info('Rescue algorithm started')
     sigma_kpkm = sum(transcript_kpkm.values())
     rescue = 0
     logging.info ('Number of transcripts quatified:' +str(len(transcript_index.keys())))
     for count,kmers in enumerate(transcript_rescue):	#loop through all kmers that have more than one transcript index
         if len(transcript_rescue[kmers][0]) <= 10:
             for indices in transcript_rescue[kmers][0]:	#loop through transcript index of kmer
                 try:
                     opt = transcript_rescue[kmers][1] * transcript_kpkm[indices]/float(sigma_kpkm)
                     transcript_index[indices] += opt
                     rescued += 1
                 except KeyError:
                     continue
     logging.info('Number of transcripts rescued :' +str(len(transcript_index.keys())))
     logging.info('Rescue completed')
     return (count)

def merge(file1,file2,database,index,output,klen):
     result = csv.writer(open(output+str(index)+'.mkx','w'),delimiter='\t')
     with open(file1) as f1, open(file2) as f2:
         line1 = get_next(f1) 
         line2 = get_next(f2) 
         while line1 != None and line2 != None:
             if int(line1[0]) < int (line2[0]):
                 result.writerow([line1[0],line1[1],line1[2]])
                 line1 = get_next(f1) 
             elif int(line1[0]) == int(line2[0]):
                 merge_list = line1[2].split(',') + line2[2].split(',')
                 result.writerow([line1[0],str(len(merge_list)),','.join(merge_list)])
                 line1 = get_next(f1) 
                 line2 = get_next(f2) 
             else:
                 result.writerow([line2[0], line2[1], line2[2]])
                 line2 = get_next(f2) 
         while line1 != None:
             result.writerow([line1[0],line1[1],line1[2]])
             line1 = get_next(f1) 
         while line2 != None:
             result.writerow([line2[0],line2[1],line2[2]])
             line2 = get_next(f2) 
     return(output+str(index)+'.mkx')

def lazy_function(fasta_file):
    """Function performs chunking of fasta file"""
    seq = True
    count = 0
    chunk = list()
    part = list()
    while fasta_file:
        if count < 5000:
            try:
                seq = next(fasta_file)
                part.append([seq.id,str(seq.seq).upper()])
                count += 1
            except StopIteration:
                chunk.append(part)
                break
        else:
            count = 0 
            chunk.append(part)
            part = list()
    return(chunk)

def kmerize(arguments):
    """Implements a quick k-mer algorithm, to k-merize a given sequence and
    store them as base of four integers. It also maps the k-mers generated
    to the respective transcripts."""
    fasta_file = arguments[0]
    index = arguments[1]
    klen = arguments[2]
    database = arguments[3]
    power = arguments[4]
    output = arguments[5]
    order = index * (10**power)
    transcript_length = dict()
    transcript_kmers = dict()
    transcript_order = dict()
    mask = (1 << (klen*2)) -1         #Create bitmask of 1's
    kmers = 0
    for count, lines in enumerate(fasta_file):
        tindex = count + order
        sequence = lines[1]
        header = lines[0]
        kmer = 0
        bit_counter = 1
        order_counter = 0
        for nuc in sequence:
            kmer =  kmer << 2         #left shift k-kmer 
            if nuc == 'A':            #add the value of the character using bitwise OR
                kmer = kmer | 0x00   
            elif nuc == 'C':
                kmer = kmer | 0x01
            elif nuc == 'G':
                kmer = kmer | 0x02
            elif nuc == 'T':
                kmer = kmer | 0x03
            else:
                bit_counter = 0
            if bit_counter == klen:   #if length equals k-mer length, store k-mer
                kmers += 1
                order_counter += 1
                try:
                    transcript_kmers[kmer & mask].append([tindex,order_counter])  #k-mer to transcript index mapping
                except KeyError:
                    transcript_kmers[kmer & mask] = [[tindex,order_counter]]
            else:
                bit_counter += 1
        logging.debug('Indexing '+header+ '; Line number '+str(tindex))
        transcript_order[tindex] = header         #transcript index to transcript name mapping
        transcript_length[header] = len(sequence) #transcript to transcript length mapping
    temp_file = csv.writer(open(output+str(index)+'.tkx','w'),delimiter='\t')
    for values in natsorted(transcript_kmers.keys()):
        temp_file.writerow([str(values),str(len(transcript_kmers[values])),str(','.join(':'.join([str(val) for val in x]) for x in transcript_kmers[values]))])
    return (transcript_order, transcript_length,kmers)

def probe_list(database, line_count, klen):
    """Method to create index from a given probe list"""
    #must be updated 
    logging.basicConfig(level=logging.INFO)
    dbfile = csv.reader(open(database),delimiter='\t')
    transcript_index = dict()
    transcript_length = dict()
    kmers = 0
    perc = [float(i) for i in range(10,110,10)]
    logging.info('Indexing reference probes file')
    for count, lines in enumerate(dbfile):
        if len(lines) > 20:
            kindex = 0
            for i in range(len(lines[17])):
                if len(lines[17][i:i+klen]) == klen:
                    if dec(lines[17][i:i+klen],klen) not in transcript_index:
                        transcript_index[dec(lines[17][i:i+klen],klen)] = [[lines[2],kindex,0]]
                    else:
                        transcript_index[dec(lines[17][i:i+klen])].append([lines[2],kindex,0])
                    kindex += 1
                    kmers += 1
            transcript_length[lines[2]] = [len(lines[17]),kindex]
    kmer_str = str(kmers)
    logging.info('Reference indexing completed; %s kmers indexed to transcripts' %kmers)
    return (transcript_index, transcript_length)

#Create transcript index, maintains refseq id, kmer and kmer index
def transcript_list(database, line_count,klen,output,threads,gene):
    """Method to create index from a fasta file"""
    logging.basicConfig(level=logging.INFO)
    transcript_length = dict()
    transcript_index = dict()
    transcript_kmers = dict()
    transcript_order = dict()
    logging.info('Indexing reference fasta file')
    fasta_file = SeqIO.parse(open(database),'fasta')
    fasta_list = lazy_function(fasta_file)
    power = 5
    count = 0
    pool = multiprocessing.Pool(processes=int(threads))
    results = pool.map(kmerize,zip(fasta_list,range(len(fasta_list)),repeat(klen),
                       repeat(database),repeat(power),repeat(output)))
    for i in range(len(results)):            #redundancy to be rectified
        for k, v in results[i][1].items():
            transcript_length[k] = v
        for k, v in results[i][0].items():
            transcript_order[k] = v
            transcript_index[k] = [0 for i in range(transcript_length[transcript_order[k]] - int(klen) +1)]
        count += results[i][2]
    logging.info('Number of transcripts read : %s ' %len(transcript_length))
    logging.info('Reference indexing completed; \
                  %s kmers indexed to transcripts' %count)
    temp_list = glob.glob(output+'*.tkx')
    i = 0
    logging.info('Merging index files')
    while len(temp_list) > 1:
        merged_file = merge(temp_list[0],temp_list[1],database, i,output,klen)
        file = temp_list.pop(0)
        os.remove(file)
        file = temp_list.pop(0)
        os.remove(file)
        temp_list.insert(0,merged_file)
        i += i
    basename = os.path.basename(output)
    shutil.move(glob.glob(output+'*.*kx')[0],output+'.mkx')
    logging.info('Merging complete')
    return (transcript_length, transcript_order, transcript_index)

def reporter(transcript_index,transcript_length,klen,output_file,seq_number,seq_len,transcript_order):
    logging.info('Generating report')
    output_csv = csv.writer(open(output_file,'w'),delimiter='\t')
    output_csv.writerow(['RefseqID','Length','KPKM','RPKM'])
    for refseqs in sorted(transcript_index):
        if refseqs != '' and transcript_length[transcript_order[refseqs]] > int(klen):
            kpkm = (sum(transcript_index[refseqs])*(10**9))/ float((transcript_length[transcript_order[refseqs]]-int(klen)+1) * seq_number) # * seq_len)
            rpkm = kpkm / float(seq_len - int(klen) +1)
            output_csv.writerow([transcript_order[refseqs], str(transcript_length[transcript_order[refseqs]]),str(kpkm),str(rpkm)])
        else:
            continue
    logging.info('Report generation complete')
    return

def kxpress(refer, seqfile, loglevel,klen,output,mode,index,db_type,type,seq_len,threads):
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
         raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level, format='%(levelname)s:%(asctime)s:%(message)s', datefmt='%m/%d/%Y;%I:%M:%S')
    db_count = 0
    seq_file_count = 0
    if mode == 'index':
        logging.info('Indexing started')
        gene = list()
        if db_type == 'fasta':
            with open(refer) as dbfile:
                for line in SeqIO.parse(dbfile,"fasta"):
                    db_count += 1
                    gene.append(line.id)
            logging.info(str(db_count)+' transcripts found')
            flanking_kmers = transcript_list(refer, db_count, int(klen),output,threads,gene)
            pickle.dump(flanking_kmers, open(output+'.kxp','wb'))
            logging.info('Indexing complete.')
        elif db_type == 'probes':
            flanking_kmers = probe_list(refer, db_count, int(klen))
            pickle.dump(flanking_kmers, open(output+'.kxp','wb'))
            logging.info('Indexing complete.')
    elif mode == 'express':
        logging.info('Parsing input files')
        logging.info('Loading index')
        transcript_length, transcript_order, transcript_index = pickle.load(open(index+'.kxp','rb'))
        basename = os.path.splitext(os.path.basename(seqfile[0]))[0].split('_')[0]
        if type == 'fastq' or type == 'fastqgz':
            if type == 'fastq':
                fastq_parser = open(seqfile[0]) 
                seq_id = fastq_parser.readline()
                seq_len = len(fastq_parser.readline()) -1
                fastq_parser.close()
            else:
                fastq_parser = gzip.open(seqfile[0],mode='rt')
                seq_id = fastq_parser.readline()
                seq_len = len(fastq_parser.readline()) -1 
                fastq_parser.close()
            logging.info('Kmerizing seqeunce files')
            kmer_file = os.path.abspath(glob.glob(index+'*.mkx')[0])
            kmerizes = subprocess.Popen([kanalyze+'count','-rcanonical','-m','dec','-d',threads,'-l',threads,'--countfilter=kmercount:c>1','-k',klen,'-f',type,'-o',os.path.abspath(os.path.dirname(__file__))+'/'+basename+'.kc',' '.join(seqfile)], stdout=subprocess.PIPE, shell=False)
            msg = kmerizes.communicate()
            logging.info('Kmers counted')
            kmer_file = os.path.abspath(os.path.dirname(__file__))+'/'+basename+'.kc'
        elif type == 'kc':
            kmer_file = seqfile
        logging.info('Retriving index')
        kmer_index = glob.glob(index+'*.mkx')[0]
        logging.info('Merging transcritps and k-mers')
        kmer_count = 0
        kmer_dist = list()
        for lines in open(kmer_file):
            kmer_count += int(lines.split('\t')[1])
            kmer_dist.append(int(lines.split('\t')[1]))
        kmer_med = np.median(kmer_dist)
        kmer_mad = np.median([(vals - kmer_med) for vals in kmer_dist])
        kmer_dist = list()
        logging.info('Total kmers found : '+str(kmer_count))
        kmer_csv = open(kmer_file)
        kmer = next(kmer_csv).split('\t')
        kmer_count = 0
        kc_chunk = dict()
        while kmer_csv:
            if len(kc_chunk) <= 50000000:
                kc_chunk[int(kmer[0])] = int(kmer[1])
                try:
                    kmer = next(kmer_csv).split('\t')
                except StopIteration:
                    kmer_count += call_expression(kc_chunk,kmer_index,transcript_index,kmer_med,kmer_mad)
                    break
            else:
                kc_chunk[int(kmer[0])] = int(kmer[1])
                kmer_count += call_expression(kc_chunk,kmer_index,transcript_index,kmer_med,kmer_mad)
                kc_chunk = dict()
        logging.info('Merging complete; Uniquely mapped k-mers:'+str(kmer_count))
        logging.info('Calculating initial kpkm')
        full_transcripts = 0
        incomplete_transcripts = 0
        unknown_kmers = 0
        for transcripts in transcript_index:
            if transcripts != '':
                full_transcripts += 1
            elif transcripts == '':
                unknown_kmers = len(transcript_index[transcripts])
        reporter(transcript_index,transcript_length,klen,output,kmer_count,seq_len,transcript_order)
        logging.info('Number of complete transcripts found = %s' %full_transcripts)
        logging.info('Number of kmers with no refseq id = %s' %unknown_kmers)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Hummingbird', 
            description='''Hummingbird is a k-mer based rna expression quantification tool developed by
Shashidhar Ravishankar, at the Vannberg Lab, Georgia Institute of Technology.
    There are two modules to Hummingbird, namely the index and the express module.
The index module, takes as input the reference cDNA fasta file or probes list,
and creates a hash based index, for quick k-mer look up. Indexing needs to be
done once for a given k-mer size.
    The express module, takes as input the index file, and a sequence file in
fastq or fastqgz format and using KAnalyze, performs a quick k-mer count and 
calculates gene expression in term of k-mers per kilobase per million k-mers 
mapped.''')
    parser.add_argument('-r','--reference',dest='refer', type=str, help='Transcript reference file')
    parser.add_argument('-f','--seqfile',dest='seqfile', type=str,nargs='+', help='Fastq file')
    parser.add_argument('-k','--kmer_length',dest='klen', type=str, help='Kmer length',default='20')
    parser.add_argument('-o','--output_file',dest='output',type=str, help='Output file path',default='transcripts.kx')
    parser.add_argument('-m','--mode',dest='mode',type=str,choices=['index','express'],help='Enter mode of execution, run index first if the database hasnt been indexed yet')
    parser.add_argument('-d','--db_type',dest='db_type',type=str,choices=['probes','fasta'],help='Database type being used for indexing')
    parser.add_argument('-t','--type', dest='file_type',type=str,choices=['fastq','kc','fastqgz'],help='Input file type') 
    parser.add_argument('-i','--index',dest='index',type=str,help='Path to indexed file')
    parser.add_argument('-l','--log',type=str,default='info',choices=['info','debug','error'],help='Verbosity parameter')
    parser.add_argument('-s','--seqlen',type=int,default=36,help='Read length, must be provided when runnnig using kc file')
    parser.add_argument('-v','--version',action='version',version='%(prog)s 0.9.7')
    parser.add_argument('-n','--num_thread',dest='threads',type=str,default='2',help='Number of threads')
    args = parser.parse_args()
    kxpress(args.refer, args.seqfile,args.log,args.klen,args.output,args.mode,args.index,args.db_type,args.file_type,args.seqlen,args.threads)
