#! /usr/bin/python3
import os
import sys
import csv
import json
import time
import logging
import argparse
import pickle
from collections import namedtuple
from collections import defaultdict
from collections import Counter
from difflib import SequenceMatcher
import numpy as np
from sklearn.metrics import jaccard_similarity_score
class Index:

    def __init__(self):
        return

    def fasta_parser(self, fasta):
        '''Parse fasta and return iterator with header and sequence'''
        fasta_handle = open(fasta)
        fasta_content = fasta_handle.read().split('>')[1:]
        reference = list()
        transcriptid = 0
        Fasta = namedtuple('Fasta',['header','seq','tid','length'])
        for lines in fasta_content:
            header = lines.split('\n')[0]
            sequence = ''.join(lines.split('\n')[1:]).upper()
            length = len(sequence)
            record = Fasta(header, sequence, transcriptid, length)
            transcriptid += 1
            yield record
    
    def split_fasta(self, fasta):
        '''Split fasta into chunks and return list of chunks'''
        record = Index.fasta_parser(self, fasta)
        pitcher = list()
        mugs = list()
        count = 0
        while True:
            try:
                line = next(record)
                mugs.append(line)
                count += 1 
                if count >0 and count % 10000 == 0:
                    pitcher.append(mugs)
                    mugs = list()
            except StopIteration:
                pitcher.append(mugs)
                break
        return(pitcher)

    def kmerize(self, sequence, klen):
        '''K-mer creator'''
        mask = (1 << (klen*2)) -1 
        kmer = 0
        bit_counter = 1
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
                frag = kmer & mask
                yield frag 
            else:
                bit_counter += 1
    
    def kmetric(self, frag):
        '''Metrics for a k-mer'''
        gc = frag.count('G') + frag.count('C')
        at = len(frag) - gc
        return (at) 
    
    def indexer(self, fasta, klen):
        '''Create index dictionary'''
        kmer_transcript = dict()
        transcript_index = list()
        lineno = 0
        for sequence in Index.fasta_parser(self, fasta):
            transcript_index.append([sequence.header,sequence.length])
            for kmers in Index.kmerize(self, sequence.seq, klen):
                if kmers == 0:
                    continue
                try:
                    kmer_transcript[kmers].append(sequence.tid)
                except KeyError:
                    kmer_transcript[kmers] = [sequence.tid]
        return(kmer_transcript, transcript_index)    


    def index_writer(self, kmers, genes):
        '''Writer module'''
        klist = open('klist.kx','w')
        glist = open('glist.gx','w')
        for lines in sorted(kmers):
            index = ','.join('{0}'.format(gene) for gene in kmers[lines])
            klist.write('{0}\t{1}\n'.format(lines, index))
        for lines in genes:
            glist.write('{0}\t{1}\n'.format(lines[0], lines[1]))
        
        return

class Quant:
    
    def __init__(self):
        return
    
    def readkc(self, kcfile):
        '''kc file reader module; modify to read ikc files'''
        kchandle = open(kcfile)
        kcount = namedtuple('kcount',['kmer', 'count'])
        for lines in kchandle:
            fields = lines.strip().split('\t')
            kcinfo =  kcount(int(fields[0]), int(fields[1]))
            yield kcinfo

    def readkx(self, kxfile):
        '''kx file reader module; modify to read ikx files'''
        kxhandle = open(kxfile)
        kindex = namedtuple('kindex', ['kmer', 'gene'])
        for lines in kxhandle:
            fields = lines.strip().split('\t')
            kmer = int(fields[0])
            gene = [int(val) for val in fields[1].split(',')]
            index = kindex(kmer, gene)
            yield index

    def readgx(self, gxfile):
        '''gx file reader module; will be removed when ikx files'''
        gxhandle = open(gxfile)
        gindex = namedtuple('gindex', ['gene','length'])
        index = list()
        for lines in gxhandle:
            ginfo = lines.strip().split('\t')
            index.append(gindex(lines[0],lines[1]))
        return index
    
    def writekc(self, kcinfo):
        '''write kc file'''
        output = open('unmatched.kc', 'w')
        for lines in kcinfo:
            output.write('{0}\t{1}\n'.format(lines.kmer, lines.count))
        return

    def expression(self, ksize, kcount, glen, gcount):
        elen = glen - ksize + 1
        kpkm = (gcount * (10**9) )/ float(elen * kcount)
        return kpkm

    def writekx(self, gcdict, gindex):
        output = open('expression.kx', 'w')
        output.write('TranscriptID\tLength\tk-mer_count\tkpkm\n')
        for lines in gcdict:
            gid = gindex[lines][0]
            glen = gindex[lines][1]
            gcount = gcdict[lines]
            kpkm = expression(ksize, kcount, tlen, gcount)
            output.write('{0}\t{1}\t{2}\t{3}\n'.format(gid, glen, gcount, kpkm))
        return

    #add maximum likelihood module
    #add ikc reader and writer

    def genetocount(self, kcinfo, kindex):
        '''merging counts and transcript list on kmers'''
        gcdict = dict()
        gene = next(kindex)
        count = next(kcinfo)
        kcount = 0
        gdata = open('gdata.ind','w')
        while True:
            try:
                if gene.kmer < count.kmer:
                    gene = next(kindex)
                elif gene.kmer > count.kmer:
                    count = next(kcinfo)
                else:
                    genes = gene.gene
                    glen = len(genes)
                    gdata.write('{0}\t{1}\n'.format(gene.gene, count.count))
                    for names in genes:
                        try:
                            gcdict[names] += count.count/float(glen)
                        except KeyError:
                            gcdict[names] = count.count/float(glen)
                    kcount += count.count
                    gene = next(kindex)
                    count = next(kcinfo)
            except StopIteration:
                break  

        while gene != None:
            try:
                genes  = gene.gene
                for names in genes:
                    try:
                        gcdict[names] += 0
                    except KeyError:
                        gcdict[names] = 0
                gene = next(kindex)
            except StopIteration:
                break
        
        while count != None:
            try:
                Quant.writekc(count)
                count = next(kcinfo)
            except StopIteration:
                break
        return(gcdict)
                    
if __name__ == '__main__':
    troch =  argparse.ArgumentParser(prog='Hummingbird')
    troch.add_argument('-r', '--reference', type=str, dest='fasta',
        help='Reference fasta file path')
    troch.add_argument('-l', '--klim', type=int, dest='klim', default=None,
        help='k-mer limit')
    troch.add_argument('-k', '--klen', type=int, dest='klen', default=27,
        help='k-mer length')
    troch.add_argument('-m', '--mem', type=int, dest='memoc', default=10000,
        help='Number of sequences to be loaded in memory')
    troch.add_argument('-v', '--version', action='version', version='%(prog)s 0.9.8')
    troch.add_argument('-c', '--kcfile', type=str, help='kc file path')
    troch.add_argument('-g', '--gxfile', type=str, help='gx file path')
    opts = troch.parse_args()
    record = Index()
    quant = Quant()
    #kmer, gene = record.indexer(opts.fasta, opts.klen)
    #record.index_writer(kmer, gene)
    #pitcher = record.split_fasta(opts.fasta)
    #for lines in pitcher:
    #    print(len(liines))
    kc = quant.readkc(opts.kcfile)
    print('Kc generator created')
    gx = quant.readkx(opts.gxfile)
    print('Kx generator created')
    gcount = quant.genetocount(kc, gx)
    print('Quant completed')
