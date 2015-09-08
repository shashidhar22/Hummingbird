Hummingbird (version 0.9.7)
===========================
>Hummingbird is a k-mer based rna expression quantification tool developed by Shashidhar Ravishankar, at the Vannberg Lab, Georgia Institute of Technology. There are two modules to Hummingbird, namely the index and the express module. The index module, takes as input the reference cDNA fasta file or probes list, and creates a hash based index, for quick k-mer look up. Indexing needs to be done once for a given k-mer size. The express module, takes as input the index file, and a sequence file in fastq or fastqgz format and using [KAnalyze](http://bioinformatics.oxfordjournals.org/content/30/14/2070) ,performs a quick k-mer count and calculates gene expression in term of k-mers per kilobase per million k-mers mapped.

Dependencies:
-------------
* [Python3.4.3](https://www.python.org/downloads/)
* [KAnalyze 0.9.7](http://sourceforge.net/projects/kanalyze/)
* [Biopython 1.65](http://biopython.org/wiki/Download)
* [Scipy 0.15](http://www.scipy.org/install.html)
* [natsort 4.0.3](https://pypi.python.org/pypi/natsort)


Example:
========
Runtime Options:
----------------
```python
usage: Hummingbird [-h] [-r REFER] [-f SEQFILE [SEQFILE ...]] [-k KLEN]
                   [-o OUTPUT] [-m {index,express}] [-d {probes,fasta}]
                   [-t {fastq,kc,fastqgz}] [-i INDEX] [-l {info,debug,error}]
                   [-s SEQLEN] [-v] [-n THREADS]

Hummingbird is a k-mer based rna expression quantification tool developed by
Shashidhar Ravishankar, at the Vannberg Lab, Georgia Institute of Technology.
There are two modules to Hummingbird, namely the index and the express module.
The index module, takes as input the reference cDNA fasta file or probes list,
and creates a hash based index, for quick k-mer look up. Indexing needs to be
done once for a given k-mer size. The express module, takes as input the index
file, and a sequence file in fastq or fastqgz format and using KAnalyze,
performs a quick k-mer count and calculates gene expression in term of k-mers
per kilobase per million k-mers mapped.

optional arguments:
  -h, --help            show this help message and exit
  -r REFER, --reference REFER
                        Transcript reference file
  -f SEQFILE [SEQFILE ...], --seqfile SEQFILE [SEQFILE ...]
                        Fastq file
  -k KLEN, --kmer_length KLEN
                        Kmer length
  -o OUTPUT, --output_file OUTPUT
                        Output file path
  -m {index,express}, --mode {index,express}
                        Enter mode of execution, run index first if the
                        database hasnt been indexed yet
  -d {probes,fasta}, --db_type {probes,fasta}
                        Database type being used for indexing
  -t {fastq,kc,fastqgz}, --type {fastq,kc,fastqgz}
                        Input file type
  -i INDEX, --index INDEX
                        Path to indexed file
  -l {info,debug,error}, --log {info,debug,error}
                        Verbosity parameter
  -s SEQLEN, --seqlen SEQLEN
                        Read length, must be provided when runnnig using kc
                        file
  -v, --version         show program's version number and exit
  -n THREADS, --num_thread THREADS
                        Number of threads
```

Indexing:
---------

The following code excerpt shows the command to execute the indexing step.
```python
 ./hummingbird.py -r reference.fa -k 27 -o testindex -m index -n 8
 ```
 
Transcript quantification:
--------------------------
 
 The following code excerpt shows the command to execute the quatification step.
 ```python
  ./hummingbird.py -f test.fastq -i testindex -k 27 -o test.kx -t fastq -m express -n 8
  ```
  
Realease Notes:
---------------
 
###version 0.9.7 (9/8/2015):
 
  * Accepts single and paired end fastq(/.gz) files.
  * Accepts .kc files as input.
  * Reference for indexing can be a fasta file or a tab delimited file with probe names and sequence.
  * Rescue algorithm to quantify isoforms does not perform well with isoforms of high overlap. (Fix coming soon in version 0.9.7).
  * Using a cleaned human reference from [DNASU](https://dnasu.org/DNASU/GetCollection.do?collectionName=Human%20hORFeome%20V8.1%20Lentiviral%20collection) will ensure better transcript quantification.
