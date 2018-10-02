# ARGDIT
## Overview
The Antimicrobial Resistance Gene Data Integration Toolkit (ARGDIT) consists of two main tools and three utilities for users to perform data validation and integration on antimicrobial resistance gene (ARG) databases. Basically it allows users to validate the fidelity of an ARG (or other coding sequence/protein) database to the coding sequence/protein information from the NCBI repositories [1], and to merge multiple validated databases into a single ARG database. It also supports automated re-annotating the output ARG sequences with NCBI sequence information, as well as sequence classification (i.e. predicting the class labels of the database sequences) according to a schema database, which is another ARG database containing classified sequences.

Note that although the default translation table used is for bacteria, other translation tables can also be used for non-bacterial coding sequence databases by specifying a different genetic code. Refer to the usage details of the two tools below.

## Main tools
* ARG database validation tool (check_arg_db.py)
* ARG database integration tool (merge_arg_db.py)

## Utilities
* Database sequence replacement utility (replace_db_seqs.py)
* UniProt identifier to NCBI protein accession number conversion utility (convert_id_uniprot_to_ncbi.py)
* Database diff utility (seq_db_diff.py)

## Sample consolidated databases
MEGARes [2] is a consolidated ARG database created from four ARG databases: ARG-ANNOT [3], CARD [4-6], ResFinder [7], and Lahey beta-lactamase archive [8] (NCBI BioProject PRJNA305729 [1]). However, its latest version was released in Dec 2016 since which its source databases have been updated. Therefore, to demonstrate the benefit of ARGDIT in highly automating the ARG database update process, a consolidated database that emsembles MEGARes database has been created from its sources. The source databases were first validated using ARGDIT validation tool, and were then corrected manually according to the validation results. The validated databases were then merged by ARGDIT integration tool into a single database, which was subsequently annotated using NCBI information [1]. The consolidated database is available in two versions:

<<<<<<< HEAD
* argdit_nt_db.fa
* argdit_nt_db_no_class.fa

The annotations (i.e. sequence headers) of the first database include predicted (but unverified) ARG class labels defined by MEGARes and the second database does not include any sequence class information. Both versions can be found in the directory "sample_integrated_dbs".
=======
* ARGDIT_integrated_DB_with_annotation_and_ontology_class.fa
* ARGDIT_integrated_DB_with_annotation.fa

The first database contains predicted MEGARes ARG ontology labels in its annotation and the second does not. Both versions and their creation logs can be found in the directory "sample_integrated_dbs".
>>>>>>> a88c3e719fe634c63395d871cab2c1bec2f35e59

## Database eligibility
In order to use the data validation and integration tool, the ARG (or other coding/protein sequence) database must be

1. In FASTA format
2. Every FASTA sequence header must contain an NCBI nucleotide/protein accession number. Uniprot ID is an alternative for protein accession number for protein sequence database (by converting the Uniprot IDs to protein NCBI accession numbers with the conversion utility provided)
3. Every sequence must either be a coding sequence (i.e. will be translated to protein product) or a protein sequence
4. ARG/Sequence class information, if any, must occupy at least one individual field in the sequence headers, in which all the fields are separated by the "|" symbol

## Important notice
All data retrieval of NCBI repositories are performed through NCBI Entrez Programming Utilities. Before using ARGDIT it is very important for every user to read its [guidelines and requirements](https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen) and avoid overwhelming the NCBI servers according to the guidelines. Based on these requirements, users are required to provide their contact email addresses (see the Installation section) so that NCBI may attempt contact before **blocking the abusing access**. Although this email address is intended for the software developers, it is more appropriate for the users to fill in their own so that they can be notified by NCBI.

## Pre-requisites
The followings must be installed for the core ARGDIT operations:
1. Python version 3.5 or higher
2. BioPython version 1.70 or higher

If sequence classification or class outlier sequence detection is required, then the followings must be installed:
1. MUSCLE version 3.8.31 or higher ([Link](https://www.drive5.com/muscle/downloads.htm))
2. OD-Seq ([Link](http://www.bioinf.ucd.ie/download/od-seq.tar.gz))
3. HMMER3 version 3.1b2 or higher ([Link](http://hmmer.org/download.html))

## Installation
No installation is required. Make sure all the third-party software in pre-requisites are in the system path.

In order to access the NCBI repositories, users must provide their own contact email addresses along with their access requests. Fill in your contact email address under the "Entrez" section in the configuration file (config.ini):

[Entrez]
Email = (your contact email address)

## Usage
### **Database validation tool**
**Command**  
./check_arg_db.py [optional arguments] *seq_db_path*

**Mandatory argument**  

| Argument name | Description                                 |
| ------------- | ------------------------------------------- |
| seq_db_path   | nucleotide/protein database FASTA file path |

**Optional arguments**  

| Argument name                 | Description                                                |
| ----------------------------- | ---------------------------------------------------------- |
| -f/--fields FIELD_NUMS        | sequence class label field numbers FIELD_NUMS for class outlier sequence detection, e.g. -f 4-5, -f ~1-~3                                            |
| -r/--refine                   | export refined DNA sequences                               |
| -c/--geneticcode GENETIC_CODE | genetic code to specify which translation table to be used |
| -e/--exportlog                | export validation results and process log                  |
| -h/--help                     | show help message and exit                                 |

**Description**  
check_arg_db.py performs ARG database validation. The --refine option allows the tool to trim at most two spurious bases before the start codon or after the stop codon, and export the trimmed sequences into an individual file specified by the tool. To perform ARG class outlier sequence detection, specify the ARG class fields after the --fields option. For example, the hierarchical ARG class information can be extracted from MEGARes database by "-f ~1-~3". The --geneticcode option overrides the default genetic code specified in the configuration file. The genetic code represents the translation table to be used when translating the DNA sequences. As a result, sequence databases for organisms other than bacteria can also be validated. Refer to [here](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for the genetic codes representing different translation tables. The validation results and process log are printed to stdout (i.e. screen) by default, and by specifying the --exportlog option they will be sent to a .log file in the same directory as the database file.

### **Database integration tool**
**Command**  
./merge_arg_db.py [optional arguments] -o *OUTPUT_SEQ_DB_PATH seq_db_paths*

**Mandatory arguments**  

| Argument name         | Description                                              |
| --------------------- | -------------------------------------------------------- |
| -o OUTPUT_SEQ_DB_PATH | specify the output database file path OUTPUT_SEQ_DB_PATH |
| seq_db_paths          | nucleotide/protein database FASTA file paths             |
	
**Optional arguments**  

| Argument name                         | Description                                         |
| ------------------------------------- | --------------------------------------------------- |
| -s/--schema SCHEMA_DB_PATH FIELD_NUMS | specify the schema database SCHEMA_DB_PATH and class label field numbers FIELD_NUMS to perform sequence class prediction                           |
| -a/--annotate                         | perform automated re-annotation using NCBI repository information                                                                                   |
| -p/--protein                          | export protein sequences                            |
| -r/--redundant                        | allow redundant sequences                           |
| -c/--geneticcode GENETIC_CODE         | genetic code to specify which translation table to be used                                                                                          |
| -e/--exportlog                        | export integration results and process log          |
| -h/--help                             | show help message and exit                          |

**Description**  
merge_arg_db.py performs integration of multiple ARG databases. The --annotate option performs re-annotation of the sequences in the output database. By specifying the schema database file path and the ARG class fields after the --schema option, the class labels of the output sequences will be predicted. However, note that the protein sequences of the schema database are not validated here, so it is advised to validate them using the validation tool. The --geneticcode option overrides the default genetic code specified in the configuration file. The genetic code represents the translation table to be used when translating the DNA sequences. As a result, sequence databases for organisms other than bacteria can also be consolidated. Refer to [here](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for the genetic codes representing different translation tables. By default only non-redundant sequences are exported, but this can be overridden with the --redundant option. The tool provides the --protein option to translate all DNA sequences to protein sequences.

### **Sequence replacement utility**
**Command**  
./replace_db_seqs.py [optional argument] *seq_db_path replace_seq_file_path output_seq_db_path*

**Mandatory arguments**  

| Argument name         | Description                                 |
| --------------------- | ------------------------------------------- |
| seq_db_path           | nucleotide/protein database FASTA file path |
| replace_seq_file_path | FASTA file path for replacement sequences   |
| output_seq_db_path	| output database file path                   |
	
**Optional argument**  

| Argument name | Description                |
| ------------- | -------------------------- |
| -h/--help	    | show help message and exit |

**Description**  
By matching identical FASTA headers of the sequences, this utility replaces the sequences in the database FASTA file with those in the replacement sequence file. The database sequences, no matter replaced or not, are exported to the output database file specified by the user.

### **Uniprot ID conversion utility**
**Command**  
./convert_id_uniprot_to_ncbi.py [optional argument] *seq_db_path output_seq_db_path*

**Mandatory arguments**  

| Argument name      | Description                                                       |
| ------------------ | ----------------------------------------------------------------- |
| seq_db_path	     | FASTA file path for protein database with Uniprot IDs             |
| output_seq_db_path | output file path for protein database with converted NCBI protein accession no.                                                                            |
	
**Optional argument**  

| Argument name | Description                |
| ------------- | -------------------------- |
| -h/--help	    | show help message and exit |

**Description**  
This tool queries the UniProt database for the Uniprot ID to NCBI protein accession number mappings, and then replaces the UniProt IDs in the sequence headers by the mapped protein accession numbers. The processed sequences are exported to the output database file specified by the user.

### **Database diff utility**
**Command**  
./seq_db_diff.py [optional arguments] *seq_db_path old_seq_db_path*

**Mandatory arguments**  

| Argument name   | Description                                       |
| --------------- | ------------------------------------------------- |
| seq_db_path	  | nucleotide/protein database FASTA file path       |
| old_seq_db_path | database FASTA file path for the previous release |
	
**Optional argument**  

| Argument name     | Description                                          |
| ----------------- | ---------------------------------------------------- |
| -m/--mode {h,s,b} | diff mode (h:header; s:sequence; b:both [default])   |
| -o/--old          | export also sequences found in old database only     |
| -f/--forward      | do not check reverse complement nucleotide sequences |
| -h/--help	        | show help message and exit                           |

**Description**  
This utility compares the sequence database with its previous version according to the diff mode specified. Header (h) mode performs diff by comparing sequence annotation header only (hence ignoring residue sequence); sequence (s) mode performs diff by comparing residue sequence only (hence ignoring sequence annotation header); both (b) mode performs diff by comparing both annotation header and residue sequence. A FASTA file storing all identical sequences (according the diff mode) is generated, as well as another FASTA file storing those that appear in the current database only. The --old option can export a FASTA file storing sequences that appear in the old database only. By default, when the diff mode is other than header mode and the database sequences are DNA sequences, residue sequence comparison is also performed with reverse complement (i.e. a DNA sequence S1 is compared with another DNA sequence S2 as well as the reverse complement of S2). However, the reverse complement comparison can be suspended by the --forward option.

## Known issues
It is sometimes (but not often) possible to have incomplete data retrieval from NCBI repositories due to server-side issues such as heavy workload. This means information for some nucleotide/protein accession numbers cannot be retrieved at the moment; the outcome is like these accession numbers are not present in the NCBI repositories. When many sequences are spuriously reported as having their accession numbers not found and/or sequence mismatches, it is advised to try using the database validation and the integration tools later.

## Citation
TBC

## References
[1] NCBI Resource Coordinators. Database Resources of the National Center for Biotechnology Information. Nucleic Acids Research 2017;45(D1):D12-D17.  
[2] Lakin, S.M., et al. MEGARes: an antimicrobial resistance database for high throughput sequencing. Nucleic Acids Research 2017;45(D1):D574-D580.  
[3] Gupta, S.K., et al. ARG-ANNOT, a New Bioinformatic Tool To Discover Antibiotic Resistance Genes in Bacterial Genomes. Antimicrobial Agents and Chemotherapy 2014;58(1):212-220.  
[4] Jia, B., et al. CARD 2017: expansion and model-centric curation of the comprehensive antibiotic resistance database. Nucleic Acids Research 2017;45(D1):D566-D573.  
[5] McArthur, A.G., et al. The Comprehensive Antibiotic Resistance Database. Antimicrobial Agents and Chemotherapy 2013;57(7):3348-3357.  
[6] McArthur, A.G. and Wright, G.D. Bioinformatics of antimicrobial resistance in the age of molecular epidemiology. Current Opinion in Microbiology 2015;27(Supplement C):45-50.  
[7] Zankari, E., et al. Identification of acquired antimicrobial resistance genes. Journal of Antimicrobial Chemotherapy 2012;67(11):2640-2644.  
[8] Bush, K. and Jacoby, G.A. Updated Functional Classification of β-Lactamases. Antimicrobial Agents and Chemotherapy 2010;54(3):969-976.  