# Loyola University Chicago COMP483 Mini Project

This is a mini project at graduate level to analyze four strains of Escherichia Coli from an RNA-seq experiment and
compare their gene expression data against the RefSeq genomes from NCBI.


## Software Dependencies/Requirements
* **Linux/Unix**
* **[Python 3+](https://www.python.org/downloads/)**
    * Developed using 3.6.7, requires 3.5+ due to subprocess python module
* **[Biopython](https://biopython.org/) 1.73**
* **[Prokka](https://github.com/tseemann/prokka) 1.13.3**
* **[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 2.3.4.1**
* **[Tophat2](https://ccb.jhu.edu/software/tophat/index.shtml) 2.1.1**
* **[Samtools](http://samtools.sourceforge.net/) 1.7**
* **[Cufflinks](https://cole-trapnell-lab.github.io/cufflinks/) 2.2.1**
* **[Cuffmerge](https://cole-trapnell-lab.github.io/cufflinks/) 1.0.0**
* **[Cuffnorm](https://cole-trapnell-lab.github.io/cufflinks/) 2.2.1**

## Installing and Input Format

### Installing

To run this wrapper, open a terminal session and clone the repository with the command below into a directory of your choosing.
```
git clone https://github.com/tbeaux18/comp483_miniproject.git
```

### Input Files and Format
Text files are only accepted, and must be constructed exactly to specification. Below are examples of the two files required and their structure.
*file containing FTP links that include both the gunzipped FASTA and FEATURE*
*file containing SRA accession IDs*

**SRA Input File**:
The format is very important since the parsers need this exact format and there is no error handling around this.
```
<Strain_name>,<SRA_Accession_ID>
```
**Example**:
```
HM27,SRR1278956\n
```

**FTP Input File**:
The format is very important since the parsers need this exact format and there is no error handling around this.
```
<Strain_name>,<FTP link for fasta .fna.gz>,<FTP link for feature file .txt.gz>\n
```
**Example**:
```
HM27,ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz,ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_feature_count.txt.gz
```

## Main Application Arguments
* `-p/--ftp_links` _text file_
* `-s/--sra_file` _text file_
* `-t/--threads` _integer_

To run the wrapper, clone the repository as described in Installing section and run the wrapper command as below:

**Example**
```
python3 run_wrapper.py -p <path/to/ftpfile.txt> -s <path/to/sra_file.txt> -t 10
```

To make an executable, run
```
chmod a+x run_wrapper.py
```
and then you can run the wrapper as so:
```
./run_wrapper.py -p <path/to/ftpfile.txt> -s <path/to/sra_file.txt> -t 10
```

### Scripts
* `run_wrapper.py`
  * Executed from the top level directory, calls the various scripts below and runs them in linear order.
* `parse_fasta.py`
  * Executed from the scripts directory, one level below the top level directory.
  * Takes the FTP links, downloads the FTP files and parses the FASTA file and outputs to UPEC.log.
* `prokka.py`
  * Executed from the scripts directory, one level below the top level directory.
  * Runs the prokka software with set arguments, and parses the output file and includes them in the UPEC.log.
  * These output files are required for downstream analysis.
* `fastq_dump.py`
  * Executed from the scripts directory, one level below the top level directory.
  * Takes the supplied SRA file and begins to fetch the SRA files and puts them in the *$HOME/ncbi/public/sra* directory on the respective machine. **DO NOT MOVE**
  * fastq_dump is executed and looks for the SRA files in the *$HOME/ncbi/public/sra* directory.
* `tophat2.py`
  * Executed from the scripts directory, one level below the top level directory.
  * Begins with using bowtie2 to construct index files from the FASTA files downloaded. The index files reside in the working directory of where this script is executed.
  * Each alignment takes roughly 2-4 hours depending on the number of threads specified.
* `cufflinks.py`
  * Executed from the scripts directory, one level below the top level directory.
  * Runs the cufflinks, cuffmerge and cuff norm software, and relies on the gtf files produced from prokka.
  * **cuffnorm will fail, there is a sorting error that is unresolved**


## Running the wrapper

After cloning the repository, run the `run_wrapper.py` command with the arguments as specified above under **Main Application Arguments**.

Pipeline will take roughly 12-24 hours in total to run depending on threading.

Place the 2 input files in the same directory as the `run_wrapper.py` script.

The following file tree will be present in the directory during and after running.
  * README.md
  * USER_DESIGNATED_FTP_INPUT_FILE.TXT
  * USER_DESIGNATED_SRA_INPUT_FILE.TXT
  * run_wrapper.py
  * ./test_data
    * ftp_input.txt
    * sra_files.txt
  * ./Timothy_Baker
    * assemblies.txt
    * ftp_links.txt
    * sra_files.txt
    * UPEC.log
    * ./name_prokout
      * name_index.gff
    * ./name_sra
      * SRAname_1.fastq
      * SRAname_2.fastq
    * ./name_tophat
      * accepted_hits.bam
    * ./name_cuff
      * transcripts.gtf
    * ./ncbi_fasta
      * name_FASTA.fna
      * name_FEAT.fna
    * name_index.bt2
    * name.fa
    * ./merged_ecoli
      * merged.gtf
  * cufflinks.py
  * fastq_dump.py
  * parse_fasta.py
  * prokka.py
  * tophat2.py
  * working.log

## Authors

* **Timothy Baker**

## Acknowledgments

* Subashchandrabose for supplying RNA-seq assembly data
* NCBI
* Loyola University Chicago
