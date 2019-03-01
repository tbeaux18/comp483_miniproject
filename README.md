# Loyola University Chicago COMP483 Mini Project

This is a mini project at graduate level to analyze four strains of Escherichia Coli from an RNA-seq experiment and
compare their gene expression data against the RefSeq genomes from NCBI.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

SRA Input File:
The format is very important since the parsers need this exact format and there is no error handling around this.
```
<Strain_name>,<SRA_Accession_ID>
```

Example:
```
HM27,SRR1278956
```

FTP Input File:
The format is very important since the parsers need this exact format and there is no error handling around this.
```
<Strain_name>,<FTP link for fasta .fna.gz>,<FTP link for feature file .txt.gz>\n
```

Example:
```
HM27,ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz,ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_feature_count.txt.gz
```


### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Authors

* **Timothy Baker**

## Acknowledgments

* Subashchandrabose for supplying RNA-seq assembly data
* NCBI
* Loyola University Chicago
