# VIRA: By-Reference Exon and CDS Viral Genome Annotation

[![PyPI version](https://badge.fury.io/py/vira-av.svg)](https://pypi.org/project/vira-av/)
[![GitHub Downloads](https://img.shields.io/github/downloads/alevar/vira/total.svg)](https://github.com/alevar/VIRA/releases/latest)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/licenses/GPL-3.0)

## Introduction

VIRA is a fully-automated protocol for transfering annotations from reference to target genomes, optimized for viral genomes and primarily developed and tested on HIV and SIV genomes.

The method uses both nucleotide and protein sequence information to search for correct alignments between genomes with high degrees of sequence divergence. VIRA is tailored to take advantage of guide protein annotations to further improve the accuracy of alignments and final annotations.

## Publications

Coming soon...

## Documentation

## Installation

### Via PyPI

The easiest way to install VIRA is through PyPI:

```bash
$ pip install vira-av
$ vira --help
```

To uninstall VIRA:

```bash
$ pip uninstall vira-av
```

### Building from source

To build from source, clone the git repository:

```bash
$ git clone https://github.com/alevar/vira.git --recursive
$ cd vira
$ pip install -r requirements.txt
$ pip install .
```

### Requirements

| Requirement | Details |
| ----------- | ------- |
| Language support | Python ≥ 3.6 |
| Dependencies | - [gffread](https://github.com/gpertea/gffread)<br>- [minimap2](https://github.com/lh3/minimap2)<br>- [miniprot](https://github.com/lh3/miniprot)<br>- [snapper](https://github.com/alevar/snapper) |

## Getting started

### Usage

```bash
vira [-h] -a ANNOTATION -g GENOME -t TARGET [-q GUIDE] [-o OUTPUT] [--force-cds] 
     [--gffread GFFREAD] [--minimap2 MINIMAP2] [--miniprot MINIPROT] [--snapper SNAPPER] 
     [--keep-tmp] [--tmp-dir TMP_DIR]
```

### Options

| Option | Description |
| ------ | ----------- |
| `-a, --annotation` | Path to the query GTF/GFF annotation file. |
| `-g, --genome` | Path to the query genome FASTA file. |
| `-t, --target` | Path to the target genome FASTA file. |
| `-q, --guide` | Optional path to the guide annotation file for the target genome. Transcripts and CDS from the guide will be used to validate the annotation. |
| `-o, --output` | Path to the output GTF file. |
| `--force-cds` | Force the CDS from the guide onto the transcript chain, even if that means merging adjacent exons together (can fix alignment artifacts such as spurious introns). If the CDS does not fit the transcript chain, the transcript will be skipped. |
| `--gffread` | Path to the gffread executable. |
| `--minimap2` | Path to the minimap2 executable. |
| `--miniprot` | Path to the miniprot executable. If not set - minimap2 will be used to align nucleotide sequence of the CDS instead. |
| `--snapper` | Path to the snapper executable. |
| `--keep-tmp` | Keep temporary files. |
| `--tmp-dir` | Directory to store temporary files. |

### Help Options

| Option | Description |
| ------ | ----------- |
| `-h, --help` | Prints help message. |

## Example Data

Sample datasets are provided in the "example" directory to test and get familiar with VIRA.

The included example can be run with the following command from the root directory of the repository:

```bash
vira --annotation example/query.gtf --output example/output.gtf --genome example/query.fasta --target example/target.fasta --guide example/guide.gtf
```