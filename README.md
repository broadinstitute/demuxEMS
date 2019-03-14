# scan_snp

Generate cell barcode by SNP matrices.

Bo Li.

* * *

Table of Contents
-----------------

* [Introduction](#introduction)
* [Installation](#install)
	* [Request RedHat7 Server](#server)
	* [Direct use scan_snp](#direct)
	* [Install conda](#conda)
	* [Compile scan_snp](#compile)
* [Usage](#usage)
	* [Inputs](#inputs)
	* [Outputs](#outputs)
	* [Example](#example)

* * *

## <a name="introduction"></a> Introduction

scan_snp generates cell barcode by SNP matrices from cellranger outputs.

## <a name="install"></a> Installation

First, you need to request a RedHat7 Server.

### <a name="server"></a> Request a RedHat7 Server

```
qrsh -q interactive -l h_vmem=4g -l os=RedHat7 -P regevlab
```

or 

```
qrsh -q interactive -l h_vmem=4g -l os=RedHat7 -pe smp 4 -binding linear:4 -P regevlab 
```

to request multiple threads.

### <a name="direct"></a> Direct use scan_snp

```
/ahg/regevdata/users/libo/software/scan_snps/scan_snp
```

### <a name="conda"></a> Install conda environment

This installation instruction has been tested on Ubuntu Linux 18.04.

Suppose your Linux user directory is /users/foo. We will create two folders /users/foo/miniconda3 and /users/foo/software.

Please use the commands below to install scCloud locally via Miniconda:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh .
bash Miniconda3-latest-Linux-x86_64.sh -p /users/foo/miniconda3
mv Miniconda3-latest-Linux-x86_64.sh /users/foo/miniconda3
source ~/.bashrc
conda install -y -c anaconda gcc // install a new version of gcc to use c++11 features 
conda install -y -c anaconda libgcc // install a new version of libgcc
conda install -y -c anaconda libstdcxx-ng // update stdlibc++ symbols
conda install -y -c anaconda cmake // install cmake
conda install -y -c bioconda htslib // install htslib
conda install -y -c anaconda boost // install boost library
```

### <a name="compile"></a> Compile scan_snp

```
git clone --single-branch --branch scan_snps https://github.com/broadinstitute/demuxEMS.git /users/foo/software/scan_snps
cd /users/foo/software/scan_snps
mkdir -p build
cd build
cmake ..
make
```

You should be able to find the executable as ``/users/foo/software/scan_snps/build/src/scan_snp``.

## <a name="usage"></a> Usage

```
Usage: scan_snp input.vcf.gz input.bam output_name [-p number_of_threads]
```

### <a name="inputs"></a> Inputs

File | Description
-----|------------
input.vcf.gz | VCF files contained SNPs. Only SNPs with 'PASS' field will be kept.
input.bam | 10x genomics sorted bam file.
output_name | output file name prefix. 

### Options

Option | Description
-------|------------
-p number_of_threads | Set number of threads used for parsing the BAM file.

### <a name="outputs"></a> Outputs

File | Description
-----|------------
output_name.barcodes.tsv | List all cell barcodes with at least 1 SNP.
output_name.snps.tsv | List all SNPs appeared in at least one read.
output_name.matrix.ref.mtx | Market format matrix, barcode by number of SNP. Counts are the UMIs for the reference allele.
output_name.matrix.alt.mtx | Market format matrix, barcode by number of SNP. Counts are the UMIs for the alternative allele.
output_name.matrix.ref.mtx | Market format matrix, barcode by number of SNP. Counts are the UMIs for errors, the nucleotides that match neither reference nor alternative allele.

### <a name="example"></a> Example:

```
/users/foo/software/scan_snps/build/src/scan_snp MantonBL_demuxlet.vcf.gz possorted_genome_bam.bam test -p 4
```
