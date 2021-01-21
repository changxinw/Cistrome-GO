# Cistrome-GO 1.0.0

Functional enrichment analysis of transcription factor ChIP-seq peaks. The package of http://go.cistrome.org.

## Requirements
##### This software is based on **python2.7**, some python packages are required.

- xlmhg==2.4.9
- mne==0.17.0
- scipy

## Installation

```
git clone https://github.com/WChangson/Cistrome-GO.git
cd Cistrome-GO
python setup.py install
```

## Usage
For all the help information, please input

```
cistromego --help
```
#### Solo mode of Cistrome-GO

```
cd test_data
cistromego -g hg38 -b MYOD1.bed
```
#### Ensemble mode of Cistrome-GO

```
cd test_data
cistromego -g hg38 -b peak.bed -e de.txt
```

- -v/--version

Check the version of Cistrome-GO you are using.

- -g/--genome

Specify your genome assembly, select from hg38, hg19, mm10, mm9.

- -b/--bed

Peak bed file path.

- -pn/--peaknumber

Top peak number to use for analysis. Set all to use all the peaks. By default is 10000.

- -d/--decay

Half decay distance to use. Select from auto, 0.5k, 1k, 10k, 50k, 100k, 500k. By default is auto.

- -n/--name

Specify the prefix of your output file.

- -o/--output

Specify the output directory of your result. The program will create it if it does not exists. By default is current directory.

- -dg/--dego

upgenes, downgenes, allgenes for gene ontology. By default is allgenes.

- -e/--expr

Expression file path. If users input this parameter, Cistrome-GO will run as ensemble mode.

- -ei/--exprinfo

Expression file format setting. Column index of geneID, logFoldChange, FDR in differential expression file. By default: 1,3,7.

- -max/--maxgenenumber

Maximum gene number in the GO or kegg terms. By default: 2000.

- -min/--mingenenumber

Minimum gene number in the GO or kegg terms. By default: 10.

- -logfc/--logfccut

logFoldChange of differential expression genes. By default: 1.0.

- -fdr/--fdrcut

FDR cutoff of differential expression genes. By default: 0.05.
