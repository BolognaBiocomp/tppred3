## tppred3 - Organelle-targeting peptide prediction

#### Publication

Savojardo C., Martelli P.L., Fariselli P., Casadio R. [TPpred3 detects and discriminates mitochondrial and chloroplastic targeting peptides in Eukaryotic proteins](http://bioinformatics.oxfordjournals.org/content/31/20/3269), *Bioinformatics* (2015) **31**(20): 3269-3275.

### The TPpred3 Docker image

Image availbale on DockerHub [https://hub.docker.com/r/bolognabiocomp/tppred3](https://hub.docker.com/r/bolognabiocomp/tppred3)

#### Usage of the image

The first step to run TPpred3 Docker container is the pull the container image. To do so, run:

```
$ docker pull bolognabiocomp/tppred3
```

Now the TPpred3 Docker image is installed in your local Docker environment and ready to be used. To show TPpred3 help page run:

```
$ docker run bolognabiocomp/tppred3 -h

usage: tppred3.py -f FASTA_File [-k P,N] [-o out_file]

 tppred3.py: Prediction of organelle targeting peptides in proteins.

   Copyright (C) 2015 Castrense Savojardo
   Bologna Biocomputing Group
   University of Bologna, Italy.
   savojard@biocomp.unibo.it

optional arguments:
  -h, --help  show this help message and exit

OPTIONS:
  -f FILE     Protein sequences in FASTA format. Required.
  -o FILE     Output prediction file. Optional, default: STDOUT.
  -k {P,N}    Protein kingdom: P="Plant" or N="non-Plant"
```
The program accepts three arguments:
- The full path of the input FASTA file containing protein sequences to be predicted;
- The output file where predictions will be stored;
- The kingdom the sequences belong to. You must specify "P" for plant species or "N" for non-plant species. When plant (P) is specified, TPpred3 predicts also chloroplastic targeting signals whereas only mitochondrial signals are predicted for non-plant sequences. The defualt value is plant mode (P).
-
Let's now try a concrete example. First of all, let's downlaod an example sequence from UniProtKB, e.g. the human Enoyl-CoA hydratase, mitochondrial (accession P30084):

```
$ wget http://www.uniprot.org/uniprot/P30084.fasta
```

Now, we are ready to predict the targeting peptide of our input protein. Run:

```
$ docker run -v $(pwd):/data/ bolognabiocomp/tppred3 -f P30084.fasta -o P30084.out -k N
```

In the example above, we are mapping the current program working directory ($(pwd)) to the /data/ folder inside the container. This will allow the container to see the external FASTA file P30084.fasta.
The file P30084.out now contains the TPpred3 prediction:
```
$ cat P30084.out

sp|P30084|ECHM_HUMAN	mTP	0.61 27
```
The first column is the protein accession/id as reported in the input fasta file. The second column report the result of targeting peptide detection. Three different outcomes are possible: mTP (i.e. a mitochondrial-targeting peptide is detected at the N-terminus of the protein), cTP (i.e. a chloroplast-targeting peptide is detected at the N-terminus of the protein, only in plant mode) or Other (i.e. no targeting peptide is detected). The third column is a reliabilty score attached to the previous prediction. Finally, the fourth column reports, for those proteins predicted as having a targeting peptide (eiher mTP or cTP), the predicted position of the cleavage site.

### Install and use TPpred3 from source

Source code available on GitHub at [https://github.com/BolognaBiocomp/tppred3](https://github.com/BolognaBiocomp/tppred3).

#### Installation and configuration

TPpred3 is designed to run on Unix/Linux platforms. The software was written using the Python programming language and it was tested under the Python version 3.

To obtain TPpred3, clone the repository from GitHub:

```
$ git clone https://github.com/BolognaBiocomp/tppred3
```

This will produce a directory “tppred3”. Before running deepsig you need to set and export a variable named TPPRED_ROOT to point to the tppred3 installation dir:
```
$ export TPPRED_ROOT='/path/to/tppred3'
```

Before running the program, you need to install TPpred3 dependencies. We suggest to use Conda (we suggest [Miniconda3](https://docs.conda.io/en/latest/miniconda.html)) create a Python virtual environment and activate it.

To create a conda env for tppred3:

```
$ conda create -n tppred3
```
To activate the environment:

```
$ conda activate tppred3
```

The following Python libraries/tools are required:

- biopython
- emboss
- meme
- libsvm

To install all requirements run the followgin commands:

```
$ conda install meme emboss -c bioconda
$ conda install libsvm -c conda-forge
$ conda install libiconv biopython
```

Now you are able to use tppred3 (see next Section). Remember to keep the environment active.
If you whish, you can copy the “tppred3.py” script to a directory in the users' PATH.

#### Usage

The program accepts three arguments:
- The full path of the input FASTA file containing protein sequences to be predicted;
- The output file where predictions will be stored;
- The kingdom the sequences belong to. You must specify "P" for plant species or "N" for non-plant species. When plant (P) is specified, TPpred3 predicts also chloroplastic targeting signals whereas only mitochondrial signals are predicted for non-plant sequences. The defualt value is plant mode (P).

As an example, run the program on the example FASTA file contained in the folder "example":

```
$ ./tppred3.py -f example/example.fasta -k P -o example/example.out
```

This will run tppred3 on sequences contained in the "example/example.fasta" file, in plant mode and storing the output in the "example/example.out" file.

Once the prediction is done, the output should look like the following:

```
$ cat example/example.out

sp|P30084|ECHM_HUMAN	mTP	0.61	27
sp|P36957|ODO2_HUMAN	mTP	0.91	67
sp|Q86U06|RBM23_HUMAN	Other	1.0	-
sp|P38646|GRP75_HUMAN	cTP	0.83	46
sp|P48047|ATPO_HUMAN	mTP	0.84	23
sp|P48201|AT5G3_HUMAN	mTP	0.73	67
sp|P49411|EFTU_HUMAN	Other	0.72	-
sp|P49448|DHE4_HUMAN	mTP	0.81	53
sp|P55084|ECHB_HUMAN	mTP	0.72	34
```
The first column is the protein accession/id as reported in the input fasta file. The second column report the result of targeting peptide detection. Three different outcomes are possible: mTP (i.e. a mitochondrial-targeting peptide is detected at the N-terminus of the protein), cTP (i.e. a chloroplast-targeting peptide is detected at the N-terminus of the protein, only in plant mode) or Other (i.e. no targeting peptide is detected). The third column is a reliabilty score attached to the previous prediction. Finally, the fourth column reports, for those proteins predicted as having a targeting peptide (eiher mTP or cTP), the predicted position of the cleavage site.

Please, reports bugs to: castrense.savojardo2@unibo.it
