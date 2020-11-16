tppred2
=======

Mitochondrial targeting peptide prediction

Requirements
============

tppred2 is entirely written in Python. The following packages are required:

 - EMBOSS hmoment executable (http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/hmoment.html)
 - python argparse library (http://docs.python.org/library/argparse.html#module-argparse)
 - python numpy and scipy libraries (http://numpy.scipy.org/)
 - python Biopython library (http://biopython.org/wiki/Main_Page)

To install these packages under Linux debian/ubuntu (you need to be a 
superuser):

sudo apt-get install emboss python-argparse python-numpy python-scipy python-biopython

Configuration
=============

In order to run tppred2 you need to edit the modules/config.py file by setting 
the TPPRED_ROOT variable to point to the tppred2 installation dir. For example, 
if tppred2 has been uncompressed into /home/cas/tppred2 you should set:

TPPRED_ROOT = '/home/cas/tppred2'

Copy the 'bin/tppred2.py' script to a directory in the users' path.

Run the program
===============

Run tppred2 on a multi-FASTA input file as follows:

tppred2.py -f input_fasta_file.fa

By default tppred2 outputs the result to stdout. If you want to write output to file run:

tppred2.py -f input_fasta_file.fa -o out_file.txt

Enjoy!!

Please, report problems or bugs to:  savojard@biocomp.unibo.it


 


