#!/usr/bin/env python

# Copyright (C) 2015 Castrense Savojardo
#     Bologna Biocomputing Group
#     University of Bologna, Italy
#     savojard@biocomp.unibo.it
#
#  tppred3.py - This file is part of TPPRED2
#
#  TPPRED3 is a free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation; either version 3 of the License,
#  or (at your option) any later version.
#
#  TPPRED3 is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with TPPRED2; if not, write to the Free Software Foundation,
#  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.


# This variable should be set to the tppred installation folder


TPPRED_ROOT = os.environ.get('TPPRED_ROOT')
import sys
sys.path.append(TPPRED_ROOT)

import os
import re
try:
    import tempfile
except ImportError:
    sys.stderr.write("Error: Python module tempfile not found.\n")
    sys.exit(1)

try:
    import argparse
except ImportError:
    sys.stderr.write("Error: Python module argparse not found.\n")
    sys.exit(1)

try:
    from Bio import SeqIO
except ImportError:
    sys.stderr.write("Error: Python module biopython not found.\n")
    sys.exit(1)

try:
    import numpy
except ImportError:
    sys.stderr.write("Error: Python module numpy not found.\n")
    sys.exit(1)

import modules.config as config
import modules.utils as utils
import modules.workenv as workenv

def parse_arguments():
    parser = argparse.ArgumentParser(prog = 'tppred3.py',
                                     formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description = config.DESCRIPTION,
                                     epilog = config.EPILOG,
                                     usage = '%(prog)s -f FASTA_File [-k P,N] [-o out_file]')
    options = parser.add_argument_group('OPTIONS')
    options.add_argument('-f',
                         help = 'Protein sequences in FASTA format. Required.',
                         type = argparse.FileType('r'),
                         dest = 'fasta',
                         metavar = 'FILE',
                         required = True)
    options.add_argument('-o',
                         help = 'Output prediction file. Optional, default: STDOUT.',
                         type = argparse.FileType('w'),
                         dest = 'outFile',
                         default = sys.stdout,
                         metavar = 'FILE',
                         required = True)
    options.add_argument('-k',
                         help = 'Protein kingdom: P="Plant" or N="non-Plant"',
                         choices = ['P', 'N'],
                         default = 'P',
                         dest = 'kingdom')

    ns = parser.parse_args()
    return ns

def main():
    args = parse_arguments()
    we = workenv.TemporaryEnv()
    for fasta in SeqIO.parse(args.fasta, 'fasta'):
        l = len(str(fasta.seq))
        seq = str(fasta.seq).replace("U", "C")[:min(l, 160)]
        fastarecfile = we.createFile("seq.", ".fasta")
        fsofs=open(fastarecfile,'w')
        SeqIO.write([fasta], fsofs, 'fasta')
        fsofs.close()
        m100, m160 = utils.compute_hmom(fastarecfile, we)
        inputdat = utils.encode_protein(seq, m100, m160,
                                        config.KDSCALE,
                                        config.HWINDOW,
                                        config.AA_ORDER)
        crf_cleavage, label_prob, state_prob = utils.crf_predict2(inputdat,
                                                                  config.CRF_MODEL_FILE,
                                                                  config.CRFBIN,
                                                                  we)
        if crf_cleavage > 0:
            crf_prob = numpy.mean(label_prob[:crf_cleavage+1])
            if args.kingdom == "P":
                seqdat = inputdat[0:crf_cleavage+1, 0:20]
                organelle = utils.elm_predict(seqdat, config.ELM_MODEL_FILE,
                                              config.ELMBIN, we)
            else:
                organelle = "M"
            motifoccs = utils.find_motifs(fastarecfilename,
                                          config.FIMO_MOTIF_FILE[organelle],
                                          config.FIMOBIN,
                                          config.FIMO_TH[organelle],
                                          we)
            svmdatfile = utils.svm_encode_protein(seq,
                                                  crf_cleavage - 1,
                                                  config.SVMWINDOW[organelle]/2,
                                                  motifoccs,
                                                  state_prob,
                                                  label_prob,
                                                  config.MDWINDOW[organelle],
                                                  config.OCCDISTR[organelle],
                                                  we)
            svmpredfile = utils.svm_predict(svmdatfile,
                                            config.SVMMODEL[organelle],
                                            config.SVMBIN,
                                            we)
            pred = []
            dlines = open(svmdatfile).readlines()
            plines = open(svmpredfile).readlines()[1:]
            for i in range(len(dlines)):
                lined = dlines[i].split()
                linep = plines[i].split()
                if int(linep[0]) == 1:
                    pred.append((int(lined[5].split("=")[1]), float(linep[1])))
            if len(pred) > 0:
                cleavage, prob = sorted(pred, key = lambda x: x[1])[-1]
                prob = crf_prob
                cleavage = str(cleavage)
                source = "SVM"
            else:
                cleavage, prob = str(crf_cleavage), crf_prob
                source = "CRF"
        else:
            cleavage = "No targeting signal detected"
            organelle = "N"
            prob = numpy.mean(1.0 - label_prob[:min(30, len(seq))])
            motifoccs = []
            source = "CRF"
        print(fasta.id, organelle, crf_cleavage, cleavage, prob, source,
              sep="\t",file=ns.outFile)
        try:
            c = int(cleavage)
            occs = filter(lambda x: x[1]>=c-config.MDWINDOW[organelle]/2 and x[1]<=c+config.MDWINDOW[organelle]/2, motifoccs)
        except:
            occs = []
    #we.destroy()
    return 0

if __name__ == "__main__":
    r = main()
    sys.exit(r)
