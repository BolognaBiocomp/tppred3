# Copyright (C) 2014 Castrense Savojardo
#     Bologna Biocomputing Group
#     University of Bologna, Italy
#     savojard@biocomp.unibo.it
#
#  utils.py - This file is part of TPPRED2
#
#  TPPRED2 is a free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation; either version 3 of the License,
#  or (at your option) any later version.
#
#  TPPRED2 is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with TPPRED2; if not, write to the Free Software Foundation,
#  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.

__all__ = ['compute_hmom',
           'encode_protein',
           'svm_encode_protein',
           'crf_predict',
           'svm_predict']

import os
import tempfile
import shutil
import subprocess as sp
import numpy
import re
import sys
import biocompy.crf as crf

def compute_hmom(fastafile, we):
    m100, m160 = [], []
    hmomfilename = we.createFile("hmom.", ".tsv")
    sp.check_call(['hmoment',fastafile,'-double','-outfile', hmomfilename],
                   stdout = open("/dev/null", 'w'),
                   stderr = open("/dev/null", 'w'))
    for line in open(hmomfilename).readlines()[4:]:
        line = line.split()
        m100.append(float(line[1]))
        m160.append(float(line[2]))
    return m100, m160

def encode_protein(sequence, m100, m160, hscale, window, alphabet):
    dat = []
    pos = 'RK'
    neg = 'ED'
    for i in range(len(sequence)):
        h = p = n = 0.0
        subseq = sequence[max(0, i-window/2):min(len(sequence),i+window/2 + 1)]
        for aa in subseq:
            h += hscale[aa]
            p += float(aa in pos)
            n += float(aa in neg)
        h = h / len(subseq)
        p = p / len(subseq)
        n = n / len(subseq)
        v = [0.0] * len(alphabet)
        try:
            v[alphabet.index(sequence[i])] = 1.0
        except ValueError:
            pass
        if i >= 5 and i < len(sequence) - 5:
            v.extend([h, p, n, m100[i-5], m160[i-5]])
        else:
            v.extend([h, p, n, 0.0, 0.0])
        dat.append(v)
    return numpy.array(dat)

def crf_predict2(inputdat, crfmodel, crfbin, we):
    # ~/biocrf/bin/biocrf -test -w 5 -m crf.model -a 32 -d posterior-viterbi-sum -q post/posterior_label -o crf.positive.pred crf.positive.pred.dat
    crfdat = we.createFile("crf.", ".dat")
    cdofs=open(crfdat,'w')
    for i in range(inputdat.shape[0]):
        for j in range(inputdat.shape[1]):
            cdofs.write("%f " % inputdat[i][j])
        cdofs.write("L\n")
    cdofs.write("\n")
    cdofs.close()
    crfpred = we.createFile("crf.", ".pred")
    crfplabel = we.createFile("crf.",".plabel")
    crfpstate = we.createFile("crf.",".pstate")
    # Run prediction to get P(State)
    sp.check_call([crfbin,
                   '-test', '-w', '5', '-m',
                   crfmodel, '-a', '1', '-d',
                   'posterior-viterbi-max', '-q',
                   crfpstate, '-o', crfpred, crfdat],
                  stderr = open("/dev/null", 'w'),
                  stdout = open("/dev/null", 'w'))
    state_prob = [float(line.split()[45]) for line in open(crfpstate+"_0").readlines()]
    #print state_prob
    # Run prediction to get P(Label) and actual predictions
    sp.check_call([crfbin,
                   '-test', '-w', '5', '-m',
                   crfmodel, '-a', '1', '-d',
                   'posterior-viterbi-sum', '-q',
                   crfplabel, '-o', crfpred, crfdat],
                  stderr = open("/dev/null", 'w'),
                  stdout = open("/dev/null", 'w'))
    label_prob = [float(line.split()[0]) for line in open(crfplabel+"_0").readlines()]
    #print crfpredname
    target = "".join([x.split()[1] for x in open(crfpred).readlines()[:-1]])
    m = re.match('P+', target)
    if not m == None:
        cleavage = m.end()
    else:
        cleavage = 0
    return cleavage, numpy.array(label_prob), numpy.array(state_prob)

def crf_predict(inputdat, crfmodel):
    crfm = crf.CRF()
    crfm.parse(crfmodel)
    target = "".join(crfm.predict([(inputdat, None)], algo = 'posterior-viterbi-sum')[0])
    label_prob = []
    state_prob = []
    for i in range(inputdat.shape[0]):
        label_prob.append(crfm.prob('P', i, label = True))
        state_prob.append(crfm.prob('F3', i, label = False))
    m = re.match('P+', target)
    #print target
    if not m == None:
        cleavage = m.end()
    else:
        cleavage = 0
    return cleavage, numpy.array(label_prob), numpy.array(state_prob)

def elm_predict(seqdat, elmmodel, elmbin, we):
    elmproffile = we.createFile("elm.",".prof")
    elmdatlistfile = we.createFile("elm.",".list.txt")
    elmoutputfile = we.createFile("elm.",".out.txt")
    oprof = open(elmproffile, 'w')
    for i in range(seqdat.shape[0]):
        for j in range(seqdat.shape[1]):
            oprof.write("%1.1f " % seqdat[i][j])
        oprof.write("\n" % seqdat[i][j])
    oprof.close()
    elmprofdir = os.path.dirname(elmproffile)
    olist = open(elmdatlistfile)
    olist.write("%s 0.0\n" % os.path.basename(elmproffile))
    olist.close()
    sp.check_call([elmbin,
                   '-p',
                   elmprofdir,
                   '-i',
                   elmdatlistfile,
                   '-m',
                   elmmodel,
                   '-o',
                   elmoutputfile,
                   '-l', '-t', '-w', '27', '-d', '20', '-b'], stderr = open("/dev/null", 'w'), stdout = open("/dev/null", 'w'))
    o = float(open(elmoutputfile).read().strip().split()[2])
    ret = 'C'
    if o < 0.5:
        ret = 'M'
    return ret

def find_motifs(fastafile, fimoinput, fimobin, fimoth, we):
    fimodir = we.createDir("fimo.")
    sp.check_call([fimobin,
                  "--thresh", str(fimoth),
                  "--oc", fimodir,
                  fimoinput, fastafile], stdout = open("/dev/null", 'w'), stderr = open("/dev/null", 'w'))
    occ = [(int(line.split()[0]),
            int(line.split()[2]) - 1,
            int(line.split()[3]) - 1, float(line.split()[5]),
            float(line.split()[6])) for line in open("%s/fimo.txt" % fimodir).readlines()[1:]]
    return occ

def read_occ_distribution(handle, window):
    mws = {}
    for line in handle.readlines():
        line = line.split()
        M = mws.get(int(line[0]), numpy.zeros(window))
        M[int(line[1])+window/2] = float(line[2])
        mws[int(line[0])] = M
    return mws

def svm_encode_protein(seq, c, w, occs, sp, lp, wdistr, distrfile, we):
    outfile = we.createFile("svm.",".dat")
    ofs=open(outfile,'w')
    dis = read_occ_distribution(open(distrfile), wdistr)
    for i in range(max(0, c - w), min(len(seq), c + w + 1)):
        sc = 0.0
        for o in occs:
            if o[1] >= i - wdistr/2 and o[1] <= i + wdistr/2:
                sc += dis.get(o[0], numpy.zeros(wdistr))[o[1]-i+wdistr/2] * o[3]
        ofs.write("1 1:%f 2:%f 3:%f # pos=%d cleavage=%d\n" % (sc, lp[i], sp[i], i+1, c+1))
    ofs.close()
    return outfile

def svm_predict(svmdatfile, svmmodelfile, svmbin, we):
    outfile = we.createFile("svm.",".pred")
    sp.check_call([svmbin, "-b", "1", svmdatfile, svmmodelfile, outfile],
                   stderr = open("/dev/null", 'w'),
                   stdout = open("/dev/null", 'w'))
    return outfilename
