# Copyright (C) 2014 Castrense Savojardo
#     Bologna Biocomputing Group
#     University of Bologna, Italy
#     savojard@biocomp.unibo.it
#
#  config.py - This file is part of TPPRED2
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

import os

TPPRED_ROOT = os.environ.get('TPPRED_ROOT')

__all__ = ['DESCRIPTION',
           'EPILOG',
           'SVM_MODEL_FILE',
           'CRF_MODEL_FILE',
           'CRF_DECODING',
           'AA_ORDER',
           'HWINDOW',
           'CRFWINDOW',
           '']


DESCRIPTION = """
 tppred3.py: Prediction of organelle targeting peptides in proteins.

   Copyright (C) 2015 Castrense Savojardo
   Bologna Biocomputing Group
   University of Bologna, Italy.
   savojard@biocomp.unibo.it
"""

EPILOG = """

"""


SVMMODEL = {'C': os.path.join(TPPRED_ROOT, 'data', 'chloro.svm.model.txt'),
            'M': os.path.join(TPPRED_ROOT, 'data', 'mito.svm.model.txt')}

#CRF_MODEL_FILE = 'data/CRF_w11_s0.05.new.model'
CRF_MODEL_FILE = os.path.join(TPPRED_ROOT, 'data', 'CRF_w11_s0.05.model')

ELM_MODEL_FILE = os.path.join(TPPRED_ROOT, 'data', 'ELM.w27.h54.model.txt')

AA_ORDER = 'VLIMFWYGAPSTCHRKQEND'

FIMO_MOTIF_FILE = {'C': os.path.join(TPPRED_ROOT, 'data', 'chloro.fimo.motif.nr.txt'),
                   'M': os.path.join(TPPRED_ROOT, 'data', 'mito.fimo.motif.nr.txt')}

HWINDOW = 7

CRFWINDOW = 11

SVMWINDOW = {'C': 45, 'M': 29}

FIMO_TH = {'C': '1e-3', 'M': '1e-2'}

MDWINDOW = {'C': 9, 'M' : 9}

KDSCALE = {'A': 1.8, 'C': 2.5, 'U': 2.5, 'D': -3.5,
           'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2,
           'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9,
           'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
           'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9,
           'Y': -1.3, 'X': 0, 'B': 0, 'Z': 0, 'O': 0}

OCCDISTR = {'C': os.path.join(TPPRED_ROOT, 'data', 'chloro.occ.distr.txt'),
            'M': os.path.join(TPPRED_ROOT, 'data', 'mito.occ.distr.txt')}

ELMBIN = os.path.join(TPPRED_ROOT, 'tools', 'elm_predict.py')

FIMOBIN = 'fimo'

SVMBIN = 'svm-predict'

CRFBIN = os.path.join(TPPRED_ROOT, 'tools', 'biocrf-static')

locmap = {"M": ("Mitochondrion", "GO:0005739"),
          "C": ("Chloroplast", "GO:0009507"),
          "N": ("Other", "")}

motifmap = {'M': {1: "RC|[YF][AS]", 2: "SVRx|Y[SA][TS]G"},
            'C': {1: "[VI][RA]|[AC]AAE", 2: "S[VI][RSV]|[CA]A", 3: "[AV]N|A[AM]AG[ED]"}}
