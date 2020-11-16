#!/usr/bin/env python

import getopt
import numpy
import sys

sys.path.append('/home/savojard/BUSCA/tools/tppred3/tppred3')

from modules.biocompy.slfn import *

def printVerbose(level, infolevel, message, nl = True):
    if level >= infolevel:
        if nl:
            print message
        else:
            print message,



def usage():
    u = '''
Copyright (C) 2011 Castrense Savojardo.
savojard@biocomp.unibo.it
elm_predict.py: Extreme Learning Machine for N-to-1 Neural Networks,
                Prediction module.

Usage: elm_predict.py OPTIONS
I/O options:
 -p dir     -> Directory containing examples files. MANDATORY.
 -i file    -> File containing the list of examples. MANDATORY.
 -m file    -> Trained model file. MANDATORY."
 -o file    -> Output prediction file. DEFAULT: examples_file + ".pred" 
 -l         -> Add example names to output. DEFAULT: False.
 -t         -> Transform both targets (if available) and outputs by means
               of a sigmoidal function

Example options:
 -w integer -> Size of the sliding window. DEFAULT: 1.
 -d integer -> Dimensionality of the input vector. DEFAULT: 20.

Network options:
 -k act     -> Type of activation function (DEFAULT: logistic):
                linearr:  one-layer linear network.
                logistic: sigmoidal activation function 1 / (1+exp(-(w*x + b)).
                rbf1:     radial basis activation function exp(-||x-w||/b) 
                          with L1-norm.
                rbf2:     radial basis activation function exp(-||x-w||^2/b^2) 
                          with L2-norm.
 -b         -> Use biased output layer. DEFAULT: 0.
General options:
 -v [0:2]   -> Verbosity level. DEFAULT: 0.
 -h         -> Print this usage.
'''
    print u

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "p:d:o:w:i:m:k:hv:ltb")
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    window = None
    model = None
    ifileName = None
    ofileName = None
    dim = 20
    verbose = 0
    netType = 'logistic'
    datadir = None
    outName = False
    transform = False
    biasedOut = False
    for o, a in opts:
        if o in ("-i"):
            ifileName = a
        elif o in ("-m"):
            model = a
        elif o in ("-p"):
            datadir = a
        elif o in ("-h"):
            usage()
            sys.exit(0)
        elif o in ("-w"):
            window = int(a)
        elif o in("-o"):
            ofileName = a
        elif o in ("-d"):
            dim = int(a)
        elif o in ("-v"):
            verbose = int(a)
        elif o in ("-l"):
            outName = True
        elif o in ("-t"):
            transform = True
        elif o in ("-b"):
            biasedOut = True
        elif o in ("-k"):
            if a == 'logistic':
                pass
            elif a == 'rbf1':
                netType = a
            elif a == 'rbf2':
                netType = a
            elif a == 'linear':
                netType = a
            else:
                print "Unknown activation function", a, ". Using logistic function."
        else:
            assert False, "unhandled option"
    if ifileName == None:
        print 'Please specify an input file'
        usage()
        sys.exit(2)
    
    if datadir == None:
        print 'Please specify the directory that contains data'
        usage()
        sys.exit(2)

    if window == None:
        window = 1
    
    if model == None:
        print 'Please specify a model'
        usage()
        sys.exit(2)
    
    if ofileName == None:
        ofileName = ifileName + '.pred'
    
    try:
        ifile = open(ifileName)
    except IOError:
        print 'Cannot open input file %s' % (ifileName,)
        raise
    
    #printVerbose(verbose, 0, 'Reading the network from the model file ' + model)
    printEvery = 400
    ffn = None
    if netType == 'logistic':
        printVerbose(verbose, 0, 'Sigmoidal network.')
        ffn = SigmoidNetwork(biasedOutput = biasedOut)
    elif netType == 'rbf1':
        printVerbose(verbose, 0, 'Radial Basis Function netwrok with L1-norm.')
        ffn = RBFL1Network(biasedOutput = biasedOut)
    elif netType == 'rbf2':
        printVerbose(verbose, 0, 'Radial Basis Function netwrok with L2-norm.')
        ffn = RBFL2Network(biasedOutput = biasedOut)
    elif netType == 'linear':
        printVerbose(verbose, 0, 'Linear network.')
        ffn = Perceptron(biasedOutput = biasedOut)
    printVerbose(verbose, 0, 'Reading the network from the model file ' + model + ' ...', nl = False)
    ffn.parse(model)
    printVerbose(verbose, 0, 'done.')
    targets = []
    examples = []
    i = 0
    d = (window - 1) / 2
    printVerbose(verbose, 0, 'Starting reading input file ' + ifileName + '...', nl = False)
    printVerbose(verbose, 1, '')
    sys.stdout.flush()
    for line in ifile:
        line = line.split()
        printVerbose(verbose, 2, line[0] + ' ' + line[1])
        examples.append(line[0])
        try:
            prot = open(datadir + '/' + line[0])
        except IOError:
            print 'Cannot open input file', line[0]
            raise
        lines = prot.readlines()
        seq = [numpy.array(x.split()[0:dim], dtype = numpy.float) for x in lines]
        prot.close()
        if len(line) > 1:
            t = float(line[1])
            if transform:
                t = 1. / (1. + numpy.exp(-t))
            targets.append(t)
        padding = [numpy.zeros(dim) for x in range(d)]
        seq = padding + seq + padding
        ex = []
        for j in range(d, len(seq) - d):
            ex.append(numpy.ravel(seq[j - d:j + d + 1]))
        ex = numpy.matrix(ex)
        ffn.compute_H_row(ex)
        if i % printEvery == 0 and i > 0:
            printVerbose(verbose, 1, 'Scanned ' + str(i) + ' sequences')
        i += 1
    printVerbose(verbose, 0, 'done.')
    printVerbose(verbose, 0, 'Number of sequences: ' + str(i))
    printVerbose(verbose, 0, 'Predicting the sample...', nl = False)
    sys.stdout.flush()
    preds = ffn.run()
    if transform:
        preds = 0.5 + numpy.tanh(preds / 2.) / 2.
    printVerbose(verbose, 0, 'done.')
    try:
        ofile = open(ofileName, 'w')
    except IOError:
        print 'Cannot open out file', ofileName, 'for writing'
        sys.exit(2)
    printVerbose(verbose, 0, 'Writing prediction to file ' + ofileName + '...') 
    for j in range(len(preds)):
        ostr = ''
        if outName:
            ostr += '%s\t' % (examples[j],)
        if len(targets) > 0:
            ostr += '%f\t%f\n' % (targets[j], preds[j])
        else:
            ostr += '%f\n' % (preds[j],)
        ofile.write(ostr)
    ofile.close()
    printVerbose(verbose, 0 , 'done.')

if __name__ == "__main__":
    main()
