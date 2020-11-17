'''
Biocompy: Python Library for Biocomputing
Single Hidden Layer Feed-Forward Network module

@author    : Castrense Savojardo
@copyright : 2011 Castrense Savojardo
@contact   : savojard@biocomp.unibo.it

'''
import numpy

class SLFNSigmoidSKLWrapper(object):
    def __init__(self, n_features = 1, n_hidden = 10, biasedOutput = True):
        self.n_features = n_features
        self.slfn = SigmoidNetwork(d = n_features, h = 10, biasedOutput = biasedOutput)

    def fit(self, X, y):
        for i in range(X.shape[0]):
            xrow = numpy.matrix(X[i])
            self.slfn.compute_H_row(xrow)
        self.slfn.trainELM(numpy.transpose(numpy.matrix(y)))

    def predict(self, X):
        self.slfn.H = numpy.matrix([])
        for i in range(X.shape[0]):
            xrow = numpy.matrix(X[i])
            self.slfn.compute_H_row(xrow)
        p = self.slfn.run()
        return numpy.ravel(numpy.array(p))

    def get_params(self, deep=True):
        return {'n_hidden': self.slfn.h, 'n_features': self.slfn.d, 'biasedOutput': self.slfn.biasedOutput}

    def set_params(self, **params):
        self.slfn = SigmoidNetwork(d = self.n_features, h = params.get('n_hidden', 10), biasedOutput = params.get('biasedOutput', True))
        return self

class SLFNRBFL1SKLWrapper(object):
    def __init__(self, n_features = 1, n_hidden = 10, biasedOutput = True):
        self.n_features = n_features
        self.slfn = RBFL1Network(d = n_features, h = 10, biasedOutput = biasedOutput)

    def fit(self, X, y):
        for i in range(X.shape[0]):
            xrow = numpy.matrix(X[i])
            self.slfn.compute_H_row(xrow)
        self.slfn.trainELM(numpy.transpose(numpy.matrix(y)))

    def predict(self, X):
        self.slfn.H = numpy.matrix([])
        for i in range(X.shape[0]):
            xrow = numpy.matrix(X[i])
            self.slfn.compute_H_row(xrow)
        p = self.slfn.run()
        return numpy.ravel(numpy.array(p))

    def get_params(self, deep=True):
        return {'n_hidden':self.slfn.h, 'n_features':self.slfn.d, 'biasedOutput': self.slfn.biasedOutput}

    def set_params(self, **params):
        self.slfn = RBFL1Network(d = self.n_features, h = params.get('n_hidden', 10), biasedOutput = params.get('biasedOutput', True))
        return self

class SLFNRBFL2SKLWrapper(object):
    def __init__(self, n_features = 1, n_hidden = 10, biasedOutput = True):
        self.n_features = n_features
        self.slfn = RBFL2Network(d = n_features, h = 10, biasedOutput = biasedOutput)

    def fit(self, X, y):
        for i in range(X.shape[0]):
            xrow = numpy.matrix(X[i])
            self.slfn.compute_H_row(xrow)
        self.slfn.trainELM(numpy.transpose(numpy.matrix(y)))

    def predict(self, X):
        self.slfn.H = numpy.matrix([])
        for i in range(X.shape[0]):
            xrow = numpy.matrix(X[i])
            self.slfn.compute_H_row(xrow)
        p = self.slfn.run()
        return numpy.ravel(numpy.array(p))

    def get_params(self, deep=True):
        return {'n_hidden':self.slfn.h, 'n_features':self.slfn.d, 'biasedOutput': self.slfn.biasedOutput}

    def set_params(self, **params):
        self.slfn = RBFL2Network(d = self.n_features, h = params.get('n_hidden', 10), biasedOutput = params.get('biasedOutput', True))
        return self

class SLFNError(Exception):
    """
    Generic SLFN exception.
    """
    pass

class SLFN(object):
    """
    A class for Single hidden-Layer Feed-forward neural Network.
    Both 1-to-1 and 1-to-N SLFNs are supported.
    """

    N_TO_ONE = 0
    ONE_TO_ONE = 1
    LEFT_RIGHT = 2
    ALL_PAIRS = 3

    def __init__(self,
                 d = None,
                 h = None,
                 wrange = (-1, 1),
                 biasrange = (-1, 1),
                 biasedOutput = False,
                 mapping = 0,
                 swap = False,
                 swapFile = './slfn.H.swap'):
        """
        Construct a new SLFN object.

        @type d: int
        @param d: The dimensionality of input layer.
        @type h: int
        @param h: The number of neurons in the hidden layer.
        @type wrange: tuple or list
        @param wrange: The interval of choice of input weights.
        @type biasrange: tuple or list
        @param biasrange: The interval of choice of hidden node biases.
        @type biasedoutput: boolean
        @param biasedoutput: If True consider biased output layers.
        """
        self.o = None
        self.biasedOutput = biasedOutput
        self.betas = None
        self.H = numpy.matrix([])
        self.mapping = mapping
        if self.mapping == SLFN.N_TO_ONE:
            self.hidden_layer_output = self.n_to_one_hidden_layer
        elif self.mapping == SLFN.ONE_TO_ONE:
            self.hidden_layer_output = self.one_to_one_hidden_layer
        elif self.mapping == SLFN.LEFT_RIGHT:
            self.hidden_layer_output = self.left_right_hidden_layer
        elif self.mapping == SLFN.ALL_PAIRS:
            self.hidden_layer_output = self.all_pairs_hidden_layer

        if not d == None and not h == None:
            assert isinstance(d, int)
            assert isinstance(h, int)
            assert d > 0
            assert h > 0
            assert (isinstance(wrange, tuple) or
                    isinstance(wrange, list))
            assert (isinstance(biasrange, tuple) or
                    isinstance(biasrange, list))
            assert (wrange[1] > wrange[0])
            assert (biasrange[1] > biasrange[0])

            self.d = d
            self.h = h
            # Generate random hidden node weights
            weights = numpy.random.uniform(wrange[0],
                                           wrange[1],
                                           (self.d, self.h))
            biases = numpy.random.uniform(biasrange[0],
                                          biasrange[1],
                                          (1, self.h))
            self.w = numpy.matrix(numpy.concatenate((weights, biases)))
        else:
            self.d = None
            self.h = None

        self.swap = swap
        if self.swap:
            try:
                self.swapFile = open(swapFile, 'w')
                self.swapFileName = swapFile
            except IOError:
                raise SLFNError('Error opening swap file ' + swapFile)

    def __str__(self):
        s = 'mapping ' + str(self.mapping) + '\n'
        s += 'h ' + str(self.h) + '\n'
        s += 'd ' + str(self.d) + '\n'
        s += 'o ' + str(self.o) + '\n'
        s += 'b'
        for b in numpy.ravel(self.betas):
            s += ' %1.15f' % (b,)
        s += '\n'
        s += 'w'
        for we in numpy.ravel(self.w):
            s += ' %1.15f' % (we,)
        s += '\n'
        return s

    def set_first_layer_weights(self, w):
        assert (self.w.shape[0] == self.d + 1 and self.w.shape[1] == self.h)
        self.w = w


    def parse(self, filename):
        '''
        Initialize the network by parsing a model file

        @type filename: str
        @param filename: A file containing a trained model.
        '''

        self.mapping = SLFN.N_TO_ONE

        try:
            ifile = open(filename)
        except IOError:
            raise SLFNError('Cannot open model file')
        for line in ifile:
            line = line.split()
            if line[0] == 'b':
                assert isinstance(self.o, int)
                assert isinstance(self.h, int)
                assert self.o > 0
                assert self.h > 0
                try:
                    b = [float(x) for x in line[1:]]
                except ValueError:
                    raise SLFNError('Invalid beta in the model file')

                if (self.mapping == SLFN.LEFT_RIGHT or
                    self.mapping == SLFN.ALL_PAIRS):
                    try:
                        self.betas = numpy.reshape(numpy.array(b), (self.h * 3, self.o))
                    except ValueError:
                        raise SLFNError('Error while reshaping beta')
                else:
                    try:
                        self.betas = numpy.reshape(numpy.array(b), (self.h, self.o))
                    except ValueError:
                        raise SLFNError('Error while reshaping beta')
                self.betas = numpy.matrix(self.betas)
            elif line[0] == 'mapping':
                self.mapping = int(line[1])
                try:
                    assert (self.mapping == SLFN.N_TO_ONE or
                            self.mapping == SLFN.ONE_TO_ONE or
                            self.mapping == SLFN.LEFT_RIGHT or
                            self.mapping == SLFN.ALL_PAIRS)
                except AssertionError:
                    self.mapping = SLFN.N_TO_ONE

            elif line[0] == 'h':
                try:
                    self.h = int(line[1])
                except ValueError:
                    raise SLFNError('Invalid h value in the model file')
                assert self.h > 0
            elif line[0] == 'd':
                try:
                    self.d = int(line[1])
                except ValueError:
                    raise SLFNError('Invalid d value in the model file')
                assert self.d > 0
            elif line[0] == 'w':
                assert isinstance(self.h, int)
                assert isinstance(self.d, int)
                assert self.h > 0
                assert self.d > 0
                try:
                    w = [float(x) for x in line[1:]]
                except ValueError:
                    raise SLFNError('Invalid weight in the model file')
                try:
                    self.w = numpy.reshape(numpy.array(w), (self.d + 1, self.h))
                except ValueError:
                    raise SLFNError('Error while reshaping w')
                self.w = numpy.matrix(self.w)
            elif line[0] == 'o':
                try:
                    self.o = int(line[1])
                except ValueError:
                    raise SLFNError('Invalid o value in the model file')
                assert self.o > 0

    def activation(self, x):
        """
        Compute the default sigmoidal hidden layer activation function.

        @type x: numpy.matrix of L rows and d columns
        @param x: The input matrix.

        @return: sig(x * w + b)
        @rtype: numpy.matrix L x h
        """
        assert isinstance(self.d, int)
        assert isinstance(self.h, int)
        assert (self.d > 0)
        assert (self.h > 0)
        assert (self.w.shape[0] == self.d + 1 and self.w.shape[1] == self.h)
        assert (x.shape[1] + 1 == self.w.shape[0])
        try:
            L = numpy.shape(x)[0]
        except IndexError:
            L = 1
        bias = numpy.matrix(numpy.ones(L)).T
        row = numpy.array(numpy.concatenate((x, bias), axis = 1) * self.w)
        # row : L x h
        row = 0.5 + numpy.tanh(row / 2.) / 2.
        return row

    def _compile(self, mat):
        assert (mat.shape[0] > 0 and mat.shape[1] > 0)
        return numpy.mean(mat, axis = 0)

    def n_to_one_hidden_layer(self, ex):
        assert isinstance(self.d, int)
        assert isinstance(self.h, int)
        assert (self.d > 0)
        assert (self.h > 0)
        assert (self.w.shape[0] == self.d + 1 and self.w.shape[1] == self.h)
        assert(ex.shape[1] + 1 == self.w.shape[0])
        row = self.activation(ex)
        # computing the mean vector m : 1 x h
        row = self._compile(row)
        if self.biasedOutput:
            row[-1] = 1.
        return row

    def one_to_one_hidden_layer(self, ex):
        assert isinstance(self.d, int)
        assert isinstance(self.h, int)
        assert (self.d > 0)
        assert (self.h > 0)
        assert (self.w.shape[0] == self.d + 1 and self.w.shape[1] == self.h)
        assert(ex.shape[1] + 1 == self.w.shape[0])
        row = self.activation(ex)
        if self.biasedOutput:
            #TODO
            pass
        return row

    def _forward(self, row):
        fw = numpy.zeros_like(row)
        for j in range(1, row.shape[0]):
            fw[j] = (fw[j - 1] * (j - 1) + row[j - 1]) / j
        return fw


    def _backward(self, row):
        bw = numpy.zeros_like(row)
        l = row.shape[0]
        for j in range(l - 2, -1, -1):
            bw[j] = (bw[j + 1] * (l - 1 - (j + 1)) + row[j + 1]) / (l - 1 - j)
        return bw

    def left_right_hidden_layer(self, ex):

        row = self.activation(ex)
        fw = self._forward(row)
        bw = self._backward(row)
        res = numpy.hstack((fw, row, bw))
        if self.biasedOutput:
            #TODO
            pass
        return res

    def all_pairs_hidden_layer(self, ex):

        nobj = ex.shape[0]
        pairs = {}
        for i in range(nobj - 1):
            for j in range(i + 1, nobj):
                #print numpy.matrix(numpy.hstack((ex[i], ex[j]))).shape
                center = self.activation(numpy.matrix(numpy.hstack((ex[i], ex[j]))))
                pairs[(i, j)] = center

        res = None


        for p in pairs:
            r = None
            l = None
            c = pairs[p]
            for i in range(nobj - 3):
                c = numpy.vstack((c, pairs[p]))
            for i in range(nobj):
                if not (i == p[1] or i == p[0]):
                    if l == None:
                        l = pairs[(min(p[0], i), max(p[0], i))]
                    else:
                        l = numpy.vstack((l, pairs[(min(p[0], i), max(p[0], i))]))
                    if r == None:
                        r = pairs[(min(p[1], i), max(p[1], i))]
                    else:
                        r = numpy.vstack((r, pairs[(min(p[1], i), max(p[1], i))]))
            if res == None:
                res = self._compile(numpy.matrix(numpy.hstack((l, c, r))))
            else:
                res = numpy.vstack((res, self._compile(numpy.matrix(numpy.hstack((l, c, r))))))
        return res

    def compute_H_row(self, ex):
        """
        Compute a row of the hidden layer output matrix for a given example.

        @type ex: numpy.matrix
        @param x: The input example.
        """
        row = self.hidden_layer_output(ex)
        # appending the row to the H matrix
        if not self.swap:
            try:
                self.H = numpy.vstack((self.H, row))
            except ValueError:
                self.H = numpy.matrix(row)
        else:
            try:
                for r in row:
                    for el in r:
                        self.swapFile.write('%f ' % (el,))
                    self.swapFile.write('\n')
            except IOError:
                raise SLFNError('Error writing swap file ' + self.swapFileName)


    def trainELM(self, targets):
        """
        Train the network with the Extreme Learning Machine algorithm.

        @type targets: numpy.matrix N x o
        @param targets: The targets values to train the network.
        """
        assert isinstance(self.h, int)
        assert (targets.shape[0] > 0 and targets.shape[1] > 0)
        assert (self.h > 0)
        if not self.swap:
            assert (self.H.shape[1] == self.h or
                    (self.H.shape[1] == (self.h * 3) and
                     (self.mapping == SLFN.LEFT_RIGHT or
                      self.mapping == SLFN.ALL_PAIRS)))
            assert (self.H.shape[0] > 0)
            assert (targets.shape[0] == self.H.shape[0])

        try:
            self.o = numpy.shape(targets)[1]
        except IndexError:
            self.o = 1
        # Compute the output weights beta = Hdagger * T
        if not self.swap:
            self.betas = (numpy.linalg.pinv(self.H.transpose() * self.H)
                          * self.H.transpose() * targets)
        else:
            try:
                self.swapFile.close()
            except:
                pass
            try:
                swf = open(self.swapFileName, 'r')
            except IOError:
                raise SLFNError('Error reading swap file ' + self.swapFileName)
            # Directly compute Hsq = H^transpose * H
            dim = self.h
            if (self.mapping == SLFN.LEFT_RIGHT or
                self.mapping == SLFN.ALL_PAIRS):
                dim = dim * 3

            Hsq = numpy.matrix(numpy.zeros((dim, dim)))
            t = numpy.matrix(numpy.zeros((dim, self.o)))
            j = 0
            for line in swf.xreadlines():
                # Compute H^T * H
                vect = numpy.matrix(numpy.array([float(x) for x in line.split()]))
                temp = vect.transpose() * vect
                Hsq += temp
                # Compute H^T * targets
                temp = vect.transpose() * targets[j]
                t += temp
                j += 1

            self.betas = numpy.linalg.pinv(Hsq) * t
            self.swapFile.close()

    def run(self):
        """
        Compute network outputs from the (pre)computed H matrix.

        @rtype: numpy.matrix N x o
        @return: the network outputs.
        """
        assert isinstance(self.o, int)
        assert isinstance(self.h, int)
        assert (self.o > 0)
        assert (self.h > 0)
        if not self.swap:
            assert (self.H.shape[1] == self.h or
                    (self.H.shape[1] == (self.h * 3) and
                     (self.mapping == SLFN.LEFT_RIGHT or
                      self.mapping == SLFN.ALL_PAIRS)))
            assert (self.H.shape[0] > 0)
        assert (self.betas.shape[0] == self.h or
                (self.betas.shape[0] == (self.h * 3) and
                 (self.mapping == SLFN.LEFT_RIGHT or
                  self.mapping == SLFN.ALL_PAIRS)))
        assert (self.betas.shape[1] == self.o)

        if not self.swap:
            o = self.H * self.betas
        else:
            try:
                self.swapFile.close()
            except:
                pass
            try:
                swf = open(self.swapFileName, 'r')
            except IOError:
                raise SLFNError('Error reading swap file ' + self.swapFileName)
            o = numpy.matrix([])
            for line in swf.xreadlines():
                vect = numpy.matrix(numpy.array([float(x) for x in line.split()]))
                temp = vect * self.betas
                try:
                    o = numpy.vstack((o, temp))
                except ValueError:
                    o = temp
            self.swapFile.close()
        return o

class Perceptron(SLFN):

    def __init__(self,
                 d = None,
                 h = None,
                 biasedOutput = False,
                 mapping = 0,
                 swap = False,
                 swapFile = './slfn.H.swap'):
        super(Perceptron, self).__init__(d,
                                         d,
                                         biasedOutput = biasedOutput,
                                         mapping = mapping,
                                         swap = swap,
                                         swapFile = swapFile)

    def activation(self, x):
        assert isinstance(self.d, int)
        assert isinstance(self.h, int)
        assert (self.d > 0)
        assert (self.h > 0)
        assert (self.h == self.d)
        row = numpy.array(x)
        return row

class SigmoidNetwork(SLFN):

    def __init__(self,
                 d = None,
                 h = None,
                 biasedOutput = False,
                 mapping = 0,
                 swap = False,
                 swapFile = './slfn.H.swap'):
        super(SigmoidNetwork, self).__init__(d,
                                             h,
                                             biasedOutput = biasedOutput,
                                             mapping = mapping,
                                             swap = swap,
                                             swapFile = swapFile)

    def activation(self, x):
        assert isinstance(self.d, int)
        assert isinstance(self.h, int)
        assert (self.d > 0)
        assert (self.h > 0)
        assert (self.w.shape[0] == self.d + 1 and self.w.shape[1] == self.h)
        assert (x.shape[1] + 1 == self.w.shape[0])
        try:
            L = numpy.shape(x)[0]
        except IndexError:
            L = 1
        bias = numpy.matrix(numpy.ones(L)).T
        row = numpy.array(numpy.concatenate((x, bias), axis = 1) * self.w)
        # row : L x h
        row = 0.5 + numpy.tanh(row / 2.) / 2.
        return row

class RBFL1Network(SLFN):
    """
    A class, derived from SLFN, for Radial Basis Function Networks
    with L1-norm.
    """
    def __init__(self,
                 d = None,
                 h = None,
                 biasedOutput = False,
                 mapping = 0,
                 swap = False,
                 swapFile = './slfn.H.swap'):
        """
        Construct a RBFL1Network object.
        """
        super(RBFL1Network, self).__init__(d,
                                           h,
                                           wrange = (0., 0.5),
                                           biasrange = (0., 1.),
                                           biasedOutput = biasedOutput,
                                           mapping = mapping,
                                           swap = swap,
                                           swapFile = swapFile)

    def activation(self, x):
        """
        RBF activation function exp(-|x - c| / sigma)
        """
        row = []
        try:
            L = numpy.shape(x)[0]
        except IndexError:
            L = 1
        for j in range(L):
            col = x[j] - self.w[:-1].T
            col = numpy.abs(col)
            col = numpy.sum(col, axis = 1)
            col = numpy.exp(-col / self.w[-1].T)
            row.append(numpy.array(numpy.transpose(col))[0])
        row = numpy.array(row)
        return row

class RBFL2Network(SLFN):
    """
    A class, derived from SLFN, for Radial Basis Function Networks
    with L2-norm.
    """
    def __init__(self,
                 d = None,
                 h = None,
                 biasedOutput = False,
                 mapping = 0,
                 swap = False,
                 swapFile = './slfn.H.swap'):
        """
        Construct a RBFL2Network object.
        """
        super(RBFL2Network, self).__init__(d,
                                           h,
                                           wrange = (0., 0.5),
                                           biasrange = (0., 1.),
                                           biasedOutput = biasedOutput,
                                           mapping = mapping,
                                           swap = swap,
                                           swapFile = swapFile)

    def activation(self, x):
        """
        RBF activation function exp(-||x - c||^2 / sigma^2)
        """
        row = []
        try:
            L = numpy.shape(x)[0]
        except IndexError:
            L = 1
        for j in range(L):
            # computing x[j] - c[0:h] : h x d
            col = x[j] - self.w[:-1].T
            # computing r_k = || x[j] - c[k] ||^2 for each k=0,...,h
            col = numpy.power(col, 2)
            col = numpy.sum(col, axis = 1)
            # computing r_k / sigma^2 for k=0,...,h
            col = numpy.exp(-col / numpy.power(self.w[-1].T, 2))
            # appending the jth row
            row.append(numpy.array(numpy.transpose(col))[0])
        # row : L x h
        row = numpy.array(row)
        return row
