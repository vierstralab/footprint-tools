"""
"""

import numpy as np
import itertools
import random

class bias_model(object):

    def __init__(self):

        self.model = {}
        
        self.k = 6
        self.mid = 3

    def __getitem__(self, key):
        
        return self.model.get(key, 1e-6)

    def __setitem__(self, key, value):

        self.model[key] = value

    def offset(self):

        return max(self.k-self.mid, self.mid)

    def shuffle(self):
        """Randomly shuffle the bias model

        Returns
        -------
        model : :class:`bias_model`
            A shuffled bias model
        """
        ret = bias_model()
        ret.model = { x: y for (x, y) in zip(list(self.model.keys()), sorted(list(self.model.values()), key = lambda k: random.random())) }
        ret.offset = self.offset
        return ret

    def predict(self, probs, n = 100):
        """Compute cleavage propensities from sequence
        
        Parameters
        ----------
        probs : :class:`numpy.ndarray`
            An array of probilities (relative values)
        n : int
            Number of total tags to distrbute
        """

        return np.around( probs / np.sum(probs) * n )

class kmer_model(bias_model):

    def __init__(self, filepath):
        
        bias_model.__init__(self)		
        self.read_model(filepath)
    
    def read_model(self, filepath):
        """Read the k-mer model from a file.

        Parameters
        ----------
        filepath : str
            Path to a K-mer model file

        """

        import urllib.request as request

        try:

            if filepath.startswith('http'):
                file = request.urlopen(filepath)
            else:
                file = open(filepath, 'r')

            for line in file:
                (seq, prob) = line.strip().split('\t')
                self.model[seq.upper()] = float(prob)

        except IOError:
            
            raise IOError("Cannot open file: %s" % filepath)
    
    def probs(self, seq):
        """Generate cleavage preference array from DNA sequence

        Parameters
        ----------
        seq: str
            A DNA sequence to compute relative sequence preference

        Returns
        -------
        out : :class:`numpy.ndarray`
            Relate sequence preferencew
        """
        k = self.k
        mid = self.mid
        offset = self.offset()

        return np.array([ self.__getitem__( seq[(i-mid):(i+(k-mid))] ) for i in range(offset, len(seq)-offset) ], dtype = np.float64)

class uniform_model(bias_model):

    def __init__(self):

        bias_model.__init__(self)

        for seq in itertools.product('ATCG', repeat = self.k):
            self.model[''.join(seq)] = 1.0

    def probs(self, seq):
        return np.ones( len(seq) )
