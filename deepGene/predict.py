from keras.models import Sequential, load_model
from keras.layers import Dense, Activation, Dropout
import h5py
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
from keras.utils import np_utils

from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.feature_extraction import DictVectorizer

import Kmers

class Predict():
    def __init__(self, args):
        ''' '''
        self.args = args
        self.model_name = args.model
        self.model = load_model(args.model)
        self.fk = args.kmers
        self.fi = args.input

    def pred(self):
        print("Loading input file and feature extraction")
        self.X = Kmers.extract_features(fi=self.fi, kf=self.fk)
        print(json.dumps(self.X, indent=10))
        # V = DictVectorizer(sparse=False)
        # self.X = V.fit_transform( self.X.values() )
        # self.N = X.keys()
        # self.FT = np.array( [str(i) for i in V.get_feature_names()] )

def main(args):
    predictor = Predict(args);
    predictor.pred()
