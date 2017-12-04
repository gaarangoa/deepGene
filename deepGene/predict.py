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
import json

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
        rd = Kmers.extract_features(fi=self.fi, kf=self.fk)
        X = []
        L = []
        for i in rd:
            X.append(i['features'])
            L.append(i['gene_name'])
        # print(json.dumps(self.X, indent=10))
        V = DictVectorizer(sparse=False)
        self.X = V.fit_transform( self.X )
        self.FT = np.array( [str(i) for i in V.get_feature_names()] )
        print(self.X, self.FT)

def main(args):
    predictor = Predict(args);
    predictor.pred()
