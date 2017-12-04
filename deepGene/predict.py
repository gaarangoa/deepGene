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
        self.X = []
        self.L = []
        for i in rd:
            self.X.append(i['features'])
            self.L.append(i['gene_id'])
        # print(json.dumps(self.X, indent=10))
        V = DictVectorizer(sparse=False)
        self.X = V.fit_transform( self.X )[:-1] # I need to remove the last element, because it corersponds to the synthetic gene 
        self.L = self.L[:-1] # remove the last position of the array that corresponds to the synthetic gene
        
        pred = self.model.predict(self.X)
        print("query\t%Low\t%High")
        for ix,i in enumerate(pred):
            print("\t".join([
                    self.L[ix],
                    str(100*i[0]),
                    str(100*i[1])
                ]))

def main(args):
    predictor = Predict(args);
    predictor.pred()
