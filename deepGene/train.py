from keras.models import Sequential, load_model
from keras.layers import Dense, Activation, Dropout
import h5py
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
from keras.utils import np_utils

from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectFromModel

import Kmers

class DLModel():
    def __init__(self, args):
        self.args = args
        self.dataset = args.dataset
        self.modelName = args.model
        try:
            self.epochs = args.epochs
        except:
            self.epochs = 100
        try:
            self.testFraction = args.test
        except:
            self.testFraction = 0.3
        try:
            self.batch_size = args.batch_size
        except:
            self.batch_size = 64
        # self.validation = args.validation
        # self.fullModel = args.fullmodel
    
    def loadDataset(self, r=False):
        self.f = h5py.File(self.dataset, 'r')
        self.min_max_scaler = preprocessing.MinMaxScaler()
        self.X = self.min_max_scaler.fit_transform( np.array( [ i for i in self.f['dataset/X'] ] ) )
        self.F = [ i for i in self.f['dataset/F'] ]
        self.Y = [ i for i in self.f['dataset/Y'] ]
        
        print("Extract features using random forest")
        clf = ExtraTreesClassifier()
        clf = clf.fit(self.X, self.Y)
        rfModel = SelectFromModel(clf, prefit=True)
        self.X = rfModel.transform(self.X)

        # self.selected_features = np.asarray(clf.get_feature_names())[rfModel.get_support()]

        print("build the Y labels (numerical values for regression or categorical for classification")

        if r:
            self.Y = self.min_max_scaler.fit_transform( np.array( [float(i) for i in self.f['dataset/Y']] ).reshape(-1, 1) )
        else:
            self.YEncoder = preprocessing.LabelEncoder()
            self.YEncoder.fit( self.Y )
            self.EncodedY = self.YEncoder.transform( [ i for i in self.f['dataset/Y'] ] )
            self.Y = np_utils.to_categorical( self.EncodedY )

    def loadModel(self):
        self.model = load_model(self.modelName)

    def createModel(self):
        print("number of features: ", len(self.F))
        print("number of selected features: ", self.X.shape[1])
        print("Input dataset: ", self.X.shape)
        self.nLabels = self.YEncoder.classes_
        self.model = Sequential()
        self.model.add( Dense(units=500, input_dim=self.X.shape[1]) )
        self.model.add( Activation('relu') )
        self.model.add( Dropout(0.5) )
        self.model.add( Dense(units=250) )
        self.model.add( Activation('relu') )
        self.model.add( Dropout(0.5) )
        self.model.add( Dense(units=100) )
        self.model.add( Activation('relu') )
        self.model.add( Dropout(0.5) )
        self.model.add( Dense(units=50) )
        self.model.add( Activation('relu') )
        self.model.add( Dense(units=len(self.nLabels)) )
        self.model.add( Activation('softmax') )

        print( self.model.summary() )

        self.model.compile(
            loss = "categorical_crossentropy",
            optimizer = "adam",
            metrics = ["accuracy"]
        )
    
    def createRegressionModel(self):
        print("number of features: ", len(self.F))
        self.nLabels = 1
        self.model = Sequential()
        self.model.add( Dense(units=1000, input_dim=self.X.shape[1]) )
        self.model.add( Activation('relu') )
        self.model.add( Dropout(0.5) )
        self.model.add( Dense(units=500) )
        self.model.add( Activation('relu') )
        self.model.add( Dropout(0.5) )
        self.model.add( Dense(units=100) )
        self.model.add( Activation('relu') )
        self.model.add( Dropout(0.5) )
        self.model.add( Dense(units=50) )
        self.model.add( Activation('relu') )
        self.model.add( Dense(units=1) )
        self.model.add( Activation('linear') )

        print( self.model.summary() )

        self.model.compile(
            loss = "mean_squared_error",
            optimizer = "rmsprop",
            metrics = ["accuracy"]
        )
        return self.model

    def trainModel(self):
        self.x_train, self.x_test, self.y_train, self.y_test = train_test_split(
            self.X, 
            self.Y, 
            test_size = self.testFraction
        )
        self.model.fit(
            self.x_train,
            self.y_train,
            epochs = self.epochs,
            batch_size = self.batch_size,
            verbose = 1,
            validation_split=0.3
        )

    def testModel(self):
        score = self.model.evaluate(self.x_test, self.y_test, batch_size=self.batch_size)
        print(score)
    
    def saveModel(self):
        self.model.save(self.modelName)
        # f = h5py.File(self.modelName+".selected_features.hdf5", "w")
        # dataset = f.create_group('selected_features')
        # dataset.create_dataset(self.selected_features, 'F')
        # f.close()
    
    def modelWeights(self):
        w0 = [ max(i) for i in self.model.layers[0].get_weights()[0] ]
        b0 = [ max(i) for i in self.model.layers[0].get_weights()[0] ]
        fo = open(self.modelName+'.wg','w')
        for ix,i in enumerate(w0):
            fo.write( "\t".join([str(i),str(self.F[ix])])+"\n" )
      

def main(args):
    print("create the object model")
    ML = DLModel(args)

    print("load the dataset")
    ML.loadDataset()

    print("create the deepLearning model")
    ML.createModel()

    print("train the deep learning model")
    ML.trainModel()

    print("test deep L model")
    ML.testModel()

    print("save model")
    ML.saveModel()

    print("save weights")
    ML.modelWeights()

def weights(args):
    ML = DLModel(args)
    
    print('loading dataset ...')
    ML.loadDataset()

    print('loading model ...')
    ML.loadModel()

    print('retrieving weights ...')
    ML.modelWeights()

def regression(args):
    ML = DLModel(args)

    # load the dataset
    ML.loadDataset(r=True)

    # create the deepLearning model
    ML.createRegressionModel()

    # train the deep learning model
    ML.trainModel()

    # test deep L model
    ML.testModel()

    # save model
    ML.saveModel()