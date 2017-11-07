from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
import h5py
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy as np
from keras.utils import np_utils

class DLModel():
    def __init__(self, args):
        self.dataset = args.dataset
        self.saveDir = args.model
        self.epochs = args.epochs
        self.testFraction = args.test
        self.batch_size = args.batch_size
        # self.validation = args.validation
        # self.fullModel = args.fullmodel
    
    def loadDataset(self):
        self.f = h5py.File(self.dataset, 'r')
        self.min_max_scaler = preprocessing.MinMaxScaler()
        self.YEncoder = preprocessing.LabelEncoder()
        self.YEncoder.fit( [ i for i in self.f['dataset/Y'] ] )
        self.X = self.min_max_scaler.fit_transform( np.array( [ i for i in self.f['dataset/X'] ] ) )
        self.EncodedY = self.YEncoder.transform( [ i for i in self.f['dataset/Y'] ] )
        self.Y = np_utils.to_categorical( self.EncodedY )
        self.F = [ i for i in self.f['dataset/F'] ]

    def createModel(self):
        self.nLabels = self.YEncoder.classes_
        self.model = Sequential()
        self.model.add( Dense(units=1000, input_dim=len(self.F)) )
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
        self.model.add( Dense(units=len(self.nLabels)) )
        self.model.add( Activation('softmax') )

        print( self.model.summary() )

        self.model.compile(
            loss = "categorical_crossentropy",
            optimizer = "adam",
            metrics = ["accuracy"]
        )
    
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
            validation_split=0.1
        )

    def testModel(self):
        score = self.model.evaluate(self.x_test, self.y_test, batch_size=self.batch_size)
        print(score)
    
    def saveModel(self):
        self.model.save(self.saveDir+"model.hdf5")
    
def main(args):
    # create the object model
    ML = DLModel(args)

    # load the dataset
    ML.loadDataset()

    # create the deepLearning model
    ML.createModel()

    # train the deep learning model
    ML.trainModel()

    # test deep L model
    ML.testModel()

    # save model
    ML.saveModel()