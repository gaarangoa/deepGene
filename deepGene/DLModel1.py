# deep learning model using LSTMs architecture
# Input are sequences in fasta format and labels in a separate file
# Output the model architecture

__dev__='gustavo arango'


import json
from keras.preprocessing.text import one_hot
from keras.preprocessing.sequence import pad_sequences
from keras.models import Sequential
from keras.layers import Activation, Dense, Flatten, Dropout, LSTM, concatenate, Input, Reshape
from keras.layers.embeddings import Embedding
from keras.utils import np_utils
from keras.models import Model

from sklearn import preprocessing
import numpy as np
import sys
from Bio import SeqIO

from functools import partial
from StringIO import StringIO

from tqdm import tqdm as bar

class Train():
    def __init__(self, args):
        self.info = ''
        self.input_genes_fasta = args.input
        self.input_genes_labels = args.labels
        self.output_model = args.output
        self.classifier_params = {
            "vocab_size": 500,
            "max_length": 5000,
            "embedding_size": 500,
            "output": 2
        }
    
    def run(self):
        # make the dataset
        labels_raw = {i.strip().split()[0]:i.strip().split()[1] for i in open(self.input_genes_labels)}
        docs = []
        labels = []
        data = SeqIO.parse(self.input_genes_fasta, "fasta")
        for record in bar(data):
            try:
                labels.append( float(labels_raw[record.id]) )
                docs.append(
                    " ".join([l for l in iter(partial(StringIO(str(record.seq)).read, 3), '')])
                )
            except Exception as e:
                pass
        
        labels = np.array(labels)
        
        # Encode documents
        encoded_docs = [one_hot(d, self.classifier_params['vocab_size']) for d in docs] #uses a hash function to represent words
        padded_docs = pad_sequences(encoded_docs, maxlen=self.classifier_params['max_length'], padding='post')

        # Build model
        text_model_input = Input(shape = (self.classifier_params['max_length'],), dtype="int32", name = 'text_model_input')
        text_model = Embedding(input_dim = self.classifier_params['vocab_size'], mask_zero=True, output_dim = self.classifier_params['embedding_size'], input_length = self.classifier_params['max_length'], name="embedding" )(text_model_input)
        text_model = LSTM(512, name = "text-lstm-1", return_sequences=True)(text_model)
        text_model_output = LSTM(256, name = 'text-lstm-2')(text_model)

        # Merge model
        merged_model = concatenate([text_model_output, text_model_output], axis=1)
        merged_model = Dense(1200, activation="relu")(merged_model)
        merged_model = Dropout(0.5)(merged_model)
        merged_model = Dense(640, activation="relu")(merged_model)
        merged_model = Dropout(0.5)(merged_model)
        merged_model = Dense(420, activation="relu")(merged_model)
        merged_model_output = Dense(1, name = 'merged_model_output')(merged_model)

        model = Model(inputs = [text_model_input], outputs = [merged_model_output])
        model.compile(optimizer='adam', loss='mean_squared_error', metrics=['acc'])
        print(model.summary())

        try:
            from keras.utils import plot_model
            plot_model(model, to_file='model.png')
        except:
            pass

        # Train the model
        model.fit([padded_docs], [labels], batch_size=128, epochs=100)

        # save model
        model.save(self.output_model+'.hdf5')


class test():
    def __init__(self):
        self.input = '../data/gene_sequence.fasta'
        self.output = '../data/model'
        self.labels = '../data/regressionLabels'

args = test()
train = Train(args)
train.run()


