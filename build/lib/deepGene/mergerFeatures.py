from sklearn.feature_extraction import DictVectorizer
import json
import h5py
import numpy as np
import sys

def main(args):
    genes = {}
    # hash gene_ids
    selectedGenes = {i.split()[0]:True for i in open(args.genes)}
    
    # hash the gene labels
    geneLabels = {i.split()[0]:i.split()[1] for i in open(args.genes)}

    for file in args.features:
        data = json.load(open(file))
        for ix,gene in enumerate(data):
            geneId = gene['gene_id']
            try:
                assert(selectedGenes[geneId])
                for feature in gene['features']:
                    try:
                        genes[geneId].update({feature['feature_name']:feature['value']})
                    except:
                        genes[geneId] = {feature['feature_name']:feature['value']}
            except Exception as inst:
                # print(inst)
                pass
    
    # print genes

    V = DictVectorizer(sparse=False)
    X = V.fit_transform( genes.values() )
    Y = [geneLabels[i] for i in genes]
    FT = np.array( [str(i) for i in V.get_feature_names()] )
    
    f = h5py.File(args.output, 'w')
    
    dataset = f.create_group('dataset')
    x = dataset.create_dataset('X', data=X)
    y = dataset.create_dataset('Y', data = Y)
    z = dataset.create_dataset('F', data = FT)
    f.close()

    # X = f['dataset/X']
    # Y = f['dataset/Y']