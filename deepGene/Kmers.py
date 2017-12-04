from Bio import SeqIO
import json

def match_kmers(sequence='', kmers={}, ks = []):
    ''' Given a list of kmers, see if those kmers are within the sequence '''
    cmpx = {}
    for i in range(len(sequence)):
        for k in ks:
            try:
                kmer = sequence[i:i+k]
                assert(kmers[kmer])
                cmpx[kmer] += 1
            except Exception as e:
                try:
                    assert(kmers[kmer])
                    cmpx[sequence[i:i+k]] = 1
                except:
                    pass
    fcmpx = []
    for i in cmpx:
        fcmpx.append(
            {
                "feature_name": i, 
                "value": cmpx[i]
            }
        )
    return [fcmpx, cmpx]


def get_kmer_list(fi=''):
    ''' from the input file [weigth, kmer], subtract the kmer and its length ''' 
    kmers = {i.strip().split()[1]:len(i.strip().split()[1]) for i in open(fi)}
    ks = list(set(kmers.values()))
    return [kmers,ks]


def process_input_fasta(fi='', kmers={}, ks=[]):
    data = []
    for record in SeqIO.parse(fi, "fasta"):
        gene_name = record.id 
        sequence = str(record.seq)
        ffeatures, features = match_kmers(sequence=sequence, kmers = kmers, ks = ks)
        data.append({
            "gene_id": gene_name,
            "features": features
        })
    
    # add a synthetic sample to get all the features, other case the matrix will not be complete

    sgene = {i:0 for i in kmers}
    sgname = 'synthetic_gene'
    data.append([{"gene_id": sgname, "features": sgene}])
    # last position of the data is a synthetic gene, the purpose is to get all the features
    return data

def extract_features(fi='', kf=''):
    ''' Extract features from a list of kmer ids
            fi: fasta file as input.
            kf: file that contains the kmers. This file is provided by the model after training
    '''
    print("Extract features: ")
    kmers, ks = get_kmer_list(fi=kf)
    features = process_input_fasta(fi=fi, kmers=kmers, ks = ks)
    # json.dump(features, open(fi+'.ft', 'w'))
    return features

    



# sequence = 'abcdefghijklmnopqrstabcde'
# kmers = {'abcdef': True, 'ab':True, 'abc':True}

# ks = [2,3,6]

# traverse(sequence=sequence, kmers=kmers, ks=ks)