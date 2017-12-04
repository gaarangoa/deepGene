

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
    fcmpx = {}
    for i in cmpx:
        fcmpx.update(
            {
                "feature_name": i, 
                "value": cmpx[i]
            }
        )
    return fcmpx


def get_kmer_list(fi):
    ''' from the input file [weigth, kmer], subtract the kmer and its length ''' 
    kmers = {i.strip().split()[1]:len(i.strip().split()[1]) for i in open(fi)}
    ks = list(set(kmers.values()))
    return [kmers,ks]


def process_input_fasta(fi):
    


# sequence = 'abcdefghijklmnopqrstabcde'
# kmers = {'abcdef': True, 'ab':True, 'abc':True}

# ks = [2,3,6]

# traverse(sequence=sequence, kmers=kmers, ks=ks)