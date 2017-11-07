import sys
import argparse
import json

def main(fi, fo):
    data = {}
    for i in open(fi):
        if i[0] == "#": continue
        i = i.split()
        try:
            if data[ (i[0],i[2]) ] < float(i[6]):
                data[ (i[0],i[2]) ] = float(i[6])
        except:
            data[ (i[0],i[2]) ] = float(i[6])

    genes = {}
    for i in data:
        item = {
            "feature_name": i[0],
            "value": data[i]
        } 
        try:
            genes[i[1]].append(item)
        except:
            genes[i[1]] = [item]
    
    del data
    gns = []
    for i in genes:
        gns.append({
            "gene_id": i, 
            "features": genes[i]
        })

    json.dump(gns, open(fo, "w"))