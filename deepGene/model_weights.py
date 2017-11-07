import train as tr

def main(args):
    ML = tr.DLModel(args)
    
    # Load model
    ML.loadModel()
    
    # Load dataset
    ML.loadDataset()

    # 
