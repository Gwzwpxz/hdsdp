import numpy as np
from scipy.io import loadmat

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--f", type=str, default="pot-test.mat")

def main():
    
    args = parser.parse_args()
    file = args.f
        
    data = loadmat(file)
    A = data['data'][0][0][0]
    b = data['data'][0][0][1].flatten()
    c = data['data'][0][0][2].flatten()
    
    Abeg = ", ".join([str(x) for x in A.indptr])
    Aidx = ", ".join([str(x) for x in A.indices])
    Adata = ", ".join([str(x) for x in A.data])
    bstr = ", ".join([str(x) for x in b])
    cstr = ", ".join([str(x) for x in c])
    
    
    with open("data.h", "w") as f:
    
        f.write("int nCol = {0}; \n".format(len(c)));
        f.write("int nRow = {0}; \n ".format(len(b)));
        f.write("int Ap[] = {{{0}}}; \n".format(Abeg))
        f.write("int Ai[] = {{{0}}}; \n".format(Aidx))
        f.write("double Ax[] = {{{0}}}; \n".format(Adata))
        f.write("double rhs[] = {{{0}}}; \n".format(bstr))
        f.write("double obj[] = {{{0}}}; \n".format(cstr))
    
    f.close()
    
    print("Generated {0}".format("{0}.h".format(file[file.find('./') + 1:-4])))
    
    
if __name__ == '__main__':
    
    main()