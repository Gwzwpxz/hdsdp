from scipy.io import loadmat
from scipy.sparse import csc_matrix
import argparse

parser = argparse.ArgumentParser(description='Read .mat matrix in CSC format and convert to .npz file.')
parser.add_argument('--input', type=str, default='/Users/gaowenzhi/Desktop/gwz/hdsdp/gpu/tmp.mat')

if __name__ == '__main__':
    
    args = parser.parse_args()
    mat = loadmat(args.input)
    mat = csc_matrix(mat['A'])

    Abeg = mat.indptr
    Aidx = mat.indices
    Adata = mat.data

    # Write to csv Ap.csv, Ai.csv, Ax.csv
    with open('Ap.csv', 'w') as f:
        f.write('\n'.join(str(i) for i in Abeg))
    with open('Ai.csv', 'w') as f:
        f.write('\n'.join(str(i) for i in Aidx))
    with open('Ax.csv', 'w') as f:
        f.write('\n'.join(str(i) for i in Adata))
    with open('dim.csv', 'w') as f:
        f.write(str(mat.shape[0]))
