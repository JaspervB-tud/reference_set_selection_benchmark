import numpy as np
import os
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster
import argparse

def read_matrix(matrix_path):
    names = []
    with open(matrix_path, "rt") as f_in:
        header = next(f_in).strip()
        num_sequences = int(header) #header stores the number of sequences
        distance_matrix = np.zeros((num_sequences, num_sequences), dtype=np.float64)
        i = 0
        for line in f_in:
            line = line.strip().split("\t")
            names.append(line[0].split("/")[-1].replace(".fna.gz", ""))
            for j in range(1, len(line)):
                distance_matrix[i,j-1] = float(line[j])
                distance_matrix[j-1,i] = float(line[j])
            i += 1
                
    return (distance_matrix, names)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=str, required=True)
    parser.add_argument("--threshold", type=float, required=True) #note that this is a DISTANCE threshold (1 - similarity threshold)
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()
    
    distance_matrix, names = read_matrix(args.matrix)
    nonzeros = np.array([distance_matrix[i,j] for i in range(distance_matrix.shape[0]) for j in range(i) if distance_matrix[i,j] != 0])
    for linkage_type in ["single", "complete"]:
        if len(nonzeros) > 0: #if there are nonzeros, perform clustering, otherwise output first genome
            Z = linkage(squareform(distance_matrix), method=linkage_type)    
            clusters = fcluster(Z, t=args.threshold, criterion="distance")
            representatives = {}
            for cluster in np.unique(clusters):
                indices = np.where(clusters == cluster)[0]
                if len(indices) == 1:
                    representatives[cluster] = names[indices[0]]
                else:
                    intra_distances = distance_matrix[np.ix_(indices, indices)]
                    centroid = np.argmin(np.sum(intra_distances, axis=1))
                    representatives[cluster] = names[indices[centroid]]
        else:
            representatives = {0: names[0]}
        
        os.makedirs(args.output, exist_ok=True)
        with open(f"{args.output}/{linkage_type}-linkage_{args.threshold}", "w") as f_out: #write output
            for cluster in representatives:
                f_out.write(representatives[cluster] + "\n")
    
if __name__ == "__main__":
    main()