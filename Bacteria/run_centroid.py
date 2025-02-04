import numpy as np

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
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()
    
    distance_matrix, names = read_matrix(args.matrix)
    centroid = np.argmin(np.sum(distance_matrix, axis=1)) #select centroid
    
    with open(args.output, "w") as f_out: #write output
        f_out.write(names[centroid] + "\n")
    
if __name__ == "__main__":
    main()