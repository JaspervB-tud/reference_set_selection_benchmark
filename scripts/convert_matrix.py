import numpy as np
import os
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
    
    os.makedirs(args.output, exist_ok=True)
    with open(f"{args.output}/converted_matrix.mat", "w") as f_out: #write output
        f_out.write("\t".join(names) + "\n")
        for i in range(len(names)):
            cur_numbers = [str(d) for d in distance_matrix[i,:]]
            cur_line = "{}\t{}".format(names[i], "\t".join(cur_numbers))
            f_out.write(cur_line + "\n")

if __name__ == "__main__":
    main()