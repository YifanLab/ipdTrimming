import numpy as np
import sys
import pandas as pd

def cluster_array(input_array, distance_threshold):
    if len(input_array) == 0:
        return []
    input_array = np.array(input_array)
    clusters = []
    current_cluster = [input_array[0]]

    for i in range(1, len(input_array)):
        if input_array[i] - input_array[i - 1] < distance_threshold:
            current_cluster.append(input_array[i])
        else:
            if len(current_cluster) >= 3:
                clusters.append((current_cluster[0], current_cluster[-1]))
            current_cluster = [input_array[i]]

    if len(current_cluster) >= 3:
        clusters.append((current_cluster[0], current_cluster[-1]))

    return clusters

def calculate_adjacent_distances(clusters):
    distances = []
    size = [int(clusters[i][1] - clusters[i][0]+1) for  i in range(len(clusters))]
    for i in range(1, len(clusters)):
        previous_cluster_end = clusters[i - 1][1]
        current_cluster_start = clusters[i][0]
        distance = current_cluster_start - previous_cluster_end
        distances.append(distance)
    return distances,size

# Example usage
def main(ipdcsv, distance_threshold):
    ipdout = sys.argv[3] + '_clusterallfmt.xls'
    sm6mA_distri = np.zeros(600, dtype = int)
    smsize_distri = np.zeros(600, dtype = int)
    smadj_distri = np.zeros(600, dtype = int)
    bin_edges = np.linspace(0, 600, 601)
    with open(ipdcsv, 'r') as f, open(ipdout, 'w') as fout:
        for line in f:
            row = line.strip().split('\t')
            if len(row) < 9:  # Ensure there are enough columns
                print(line.strip(), "Missing columns", sep="\t")
                continue

            sm6mA = row[7]
            ref6mA = row[8]
            
            try:
                sm6mA_array = np.fromstring(sm6mA[1:-1], sep=', ')
                ref6mA_array = np.fromstring(ref6mA[1:-1], sep=', ')
            except ValueError:
                print(line.strip(), "Invalid array format", sep="\t")
                continue

            if sm6mA_array.size == 0 or ref6mA_array.size == 0:
                print(line.strip(), "Empty array", sep="\t")
                continue

            sm6mA_arraydf = np.diff(sm6mA_array)
            sm6mA_arraydf = [int(x) for x in sm6mA_arraydf]
            output_arr = cluster_array(sm6mA_array, distance_threshold)
            adjdis_arr, sizes = calculate_adjacent_distances(output_arr)
            smadj_distri += np.histogram(adjdis_arr, bins = bin_edges)[0]
            smsize_distri += np.histogram(sizes, bins = bin_edges)[0]
            sm6mA_distri += np.histogram(sm6mA_arraydf, bins = bin_edges)[0]
            fout.write(f"{line.strip()}\t{output_arr}\t{adjdis_arr}\t{sizes}\t{sm6mA_arraydf}\n")
    smcluadjcfmt = pd.DataFrame({'Distance': bin_edges[:-1], 'Count': smadj_distri})
    smcluadjcfmt.to_csv(sys.argv[3]+'_cluadjdistance.xls', header = False, sep = "\t", index = False)
    smadjcfmt = pd.DataFrame({'Distance': bin_edges[:-1], 'Count': sm6mA_distri})
    smadjcfmt.to_csv(sys.argv[3]+'_6mAadjdistance.xls', header = False, sep = "\t", index = False)
    smsizecfmt = pd.DataFrame({'Distance': bin_edges[:-1], 'Count': smsize_distri})
    smsizecfmt.to_csv(sys.argv[3]+'_smsizesdistance.xls', header = False, sep = "\t", index = False)
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python array_distance.py 6mAvs5mC.csv cluster_threshold")
        sys.exit(1)
    main(sys.argv[1], int(sys.argv[2]))
