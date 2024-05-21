import pysam
import numpy as np
import multiprocessing
import sys
import os
import tempfile
import shutil

os.environ['TMPDIR'] = os.getcwd()

def calculate_statis(data):
    data_array = np.array(data)
    mean = np.mean(data_array)
    std = np.std(data_array)
    variance = np.std(data_array) / np.mean(data_array) if np.mean(data_array) != 0 else float('inf')
    return mean, std, variance

def calbase_statis(arr1, arr2, arr3):
    (arr1_m, arr1_std, arr1_var) = calculate_statis(arr1)
    (arr2_m, arr2_std, arr2_var) = calculate_statis(arr2)
    (arr3_m, arr3_std, arr3_var) = calculate_statis(arr3)
    return arr1_m, arr1_std, arr1_var, arr2_m, arr2_std, arr2_var, arr3_m, arr3_std, arr3_var

def reverse_complement(nucleotide):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement[base] for base in nucleotide[::-1])

def read_fetcher(bam_path, read_queue, chunk_size):
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam_file:
        buffer = []
        for read in bam_file.fetch(until_eof=True):
            buffer.append(read.to_dict())  # Convert read to a dictionary to ensure it's pickle-able
            if len(buffer) >= chunk_size:
                read_queue.put(buffer)
                buffer = []
        if buffer:
            read_queue.put(buffer)  # Put the last chunk in the queue
    read_queue.put(None)  # Signal that reading is complete

def process_reads(reads, result_file):
    keys = ['aipd', 'tipd', 'cipd', 'gipd', 'afipd', 'tfipd', 'cfipd', 'dfipd', 'aipdr', 'tipdr', 'cipdr', 'gipdr', 'gfipd', 'afipdr', 'tfipdr', 'cfipdr', 'gfipdr', 'nonafipdr']
    strand = [0, 1, 2]
    list_of_dict = [{key: [] for key in keys} for _ in strand]
    
    with open(result_file, 'a') as f:
        for read_dict in reads:
            read = pysam.AlignedSegment.from_dict(read_dict, None)  # Convert dictionary back to read
            seq = read.query_sequence
            try:
                ec = read.get_tag('ec')
                fcontrol =  read.get_tag('FC')
                fnormal =  read.get_tag('FZ')
                fframenorm =  read.get_tag('FF')
                fnipdr =  read.get_tag('FN')
                fripdr =  read.get_tag('FR')
                icontrol =  read.get_tag('IC')[::-1]
                inormal =  read.get_tag('IZ')[::-1]
                iframenorm  = read.get_tag('IF')[::-1]
                inipdr = read.get_tag('IN')[::-1]
                iripdr = read.get_tag('IR')[::-1]
            except KeyError as e:
                print(f"Missing tag in read {read.query_name}:{e}")
                continue
            
            if len(fcontrol) == len(fnormal) == len(fframenorm) == len(fnipdr) == len(fripdr) == len(seq):
                fnormal_m, fnormal_std, fnormal_var = calculate_statis(fnormal)
                fframenorm_m, fframenorm_std, fframenorm_var = calculate_statis(fframenorm)
                fnipdr_m, fnipdr_std, fnipdr_var = calculate_statis(fnipdr)
                fripdr_m, fripdr_std, fripdr_var = calculate_statis(fripdr)
                inormal_m, inormal_std, inormal_var = calculate_statis(inormal)
                iframenorm_m, iframenorm_std, iframenorm_var = calculate_statis(iframenorm)
                inipdr_m, inipdr_std, inipdr_var = calculate_statis(inipdr)
                iripdr_m, iripdr_std, iripdr_var = calculate_statis(iripdr)
                tnormal_m, tnormal_std, tnormal_var = calculate_statis(np.concatenate((fnormal, inormal), axis=None))
                tframenorm_m, tframenorm_std, tframenorm_var = calculate_statis(np.concatenate((fframenorm, iframenorm), axis=None))
                tnipdr_m, tnipdr_std, tnipdr_var = calculate_statis(np.concatenate((fnipdr, inipdr), axis=None))
                tripdr_m, tripdr_std, tripdr_var = calculate_statis(np.concatenate((fripdr, iripdr), axis=None))
                for i in range(len(seq)):
                    if seq[i] == 'A':
                        list_of_dict[0]['aipd'].append(inormal[i])
                        list_of_dict[1]['tipd'].append(fnormal[i])
                        list_of_dict[0]['afipd'].append(iframenorm[i])
                        list_of_dict[1]['tfipd'].append(fframenorm[i])
                        list_of_dict[0]['aipdr'].append(inipdr[i])
                        list_of_dict[1]['tipdr'].append(fnipdr[i])
                        list_of_dict[0]['afipdr'].append(iripdr[i])
                        list_of_dict[1]['tfipdr'].append(fripdr[i])
                    if seq[i] == 'T':
                        list_of_dict[0]['tipd'].append(inormal[i])
                        list_of_dict[1]['aipd'].append(fnormal[i])
                        list_of_dict[0]['tfipd'].append(iframenorm[i])
                        list_of_dict[1]['afipd'].append(fframenorm[i])
                        list_of_dict[0]['tipdr'].append(inipdr[i])
                        list_of_dict[1]['aipdr'].append(fnipdr[i])
                        list_of_dict[0]['tfipdr'].append(iripdr[i])
                        list_of_dict[1]['afipdr'].append(fripdr[i])
                    if seq[i] == 'C':
                        list_of_dict[0]['cipd'].append(inormal[i])
                        list_of_dict[1]['gipd'].append(fnormal[i])
                        list_of_dict[0]['cfipd'].append(iframenorm[i])
                        list_of_dict[1]['gfipd'].append(fframenorm[i])
                        list_of_dict[0]['cipdr'].append(inipdr[i])
                        list_of_dict[1]['gipdr'].append(fnipdr[i])
                        list_of_dict[0]['cfipdr'].append(iripdr[i])
                        list_of_dict[1]['gfipdr'].append(fripdr[i])
                    if seq[i] == 'G':
                        list_of_dict[0]['gipd'].append(inormal[i])
                        list_of_dict[1]['cipd'].append(fnormal[i])
                        list_of_dict[0]['gfipd'].append(iframenorm[i])
                        list_of_dict[1]['cfipd'].append(fframenorm[i])
                        list_of_dict[0]['gipdr'].append(inipdr[i])
                        list_of_dict[1]['cipdr'].append(fnipdr[i])
                        list_of_dict[0]['gfipdr'].append(iripdr[i])
                        list_of_dict[1]['cfipdr'].append(fripdr[i])
                
                list_of_dict[2]['aipd'] = list_of_dict[0]['aipd'] + list_of_dict[1]['aipd']
                list_of_dict[2]['tipd'] = list_of_dict[0]['tipd'] + list_of_dict[1]['tipd']
                list_of_dict[2]['cipd'] = list_of_dict[0]['cipd'] + list_of_dict[1]['cipd']
                list_of_dict[2]['gipd'] = list_of_dict[0]['gipd'] + list_of_dict[1]['gipd']
                list_of_dict[2]['afipd'] = list_of_dict[0]['afipd'] + list_of_dict[1]['afipd']
                list_of_dict[2]['tfipd'] = list_of_dict[0]['tfipd'] + list_of_dict[1]['tfipd']
                list_of_dict[2]['cfipd'] = list_of_dict[0]['cfipd'] + list_of_dict[1]['cfipd']
                list_of_dict[2]['gfipd'] = list_of_dict[0]['gfipd'] + list_of_dict[1]['gfipd']
                list_of_dict[2]['aipdr'] = list_of_dict[0]['aipdr'] + list_of_dict[1]['aipdr']
                list_of_dict[2]['tipdr'] = list_of_dict[0]['tipdr'] + list_of_dict[1]['tipdr']
                list_of_dict[2]['cipdr'] = list_of_dict[0]['cipdr'] + list_of_dict[1]['cipdr']
                list_of_dict[2]['gipdr'] = list_of_dict[0]['gipdr'] + list_of_dict[1]['gipdr']
                list_of_dict[2]['afipdr'] = list_of_dict[0]['afipdr'] + list_of_dict[1]['afipdr']
                list_of_dict[2]['tfipdr'] = list_of_dict[0]['tfipdr'] + list_of_dict[1]['tfipdr']
                list_of_dict[2]['cfipdr'] = list_of_dict[0]['cfipdr'] + list_of_dict[1]['cfipdr']
                list_of_dict[2]['gfipdr'] = list_of_dict[0]['gfipdr'] + list_of_dict[1]['gfipdr']
                list_of_dict[0]['nonafipdr'] = list_of_dict[0]['tfipdr'] + list_of_dict[0]['cfipdr'] + list_of_dict[0]['gfipdr']
                list_of_dict[1]['nonafipdr'] = list_of_dict[1]['tfipdr'] + list_of_dict[1]['cfipdr'] + list_of_dict[1]['gfipdr']
                list_of_dict[2]['nonafipdr'] = list_of_dict[0]['nonafipdr'] + list_of_dict[1]['nonafipdr']
                
                aipdall_m, aipdall_std, aipdall_var, aipds0_m, aipds0_std, aipds0_var, aipds1_m, aipds1_std, aipds1_var = calbase_statis(list_of_dict[2]['aipd'], list_of_dict[0]['aipd'], list_of_dict[1]['aipd'])
                afipdall_m, afipdall_std, afipdall_var, afipds0_m, afipds0_std, afipds0_var, afipds1_m, afipds1_std, afipds1_var = calbase_statis(list_of_dict[2]['afipd'], list_of_dict[0]['afipd'], list_of_dict[1]['afipd'])
                aipdrall_m, aipdrall_std, aipdrall_var, aipdrs0_m, aipdrs0_std, aipdrs0_var, aipdrs1_m, aipdrs1_std, aipdrs1_var = calbase_statis(list_of_dict[2]['aipdr'], list_of_dict[0]['aipdr'], list_of_dict[1]['aipdr'])
                afipdrall_m, afipdrall_std, afipdrall_var, afipdrs0_m, afipdrs0_std, afipdrs0_var, afipdrs1_m, afipdrs1_std, afipdrs1_var = calbase_statis(list_of_dict[2]['afipdr'], list_of_dict[0]['afipdr'], list_of_dict[1]['afipdr'])
                
                tipdall_m, tipdall_std, tipdall_var, tipds0_m, tipds0_std, tipds0_var, tipds1_m, tipds1_std, tipds1_var = calbase_statis(list_of_dict[2]['tipd'], list_of_dict[0]['tipd'], list_of_dict[1]['tipd'])
                tfipdall_m, tfipdall_std, tfipdall_var, tfipds0_m, tfipds0_std, tfipds0_var, tfipds1_m, tfipds1_std, tfipds1_var = calbase_statis(list_of_dict[2]['tfipd'], list_of_dict[0]['tfipd'], list_of_dict[1]['tfipd'])
                tipdrall_m, tipdrall_std, tipdrall_var, tipdrs0_m, tipdrs0_std, tipdrs0_var, tipdrs1_m, tipdrs1_std, tipdrs1_var = calbase_statis(list_of_dict[2]['tipdr'], list_of_dict[0]['tipdr'], list_of_dict[1]['tipdr'])
                tfipdrall_m, tfipdrall_std, tfipdrall_var, tfipdrs0_m, tfipdrs0_std, tfipdrs0_var, tfipdrs1_m, tfipdrs1_std, tfipdrs1_var = calbase_statis(list_of_dict[2]['tfipdr'], list_of_dict[0]['tfipdr'], list_of_dict[1]['tfipdr'])
                
                cipdall_m, cipdall_std, cipdall_var, cipds0_m, cipds0_std, cipds0_var, cipds1_m, cipds1_std, cipds1_var = calbase_statis(list_of_dict[2]['cipd'], list_of_dict[0]['cipd'], list_of_dict[1]['cipd'])
                cfipdall_m, cfipdall_std, cfipdall_var, cfipds0_m, cfipds0_std, cfipds0_var, cfipds1_m, cfipds1_std, cfipds1_var = calbase_statis(list_of_dict[2]['cfipd'], list_of_dict[0]['cfipd'], list_of_dict[1]['cfipd'])
                cipdrall_m, cipdrall_std, cipdrall_var, cipdrs0_m, cipdrs0_std, cipdrs0_var, cipdrs1_m, cipdrs1_std, cipdrs1_var = calbase_statis(list_of_dict[2]['cipdr'], list_of_dict[0]['cipdr'], list_of_dict[1]['cipdr'])
                cfipdrall_m, cfipdrall_std, cfipdrall_var, cfipdrs0_m, cfipdrs0_std, cfipdrs0_var, cfipdrs1_m, cfipdrs1_std, cfipdrs1_var = calbase_statis(list_of_dict[2]['cfipdr'], list_of_dict[0]['cfipdr'], list_of_dict[1]['cfipdr'])
                
                gipdall_m, gipdall_std, gipdall_var, gipds0_m, gipds0_std, gipds0_var, gipds1_m, gipds1_std, gipds1_var = calbase_statis(list_of_dict[2]['gipd'], list_of_dict[0]['gipd'], list_of_dict[1]['gipd'])
                gfipdall_m, gfipdall_std, gfipdall_var, gfipds0_m, gfipds0_std, gfipds0_var, gfipds1_m, gfipds1_std, gfipds1_var = calbase_statis(list_of_dict[2]['gfipd'], list_of_dict[0]['gfipd'], list_of_dict[1]['gfipd'])
                gipdrall_m, gipdrall_std, gipdrall_var, gipdrs0_m, gipdrs0_std, gipdrs0_var, gipdrs1_m, gipdrs1_std, gipdrs1_var = calbase_statis(list_of_dict[2]['gipdr'], list_of_dict[0]['gipdr'], list_of_dict[1]['gipdr'])
                gfipdrall_m, gfipdrall_std, gfipdrall_var, gfipdrs0_m, gfipdrs0_std, gfipdrs0_var, gfipdrs1_m, gfipdrs1_std, gfipdrs1_var = calbase_statis(list_of_dict[2]['gfipdr'], list_of_dict[0]['gfipdr'], list_of_dict[1]['gfipdr'])
                nonafipdrall_m, nonafipdrall_std, nonafipdrall_var, nonafipdrs0_m, nonafipdrs0_std, nonafipdrs0_var, nonafipdrs1_m, nonafipdrs1_std, nonafipdrs1_var = calbase_statis(list_of_dict[2]['nonafipdr'], list_of_dict[0]['nonafipdr'], list_of_dict[1]['nonafipdr'])
                
                f.write(f"{read.query_name},{ec},{'all'},{tnormal_m},{tnormal_std},{tnormal_var},{tframenorm_m},{tframenorm_std},{tframenorm_var},{tnipdr_m},{tnipdr_std},{tnipdr_var},{tripdr_m},{tripdr_std},{tripdr_var},{'0'},{inormal_m},{inormal_std},{inormal_var},{iframenorm_m},{iframenorm_std},{iframenorm_var},{inipdr_m},{inipdr_std},{inipdr_var},{iripdr_m},{iripdr_std},{iripdr_var},{'1'},{fnormal_m},{fnormal_std},{fnormal_var},{fframenorm_m},{fframenorm_std},{fframenorm_var},{fnipdr_m},{fnipdr_std},{fnipdr_var},{fripdr_m},{fripdr_std},{fripdr_var},{'all_A'},{aipdall_m},{aipdall_std},{aipdall_var},{afipdall_m},{afipdall_std},{afipdall_var},{aipdrall_m},{aipdrall_std},{aipdrall_var},{afipdrall_m},{afipdrall_std},{afipdrall_var},{'s0_A'},{aipds0_m},{aipds0_std},{aipds0_var},{afipds0_m},{afipds0_std},{afipds0_var},{aipdrs0_m},{aipdrs0_std},{aipdrs0_var},{afipdrs0_m},{afipdrs0_std},{afipdrs0_var},{'s1_A'},{aipds1_m},{aipds1_std},{aipds1_var},{afipds1_m},{afipds1_std},{afipds1_var},{aipdrs1_m},{aipdrs1_std},{aipdrs1_var},{afipdrs1_m},{afipdrs1_std},{afipdrs1_var},{'all_T'},{tipdall_m},{tipdall_std},{tipdall_var},{tfipdall_m},{tfipdall_std},{tfipdall_var},{tipdrall_m},{tipdrall_std},{tipdrall_var},{tfipdrall_m},{tfipdrall_std},{tfipdrall_var},{'s0_T'},{tipds0_m},{tipds0_std},{tipds0_var},{tfipds0_m},{tfipds0_std},{tfipds0_var},{tipdrs0_m},{tipdrs0_std},{tipdrs0_var},{tfipdrs0_m},{tfipdrs0_std},{tfipdrs0_var},{'s1_T'},{tipds1_m},{tipds1_std},{tipds1_var},{tfipds1_m},{tfipds1_std},{tfipds1_var},{tipdrs1_m},{tipdrs1_std},{tipdrs1_var},{tfipdrs1_m},{tfipdrs1_std},{tfipdrs1_var},{'all_C'},{cipdall_m},{cipdall_std},{cipdall_var},{cfipdall_m},{cfipdall_std},{cfipdall_var},{cipdrall_m},{cipdrall_std},{cipdrall_var},{cfipdrall_m},{cfipdrall_std},{cfipdrall_var},{'s0_C'},{cipds0_m},{cipds0_std},{cipds0_var},{cfipds0_m},{cfipds0_std},{cfipds0_var},{cipdrs0_m},{cipdrs0_std},{cipdrs0_var},{cfipdrs0_m},{cfipdrs0_std},{cfipdrs0_var},{'s1_C'},{cipds1_m},{cipds1_std},{cipds1_var},{cfipds1_m},{cfipds1_std},{cfipds1_var},{cipdrs1_m},{cipdrs1_std},{cipdrs1_var},{cfipdrs1_m},{cfipdrs1_std},{cfipdrs1_var},{'all_G'},{gipdall_m},{gipdall_std},{gipdall_var},{gfipdall_m},{gfipdall_std},{gfipdall_var},{gipdrall_m},{gipdrall_std},{gipdrall_var},{gfipdrall_m},{gfipdrall_std},{gfipdrall_var},{'s0_G'},{gipds0_m},{gipds0_std},{gipds0_var},{gfipds0_m},{gfipds0_std},{gfipds0_var},{gipdrs0_m},{gipdrs0_std},{gipdrs0_var},{gfipdrs0_m},{gfipdrs0_std},{gfipdrs0_var},{'s1_G'},{gipds1_m},{gipds1_std},{gipds1_var},{gfipds1_m},{gfipds1_std},{gfipds1_var},{gipdrs1_m},{gipdrs1_std},{gipdrs1_var},{gfipdrs1_m},{gfipdrs1_std},{gfipdrs1_var},{'nonA_all'},{nonafipdrall_m},{nonafipdrall_std},{nonafipdrall_var},{'nonA_s0'},{nonafipdrs0_m},{nonafipdrs0_std},{nonafipdrs0_var},{'nonA_s1'},{nonafipdrs1_m},{nonafipdrs1_std},{nonafipdrs1_var}\n")

def worker(read_queue, temp_dir):
    while True:
        reads = read_queue.get()
        if reads is None:
            read_queue.put(None)  # Signal to other workers that reading is done
            break
        temp_file = os.path.join(temp_dir, f"temp_{multiprocessing.current_process().pid}.txt")
        process_reads(reads, temp_file)

def main(bam_path, result_file, num_processes=10):
    read_queue = multiprocessing.Queue(maxsize=20)
    temp_dir = tempfile.mkdtemp()

    read_process = multiprocessing.Process(target=read_fetcher, args=(bam_path, read_queue, 30))
    read_process.start()
    
    workers = []
    for _ in range(num_processes):
        p = multiprocessing.Process(target=worker, args=(read_queue, temp_dir))
        p.start()
        workers.append(p)
    
    read_process.join()
    for _ in range(num_processes):
        read_queue.put(None)
    
    for p in workers:
        p.join()

    with open(result_file, 'w') as f_out:
        for temp_file in os.listdir(temp_dir):
            temp_file_path = os.path.join(temp_dir, temp_file)
            with open(temp_file_path, 'r') as f_temp:
                shutil.copyfileobj(f_temp, f_out)
    
    shutil.rmtree(temp_dir)

if __name__ == "__main__":
    bam_path = sys.argv[1]  # Path to your unaligned BAM file
    result_file = sys.argv[2]  # Output file for results
    main(bam_path, result_file, num_processes=10)
