import pysam
import multiprocessing
import queue
import threading
import sys
import os
import tempfile
import shutil


os.environ['TMPDIR'] = os.getcwd()

def reverse_complement(nucleotide):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement[base] for base in nucleotide[::-1])

def read_fetcher(bam_path, read_queue, chunk_size):
    """
    Fetch reads from the BAM file and place them in a queue for processing.
    """
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam_file:
        buffer = []
        for read in bam_file.fetch(until_eof=True):
            buffer.append(read.to_dict())
            if len(buffer) >= chunk_size:
                read_queue.put(buffer)
                buffer = []
        if buffer:
            read_queue.put(buffer)  # Put the last chunk in the queue
    read_queue.put(None)  # Signal that reading is complete


def process_reads(reads, result_file):
    """
    Process a list of reads to find 'Ax' dinucleotides and their qualities, and write to file.
    """
    with open(result_file, 'a') as f:  # Open file in append mode
        for read_dict in reads:
            read = pysam.AlignedSegment.from_dict(read_dict, None)
            seq = read.query_sequence
            try:
                ec = read.get_tag('ec')
                fcontrol =  read.get_tag('FC')
                #fnormal =  read.get_tag('FZ')
                fframenorm =  read.get_tag('FF')
                #fnipdr =  read.get_tag('FN')
                fripdr =  read.get_tag('FR')
                #fripdr_trim = read.get_tag('FT')
                icontrol =  read.get_tag('IC')
                icontrol = icontrol[::-1]
                #inormal =  read.get_tag('IZ')
                #inormal = inormal[::-1]
                iframenorm  = read.get_tag('IF')
                iframenorm = iframenorm[::-1]
                #inipdr = read.get_tag('IN')
                #inipdr = inipdr[::-1]
                iripdr = read.get_tag('IR')
                iripdr = iripdr[::-1]
                #iripdr_trim = read.get_tag('IT')
                #iripdr_trim = iripdr_trim[::-1]
            except KeyError as e:
                print(f"Missing tag in read {read.query_name}:{e}")
                continue
 
            #print(len(seq), len(ripd), len(rpw))
            if len(fcontrol) == len(fframenorm) ==  len(fripdr) == len(seq) and int(ec+1) >= 30:
                for i in range(len(seq)):
                    if seq[i] == 'A':
                        if i+1 < len(seq) and seq[i+1] in 'GATC':
                            dinucleotide = seq[i] + seq[i+1]
                            rev_comp_dinuc = reverse_complement(dinucleotide)
                            f.write(f"{read.query_name},{0},{i+1},{'A'},{iframenorm[i]},{icontrol[i]},{iripdr[i]},{dinucleotide},{ec}\n")
                            #f.write(f"{read.query_name},{1},{i+1},{'T'},{fnormal[i]},{fframenorm[i]},{fcontrol[i]},{fnipdr[i]},{fripdr[i]},{rev_comp_dinuc},{ec}\n")
                        else:
                           f.write(f"{read.query_name},{0},{i+1},{'A'},{iframenorm[i]},{icontrol[i]},{iripdr[i]},{'A'},{ec}\n")
                           #f.write(f"{read.query_name},{1},{i+1},{'T'},{fnormal[i]},{fframenorm[i]},{fcontrol[i]},{fnipdr[i]},{fripdr[i]},{'T'},{ec}\n")
                    if seq[i] == 'T':
                        if i>0 and  seq[i-1] in 'GATC':
                           dinucleotide = seq[i-1] + seq[i]
                           rev_comp_dinuc = reverse_complement(dinucleotide)
                           #f.write(f"{read.query_name},{0},{i+1},{'T'},{inormal[i]},{iframenorm[i]},{icontrol[i]},{inipdr[i]},{iripdr[i]},{dinucleotide},{ec}\n")
                           f.write(f"{read.query_name},{1},{i+1},{'A'},{fframenorm[i]},{fcontrol[i]},{fripdr[i]},{rev_comp_dinuc},{ec}\n")
                        else:
                           #f.write(f"{read.query_name},{0},{i+1},{'T'},{inormal[i]},{iframenorm[i]},{icontrol[i]},{inipdr[i]},{iripdr[i]},{'T'},{ec}\n")
                           f.write(f"{read.query_name},{1},{i+1},{'A'},{fframenorm[i]},{fcontrol[i]},{fripdr[i]},{'A'},{ec}\n")

def worker(read_queue, temp_dir):
    while True:
        reads = read_queue.get()
        if reads is None:
            read_queue.put(None)  # Signal to other workers that reading is done
            break
        temp_file = os.path.join(temp_dir, f"temp_{multiprocessing.current_process().pid}.txt")
        process_reads(reads, temp_file)

def main(bam_path, result_file, num_processes=10):
    # Queue for thread communication
    read_queue = multiprocessing.Queue(maxsize=100)
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

# Example usage
if __name__ == "__main__":
    bam_path = sys.argv[1]  # Path to your unaligned BAM file
    result_file = sys.argv[2]  # Output file for results
    main(bam_path, result_file, num_processes=10)
