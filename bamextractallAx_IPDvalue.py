import pysam
import concurrent.futures
import queue
import threading
import sys

def reverse_complement(nucleotide):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement[base] for base in nucleotide[::-1])

def read_fetcher(bam_path, read_queue, chunk_size, num_threads):
    """
    Fetch reads from the BAM file and place them in a queue for processing.
    """
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam_file:
        buffer = []
        for read in bam_file.fetch(until_eof=True):
            buffer.append(read)
            if len(buffer) >= chunk_size:
                read_queue.put(buffer)
                buffer = []
        if buffer:
            read_queue.put(buffer)  # Put the last chunk in the queue
    for _ in range(num_threads):
        read_queue.put(None)  # Signal that reading is complete


def process_reads(reads, result_file):
    """
    Process a list of reads to find 'Ax' dinucleotides and their qualities, and write to file.
    """
    with open(result_file, 'a') as f:  # Open file in append mode
        for read in reads:
            seq = read.query_sequence
            try:
                ec = read.get_tag('ec')
                fcontrol =  read.get_tag('FC')
                fnormal =  read.get_tag('FZ')
                fframenorm =  read.get_tag('FF')
                fnipdr =  read.get_tag('FN')
                fripdr =  read.get_tag('FR')
                icontrol =  read.get_tag('IC')
                icontrol = icontrol[::-1]
                inormal =  read.get_tag('IZ')
                inormal = inormal[::-1]
                iframenorm  = read.get_tag('IF')
                iframenorm = iframenorm[::-1]
                inipdr = read.get_tag('IN')
                inipdr = inipdr[::-1]
                iripdr = read.get_tag('IR')
                iripdr = iripdr[::-1]
            except KeyError as e:
                print(f"Missing tag in read {read.query_name}:{e}")
                continue
 
            #print(len(seq), len(ripd), len(rpw))
            if len(fcontrol) == len(fnormal) == len(fframenorm) == len(fnipdr) == len(fripdr) == len(seq):
                for i in range(len(seq)):
                    if seq[i] == 'A':
                        if i+1 < len(seq) and seq[i+1] in 'GATC':
                            dinucleotide = seq[i] + seq[i+1]
                            rev_comp_dinuc = reverse_complement(dinucleotide)
                            f.write(f"{read.query_name},{0},{i+1},{'A'},{inormal[i]},{iframenorm[i]},{icontrol[i]},{inipdr[i]},{iripdr[i]},{dinucleotide},{ec}\n")
                            f.write(f"{read.query_name},{1},{i+1},{'T'},{fnormal[i]},{fframenorm[i]},{fcontrol[i]},{fnipdr[i]},{fripdr[i]},{rev_comp_dinuc},{ec}\n")
                        else:
                           f.write(f"{read.query_name},{0},{i+1},{'A'},{inormal[i]},{iframenorm[i]},{icontrol[i]},{inipdr[i]},{iripdr[i]},{'A'},{ec}\n")
                           f.write(f"{read.query_name},{1},{i+1},{'T'},{fnormal[i]},{fframenorm[i]},{fcontrol[i]},{fnipdr[i]},{fripdr[i]},{'T'},{ec}\n")
                    if seq[i] == 'T':
                        if i>0 and  seq[i-1] in 'GATC':
                           dinucleotide = seq[i-1] + seq[i]
                           rev_comp_dinuc = reverse_complement(dinucleotide)
                           f.write(f"{read.query_name},{0},{i+1},{'T'},{inormal[i]},{iframenorm[i]},{icontrol[i]},{inipdr[i]},{iripdr[i]},{dinucleotide},{ec}\n")
                           f.write(f"{read.query_name},{1},{i+1},{'A'},{fnormal[i]},{fframenorm[i]},{fcontrol[i]},{fnipdr[i]},{fripdr[i]},{rev_comp_dinuc},{ec}\n")
                        else:
                           f.write(f"{read.query_name},{0},{i+1},{'T'},{inormal[i]},{iframenorm[i]},{icontrol[i]},{inipdr[i]},{iripdr[i]},{'T'},{ec}\n")
                           f.write(f"{read.query_name},{1},{i+1},{'A'},{fnormal[i]},{fframenorm[i]},{fcontrol[i]},{fnipdr[i]},{fripdr[i]},{'A'},{ec}\n")

                           
                        

def main(bam_path, result_file, num_threads=4):
    # Queue for thread communication
    read_queue = queue.Queue(maxsize=10)  # Limit the queue size to control memory usage

    # Start the read fetcher thread
    threading.Thread(target=read_fetcher, args=(bam_path, read_queue, 20 , num_threads), daemon=True).start()

    # Process reads in parallel
    def worker():
        while True:
            reads = read_queue.get()
            if reads is None:
                read_queue.task_done()
                break  # If None, signal to stop processing
            process_reads(reads, result_file)
            read_queue.task_done()

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        for _ in range(num_threads):
            executor.submit(worker)

        read_queue.join()  # Ensure all items in the queue have been processed

# Example usage
if __name__ == "__main__":
    bam_path = sys.argv[1]  # Path to your unaligned BAM file
    result_file = sys.argv[2]  # Output file for results
    main(bam_path, result_file, num_threads=4)
