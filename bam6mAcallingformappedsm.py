import pysam
import numpy as np
import multiprocessing
import sys
import os
import tempfile
import shutil
import re

os.environ['TMPDIR'] = os.getcwd()

def reverse_complement(nucleotide):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement[base] for base in nucleotide[::-1])

def qlencigar(cigar):
    return sum(int(num) for num, type in re.findall(r'(\d+)([MIDNSHPX=])', cigar) if type in "=MXSI")

def hasharr(qs, qe, rs, re):
    if qs > qe:
        keys = range(qs, qe, -1)  # decrement if qs > qe
    else:
        keys = range(qs, qe)  # increment if qs <= qe

    if rs > re:
        values = range(rs, re , -1)  # decrement if rs > re
    else:
        values = range(rs, re)  # increment if rs <= re

    # Create a dictionary that maps keys to values
    hs_arr = {key: value for key, value in zip(keys, values)}

    return hs_arr

def cigarparse(cigar, qlen, refstart, direct):
    r_s = refstart
    q_s = 1 if direct == 0 else qlen
    hs_q2r = {}

    for num, type in re.findall(r'(\d+)([MIDNSHPX=])', cigar):
        num = int(num)

        if direct == 0:
            if type == "=":
                q_e = q_s + num
                r_e = r_s + num
                hs_q2r.update(hasharr(q_s, q_e , r_s, r_e ))
            elif type == 'D':
                r_e = r_e + num
            elif type == 'X':
                q_e = q_s + num
                r_e = r_s + num
            elif type == 'I':
                q_e = q_e + num
            elif type == 'S':
                q_s = q_s + num
                continue
            q_s = q_e
            r_s = r_e
        else:
            if type in "=":
                q_e = q_s - num
                r_e = r_s + num
                hs_q2r.update(hasharr(q_s, q_e , r_s, r_e ))
            elif type == 'X':
                r_e = r_s + num
                q_e = q_s - num
            elif type == 'D':
                r_e = r_e + num
            elif type == 'I':
                q_e = q_e - num
            elif type == 'S':
                q_s = q_s - num
                continue
            q_s = q_e
            r_s = r_e
    return r_s, hs_q2r


def read_fetcher(bam_path, read_queue, chunk_size):
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam_file:
        header_dict = bam_file.header.to_dict()
        buffer = []
        for read in bam_file.fetch(until_eof=True):
            buffer.append(read.to_dict())  # Convert read to a dictionary to ensure it's pickle-able
            if len(buffer) >= chunk_size:
                read_queue.put((buffer, header_dict))
                buffer = []
        if buffer:
            read_queue.put((buffer, header_dict))  # Put the last chunk in the queue
    read_queue.put(None)  # Signal that reading is complete

def process_reads(reads, result_file, cutoff, header_dict):
    header = pysam.AlignmentHeader.from_dict(header_dict)
    with open(result_file, 'a') as f:
        for read_dict in reads:
            methysite = []
            rs = ''
            qp_len = ''
            refat_pos = ''
            read = pysam.AlignedSegment.from_dict(read_dict, header)  # Convert dictionary back to read
            try:
                ec = read.get_tag('ec')
                qp_direct = 1 if read.is_reverse else 0
                seq = reverse_complement(read.query_sequence) if read.is_reverse else read.query_sequence
                cigarstring = read.cigarstring
                fripdr = read.get_tag('FR')
                iripdr = read.get_tag('IR')[::-1]
                qp_len = qlencigar(cigarstring)
                qp_refstart = read.reference_start
                ref_id = read.reference_id
                ref_name = header.get_reference_name(ref_id)
                rs, hs_qrc = cigarparse(cigarstring, qp_len, qp_refstart, qp_direct)
            except KeyError as e:
                print(f"Missing tag in read {read.query_name}:{e}")
                continue
            if len(seq) == len(fripdr) and int(ec + 1) >= 20:
                for i in range(len(seq)):
                    if seq[i] == 'A' and iripdr[i] >= float(cutoff):
                        methysite.append(i+1)
                    if seq[i] == 'T' and fripdr[i] >= float(cutoff):
                        methysite.append(i+1)
                refat_pos = [hs_qrc.get(m, -1) for m in methysite]
                f.write(f"{read.query_name}\t{ref_name}\t{ec}\t{qp_direct}\t{qp_len}\t{qp_refstart}\t{rs}\t{methysite}\t{refat_pos}\n")

def worker(read_queue, temp_dir, cutoff):
    while True:
        reads = read_queue.get()
        if reads is None:
            read_queue.put(None)  # Signal to other workers that reading is done
            break
        temp_file = os.path.join(temp_dir, f"temp_{multiprocessing.current_process().pid}.txt")
        reads_chunk, header_dict = reads
        process_reads(reads_chunk, temp_file, cutoff, header_dict)

def main(bam_path, result_file, num_processes=10, cutoff=1.932):
    read_queue = multiprocessing.Queue(maxsize=20)
    temp_dir = tempfile.mkdtemp()

    read_process = multiprocessing.Process(target=read_fetcher, args=(bam_path, read_queue, 30))
    read_process.start()
    
    workers = []
    for _ in range(num_processes):
        p = multiprocessing.Process(target=worker, args=(read_queue, temp_dir, cutoff))
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
    cutoff = float(sys.argv[3])
    num_processes = int(sys.argv[4])
    main(bam_path, result_file, num_processes, cutoff)
