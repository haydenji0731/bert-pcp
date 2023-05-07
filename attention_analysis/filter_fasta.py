import argparse
import os
from motif_utils import seq2kmer
import numpy as np

transcript = list()
OUTPUT_TSV_FILE = "dev.tsv"
OUTPUT_FA_FILE = "dev.fa"
OUTPUT_NPY_FILE = "annotation.npy"


def retrieve_transcript(args):
    global transcript
    gtf_fh = open(args.in_gtf, 'r')
    target = False
    for line in gtf_fh:
        if line.split("\t")[2] == "transcript":
            if target:
                return tid, tstart
            infos = line.split("\t")[8].split(";")
            tstart = int(line.split("\t")[3])
            for info in infos:
                key = info.strip().split(" ")[0].strip()
                val = info.strip().split(" ")[1].strip().replace('"', '')
                if key == "transcript_id":
                    tid = val
                elif key == "gene_id":
                    gid = val
                    if gid == args.gene_name:
                        target = True
                    break
        else:
            if target:
                region = line.split("\t")[2]
                if region != "exon":
                    start = int(line.split("\t")[3])
                    end = int(line.split("\t")[4])
                    transcript.append((region, start, end))


def find_region(i):
    global transcript
    # 5' UTR / 3' UTR / CDS
    res = [0] * 3

    for i in range(i, i + 100):
        for j in range(len(transcript)):
            five_prime = True
            feature = transcript[j]
            start = feature[1]
            end = feature[2]
            if five_prime:
                if i < start:
                    res[0] = 1
                    break
                five_prime = False
            if start <= i < end:
                res[2] = 1
                break
            elif j == len(transcript) - 1:
                if end <= 1:
                    res[1] = 1
                    break
    return res


def kmerize(args, target, target_start):
    fa_fh = open(args.in_fasta, 'r')
    is_pc = args.pc
    annotation = np.empty((0,3), int)
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    out_tsv_fn = args.out_dir + OUTPUT_TSV_FILE
    out_tsv_fh = open(out_tsv_fn, 'w')
    out_tsv_fh.write("sequence\tlabel\n")
    out_fa_fn = args.out_dir + OUTPUT_FA_FILE
    out_fa_fh = open(out_fa_fn, 'w')
    out_npy_fn = args.out_dir + OUTPUT_NPY_FILE
    seq = ""
    terminate = False
    for line in fa_fh:
        if line[0] == '>':
            if terminate:
                break
            if len(line.split(" ")) > 1:
                tid = line.replace('>', '').split(" ")[0]
            else:
                tid = line.replace('>', '')[:-1]
        else:
            if tid == target:
                terminate = True
                seq += line[:-1]
    fa_fh.close()
    cnt = 1
    for i in range(0, len(seq), 100):
        out_tsv_fh.write(seq2kmer(seq[i:i+100], args.kmer_size) + "\t" + str(is_pc) + "\n")
        out_fa_fh.write(">" + tid + "-" + str(cnt) + "\n")
        cnt += 1
        out_fa_fh.write(seq[i:i+100] + "\n")
        res = find_region(i + target_start)
        annotation = np.append(annotation, np.array([res]), axis=0)
    np.save(out_npy_fn, annotation)
    out_tsv_fh.close()
    out_fa_fh.close()


def main():
    parser = argparse.ArgumentParser(description="generation annotation.npy file for attention analysis")
    parser.add_argument("-in_fasta", type=str, help="full path please", required=True)
    parser.add_argument("-in_gtf", type=str, help="", required=True)
    parser.add_argument("-gene_name", type=str, help="", required=True)
    parser.add_argument("-out_dir", type=str, help="", required=True)
    parser.add_argument("-pc", type=int, help="", default=1, required=False)
    parser.add_argument("--kmer_size", type=int, help="", default=3, required=False)
    args = parser.parse_args()
    target, target_start = retrieve_transcript(args)
    kmerize(args, target, target_start)


if __name__ == "__main__":
    main()