import argparse
import sys
from time import strftime
import os
from motif_utils import seq2kmer
import pandas as pd

OUTPUT_FILE = "original.tsv"
label_d = dict()
gtf_d = dict()


def load_gtf(args):
    global gtf_d
    gtf_fh = open(args.in_gtf, 'r')
    for line in gtf_fh:
        if line.split("\t")[2] == "transcript":
            infos = line.split("\t")[8].split(";")
            for info in infos:
                key = info.strip().split(" ")[0].strip()
                val = info.strip().split(" ")[1].strip().replace('"', '')
                if key == "transcript_id":
                    tid = val
                elif key == "gene_id":
                    gid = val
                    break
        gtf_d[tid] = gid


def kmerize(args):
    fa_fh = open(args.in_fasta, 'r')
    label_fh = open(args.label, 'r')
    global label_d
    global gtf_d
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    out_fn = args.out_dir + OUTPUT_FILE
    out_fh = open(out_fn, 'w')
    out_fh.write("gene_id\tsequence\tlabel\n")
    for line in label_fh:
        tid = line.split("\t")[0]
        label = int(line.split("\t")[1])
        label_d[tid] = label
    label_fh.close()
    seq = ""
    # seq_mu = 0
    total_seq = 0
    for line in fa_fh:
        if line[0] == '>':
            if seq != "":
                # seq_mu += len(seq)
                # total_seq += 1
                if len(seq) < 522:
                    total_seq += 1
                    gid = gtf_d[tid]
                    out_fh.write(gid + "\t" + seq2kmer(seq, args.kmer_size).upper() + "\t" + str(label) + "\n")
                seq = ""
            if len(line.split(" ")) > 1:
                tid = line.replace('>', '').split(" ")[0]
            else:
                tid = line.replace('>', '')[:-1]
            label = label_d[tid]
        else:
            seq += line[:-1]
    if len(seq) < 522:
        gid = gtf_d[tid]
        out_fh.write(gid + "\t" + seq2kmer(seq[:-1], args.kmer_size).upper() + "\t" + str(label) + "\n")
    # seq_mu += len(seq)
    # total_seq += 1
    # if len(seq) > 512:
    #     over512 += 1
    # print(over512)
    # seq_mu /= total_seq
    print(total_seq)
    fa_fh.close()
    out_fh.close()


def split_data(args):
    seq_fn = args.out_dir + OUTPUT_FILE
    df_seq = pd.read_csv(seq_fn, sep='\t')
    df_seq = df_seq.dropna(axis=0)
    pos = df_seq[df_seq['label'] == 1]
    neg = df_seq[df_seq['label'] == 0]
    test_pos = pos.sample(frac=args.test_frac, random_state=args.seed)
    test_neg = neg.sample(frac=args.test_frac, random_state=args.seed)
    test_pos_sameG = pos[pos.gene_id.isin(test_pos.gene_id)].dropna()
    test_neg_sameG = neg[neg.gene_id.isin(test_neg.gene_id)].dropna()
    test_pos_comb = pd.merge(test_pos, test_pos_sameG, how='outer')
    test_neg_comb = pd.merge(test_neg, test_neg_sameG, how='outer')
    df_test = pd.merge(test_pos_comb, test_neg_comb, how='outer').sample(frac=1, random_state=args.seed)
    pos_re = pos[~pos.gene_id.isin(test_pos_comb.gene_id)].dropna()
    neg_re = neg[~neg.gene_id.isin(test_neg_comb.gene_id)].dropna()
    val_pos = pos_re.sample(frac=args.val_frac, random_state=args.seed)
    val_neg = neg_re.sample(frac=args.val_frac, random_state=args.seed)
    val_pos_sameG = pos_re[pos_re.gene_id.isin(val_pos.gene_id)].dropna()
    val_neg_sameG = neg_re[neg_re.gene_id.isin(val_pos.gene_id)].dropna()
    val_pos_comb = pd.merge(val_pos, val_pos_sameG, how='outer')
    val_neg_comb = pd.merge(val_neg, val_neg_sameG, how='outer')
    df_val = pd.merge(val_pos_comb, val_neg_comb, how='outer').sample(frac=1, random_state=args.seed)
    train_pos = pos_re[~pos_re.gene_id.isin(val_pos_comb.gene_id)].dropna()
    train_neg = neg_re[~neg_re.gene_id.isin(val_neg_comb.gene_id)].dropna()
    train_pos_len = len(train_pos.index)
    train_neg_len = len(train_neg.index)
    # train_neg_under = train_neg.sample(train_pos_len, random_state=args.seed)
    under_count = train_neg_len - train_pos_len
    # print(len(train_pos_under.index))
    # print(len(train_neg.index))
    train_neg_under = train_neg.sample(under_count, random_state=args.seed)
    # train_neg_comb = pd.concat([train_neg, train_neg_over], axis=0).sample(frac=1, random_state=args.seed)
    df_train = pd.merge(train_pos, train_neg_under, how='outer').sample(frac=1, random_state=args.seed)

    test_dir = os.path.join(args.out_dir, "test")
    train_val_dir = os.path.join(args.out_dir, "train")

    if not os.path.isdir(test_dir):
        os.makedirs(test_dir)
    if not os.path.isdir(train_val_dir):
        os.makedirs(train_val_dir)

    test_df_fn = os.path.join(test_dir, "dev.tsv")
    df_test = df_test.drop('gene_id', axis=1)
    df_test.to_csv(test_df_fn, sep='\t', index=False)
    val_df_fn = os.path.join(train_val_dir, "dev.tsv")
    df_val = df_val.drop('gene_id', axis=1)
    df_val.to_csv(val_df_fn, sep='\t', index=False)
    train_df_fn = os.path.join(train_val_dir, "train.tsv")
    df_train = df_train.drop('gene_id', axis=1)
    df_train.to_csv(train_df_fn, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(description="kmer-ize input FASTA for training/testing")
    parser.add_argument("-in_fasta", type=str, help="full path please", required=True)
    parser.add_argument("-in_gtf", type=str, help="", required=True)
    parser.add_argument("-label", type=str, help="", required=True)
    parser.add_argument("-out_dir", type=str, help="", required=True)
    # consider changing this to 6; DNABERT used 6 instead of 3
    parser.add_argument("--kmer_size", type=int, help="", default=3, required=False)
    parser.add_argument("--seed", type=int, help="", default=0, required=False)
    # parser.add_argument("--train_frac", type=float, help="", default=0.999, required=False)
    # parser.add_argument("--val_frac", type=float, help="", default=0.0005, required=False)
    # parser.add_argument("--test_frac", type=float, help="", default=0.0005, required=False)
    parser.add_argument("--train_frac", type=float, help="", default=0.8, required=False)
    parser.add_argument("--val_frac", type=float, help="", default=0.1, required=False)
    parser.add_argument("--test_frac", type=float, help="", default=0.1, required=False)
    args = parser.parse_args()
    train_frac = args.train_frac
    val_frac = args.val_frac
    test_frac = args.test_frac
    if round(train_frac + val_frac + test_frac) != 1.0:
        sys.exit("")
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Loading annotation file.")
    load_gtf(args)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Starting kmer-izing sequences")
    kmerize(args)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Kmer-ization completed. Extracting train/val/test sets.")
    split_data(args)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Pre-processing completed. Proceed with model training.")


if __name__ == "__main__":
    main()