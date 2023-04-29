import argparse
import sys
from time import strftime
import os
from motif_utils import seq2kmer
import pandas as pd

OUTPUT_FILE = "original.tsv"
label_d = dict()


def kmerize(args):
    fa_fh = open(args.in_fasta, 'r')
    label_fh = open(args.label, 'r')
    global label_d
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    out_fn = args.out_dir + OUTPUT_FILE
    out_fh = open(out_fn, 'w')
    out_fh.write("sequence\tlabel\n")
    for line in label_fh:
        tid = line.split("\t")[0]
        label = int(line.split("\t")[1])
        label_d[tid] = label
    label_fh.close()
    for line in fa_fh:
        if line[0] == '>':
            if len(line.split(" ")) > 1:
                tid = line.replace('>', '').split(" ")[0]
            else:
                tid = line.replace('>', '')[:-1]
            label = label_d[tid]
        else:
            out_fh.write(seq2kmer(line[:-1], args.kmer_size) + "\t" + str(label) + "\n")
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
    # TODO: is this correct?
    df_test = pd.merge(test_pos, test_neg, how='outer').sample(frac=1, random_state=args.seed)
    pos_re = pos[~pos.sequence.isin(test_pos.sequence)].dropna()
    neg_re = neg[~neg.sequence.isin(test_neg.sequence)].dropna()
    val_pos = pos_re.sample(frac=args.val_frac, random_state=args.seed)
    val_neg = neg_re.sample(frac=args.val_frac, random_state=args.seed)
    df_val = pd.merge(val_pos, val_neg, how='outer').sample(frac=1, random_state=args.seed)
    train_pos = pos_re[~pos_re.sequence.isin(val_pos.sequence)].dropna()
    train_neg = neg_re[~neg_re.sequence.isin(val_neg.sequence)].dropna()
    df_train = pd.merge(train_pos, train_neg, how='outer').sample(frac=1, random_state=args.seed)

    test_dir = os.path.join(args.out_dir, "test")
    train_val_dir = os.path.join(args.out_dir, "train")

    if not os.path.isdir(test_dir):
        os.makedirs(test_dir)
    if not os.path.isdir(train_val_dir):
        os.makedirs(train_val_dir)

    test_df_fn = os.path.join(test_dir, "dev.tsv")
    df_test.to_csv(test_df_fn, sep='\t', index=False)
    val_df_fn = os.path.join(train_val_dir, "dev.tsv")
    df_val.to_csv(val_df_fn, sep='\t', index=False)
    train_df_fn = os.path.join(train_val_dir, "train.tsv")
    df_train.to_csv(train_df_fn, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(description="kmer-ize input FASTA for training/testing")
    parser.add_argument("-in_fasta", type=str, help="full path please", required=True)
    parser.add_argument("-label", type=str, help="", required=True)
    parser.add_argument("-out_dir", type=str, help="", required=True)
    # consider changing this to 6; DNABERT used 6 instead of 3
    parser.add_argument("--kmer_size", type=int, help="", default=3, required=False)
    parser.add_argument("--seed", type=int, help="", default=0, required=False)
    parser.add_argument("--train_frac", type=float, help="", default=0.8, required=False)
    parser.add_argument("--val_frac", type=float, help="", default=0.1, required=False)
    parser.add_argument("--test_frac", type=float, help="", default=0.1, required=False)
    args = parser.parse_args()
    train_frac = args.train_frac
    val_frac = args.val_frac
    test_frac = args.test_frac
    if train_frac + val_frac + test_frac != 1.0:
        sys.exit("")
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Starting kmer-izing sequences")
    kmerize(args)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Kmer-ization completed. Extracting train/val/test sets.")
    split_data(args)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Pre-processing completed. Proceed with model training.")


if __name__ == "__main__":
    main()
