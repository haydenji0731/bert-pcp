#!/bin/sh

SCRIPT_PATH=$1
MODEL_PATH="./ft/kmer$KMER/"
DATA_PATH="./data/kmer$KMER/test/"
PRED_PATH="./results/kmer$KMER"

# TODO: consider changing kmer size to 6
KMER=3
BSIZE=128
EPOCH=3

if [ -d $PRED_PATH ]
then
  mkdir $PRED_PATH
fi

# TODO: consider ensemble prediction
# TODO: check max_seq_length parameter
# TODO: check overwrite_output parameter
python $SCRIPT_PATH --model_type dna --tokenizer_name=dna$KMER \
  --model_name_or_path $MODEL_PATH --task_name dnaprom --do_predict \
  --data_dir $DATA_PATH --max_seq_length 100 --per_gpu_pred_batch_size=$BSIZE \
  --output_dir $MODEL_PATH --predict_dir $PRED_PATH --n_process 32