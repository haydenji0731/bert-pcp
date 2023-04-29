#!/bin/sh

SCRIPT_PATH=$1
MODEL_PATH=$2
DATA_PATH="./data/kmer$KMER/train/"
OUTPUT_PATH="./ft/kmer$KMER/"

# TODO: consider changing kmer size to 6
KMER=3
BSIZE=128
EPOCH=3

if [ -d $OUTPUT_PATH ]
then
  mkdir $OUTPUT_PATH
fi

# TODO: check bsize
python $SCRIPT_PATH --model_type dna --tokenizer_name dna$KMER \
  --model_name_or_path $MODEL_PATH --task_name dnaprom \
  --do_train --do_eval --data_dir $DATA_PATH --max_seq_length 100 \
  --per_gpu_eval_batch_size=$BSIZE --per_gpu_train_batch_size=$BSIZE \
  --learning_rate 2e-4 --num_train_epochs $EPOCH --output_dir $OUTPUT_PATH \
  --evaluate_during_training --logging_steps 100 --save_steps 4000 \
  --warmup_percent 0.1 --hidden_dropout_prob 0.1 --overwrite_output \
  --weight_decay 0.01 --n_process 8
