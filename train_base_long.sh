#!/bin/sh

KMER=3
BSIZE=4
EPOCH=1

SCRIPT_PATH=$1
MODEL_PATH=$2
DATA_PATH="./data/kmer$KMER/train/"
OUTPUT_PATH="./ft/base_kmer$KMER/"

echo "checking if the output directory is present"
if [ -d $OUTPUT_PATH ]
then
  echo "output directory found"
else
  echo "output directory not found; making it"
  mkdir $OUTPUT_PATH
fi

python $SCRIPT_PATH --model_type dnalong --tokenizer_name dna$KMER \
  --model_name_or_path $MODEL_PATH --task_name dnaprom \
  --do_train --data_dir $DATA_PATH --max_seq_length 2048 \
  --per_gpu_eval_batch_size=$BSIZE --per_gpu_train_batch_size=$BSIZE \
  --learning_rate 1e-5 --output_dir $OUTPUT_PATH \
  --logging_steps 1600 --save_steps 100000000000000000000 \
  --warmup_percent 0.1 --hidden_dropout_prob 0.1 --overwrite_output \
  --weight_decay 0.01 --n_process 4 --max_steps 1