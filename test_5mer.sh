#!/bin/sh

KMER=3
BSIZE=128
EPOCH=3

SCRIPT_PATH=$1
IS_BASE=$2
if [[ $IS_BASE == 1 ]]
then
  MODEL_PATH="../3-new-12w-0/"
  PRED_PATH="./results/base_kmer$KMER"
else
  MODEL_PATH="./ft/kmer$KMER/"
  PRED_PATH="./results/kmer$KMER"
fi

DATA_PATH="./data/kmer$KMER/test/"

echo "checking if the prediction directory is present"
if [ -d $PRED_PATH ]
then
  echo "prediction directory found"
else
  echo "prediction directory not found; making it"
  mkdir $PRED_PATH
fi

python $SCRIPT_PATH --model_type dna --tokenizer_name=dna$KMER \
  --model_name_or_path $MODEL_PATH --task_name dnaprom --do_predict \
  --data_dir $DATA_PATH --max_seq_length 100 --per_gpu_pred_batch_size=$BSIZE \
  --output_dir $MODEL_PATH --predict_dir $PRED_PATH --n_process 32