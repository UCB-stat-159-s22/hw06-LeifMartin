#!/bin/bash -i
conda init bash

mamba env create -f environment.yml -p $HOME/envs/ligo

echo ello
echo $HOME
conda env list

conda activate ligo

conda env list


pip install pytest

python2 -m ipykernel install --user --name ligo --display_name "Python2"