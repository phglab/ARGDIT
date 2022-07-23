#!/bin/bash

target_path=''

function get_conda_env_bin_path()
{
    IFS=':' read -r -a env_paths <<< $PATH

    for env_path in $env_paths
    do
        if [[ $env_path =~ ^.+/.*conda.*/.+$ ]] && [[ $env_path =~ ^.+/envs/.+$ ]]; then
            target_path=$env_path
            break
        fi
    done
}

get_conda_env_bin_path

wget http://www.bioinf.ucd.ie/download/od-seq.tar.gz
tar -xvf od-seq.tar.gz
mv OD-Seq OD-Seq-src
cd OD-Seq-src
clang++ -fopenmp -o OD-seq \
    AliReader.cpp \
    Bootstrap.cpp \
    DistCalc.cpp \
    DistMatReader.cpp \
    DistWriter.cpp \
    FastaWriter.cpp \
    IQR.cpp \
    ODseq.cpp \
    PairwiseAl.cpp \
    Protein.cpp \
    ResultWriter.cpp \
    runtimeargs.cpp \
    util.cpp
cd ..

if [[ $target_path == '' ]]; then
    mv OD-seq-src/OD-seq .
else
    mv OD-seq-src/OD-seq $target_path
fi

rm -R OD-Seq-src
rm od-seq.tar.gz
