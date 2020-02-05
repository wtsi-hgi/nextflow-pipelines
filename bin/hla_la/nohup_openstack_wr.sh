#!/usr/bin/env bash

export PATH=/home/ubuntu:$PATH

rm -f hs_err_pid* && rm -f timeline* && rm -f trace* && rm -rf report* && rm -f bsub.o && rm -f bsub.e && rm -f .nextflow.log && rm -f nohup.log && rm -f nohup.out && nohup ./nextflow-pipelines/bin/ha_la/nextflow.sh > nohup.log 2>&1 &
