#!/bin/bash

macs2 callpeak -t ${sample2}.bwa.sort.filter.bam -f BAMPE --nomodel --shift --100 --extsize 200 -f BAM -g 2.2e9 --outdir ./
