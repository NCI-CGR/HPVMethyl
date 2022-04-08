#!/bin/bash
set -x
# HPV Methyl
VERSION="1.0.0"
#autorundisable
echo Pipeline version $VERSION

ln ../../*.bam ./

docker run -i -v $(pwd):/mnt -v /mnt:/user_files \
    cgrlab/hpvmethyl:dev_v1.0.0 \
        Rscript /HPVMethyl/workflows/ion_methyl_workflow.R \
        --is_torrent_server yes \
        --config_file config_file.csv \
        --barcode_file barcodes.csv \
        --control_definitions control_defs.csv \
        --control_freq control_freq.csv \
        --cores 22 \
        --hotspot_vcf TS2-T52_v1-HOTSPOT.hotspot.vcf \
        --manifest manifest.csv \
        --ram 80G \
        --reference reference.fasta \
        --region_bed region.bed \
        --tvc_cores 4 \
        --tvc_parameters local_parameters.json

rm *rawlib.bam
