#!/usr/bin/env bash

### inst/methylation/create_plugin_zip.sh TypeSeqHPV-Methyl.1_0_0.zip 
zip_fn=$1

rm -fr TypeSeqHPV-Methyl
mkdir -p TypeSeqHPV-Methyl

# ion torrent plugin specific files
cp inst/methylation/instance.html TypeSeqHPV-Methyl/
cp inst/methylation/launch.sh TypeSeqHPV-Methyl/
cp inst/methylation/plan.html TypeSeqHPV-Methylv/
cp inst/methylation/pluginsettings.json TypeSeqHPV-Methyl/

# zip plugin package (change the target zip file name accordinlgy)
zip -r $zip_fn  TypeSeqHPV-Methyl
rm -fr TypeSeqHPV-Methyl
