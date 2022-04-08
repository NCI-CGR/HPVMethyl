#!/usr/bin/env bash

# Example usage:
# inst/methylation/create_plugin_zip.sh  HPVMethyl HPVMethyl.1_0_0.zip
# inst/methylation/create_plugin_zip.sh  HPVMethyl-dev HPVMethyl-dev.1_0_0.zip

zip_dir=$1
zip_fn=$2

rm -fr $zip_dir
mkdir -p $zip_dir

# ion torrent plugin specific files
cp inst/methylation/{instance.html,launch.sh,plan.html,pluginsettings.json} $zip_dir

# zip plugin package (change the target zip file name accordinlgy)
zip -r $zip_fn  $zip_dir
rm -fr $zip_dir
