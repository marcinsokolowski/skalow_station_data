#!/bin/bash

fits2fil=0
if [[ -n "$1" && "$1" != "-" ]]; then
   fits2fil=$1
fi


export PATH=$HOME/github/hdf5_correlator/scripts/:$PATH

done_file=processed.txt

echo "############################################"
echo "PARAMETERS of script skalow_station_fold_allcc.sh :"
echo "############################################"
echo "datadir       = $datadir"
echo "template      = $template"
echo "file template = $file_template"
echo "freq_start    = $freq_start"
echo "objects       = $objects"
echo "n_cc_channels = $n_cc_channels"
echo "n_channels    = $n_channels"
echo "force         = $force"
echo "period        = $period"
echo "file format   = $ext (skip $skip_bytes [bytes])"
echo "options       = $options"
echo "fits2fil      = $fits2fil"
echo "single_file   = $single_file"
echo "############################################"
date

# added now stitching :
ls ch???/dynspec_avg*_i.fits > fits_list
ls ch???/dynspec_avg*_q.fits > fits_list_q
ls ch???/dynspec_avg*_u.fits > fits_list_u
ls ch???/dynspec_avg*_v.fits > fits_list_v

echo "skalow_stitch_data fits_list all_channels_i.fits"
skalow_stitch_data fits_list all_channels_i.fits

echo "skalow_stitch_data fits_list_q all_channels_q.fits"
skalow_stitch_data fits_list_q all_channels_q.fits

echo "skalow_stitch_data fits_list_u all_channels_u.fits"
skalow_stitch_data fits_list_u all_channels_u.fits

echo "skalow_stitch_data fits_list_v all_channels_v.fits"
skalow_stitch_data fits_list_v all_channels_v.fits


echo "total_power all_channels_i.fits all_channels_i_total_power.txt -c 12"
total_power all_channels_i.fits all_channels_i_total_power.txt -c 12

if [[ $fits2fil -gt 0 ]]; then
   echo "fits2fil2 all_channels_i.fits  all_channels_i.fil -u"
   fits2fil2 all_channels_i.fits  all_channels_i.fil -u
else
   echo "WARNING : fits2fil conversion is not required"
fi
