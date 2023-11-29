#!/bin/bash

datadir=/data_archive/
if [[ -n "$1" && "$1" != "-" ]]; then
   datadir=$1
fi

template="2021*pulsar*"
if [[ -n "$2" && "$2" != "-" ]]; then
   template=$2
fi

freq_start=410
if [[ -n "$3" && "$3" != "-" ]]; then
   freq_start=$3
fi

objects="J* B*"
if [[ -n "$4" && "$4" != "-" ]]; then
   objects="$4"
fi

n_cc_channels=16
if [[ -n "$5" && "$5" != "-" ]]; then
   n_cc_channels=$5
fi

force=1
if [[ -n "$6" && "$6" != "-" ]]; then
   force=$6
fi

conversion_options=""
if [[ -n "$7" && "$7" != "-" ]]; then
   conversion_options="$7"   
fi

n_channels=128
if [[ -n "$8" && "$8" != "-" ]]; then
   n_channels=$8
fi

period=0.0894189988
if [[ -n "$9" && "$9" != "-" ]]; then
   period=$9
fi

skip_bytes=4096
ext="dada"
if [[ -n "${10}" && "${10}" != "-" ]]; then
   ext=${10}
fi
if [[ ${ext} != "dada" ]]; then
   skip_bytes=0
fi

file_template="*"
if [[ -n "${11}" && "${11}" != "-" ]]; then
   file_template="${11}"
fi

options=""
if [[ -n "${12}" && "${12}" != "-" ]]; then
   options=${12}
else 
   options="-"
fi

fits2fil=0
if [[ -n "${13}" && "${13}" != "-" ]]; then
   fits2fil=${13}
fi

single_file=0 # for a single big file the loop has to be where, but if each channel goes to a separate .dada/.dat file the loop is inside skalow_station_fold.sh
if [[ -n "${14}" && "${14}" != "-" ]]; then
   single_file=${14}
fi

transposed=0
if [[ -n "${15}" && "${15}" != "-" ]]; then
   transposed=${15}
fi

n_avg=7
if [[ -n "${16}" && "${16}" != "-" ]]; then
   n_avg=${16}
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
echo "transposed    = $transposed"
echo "n_avg         = $n_avg"
echo "############################################"
date

echo "DEBUG : using script skalow_station_fold.sh from the location:"
which skalow_station_fold.sh

cc=0

if [[ $single_file -gt 0 ]]; then
   echo "DEBUG : single_file = $single_file -> running loop over channels here (as only one BIG file)"
   
   # 1 BIG FILE : all freq. channels in a single .dada/.dat file :
   while [[ $cc -lt $n_cc_channels ]]; 
   do
      echo "skalow_station_fold.sh $datadir $template $freq_start \"$objects\" $n_cc_channels $force - $n_channels $period $ext \"$file_template\" \"${options}\" $cc"
      skalow_station_fold.sh $datadir $template $freq_start "$objects" $n_cc_channels $force - $n_channels $period $ext "$file_template" "${options}" $cc
   
      cc=$(($cc+1))
   done
else
   echo "DEBUG : single_file = $single_file -> multiple files (one per .dada/.dat file) -> no loop here, loop over .dada/.dat files inside skalow_station_fold.sh will handle this"

   # multiple files .dada/.dat file per frequency channel :
   # 5th parameter = 1 is because there is a single coarse channel per dada file in this path (else)
   echo "skalow_station_fold.sh $datadir $template $freq_start \"$objects\" 1 $force - $n_channels $period $ext \"$file_template\" \"${options}\" $cc"
   skalow_station_fold.sh $datadir $template $freq_start "$objects" 1 $force - $n_channels $period $ext "$file_template" "${options}" $cc
fi

echo "skalow_station_stitch_channels.sh ${datadir} ${template} "${objects}" - ${fits2fil} ${transposed} ${n_avg}"
skalow_station_stitch_channels.sh ${datadir} ${template} "${objects}" - ${fits2fil} ${transposed} ${n_avg}

# added now stitching :
# ls ch???/dynspec_avg*_i.fits > fits_list
# ls ch???/dynspec_avg*_q.fits > fits_list_q
# ls ch???/dynspec_avg*_u.fits > fits_list_u
# ls ch???/dynspec_avg*_v.fits > fits_list_v

# echo "skalow_stitch_data fits_list all_channels_i.fits"
# skalow_stitch_data fits_list all_channels_i.fits

# echo "skalow_stitch_data fits_list_q all_channels_q.fits"
# skalow_stitch_data fits_list_q all_channels_q.fits

# echo "skalow_stitch_data fits_list_u all_channels_u.fits"
# skalow_stitch_data fits_list_u all_channels_u.fits

# echo "skalow_stitch_data fits_list_v all_channels_v.fits"
# skalow_stitch_data fits_list_v all_channels_v.fits


# echo "total_power all_channels_i.fits all_channels_i_total_power.txt -c 12"
# total_power all_channels_i.fits all_channels_i_total_power.txt -c 12

# if [[ $fits2fil -gt 0 ]]; then
#   echo "fits2fil2 all_channels_i.fits  all_channels_i.fil -u"
#   fits2fil2 all_channels_i.fits  all_channels_i.fil -u
#else
#   echo "WARNING : fits2fil conversion is not required"
#fi
