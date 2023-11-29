// based on merge_fits_hor
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include <bg_globals.h>
#include <bg_fits.h>
#include <bg_array.h>

#include <mystring.h>
#include <myfile.h>

#include <vector>
using namespace std;

double gTimesampleInSec=1.08/1000000.00;
double gFullBW = (400.00/512.00)*(32.00/27.00);
double gChannel2FreqMultiplier = (400.00/512.00);
double gOversamplingRatio = (32.00/27.00);
double alpha = (400.00/512.00);
double beta  = (32.00/27.00);

string fits_list="fits_list";
string out_fits_base="out.fits";
string out_fits="out.fits";
int gOutFileCounter=0;
string outdir;
int gVerb=0;
int gInitialSizeY=-1;
int gMaxYSize=10000;

void usage()
{
   printf("skalow_stitch_data FITS_LIST FITS_OUT -o OUTDIR\n");
   printf("\n");
   printf("Options :\n");
   exit(-1);
}

void print_parameters()
{
   printf("#####################################\n");
   printf("PARAMETERS :\n");
   printf("#####################################\n");
   printf("Outdir            = %s\n",outdir.c_str());
   printf("verb level        = %d\n",gVerb);
   printf("#####################################\n");   
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "vo:f:y:m:";
   int opt,opt_param,i;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
         case 'o':
            if( optarg ){
               outdir = optarg;
            }
            break;
         case 'y':
            if( optarg ){
               gInitialSizeY = atol(optarg);
            }
            break;
/*         case 'f':
            if( optarg ){
               out_file_base = optarg;
            }
            break;*/
         case 'm':
            if( optarg ){
               gMaxYSize = atol(optarg);
            }
            break;
         case 'v':
            gVerb++;
            break;
         case 'h':
            usage();
            break;
         default:  
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();
      }
   }
}
                                                        
                                                                 
int main(int argc,char* argv[])
{
  if( argc < 3 ){
     usage();
  }

  fits_list = argv[1];
  out_fits_base = argv[2];
  out_fits = out_fits_base.c_str();
  
  // parse command line :
  parse_cmdline(argc-3,argv+3);
  print_parameters();
  
  if( strlen(outdir.c_str()) > 0 ){
     char szTmp[1024];
     sprintf(szTmp,"%s/%s",outdir.c_str(),out_fits_base.c_str());     
     out_fits = szTmp;
  }

  if( strlen(outdir.c_str()) ){
     mkdir(outdir.c_str());
  }
  
  // int bg_read_list( const char* file, vector<string>& out_list );
  vector<string> fits_list_vec;
  bg_read_list( fits_list.c_str(), fits_list_vec );

  CBgFits fits( fits_list_vec[0].c_str() );
  if( fits.ReadFits( fits_list_vec[0].c_str() ) ){
     printf("ERROR : could not read fits file %s\n",fits_list_vec[0].c_str());
     exit(-1);
  }

  int n_fine_ch = fits.GetXSize();
  int n_coarse_ch = fits_list_vec.size();
  int n_cc_center = n_coarse_ch -2;
  int nc = round( n_fine_ch*(beta-1)/beta );   
//  int n_out_channels = (n_coarse_ch-2)*(n_fine_ch-nc-1) + (n_fine_ch-nc/2)*2 - 1;
  int n_out_channels = (n_coarse_ch-2)*(n_fine_ch-nc) + (n_fine_ch-nc/2)*2;
  if( (n_out_channels % 2) == 1 ){
     // make it even :
     printf("WARNING : number of output channels = %d (odd number), will skip last channel to make it even = %d\n",n_out_channels,(n_out_channels-1));
//     n_out_channels = n_out_channels - 1;
  }  
  printf("INFO : number of fine channels = %d -> overlapping channels = %d (%.8f)\n",n_fine_ch,nc,(n_fine_ch*(beta-1)/beta));
  printf("INFO : final number of channels in the stitched FITS file = %d\n",n_out_channels);
  
  CBgFits out( fits.GetXSize(), n_out_channels );
  out.SetValue(0.00);
  // only using keywords from the LEFT file :
  out.SetKeysWithoutStates( fits.GetKeys() );
  double start_freq = fits.start_freq; // frequency starts from the first channel (nothing skipped at the start) + (nc/2)*fits.delta_freq;
  out.PrepareBigHornsHeaderTransposed( fits.dtime_fs + fits.dtime_fu/1000000.00, fits.inttime, start_freq, fits.delta_freq );
  // void PrepareBigHornsHeader( double ux_start, double _inttime, double freq_start, double delta_freq_mhz );
  
  int out_channel_index = 0;

/*
  - Fixed -> there should be 1748 fine channels !!! 
  - Proper equation is : n_fine_total = 2*(n_fine-nc/2) + (ncc-2)(n_fine-nc)
  - where nc=20 , from n_fine*(1-1/beta)


  - This is because lower branch has 10 points below crossing and upper branch has 9
    - So we count
      - Lowest-freq coarse channel : (n_fine - nc/2) , as nc/2 = 10 = 9 +1 we remove 9 points above the crossing and the crossing point it will be added as part of next coarse channel 
      - Next ncc-2 coarse channels (middle without edges) we count :
         - remove lowest nc/2=10  fine channels and keep the 11th point (crossing which is common to both but was removed previously)
         - remove 10 due to :  9 points above the crossing and the crossing point it will be added as part of next coarse channel 
         - (ncc-2)	*(n_fine-nc)
      - highest frequency channel :
        - Remove 10 lowest freq channels and keep the 11th point 
        - (n_fine-nc/2)
   Hence total number of fine channels in the stitched spectrum :
     -  N_total = 2*(n_fine - nc/2) + (ncc-2)*(n_fine-nc)
   If ncc is even -> N_total is also even !!!
*/   
   
  // add lowest frequency channel (skip only nc/2 channels at the top of the band but keep the overlapping point):
  for(int ch=0;ch<(n_fine_ch-(nc/2-1));ch++){
     for(int t=0;t<fits.GetXSize();t++){
        double val = fits.getXY(t,ch);
        out.setXY(t,ch,val);
     }
     out_channel_index++;     
  }
  // out filled up to channel 118 (inclusive), out_channel_index ends up to be 128-9 = 119 
  
//  out_channel_index--; // step back as the next channel needs to be averaged with last from previous
  // out_channel_index += (n_fine_ch-nc/2) - 1; // we need to average the first channel in the next coarse channel with the last in previous
  printf("Add channel 1 (skipped last %d fine channels), next fine channel to fill is %d\n",nc/2,out_channel_index);
  
  for(int file=1;file<fits_list_vec.size();file++){
     printf("PROGRESS : adding file %s\n",fits_list_vec[file].c_str());
     if( fits.ReadFits( fits_list_vec[file].c_str() ) ){
        printf("ERROR : could not read fits file %s\n",fits_list_vec[file].c_str());
        exit(-1);
     }
     
     // average channel nc/2+1 with last fine channel in previous coarse channel :
     for(int t=0;t<fits.GetXSize();t++){
        double val1 = out.getXY(t,out_channel_index-1); // adding to previously filled channel, for the 2nd coarse channel it will be 119-1=118 to be averaged with 10th (nc/2) from 2nd coarse channel
        double val2 = fits.getXY(t,nc/2); // to be averaged with 10th (nc/2) from 2nd coarse channel
        
        out.setXY( t,out_channel_index-1,(val1+val2)/2.00);
     }
//     out_channel_index++;

     // add next fine channels 
     int last_ch = (n_fine_ch-(nc/2-1)); // again skip the last 9 channels from the upper end, but keep 118th to be averaged with 10th fine channel from the next coarse cc.
     if( file == (fits_list_vec.size()-1) ){
        // add channels up to the end from the last file :
        last_ch = n_fine_ch;
     }     
     for(int ch=(nc/2+1);ch<last_ch;ch++){ // channel nc/2 was already added and averaged with the channel from previous Coarse Channel, so we need to start from nc/2+1 = 11th channel (counting from ZERO)
        for(int t=0;t<fits.GetXSize();t++){
           double val = fits.getXY(t,ch);
           out.setXY( t,out_channel_index,val); 
        }
      
        out_channel_index++;
     }
     
  }
  printf("CHECK : filled %d fine channels\n",out_channel_index);
  
  
  if( out.WriteFits(out_fits.c_str()) ){
     printf("ERROR : could not write output fits file %s\n",out_fits.c_str());
     exit(-1);
  }  
  printf("SUCCESS : written output fits file %s\n",out_fits.c_str());
}

