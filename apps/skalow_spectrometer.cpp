#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <complex>
#include <vector>
#include <spectrometer.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <bg_fits.h>
#include <myfile.h>
#include <mystring.h>

#include "SigprocFile.h" // to save .fil files too

enum eBandPassCalMethod_T { eBWCal_donothing=0, eBWCal_divide=1, eBWCal_subtract=2 };
enum eOutputType_T { eUnknown=0, eI=1, eQ=2, eU=3, eV=4, eIQUV=5, ePhaseXY=6, eXY=7 };
eOutputType_T gOutPutType = eI;
bool gShiftReIm=true; // for now always do XY phase 

// saving filterbank files:
bool gSaveFilterbankFiles=false;
char gOutFilFile[2048];

eOutputType_T parse_output_type( const char* szOutputType ){
   if( strcasecmp(szOutputType,"iquv") == 0 ){ 
      return eIQUV;
   }
   if( strcasecmp(szOutputType,"pxy") == 0 || strcasecmp(szOutputType,"phase_xy") == 0 || strcasecmp(szOutputType,"phasexy") == 0 ){ 
      return ePhaseXY;
   }
   if( strcasecmp(szOutputType,"xy") == 0 ){
      return eXY;
   }
   
   if( szOutputType[0] == 'I' || szOutputType[0] == 'i' ){
      return eI;
   }
   if( szOutputType[0] == 'Q' || szOutputType[0] == 'q' ){
      return eQ;
   }
   if( szOutputType[0] == 'U' || szOutputType[0] == 'u' ){
      return eU;
   }
   if( szOutputType[0] == 'V' || szOutputType[0] == 'v' ){
      return eV;
   }
   
   return eUnknown;
}

// 8192 -> 8200 (for timestamp)
#define N_CHANNELS (24*32*4) // 10kHz channels 24*32*4 = 3072
#define MAX_N_SAMPLES 2*N_CHANNELS // Re/Im of 24*32 fine channels 
#define DATA_TYPE     short

double gTimesampleInSec=1.08/1000000.00;
double gFullBW = (400.00/512.00)*(32.00/27.00);
double gChannel2FreqMultiplier = (400.00/512.00);
int gFreqChannel = 410;

int gVerb=0;
int gDumpChannel=0;
char filename[1024];
char gOutFile[1024],gOutFileAvg[2048],gOutFileDat[2048],gOutFileTotPower[1024];
mystring gOutDir;
int gAddIndex=0;
int gPol=1; // 0 - X , 1 - Y 
int gNPols=2;
int gNChannels=1;
int gSaveNChannels=1;
int gOutFineChannels=32;

#define BYTES_PER_SAMPLE 4 // this is re/im for X and re/im for Y -> 2 x 2 bytes = 4 bytes 

int gNBytesPerSample = BYTES_PER_SAMPLE;
int n_samples = 1024 * BYTES_PER_SAMPLE; // number of samples to read, used to be N_CHANNELS; 
int gStartIndex=-1;
int gSkipNBytes=-1;
int gProcessNSamples = -1;
int gMultiplier=1;
int gCheckZeros=0;

double gCorrWithInvJones=-1; // normalisation factor

// maximum time to process
double gMaxTimeToProcessInSec=-1; 

// output FITS files :
string gOutFITSBase; // empty -> do not save FITS
int gAvgNSpectraToFITS=29; // 1ms / (1.08*32/1000.00) = 28.935185185 , 0.1ms = 2.8935185185 ~= 3 spectra 
double gIntegrationTimeSec=-1000; // 1ms <1000 -> no specified
bool gTransposedFITS=false; 

// folding :
int gNFoldingBinsCount=128;
double gFoldingPulsarPeriod=0.089328385024; // in seconds 
eBandPassCalMethod_T gBandpassCalMethod=eBWCal_donothing;
string gOutFoldedFits;
string gOutFoldedMaxVsFreq;

// observation direction:
string gTarget  = "B0950+08";
double gRA_deg  = 148.28750000;
double gDEC_deg = 7.92638889;

// # /home/msok/Desktop/Curtin/students/2022/Chris_Lee_Honours_Project/logbook/20220412_beam-correcting_VELA_data.odt
// std::complex<double> gInvJones00( -0.19391154 , -0.00342333 ); // -0.133340316282896 , 0.039036116382880 );
// std::complex<double> gInvJones01( 1.1539054   , +0.37828104 ); // -0.860154030652396 , 0.187848858551529 );
// std::complex<double> gInvJones10( -1.07909468 , -0.24393281 ); // 0.764011070552511 , -0.244298897390238 );
// std::complex<double> gInvJones11( -0.18581259 , -0.04685302 ); // -0.122434290299445 , 0.069153015019698 );

// # page 5 of /home/msok/Desktop/Curtin/students/2022/Chris_Lee_Honours_Project/logbook/20220509_VELA_obs_amps_and_InvAmps_BeamCorr.odt
std::complex<double> gInvJones00( 0.01824534 , 0.00622549 );
std::complex<double> gInvJones01( 0.07101837 , 0.02065822 );
std::complex<double> gInvJones10( -0.06648148 , -0.01180349 );
std::complex<double> gInvJones11(  0.01690944 , 0.0054325 );
bool gAddXYPhase=false;

double fine_channel_freq( int cc, int fine_ch, int n_fine_ch_total )
{
   double center_freq = cc*gChannel2FreqMultiplier;
   double start_freq  = center_freq - gFullBW/2.00;
   
   double fine_ch_bw = gFullBW/n_fine_ch_total;
   
   double fine_ch_freq = start_freq + fine_ch_bw/2.00 + fine_ch_bw*fine_ch;
   
   return fine_ch_freq;   
}

double GetFineChanFreqMHz( int fine_channel )
{
   double fine_channel_bw = gFullBW/gOutFineChannels;
   double freq_center = gFreqChannel * gChannel2FreqMultiplier;
   double freq_start = freq_center - gFullBW/2.00 + fine_channel_bw/2.00;

   double freq_mhz = freq_start + fine_channel*fine_channel_bw;
   return freq_mhz;
}


bool DoesFileExist(const char* fname)
{
        bool bRet=false;
        if(fname && fname[0] ){
                std::string szFile=fname;
                // szFile.env2str();

                if( access( szFile.c_str(), F_OK ) == 0 ) {
                        bRet = true;
                }
        }
        return bRet;
}


long int GetFileSize( const char* filename )
{
        if( !DoesFileExist( filename ) ){
                return -1;
        }

        struct stat buf;
        stat( filename, &buf );

        return (long int)buf.st_size;
}

void save_spectrum( const char* filename, std::vector<double>& spectrum )
{
   FILE* outf = fopen(filename,"w");
   for(int i=0;i<spectrum.size();i++){
      fprintf(outf,"%d %.8f\n",i,spectrum[i]);
   }
   fclose(outf);
}

void normalise_folded( CBgFits& folded )
{

   if( gBandpassCalMethod != eBWCal_donothing ){
      CBgArray mean_lines, rms_lines;
      folded.MeanLines( mean_lines, rms_lines );

      // folded :
//    pFoldedI->Divide( *pFoldedCounter ); // divide by counter of averaged spectra in each phase bin
      // normalise by the mean spectrum (WARNING : channels are in vertical direction - not horizontal)
      for(int fch=0;fch<folded.GetYSize();fch++){
         for(int t=0;t<folded.GetXSize();t++){
            double val = folded.getXY(t,fch);

            if( gBandpassCalMethod == eBWCal_divide ){
               folded.setXY(t,fch,val / mean_lines[fch]); // or - or / - divide is correct otherwise bandpass ripple stays there 
            }else{
               folded.setXY(t,fch,val - mean_lines[fch]); // may have to be running mean subtracted earlier 
            }
         }
      }
   }
}

void usage()
{
   printf("read_binary_station_beam_test_order1_2pol voltages_r.dat -d -u -c -f outfile -i integration_index -v -a AVERAGE_N_POINTS -p POL -s START_INDEX -n PROCESS_TIMESTAMPS -M multiplier -Z -S SAVE_N_CHANNELS\n");
   printf("\n");
   printf("Options :\n");
//   printf("\t-u - assumes that UNIXTIME timestamp is added to line as extra record of DATA_TYPE structure\n");
   printf("\t-d - increases debug level\n");
   printf("\t-v - increases debug level\n");   
   printf("\t-c CHANNEL : specifies channel to be dumped [default 0]\n");      
   printf("\t-C N_CHANNELS : number of channels [default %d]\n",gNChannels);
   printf("\t-f OUTFILENAME : name of file to which data is dumped\n");
   printf("\t-i INTEGRATION_TIME[sec] : integration time in seconds [default %.2f sec], <0 -> not specified\n",gIntegrationTimeSec);   
   printf("\t-a AVERAGE_N_POINTS : average N data points [default %d], <=1 -> do not average\n",gAvgNSpectraToFITS);
   printf("\t-p POL : POL = 0 = X , POL = 1 = Y , POL = -1 X+Y (power)\n");
   printf("\t-s START_INDEX : where to set start index [default %d] , <0 -> starting from the start of the file (and save these many bytes to file header_saved.txt\n",gStartIndex);   
   printf("\t-S SAVE_N_CHANNELS : how many channels to save into the output file [default %d]\n",gSaveNChannels);
   printf("\t-n PROCESS_TIMESTAMPS : number of timestamps to process [default %d and <0 is all]\n",gProcessNSamples);
   printf("\t-M multiplier : multiply timestep index (for example 128 fine channels * 7 averages to get ~1ms)\n");
   printf("\t-Z : check statistics of ZEROS\n");
   printf("\t-X N_BYTES : skip these many bytes [default %d]\n",gSkipNBytes);
   printf("\t-I : if >0 correct with inverse Jones and multiply by this factor [default %.4f]\n",gCorrWithInvJones);
   printf("\t-m time interval to process in seconds [default %.6f] , <=0 -> all\n",gMaxTimeToProcessInSec);
   printf("\t-O OUTPUT_FITSFILE : name of output FITS file basename (no extension .fits required) [default not set]\n");
   printf("\t-F frequency channel [default %d]\n",gFreqChannel);
   printf("\t-N N_FINE_CHANNELS [default %d]\n",gOutFineChannels);
   printf("\t-P FOLDING_PERIOD - period of the pulsar [default %.6f sec]\n",gFoldingPulsarPeriod);
   printf("\t-B N_PHASE_BINS   - number of folding phase bins [default %d]\n",gNFoldingBinsCount);
   printf("\t-T 0 or 1: transposed FITS file (time on vertical axis and frequency on horizontal) [default other way around]\n");
   printf("\t-D : divide by mean power spectrum\n");
   printf("\t-t OUTPUT_TYPE : iquv, i, q, u, v, or pxy - is Phase(xy*), xy* is correlation XY* etc\n");
   printf("\t-A OUTPUT_DIRECTORY : [default %s]\n",gOutDir.c_str());
   printf("\t-Y : add XY phase\n");
   printf("\t-b : save filterbank files [default disabled]\n");
   
   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "vudc:C:f:i:r:m:p:s:n:a:M:ZS:X:I:O:F:N:P:T:B:D:o:t:A:Yb";
   int opt,opt_param,i;
   
//   strcpy(filename,"");
   sprintf(gOutFile,"%s/out_file",gOutDir.c_str());
      
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
         case 'A':
            if( optarg ){
               gOutDir = optarg;
            }
            break;
         
         case 'b':
            gSaveFilterbankFiles = true;
            break;   
            
         case 'c':
            if( optarg ){   
               if( atol(optarg)>=0 && atol(optarg)<=512 ){
                 gDumpChannel = atol(optarg);
               }else{
                 printf("ERROR : wrong number of channels %d , not in range [0-512]!\n",(int)atol(optarg));
               }
            }
            break;

         case 'f':
            if( optarg ){   
               strcpy(gOutFile,optarg);
            }
            break;

         case 'D':
            gBandpassCalMethod = (eBandPassCalMethod_T)(atol(optarg));
            break;

         case 'F':
            if( optarg ){   
               gFreqChannel = atol(optarg);
            }
            break;

         case 'O':
            if( optarg ){   
               gOutFITSBase = optarg;
            }
            break;

         case 'o':
            if( optarg ){   
               gOutFoldedFits = optarg;
            }
            break;

         case 'i':
            if( optarg ){   
               gIntegrationTimeSec = atof(optarg);
            }
            break;

         case 'I':
            gCorrWithInvJones = atof(optarg);
            break;

         case 'a':
            if( optarg ){   
               gAvgNSpectraToFITS = atol(optarg);
            }
            break;

         case 'm':
            if( optarg ){   
               gMaxTimeToProcessInSec = atof(optarg);
            }
            break;

         case 'M':
            if( optarg ){   
               gMultiplier = atol(optarg);
            }
            break;

         case 'C':
            if( optarg ){   
               gNChannels = atol(optarg);
            }
            break;

         case 'p':
            if( optarg ){   
               gPol = atol(optarg);
               if( gPol >= 2 ){
                  gPol = 1;
               }
               if( gPol <= -1 ){
                  gPol = 0;
               }
            }
            break;
            
         case 'P':
            if( optarg ){   
               gFoldingPulsarPeriod = atof(optarg);
            }
            break;
            
         case 'B':
            if( optarg ){   
               gNFoldingBinsCount = atol(optarg);
            }
            break;
            
         case 'd':
            gVerb++;
            break;

         case 's':
            if( optarg ){
               gStartIndex = atol( optarg );
            }
            break;
            
         case 'S':
            if( optarg ){
               gSaveNChannels = atol( optarg );
            }
            break;
            
         case 'n':
            if( optarg ){
               gProcessNSamples = atol( optarg );
            }
            break;
            
         case 'N':
            if( optarg ){
               gOutFineChannels = atol( optarg );
            }
            break;
            
         case 'v':
            gVerb++;
            break;
            
         case 'X':
            if( optarg ){
               gSkipNBytes = atol( optarg );
            }
            break;
            
         case 'T':            
            if( optarg ){
               gTransposedFITS = ( atol( optarg ) > 0 );
            }else{            
               gTransposedFITS = true;
            }
            break;
            
         case 'Z':
            gCheckZeros = 1;
            break;
            
         case 't' :
            gOutPutType = parse_output_type( optarg );
            break;

         case 'Y' :
            gAddXYPhase = true;
            break;
         
         default:
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();                                          
      }
   }
   
   printf("INFO : frequency channel %d -> %d\n",gFreqChannel,(gFreqChannel + gDumpChannel));
   gFreqChannel = gFreqChannel + gDumpChannel;
   
/*   if( strlen(filename) == 0 ){
      if( gStartIndex > 0 && gProcessNSamples > 0 ){
         sprintf(filename,"%ld_%ld.dat",gStartIndex,(gStartIndex+gProcessNSamples));
      }else{
         strcpy(filename,"test.dat");
      }
   }*/

   if( strcmp(gOutDir.c_str(),"./") ){
      if( !strstr(gOutDir.c_str(),"/") ){
         gOutDir += "//";
      }
      MyFile::CreateDir(gOutDir.c_str());
      printf("INFO : created output directory %s\n",gOutDir.c_str());
   }

   if( (strlen(gOutDir.c_str()) + strlen(gOutFile) + 8) < 1024 ){
      sprintf(gOutFileAvg,"%s/%s_avg.txt",gOutDir.c_str(),gOutFile);
   }
   sprintf(gOutFileDat,"%s/%s.dat",gOutDir.c_str(),gOutFile);
   sprintf(gOutFileTotPower,"%s/total_power_out.txt",gOutDir.c_str());

   if( strlen(gOutFile) == 0 ){
      if( gStartIndex > 0 && gProcessNSamples > 0 ){
         sprintf(gOutFile,"%s/%d_%d.txt",gOutDir.c_str(),gStartIndex,(gStartIndex+gProcessNSamples));
      }else{
         sprintf(gOutFile,"%s/test.txt",gOutDir.c_str());
      }
      
      printf("DEBUG : output file name was empty -> set to %s\n",gOutFile);
   }   
   
   gStartIndex = gStartIndex * gMultiplier;
   gProcessNSamples = gProcessNSamples * gMultiplier;
   
   if( strlen(gOutFoldedFits.c_str())<=0 ){
      char szTmp[1024];
      sprintf(szTmp,"%s/folded_ch%d",gOutDir.c_str(),gFreqChannel);
      gOutFoldedFits=szTmp;
      
      sprintf(szTmp,"%s/folded_ch%d.txt",gOutDir.c_str(),gFreqChannel);
      gOutFoldedMaxVsFreq=szTmp;      
   }
   
   if( gIntegrationTimeSec > 0 ){
      // integration time overwrittes option -a 
      double n_integr = (gIntegrationTimeSec)/(gOutFineChannels*gTimesampleInSec);
      printf("INFO : integration time specified to be %.6f [sec] -> setting number of averaged integrations to %d (exact %.6f)\n",gIntegrationTimeSec,int(n_integr),n_integr);
      gAvgNSpectraToFITS = int(n_integr);
   }
   
   if( gOutPutType == eU || gOutPutType == eV || gOutPutType ==  ePhaseXY || gOutPutType == eXY )
   {
      gShiftReIm = true;
   }
}                                                                           

void printf_parameters()
{
  printf("#####################################\n");
  printf("PARAMETERS :\n");
  printf("#####################################\n");
  printf("Input file      = %s\n",filename);
  if( strlen(gOutFile) ){
     printf("Output file        = %s\n",gOutFile);
  }else{
     printf("Output file        = -\n");
  }  
  printf("Output type         = %d\n",gOutPutType);
  printf("Output folded FITS  = %s\n",gOutFoldedFits.c_str());
  printf("Output directory    = %s\n",gOutDir.c_str());
  printf("Save filterbank     = %d\n",gSaveFilterbankFiles);
  printf("Max time to process = %.6f [sec]\n",gMaxTimeToProcessInSec);
  printf("Dump channeld idx  = %d\n",gDumpChannel);
  printf("Start Index        = %d (save to header file -s option)\n",gStartIndex);
  printf("Skip N bytes       = %d (no save to header file -X option)\n",gSkipNBytes);
  printf("# of timesteps to process = %d\n",gProcessNSamples);
  printf("# averaged points  = %d\n",gAvgNSpectraToFITS);
  printf("Multiplier         = %d\n",gMultiplier);
  printf("Verb level         = %d\n",gVerb);
  printf("Check zeros        = %d\n",gCheckZeros);
  printf("Number of channels = %d\n",gNChannels);
  printf("Save number of channels = %d\n",gSaveNChannels);
  printf("Correct with inverse Jones = %.4f\n",gCorrWithInvJones);
  printf("Output FITS file base = %s\n",gOutFITSBase.c_str());
  printf("\tTransposed FITS = %d\n",gTransposedFITS);
  printf("Pulsar folding period = %.6f [sec]\n",gFoldingPulsarPeriod);
  printf("Pulsar folding phase bins number = %d\n",gNFoldingBinsCount);
  printf("Divide by mean power spectrum = %d\n",gBandpassCalMethod);
  printf("Add XY phase = %d\n",gAddXYPhase);
  printf("#####################################\n");       
}

void init_defaults()
{
   gOutDir = "./";
}

// channel_0_16_1646062557.497930.dada
double parse_filename( const char* filename )
{
   double uxtime=0.00;
   int ch,n_ch;
   if( sscanf( filename, "channel_%d_%d_%lf.dada", &ch, &n_ch, &uxtime) == 3 ){
      return uxtime;
   }
   if( sscanf( filename, "channel_%d_%d_%lf.dat", &ch, &n_ch, &uxtime) == 3 ){
      return uxtime;
   }
   
   printf("WARNING : unknown format of input datafile %s\n",filename);   
   
   return 0.00;
}

int main(int argc,char* argv[])
{
  init_defaults();
  if( argc<2 || (argc==2 && strncmp(argv[1],"-h",2)==0) ){
    usage();
  }

  if( argc>=2 ){
    strcpy(filename,argv[1]);
  }
  
  // parsing and printing paramters :
  parse_cmdline(argc-1,argv+1);
  printf_parameters();

  n_samples=gNBytesPerSample*gNChannels*1024;  
  FILE *f=NULL;  
  FILE* outf=NULL;
  char* buffer = new char[n_samples];
  int line_size = n_samples;
  int avg_count=0;
  int n = n_samples;

  
  if( gVerb ){
    printf("Reading file %s\n",filename);fflush(stdout);
  }
  f = fopen(filename, "rb");  
  if( strlen(gOutFileAvg) ){
     outf = fopen(gOutFileAvg,"w");
  }
  
  long int n_filesize_bytes = GetFileSize( filename );
  long int n_file_timesamples = n_filesize_bytes / gNChannels / gNPols / 2; // the last 2 for RE/IM
  if( gMaxTimeToProcessInSec > 0 ){
     // double time_processed = n_spectra*(1.08/1e6)*gOutFineChannels;
     int new_file_timesamples = int( gMaxTimeToProcessInSec/gTimesampleInSec );
     printf("INFO : file samples = %ld overwritten with %d as only %.6f [sec] required to be processed\n",n_file_timesamples,new_file_timesamples,gMaxTimeToProcessInSec);
     n_file_timesamples = new_file_timesamples;
  }
  double file_duration_sec = double(n_file_timesamples)*gTimesampleInSec;
  printf("File size = %.0f (%ld) bytes -> %.0f timesamples -> %.6f [sec]\n",double(n_filesize_bytes),n_filesize_bytes,double(n_file_timesamples),file_duration_sec);
  

  double uxtime_start = parse_filename( filename );  
  printf("INFO : parsed input filename %s to get unixtime = %.8f\n",filename,uxtime_start);
  
  FILE* out_dat_f = fopen( gOutFileDat , "w" );
  FILE* out_totpower_f = fopen( gOutFileTotPower ,"w" );
  if( !out_totpower_f ){
     printf("ERROR : could not create file %s\n",gOutFileTotPower);
     exit(-1);
  }else{
     printf("INFO : file %s created ok\n",gOutFileTotPower);
  }

  
  if (f)
  {
      int i=0,line=0;
      int line_raw=0;
      long int total_processed_samples = 0;
      long int zeros_count = 0;
      long int minus128_count = 0;
      
      double accum = 0;
      int accum_count = 0;
      long int avg_index = 0;
      
      if( gSkipNBytes > 0 ){
         fseek( f , gSkipNBytes*1, SEEK_SET ); // was *BYTES_PER_SAMPLE
         printf("DEBUG : set file position to index = %d (option -X)\n",gSkipNBytes*1); // was *BYTES_PER_SAMPLE
      }else{      
         printf("INFO : skiping of bytes with fseek (-X option) is not required, checking if -s to skip and save (.dada header is required)\n");
         // add fseek to find some specific pulses at a given index :                  
         if ( gStartIndex > 0 ){
              printf("INFO skip and save (option -c) is required (gStartIndex = %d bytes)\n",gStartIndex);
              char* buffer_tmp = new char[gStartIndex+10];
              int n_tmp = fread(buffer_tmp, 1, gStartIndex, f);
              if ( n_tmp != gStartIndex ){
                 printf("ERROR : could not read %d bytes of header only read %d\n",gStartIndex,n_tmp);
              }else{
                 printf("INFO : read %d bytes, same as expected %d bytes\n",n_tmp,gStartIndex);
              }
              FILE* out_header_f = fopen("header_saved.txt","wb");
//              fprintf(out_header_f,"%s",buffer_tmp);
              int fd = fileno(out_header_f);
              size_t n_written = write(fd, buffer_tmp, n_tmp);
              fclose(out_header_f);
              delete [] buffer_tmp;
         }else{
            printf("INFO skip and save (option -c) is not required (gStartIndex = %d bytes)\n",gStartIndex);
         }
      }
      printf("PROGRESS : starting to read %d bytes in each iteration\n",gNBytesPerSample*line_size);fflush(stdout);
      
//      char* out_buffer = new char[line_size/gNChannels];
  
      int n_fine_fft_samples = gOutFineChannels;    
      int counter=0;
      std::complex<float>* in_x = new std::complex<float>[gOutFineChannels];
      std::complex<float>* in_y = new std::complex<float>[gOutFineChannels];
      double* out_spectrum_i = new double[gOutFineChannels];      
      double* out_spectrum_x = new double[gOutFineChannels];      
      std::complex<float>* out_spectrum_reim_x = new std::complex<float>[gOutFineChannels];
      std::vector<double> out_spectrum_x_shifted, out_spectrum_y_shifted, out_xy_phase, out_stokes_q, out_stokes_u, out_stokes_v, out_stokes_i;
//      std::vector< std::complex<float> > out_stokes_u, out_stokes_v;
      std::vector< std::complex<float> > out_spectrum_reim_x_shifted, out_spectrum_reim_y_shifted;
      int out_count_x=0;
      double* out_spectrum_y = new double[gOutFineChannels];      
      std::complex<float>* out_spectrum_reim_y = new std::complex<float>[gOutFineChannels];
      int out_count_y=0;
      
      std::vector<double> avg_spectrum2fits_x, avg_spectrum2fits_y, avg_spectrum2fits_i, avg_spectrum2fits_q, avg_xyphase, avg_spectrum2fits_u, avg_spectrum2fits_v; 
//      std::vector< std::complex<double> > avg_spectrum2fits_u, avg_spectrum2fits_v;
      out_xy_phase.assign(gOutFineChannels,0);
      out_stokes_q.assign(gOutFineChannels,0);
      out_stokes_u.assign(gOutFineChannels,0);
      out_stokes_v.assign(gOutFineChannels,0);
      out_stokes_i.assign(gOutFineChannels,0);
      avg_spectrum2fits_x.assign(gOutFineChannels,0);
      avg_spectrum2fits_y.assign(gOutFineChannels,0);
      avg_spectrum2fits_i.assign(gOutFineChannels,0);
      avg_spectrum2fits_q.assign(gOutFineChannels,0);
      avg_spectrum2fits_u.assign(gOutFineChannels,0);
      avg_spectrum2fits_v.assign(gOutFineChannels,0);
      avg_xyphase.assign(gOutFineChannels,0);
      int n_avg_spectrum2fits = 0;
      
      // output variables for:
      std::vector<double> mean_power_x;
      std::vector<double> mean_power_y;
      mean_power_x.assign(gOutFineChannels,0.00);
      mean_power_y.assign(gOutFineChannels,0.00);
      double norm = 1.0/sqrt(gOutFineChannels);
      int n_spectra=0;      
      int timesample=0;
      
      CBgFits* pOutFitsX = NULL;
      CBgFits* pOutFitsY = NULL;
      CBgFits* pOutFitsXYPhase = NULL;
      CBgFits* pOutFitsU = NULL;
      CBgFits* pOutFitsV = NULL;
      CBgFits* pFoldedI = NULL;
      CBgFits* pFoldedQ = NULL;
      CBgFits* pFoldedU = NULL;
      CBgFits* pFoldedV = NULL;
      CBgFits* pFoldedCounter = NULL;
      SigprocFile* pOutFilterbankI = NULL;
      float* filterbank_buffer=NULL;
      int out_fits_line_index=0;
      int nTimesteps = ((n_file_timesamples / gOutFineChannels)/gAvgNSpectraToFITS); // ignore partial 
      printf("INFO : calculated %d timesteps to be saved as averaged spectra\n",nTimesteps);
      if( strlen(gOutFITSBase.c_str()) > 0 ){
         if( gTransposedFITS ){
            pOutFitsX = new CBgFits( nTimesteps, gOutFineChannels );
            pOutFitsY = new CBgFits( nTimesteps, gOutFineChannels );
            pOutFitsU = new CBgFits( nTimesteps, gOutFineChannels );
            pOutFitsV = new CBgFits( nTimesteps, gOutFineChannels );
            pOutFitsXYPhase = new CBgFits( nTimesteps, gOutFineChannels );
         }else{
            pOutFitsX = new CBgFits( gOutFineChannels , nTimesteps );
            pOutFitsY = new CBgFits( gOutFineChannels , nTimesteps );
            pOutFitsU = new CBgFits( gOutFineChannels , nTimesteps );
            pOutFitsV = new CBgFits( gOutFineChannels , nTimesteps );
            pOutFitsXYPhase = new CBgFits( gOutFineChannels , nTimesteps );
         }
         pFoldedI = new CBgFits( gNFoldingBinsCount, gOutFineChannels );
         pFoldedI->SetValue(0.00);
         pFoldedQ = new CBgFits( gNFoldingBinsCount, gOutFineChannels );
         pFoldedQ->SetValue(0.00);
         pFoldedU = new CBgFits( gNFoldingBinsCount, gOutFineChannels );
         pFoldedU->SetValue(0.00);
         pFoldedV = new CBgFits( gNFoldingBinsCount, gOutFineChannels );
         pFoldedV->SetValue(0.00);
         pFoldedCounter = new CBgFits( gNFoldingBinsCount, gOutFineChannels );
         pFoldedCounter->SetValue(0.00);
         
         if( gSaveFilterbankFiles ){
            double freq_higher_end = gFreqChannel*gChannel2FreqMultiplier + gFullBW/2.00;
            double fine_channel = gFullBW/gOutFineChannels;
            double foff = -fine_channel;
            double tstart = uxtime_start;
            double tsamp  = (gTimesampleInSec*gOutFineChannels*gAvgNSpectraToFITS); // in seconds
            int nbits = 8*int(sizeof(float));
            pOutFilterbankI = new SigprocFile( nbits, 1, gOutFineChannels, freq_higher_end, foff, tstart, tsamp );
            printf("INFO : creating filterbank filewith the following parameters : freq_higher_end = %.6f MHz, fine_channel_bw = %.6f MHz, foff = %.6f MHz, tstart = %.6f ux, tsamp = %.6f sec\n",freq_higher_end,fine_channel,foff,tstart,tsamp);
            filterbank_buffer = new float[gOutFineChannels];
            
            sprintf(gOutFilFile,"%s/%s_avg%d_i.fil",gOutDir.c_str(),gOutFITSBase.c_str(),gAvgNSpectraToFITS);
            printf("DEBUG : writting header ...\n");
            pOutFilterbankI->name( gOutFilFile );
            pOutFilterbankI->sourcename("J0835-4510");
            pOutFilterbankI->telescope_id(30);// MWA
            pOutFilterbankI->src_raj(83520.61149);
            pOutFilterbankI->src_dej(-451034.8751);
            pOutFilterbankI->FillHeader( true, false );
            int ret = pOutFilterbankI->WriteHeader( gOutFilFile , false, true );
            printf("DEBUG : header written to fil file %s OK ( %d bytes written )\n",gOutFilFile,ret);
         }         
         
         // init keywords :
         pOutFitsX->SetKeywordFloat( "RA_deg", gRA_deg );
         pOutFitsX->SetKeywordFloat( "DEC_deg", gDEC_deg );
         pOutFitsY->SetKeywordFloat( "RA_deg", gRA_deg );
         pOutFitsY->SetKeywordFloat( "DEC_deg", gDEC_deg );
         pOutFitsU->SetKeywordFloat( "RA_deg", gRA_deg );
         pOutFitsU->SetKeywordFloat( "DEC_deg", gDEC_deg );
         pOutFitsU->SetKeywordFloat( "RA_deg", gRA_deg );
         pOutFitsV->SetKeywordFloat( "DEC_deg", gDEC_deg );
         
         pFoldedI->SetKeywordFloat( "RA_deg", gRA_deg );
         pFoldedI->SetKeywordFloat( "DEC_deg", gDEC_deg );
         pFoldedQ->SetKeywordFloat( "RA_deg", gRA_deg );
         pFoldedQ->SetKeywordFloat( "DEC_deg", gDEC_deg );
         pFoldedU->SetKeywordFloat( "RA_deg", gRA_deg );
         pFoldedU->SetKeywordFloat( "DEC_deg", gDEC_deg );
         pFoldedV->SetKeywordFloat( "RA_deg", gRA_deg );
         pFoldedV->SetKeywordFloat( "DEC_deg", gDEC_deg );
         
         double memory_size_gb = (gOutFineChannels*nTimesteps)*sizeof(float)/1e9;
         printf("INFO : allocated memory for dynamic spectra for X and Y pols. of size (%d,%d) -> approx. %.4f [GB] each\n",gOutFineChannels , nTimesteps, memory_size_gb );
      }
      
      // total power :
      double total_power_x=0.00, total_power_y=0.00;
      long int total_power_count=0;
      
      // gNBytesPerSample = BYTES_PER_SAMPLE = 4        
      bool bContinue = true;
      while( (n = fread(buffer, 1, n_samples, f)) > 0 && bContinue ){ // reads a line of all fine channels 
      
        int bytes_written = 0;
        int all_zeros_count=0;
        for(int i=(gDumpChannel*gNBytesPerSample);i<line_size;i+=(gNBytesPerSample*gNChannels)){
           // ordering 1 : X Y | X Y | X Y | X Y ...
           //            CH1    CH2  CH1   CH2 ...
           char x_re = buffer[i];
           char x_im = buffer[i+1];
           char y_re = buffer[i+2];        
           char y_im = buffer[i+3];
           
           if( gVerb > 0 ){
              if( x_re == 0 && x_im == 0 && y_re == 0 && y_im == 0 ){
                 all_zeros_count++;
              }
           }
           
           if( gCheckZeros ){
              for(int ttt=0;ttt<4;ttt++){
                 if( buffer[i+ttt] == 0 ){ 
                    zeros_count++;          
                 }
                 if( buffer[i+ttt] == -128 ){ 
                    minus128_count++;          
                 }
              }
              total_processed_samples += 4;
           }
           
           if( gAddXYPhase && 0 ){
              double x = gFreqChannel*(400.00/512.00);// freq in MHz 
//              double xy_phase = (3.8225454867455517238283846381818*(x*x*x*(-0.00000801329232402725466499759)+x*x*(0.00249363049704243167070)+x-323.5+3.1415/2.00+6+3.13792111655815419)); //  + 380*3.141592654);
              double diff = (x-323.75950000);
              double xy_phase = 2.00*M_PI*( diff*0.10675035 + diff*diff*(-0.00649167) + diff*diff*diff*0.00001584442515310631 -0.00010545467580092109 ) + M_PI/2.00;

              std::complex<float> tmp( y_re, y_im );
              double arg = std::arg( tmp ) + xy_phase;
              double r = sqrt(y_re*y_re + y_im*y_im);
              std::complex<float> tmp2( r*cos(arg) , r*sin(arg) );
              
              y_re = tmp2.real();
              y_im = tmp2.imag();

              // Perhaps better to set Beam_x = 1 -> Beam_y = (I-Q)/(I+Q) and divide y_re and y_im by SQRT(Beam_y)              
              // x_re = x_re / sqrt(Beam_x) : add division by sqrt(Beam_x) - as fitted to (I+Q)/)I-Q)
              // y_re = y_re / sqrt(Beam_x) : add division by sqrt(Beam_x) - as fitted to (I+Q)/)I-Q)
           }
           
           // total power X and Y for a given integration (n_fine_fft_samples)
           total_power_x += x_re*x_re + x_im*x_im;
           total_power_y += y_re*y_re + y_im*y_im;
           total_power_count += 1;
           
           in_x[counter] = std::complex<float>( x_re, x_im );
           in_y[counter] = std::complex<float>( y_re, y_im );
           counter++;
           
           if( counter == n_fine_fft_samples )
           {
              // do FFT , accumlate 
              // int CSpectrometer::doFFT( std::complex<float>* in, int in_count, double* spectrum, std::complex<float>* spectrum_reim, int& out_count, double norm )
              // norm = 1.0/sqrt(context->n_chan)
              CSpectrometer::doFFT( in_x, n_fine_fft_samples, out_spectrum_x, out_spectrum_reim_x, out_count_x, norm );
              CSpectrometer::doFFT( in_y, n_fine_fft_samples, out_spectrum_y, out_spectrum_reim_y, out_count_y, norm );
              
              // do fft shift !!! to have DC in bin 0 !
              // TODO :
              if( pOutFitsX && pOutFitsY ){
                 // only required if FITS has to be saved :
                 CSpectrometer::fft_shift( out_spectrum_x, out_count_x, out_spectrum_x_shifted );
                 CSpectrometer::fft_shift( out_spectrum_y, out_count_y, out_spectrum_y_shifted );
                 
                 
                 if( gShiftReIm ){
                    CSpectrometer::fft_shift( out_spectrum_reim_x, out_count_x, out_spectrum_reim_x_shifted );
                    CSpectrometer::fft_shift( out_spectrum_reim_y, out_count_y, out_spectrum_reim_y_shifted );
                    
                    if( gAddXYPhase ){
                       for(int jjj=0;jjj<out_spectrum_reim_x_shifted.size();jjj++){
                          // TODO 
                          double freq_fine_mhz = GetFineChanFreqMHz(jjj);
                          double diff = (freq_fine_mhz-323.75950000);
                          double xy_phase = 2.00*M_PI*( diff*0.10675035 + diff*diff*(-0.00649167) + diff*diff*diff*0.00001584442515310631 -0.00010545467580092109 ) + M_PI/2.00;

                          // X pol :
                          std::complex<float> tmp( out_spectrum_reim_y_shifted[jjj].real(), out_spectrum_reim_y_shifted[jjj].imag() );
                          double arg = std::arg( tmp ) + xy_phase;
                          double r = sqrt( tmp.real()*tmp.real() + tmp.imag()*tmp.imag() );
                          std::complex<float> tmp2( r*cos(arg) , r*sin(arg) );
                          out_spectrum_reim_y_shifted[jjj] = tmp2;
                          out_spectrum_y_shifted[jjj] = tmp2.real()*tmp2.real() + tmp2.imag()*tmp2.imag();
                      }
                    }

                    
                    // TODO : calculate XY phase and put into : out_xy_phase
                    for(int fch=0;fch<out_count_x;fch++){
//                       std::complex<float> y_conj = (out_spectrum_reim_y_shifted[fch]).conj();
                       std::complex<float> xx = out_spectrum_reim_x_shifted[fch]*std::conj( out_spectrum_reim_x_shifted[fch] );
                       std::complex<float> yy = out_spectrum_reim_y_shifted[fch]*std::conj( out_spectrum_reim_y_shifted[fch] );
                       std::complex<float> xy = out_spectrum_reim_x_shifted[fch]*std::conj( out_spectrum_reim_y_shifted[fch] );
                       out_xy_phase[fch] = std::arg( xy ) *(180.00/M_PI); // phase in degree
                       
                       std::complex<float> yx = std::conj( out_spectrum_reim_x_shifted[fch] )*out_spectrum_reim_y_shifted[fch];
                       out_stokes_u[fch] = 0.5*(xy + yx).real(); // u = (xy+yx)/2 // was std::complex<float>(0.5,0.00)*(xy + yx)
                       out_stokes_v[fch] = 0.5*(xy - yx).imag(); // V = (xy-yx)/2i // was std::complex<float>(0.5,0.00)*(xy - yx)
                       out_stokes_q[fch] = xx.real() - yy.real();
                       out_stokes_i[fch] = xx.real() + yy.real();
                    }
                 }
                 
                 for(int fch=0;fch<out_count_x;fch++){ // was n_fine_fft_samples
                     avg_spectrum2fits_x[fch] += out_spectrum_x_shifted[fch];
                     avg_spectrum2fits_y[fch] += out_spectrum_y_shifted[fch];
                     avg_xyphase[fch] += out_xy_phase[fch];
                     avg_spectrum2fits_u[fch] += out_stokes_u[fch];
                     avg_spectrum2fits_v[fch] += out_stokes_v[fch];
                     avg_spectrum2fits_q[fch] += out_stokes_q[fch];
                 }
                 // center_time += ((double(n_avg_spectrum2fits+0.5)*double(n_fine_fft_samples)*1.08)/1000000.00);
                 n_avg_spectrum2fits++;
                 
                 if( n_avg_spectrum2fits == gAvgNSpectraToFITS){
                    double center_time = (out_fits_line_index+0.5)*(gAvgNSpectraToFITS*n_fine_fft_samples*gTimesampleInSec);
                    if( gVerb>=2 ){
                       printf("DEBUG : Dumping %d-th spectrum, centre time = %.8f [sec] ( = (%d+0.5)*%d*%d*%.8f\n",out_fits_line_index,center_time,out_fits_line_index,gAvgNSpectraToFITS,n_fine_fft_samples,gTimesampleInSec);
                    }
                    
                    
                    // average all (x,y,xyphase,u,v,) :
                    for(int fch=0;fch<avg_spectrum2fits_x.size();fch++){
                       avg_spectrum2fits_x[fch] = avg_spectrum2fits_x[fch] / n_avg_spectrum2fits;
                    }
                    for(int fch=0;fch<avg_spectrum2fits_y.size();fch++){
                       avg_spectrum2fits_y[fch] = avg_spectrum2fits_y[fch] / n_avg_spectrum2fits;
                    }
                    for(int fch=0;fch<avg_xyphase.size();fch++){
                       avg_xyphase[fch] = avg_xyphase[fch] / n_avg_spectrum2fits;
                       avg_spectrum2fits_u[fch] = avg_spectrum2fits_u[fch] / n_avg_spectrum2fits; // was std::complex<double>( 1.00/n_avg_spectrum2fits , 0.00) *
                       avg_spectrum2fits_v[fch] = avg_spectrum2fits_v[fch] / n_avg_spectrum2fits; // was std::complex<double>( 1.00/n_avg_spectrum2fits , 0.00) * 
                       avg_spectrum2fits_q[fch] = avg_spectrum2fits_q[fch] / n_avg_spectrum2fits;
                    }
                    
                    
                    // folding :
                    double period_index = double(center_time)/double(gFoldingPulsarPeriod);
                    int folding_bin = int((period_index-(int(period_index)))*gNFoldingBinsCount); // fractional number of periods times number of bins 
                    for(int fch=0;fch<avg_spectrum2fits_x.size();fch++){
                       avg_spectrum2fits_i[fch] = avg_spectrum2fits_x[fch] + avg_spectrum2fits_y[fch];
                       
                       int    corrent_counter = pFoldedCounter->getXY( folding_bin, fch );

                       // StokeS I folding :
                       double current_value = pFoldedI->getXY( folding_bin, fch );
                       pFoldedI->setXY( folding_bin, fch, (current_value + avg_spectrum2fits_i[fch]) );

                       // Stokes Q folding :
                       current_value = pFoldedQ->getXY( folding_bin, fch );
                       pFoldedQ->setXY( folding_bin, fch, (current_value + avg_spectrum2fits_q[fch]) );

                       // Stokes U folding :
                       current_value = pFoldedU->getXY( folding_bin, fch );
                       pFoldedU->setXY( folding_bin, fch, (current_value + avg_spectrum2fits_u[fch]) );

                       // Stokes V folding :
                       current_value = pFoldedV->getXY( folding_bin, fch );
                       pFoldedV->setXY( folding_bin, fch, (current_value + avg_spectrum2fits_v[fch]) );
                       
                       pFoldedCounter->setXY( folding_bin, fch, corrent_counter+1 );
                    }
                    
                 
                    // copy to output FITS 
                    bool isTimestepOK = false;
                    if( gTransposedFITS ){
                       if( out_fits_line_index < pOutFitsX->GetXSize() && out_fits_line_index < pOutFitsY->GetXSize() ){
                          isTimestepOK = true;
                       }
                    }else{
                       if( out_fits_line_index < pOutFitsX->GetYSize() && out_fits_line_index < pOutFitsY->GetYSize() ){
                          isTimestepOK = true;
                       }
                    }                   
                    
                    if( isTimestepOK ){
                       if( gTransposedFITS ){
                          for(int fch=0;fch<avg_spectrum2fits_x.size();fch++){
                             pOutFitsX->setXY(out_fits_line_index,fch,avg_spectrum2fits_x[fch]);
                             pOutFitsY->setXY(out_fits_line_index,fch,avg_spectrum2fits_y[fch]);
                             pOutFitsXYPhase->setXY(out_fits_line_index,fch,avg_xyphase[fch]);
                             pOutFitsU->setXY(out_fits_line_index,fch,avg_spectrum2fits_u[fch]); // TODO: should really save both RE/IM as a test to check if IM=0 ?
                             pOutFitsV->setXY(out_fits_line_index,fch,avg_spectrum2fits_v[fch]); // TODO: should really save both RE/IM as a test to check if RE=0 ? actually I did one save RE and it was 0 !!!
                          }
                       }else{
                          pOutFitsX->set_line( out_fits_line_index , avg_spectrum2fits_x );
                          pOutFitsY->set_line( out_fits_line_index , avg_spectrum2fits_y );
                          pOutFitsXYPhase->set_line( out_fits_line_index , avg_xyphase );
                          
                          
                          for(int fch=0;fch<avg_spectrum2fits_x.size();fch++){
                              pOutFitsU->setXY( fch, out_fits_line_index, avg_spectrum2fits_u[fch]); // TODO: should really save both RE/IM as a test to check if IM=0 ?
                              pOutFitsV->setXY( fch, out_fits_line_index, avg_spectrum2fits_v[fch]); // TODO: should really save both RE/IM as a test to check if RE=0 ? actually I did one save RE and it was 0 !!!
                          }
                          // pOutFitsU->set_line( out_fits_line_index , avg_spectrum2fits_u );
                          // pOutFitsV->set_line( out_fits_line_index , avg_spectrum2fits_v );
                       }
                       
                       // write to filterbank (if required)
                       if( pOutFilterbankI && filterbank_buffer ){
                          for(int ch=0;ch<avg_spectrum2fits_i.size();ch++){                           
                             filterbank_buffer[ch] = avg_spectrum2fits_i[ch];
                          }
                          pOutFilterbankI->WriteData( filterbank_buffer , avg_spectrum2fits_i.size() ); 
                          // printf("DEBUG : wrote %d channels to filterbank file\n",avg_spectrum2fits_i.size() );
                       }
                       
                       
                       if( (out_fits_line_index % 1000)==0 || nTimesteps < 1000 ){
                          printf("PROGRESS : averaged spectrum %d / %d\n",out_fits_line_index,nTimesteps);
                       }
                       
                       out_fits_line_index++;
                    }else{
                       printf("WARNING : timestep outside bounderies !!!\n");
                    }
                    
                    if( (out_fits_line_index % 100) == 0 && 0 ){
                       char szControlSpectrum[1024];
                       sprintf(szControlSpectrum,"spectrum_%08d.txt",out_fits_line_index);
                       
                       save_spectrum(szControlSpectrum,avg_spectrum2fits_x);
                    }
                 
                    double total_power_i = total_power_x+total_power_y;
                    fprintf(out_totpower_f,"%.8f %.8f %.8f %.8f\n",center_time,total_power_x/double(n_avg_spectrum2fits),total_power_y/double(n_avg_spectrum2fits),total_power_i/double(n_avg_spectrum2fits));
                    total_power_x=0.00;
                    total_power_y=0.00;
                    total_power_count=0;

                    avg_spectrum2fits_x.assign(gOutFineChannels,0.00);
                    avg_spectrum2fits_y.assign(gOutFineChannels,0.00);
                    avg_xyphase.assign(gOutFineChannels,0.00);
                    n_avg_spectrum2fits = 0;
                 }                 
              }
              
              
              for(int fch=0;fch<out_count_x;fch++){ // was n_fine_fft_samples
                  mean_power_x[fch] += out_spectrum_x[fch];
                  mean_power_y[fch] += out_spectrum_y[fch];
              }
              
              // clean :
              counter = 0;
              
              n_spectra++;
           }
           
 
           if( (n_spectra % 100)==0 ){
              double time_processed = n_spectra*(gTimesampleInSec)*gOutFineChannels;
              if( gVerb>=2 ){
                 printf("PROGRESS : generated %d spectra (processed interval = %.6f sec)\n",n_spectra,time_processed);
                 printf("DEBUG : out_count_x = %d , out_count_y = %d\n",out_count_x,out_count_y);
              }
              
              if( gMaxTimeToProcessInSec > 0  ){
                 if( time_processed >= gMaxTimeToProcessInSec ){
                    printf("Processed %.6f [seconds] as required (%.6f [sec])\n",time_processed,gMaxTimeToProcessInSec);
                    bContinue = false;
                    break;
                 }
              }
              
           }
           timesample++;
        }
        
        if( gVerb > 0 ){
           printf("read line = %d, all zero count = %d\n",line,all_zeros_count);
        }
        line++;
      }      
      
      if( gCheckZeros ){
         printf("################################## STATISTICS ##################################\n");
         printf("Number of zeros          = %ld ( %.3f %% )\n",zeros_count,(double(zeros_count)/double(total_processed_samples))*100.00);
         printf("Number of -128           = %ld ( %.3f %% )\n",minus128_count,(double(minus128_count)/double(total_processed_samples))*100.00);
         printf("Total number of samples  = %ld\n",total_processed_samples);
         printf("################################################################################\n");
      }
      
      vector<double> mean_x_shifted, mean_y_shifted;
      CSpectrometer::fft_shift( mean_power_x, mean_x_shifted );
      CSpectrometer::fft_shift( mean_power_y, mean_y_shifted );
      
      // write to output text file :
      for(int fch=0;fch<mean_x_shifted.size();fch++){
         fprintf(outf,"%d %.8f %.8f\n",fch,mean_x_shifted[fch]/n_spectra,mean_y_shifted[fch]/n_spectra);
      }
      printf("INFO : mean of %d spectra saved\n",n_spectra);
      
      if( pOutFitsX && pOutFitsY ){
         char szOutFits[1024];
         
         // X pol :
         sprintf(szOutFits,"%s/%s_avg%d_x.fits",gOutDir.c_str(),gOutFITSBase.c_str(),gAvgNSpectraToFITS);         
         printf("PROGRESS : saving output FITS file %s ...\n",szOutFits);fflush(stdout);
         
         double inttime = gTimesampleInSec*gOutFineChannels*gAvgNSpectraToFITS;
         double uxtime  = uxtime_start;
         double fine_channel_bw = gFullBW/gOutFineChannels;
         double freq_center = gFreqChannel * gChannel2FreqMultiplier;
         double freq_start = freq_center - gFullBW/2.00 + fine_channel_bw/2.00;
         double freq_header = freq_start;
         printf("INFO : saving FITS files, freq. start = %.6f MHz, freq. center = %.6f , bw = %.6f MHz , BW = %.6f MHz -> Header freq = %.6f MHz\n",freq_start,freq_center,fine_channel_bw,gFullBW,freq_header);
         if( gTransposedFITS ){
            pOutFitsX->PrepareBigHornsHeaderTransposed( uxtime, inttime, freq_header, fine_channel_bw );
         }else{
            pOutFitsX->PrepareBigHornsHeader( uxtime, inttime, freq_header, fine_channel_bw );
         }
         if( pOutFitsX->WriteFits( szOutFits ) ){
            printf("ERROR : when writting output FITS %s\n",szOutFits);
         }
         
         // Y pol :
         sprintf(szOutFits,"%s/%s_avg%d_y.fits",gOutDir.c_str(),gOutFITSBase.c_str(),gAvgNSpectraToFITS);         
         printf("PROGRESS : saving output FITS file %s ...\n",szOutFits);fflush(stdout);
         if( gTransposedFITS ){
            pOutFitsY->PrepareBigHornsHeaderTransposed( uxtime, inttime, freq_header, fine_channel_bw );
         }else{
            pOutFitsY->PrepareBigHornsHeader( uxtime, inttime, freq_header, fine_channel_bw );
         }
         if( pOutFitsY->WriteFits( szOutFits ) ){
            printf("ERROR : when writting output FITS %s\n",szOutFits);
         }
         
         // Stokes U :
         sprintf(szOutFits,"%s/%s_avg%d_u.fits",gOutDir.c_str(),gOutFITSBase.c_str(),gAvgNSpectraToFITS);         
         printf("PROGRESS : saving output FITS file %s ...\n",szOutFits);fflush(stdout);
         if( gTransposedFITS ){
            pOutFitsU->PrepareBigHornsHeaderTransposed( uxtime, inttime, freq_header, fine_channel_bw );
         }else{
            pOutFitsU->PrepareBigHornsHeader( uxtime, inttime, freq_header, fine_channel_bw );
         }
         if( pOutFitsU->WriteFits( szOutFits ) ){
            printf("ERROR : when writting output FITS %s\n",szOutFits);
         }

         // Stokes V :
         sprintf(szOutFits,"%s/%s_avg%d_v.fits",gOutDir.c_str(),gOutFITSBase.c_str(),gAvgNSpectraToFITS);         
         printf("PROGRESS : saving output FITS file %s ...\n",szOutFits);fflush(stdout);
         if( gTransposedFITS ){
            pOutFitsV->PrepareBigHornsHeaderTransposed( uxtime, inttime, freq_header, fine_channel_bw );
         }else{
            pOutFitsV->PrepareBigHornsHeader( uxtime, inttime, freq_header, fine_channel_bw );
         }
         if( pOutFitsV->WriteFits( szOutFits ) ){
            printf("ERROR : when writting output FITS %s\n",szOutFits);
         }
         
         // XY phase :
         sprintf(szOutFits,"%s/%s_avg%d_xyphase.fits",gOutDir.c_str(),gOutFITSBase.c_str(),gAvgNSpectraToFITS);         
         printf("PROGRESS : saving output FITS file %s ...\n",szOutFits);fflush(stdout);
         if( gTransposedFITS ){
            pOutFitsXYPhase->PrepareBigHornsHeaderTransposed( uxtime, inttime, freq_header, fine_channel_bw );
         }else{
            pOutFitsXYPhase->PrepareBigHornsHeader( uxtime, inttime, freq_header, fine_channel_bw );
         }
         if( pOutFitsXYPhase->WriteFits( szOutFits ) ){
            printf("ERROR : when writting output FITS %s\n",szOutFits);
         }

         // normalise if required :         
         normalise_folded( *pFoldedI );
         normalise_folded( *pFoldedQ );
         normalise_folded( *pFoldedU );
         normalise_folded( *pFoldedV );
         
         // save output FITS files:
         string szFoldedFits;
         szFoldedFits = gOutFoldedFits + "_I.fits";
         if( pFoldedI->WriteFits( szFoldedFits.c_str() ) ){
            printf("ERROR : when writting output FITS folded fits %s\n",szFoldedFits.c_str());
         }
         szFoldedFits = gOutFoldedFits + "_Q.fits";
         if( pFoldedQ->WriteFits( szFoldedFits.c_str() ) ){
            printf("ERROR : when writting output FITS folded %s\n",szFoldedFits.c_str());
         }
         szFoldedFits = gOutFoldedFits + "_U.fits";
         if( pFoldedU->WriteFits( szFoldedFits.c_str() ) ){
            printf("ERROR : when writting output FITS folded %s\n",szFoldedFits.c_str());
         }
         szFoldedFits = gOutFoldedFits + "_V.fits";
         if( pFoldedV->WriteFits( szFoldedFits.c_str() ) ){
            printf("ERROR : when writting output FITS folded %s\n",szFoldedFits.c_str());
         }
         
         // Stokes I : use the same memory buffer as for Stokes V :
         (*pOutFitsV) = (*pOutFitsX); // this is causing later core dump in delete pOutFitsV
         pOutFitsV->AddImages( *pOutFitsY );
         sprintf(szOutFits,"%s/%s_avg%d_i.fits",gOutDir.c_str(),gOutFITSBase.c_str(),gAvgNSpectraToFITS);
         printf("PROGRESS : saving output FITS file %s ...\n",szOutFits);fflush(stdout);
         if( gTransposedFITS ){
            pOutFitsV->PrepareBigHornsHeaderTransposed( uxtime, inttime, freq_header, fine_channel_bw );
         }else{
            pOutFitsV->PrepareBigHornsHeader( uxtime, inttime, freq_header, fine_channel_bw );
         }
         if( pOutFitsV->WriteFits( szOutFits ) ){
            printf("ERROR : when writting output FITS %s\n",szOutFits);
         }
         
         // Stokes Q : use the same memory buffer as for Stokes V :
         (*pOutFitsV) = (*pOutFitsX); // this is causing later core dump in delete pOutFitsV
         pOutFitsV->Subtract( *pOutFitsY );
         sprintf(szOutFits,"%s/%s_avg%d_q.fits",gOutDir.c_str(),gOutFITSBase.c_str(),gAvgNSpectraToFITS);
         printf("PROGRESS : saving output FITS file %s ...\n",szOutFits);fflush(stdout);
         if( gTransposedFITS ){
            pOutFitsV->PrepareBigHornsHeaderTransposed( uxtime, inttime, freq_header, fine_channel_bw );
         }else{
            pOutFitsV->PrepareBigHornsHeader( uxtime, inttime, freq_header, fine_channel_bw );
         }
         if( pOutFitsV->WriteFits( szOutFits ) ){
            printf("ERROR : when writting output FITS %s\n",szOutFits);
         }
         
         // save maximum vs. frequency :
         if( strlen(gOutFoldedMaxVsFreq.c_str())>0 ){
            double max_freq=0;
            FILE* out_max_f = fopen(gOutFoldedMaxVsFreq.c_str(),"w");
            for(int ch=0;ch<pFoldedI->GetYSize();ch++){
               double max_val_i = pFoldedI->GetMaxPower(ch,max_freq);               
               double freq_mhz = GetFineChanFreqMHz(ch);
               fprintf(out_max_f,"%.6f %.8f %.4f\n",freq_mhz,max_val_i,max_freq);
            }
            fclose(out_max_f);
         }
         
         delete pOutFitsX;
         delete pOutFitsY;
         delete pOutFitsU;
         delete pOutFitsV; // TODO : this is causing core dump, possibly due to operator= doing something weired !???
         delete pOutFitsXYPhase;
         
         delete pFoldedI;
         delete pFoldedQ;
         delete pFoldedU;
         delete pFoldedV;
         delete pFoldedCounter;
         
         if(pOutFilterbankI){
            pOutFilterbankI->Close();
            printf("INFO : output filterbank file closed\n");
            delete pOutFilterbankI;
         }         
      }
      
      if( in_x ){
        delete [] in_x;
      }
      if( in_y ){
         delete [] in_y;
      }
      if( out_spectrum_x ){
         delete [] out_spectrum_x;
      }
      if( out_spectrum_i ){
         delete [] out_spectrum_i;
      }
      if( out_spectrum_y ){
         delete [] out_spectrum_y;
      }
      if( out_spectrum_reim_x ){
         delete [] out_spectrum_reim_x;
      }
      if( out_spectrum_reim_y ){
         delete [] out_spectrum_reim_y;
      }
      if( filterbank_buffer ){
         delete [] filterbank_buffer;
      }
  }else{
     printf("ERROR : could not open file %s\n",filename);
  }

  fclose(f);
  
  if( out_totpower_f ){
     fclose( out_totpower_f );
  }
  
  if( out_dat_f ){
     fclose( out_dat_f );
  }
  
  if( outf ){
     fclose(outf);
  }  
}

