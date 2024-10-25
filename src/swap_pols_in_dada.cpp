#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <complex>

// 8192 -> 8200 (for timestamp)
#define N_CHANNELS (24*32*4) // 10kHz channels 24*32*4 = 3072
#define MAX_N_SAMPLES 2*N_CHANNELS // Re/Im of 24*32 fine channels 
#define DATA_TYPE     short
#define DADA_HEADER_SIZE 4096

int gVerb=0;
int gSnapshotIdx=-1;
int gDumpChannel=0;
char filename[1024];
char gOutFile[1024],gOutFileAvg[1024],gOutFileDat[1024];
int gAddIndex=0;
int gAverageNSpectra=0;
int gPol=1; // 0 - X , 1 - Y 
int gNChannels=1;
int gSaveNChannels=1;

#define BYTES_PER_SAMPLE 4 // this is re/im for X and re/im for Y -> 2 x 2 bytes = 4 bytes 

int gNBytesPerSample = BYTES_PER_SAMPLE;
int n_samples = 1024 * BYTES_PER_SAMPLE; // number of samples to read, used to be N_CHANNELS; 
int gStartIndex=-1;
int gSkipNBytes=-1;
int gProcessNSamples = -1;
int gMultiplier=1;
int gCheckZeros=0;

double gCorrWithInvJones=-1; // normalisation factor
bool gXYPhase = false;

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


void usage()
{
   printf("swap_pols_in_dada.c in.dada out.dada\n");
   printf("\n");
   printf("Options :\n");
//   printf("\t-u - assumes that UNIXTIME timestamp is added to line as extra record of DATA_TYPE structure\n");
/*   printf("\t-d - increases debug level\n");
   printf("\t-v - increases debug level\n");   
   printf("\t-c CHANNEL : specifies channel to be dumped [default 0]\n");      
   printf("\t-C N_CHANNELS : number of channels [default %d]\n",gNChannels);
   printf("\t-f OUTFILENAME : name of file to which data is dumped\n");
   printf("\t-i INTEGRATION_INDEX : index of integration (= adcv dump) to be dumped, by default all are dumped\n");   
   printf("\t-a AVERAGE_N_POINTS : average N data points [default %d], <=1 -> do not average\n",gAverageNSpectra);
   printf("\t-p POL : POL = 0 = X , POL = 1 = Y , POL = -1 X+Y (power)\n");
   printf("\t-s START_INDEX : where to set start index [default %d] , <0 -> starting from the start of the file (and save these many bytes to file header_saved.txt\n",gStartIndex);   
   printf("\t-S SAVE_N_CHANNELS : how many channels to save into the output file [default %d]\n",gSaveNChannels);
   printf("\t-n PROCESS_TIMESTAMPS : number of timestamps to process [default %d and <0 is all]\n",gProcessNSamples);
   printf("\t-M multiplier : multiply timestep index (for example 128 fine channels * 7 averages to get ~1ms)\n");
   printf("\t-Z : check statistics of ZEROS\n");
   printf("\t-X N_BYTES : skip these many bytes [default %d]\n",gSkipNBytes);
   printf("\t-I : if >0 correct with inverse Jones and multiply by this factor [default %.4f]\n",gCorrWithInvJones);
   printf("\t-Y : enable XY phase correction\n");*/
   
   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "vudc:C:f:i:r:m:p:s:n:a:M:ZS:X:I:Y";
   int opt,opt_param,i;
   
//   strcpy(filename,"");
   strcpy(gOutFile,"out_file");
      
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
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

         case 'i':
            if( optarg ){   
               gSnapshotIdx = atol(optarg);
            }
            break;

         case 'I':
            gCorrWithInvJones = atof(optarg);
            break;

         case 'a':
            if( optarg ){   
               gAverageNSpectra = atol(optarg);
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
            
         case 'v':
            gVerb++;
            break;
            
         case 'X':
            if( optarg ){
               gSkipNBytes = atol( optarg );
            }
            break;
            
         case 'Z':
            gCheckZeros = 1;
            break;

         case 'Y':
            gXYPhase = true;
            break;
            
         default:
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();                                          
      }
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
/*  printf("Dump channeld idx  = %d\n",gDumpChannel);
  printf("Start Index        = %d (save to header file -s option)\n",gStartIndex);
  printf("Skip N bytes       = %d (no save to header file -X option)\n",gSkipNBytes);
  printf("# of timesteps to process = %d\n",gProcessNSamples);
  printf("# averaged points  = %d\n",gAverageNSpectra);
  printf("Multiplier         = %d\n",gMultiplier);
  printf("Verb level         = %d\n",gVerb);
  printf("Check zeros        = %d\n",gCheckZeros);
  printf("Number of channels = %d\n",gNChannels);
  printf("Save number of channels = %d\n",gSaveNChannels);
  printf("Correct with inverse Jones = %.4f\n",gCorrWithInvJones);
  printf("XY phase correction = %d\n",gXYPhase);*/
  printf("#####################################\n");       
}

int main(int argc,char* argv[])
{
  if( argc<2 || (argc==2 && strncmp(argv[1],"-h",2)==0) ){
    usage();
  }

  if( argc>=2 ){
    strcpy(filename,argv[1]);
  }
  if( argc>=3 ){
    strcpy(gOutFile,argv[2]);
  }
  
  // parsing and printing paramters :
//  parse_cmdline(argc-2,argv+2);
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
  outf = fopen(gOutFile,"w");
  
  if (f)
  {
      int i=0,line=0;
      int line_raw=0;
      long int total_processed_samples = 0;
      long int zeros_count = 0;
      
      double accum = 0;
      int accum_count = 0;
      long int avg_index = 0;
      
      char dada_header[DADA_HEADER_SIZE];
      if( fread(buffer, 1, DADA_HEADER_SIZE, f) != DADA_HEADER_SIZE ){
         printf("ERROR : could not read .dada header\n");
         exit(-1);
      }
      if( fwrite( buffer, DADA_HEADER_SIZE, 1, outf ) != DADA_HEADER_SIZE ){
         printf("ERROR : could not write .dada header\n");
         exit(-1);
      }

      printf("PROGRESS : starting to read %d bytes in each iteration\n",gNBytesPerSample*line_size);fflush(stdout);
      
      char* out_buffer = new char[n_samples];

      // gNBytesPerSample = BYTES_PER_SAMPLE = 4        
      bool bContinue = true;
      while( (n = fread(buffer, 1, n_samples, f)) > 0 && bContinue ){ // reads a line of all fine channels 
      
        int timesample=0;
        int bytes_written = 0;
        for(int i=(gDumpChannel*gNBytesPerSample);i<line_size;i+=(gNBytesPerSample*gNChannels)){
           // ordering 1 : X Y | X Y | X Y | X Y ...
           //            CH1    CH2  CH1   CH2 ...
           char x_re = buffer[i];
           char x_im = buffer[i+1];
           char y_re = buffer[i+2];        
           char y_im = buffer[i+3];
           
           
           out_buffer[i] = y_re;
           out_buffer[i+1] = y_im;
           out_buffer[i+2] = x_re;
           out_buffer[i+3] = x_im;
           
           timesample++;
        }
        
        if( fwrite( out_buffer, n_samples, 1, outf ) != n_samples ){
           printf("ERROR : when writting output .dada file\n");
           exit(-1);
        }
        if( gVerb > 0 ){
           printf("read line = %d\n",line);
        }
        line++;
      }      
      
  }else{
     printf("ERROR : could not open file %s\n",filename);
  }

  fclose(f);
  
  if( outf ){
     fclose(outf);
  }
}

