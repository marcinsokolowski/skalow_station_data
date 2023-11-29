#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>

// 8192 -> 8200 (for timestamp)
#define N_CHANNELS (24*32*4) // 10kHz channels 24*32*4 = 3072
#define MAX_N_SAMPLES 2*N_CHANNELS // Re/Im of 24*32 fine channels 
#define DATA_TYPE     short

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


// start value of the generator :
int gStartValue=1; // was 0, but 1 from 2022-02-18

#define BYTES_PER_SAMPLE 4 // this is re/im for X and re/im for Y -> 2 x 2 bytes = 4 bytes 

int gNBytesPerSample = BYTES_PER_SAMPLE;
int n_samples = 1024 * BYTES_PER_SAMPLE; // number of samples to read, used to be N_CHANNELS; 
int gStartIndex=-1;
int gProcessNSamples = -1;
int gMultiplier=1;
int gCheckZeros=0;

void usage()
{
   printf("test_pattern_generated_data DATA_FILE.dat -d -u -c -f outfile -i integration_index -v -a AVERAGE_N_POINTS -p POL -s START_INDEX -n PROCESS_TIMESTAMPS -M multiplier -Z -S SAVE_N_CHANNELS\n");
   printf("\n");
   printf("Options :\n");
//   printf("\t-u - assumes that UNIXTIME timestamp is added to line as extra record of DATA_TYPE structure\n");
   printf("\t-d - increases debug level\n");
   printf("\t-v - increases debug level\n");   
   printf("\t-c CHANNEL : specifies channel to be dumped [default 0]\n");      
   printf("\t-C N_CHANNELS : number of channels [default %d]\n",gNChannels);
   printf("\t-f OUTFILENAME : name of file to which data is dumped\n");
   printf("\t-i INTEGRATION_INDEX : index of integration (= adcv dump) to be dumped, by default all are dumped\n");   
   printf("\t-a AVERAGE_N_POINTS : average N data points [default %d], <=1 -> do not average\n",gAverageNSpectra);
   printf("\t-p POL : POL = 0 = X , POL = 1 = Y , POL = -1 X+Y (power)\n");
   printf("\t-s START_INDEX : where to set start index [default %d] , <0 -> starting from the start of the file\n",gStartIndex);
   printf("\t-S start expected value [default %d]\n",gStartValue);
   printf("\t-n PROCESS_TIMESTAMPS : number of timestamps to process [default %d and <0 is all]\n",gProcessNSamples);
   printf("\t-M multiplier : multiply timestep index (for example 128 fine channels * 7 averages to get ~1ms)\n");
   printf("\t-Z : check statistics of ZEROS\n");
   
   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "vudc:C:f:i:r:m:p:s:n:a:M:ZS:";
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
               gStartValue = atol( optarg );
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
            
         case 'Z':
            gCheckZeros = 1;
            break;
            
         default:
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();                                          
      }
   }
   
/*   if( strlen(filename) == 0 ){
      if( gStartIndex > 0 && gProcessNSamples > 0 ){
         sprintf(filename,"%ld_%ld.dat",gStartIndex,(gStartIndex+gProcessNSamples));
      }else{
         strcpy(filename,"test.dat");
      }
   }*/

   sprintf(gOutFileAvg,"%s_avg.txt",gOutFile);
   sprintf(gOutFileDat,"%s.dat",gOutFile);

   if( strlen(gOutFile) == 0 ){
      if( gStartIndex > 0 && gProcessNSamples > 0 ){
         sprintf(gOutFile,"%d_%d.txt",gStartIndex,(gStartIndex+gProcessNSamples));
      }else{
         strcpy(gOutFile,"test.txt");
      }
      
      printf("DEBUG : output file name was empty -> set to %s\n",gOutFile);
   }   
   
   gStartIndex = gStartIndex * gMultiplier;
   gProcessNSamples = gProcessNSamples * gMultiplier;
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
  printf("Dump channeld idx  = %d\n",gDumpChannel);
  printf("Start Index        = %d\n",gStartIndex);
  printf("# of timesteps to process = %d\n",gProcessNSamples);
  printf("# averaged points  = %d\n",gAverageNSpectra);
  printf("Multiplier         = %d\n",gMultiplier);
  printf("Verb level         = %d\n",gVerb);
  printf("Check zeros        = %d\n",gCheckZeros);
  printf("Number of channels = %d\n",gNChannels);
  printf("Start value        = %d\n",gStartValue);
  printf("#####################################\n");       
}

void DumpBuffer( char* data, int len )
{
   printf("\t\tBuffer :");
   for(int i=0;i<len;i++){
      printf("%d ",data[i]);
   }
   printf("\n\n");   
}

int main(int argc,char* argv[])
{
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
  
  FILE* out_dat_f = fopen( gOutFileDat , "w" );
  
  if (f)
  {
      int i=0,line=0;
      int line_raw=0;
      long int total_processed_samples = 0;
      long int zeros_count = 0;
      
      double accum = 0;
      int accum_count = 0;
      long int avg_index = 0;
      
      // add fseek to find some specific pulses at a given index :
//      if ( gStartIndex > 0 ){
//         fseek( f , gStartIndex*1, SEEK_SET ); // was *BYTES_PER_SAMPLE
//         printf("DEBUG : set file position to index = %d\n",gStartIndex*1); // was *BYTES_PER_SAMPLE
//      }
      
      printf("PROGRESS : starting to read %d bytes in each iteration\n",gNBytesPerSample*line_size);fflush(stdout);
      
      // gNBytesPerSample = BYTES_PER_SAMPLE = 4        
      bool bContinue = true;
      int n_bytes = gNChannels*4;
      int n_errors = 0;
      while( (n = fread(buffer, 1, n_bytes, f)) > 0 && bContinue ){ // reads a line of all fine channels 
        if( n != n_bytes ){
           printf("ERROR : could not read %d bytes from file read only %d\n",n_bytes,n);
        }
      
        int timesample=0;
        
        if( total_processed_samples <= 0 ){
           printf("DEBUG : %d %d %d %d\n",buffer[0],buffer[1],buffer[2],buffer[3]);
        }
        
                
        int expected_value = gStartValue;        
        int byte=0;
        while( byte < n ){
           int value = buffer[byte];
           if( gVerb >= 2 ){ 
              printf("DEBUG : %d : %d vs. %d\n",byte,value,expected_value);
           }

           if( value == 0 ){
              zeros_count++;
           }
           
           if( value != expected_value ){
              printf("ERROR : wrong value at byte = %ld , expected value = %d , actual value = %d\n",total_processed_samples,expected_value,value);
              DumpBuffer( buffer, n );
              n_errors++;
           }
           
           byte++;
           total_processed_samples++;
           expected_value++;
        }
        
        if( gVerb > 0 ){
           printf("read line = %d\n",line);
        }
        line++;
      }      
      
      printf("################################## STATISTICS ##################################\n");
      printf("Number of zeros          = %ld ( %.3f %% )\n",zeros_count,(double(zeros_count)/double(total_processed_samples))*100.00);
      printf("Total number of samples  = %ld\n",total_processed_samples);
      printf("Total number of errors   = %d\n",n_errors);
      printf("################################################################################\n");
  }else{
     printf("ERROR : could not open file %s\n",filename);
  }

  fclose(f);  
}

