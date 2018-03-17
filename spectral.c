#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <complex.h>
#include <sndfile.h>
#include <fftw3.h>

#include <math.h>

#include "gnuplot_i.h"

/* taille de la fenetre */
// #define	FRAME_SIZE 4096
#define	FRAME_SIZE 1764 // 40ms à 44100 Hz
/* avancement */
// #define HOP_SIZE 1024
#define HOP_SIZE 1323 // 30ms à 44100 Hz

static gnuplot_ctrl *h;


void dft (double s[FRAME_SIZE], complex S[FRAME_SIZE]){
    for (size_t m = 0; m < FRAME_SIZE; m++) {
        S[m] = 0;
        for (size_t n = 0; n < FRAME_SIZE; n++) {
            S[m] += s[n] * cexp(-I * 2 * M_PI / FRAME_SIZE * n * m);
        }
    }
}



static void
print_usage (char *progname)
{
    printf ("\nUsage : %s <input file> \n", progname) ;
    puts ("\n");

}
static void
fill_buffer(double *buffer, double *new_buffer)
{
    int i;
    double tmp[FRAME_SIZE-HOP_SIZE];

    /* save */
    for (i=0;i<FRAME_SIZE-HOP_SIZE;i++)
    tmp[i] = buffer[i+HOP_SIZE];

    /* save offset */
    for (i=0;i<(FRAME_SIZE-HOP_SIZE);i++)
    {
        buffer[i] = tmp[i];
    }

    for (i=0;i<HOP_SIZE;i++)
    {
        buffer[FRAME_SIZE-HOP_SIZE+i] = new_buffer[i];
    }
}

static int
read_n_samples (SNDFILE * infile, double * buffer, int channels, int n)
{

    if (channels == 1)
    {
        /* MONO */
        int readcount ;

        readcount = sf_readf_double (infile, buffer, n);

        return readcount==n;
    }
    else if (channels == 2)
    {
        /* STEREO */
        double buf [2 * n] ;
        int readcount, k ;
        readcount = sf_readf_double (infile, buf, n);
        for (k = 0 ; k < readcount ; k++)
        buffer[k] = (buf [k * 2]+buf [k * 2+1])/2.0 ;

        return readcount==n;
    }
    else
    {
        /* FORMAT ERROR */
        printf ("Channel format error.\n");
    }

    return 0;
}

static int
read_samples (SNDFILE * infile, double * buffer, int channels)
{
    return read_n_samples (infile, buffer, channels, HOP_SIZE);
}

int
main (int argc, char * argv [])
{	char 		*progname, *infilename;
    SNDFILE	 	*infile = NULL ;
    SF_INFO	 	sfinfo ;

    progname = strrchr (argv [0], '/') ;
    progname = progname ? progname + 1 : argv [0] ;

    if (argc != 2)
    {	print_usage (progname) ;
        return 1 ;
    } ;

    infilename = argv [1] ;

    if ((infile = sf_open (infilename, SFM_READ, &sfinfo)) == NULL)
    {	printf ("Not able to open input file %s.\n", infilename) ;
    puts (sf_strerror (NULL)) ;
    return 1 ;
} ;

/* Read WAV */
int nb_frames = 0;
double new_buffer[HOP_SIZE];
double buffer[FRAME_SIZE];

/* Plot Init */
h = gnuplot_init();
gnuplot_setstyle(h, "lines");

int i;
for (i=0;i<(FRAME_SIZE/HOP_SIZE-1);i++)
{
    if (read_samples (infile, new_buffer, sfinfo.channels)==1)
    fill_buffer(buffer, new_buffer);
    else
    {
        printf("not enough samples !!\n");
        return 1;
    }
}

/* Info file */
printf("sample rate %d\n", sfinfo.samplerate);
printf("channels %d\n", sfinfo.channels);
printf("size %d\n", (int)sfinfo.frames);

static fftw_plan plan;
fftw_complex data_in[sfinfo.samplerate];
fftw_complex data_out[sfinfo.samplerate];
plan = fftw_plan_dft_1d(sfinfo.samplerate, data_in, data_out, FFTW_FORWARD, FFTW_ESTIMATE);




// fftw_real s[FRAME_SIZE]; /* domaine temporel */



for (int i = FRAME_SIZE; i < sfinfo.samplerate; i++) {
    data_in[i] = 0;
}
while (read_samples (infile, new_buffer, sfinfo.channels)==1)
{
    /* Process Samples */
    printf("Processing frame %d\n",nb_frames);

    /* hop size */
    fill_buffer(buffer, new_buffer);


    /* TODO */
    for (i=0; i<FRAME_SIZE; i++){
        // data_in[i] = buffer[i];
        // Fenetre de Hann :
        data_in[i] = buffer[i] * (0.5 - 0.5 * cos(2.0* M_PI * i/FRAME_SIZE));
    }
    fftw_execute (plan);

    // complex S[FRAME_SIZE];
    double amp[FRAME_SIZE];
    double phase[FRAME_SIZE];
    int ind_max = 0;
    int ind_max2 = 0;
    double max_val = amp[0];
    // dft(buffer, S);
    for (size_t i = 0; i < FRAME_SIZE; i++) {
        amp[i] = cabs(data_out[i]);
        phase[i] = carg(data_out[i]);
        if (i < FRAME_SIZE/2) {
            if (amp[ind_max] < amp[i]) {
                ind_max2 = ind_max;
                ind_max = i;
            } else if (amp[ind_max2] < amp[i]) {
                ind_max2 = i;
            }
        }
    }

    printf("frequence = %d / %d\n", ind_max, ind_max2);

    /* PLOT */
    gnuplot_resetplot(h);
    gnuplot_plot_x(h,amp,FRAME_SIZE,"temporal frame");
    // gnuplot_plot_x(h,phase,FRAME_SIZE,"phase");
    sleep(1);

    nb_frames++;
}
fftw_destroy_plan (plan);

sf_close (infile) ;

return 0 ;
} /* main */
