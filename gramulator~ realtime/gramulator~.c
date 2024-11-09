#include "m_pd.h"
#include "math.h"
#define TWOPI 6.2831853072

static t_class *gramulator_class;

typedef struct _gramulator
{
  t_object obj;
  t_sample x_f;

  int delay_length; 
  int grain_length;
  int grain_bytes; 
  int delay_bytes;
  float maximum_delay_time;
  float *delay_line;
  float *grain_window;
  int write_index;
  float sr;

  /* -- samphold -- */
  float lastin1;
  float lastout1;
  float lastin2;
  float lastout2;

  /* buffer record */
  short acquire_sample; 

} t_gramulator;



void grain_w(t_gramulator *x) { //read half of 4096 sample sine
  for(int a = 0; a < x->grain_length * 2; a++) {
  x->grain_window[a] = sin(1 * TWOPI * (float)a / (x->grain_length * 2));
  }    
} 


void *gramulator_new(void)
{

  t_gramulator *x = (t_gramulator *) pd_new(gramulator_class);
  
  inlet_new(&x->obj, &x->obj.ob_pd, gensym("signal"), gensym("signal")); //main playhead 
  inlet_new(&x->obj, &x->obj.ob_pd, gensym("signal"), gensym("signal")); //grain playhead
  outlet_new(&x->obj, gensym("signal")); //L
  outlet_new(&x->obj, gensym("signal")); //R
  outlet_new(&x->obj, gensym("signal")); //Sync

  x->sr = sys_getsr();
  x->grain_length = 4096;
  x->delay_length = 88200; //2 second buffer
  x->delay_bytes = x->delay_length * sizeof(float);
  x->grain_bytes = (x->grain_length * 2) * sizeof(float); //for window; double size; read half
  x->delay_line = (float *)getbytes(x->delay_bytes);
  x->grain_window = (float *)getbytes(x->grain_bytes);

  /* clear the delayline */
  for(int i = 0; i < x->delay_length; i++){
    x->delay_line[i] = 0.0;
  }

  for(int i = 0; i < x->grain_length * 2; i++){
    x->grain_window[i] = 0.0;
  }

  /* compute grain window */
  grain_w(x);

  x->acquire_sample = 0;
  x->write_index = 0;
  x->lastin1 = 0.0;
  x->lastin2 = 0.0;
  x->lastout1 = 0.0;
  x->lastout2 = 0.0;

  return x;
}



void sample(t_gramulator *x)
{
  /* Turn on the flag to acquire a sample */
  x->write_index = 0;
  x->acquire_sample = 1;
}


// Precise method, which guarantees v = v1 when t = 1.
float lerp(float v0, float v1, float t) {
  return (1 - t) * v0 + t * v1;
}


t_int *gramulator_perform(t_int *w) {
  t_gramulator *x = (t_gramulator *) (w[1]);
  t_sample *in = (t_sample *)(w[2]);
  t_sample *phasor1 = (t_sample *)(w[3]);
  t_sample *phasor2 = (t_sample *)(w[4]);
  float *out1 = (float *)(w[5]);
  float *out2 = (float *)(w[6]);
  float *sync = (float *)(w[7]);
  int n = w[8];

  float *delay_line = x->delay_line;
  float *grain_window = x->grain_window;
  int delay_length = x->delay_length;
  int grain_length = x->grain_length;
  int write_index = x->write_index;
  int grain_a, grain_b;
  float ga_out, gb_out;
  float lastin1 = x->lastin1; 
  float lastout1 = x->lastout1;
  float lastin2 = x->lastin2; 
  float lastout2 = x->lastout2;
  float out_l, out_r, window1, window2, sync_val;
  short acquire_sample = x->acquire_sample;
  int maxindex = 88200;
  int i;

/* Write */
if(acquire_sample) {

    for(i = 0; i < n; i++) {
      delay_line[write_index++] = in[i];
      sync_val = (float)write_index;
      //post("%f", delay_line[i]);
      
      if(write_index < 0) {
      write_index = 0;
      } else if (write_index >= maxindex) {
      acquire_sample = 0;
      }

    *out1++ = 0.0;
    *out2++ = 0.0;
    *sync++ = sync_val;

    }

} else if (acquire_sample == 0) {     


while(n--){
    
  int index = *phasor1++ * 88200;

  if (index < 0) {
      index = 0;
      }else if (index >= maxindex) {
      index = maxindex;
      
      } 

sync_val = (float)index;

float in2  = *phasor2++ ; //grain playhead - input 3
grain_a = in2 * grain_length;
grain_b = fmod(in2  + 0.5, 1) * grain_length; //grain 2 - overlapping grain

//post("%d", grain_a);
//post("%d", grain_b);

if(grain_a < lastin1) lastout1 = index; //samphold1
lastin1 = grain_a;

if(grain_b < lastin2) lastout2 = index; //samphold2
lastin2 = grain_b;

window1 = grain_window[(int)grain_a]; //grain window 1
window2 = grain_window[(int)grain_b]; //grain window 2

ga_out = fmod(lastout1 + grain_a, delay_length); //index for samphold1 output plus grain1
   
       /* limit boundaries */
      if(ga_out >= delay_length){
      ga_out -= delay_length; }

gb_out = fmod(lastout2 + grain_b, delay_length); //index for samphold2 output plus grain2

      if(gb_out >= delay_length){
      gb_out -= delay_length; }
    
      out_l = delay_line[(int)ga_out];
      out_r = delay_line[(int)gb_out];

    *out1++ = out_l * window1;
    *out2++ = out_r * window2;
    *sync++ = sync_val;
}

} 
    x->write_index = write_index; 
    x->acquire_sample = acquire_sample;
    x->lastin1 = lastin1;
    x->lastout1 = lastout1;
    x->lastin2 = lastin2;
    x->lastout2 = lastout2;
   
    return w + 9;
}



void clear(t_gramulator *x)
{
  for(int i = 0; i < x->delay_length; i++){
    x->delay_line[i] = 0.0;
  }
}



void gramulator_dsp(t_gramulator *x, t_signal **sp)
{
    dsp_add(gramulator_perform, 8, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, sp[0]->s_n);
}


void gramulator_free(t_gramulator *x)
{
  t_freebytes(x->delay_line, x->delay_bytes);
  t_freebytes(x->grain_window, x->grain_bytes);
}


void gramulator_tilde_setup(void)
{
  gramulator_class = class_new(gensym("gramulator~"), (t_newmethod)gramulator_new, (t_method)gramulator_free, sizeof(t_gramulator), CLASS_DEFAULT, A_GIMME, 0);
  class_addmethod(gramulator_class, (t_method)gramulator_dsp, gensym("dsp"),0);
  CLASS_MAINSIGNALIN(gramulator_class, t_gramulator, x_f);
  class_addmethod(gramulator_class, (t_method)clear, gensym("clear"),0);
  class_addmethod(gramulator_class, (t_method)sample, gensym("sample"), 0);
  post("gramulator~ v1: ricky@rickygraham.net");
}
