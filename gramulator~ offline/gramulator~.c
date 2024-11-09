#include "m_pd.h"
#include "math.h"
#define TWOPI 6.2831853072
#define BUFFER_EMPTY 0
#define BUFFER_FULL 1


static t_class *gramulator_class;

typedef struct _gramulator
{
  t_object obj;
  t_sample x_f;

  long delay_length;
  long delay_bytes; 
  long grain_length;
  long grain_bytes; 
  float maximum_delay_time;
  float *delay_line;
  float *grain_window;

  /* -- samphold -- */
  float lastin1;
  float lastout1;
  float lastin2;
  float lastout2;
  float sr;

  /* --tableread -- */
   int x_npoints;
   t_word *x_vec;
   t_symbol *x_arrayname;


} t_gramulator;



void grain_w(t_gramulator *x) { //read half of 2048 sample sine
  for(int a = 0; a < x->grain_length * 2; a++) {
  x->grain_window[a] = sin(1 * TWOPI * (float)a / (x->grain_length * 2));
  }    
} 


void *gramulator_new(t_symbol *s)
{

  t_gramulator *x = (t_gramulator *) pd_new(gramulator_class);
  
  inlet_new(&x->obj, &x->obj.ob_pd, gensym("signal"), gensym("signal")); //grain playhead
  outlet_new(&x->obj, gensym("signal")); //L
  outlet_new(&x->obj, gensym("signal")); //R
  outlet_new(&x->obj, gensym("signal")); //Sync

  //sr/
  x->x_arrayname = s;
  x->x_vec = 0;
  x->x_f = 0;

  x->sr = sys_getsr();
  x->grain_length = 4096;
  x->grain_bytes = (x->grain_length * 2) * sizeof(float) ; //for window; double size; read half
  x->grain_window = (float *)getbytes(x->grain_bytes);

  for(int i = 0; i < x->grain_length; i++){
    x->grain_window[i] = 0.0;
  }

  /* compute grain window */
  grain_w(x);

  x->lastin1 = 0.0;
  x->lastin2 = 0.0;
  x->lastout1 = 0.0;
  x->lastout2 = 0.0;

  return x;
}


// Precise method, which guarantees v = v1 when t = 1.
float lerp(float v0, float v1, float t) {
  return (1 - t) * v0 + t * v1;
}


t_int *gramulator_perform(t_int *w) {
  t_gramulator *x = (t_gramulator *) (w[1]);
  float *phasor1 = (float *)(w[2]);
  float *phasor2 = (float *)(w[3]);
  float *out1 = (float *)(w[4]);
  float *out2 = (float *)(w[5]);
  float *sync = (float *)(w[6]);
  int n = w[7];

  /* table reasd */  
  int maxindex;
  t_word *buf = x->x_vec;
  int i;
  float *grain_window = x->grain_window;
  float grain_length = x->grain_length;
  float playhead1, playhead2;
  float grain_a, grain_b;
  float ga_out, gb_out;
  float lastin1 = x->lastin1; 
  float lastout1 = x->lastout1;
  float lastin2 = x->lastin2; 
  float lastout2 = x->lastout2;
  float window1, window2;
  long read_1, read_2;
  float index;

 /* simple table read */ 
    maxindex = x->x_npoints - 1;
    if(maxindex<0) goto zero;
    if(!buf) goto zero;


for (i = 0; i < n; i++) {

index = *phasor1++;

/* keep within limits - maxindex can be used for linear read out table - otherwise excise code - unnecessary */
    
    if(index < 0)
    index = 0;
    else if (index > maxindex)
    index = maxindex;

/* ------------------ */  

playhead1 = index * x->x_npoints - 1; //main playhead; convert hz to samples (x * SR); SR / HZ = original playback speed
playhead2 = *phasor2++; //grain playhead

grain_a = playhead2 * grain_length; //grain 1
grain_b = fmod(playhead2 + 0.5, 1) * grain_length; //grain 2 - overlapping grain

if(grain_a < lastin1) lastout1 = playhead1; //output main playhead when grain read is done
lastin1 = grain_a; //grain

if(grain_b < lastin2) lastout2 = playhead1; //output main playhead when grain read is done
lastin2 = grain_b; //grain

window1 = grain_window[(long)grain_a];
window2 = grain_window[(long)grain_b];

ga_out = fmod(lastout1 + grain_a, x->x_npoints - 1); //index for samphold ph1 output plus grain1
gb_out = fmod(lastout2 + grain_b, x->x_npoints - 1);

      read_1 = (long)fmod(ga_out, x->x_npoints - 1);
      read_2 = (long)fmod(gb_out, x->x_npoints - 1);

          *out1++ = buf[read_1].w_float * window1;
          *out2++ = buf[read_2].w_float * window2;
          *sync++ = index;
    }

    x->lastin1 = lastin1;
    x->lastout1 = lastout1;
    x->lastin2 = lastin2;
    x->lastout2 = lastout2;

    return w + 8;
zero: //good way to force zero
    while (n--) *out1++ = 0; *out2++ = 0; *sync++ = 0;
    return w + 8;
}



static void gramulator_set(t_gramulator *x, t_symbol *s)
{
    t_garray *a;

    x->x_arrayname = s;
    if (!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
    {
        if (*s->s_name)
            pd_error(x, "tabread~: %s: no such array", x->x_arrayname->s_name);
        x->x_vec = 0;
    }
    else if (!garray_getfloatwords(a, &x->x_npoints, &x->x_vec))
    {
        pd_error(x, "%s: bad template for tabread~", x->x_arrayname->s_name);
        x->x_vec = 0;
    }
    else garray_usedindsp(a);
}


void gramulator_dsp(t_gramulator *x, t_signal **sp)
{
    gramulator_set(x, x->x_arrayname);
    x->delay_length = x->x_npoints;
    dsp_add(gramulator_perform, 7, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[0]->s_n);
}


void gramulator_free(t_gramulator *x)
{
  t_freebytes(x->grain_window, x->grain_length * 2);
}


void gramulator_tilde_setup(void)
{
  gramulator_class = class_new(gensym("gramulator~"), (t_newmethod)gramulator_new, (t_method)gramulator_free, sizeof(t_gramulator), CLASS_DEFAULT, A_DEFSYM, 0);
  class_addmethod(gramulator_class, (t_method)gramulator_dsp, gensym("dsp"),0);
  CLASS_MAINSIGNALIN(gramulator_class, t_gramulator, x_f);
  class_addmethod(gramulator_class, (t_method)gramulator_set, gensym("set"), A_SYMBOL, 0);
  post("gramulator~ v1: ricky@rickygraham.net");
}
