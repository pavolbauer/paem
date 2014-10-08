/* aemore.c - Core AEM solver. */

/* P. Bauer and S. Engblom 2012-05-10 */

// three lines to enable pinning of threads (together with _GNU_SOURCE)
#define _GNU_SOURCE
#include <sched.h>
#include <unistd.h>
#include <sys/syscall.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <pthread.h>
#include <inttypes.h>
#include <assert.h>
#include <glib.h>
#include <stdint.h>

#include "propensities.h"
#include "aem.h"
#include "binheap.h"
#include "stack.h"
#include "report.h"
#include "spsc_cb.h"
#include "rlcg.h"

//#define ROUNDOFF 1

//#define FUJIMOTO 1

/* The length of the prefix of the diffheap to scan */
#define DIFF_SCANLENGTH 30

/* update time fossil collection */
#define UPDATETIME 20000
/* CPU frequency */
#define CPUFREQ 2700000000UL
/* ?? */
#define PSTEP 27000UL
/* ?? */
#define PTOTALSTEPS 1000
/* ?? */
#define CL_SIZE 64
/* Rollback ON/OFF */
#define ROLLBACK
/* Fossil collection ON/OFF */
#define FOSSIL
/* Reversible RNG ON/OFF */
#define RLCG
/* Rollback test - time to rollback*/
#define RLB_DBG_T1 0.1
//#define RLB_DBG_T1 0.2
/* Rollback test - rollbacking time */
#define RLB_DBG_T2 3

/* for rollback testing */
void rlb_dbg_save();
int rlb_dbg_compare();
__thread int rlb_dbg_record=0;
__thread int rlb_dbg_fwd=0;
__thread double rlb_dbg_tt=0;
evnode_t *diffnodes_dbg, *reactnodes_dbg;
evtime_t *difftimes_dbg, *reacttimes_dbg;
__thread int *xx_dbg;


void peekdiffs();
void roundoff (double *val);
int updatediffs (int subvol, int spec, double tt, int multi, int rb);
int updatereacts (int subvol, int index, double tt, int rb);
void foreach_spec (int subvol, int reaction, double tt, double xxt, int modifier);

/* kind of accurate time measurement on late intel cpus */
static inline uint64_t __attribute__((always_inline))
read_tsc_p()
{
   uint64_t tsc;
   __asm__ __volatile__ ("rdtscp\n"
	 "shl $32, %%rdx\n"
	 "or %%rdx, %%rax"
	 : "=a"(tsc)
	 :
	 : "%rcx", "%rdx");
   return tsc;
}

#define now() read_tsc_p()

/* rollback list entry */
typedef struct rlb_s {
    double tt; // time of event
    double xxt;
    double toxxt;
    int type;  // type of eventf
    int id;    // event ID
    int dom;   // should be set to sender domain if type is msg
} rlb_t;

/* Global variables for debug statistics */
__thread unsigned long updates = 0;
__thread unsigned long wakeups = 0;
__thread unsigned long react2diff = 0;
__thread unsigned long react2react = 0;
__thread unsigned long diff2react = 0;
__thread unsigned long diff2diff = 0;

__thread int *to_nbs;
__thread int to_nbs_cnt = 0;
__thread int *from_nbs_q;



__thread int rollbacked_events = 0;
__thread int cancelled_msgs = 0;
__thread int rb_anti = 0;
__thread int rb_newanti = 0;
__thread int rb_ooo = 0;
__thread int rb_ooo_local = 0;

__thread long pdex = 0;
__thread long pdupdated = 0;
__thread long pdupdatedinf = 0;


/* pinning of threads */
/* straight core allocation on sandy (i.e., threads 0 - 15 on socket 1,
 *  threads 16-31 on socket 2, usw.)
 */
#define topo32(i) ( (i)/8 + 4*((i) % 8))
#define topo32ht(i) ( ((i)/2)/8 + (4*((i)/2) % 32) + ((i)%2)*32 )

#define PAD(_n) char __pad ## _n [CL_SIZE]
#define PAD_ALLOC(_s)                                           \
    ((void *)(((unsigned long)malloc((_s)+CL_SIZE*2) +  \
        CL_SIZE - 1) & ~(CL_SIZE-1)))

pid_t
gettid(void) 
{
    return (pid_t) syscall(SYS_gettid);
}

void
pin(pid_t t, int cpu)
{
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(cpu, &cpuset);
  sched_setaffinity(t, sizeof(cpu_set_t), &cpuset);
}

double
ts_diff_sec(struct timespec *start, struct timespec *end)
{
    double seconds = end->tv_sec - start->tv_sec;
    seconds += (end->tv_nsec - start->tv_nsec)/1000000000.0;
    return seconds;
}

/*
 *void calcTimes(double * restrict time, double * restrict infTime, const double tt, 
          const double old_rate, const double new_rate, uint64_t* seed) {
    double oldtime = time[0];
    if (isinf(oldtime)) {
        if (new_rate > 0.0) {    
            if (infTime[0] == 0) {
                time[0] = -log(1.0 - rng_next(seed)) / new_rate + tt;                            
            } else { 
                time[0] = tt + (infTime[0] / new_rate);
            }
        }
        wakeups++;
    } else {
        if (new_rate >= DBL_MIN) {
            if (time[0] == tt) {
                time[0] = -log(1.0 - rng_next(seed)) / new_rate + tt;               
            } else { 
                time[0] = ((old_rate / new_rate) * (time[0] - tt)) + tt; 
            }
            updates++;
        } else { 
            infTime[0] = (oldtime - tt) * old_rate;
            time[0] = INFINITY;            
        }
    }
}
 */

void calcTimes(double * restrict time, double * restrict infTime, const double tt, 
          const double old_rate, const double new_rate, uint64_t* seed) {
    double oldtime = time[0];
    if (isinf(oldtime)) {
       time[0] = -log(1.0 - rng_next(seed)) / new_rate + tt;                            
    } else {
        if (new_rate >= DBL_MIN) {
            if (time[0] == tt) {
                time[0] = -log(1.0 - rng_next(seed)) / new_rate + tt;               
            } else { 
                time[0] = ((old_rate / new_rate) * (time[0] - tt)) + tt; 
            }
            updates++;
        } else {
            time[0] = INFINITY;            
        }
    }
}


int test_heap_prty(evtime_t* data,int N)
{
  int i;
  for(i=0; i<(N-1)/2; i++)
    if(data[i].time>data[2*i+1].time || data[i].time>data[2*i+2].time)
      return 0; 
  return 1;
}

/* thread state */
typedef struct tst_s {
    int id;
    pthread_t t; // pthread specific (for join)
    PAD(0);
    long cnt; // count something
    PAD(1);
    double localtime;
    double wait;
    PAD(2);
    double time;
    PAD(3);
    cb_t **nbs; // neighbours?
    int nbscnt;
    int *nbsfrom; // from dom to queue
    PAD(4);
} tst_t;

static tst_t *tst;
PAD(0);
volatile int initblock = 0;
volatile int startblock = 0;
volatile int block = 0;

volatile int block1 = 0;


unsigned long glob_pdex = 0, glob_pdupdated = 0, glob_pdupdatedinf = 0;

unsigned long glob_rb_ooo = 0, glob_rb_anti = 0, glob_cancelled_msgs = 0, glob_rollbacked_events;


void *run(void *);

/* Global variables read by all threads */

const double *vol;
const int *sd;
const double *data;
const double *tspan;
static size_t tlen;
static int nthreads;
static int **domainLookup;
static int *g2lLookup;
static int *sd_length;


const size_t *irN, *jcN; const int *prN;
const size_t *irG, *jcG;
const size_t *irD, *jcD; const double *prD;
const int *u0;
int *globU;


static PropensityFun *rfun;
static size_t globNdofs;
static size_t globNcells;
static size_t Mspecies;
static size_t Mreactions;
static size_t dsize;
ReportFun report;     
static int report_level;
static unsigned long long starttime;

PAD(1);

double minitime;
int settime = 0;
PAD(2);

#ifdef FUJIMOTO
  /* Flag to keep count of computing threads */
	int GVTFlag;
  /* Generation of GVT */
  int GVTGen;
  /* Minimum for every PE */
	double* PEMin;

  /* Minimum sent message timestamp */
	__thread double SendMin=INFINITY;
  /* Copy of GVTFlag */
	__thread int LocalGVTFlag;
  /* Make sure GVT compues only once */
  int GVTCnt=0;
  __thread int GVTUpdateCnt=0;
#endif

#ifdef RLB_DBG_T1
  rlb_cnt=0;
#endif

/* Thread specific global variables */
__thread double *xxt;
__thread int *xx;

__thread int diff_updated = 0;

__thread int *U, *irE, *jcE, *vxs, localdom;//, dom;

__thread evnode_t *diffnodes, *reactnodes;
__thread evtime_t *difftimes, *reacttimes;

__thread size_t Ncells, Ndofs, diffHeapSize, reactHeapSize;

__thread unsigned short rng[3];

/* the pseudo-diagonal of D is used very frequently */
__thread double *diag;

__thread tst_t *st;


/* Rollback list */
__thread GList *rblist;
__thread GList *rblistend;

__thread int antidom;

//__thread double wavefront[PTOTALSTEPS];

//double globwavefront[32][PTOTALSTEPS];



void aem_core(const size_t *_irD,const size_t *_jcD,const double *_prD,
              const int *_u0,
              const size_t *_irN,const size_t *_jcN,const int *_prN,
              const size_t *_irG,const size_t *_jcG,
              const double *_tspan,const size_t _tlen,
              int *_U,
              const double *_vol,const double *_data,const int *_sd,
              const size_t _Ncells,
              const size_t _Mspecies,const size_t _Mreactions,
              const size_t _dsize,int _report_level, int _nthreads,
		long int *timeHorizon)




/* Specification of the inputs:

 Ncells
 Number of subvolumes.

 Mspecies
 Number of species.

 Hence Ndofs = Ncells*Mspecies.

 Mreactions
 Total number of reactions.

 dsize
 Size of data vector sent to propensities.

 tlen
 Number of sampling points in time.

 report_level
 The desired degree of feedback during simulations. 0, 1, and 2 are
 currently supported options.

 Diffusion matrix D. Double sparse (Ndofs X Ndofs).
 Macroscopic diffusion matrix. D(i,j) is the diffusion rate from dof #j to
 dof #i. This matrix uses the CSR-format and not CSC because fast access to
 rows is needed.

 Initial state vector u0. Integer (Mspecies X Ncells).
 Gives the initial copy number of the species in each subvolume.

 Stochiometric matrix N. Integer sparse (Mspecies X Nreactions).
 N(:,j) describes how reaction j changes the number of species.

 Dependency graph G. Integer sparse (Mreactions X Mspecies+Mreactions).
 G(i,Mspecies+j) is non-zero if executing reaction j means that reaction i
 needs to be re-evaluated. The first Mspecies columns of G similarily cover
 diffusion events.

 tspan. Double vector.
 Output times. tspan[0] is the start time and tspan[length(tspan)-1] is the
 stop time.

 vol. Double vector (length Ncells).
 vol[i] gives the volume of cell #i.

 data. Double matrix (dsize X Ncells).
 Generalized data matrix, data(:,j) gives a data vector for cell #j.

 sd. Integer vector (length Ncells).
 Subdomain number. sd[i] is the subdomain of cell #i. The vector sd can also
 be used to separate boundaries, line segments and points.

 Format of sparse matrices:
 G, N and S are sparse matrices in compressed column format (CCS). D is sparse
 but in compressed row format (CRS), or equivalently, a transposed matrix in
 CCS format.
 jcD, irD, prD (double *)
 jcN, irN, prN (int *)
 jcG, irG (int *)

 Propensities:
 a vector of function pointers (length Mreactions) is input by
 linking with the prototypes in propensities.h and function
 definitions in a user-specified .c-file. The type of this vector is
 PropensityFun which defines the input to a property function. See
 propensities.h for more details.

 Ordering of the dofs:
 Dof #i is located in cell #(i/Mspecies), and the dofs located in
 cell #j is u0(:,j). Thus, u0 is understood as a matrix of size
 Mspecies X Ncells.

 The output is a matrix U (Ndofs X length(tspan)).
 U(:,j) contains the state of the system at tspan(j).
 */
{
  size_t i,j = 0;


  globNcells   = _Ncells;
  Mspecies     = _Mspecies;
  Mreactions   = _Mreactions;
  dsize        = _dsize;
  report_level = _report_level;
  globNdofs    = _Ncells*Mspecies;
  vol          = _vol;
  sd           = _sd;
  data         = _data;
  tspan        = _tspan;
  tlen         = _tlen;

  jcN = _jcN;
  irN = _irN;
  prN = _prN;
  irD = _irD;
  jcD = _jcD;
  prD = _prD;
  irG = _irG;
  jcG = _jcG;
  u0  = _u0;
  globU = _U;

  minitime = 0.0;
  
  /* build up seed reaction+diffusion matrix here */

  struct timespec startsim,endsim;
  clock_gettime(CLOCK_MONOTONIC, &startsim);

#ifdef FUJIMOTO
  PEMin = calloc(nthreads, sizeof (double *));
#endif

  rfun = ALLOC_propensities();

  nthreads = _nthreads;

  ReportFun report = NULL;
  if (report_level) report = &reportFun1;

  /* create subdomain lookup */
  sd_length=calloc(nthreads,sizeof(int));
  for(i = 0; i < globNcells; i++)
      sd_length[sd[i]-1]++;

  domainLookup = malloc (nthreads*sizeof(int*));

  int cnt = 0;  
  g2lLookup = malloc (globNcells*sizeof(int));
  for (i = 0; i < nthreads; i++)
  {
      domainLookup[i] =  malloc(sd_length[i]*sizeof(int));

      for(cnt = 0, j = 0; j < globNcells; j++) {
          if( (sd[j]-1) == i) {
              domainLookup[i][cnt] = j;
              g2lLookup[j] = cnt;
              cnt++;
          }
      }
  }  

  
  tst = PAD_ALLOC(nthreads * sizeof(tst_t));
  

  // first init structs
    for (int i = 0; i < nthreads; i++)
    {
        tst[i].id = i;
        tst[i].cnt = 0;
        tst[i].time = 0.0;
        tst[i].nbscnt = 0;
        tst[i].nbsfrom = PAD_ALLOC(sizeof(int) * nthreads);
        for (j = 0; j < nthreads; j++)
            tst[i].nbsfrom[j] = -1;
    }
    // then start threads
    for (int i = 0; i < nthreads; i++) {
        pthread_create (&tst[i].t, NULL, run, (void *)&tst[i]);
    }

    while (startblock < nthreads)
        ;

    starttime = now();
    __sync_synchronize();
    
    // tso, no barrier needed.
    startblock = 0;

    for (int i = 0; i < nthreads; i++)
    {
        pthread_join (tst[i].t, NULL);
    }

/*
    FILE *wffile;
    wffile = fopen("wavefront.dat","w");
    for (i = 0; i < PTOTALSTEPS; i++) {
        for (j = 0; j < nthreads; j++) {
            fprintf(wffile, "%.10f ", globwavefront[j][i]);
        }
        fprintf(wffile, "\n");
    }
    fclose(wffile);
*/
    
#ifdef TEST
    printf("%f ", tspan[tlen - 1]);
    for (j = 0; j < globNdofs; j++) {
        printf("%d ", _U[(tlen - 1)*globNdofs + j]);
    }
    printf("\n");
#endif
    
    int ooo = 0;
    int msgcnt = 0;
    unsigned long long totalwait = 0;
    
    for (int i = 0; i < nthreads; i++) {
        ooo += tst[i].cnt;
        totalwait += tst[i].wait;
        for (int d = 0; d < tst[i].nbscnt; d++) {
            msgcnt += tst[i].nbs[d]->cnt;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &endsim);


//    printf("%lu %lu %lu\n", glob_pdex, glob_pdupdated, glob_pdupdatedinf);
    printf("%lu %lu %lu %lu\n", glob_cancelled_msgs, glob_rollbacked_events, glob_rb_anti, glob_rb_ooo);
    

    printf("%d %.1f %d %d %d %.4f %.1f\n", nthreads, 
           ts_diff_sec(&startsim, &endsim), 
           glob_cancelled_msgs+glob_rollbacked_events, ooo, msgcnt, 
           (double)ooo/(double)msgcnt, 
           (double)totalwait/(double)CPUFREQ/nthreads);

#ifdef RLB_DBG_T1
    printf("Number of rollbacks: %d \n",rlb_cnt);
#endif
//}
             
//    free (tst);
//    FREE_propensities(rfun);
}



inline void
roundoff (double *val) {
#ifdef ROUNDOFF
    if (fabs(*val) < 1.0E-12) {
        *val = 0.0;
    }
#endif
}

/* returns a mapping from local voxel id to global voxel id for
 * subdomain subdom */
int *
lookupSubDomain(int subdom)
{
    return domainLookup[subdom];
}


/* input : global voxel id
 * output: subdomain id of voxel 
 */
int 
lookupVoxel (int voxel)
{
    return sd[voxel] - 1;
}

/* input : global voxel id
 * output: local voxel id
 * NOTE: the return value is only valid when 
 * the calling thread is the owner of the 
 * voxel. */
int
lookupLocalVoxel (int global_voxel) {
    return g2lLookup[global_voxel];
}


double 
readNeighbourDiffTimesMD (int id, int *dir, int *subvol) 
{
    int mindir = -1;
    int minvol = -1;
    double min = INFINITY;
    for (int i = 0; i < st->nbscnt; i++) {
        if (st->nbs[i]->nextdiff < min) {
            min = st->nbs[i]->nextdiff;
            mindir = i;
            minvol = st->nbs[i]->nextsubvol;
        }
    }
    *dir = mindir;
    *subvol = minvol;
    return min;
}


static void 
idle_wait(unsigned int ticks) 
{
    volatile int x;
    for (x = 0; x < ticks; x++);
}


#ifdef ROLLBACK

static gint
compare_rlb (gconstpointer a, gconstpointer b) 
{
    const rlb_t *rlb0 = a, *rlb1 = b;
    
    return ((rlb0->tt < rlb1->tt) ? 1 : (rlb0->tt == rlb1->tt) ? 0 : -1);
}


static void
rb_push_event (int type, double time, double xxt, int id, int dom, double toxxt)
{
   rlb_t *rb;
   rb = g_slice_new (rlb_t);
   rb->type = type;
   rb->id = id;
   rb->tt = time;
   rb->xxt = xxt;
   rb->dom = dom;
   rb->toxxt = toxxt;
   rblist = g_list_insert_sorted(rblist, (gpointer)rb,
                                     (GCompareFunc)compare_rlb);
}

static double
update_gvt(double calc_min)
{
    double new_local;
    if (calc_min == minitime && settime == 0 && !__sync_fetch_and_add(&settime, 1)) {
        // i.e., has not changed, update gvt
        new_local = DBL_MAX;
        for (int i = 0; i < nthreads; i++) {
            new_local = (tst[i].time < new_local? tst[i].time : new_local);
        }
        assert(new_local >= minitime);
        minitime = new_local;
        settime = 0; // signal updating is done
    } else {
        new_local = minitime;
    }
    return new_local;
}

void
fossil_collect()
{
    int cnt = 0;
    rlb_t *event;
    GList *prev, *node = g_list_previous(rblistend);
    if (!node) return;
    event = node->data;
    while (event && event->tt < minitime && cnt < 100) {
        cnt++;
        g_slice_free(rlb_t, node->data);
        prev = g_list_previous(node);
        rblist = g_list_delete_link(rblist, node); 
        node = prev;
        if (!node) break;
        event = node->data;
    }
}

/*-----------------------------------------------------------------------*/
/* Rollback one event */

void
rollback_event (const rlb_t * restrict event, double * restrict nbs)
{  
    int i;
    msg_t msg, msgtmp;
    
    int subvol,re,spec,to_vol,to_spec;
    int localpos, dom;
    
  if (event->type == 0) {
    subvol = event->id / Mreactions;
    re = event->id % Mreactions;
    
    /* includes updatediffs */
    foreach_spec (subvol, re, event->tt, event->xxt, -1);
    
    /* c) Update dependent reaction events (including the
     * actual one) */

    updatereacts (subvol, Mspecies + re, event->tt, 1);

  } else if (event->type == 1) {
      for (i = 0; event->id >= jcE[i]; i++)
          ;
      subvol = (i - 1) / Mspecies;
      spec = (i - 1) % Mspecies;
      to_vol = irE[event->id] / Mspecies;
      to_spec = irE[event->id] % Mspecies;
      
      /* revert local diffusion */
      xx[subvol * Mspecies + spec]++;
      assert(xx[subvol * Mspecies + spec]>=0);
      xxt[subvol] = event->xxt;
      
      localpos = lookupLocalVoxel(to_vol);
      dom = lookupVoxel(to_vol);
      assert(event->tt > 0.0);
      
      /* message will serve to determine rollback cascades */
      if (localdom != dom) {
          if (nbs[st->nbsfrom[dom]] > event->tt)
              nbs[st->nbsfrom[dom]] = event->tt;
          
          /* empty saved msgs from domain to which we send 
           * an anti message */
          while (st->nbsfrom[dom] >= 0 && 
                 bs_peek(st->nbs[st->nbsfrom[dom]]->stack, &msg) && 
                 msg.time >= nbs[st->nbsfrom[dom]]) {
              bs_pop(st->nbs[st->nbsfrom[dom]]->stack, &msg);
          }
      }
      else if (localdom == dom) {
          /* only update the incoming subvolume if it is internal to the
           * outgoing subvolume's subdomain */
        
          /* revert incoming subvolume if local */
          xx[localpos * Mspecies + to_spec]--;

          xxt[localpos] = event->toxxt;
          assert(xx[localpos * Mspecies + to_spec]>=0);
          
          /* update dependent reaction events of incoming subvolume */
          updatereacts(localpos, to_spec, event->tt, 0);
          /* Update all diffusion events of affected species in
           * the incoming subvolume */
          updatediffs(localpos, to_spec, event->tt, -1, 0);
      }
      
      /* update reaction events of outgoing subvolume */
      updatereacts(subvol, spec, event->tt, 0);
      
      /* Update all diffusion events of affected species in
       * the outgoing subvolume */
      /* In the AEM context decrease the intensity of all
       * events related to non-diagonal row entries of D */
      
      
      updatediffs(subvol, spec, event->tt, 1, 1);

  } else { // type 2
      // recvd diffusion msg from another node
      
      xx[event->id]--;
      
      assert(xx[event->id]>=0);
      subvol = event->id / Mspecies;
      spec   = event->id % Mspecies;
      xxt[subvol] = event->xxt;
      /* external msgs should be saved in order to be re-received */
      /* if they do not come from the node causing the rollback */
      /* if ooo and not anti, ooo node has not sent later diffusion messages
       * than the one that is out of order */
      if (antidom != event->dom) {
          assert(event->tt > 0.0);
          msgtmp = (msg_t) {
              .time = event->tt,
              .type = 0,
              .to_vol = vxs[subvol],
              .to_spec = spec,
          };
          assert(cb_pushback(st->nbs[st->nbsfrom[event->dom]], msgtmp));
      }
      
      updatereacts(subvol, spec, event->tt, 0);
      /* Update all diffusion events of affected species in
       * the incoming subvolume */
      updatediffs(subvol, spec, event->tt, -1, 0);
  }
}




void
rollbackToTime (const double rolltt) {

    int i, dom;
    msg_t msg;
    GList *l = g_list_first(rblist);
    rlb_t* event = l->data;
    double nbs[nthreads];

    assert(rolltt != 0.0);

    for (i = 0; i < nthreads; i++) 
        nbs[i] = DBL_MAX;

    while (event->tt >= rolltt)
    {
        GList *next = g_list_next(l);
        rollback_event(event, nbs);
        rollbacked_events++;
#ifdef RLB_DBG_T1
        printf("rollbacking to %f (type %d).\n",event->tt,event->type);
        rlb_cnt++;
        printf("reactnode 1 rate: %f\n",reactnodes[1].rate);
#endif
        g_slice_free (rlb_t, l->data);
        rblist = g_list_delete_link (rblist, l);
        l = next;
        event = l->data;
    }
    
//    printf("%d: rb to %.6f\n", st->id, event->tt);
    
    // check if necessary to send anti message
    for (i = 0; i < st->nbscnt; i++) {
        if (nbs[i] < DBL_MAX) {
            msg = (msg_t) {
                .type = 1,
                .time = nbs[i],
                .to_vol = -1,
                .to_spec = -1,
            };
            dom = from_nbs_q[i];
            while (!cb_push(tst[dom].nbs[to_nbs[dom]], msg));
            __sync_fetch_and_add(&tst[dom].nbs[to_nbs[dom]]->anti, 1);
        }
    }
    antidom = -1;
}

#endif

void *
run (void *_args)
{
    int pcnt = 0;
    double tt = tspan[0];
    int i, j, k, re, cnt;
    int event;
    int dom;
    int subvol, spec, to_vol, to_spec;
    
    short errcode = 0;
    long int total_diffusion = 0, total_reactions = 0;

    size_t it = 0;

    int local_i = -1, global_i;
    int myNnz = 0;
    int updatetimecnt = 0;
    double tn = DBL_MAX;
    msg_t msgp;
    msg_t msg;

    double tnm;
    int localpos;
    int next = -1;
    int nextdiffdir = -1;
    int rcvd = 0;
    double tempmin;
    
    unsigned long long localwait = 0;
    unsigned long long waitstart;

    double min_nbs,min_cb;

// for FOSSIL
    double calc_min = 0.0;

#ifdef ROLLBACK
    // Insert dummy element so that we have a pointer to 
    // the last element.
    rb_push_event (-1, 0.0, 0.0, 0, -1, 0.0);
    rblistend = rblist;
#endif
    
    #ifdef RLCG
    rng_init();
    #endif

    // thread state
    st = (tst_t *)_args;

    pin(gettid(), topo32ht(st->id));

    to_nbs = PAD_ALLOC(sizeof(int) * nthreads);
    memset(to_nbs, -1, sizeof(int) * nthreads);

    //initialize thread local rng
    for (i = 0; i < 3; i++)
        rng[i] = rand() % 65535;

    localdom = st->id;
    vxs = lookupSubDomain(localdom);
    
    Ncells = sd_length[st->id];
    Ndofs = Ncells*Mspecies;

    U = PAD_ALLOC(tlen * Ndofs * sizeof(int));

    /* Set xx to the initial state. */
    xx = PAD_ALLOC (Ndofs * sizeof(int));

    /* keep track of update time for each voxel */
    xxt = PAD_ALLOC (Ncells * sizeof(double));
    for (i = 0; i < Ncells; i++) 
        xxt[i] = 0.0;

    /* read in initial state */
    for (i = 0; i < Ncells; i++) 
        for (j = 0; j < Mspecies; j++) 
        {
            xx[(i*Mspecies)+j] = u0[vxs[i]*Mspecies+j]; 
        }

    /* Create Reaction-heap. */
    /* Binary heap storing all reaction events */
    reactHeapSize = Ncells * Mreactions;

    reactnodes = PAD_ALLOC(reactHeapSize * sizeof(evnode_t));
    reacttimes = PAD_ALLOC(reactHeapSize * sizeof(evtime_t));

    /* Create reaction rate matrix (Mreactions X Ncells) and total rate
       vector. In rrate we store all propensities for chemical rections,
        */

    /* Calculate the propensity for every reaction and every
       subvolume. */
    for (i = 0; i < Ncells; i++) {
        for (j = 0; j < Mreactions; j++) {
            reactnodes[i * Mreactions + j].rate = 
                (*rfun[j])(&xx[i * Mspecies], tt, vol[vxs[i]],
                           &data[vxs[i] * dsize], sd[vxs[i]]);
        }
    }

    /* make sure rand is seeded differently on every core */
    srand((long int)st->id);	
    
    /* Initialize reaction seeds */
    for (i = 0; i < reactHeapSize; i++) {
  #ifdef RLCG      
      reactnodes[i].seed = rand(); 
  #else
      for (j = 0; j < 3; j++)
        reactnodes[i].seed[j]=rand() % 65535;
  #endif
    }

    /* Calculate times to next reaction event */
    for (i = 0; i < reactHeapSize; i++) {
   #ifdef RLCG
        reacttimes[i].time = -log(1.0 - rng_next(&reactnodes[i].seed)) / 
                             reactnodes[i].rate + tspan[0];
    #else
        reacttimes[i].time = -log(1.0 - erand48(reactnodes[i].seed)) / 
                             reactnodes[i].rate + tspan[0];
    #endif
        if (reacttimes[i].time <= 0.0) 
            reacttimes[i].time = INFINITY;

        reacttimes[i].node = reactnodes[i].heapidx = i;
        reactnodes[i].inftime = 0.0;
        //printf("thr<ead %d. reacttimes[%d]: %f \n",st->id,i,reacttimes[i].time);
    }

    /* Initialize reaction heap */
    initialize_heap(reacttimes, reactnodes, reactHeapSize);
        
    
    for (i = 0; i < Ncells; i++) {
        for (k = 0; k<Mspecies; k++) {
            /* ndof_i -> mapped into global Ndofs adress */
            global_i=(vxs[i]*Mspecies)+k;
            for (j = jcD[global_i]; j < jcD[global_i+ 1]; j++) {
                if (irD[j] != global_i)
                    myNnz++;
            }
        }
    }
    
    /* diffHeapSize for unidirectional diffusion problems -> subtract
     * non-zero'd Nrows only */
    diffHeapSize = myNnz;
    
    assert(myNnz >= 0);

    difftimes = PAD_ALLOC(diffHeapSize * sizeof(evtime_t));
    diffnodes = PAD_ALLOC(diffHeapSize * sizeof(evnode_t)); 
    
    /* Initialize diffusion seeds */
    for (i = 0; i < diffHeapSize; i++) {
  #ifdef RLCG      
      diffnodes[i].seed = rand(); 
  #else
      for (j = 0; j < 3; j++)
        diffnodes[i].seed[j]=rand() % 65535;
  #endif
    }    
    
    /* Creating new sparse matrix for diffusion events */
    jcE = PAD_ALLOC ((Ndofs + 1) * sizeof(int));
    irE = PAD_ALLOC (diffHeapSize * sizeof(int));

    diag = PAD_ALLOC (diffHeapSize * sizeof(double));

    /* Queue all non-diagonal entries of D as diffusion events */
    for (cnt = 0, i = 0; i < Ncells; i++) {
        for(k=0; k<Mspecies; k++) {
            local_i=(i*Mspecies)+k;
            global_i=(vxs[i]*Mspecies)+k;
            jcE[local_i] = cnt;
            for (j = jcD[global_i]; j < jcD[global_i + 1]; j++)
                if (irD[j] != global_i) { // if not diagonal

                    // set up comm.
                    dom = lookupVoxel(irD[j]/Mspecies);
                    
                    if (dom != st->id && to_nbs[dom] < 0) {
                        to_nbs_cnt ++;
                        to_nbs[dom] = __sync_fetch_and_add(&tst[dom].nbscnt, 1);
                        // node with queues should know which domain each queue is
                        tst[dom].nbsfrom[st->id] = to_nbs[dom];
                    }
                    diffnodes[cnt].internode = dom != st->id;
                    
                    diffnodes[cnt].heapidx = difftimes[cnt].node = cnt;
                    diffnodes[cnt].rate = xx[i*Mspecies+k] * prD[j];
#ifdef RLCG                    
                    difftimes[cnt].time = -log(1.0 - rng_next(&diffnodes[cnt].seed)) / 
                                          diffnodes[cnt].rate + tspan[0];
#else
                    difftimes[cnt].time = -log(1.0 - erand48(diffnodes[cnt].seed)) / 
                                          diffnodes[cnt].rate + tspan[0];                    
#endif
                    if (difftimes[cnt].time <= 0.0) 
                        difftimes[cnt].time = INFINITY;
                    irE[cnt] = irD[j];
//                    difftimes[cnt].irE = irD[j];
                    diag[cnt] = prD[j];
                    diffnodes[cnt].inftime = 0.0;
                    cnt++;
                }
        }
    }
    jcE[local_i+1] = cnt;

    /* Initialize diffusion heap */
    initialize_heap(difftimes, diffnodes, diffHeapSize);
    
    /* Fix to run reaction-only tests */
    if (diffHeapSize == 0) {
        difftimes = malloc(sizeof(evtime_t));
        difftimes[0].time = INFINITY;
	}

    /* Fix to run diffusion-only tests */
    if (reactHeapSize == 0) {
        reacttimes = malloc(sizeof(evtime_t));
        reacttimes[0].time = INFINITY;
    }
    
    if (report_level == 2) {
        printf("React-Heap Size: %d\n", (int) reactHeapSize);
        printf("Diff-Heap Size: %d\n", (int) diffHeapSize);
    }

    /* barrier, guarantee that all threads are initialized */
    __sync_fetch_and_add(&initblock, 1);
    while (initblock < nthreads);

    // set up nbs channels

    // set up msgs channels to neighbours
    st->nbs = PAD_ALLOC(sizeof (cb_t *) * st->nbscnt);
    for (i = 0; i < st->nbscnt; i++) {
        st->nbs[i] = cb_init(2048);
    }

    from_nbs_q = PAD_ALLOC(sizeof (int) * st->nbscnt);

    for (i = 0; i < nthreads; i++){
        if (st->nbsfrom[i] >= 0) {
            from_nbs_q[st->nbsfrom[i]] = i;
        }
    }

    /* barrier, guarantee that all threads are initialized */
    __sync_fetch_and_add(&startblock, 1);
    while (startblock);

    peekdiffs();
    
    msg_t dummy;
    
    /* Main loop. */
    for (;;) {
      

    latereceive:
        /* Calculate tt with diff-time and react-time */
        tt = MIN(reacttimes[0].time,difftimes[0].time);

/*        
        if(now() > pcnt*PSTEP + starttime) {
            if (pcnt == PTOTALSTEPS) {
                for (i = 0; i < pcnt; i++) {
                    globwavefront[st->id][i] = wavefront[i];
                }
                cb_reset(st->nbs[0]);
                cb_reset(st->nbs[1]);
                st->nextdiff[0] = INFINITY;
                st->nextdiff[1] = INFINITY;

                return NULL;
            }
            wavefront[pcnt++] = tt;

        }
*/      


#ifdef RLB_DBG_T1
        if(tt>RLB_DBG_T1 && rlb_dbg_record && tt<=RLB_DBG_T2) {
            rlb_dbg_fwd++;
            printf("%d iteration: %f.\n",rlb_dbg_fwd,tt);
        }
        /* save state */
        if(tt>RLB_DBG_T1 && !rlb_dbg_record) {
            rlb_dbg_record=1;
            rlb_dbg_save();
            rlb_dbg_tt=tt;
            rlb_dbg_fwd=0;
            /*
            printf("reactnode 1 rate: %f\n",reactnodes[1].rate);
            printf("saving at %f.\n",tt);
            
                        printf("Did %d forward iterations and %d rollbacks.\n",rlb_dbg_fwd,rlb_cnt);
                        for(int c=0; c<reactHeapSize; c++) {
              printf("reacttimes; %d: %f, ",c,reacttimes[c].time);
            }
            printf("\n");
            assert(test_heap_prty(reacttimes,reactHeapSize));
            
            for(int c=0; c<diffHeapSize; c++) {
              printf("difftimes; %d: %f, ",c,difftimes[c].time);
            }
            printf("\n");
            assert(test_heap_prty(difftimes,diffHeapSize));
            */
        }
        
 
        if(tt>RLB_DBG_T1 && rlb_dbg_record) {
            if(reacttimes[0].time<difftimes[0].time)
              printf("reaction event at : %f.\n",tt);
            else
              printf("diffusion event at: %f.\n",tt);
            printf("reactnode 1 rate: %f\n",reactnodes[1].rate);
        }
        
        
        /* do rollback */
        if(rlb_dbg_fwd==RLB_DBG_T2) {
            //printf("rollbacking to %f at tt %f\n",rlb_dbg_tt,tt);
            rollbackToTime(rlb_dbg_tt-DBL_EPSILON);
          
            printf("Did %d forward iterations and %d rollbacks.\n",rlb_dbg_fwd,rlb_cnt);
            
            /*
            for(int c=0; c<reactHeapSize; c++) {
              printf("reacttimes; %d: %f, ",c,reacttimes[c].time);
            }
            printf("\n");
            assert(test_heap_prty(reacttimes,reactHeapSize));
            
            for(int c=0; c<diffHeapSize; c++) {
              printf("difftimes; %d: %f, ",c,difftimes[c].time);
            }
            printf("\n");
            assert(test_heap_prty(difftimes,diffHeapSize));
            */
            
            if(rlb_dbg_compare())
              printf("Rollback test passed.\n");
            else
              printf("Rollback test failed.\n");
            assert(0);
        }
#endif

#ifdef FUJIMOTO
        LocalGVTFlag=GVTFlag;
#endif
        diff_updated = 0;
         /* process all msgs occurring before next local event (tt) */
        while (1) {
            
            /* process earliest msg first */
            rcvd = 0;
            tnm = DBL_MAX;
            for (int dir = 0; dir < st->nbscnt; dir++) {
                if (cb_peek (st->nbs[dir], &msgp)
                    && msgp.time < tnm) {
                    next = dir;
                    tnm = msgp.time;
                    rcvd = 1;
                } 
            }
            // next is now the direction which has the message with 
            // the earliest timestamp

#ifdef ROLLBACK
            // if an anti-message will occur later in the message queue,
            // pop all messages prior to the anti-message
            if (rcvd && st->nbs[next]->anti) {
                if (!st->nbs[next]->nextanti > 0.0) {
                    st->nbs[next]->nextanti = cb_antipeek(st->nbs[next]);
                }
                if (cb_peek(st->nbs[next],&msgp) && msgp.time >= st->nbs[next]->nextanti) {
                    while (cb_pop(st->nbs[next],&msgp) && msgp.type != 1) 
                    {
                        cancelled_msgs++;
                    }
                    antidom = from_nbs_q[next];
                    rollbackToTime(msgp.time);
                    st->cnt++;
                    rb_anti++;
                    __sync_fetch_and_add(&st->nbs[next]->anti, -1);
                    st->nbs[next]->nextanti = 0.0;
                    tt = MIN(reacttimes[0].time,difftimes[0].time);
                    continue;
                }
             }
#endif

            if (rcvd && tnm <= tt) {
                cb_pop (st->nbs[next], &msgp);
            } else {
                if (rcvd && tnm > tt)
                tnm = DBL_MAX;
            }

            if (!rcvd || tt < tnm) break;
            // message is rcvd.
            to_vol = msgp.to_vol;
            subvol = lookupLocalVoxel(to_vol);
            spec = msgp.to_spec;

            assert(msgp.time != 0.0);
            
            // is msg older than last update of xx for subvol?
            // or is it an antimessage
            if (msgp.type == 1 || msgp.time < xxt[subvol]) {
                // count out of order msgs / anti
                st->cnt++;
#ifdef ROLLBACK
                if (msgp.type == 1) {
                    antidom = from_nbs_q[next];
                    rb_anti++;
                } else {
                    assert(msgp.time > minitime);
                    int outgoingisempty = 1;
                    for (i = 0; i < nthreads; i++){
                        if (to_nbs[i] >= 0) {
                            outgoingisempty &= tst[i].nbs[to_nbs[i]]->buf[tst[i].nbs[to_nbs[i]]->r].active == 0;
                            if (!outgoingisempty) break;
                        }
                        
                    }
                    if (outgoingisempty) st->time = msgp.time; //update local time value
                    antidom = -1;
                    rb_ooo++;
                }
                rollbackToTime(msgp.time);
                if (msgp.type == 1) {
                    __sync_fetch_and_add(&st->nbs[next]->anti, -1);
                    tt = MIN(reacttimes[0].time,difftimes[0].time);
                    continue;
                }
#endif
                tt = MIN(reacttimes[0].time,difftimes[0].time);
            }

#ifdef ROLLBACK
            rb_push_event (2, msgp.time, xxt[subvol], subvol * Mspecies + spec, from_nbs_q[next], 0.0);
#endif            
            xx[subvol * Mspecies + spec]++;
            xxt[subvol] = msgp.time;
                    
            //Recalculate the reaction rates using dependency graph G.
            updatereacts(subvol, spec, msgp.time, 0);
            
            // dependent diffusion events
            updatediffs(subvol, spec, msgp.time, 1, 0);
            
            tt = MIN(reacttimes[0].time,difftimes[0].time);
            
        }
        if (diff_updated)
            peekdiffs();

        diff_updated = 0;

        int nextsubvol;
        /* time of next incoming diffusion event (in the future) */
        tn = readNeighbourDiffTimesMD(st->id, &nextdiffdir, &nextsubvol);

        waitstart = now();




        /* wait for neighbour - fossil collect while waiting */
        while(tt > tn && !st->nbs[nextdiffdir]->anti && tn != 0.0 && !isinf(tt) && tt < tspan[tlen-1]) {

#ifdef FOSSIL
            updatetimecnt++;
            if (updatetimecnt > UPDATETIME) {
                updatetimecnt = 0;
#ifdef FUJIMOTO
                /* Initiate GVT computation */
                if(__sync_bool_compare_and_swap(&GVTFlag, 0, nthreads)) {
                  GVTCnt++;
                  printf("Thread %d initiating GVT computation at tt=%f.\n",st->id,tt);
                }
                
#else
                calc_min = update_gvt(calc_min);
#endif
            }
            if (!(updatetimecnt % 10))
                fossil_collect();
#endif

            for (int dir = 0; dir < st->nbscnt; dir++) {
                if (cb_peek (st->nbs[dir],  &dummy)) {
                    localwait += now() - waitstart;
                    goto latereceive;
                }
            }
            tn = readNeighbourDiffTimesMD(st->id, &nextdiffdir, &nextsubvol);
        }
        
        localwait += now() - waitstart;
        
        // we have a anti message from one neighbour, and 
        // an earlier message from the other neighbour
        if (tt > tn && st->nbs[nextdiffdir]->anti) {
            // critical?!
                goto latereceive;
        }


        if (isinf(tt) || tt >= tspan[tlen-1]) {
            for (i = 0; i < st->nbscnt; i++){
                st->nbs[i]->nextdiff = INFINITY;
            }
            
            // only signal readiness if no threads have already observed
            // finishing condition
            while (block1 > 0) ;
            __sync_fetch_and_add(&block, 1);
            
            while (1) {
                if (block >= nthreads) {
                    for (int dir = 0; dir < st->nbscnt; dir++) {
                        if (cb_peek (st->nbs[dir],  &dummy)) {
                            __sync_fetch_and_add(&block, -1);
                            if(report_level)
                              printf("%d: late receive\n", st->id);
                            goto latereceive;
                        }
                    }
                    __sync_fetch_and_add(&block1, 1);
                    while (1) {
                        if (block >= nthreads && block1 >= nthreads)
                            goto exitloops;
                        if (block < nthreads) {
                            __sync_fetch_and_add(&block1, -1);
                            break;
                        }
                    }
                }
                    
                for (int dir = 0; dir < st->nbscnt; dir++) {
                    if (cb_peek (st->nbs[dir],  &dummy)) {
                        __sync_fetch_and_add(&block, -1);
                        goto latereceive;
                    }
                }
            }
        }

        /* What happens here ? */
    exitloops:
        if (tt >= tspan[it] || isinf(tt)) {
            for (; it < tlen && (tt >= tspan[it] || isinf(tt)); it++)
                memcpy(&U[Ndofs * it], &xx[0], Ndofs * sizeof(int));
                if(report_level)
                  printf("%d: minitime: %.6f, time: %.6f\n", st->id, minitime, st->time);
            if (it >= tlen) {


                for (i = 0; i < st->nbscnt; i++){
                    assert(cb_empty(st->nbs[i]));
                }

                for (i = 0; i < tlen; i++) 
                    for (j = 0; j < Ncells; j++) 
                        for (k = 0; k < Mspecies; k++)
                            globU[(vxs[j]*Mspecies+k)+i*globNdofs] = 
                                U[j*Mspecies+k + Ndofs*i];
                for (int dir = 0; dir < st->nbscnt; dir++)
                    while (cb_pop (st->nbs[dir], &dummy)) {
                        printf("%d: WARNING\n", st->id);
                    }
                // leave simulator loop
                break;
            }
        }

#ifdef FUJIMOTO
        if(LocalGVTFlag>0 && GVTUpdateCnt<GVTCnt) {

        	/* scan for smallest value in all neighbours */
        	min_nbs=INFINITY;
          for (int dir = 0; dir < st->nbscnt; dir++) {
              min_cb=cb_scan_min(st->nbs[dir]);
              if(min_cb<min_nbs)
              	min_nbs=min_cb;
          }
          printf("Thread %d: Minimum neighbour stamp: %f.\n",st->id,min_nbs);
          /* BUG: min_nbs can be negative */
          if(min_nbs<0) min_nbs=INFINITY;
          
          tempmin=MIN(SendMin,min_nbs);
          PEMin[st->id]=MIN(tempmin,tt);
          GVTUpdateCnt++; //computed local min
          
          if(__sync_add_and_fetch(&GVTFlag,-1)==0) {
              /* compute GVT */
              double min=INFINITY;
              for(int i=0;i<nthreads;i++)
                  if(PEMin[i]<min)
                      min=PEMin[i];
              minitime=min;
              printf("Thread %d: Computed new GVT: %f.\n",st->id,minitime);
          }
        }
#endif

        if (reacttimes[0].time <= difftimes[0].time) {

            /* Reaction event. */
            event = 0;
            
            /* a) Extract subvolume and reaction id */
            subvol = reacttimes[0].node / Mreactions;
            re = reacttimes[0].node % Mreactions;
            
            /* Write into the event history list */
#ifdef ROLLBACK
            rb_push_event (event, tt, xxt[subvol], reacttimes[0].node, -1, 0.0);
#endif            
                
            /* b) Update the state of the subvolume subvol and
             * diffusion events. */
            /* Loop over species that are affected by the reaction */
            foreach_spec (subvol, re, tt, tt, 1);
            
            /* c) Update dependent reaction events (including the
             * actual one) */
            total_reactions += updatereacts(subvol, re + Mspecies, tt, 0);
        } else {
            /* Diffusion event. */
            event = 1;
            
            /* Determine parameters */

            for (i = 0; difftimes[0].node >= jcE[i]; i++)
                ;
            subvol = (i - 1) / Mspecies;
            spec = (i - 1) % Mspecies;
            to_vol = irE[difftimes[0].node] / Mspecies;
            to_spec = irE[difftimes[0].node] % Mspecies;

            localpos = lookupLocalVoxel(to_vol);
            dom = lookupVoxel(to_vol);

            assert(xxt[subvol] < tt);

#ifdef ROLLBACK            
            if (dom == localdom && xxt[localpos] > tt)
            {
                antidom = -1;
                rollbackToTime(tt);
                rb_ooo_local++;
                rb_ooo++;
            }
            /* Write into the event history list */
            rb_push_event (event, tt, xxt[subvol], difftimes[0].node, dom == localdom ? -1 : dom, xxt[localpos]);
#endif
            /* Execute the diffusion event (check for negative elements). */
            xx[subvol * Mspecies + spec]--;

            xxt[subvol] = tt;
            
            if (xx[subvol * Mspecies + spec] < 0) {
                errcode = 1;
            }

            /* send message to other domain */
            if (localdom != dom ) {
                msg = (msg_t) {
                    .type = 0,
                    .time = tt,
                    .to_vol = to_vol,
                    .to_spec = to_spec
                };
                while (!cb_push(tst[dom].nbs[to_nbs[dom]], msg));
#ifdef FUJIMOTO
                if(GVTFlag>0 && GVTUpdateCnt<GVTCnt) { 
                	/* scan for minimum in all queues */
                    min_nbs=INFINITY;
                    for (int dir = 0; dir < st->nbscnt; dir++) {
                        min_cb=cb_scan_min(st->nbs[dir]);
                        if(min_cb<min_nbs)
                            min_nbs=min_cb;
                    }
                    SendMin=MIN(SendMin,min_nbs);
                    printf("Thread %d computed SendMin: %f \n",st->id,SendMin);
                }
#endif

            } else if (localdom == dom) {
                /* only update the incoming subvolume if it is internal to the
                 * outgoing subvolume's subdomain */
                xx[localpos * Mspecies + to_spec]++;
                xxt[localpos] = tt;

                /* update reaction rates of incoming subvolume */
                updatereacts(localpos, to_spec, tt, 0);
                /* Update all diffusion events of affected species in
                 * the incoming subvolume */
                updatediffs(localpos, to_spec, tt, 1, 0);
            }

            /* Recalculate the reaction rates using dependency graph G. */
            updatereacts(subvol, spec, tt, 0);
            
            /* Update all diffusion events of affected species in
             * the outgoing subvolume */
            /* In the AEM context decrease the intensity of all
             * events related to non-diagonal row entries of D */
            updatediffs(subvol, spec, tt, -1, 0);


            total_diffusion++; /* counter */
            if (report_level == 3) {
                printf("reactTimes[0]: %f reactTimes[1]: %f \n",reacttimes[0].time,reacttimes[1].time);
                printf("diffTimes[0]: %f diffTimes[1]: %f diffTimes[2]: %f \n",difftimes[0].time,difftimes[1].time,difftimes[2].time);
                printf("-------------------------------------------------------\n");
            }
        }
        if (diff_updated)
            peekdiffs();

        if (errcode) {
            memcpy(&U[Ndofs * it], &xx[0], Ndofs * sizeof(int));
            memset(&U[Ndofs * (it + 1)], -1, Ndofs * (tlen - it - 1) * sizeof(int));
            printf("thread %d, errcode: %d\n", st->id, errcode);
            break;
        }
    }

/*
    for (i = 0; i < pcnt; i++) {
        globwavefront[st->id][i] = wavefront[i];
    }
*/

    int rbcnt = 0;
    GList *rbnext = g_list_first(rblist);
    while (rbnext)
    {
        rbcnt++;
        rbnext = g_list_next(rbnext);
    }




__sync_fetch_and_add(&glob_pdex, pdex);
__sync_fetch_and_add(&glob_pdupdated, pdupdated);
__sync_fetch_and_add(&glob_pdupdatedinf, pdupdatedinf);

    if(report_level)
printf("%d: rblist left: %d, wt: %.1f\n", st->id, rbcnt,(double)localwait/(double)CPUFREQ);
//        printf ("pd: %ld %ld %ld\n",pdex, pdupdated, pdupdatedinf);
//    }
    
//    printf("%d: rb_newanti: %d, rb_anti: %d, rb_ooo: %d, rb_ooo_local: %d, cancelled: %d, rollbacked: %d\n", st->id, rb_newanti, rb_anti, rb_ooo, rb_ooo_local, cancelled_msgs, rollbacked_events);

__sync_fetch_and_add(&glob_cancelled_msgs, cancelled_msgs);
__sync_fetch_and_add(&glob_rollbacked_events, rollbacked_events);
__sync_fetch_and_add(&glob_rb_anti, rb_anti);
__sync_fetch_and_add(&glob_rb_ooo, rb_ooo);


    if (report_level == 2) {

        printf("%d: rb_newanti: %d, rb_anti: %d, rb_ooo: %d, rb_ooo_local: %d\n", st->id, rb_newanti, rb_anti, rb_ooo, rb_ooo_local);



//        printf("%d waitcnt: %d\n", st->id, waitcnt);
//        printf("%2d localwait: %llu\n", st->id, localwait);

//        printf("%2d depcnt: %llu\n", st->id, depcnt);
//        printf("Updates: %lu Wake-ups: %lu\n", updates, wakeups);
//        printf("React2React: %lu React2diff: %lu\n", react2react, react2diff);
//        printf("Diff2diff: %lu Diff2react: %lu\n", diff2diff, diff2react);

    }
  
//    free(reactHeap);
//    free(reactNode);
//    free(reactTimes);
//    free(reactInf);
    
//    free(diffHeap);
//    free(diffNode);
//    free(diffTimes);
//    free(diffInf);
//    free(jcE);
//    free(irE);
        
//	free(sdrate);
//	free(rrate);
//	free(xx);
    st->wait = localwait;
    return NULL;
}


/* b) Update the state of the subvolume subvol and
 * diffusion events. */
/* Loop over species that are affected by the reaction */
void
foreach_spec (int subvol, int reaction, double tt, double xxtt, int modifier) 
{
    for (int i = jcN[reaction]; i < jcN[reaction + 1]; i++) {
        /* Update the state according to the stochiometry matrix */
        xx[subvol * Mspecies + irN[i]] += prN[i] * modifier;
        xxt[subvol] = xxtt;

        assert (xx[subvol * Mspecies + irN[i]] >= 0);
        
        /* Update all dependent diffusion events */
        updatediffs(subvol, irN[i], tt, prN[i] * modifier, 0);
    }
}

int 
updatediffs (const int subvol,const int spec,const double tt,
             const int multi, const int rb) 
{
    int cnt, n = 0;
    double oldrate;
    
    /* Update all diffusion events of affected species in
     * the outgoing subvolume */
    /* In the AEM context decrease the intensity of all
     * events related to non-diagonal row entries of D */
    for (cnt = jcE[subvol * Mspecies + spec];
         cnt < jcE[subvol * Mspecies + spec + 1]; cnt++) {
        
        oldrate = diffnodes[cnt].rate;
        /* Update the diffusion rate in respect to the
         * affected species & subvolume */
        diffnodes[cnt].rate += diag[cnt] * multi;
        
        roundoff(&diffnodes[cnt].rate);
        //update of the affected node
        if (rb && cnt == (subvol * Mspecies + spec)) { 
            difftimes[diffnodes[cnt].heapidx].time = tt;
            /* rewind the RNG state for the main rollbacked event */
#ifdef RLCG            
            printf("diff seed %d before revert: %lu.\n",cnt,diffnodes[cnt].seed);
            rng_prev(&diffnodes[cnt].seed);
            printf("diff seed %d after revert: %lu.\n",cnt,diffnodes[cnt].seed);
#endif
        } else {
            calcTimes(&difftimes[diffnodes[cnt].heapidx].time, 
                      &diffnodes[cnt].inftime, tt, 
                      oldrate, diffnodes[cnt].rate, &diffnodes[cnt].seed);
        }
        /* update heap */
        update(diffnodes[cnt].heapidx, difftimes, diffnodes, 
               diffHeapSize);
        diff_updated |= diffnodes[cnt].internode;
//        diff_updated |= diffnodes[cnt].internode && (diffnodes[cnt].heapidx < DIFF_SCANLENGTH);
        n++;

    }
    return n;
}



/* index : re + Mspecies or spec */
int
updatereacts (const int subvol,const int index, const double tt, const int rb)
{
    int i, j, n = 0;
    double oldrate;
    
    /* Update dependent reaction events (including the
     * actual one) */
    /* Recalculate the reaction rates using dependency graph G. */
    for (i = jcG[index]; i < jcG[index + 1]; i++) {
        j = irG[i];
                    
        oldrate = reactnodes[subvol * Mreactions + j].rate;
        assert(subvol*Mreactions + j < Mreactions * Ncells);
        
        reactnodes[subvol * Mreactions + j].rate = 
            (*rfun[j])(&xx[subvol * Mspecies], tt, vol[vxs[subvol]], 
                       &data[vxs[subvol] * dsize], sd[vxs[subvol]]);
        
        if (rb && j == (index-Mspecies)) { // update of the affected node
            reacttimes[reactnodes[subvol * Mreactions + j].heapidx].time = tt;
            /* rewind the RNG state for the main rollbacked event */
#ifdef RLCG
            printf("react seed %d before revert: %llu.\n",(subvol * Mreactions + j),reactnodes[subvol * Mreactions + j].seed);
            rng_prev(&reactnodes[subvol * Mreactions + j].seed);
            printf("react seed %d after revert: %llu.\n",(subvol * Mreactions + j),reactnodes[subvol * Mreactions + j].seed);
#endif
        } else { // Update of dependent nodes
            calcTimes(&reacttimes[reactnodes[subvol * Mreactions + j].heapidx].time, 
                      &reactnodes[subvol * Mreactions + j].inftime, tt, 
                      oldrate, reactnodes[subvol * Mreactions + j].rate,
                      &reactnodes[subvol * Mreactions + j].seed);
        }

        /* update heap */
        update(reactnodes[subvol * Mreactions + j].heapidx, 
               reacttimes, reactnodes, reactHeapSize);

        n++;
    }
    return n;
}





void 
peekdiffs() 
{
    int dom, to_vol;
    int updated[nthreads];
    int updatedcnt = 0;
    
    memset(updated, 0, sizeof updated);
    pdex++;
    for (int i = 0; i < diffHeapSize && i < DIFF_SCANLENGTH && updatedcnt < to_nbs_cnt; i++) {
        to_vol = irE[difftimes[i].node] / Mspecies;
        dom = lookupVoxel(to_vol);
        if (dom != st->id && updated[dom] == 0) {
            updated[dom] = 1;
            if (tst[dom].nbs[to_nbs[dom]]->nextdiff != difftimes[i].time) {
                pdupdated++;
                tst[dom].nbs[to_nbs[dom]]->nextdiff = difftimes[i].time;
                tst[dom].nbs[to_nbs[dom]]->nextsubvol = to_vol;
                updatedcnt++;
            }
        }
    }
    for (int i = 0; i < nthreads; i++) {
        if (!updated[i] && to_nbs[i] >= 0) {
            if (tst[i].nbs[to_nbs[i]]->nextdiff != INFINITY)
            {                
                tst[i].nbs[to_nbs[i]]->nextdiff = INFINITY;
                pdupdatedinf++;
            }
            
        }
    }
}

/* ---------------------------------------------------------------------- */

void rlb_dbg_save() {
	int i;

  //malloc
  xx_dbg = PAD_ALLOC(Ndofs * sizeof(int));
  reacttimes_dbg = PAD_ALLOC(reactHeapSize * sizeof(evtime_t));
  difftimes_dbg = PAD_ALLOC(diffHeapSize * sizeof(evtime_t));
  reactnodes_dbg = PAD_ALLOC(reactHeapSize * sizeof(evnode_t));
  diffnodes_dbg = PAD_ALLOC(diffHeapSize * sizeof(evnode_t));
  
  memcpy(xx_dbg,xx,Ndofs * sizeof(int));
  memcpy(reacttimes_dbg,reacttimes,reactHeapSize * sizeof(evtime_t));
  memcpy(difftimes_dbg,difftimes,diffHeapSize * sizeof(evtime_t));
  memcpy(reactnodes_dbg,reactnodes,reactHeapSize * sizeof(evnode_t));
  memcpy(diffnodes_dbg,diffnodes,diffHeapSize * sizeof(evnode_t));
  
  // more? msg queue, rlb_queue,...
}

/* ---------------------------------------------------------------------- */

int rlb_dbg_compare(){
	int i,ok=1;
  
  printf("****************************************************\n");

	// compare state
  for (i = 0; i < Ndofs; i++) {
  	if(xx_dbg[i]!=xx[i]) {
  		printf("state of node %d differs (%d vs %d).\n",i,xx_dbg[i],xx[i]);
  		ok=0;
  	}
  }

	// compare reacttimes
  for (i = 0; i < reactHeapSize; i++) {
  	if(fabs(reacttimes_dbg[i].time-reacttimes[i].time)>DBL_EPSILON){
  		printf("reacttimes.times %d differs (%f vs %f).\n",i,reacttimes_dbg[i].time,reacttimes[i].time);
  		ok=0;
  	}
  	if(reacttimes_dbg[i].node!=reacttimes[i].node){
  		printf("reacttimes.node %d differs (%d vs %d).\n",i,reacttimes_dbg[i].node,reacttimes[i].node);
  		ok=0;
  	}
  }

  // compare difftimes
  for (i = 0; i < diffHeapSize; i++) {
  	if(fabs(difftimes_dbg[i].time-difftimes[i].time)>DBL_EPSILON){
  		printf("difftimes.times %d differs (%f vs %f).\n",i,difftimes_dbg[i].time,difftimes[i].time);
  		ok=0;
  	}
  	if(difftimes_dbg[i].node!=difftimes[i].node){
		printf("difftimes.node %d differs (%d vs %d).\n",i,difftimes_dbg[i].node,difftimes[i].node);
		ok=0;
  	}
  }

	// compare reactnodes
  for (i = 0; i < reactHeapSize; i++) {
  	if(reactnodes_dbg[i].heapidx!=reactnodes[i].heapidx){
  			printf("reactnodes.heapidx %d differs (%d vs %d).\n",i,reactnodes_dbg[i].heapidx,reactnodes[i].heapidx);
  			ok=0;
  	}
  	if(reactnodes_dbg[i].rate!=reactnodes[i].rate){
  			printf("reactnodes.rate %d differs (%f vs %f).\n",i,reactnodes_dbg[i].rate,reactnodes[i].rate);
  			ok=0;
  	}
  	if(reactnodes_dbg[i].inftime!=reactnodes[i].inftime){
  			printf("reactnodes.inftime %d differs (%f vs %f).\n",i,reactnodes_dbg[i].inftime,reactnodes[i].inftime);
  			ok=0;
  	}
#ifdef RLCG
  	if(reactnodes_dbg[reacttimes_dbg[i].node].seed!=reactnodes[reacttimes[i].node].seed){
  			printf("reactnodes.seed %d differs (%llu vs %llu).\n",i,reactnodes_dbg[i].seed,reactnodes[i].seed);
  			ok=0;
  	}
#else
    for(int j=0; j<3; j++) {
        if(reactnodes_dbg[i].seed[j]!=reactnodes[i].seed[j]){
                printf("reactnodes[%d].seed[%d] differs (%d vs %d).\n",i,j,reactnodes_dbg[i].seed[j],reactnodes[i].seed[j]);
                ok=0;
        }
    }
#endif
  }

  // compare diffnodes
  for (i = 0; i < diffHeapSize; i++) {
  	if(diffnodes_dbg[i].heapidx!=diffnodes[i].heapidx){
			printf("diffnodes.heapidx %d differs (%d vs %d).\n",i,diffnodes_dbg[i].heapidx,diffnodes[i].heapidx);
			ok=0;
  	}
  	if(diffnodes_dbg[i].rate!=diffnodes[i].rate){
  		printf("diffnodes.rate %d differs (%f vs %f).\n",i,diffnodes_dbg[i].rate,diffnodes[i].rate);
  		ok=0;
  	}
  	if(diffnodes_dbg[i].inftime!=diffnodes[i].inftime){
  		printf("diffnodes.inftime %d differs (%f vs %f).\n",i,diffnodes_dbg[i].inftime,diffnodes[i].inftime);
  		ok=0;
  	}
#ifdef RLCG    
  	if(diffnodes_dbg[i].seed!=diffnodes[i].seed){
  		printf("diffnodes.seed %d differs (%llu vs %llu).\n",i,diffnodes_dbg[i].seed,diffnodes[i].seed);
  		ok=0;
  	} 
 #else
    for(int j=0; j<3; j++) {
        if(diffnodes_dbg[i].seed[j]!=diffnodes[i].seed[j]){
            printf("diffnodes[%d].seed[%d] differs (%d vs %d).\n",i,j,diffnodes_dbg[i].seed[j],diffnodes[i].seed[j]);
            ok=0;
        }
  	}
 #endif
  }
    
  printf("****************************************************\n");

	return(ok);
}
