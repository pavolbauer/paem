/* aemore.c - Core AEM solver. */

/* P. Bauer and S. Engblom 2012-05-10 */

// three lines to enable pinning of threads (together with _GNU_SOURCE)
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

#include "propensities.h"
#include "aem.h"
#include "../nsm/binheap.h"
#include "report.h"
#include "spsc_cb.h"
#include <glib.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define ROUNDOFF 1

#define NBS 2

#define BLOCKDIFF

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
    int type;  // type of event
    int id;    // event ID
} rlb_t;

/* Global variables for debug statistics */
unsigned long updates = 0;
unsigned long wakeups = 0;
unsigned long react2diff = 0;
unsigned long react2react = 0;
unsigned long diff2react = 0;
unsigned long diff2diff = 0;

/* pinning of threads */
/* straight core allocation on sandy (i.e., threads 0 - 15 on socket 1,
 *  threads 16-31 on socket 2, usw.)
 */
#define topo32ht(i) ( ((i)/2)/8 + (4*((i)/2) % 32) + ((i)%2)*32 )

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

/* thread state */
typedef struct tst_s {
    pthread_t t; // pthread specific (for join)
    int id; 
    long cnt; // count something

    char pad0[128];
    double localtime;
    double nextdiff[NBS];
    int nextdiffvol[NBS];
    char pad[128];
    cb_t **nbs; // neighbours?
} tst_t;

static tst_t *tst;

volatile int block = 0;


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


/* Thread specific global variables */

__thread double *rrate, *isdrate2, *xxt;
__thread double *reactTimes, *diffTimes, *reactInf, *diffInf;

__thread int *reactNode, *reactHeap, *diffNode, *diffHeap, *xx, *U;
__thread int *jcE, *irE, *jvec, *vxs, localdom;

__thread size_t Ncells, Ndofs, diffHeapSize, reactHeapSize;
__thread unsigned short rng[3];

/* Rollback list */
__thread GList *rblist;

void aem_core(const size_t *_irD,const size_t *_jcD,const double *_prD,
              const int *_u0,
              const size_t *_irN,const size_t *_jcN,const int *_prN,
              const size_t *_irG,const size_t *_jcG,
              const double *_tspan,const size_t _tlen,
              int *_U,
              const double *_vol,const double *_data,const int *_sd,
              const size_t _Ncells,
              const size_t _Mspecies,const size_t _Mreactions,
              const size_t _dsize,int _report_level, int _nthreads)

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

  const size_t Ndofs = globNcells * Mspecies;

  rfun = ALLOC_propensities();

  nthreads = _nthreads;

  ReportFun report = NULL;
  if (report_level) report = &reportFun1;


  /* create subdomain lookup */
  sd_length=calloc(nthreads,sizeof(int));
  for(i = 0; i < globNcells; i++)
      sd_length[sd[i]-1]++;

  domainLookup = malloc(nthreads*sizeof(int*));

  int cnt = 0;  
  g2lLookup = (int*) malloc(globNcells*sizeof(int));
  for (i = 0; i < nthreads; i++)
  {
      domainLookup[i] = (int*) malloc(sd_length[i]*sizeof(int));

      for(cnt = 0, j = 0; j < globNcells; j++) {
          if( (sd[j]-1) == i) {
              domainLookup[i][cnt] = j;
              g2lLookup[j] = cnt;
              cnt++;
          }
      }
  }
    
  tst = calloc(nthreads, sizeof (tst_t));
    
    for (int i = 0; i < nthreads; i++)
    {
        tst[i].id = i;
        tst[i].nextdiff[0] = 0.0;
        tst[i].nextdiff[1] = 0.0;

        pthread_create (&tst[i].t, NULL, run, (void *)&tst[i]);
    }


    while (block < nthreads)
        ;

    block = 0;
    

    for (int i = 0; i < nthreads; i++)
    {
        pthread_join (tst[i].t, NULL);
    }

/*
   for (i = 0; i < tlen; i++) {
       printf("%f ", tspan[i]);
       for (j = 0; j < globNdofs; j++) {
           printf("%d ", _U[i*globNdofs + j]);
       }
       printf("\n");
   }
*/

    int ooo = 0;
    int msgcnt = 0;
    
    for (int i = 0; i < nthreads; i++) {
        printf("%d:\n", i);
        printf("ooo: %d\n", tst[i].cnt);
        ooo += tst[i].cnt;
        for (int d = 0; d < NBS; d++) {
            printf("%d ", tst[i].nbs[d]->cnt);
            msgcnt += tst[i].nbs[d]->cnt;
        }
        printf("\n");

    }
    printf("total, ooo: %d, msg: %d\n", ooo, msgcnt);

    free (tst);
    FREE_propensities(rfun);
}


int *
lookupSubDomain(int subdom)
{
    return domainLookup[subdom];
}


int 
lookupVoxel(int voxel)
{
    return sd[voxel] - 1;
}

/* the return value is only valid when 
 * the calling thread is the owner of the 
 * voxel. */
int
lookupLocalVoxel(int global_voxel) {
    return g2lLookup[global_voxel];
}



void
broadcastTime(int id, double tt)
{
    tst[id].localtime = tt;
}

void
broadcastDiff(int id, double diff, int dom)
{
    int dir = dom - id + (dom < id);
    if (tst[id].nextdiff[dir] != diff)
        tst[id].nextdiff[dir] = diff;
}

double
readNeighbourTimes(int id) {
    if (id == 0) return tst[id + 1].localtime;
    if (id == nthreads - 1) return tst[id - 1].localtime;
    return MIN(tst[id - 1].localtime, tst[id + 1].localtime);
}

double
readNeighbourDiffTimes(int id, int *vol) {
    if (id == 0) {
        *vol = tst[id + 1].nextdiffvol[0];
        return tst[id + 1].nextdiff[0];
    }
    if (id == nthreads - 1) {
        *vol = tst[id - 1].nextdiffvol[1];
        return tst[id - 1].nextdiff[1];
    }
    if (tst[id - 1].nextdiff[1] < tst[id + 1].nextdiff[0]) {
        *vol = tst[id - 1].nextdiffvol[1];
    } else {
        *vol = tst[id + 1].nextdiffvol[0];
    }
    return MIN(tst[id - 1].nextdiff[1], tst[id + 1].nextdiff[0]);
}



int
areNeighbours3d(int sd0, int sd1) {
    return (sd0 == sd1) ||
        (sd0+1 == sd1) || (sd0+4 == sd1) || (sd0+8 == sd1) ||
        (sd0-1 == sd1) || (sd0-4 == sd1) || (sd0-8 == sd1);
}

int
areNeighbours1d(int sd0, int sd1) {
    return (sd0 == sd1) || (sd0+1 == sd1) || (sd0-1 == sd1);
}


cb_t * 
receiver(int id, int nbor) {
    int direction = abs(nbor-id)/4 + 1 - ((nbor - id)>=0);
    return tst[nbor].nbs[(direction+(NBS/2)) % (NBS)];
}

cb_t * 
receiver1d(int id, int nbor) {
    int direction = nbor - id + ((nbor - id) < 0);
    return tst[nbor].nbs[1 - direction];
}


static void 
idle_wait(unsigned int ticks) {
    volatile int x;
    for (x = 0; x < ticks; x++);
}


void 
calcTimesBack(double* time, double* infTime, double tt, double old_rate,
               double new_rate, seeds seed) {

    double oldtime = time[0];
    printf("Inftime: %f. new_rate: %f \n",infTime[0],new_rate);
    if (isinf(oldtime)) {
        
        if (infTime[0] == 0) { // Waking up first time
            if (new_rate > 0.0) {
                time[0] = -log(1.0 - erand48(seed)) / new_rate + tt;
                printf("hehe");
            }
            
        } else { // Waking up the 2nd..nth time
            if (new_rate > 0.0) {
                time[0] = tt + (infTime[0] / new_rate);
                printf("haha (using InfTime=%f).",infTime[0]);
            }
        }
        wakeups++;
    } else {
        if (new_rate >= DBL_MIN) {
            if (time[0] == tt) {// Regular update of current event 
                time[0] = -log(1.0 - erand48(seed)) / new_rate + tt;
                printf("huhu.");
            } else { // Regular update of dependent events (rescaling)
                time[0] = ((old_rate / new_rate) * (oldtime - tt)) + tt;
                printf("here.");
            }
            updates++;

        } else { // Next event time set to infinity
            infTime[0] = (oldtime - tt) * old_rate;
            time[0] = INFINITY;
            printf("huhu");
        }
    }
    printf("newly calculated difftime: %f.\n",time[0]);
}

/*-----------------------------------------------------------------------*/

msg_t 
rollback_event(rlb_t* event){  
  int i,j,k,errcode,cnt;
  double oldrate;
  msg_t msg;
  
  int subvol,re,spec,to_vol,to_spec;
  int localpos,dom;
  
  if(event->type==0) {
    subvol = event->id / Mreactions;
    re = event->id % Mreactions;
    
    for (i = jcN[re]; i < jcN[re + 1]; i++) { // foreach species?
      assert(irN[i] < Mspecies);
      
      /* Update the state according to the stochiometry matrix */
      xx[subvol * Mspecies + irN[i]] -= prN[i];
      
      //xxt[subvol] = event.tt;
      
      if (xx[subvol * Mspecies + irN[i]] < 0)
          errcode = 1;

      /* Update all dependent diffusion events */
      for (cnt = jcE[subvol * Mspecies + irN[i]];
           cnt < jcE[subvol * Mspecies + irN[i] + 1]; cnt++) {

          oldrate = isdrate2[cnt];

          // Update the diffusion rate in respect to the
          // affected species & subvolume 
          isdrate2[cnt] -= prD[jvec[cnt]] * prN[i];

          if (fabs(isdrate2[cnt]) < 1.0E-12) {
              isdrate2[cnt] = 0.0;
          }
          
          calcTimes(&diffTimes[diffHeap[cnt]], &diffInf[cnt], 
                    event->tt, oldrate, isdrate2[cnt], rng);

          update(diffHeap[cnt], diffTimes, diffNode, 
                 diffHeap, diffHeapSize);
          
        }
      }
    
      /* c) Update dependent reaction events (including the
       * actual one) */
      
      for (i = jcG[Mspecies + re]; i < jcG[Mspecies + re + 1]; i++) {
          j = irG[i];

          oldrate = rrate[subvol * Mreactions + j];

          assert(subvol*Mreactions + j < Mreactions * Ncells);

          rrate[subvol * Mreactions + j] = 
              (*rfun[j])(&xx[subvol * Mspecies], event->tt, vol[vxs[subvol]], 
                         &data[vxs[subvol] * dsize], sd[vxs[subvol]]);
          
          //assert(rrate[subvol * Mreactions + j]>0.0);
	  printf("new rate: %f.\n",rrate[subvol * Mreactions + j]);
          
          if(j==re) { // update of the affected node
            reactTimes[reactHeap[subvol * Mreactions + j]]=event->tt;
            printf("yes.");
          } else { // Update of dependent nodes
            calcTimes(&reactTimes[reactHeap[subvol * Mreactions + j]],
              &reactInf[subvol * Mreactions + j], event->tt, oldrate,
              rrate[subvol * Mreactions + j],rng);
          }

          update(reactHeap[subvol * Mreactions + j], 
                 reactTimes, reactNode,reactHeap, reactHeapSize);
      }
  }
  else{
    for (i = 0; event->id >= jcE[i]; i++)
        ;
    
    subvol = (i - 1) / Mspecies;
    spec = (i - 1) % Mspecies;
    to_vol = irE[event->id] / Mspecies;
    to_spec = irE[event->id] % Mspecies;
    
    /* revert local diffusion */
    xx[subvol * Mspecies + spec]++;
    assert(xx[subvol * Mspecies + spec]>=0);
    //xxt[subvol] = tt;

    if (xx[subvol * Mspecies + spec] < 0) {
        errcode = 1;
    }
    
    localpos = lookupLocalVoxel(to_vol);
    dom = lookupVoxel(to_vol);

    /* message will serve to determine rollback cascades */
    if (localdom != dom) {
        msg = (msg_t) {
            .time = event->tt,
            .to_vol = to_vol,
            .to_spec = to_spec
        };
    }
    else if (localdom == dom) {
        //assert(localpos < sd_length[st->id]);
        xx[localpos * Mspecies + to_spec]--;
	assert(xx[localpos * Mspecies + to_spec]>=0);
        //xxt[localpos] = event->tt;
    }
    else {
        // bugger
        assert(0);
    }

    /* Recalculate the reaction rates using dependency graph G. */
    for (i = jcG[spec]; i < jcG[spec + 1]; i++) {
        j = irG[i];

        oldrate = rrate[subvol * Mreactions + j];
        
        //assert(subvol*Mreactions + j < Mreactions * Ncells);
 
        rrate[subvol * Mreactions + j] = 
            (*rfun[j])(&xx[subvol * Mspecies], event->tt, vol[vxs[subvol]],                                
                &data[vxs[subvol] * dsize], sd[vxs[subvol]]);


        calcTimesBack(&reactTimes[reactHeap[subvol * Mreactions + j]],
                  &reactInf[subvol * Mreactions + j], event->tt, oldrate,
                  rrate[subvol * Mreactions + j], rng);

        update(reactHeap[subvol * Mreactions + j], reactTimes, 
               reactNode, reactHeap, reactHeapSize);
         
 
        if (dom == localdom) {
            //assert(localpos < sd_length[st->id]);
            oldrate = rrate[localpos * Mreactions + j];
            assert(localpos * Mreactions + j < Mreactions * Ncells);
            rrate[localpos * Mreactions + j] = 
                (*rfun[j])(&xx[localpos * Mspecies],event->tt,vol[to_vol], 
                           &data[to_vol * dsize], sd[to_vol]);

            calcTimesBack(&reactTimes[reactHeap[localpos * Mreactions + j]],
                      &reactInf[localpos * Mreactions + j], event->tt, 
                      oldrate, rrate[localpos * Mreactions + j], rng);

            update(reactHeap[localpos * Mreactions + j], reactTimes, 
                   reactNode, reactHeap, reactHeapSize);
        }
    } 

    /* Update all diffusion events of affected species in
     * the outgoing subvolume */
    /* In the AEM context decrease the intensity of all
     * events related to non-diagonal row entries of D */
    for (cnt = jcE[subvol * Mspecies + spec];
         cnt < jcE[subvol * Mspecies + spec + 1]; cnt++) {
        assert(cnt < diffHeapSize);

        oldrate = isdrate2[cnt];
        isdrate2[cnt] -= prD[jvec[cnt]];
#ifdef ROUNDOFF  
        if (fabs(isdrate2[cnt]) < 1.0E-12) {
            isdrate2[cnt] = 0.0;
        }
#endif
        if(cnt==(subvol * Mspecies + spec)) { // update of the affected node
            diffTimes[diffHeap[cnt]]=event->tt;
            printf("yes."); 
        } else {
        calcTimes(&diffTimes[diffHeap[cnt]], &diffInf[cnt], event->tt, 
                  oldrate, isdrate2[cnt], rng);
        }
        update(diffHeap[cnt], diffTimes, diffNode, 
               diffHeap, diffHeapSize);
    }

    /* only update the incoming subvolume if it is internal to the
     * outgoing subvolume's subdomain */
    if (dom == localdom) {
        /* Update all diffusion events of affected species in
         * the incoming subvolume */
        //assert(localpos < sd_length[st->id]);
        for (cnt = jcE[localpos * Mspecies + to_spec];
             cnt < jcE[localpos * Mspecies + to_spec + 1]; cnt++) {
            assert(cnt < diffHeapSize);
            oldrate = isdrate2[cnt];
            isdrate2[cnt] += prD[jvec[cnt]];
#ifdef ROUNDOFF
            if (fabs(isdrate2[cnt]) < 1.0E-12) {
                isdrate2[cnt] = 0.0;
            }
#endif
            if(cnt==(localpos * Mspecies + to_spec)) {
              diffTimes[diffHeap[cnt]]=event->tt;
              printf("yes."); 
            } else {
            calcTimes(&diffTimes[diffHeap[cnt]], &diffInf[cnt], 
                      event->tt, oldrate, isdrate2[cnt], rng);
            }
            update(diffHeap[cnt], diffTimes, diffNode, diffHeap, diffHeapSize);
        }
    }
  }
  
  printf("Rollbacked type %d at tt=%f \n",event->type,event->tt);
  printf("reactTimes[0]: %f reactTimes[1]: %f \n",reactTimes[0],reactTimes[1]);
  printf("diffTimes[0]: %f diffTimes[1]: %f diffTimes[2]: %f \n",diffTimes[0],diffTimes[1],diffTimes[2]);
  printf("-------------------------------------------------------\n");
  
  /* Return information regarding rollbacked diffusion into other domains */
  return msg;
}


/* Rollbacks domain until time>rolltt
 * Returns the last valid timestamp<rolltt
*/

double
testHistoryOut(double rolltt) {
  
  rlb_t* event; 
  GList *l = g_list_last(rblist);
  msg_t msg;
  int found=0;
  double tt;
  
  while (l != NULL && found==0)
  {
    GList *next = g_list_previous(l);
    event=l->data;
    tt=event->tt;
    if(event->tt>=rolltt) { 
      msg=rollback_event(event);
      free(l->data);
      rblist = g_list_delete_link (rblist, l);
      l = next;
    } else {
     found=1; 
    }
  }
  printf("Last valid tt=%f.\n",tt);
  return tt;
  
  // return rollback'd tt
  // return last outgoing diffusion time
}


void *
run (void *_args)
{
    double tt = tspan[0];
    double oldtime, oldrate;
    
    int waitcnt = 0;    
    
    int i, j, re, cnt;
    
    int event, subvol, spec, to_vol, to_spec;
    
    short errcode = 0;
    long int total_diffusion = 0, total_reactions = 0;

    size_t it = 0;

    tst_t *st = (tst_t *)_args;
    
    rlb_t* rb;

    pin(gettid(), topo32ht(st->id));

    // set up msgs channels to neighbours
    st->nbs = malloc( sizeof(cb_t *) * NBS);
    for (i = 0; i < NBS; i++) {
        st->nbs[i] = cb_init(512);
    }

    //initialize thread local rng
    for (i = 0; i < 3; i++)
        rng[i] = rand() % 65535;


    localdom = st->id;
    vxs = lookupSubDomain(localdom);
    
    Ncells = sd_length[st->id];
    Ndofs = Ncells*Mspecies;

    U = calloc (tlen * Ndofs,sizeof(int));

    xxt = calloc(Ncells, sizeof(double));
    
    /* Set xx to the initial state. */
    xx = (int *) malloc(Ndofs * sizeof(int));

    for (i = 0; i < Ncells; i++) 
        for (j = 0; j < Mspecies; j++) 
        {
            xx[(i*Mspecies)+j] = u0[vxs[i]*Mspecies+j]; 
        }

    /* Create reaction rate matrix (Mreactions X Ncells) and total rate
       vector. In rrate we store all propensities for chemical rections,
       and in srrate the sum of propensities in every subvolume. */
    rrate = (double *) malloc(Mreactions * Ncells * sizeof(double));

    /* Calculate the propensity for every reaction and every
       subvolume. Store the sum of the reaction intensities in each
       subvolume in srrate. */
    for (i = 0; i < Ncells; i++) {
        for (j = 0; j < Mreactions; j++) {
            rrate[i * Mreactions + j] = 
                (*rfun[j])(&xx[i * Mspecies], tt, vol[vxs[i]],
                           &data[vxs[i] * dsize], sd[vxs[i]]);
        }
    }

    /* Binary heap storing all reaction events */
    reactHeapSize = Ncells * Mreactions;
    
    /* Create Reaction-heap. */
    reactTimes = (double *) malloc(reactHeapSize * sizeof(double));
    reactNode = (int *) malloc(reactHeapSize * sizeof(int));
    reactHeap = (int *) malloc(reactHeapSize * sizeof(int));

    /* Create Array to store pre-infinity times */
    reactInf = (double *) malloc(reactHeapSize * sizeof(double));
    for (i = 0; i < reactHeapSize; i++)
        reactInf[i] = 0;

    /* Calculate times to next reaction event */
    for (i = 0; i < reactHeapSize; i++) {
        reactTimes[i] = -log(1.0 - erand48(rng)) / rrate[i] + tspan[0];
        
        if (reactTimes[i] <= 0.0)
            reactTimes[i] = INFINITY;
        
        reactHeap[i] = reactNode[i] = i;
    }

    /* Initialize reaction heap */
    initialize_heap(reactTimes, reactNode, reactHeap, reactHeapSize);
        
    int my_i, global_i;
    int myNnz = 0;
    int k;
    
    
    for (i = 0; i < Ncells; i++) {
        for (k = 0; k<Mspecies; k++) {
            /* ndof_i -> mapped into global Ndofs adress */
            my_i=(i*Mspecies)+k;
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
    
    diffTimes = (double *) malloc(diffHeapSize * sizeof(double));
    diffNode = (int *) malloc(diffHeapSize * sizeof(int));
    diffHeap = (int *) malloc(diffHeapSize * sizeof(int));
    isdrate2 = (double *) malloc(diffHeapSize * sizeof(double));

    /* Create Array to store pre-infinity times */
    diffInf = (double *) calloc(diffHeapSize,sizeof(double));
    
    /* Creating new sparse matrix for diffusion events */
    jcE = (int*) malloc((Ndofs + 1) * sizeof(int));
    irE = (int*) malloc(diffHeapSize * sizeof(int));
    jvec = (int*) malloc(diffHeapSize * sizeof(int));

    /* Queue all non-diagonal entries of D as diffusion events */
    for (cnt = 0, i = 0; i < Ncells; i++) {
        for(k=0; k<Mspecies; k++) {
            my_i=(i*Mspecies)+k;
            global_i=(vxs[i]*Mspecies)+k;
            jcE[my_i] = cnt;
            for (j = jcD[global_i]; j < jcD[global_i + 1]; j++)
                if (irD[j] != global_i) { // if not diagonal
                    // check that only neighbours comm.
                    assert(areNeighbours1d(st->id, lookupVoxel(irD[j]/Mspecies)));
                    diffNode[cnt] = diffHeap[cnt] = cnt;
                    isdrate2[cnt] = xx[i*Mspecies+k] * prD[j];
                    diffTimes[cnt] = -log(1.0 - erand48(rng)) / isdrate2[cnt] + tspan[0];
                    if (diffTimes[cnt] <= 0.0)
                        diffTimes[cnt] = INFINITY;
                    irE[cnt] = irD[j];
                    jvec[cnt] = j;
                    cnt++;
                }
        }
    }
    jcE[my_i+1] = cnt;

    /* Initialize diffusion heap */
    initialize_heap(diffTimes, diffNode, diffHeap, diffHeapSize);
    
    /* Fix to run reaction-only tests */
    if (diffHeapSize == 0) {
        diffTimes = (double *) malloc(sizeof(double));
        diffTimes[0] = INFINITY;
	}
    /* Fix to run diffusion-only tests */
    if (reactHeapSize == 0) {
        reactTimes = (double *) malloc(sizeof(double));
        reactTimes[0] = INFINITY;
    }
    
    if (report_level == 2) {
        printf("React-Heap Size: %d\n", (int) reactHeapSize);
        printf("Diff-Heap Size: %d\n", (int) diffHeapSize);
    }
    
    int dom;
    int lf = 0, rf = 0;
    if (st->id == 0) lf = 1;
    if (st->id == nthreads - 1) rf = 1;
    for (i = 0; i < diffHeapSize; i++) {
        if (lf && rf) break;
        to_vol = irE[diffNode[i]] / Mspecies;
        dom = lookupVoxel(to_vol);
        if (dom == st->id + 1 && !rf) {
            rf = 1;
            if (st->nextdiff[1] != diffTimes[i]) {
                st->nextdiff[1] = diffTimes[i];
                st->nextdiffvol[1] = to_vol;
            }
        } else if (dom == st->id - 1 && !lf) {
            lf = 1;
            if (st->nextdiff[0] != diffTimes[i]) {
                st->nextdiff[0] = diffTimes[i];
                st->nextdiffvol[0] = to_vol;
            }
        }
    }
    if (lf != 1) {
        st->nextdiff[0] = INFINITY;
    } else if (rf != 1) {
        st->nextdiff[0] = INFINITY;
    }


    /* barrier, guarantee that all threads are initialized */
    __sync_fetch_and_add(&block, 1);
    while (block);
    
    double tn;
    msg_t msgp;
    msg_t msg;

    int next;

    int nextsubvol = -1;
    int nextdiffvol;
    double told;
    /* Main loop. */
    for (;;) {
        told = tt;
    latereceive:

/*
        for (int dir = 0; dir < NBS; dir++) {
            while(cb_pop (st->nbs[dir], &msgp)) {
                if (msgp.time < told) {
                    st->cnt++;
                }
            }
        
        }
*/
        /* Calculate tt with diff-time and react-time */
        tt = MIN(reactTimes[0],diffTimes[0]);
//        if (reactTimes[0] < diffTimes[0]) {
//            nextsubvol = reactNode[0] / Mreactions;
//        } else {
//            for (i = 0; diffNode[0] >= jcE[i]; i++)
//                ;
//            nextsubvol = (i - 1) / Mspecies;
//        }
        

        // broadcastTime(st->id, tt);

         /* process all msgs occurring before next local event (tt) */
        while (1) {
            tn = DBL_MAX;
            /* process earliest msg first */
            for (int dir = 0; dir < NBS; dir++) {
                if (cb_peek (st->nbs[dir], &msgp) && 
                    msgp.time < tt && msgp.time < tn) {
                    next = dir;
                    tn = msgp.time;
                }
            }
            if (!(tn < DBL_MAX)) break;

            // next is now direction which has message with earliest timestamp
            
            cb_pop (st->nbs[next], &msgp);

            to_vol = msgp.to_vol;
            subvol = lookupLocalVoxel(to_vol);
            spec = msgp.to_spec;
            // is msg older than last update of xx for subvol?
            if (msgp.time < xxt[subvol]) {
                // count out of order msgs
                st->cnt++;
            }
            
            xx[subvol * Mspecies + spec]++;
            xxt[subvol] = MAX(msgp.time, xxt[subvol]);
                    
            //Recalculate the reaction rates using dependency graph G.
            for (i = jcG[spec]; i < jcG[spec + 1]; i++) {
                j = irG[i];
                oldrate = rrate[subvol * Mreactions + j];
                assert(subvol * Mreactions + j < Mreactions * Ncells);
                rrate[subvol * Mreactions + j] = 
                    (*rfun[j])(&xx[subvol * Mspecies],tt,vol[to_vol],
                               &data[to_vol * dsize], sd[to_vol]);
                
                /* Update reaction waiting time in incoming subvolume 
                calcTimes(&reactTimes[reactHeap[subvol * Mreactions + j]],
                          &reactInf[subvol * Mreactions + j], tt, 
                          oldrate, rrate[subvol * Mreactions + j], rng);
                 */
                
                update(reactHeap[subvol * Mreactions + j], reactTimes, 
                       reactNode, reactHeap, reactHeapSize);
            }
            
            // dependent diffusion events
            for (cnt = jcE[subvol * Mspecies + spec];
                 cnt < jcE[subvol * Mspecies + spec + 1]; cnt++) {
                oldrate = isdrate2[cnt];
                isdrate2[cnt] += prD[jvec[cnt]];
#ifdef ROUNDOFF
                if (fabs(isdrate2[cnt]) < 1.0E-12) {
                    isdrate2[cnt] = 0.0;
                }
#endif
                
                calcTimes(&diffTimes[diffHeap[cnt]], &diffInf[cnt], 
                          tt, oldrate, isdrate2[cnt], rng);
                update(diffHeap[cnt], diffTimes, diffNode, diffHeap, diffHeapSize);
            }
            /* update tt, since message may have affected */
            tt = MIN(reactTimes[0],diffTimes[0]);
        }


        
#ifdef BLOCKDIFF
        /* time of next incoming diffusion event (in the future) */
        tn = readNeighbourDiffTimes(st->id, &nextdiffvol);
//        nextdiffvol = lookupLocalVoxel(nextdiffvol);
//&& nextdiffvol == nextsubvol
        while(tt >= tn  && !isinf(tt) && tt < tspan[tlen-1]) {
            idle_wait(200);
            waitcnt++;
            if (cb_peek (st->nbs[0], &msgp) || cb_peek (st->nbs[1], &msgp))
            {
                goto latereceive;
            }
            tn = readNeighbourDiffTimes(st->id, &nextdiffvol);
//            nextdiffvol = lookupLocalVoxel(nextdiffvol);
        }
#else
        tn = readNeighbourTimes(st->id);
        while (tt > tn && !isinf(tt) && tt < tspan[tlen-1]) {
            idle_wait(200);
            tn = readNeighbourTimes(st->id);
        }
#endif
        

        if (isinf(tt) || tt >= tspan[tlen-1]) {
            __sync_fetch_and_add(&block, 1);
            printf("%d: waiting.\n", st->id);
            while (1) {
                
                if (block >= nthreads) {
                    break;
                }

                for (int dir = 0; dir < NBS; dir++) {
                    if (cb_peek (st->nbs[dir],  &msgp))
                    {
                        printf("%d: going back\n", st->id);
                        __sync_fetch_and_add(&block, -1);
                        goto latereceive;
                    }
                }
            }
        }


        if (tt >= tspan[it] || isinf(tt)) {
            if (report_level >= 2){
                printf("%d: %d\n", st->id, it);
            }

            for (; it < tlen && (tt >= tspan[it] || isinf(tt)); it++)
                memcpy(&U[Ndofs * it], &xx[0], Ndofs * sizeof(int));
            if (it >= tlen) {

                for (i = 0; i < tlen; i++) 
                    for (j = 0; j < Ncells; j++) 
                        for (k = 0; k < Mspecies; k++)
                            globU[(vxs[j]*Mspecies+k)+i*globNdofs] = U[j*Mspecies+k + Ndofs*i];

                for (int dir = 0; dir < NBS; dir++) {
                    while (cb_pop (st->nbs[dir], &msgp)) {
                        printf("%d: WARNING\n", st->id);
                    }
                }
                
                break;
            }
        }

        if (reactTimes[0] <= diffTimes[0]) {
            /* Reaction event. */
            event = 0;
            
//                broadcastTime(st->id, tt);
            
            /* a) Extract subvolume and reaction id */
            subvol = reactNode[0] / Mreactions;
            re = reactNode[0] % Mreactions;
            
            /* Write into the event history list */
            rb = (rlb_t*)malloc(sizeof(rlb_t));
            rb->tt=tt;
            rb->type=event;
            rb->id=reactNode[0];  // diffNode[0]
            
            rblist=g_list_append(rblist,(gpointer)rb);
            printf("added event %d at time %f. xx was %d and rate was %f.\n",rb->type,rb->tt, xx[subvol * Mspecies + re],rrate[subvol * Mreactions + re]);
            
                
            /* b) Update the state of the subvolume subvol and
             * diffusion events. */
            /* Loop over species that are affected by the reaction */
            for (i = jcN[re]; i < jcN[re + 1]; i++) { // foreach species?
                assert(irN[i] < Mspecies);
                /* Update the state according to the stochiometry matrix */
                xx[subvol * Mspecies + irN[i]] += prN[i];
                xxt[subvol] = tt;
                if (xx[subvol * Mspecies + irN[i]] < 0)
                    errcode = 1;

                /* Update all dependent diffusion events */
                for (cnt = jcE[subvol * Mspecies + irN[i]];
                     cnt < jcE[subvol * Mspecies + irN[i] + 1]; cnt++) {
                    oldtime = diffTimes[diffHeap[cnt]];
                    assert(cnt < diffHeapSize);
                    assert(isdrate2[cnt] >= 0.0 || isinf(oldtime));
                    /* get old diffusion rate */
                    oldrate = isdrate2[cnt];
                    
                    /* Update the diffusion rate in respect to the
                     * affected species & subvolume */
                    isdrate2[cnt] += prD[jvec[cnt]] * prN[i];

#ifdef ROUNDOFF
                    if (fabs(isdrate2[cnt]) < 1.0E-12) {
                        isdrate2[cnt] = 0.0;
                    }
#endif
                    calcTimes(&diffTimes[diffHeap[cnt]], &diffInf[cnt], 
                              tt, oldrate, isdrate2[cnt], rng);
                    update(diffHeap[cnt], diffTimes, diffNode, 
                           diffHeap, diffHeapSize);
                    react2diff++;
                }
            }

            /* c) Update dependent reaction events (including the
             * actual one) */
            for (i = jcG[Mspecies + re]; i < jcG[Mspecies + re + 1]; i++) {
                j = irG[i];
                oldtime = reactTimes[reactHeap[subvol * Mreactions + j]];
                    
                oldrate = rrate[subvol * Mreactions + j];

                assert(subvol*Mreactions + j < Mreactions * Ncells);

                rrate[subvol * Mreactions + j] = 
                    (*rfun[j])(&xx[subvol * Mspecies], tt, vol[vxs[subvol]], 
                               &data[vxs[subvol] * dsize], sd[vxs[subvol]]);
                    
                calcTimes(&reactTimes[reactHeap[subvol * Mreactions + j]],
                          &reactInf[subvol * Mreactions + j], tt, oldrate,
                          rrate[subvol * Mreactions + j],rng);
                
                update(reactHeap[subvol * Mreactions + j], 
                       reactTimes, reactNode,reactHeap, reactHeapSize);
                react2react++;
            }
            
            total_reactions++; /* counter */
        } else {
            /* Diffusion event. */
            event = 1;
            
            /* Write into the event history list */
            rb = (rlb_t*)malloc(sizeof(rlb_t));
            rb->tt=tt;
            rb->type=event;
            rb->id=diffNode[0];  // diffNode[0]
            
            rblist=g_list_append(rblist,(gpointer)rb);
            printf("added event %d at time %f. xx was %d and rate was %f.\n",rb->type,rb->tt, xx[subvol * Mspecies + re],rrate[subvol * Mreactions + re]);

            /* Determine parameters */

            for (i = 0; diffNode[0] >= jcE[i]; i++)
                ;
            subvol = (i - 1) / Mspecies;
            spec = (i - 1) % Mspecies;
            to_vol = irE[diffNode[0]] / Mspecies;
            to_spec = irE[diffNode[0]] % Mspecies;

            /* Execute the diffusion event (check for negative elements). */
            xx[subvol * Mspecies + spec]--;
            xxt[subvol] = tt;

            if (xx[subvol * Mspecies + spec] < 0) {
                errcode = 1;
            }
            
            int localpos = lookupLocalVoxel(to_vol);
            dom = lookupVoxel(to_vol);


            if (localdom != dom && areNeighbours1d(st->id, dom)) {
                msg = (msg_t) {
                    .time = tt,
                    .to_vol = to_vol,
                    .to_spec = to_spec
                };
                while (!cb_push(receiver1d(st->id, dom), msg));
            }
            else if (localdom == dom) {
                assert(localpos < sd_length[st->id]);
                xx[localpos * Mspecies + to_spec]++;
                xxt[localpos] = tt;
            }
            else {
                // bugger
                assert(0);
            }

            /* Recalculate the reaction rates using dependency graph G. */
            for (i = jcG[spec]; i < jcG[spec + 1]; i++) {
                j = irG[i];
                
                oldrate = rrate[subvol * Mreactions + j];
                assert(subvol*Mreactions + j < Mreactions * Ncells);
                rrate[subvol * Mreactions + j] = 
                    (*rfun[j])(&xx[subvol * Mspecies], tt, vol[vxs[subvol]],                                
                        &data[vxs[subvol] * dsize], sd[vxs[subvol]]);

                
                /* Update reaction waiting time in outgoing subvolume */
                calcTimes(&reactTimes[reactHeap[subvol * Mreactions + j]],
                          &reactInf[subvol * Mreactions + j], tt, oldrate,
                          rrate[subvol * Mreactions + j], rng);
                
                update(reactHeap[subvol * Mreactions + j], reactTimes, 
                       reactNode, reactHeap, reactHeapSize);

                /* only update reaction rates of incoming subvolume if
                 * it is internal to the outgoing subvolume's
                 * subdomain */
                if (dom == localdom) {
                    assert(localpos < sd_length[st->id]);
                    oldrate = rrate[localpos * Mreactions + j];
                    assert(localpos * Mreactions + j < Mreactions * Ncells);
                    rrate[localpos * Mreactions + j] = 
                        (*rfun[j])(&xx[localpos * Mspecies],tt,vol[to_vol], 
                                   &data[to_vol * dsize], sd[to_vol]);

                    /* Update reaction waiting time in incoming subvolume */
                    calcTimes(&reactTimes[reactHeap[localpos * Mreactions + j]],
                              &reactInf[localpos * Mreactions + j], tt, 
                              oldrate, rrate[localpos * Mreactions + j], rng);

                    update(reactHeap[localpos * Mreactions + j], reactTimes, 
                           reactNode, reactHeap, reactHeapSize);
                }
                // FIX: split this one?
                diff2react += 2;
            }

            /* Update all diffusion events of affected species in
             * the outgoing subvolume */
            /* In the AEM context decrease the intensity of all
             * events related to non-diagonal row entries of D */
            for (cnt = jcE[subvol * Mspecies + spec];
                 cnt < jcE[subvol * Mspecies + spec + 1]; cnt++) {
                assert(cnt < diffHeapSize);
                
                oldrate = isdrate2[cnt];
                isdrate2[cnt] -= prD[jvec[cnt]];
#ifdef ROUNDOFF  
                if (fabs(isdrate2[cnt]) < 1.0E-12) {
                    isdrate2[cnt] = 0.0;
                }
#endif

                calcTimes(&diffTimes[diffHeap[cnt]], &diffInf[cnt], tt, 
                          oldrate, isdrate2[cnt], rng);

                update(diffHeap[cnt], diffTimes, diffNode, 
                       diffHeap, diffHeapSize);
                diff2diff++;
            }

            /* only update the incoming subvolume if it is internal to the
             * outgoing subvolume's subdomain */
            if (dom == localdom) {
                /* Update all diffusion events of affected species in
                 * the incoming subvolume */
                assert(localpos < sd_length[st->id]);
                for (cnt = jcE[localpos * Mspecies + to_spec];
                     cnt < jcE[localpos * Mspecies + to_spec + 1]; cnt++) {
                    assert(cnt < diffHeapSize);
                    oldrate = isdrate2[cnt];
                    isdrate2[cnt] += prD[jvec[cnt]];
#ifdef ROUNDOFF
                    if (fabs(isdrate2[cnt]) < 1.0E-12) {
                        isdrate2[cnt] = 0.0;
                    }
#endif

                    calcTimes(&diffTimes[diffHeap[cnt]], &diffInf[cnt], 
                              tt, oldrate, isdrate2[cnt], rng);
                    update(diffHeap[cnt], diffTimes, diffNode, diffHeap, diffHeapSize);
                    diff2diff++;
                }
            }

            total_diffusion++; /* counter */
 }

        printf("reactTimes[0]: %f reactTimes[1]: %f \n",reactTimes[0],reactTimes[1]);
        printf("diffTimes[0]: %f diffTimes[1]: %f diffTimes[2]: %f \n",diffTimes[0],diffTimes[1],diffTimes[2]);
        printf("-------------------------------------------------------\n");

        lf = 0; rf = 0;
        if (st->id == 0) lf = 1;
        if (st->id == nthreads - 1) rf = 1;
        
        for (i = 0; i < diffHeapSize; i++) {
            if (lf && rf) break;
            to_vol = irE[diffNode[i]] / Mspecies;
            dom = lookupVoxel(to_vol);
            if (dom == st->id + 1 && !rf) {
                rf = 1;
                if (st->nextdiff[1] != diffTimes[i]) {
                    st->nextdiff[1] = diffTimes[i];
                }
            } else if (dom == st->id - 1 && !lf) {
                lf = 1;
                if (st->nextdiff[0] != diffTimes[i]) {
                    st->nextdiff[0] = diffTimes[i];
                }
            }
        }
        if (lf != 1) {
            st->nextdiff[0] = INFINITY;
        } else if (rf != 1) {
            st->nextdiff[0] = INFINITY;
        }
        


        if (errcode) {
            memcpy(&U[Ndofs * it], &xx[0], Ndofs * sizeof(int));
            memset(&U[Ndofs * (it + 1)], -1, Ndofs * (tlen - it - 1) * sizeof(int));
            printf("errcode: %d\n", errcode);
            break;
        }
    }

    if (report_level == 2) {
        printf("%d waitcnt: %d\n", st->id, waitcnt);
        printf("Updates: %lu Wake-ups: %lu\n", updates, wakeups);
//        printf("React2React: %lu React2diff: %lu\n", react2react, react2diff);
//        printf("Diff2diff: %lu Diff2react: %lu\n", diff2diff, diff2react);
    }
    
    printf("****************************************************\n");
    printf("****************************************************\n");
    printf("****************************************************\n");
    
    testHistoryOut(1.6);
    //testHistoryOut(0.02);
  
    free(reactHeap);
    free(reactNode);
    free(reactTimes);
    free(reactInf);
    
    free(diffHeap);
    free(diffNode);
    free(diffTimes);
    free(diffInf);
    free(jcE);
    free(irE);
        
//	free(sdrate);
//	free(rrate);
//	free(xx);
}



/*----------------------------------------------------------------------- */
/* CalcTimes: calculate update of waiting times, inclusive sleeping times */
/*----------------------------------------------------------------------- */

void calcTimes(double* time, double* infTime, double tt, double old_rate,
               double new_rate, seeds seed) {

    double oldtime = time[0];
    if (isinf(oldtime)) {

        if (infTime[0] == 0) { // Waking up first time
            if (new_rate > 0.0) {
                time[0] = -log(1.0 - erand48(seed)) / new_rate + tt;
                
                
            }
            
        } else { // Waking up the 2nd..nth time
            if (new_rate > 0.0) {
                time[0] = tt + (infTime[0] / new_rate);
            }
        }
        wakeups++;
    } else {
        if (new_rate >= DBL_MIN) {
            if (time[0] == tt) {// Regular update of current event 
                time[0] = -log(1.0 - erand48(seed)) / new_rate + tt;
                
            } else { // Regular update of dependent events (rescaling)
                time[0] = ((old_rate / new_rate) * (time[0] - tt)) + tt;
            }
            updates++;

        } else { // Next event time set to infinity
            /* BUG: this only holds for the currently non-active event*/
            infTime[0] = (oldtime - tt) * old_rate;
            time[0] = INFINITY;
            //printf("going to infinity with inftime=%f (oldtime=%f, old_rate=%f.\n",infTime[0],oldtime,old_rate);
        }
    }
}
