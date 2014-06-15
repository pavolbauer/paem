/* aemore.c - Core AEM solver. */

/* P. Bauer and S. Engblom 2012-05-10 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "propensities.h"
#include "aem.h"
#include "../nsm/binheap.h"
#include "report.h"

/* Global variables for debug statistics */
int updates = 0;
int wakeups = 0;
int react2diff = 0;
int react2react = 0;
int diff2react = 0;
int diff2diff = 0;

void aem_core(const size_t *irD, const size_t *jcD, const double *prD,
		const int *u0, const size_t *irN, const size_t *jcN, const int *prN,
		const size_t *irG, const size_t *jcG, const double *tspan,
		const size_t tlen, int *U, const double *vol, const double *data,
		const int *sd, const size_t Ncells, const size_t Mspecies,
		const size_t Mreactions, const size_t dsize, int report_level)

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

	double *rrate, old;
	double *sdrate, *Ddiag, *isdrate2;
	double *reactTimes, *diffTimes, *reactInf, *diffInf;

	double tt = tspan[0];
	double oldtime, oldrate;

	int *reactNode, *reactHeap, *diffNode, *diffHeap, *xx;
	int *jcE, *irE, *xx_diag, *jvec;

	int i, j, re, cnt, nSub, i_temp;
	int diffHeapSize, reactHeapSize, event;
	int subvol, spec, to_vol, to_spec;

	short errcode = 0;
	long int total_diffusion = 0, total_reactions = 0;

	const size_t Ndofs = Ncells * Mspecies;
	seeds *reactSeeds, *diffSeeds;
	size_t it = 0;

	PropensityFun *rfun;
	rfun = ALLOC_propensities();

	ReportFun report;
	if (report_level)
		report = &reportFun1;
	else
		report = NULL;

	//printf("subdomain of voxel 1223 is %d.\n",lookupVoxel(1223,sd));

	/* TODO
	 * 1) make heaps private -> take psiminf code
	 *
	 * 2) divide Ncells:
	 *
	 * for (i = 0; i < Ncells; i++)
	 *   Ncells_thread(sd[i]).push(i);
	 *
	 * 3) write into global xx using mapping function x->y
	 *
	 * 4) divide diffusion matrix
	 * for(i=0;i<Ndfofs;i+=Mspecies)
	 * 	for (j = jcD[i]; j < jcD[i + 1]; j++) {
	 * 	  thread[sd[i%2]].jcD[j]=jcD[j];
	 *		thread[sd[i%2]].prD[j]=prD[j];
	 * }
	 *
	 * 5) restructure diffusion matrix privately for each subdomain
	 */


	/* Set xx to the initial state. --> global */
	xx = (int *) malloc(Ndofs * sizeof(int));
	memcpy(xx, u0, Ndofs * sizeof(int));

	/* Create reaction rate matrix (Mreactions X Ncells) and total rate
	 vector. In rrate we store all propensities for chemical rections,
	 and in srrate the sum of propensities in every subvolume. */
	rrate = (double *) malloc(Mreactions * Ncells * sizeof(double));

	/* Calculate the propensity for every reaction and every
	 subvolume. Store the sum of the reaction intensities in each
	 subvolume in srrate. */
	for (i = 0; i < Ncells; i++) {
		for (j = 0; j < Mreactions; j++) {
			rrate[i * Mreactions + j] = (*rfun[j])(&xx[i * Mspecies], tt, vol[i],
					&data[i * dsize], sd[i]);
		}
	}

	/* Binary heap storing all reaction events */
	reactHeapSize = Ncells * Mreactions;

	/* Create Reaction-heap. */
	reactTimes = (double *) malloc(reactHeapSize * sizeof(double));
	reactNode = (int *) malloc(reactHeapSize * sizeof(int));
	reactHeap = (int *) malloc(reactHeapSize * sizeof(int));
	reactSeeds = (seeds *) malloc(reactHeapSize * sizeof(seeds));

	/* Create Array to store pre-infinity times */
	reactInf = (double *) malloc(reactHeapSize * sizeof(double));
	for (i = 0; i < reactHeapSize; i++)
		reactInf[i] = 0;

	/* Initialize erand48 seeds */
	for (i = 0; i < reactHeapSize; i++)
		for (j = 0; j < 3; j++)
			reactSeeds[i][j] = rand() % 65535;

	/* Calculate times to next reaction event */
	for (i = 0; i < reactHeapSize; i++) {
		reactTimes[i] = -log(1.0 - erand48(reactSeeds[i])) / rrate[i] + tspan[0];

		if (reactTimes[i] <= 0.0)
			reactTimes[i] = INFINITY;

		reactHeap[i] = reactNode[i] = i;
	}

	/* Initialize reaction heap */
	initialize_heap(reactTimes, reactNode, reactHeap, reactHeapSize);

	/* The diagonal value of the D-matrix is used frequently. For
	 efficiency, we store the negative of D's diagonal in Ddiag. */
	nSub = 0;
	Ddiag = (double *) malloc(Ndofs * sizeof(double));
	for (i = 0; i < Ndofs; i++) {
		Ddiag[i] = 0.0;
		for (j = jcD[i]; j < jcD[i + 1]; j++)
			if (irD[j] == i)
				Ddiag[i] = -prD[j];
		if (Ddiag[i] > 0)
			nSub++;
	}

	/* Store all entries in diagonal row */
	xx_diag = (int *) malloc(Ndofs * sizeof(int));

	for (i = 0; i < Ncells; i++) {
		for (j = 0; j < Mspecies; j++) {
			xx_diag[i * Mspecies + j] = xx[i * Mspecies + j];
		}
	}

	/* diffHeapSize for unidirectional diffusion problems -> subtract non-zero'd Nrows only */
	diffHeapSize = (jcD[Ndofs] - nSub);

	/* allocate space for diffusion-event heap */
	if (diffHeapSize < 0)
		diffHeapSize = 0;
	diffTimes = (double *) malloc(diffHeapSize * sizeof(double));
	diffNode = (int *) malloc(diffHeapSize * sizeof(int));
	diffHeap = (int *) malloc(diffHeapSize * sizeof(int));
	diffSeeds = (seeds *) malloc(diffHeapSize * sizeof(seeds));
	isdrate2 = (double *) malloc(diffHeapSize * sizeof(double));

	/* Create Array to store pre-infinity times */
	diffInf = (double *) malloc(diffHeapSize * sizeof(double));
	for (i = 0; i < diffHeapSize; i++)
		diffInf[i] = 0;

	/* Initialize erand48 seeds */
	for (i = 0; i < diffHeapSize; i++)
		for (j = 0; j < 3; j++)
			diffSeeds[i][j] = rand() % 65535;

	/* Creating new sparse matrix for diffusion events */
	jcE = (int*) malloc((Ndofs + 1) * sizeof(int));
	irE = (int*) malloc(diffHeapSize * sizeof(int));
	jvec = (int*) malloc(diffHeapSize * sizeof(int));

	/* Queue all non-diagonal entries of D as diffusion events */
	for (cnt = 0, i = 0; i < Ndofs; i++) {
		jcE[i] = cnt;
		for (j = jcD[i]; j < jcD[i + 1]; j++)
			if (irD[j] != i) {
				diffNode[cnt] = diffHeap[cnt] = cnt;
				isdrate2[cnt] = xx_diag[i] * prD[j];
				diffTimes[cnt] = -log(1.0 - erand48(diffSeeds[cnt])) / isdrate2[cnt]
						+ tspan[0];

				if (diffTimes[cnt] <= 0.0)
					diffTimes[cnt] = INFINITY;
				irE[cnt] = irD[j];
				jvec[cnt] = j;

				cnt++;
			}
	}
	jcE[i] = cnt;

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

	/* Main loop. */
	for (;;) {

		/* Calculate tt with diff-time and react-time */
		if (reactTimes[0] <= diffTimes[0])
			tt = reactTimes[0];
		else
			tt = diffTimes[0];

		if (tt >= tspan[it] || isinf(tt)) {
			for (; it < tlen && (tt >= tspan[it] || isinf(tt)); it++)
				memcpy(&U[Ndofs * it], &xx[0], Ndofs * sizeof(int));
			if (it >= tlen)
				break;
		}

		if (reactTimes[0] <= diffTimes[0]) {
			/* Reaction event. */
			event = 0;

			/* a) Extract subvolume and reaction id */
			subvol = reactNode[0] / Mreactions;
			re = reactNode[0] % Mreactions;

			/* b) Update the state of the subvolume subvol and diffusion events. */
			/* Loop over species that are affected by the reaction */
			for (i = jcN[re]; i < jcN[re + 1]; i++) {

				/* Update the state according to the stochiometry matrix*/
				xx[subvol * Mspecies + irN[i]] += prN[i];
				if (xx[subvol * Mspecies + irN[i]] < 0)
					errcode = 1;

				/* Update all dependent diffusion events */
				for (cnt = jcE[subvol * Mspecies + irN[i]];
						cnt < jcE[subvol * Mspecies + irN[i] + 1]; cnt++) {
					oldtime = diffTimes[diffHeap[cnt]];

					/* get old diffusion rate */
					oldrate = isdrate2[cnt];

					/* Update the diffusion rate in respect to the affected species & subvolume */
					isdrate2[cnt] += prD[jvec[cnt]] * prN[i];

					calcTimes(&diffTimes[diffHeap[cnt]], &diffInf[cnt], tt, oldrate,
							isdrate2[cnt], diffSeeds[cnt]);

					update(diffHeap[cnt], diffTimes, diffNode, diffHeap, diffHeapSize);
					react2diff++;
				}
			}

			/* c) Update dependent reaction events (including the actual one) */
			for (i = jcG[Mspecies + re]; i < jcG[Mspecies + re + 1]; i++) {
				j = irG[i];
				oldtime = reactTimes[reactHeap[subvol * Mreactions + j]];

				oldrate = rrate[subvol * Mreactions + j];
				rrate[subvol * Mreactions + j] = (*rfun[j])(&xx[subvol * Mspecies], tt,
						vol[subvol], &data[subvol * dsize], sd[subvol]);

				calcTimes(&reactTimes[reactHeap[subvol * Mreactions + j]],
						&reactInf[subvol * Mreactions + j], tt, oldrate,
						rrate[subvol * Mreactions + j],
						reactSeeds[subvol * Mreactions + j]);

				update(reactHeap[subvol * Mreactions + j], reactTimes, reactNode,
						reactHeap, reactHeapSize);
				react2react++;
			}

			total_reactions++; /* counter */
		} else {
			/* Diffusion event. */
			event = 1;

			/* Determine parameters */
			for (i = 0; diffNode[0] >= jcE[i]; i++)
				;

			i_temp = i;
			subvol = (i - 1) / Mspecies;
			spec = (i - 1) % Mspecies;
			to_vol = irE[diffNode[0]] / Mspecies;
			to_spec = irE[diffNode[0]] % Mspecies;

			/* Execute the diffusion event (check for negative elements). */
			xx[subvol * Mspecies + spec]--;
			if (xx[subvol * Mspecies + spec] < 0)
				errcode = 1;
			xx[to_vol * Mspecies + to_spec]++;

			/* Recalculate the reaction rates using dependency graph G. */
			for (i = jcG[spec]; i < jcG[spec + 1]; i++) {
				old = rrate[subvol * Mreactions + irG[i]];
				j = irG[i];

				oldrate = rrate[subvol * Mreactions + j];
				rrate[subvol * Mreactions + j] = (*rfun[j])(&xx[subvol * Mspecies], tt,
						vol[subvol], &data[subvol * dsize], sd[subvol]);

				/* Update reaction waiting time in outgoing subvolume */

				calcTimes(&reactTimes[reactHeap[subvol * Mreactions + j]],
						&reactInf[subvol * Mreactions + j], tt, oldrate,
						rrate[subvol * Mreactions + j],
						reactSeeds[subvol * Mreactions + j]);

				update(reactHeap[subvol * Mreactions + j], reactTimes, reactNode,
						reactHeap, reactHeapSize);

				oldrate = rrate[to_vol * Mreactions + j];
				rrate[to_vol * Mreactions + j] = (*rfun[j])(&xx[to_vol * Mspecies], tt,
						vol[to_vol], &data[to_vol * dsize], sd[to_vol]);

				/* Update reaction waiting time in incoming subvolume */

				calcTimes(&reactTimes[reactHeap[to_vol * Mreactions + j]],
						&reactInf[to_vol * Mreactions + j], tt, oldrate,
						rrate[to_vol * Mreactions + j],
						reactSeeds[to_vol * Mreactions + j]);

				update(reactHeap[to_vol * Mreactions + j], reactTimes, reactNode,
						reactHeap, reactHeapSize);
				diff2react += 2;
			}

			/* Update all diffusion events of affected species in the outgoing subvolume */
			/* In the AEM context decrease the intensity of all events related to non-diagonal row entries of D */
			for (cnt = jcE[subvol * Mspecies + spec];
					cnt < jcE[subvol * Mspecies + spec + 1]; cnt++) {
				oldrate = isdrate2[cnt];
				isdrate2[cnt] -= prD[jvec[cnt]];

				calcTimes(&diffTimes[diffHeap[cnt]], &diffInf[cnt], tt, oldrate,
						isdrate2[cnt], diffSeeds[cnt]);

				update(diffHeap[cnt], diffTimes, diffNode, diffHeap, diffHeapSize);
				diff2diff++;
			}

			/* Update all diffusion events of affected species in the incoming subvolume */
			for (cnt = jcE[to_vol * Mspecies + spec];
					cnt < jcE[to_vol * Mspecies + spec + 1]; cnt++) {
				oldrate = isdrate2[cnt];

				isdrate2[cnt] += prD[jvec[cnt]];

				calcTimes(&diffTimes[diffHeap[cnt]], &diffInf[cnt], tt, oldrate,
						isdrate2[cnt], diffSeeds[cnt]);

				update(diffHeap[cnt], diffTimes, diffNode, diffHeap, diffHeapSize);
				diff2diff++;
			}

			total_diffusion++; /* counter */
		}

		if (errcode) {
			memcpy(&U[Ndofs * it], &xx[0], Ndofs * sizeof(int));
			memset(&U[Ndofs * (it + 1)], -1, Ndofs * (tlen - it - 1) * sizeof(int));
			break;
			break;
		}
	}

	if (report_level == 2) {
		printf("Updates: %d Wake-ups: %d \n", updates, wakeups);
		printf("React2React: %d React2diff: %d \n", react2react, react2diff);
		printf("Diff2diff: %d Diff2react: %d \n", diff2diff, diff2react);
	}

	free(reactHeap);
	free(reactNode);
	free(reactTimes);
	free(reactSeeds);
	free(reactInf);

	free(diffHeap);
	free(diffNode);
	free(diffTimes);
	free(diffSeeds);
	free(diffInf);
	free(jcE);
	free(irE);

	FREE_propensities(rfun);
	free(Ddiag);
	free(sdrate);
	free(rrate);
	free(xx);
}

/*----------------------------------------------------------------------- */
/* CalcTimes: calculate update of waiting times, inclusive sleeping times */
/*----------------------------------------------------------------------- */

void calcTimes(double* time, double* infTime, double tt, double old_rate,
		double new_rate, seeds seed) {
	double oldtime = time[0];

	if (isinf(oldtime)) {

		if (infTime[0] == 0) { // Waking up first time
			time[0] = -log(1.0 - erand48(seed)) / new_rate + tt;
		} else { // Waking up the 2nd..nth time
			if (new_rate > 0.0) {
				time[0] = tt + (infTime[0] / new_rate);
			}
		}
		wakeups++;
	} else {
		if (new_rate >= DBL_MIN) {
			if (time[0] == tt) // Regular update of current event
				time[0] = -log(1.0 - erand48(seed)) / new_rate + tt;
			else
				// Regular update of dependent events (rescaling)
				time[0] = ((old_rate / new_rate) * (time[0] - tt)) + tt;
			updates++;
		} else { // Next event time set to infinity

			infTime[0] = (oldtime - tt) * old_rate;
			time[0] = INFINITY;
		}
	}
}

/*-----------------------------------------------------------------------

int
lookupVoxel(int voxel,const int *sd)
{
    return sd[voxel];
}

----------------------------------------------------------------------*/
