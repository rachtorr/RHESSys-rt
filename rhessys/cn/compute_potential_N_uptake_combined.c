/*--------------------------------------------------------------*/
/* 								*/
/*		compute_potential_N_uptake_combined		*/
/*								*/
/*								*/
/*	NAME							*/
/*	compute_potential_N_uptake_combined -			*/
/*		integrates dickenson allocation within		*/
/*		original Waring potential N uptake file		*/
/*								*/
/*		computes potential N uptake from soil		*/
/*		for this strata without mineralization		*/
/*		limitation					*/
/*								*/
/*	SYNOPSIS						*/
/*	int compute_potential_N_uptake_combined(		*/
/*                          struct epconst_struct,              */
/*			    struct epvar_struct *epv,		*/
/*                          struct cstate_struct *,             */
/*                          struct nstate_struct *,             */
/*                          struct cdayflux_struct *)           */
/*								*/
/*								*/
/*	returns int:						*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*								*/
/*              modified from Peter Thornton (1998)             */
/*                      dynamic - 1d-bgc ver4.0                 */
/*--------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"
double compute_potential_N_uptake_combined(
								  struct	epconst_struct epc,
								  struct	epvar_struct *epv,
								  struct cstate_struct *cs,
								  struct nstate_struct *ns,
								  struct cdayflux_struct *cdf)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	double day_gpp;     /* daily gross production */
	double day_mresp;   /* daily total maintenance respiration */
	double fbroot, fcroot, fstem, froot, fleaf, fwood;          /* RATIO   new fine root C : new leaf C     */
	double f2;          /* RATIO   fraction to leaf and fraction to root*/
	double f4;          /* RATIO   new live wood C : new wood C     */
	double f3;          /* RATIO   new leaf C : new wood C     */
	double g1;          /* RATIO   C respired for growth : C grown  */
	double cnl;         /* RATIO   leaf C:N      */
	double cnfr;        /* RATIO   fine root C:N */
	double cnlw;        /* RATIO   live wood C:N */
	double cndw;        /* RATIO   dead wood C:N */
	double cnmax;       /* RATIO   max of root and leaf C:N      */
	double mean_cn, transfer, ratio;
	double plant_calloc, plant_ndemand;
	double k2, c; /* working variables */
	double dickenson_k; /* working variable for LAI exponential decay constant */
	/*---------------------------------------------------------------
	Assess the carbon availability on the basis of this day's
	gross production and maintenance respiration costs
	----------------------------------------------------------------*/
	cs->availc = cdf->psn_to_cpool-cdf->total_mr;
	/* no allocation when the daily C balance is negative */
	if (cs->availc < 0.0) cs->availc = 0.0;
	/* test for cpool deficit */
	if (cs->cpool < 0.0){
	/*--------------------------------------------------------------
	running a deficit in cpool, so the first priority
	is to let today's available C accumulate in cpool.  The actual
	accumulation in the cpool is resolved in day_carbon_state().
		--------------------------------------------------------------*/
		/*------------------------------------------------
		cpool deficit is less than the available
		carbon for the day, so aleviate cpool deficit
		and use the rest of the available carbon for
		new growth and storage.
			-----------------------------------------------*/
			transfer = min(cs->availc, -cs->cpool);
			cs->availc -= transfer;
			cs->cpool += transfer;
	} /* end if negative cpool */
	/* assign local values for the allocation control parameters */
	f2 = epc.alloc_crootc_stemc;
	f3 = epc.alloc_stemc_leafc;
	f4 = epc.alloc_livewoodc_woodc;
	g1 = epc.gr_perc;
	cnl = epc.leaf_cn;
	cnfr = epc.froot_cn;
	cnlw = epc.livewood_cn;
	cndw = epc.deadwood_cn;
	/*--------------------------------------------------------------- */
	/*	given the available C, use Waring allometric relationships to */
	/*	estimate N requirements -					*/ 
	/*	constants a and b are taken from Landsberg and Waring, 1997 */
	/*--------------------------------------------------------------*/
	
	if (((cdf->potential_psn_to_cpool) > ZERO) && (cdf->psn_to_cpool > ZERO)) {
	c = max(cdf->potential_psn_to_cpool, cdf->psn_to_cpool);
	plant_calloc = cs->availc;
	if ((plant_calloc > ZERO) && (c > 0)){
		
		/*--------------------------------------------------------------- 
		combined allocation 
		waring_pa and waring_pb control the exponential decay constant (k)
		in the dickenson allocation (i.e. old dickenson_pa) 
		original dickenson allocation was fleaf=exp(-1*k*lai) use  
		froot = (1-leaf) for simplicity  
		/* --------------------------------------------------------------- */

			/* --------------------------------------------------------------- */
			/* uses approach published in Reyes et al., (2017) Assessing the Impact of...JAMES  */
			/* adapted for trees here (unpublished) */
			/* --------------------------------------------------------------- */

		if (epc.veg_type == TREE) {


		fleaf = exp(-1.0*epc.waring_pa * epv->proj_lai);
		fbroot = fleaf * epc.alloc_frootc_leafc *  (1+epc.waring_pb )/ (1.0 + epc.waring_pb * (cdf->psn_to_cpool)/c);
		ratio = fbroot/fleaf;

		if ((epc.veg_type == TREE)) {
		if (fbroot+fleaf > 0.95) {
			fleaf = 0.95/(1+ratio);
			fbroot = fleaf*ratio;
			}
		
		fcroot = fbroot/(1+epc.alloc_frootc_crootc);
		froot = fbroot-fcroot;

		fstem = 1.0-(froot+fcroot+fleaf);
		fwood = fstem+fcroot;
		}
		else {
			fleaf = 1.0-fbroot;
			froot = fbroot;
		}
		}

		else {
		dickenson_k = (epc.waring_pa / (1.0 + epc.waring_pb * (cdf->psn_to_cpool) / c));
		froot = (1-exp(-1.0*dickenson_k * epv->proj_lai));
		fleaf = 1.0-froot;
		fcroot=0.0;
		fwood=0.0;
		}

	}
	else {
		froot = 0.0;
		fleaf = 0.0;
		fstem=0.0;
		fcroot=0.0;
		fwood=0.0;
		f3 = 0.0;
	}




	if (epc.veg_type == TREE){
		if ((fleaf + froot + fwood) > ZERO) 
			mean_cn = 1.0 / (fleaf / cnl + froot / cnfr + f4 * fwood / cnlw + fwood * (1.0-f4) / cndw);
		else mean_cn = 1.0;
	}
        else{
	if ((fleaf + froot) > ZERO) 	
           mean_cn = 1.0 / (fleaf / cnl + froot / cnfr);
	else
		mean_cn=1.0;
        }

	if (mean_cn > ZERO)
		plant_ndemand = cs->availc / (1.0+epc.gr_perc) / mean_cn;
	else
		plant_ndemand = 0.0;

	}
	else {
		plant_ndemand = 0.0;
		fleaf = 0.0;
		froot = 0.0;
		fwood = 0.0;
		fcroot = 0.0;
	}

	cdf->fleaf = fleaf;
	cdf->fwood = fwood;
	cdf->froot = froot;
	cdf->fcroot = fcroot;

/*
	if (c > ZERO)
		printf("\n lai %lf  fleaf %lf froot %lf fcroot %lf fstem %lf stress %lf", epv->proj_lai, fleaf, froot, fcroot, fstem, (cdf->psn_to_cpool)/c);
	else
		printf("\n lai %lf  fleaf %lf froot %lf fcroot %lf fstem %lf nopsn", epv->proj_lai, fleaf, froot, fcroot, fstem );
*/

	return(plant_ndemand);
} /* 	end compute_potential_N_uptake_Waring */
