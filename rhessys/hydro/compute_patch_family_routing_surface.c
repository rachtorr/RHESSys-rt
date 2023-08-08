/*--------------------------------------------------------------*/
/*					                                            */
/*		compute_patch_family_routing		                    */
/*					                                            */
/*	compute_patch_family_routing.c - routes within patch family	*/
/*					                                            */
/*	NAME						                                */
/*	compute_patch_family_routing.c - routes within patch family	*/
/*					                                            */
/*	SYNOPSIS			                                        */
/*	void compute_patch_family_routing( 	                        */
/*						    struct zone_object *zone)           */
/*										                        */
/*	OPTIONS								                    	*/
/*										                        */
/*										                        */
/*	DESCRIPTION								                    */
/*  For all patch families in a zone, routes surface water between      */
/*	patches in each patch family. Routing is based on detention  */
/*	store available water and ability to receive water      */
/*  modified by surf_g coefficients.                     */
/*                                                              */
/*	PROGRAMMER NOTES							                */
/*										                        */
/*	July, 2023 RT   						            */
/*											                    */
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h> // i think not needed?
#include <math.h>
#include "../include/rhessys.h"

void compute_patch_family_routing_surface(struct zone_object *zone,
                                  struct command_line_object *command_line,
		                          struct date current_date,
                                  int n_timesteps)

{

double compute_infiltration(int, double, double, double, double, double,
			double, double, double, double, double);

    /*--------------------------------------------------------------*/
    /*	Local variable definition.			                        */
    /*--------------------------------------------------------------*/
    int verbose_flag; 
    int grow_flag; 
    int pf;
    int i;
    int p_ct;            // number of patches in patch family
    int p_ct_incl_surf;       // number of patches in patch family without skipped patches
    double wet_mean_surf;     // mean wetness for surface above detention store size, meters water
    double area_sum_loss;     // sum of areas in patch family
    double area_sum_g;       // sum of gaining patches area 
    double dG_pot;       // sum of (potential) surface gains over entire family 
    double dL_pot;      // sum of surface losses over entire family 
    double surf_trans; //intermediate vars for surface transfer
    double time_int;  //
    double infiltration; 


    /* initializations */
    time_int = 1.0 / n_timesteps;
    grow_flag = command_line[0].grow_flag;
    verbose_flag = command_line[0].verbose_flag;

    // for testing
    // if (current_date.year == 2011 && current_date.month == 3 && current_date.day == 22) {
    //     command_line[0].verbose_flag = -6;
    // } else {
    //     command_line[0].verbose_flag = 0;
    // }

    /*--------------------------------------------------------------*/
    /*	Loop through patch families in the zone   	                */
    /*--------------------------------------------------------------*/
    for (pf = 0; pf < zone[0].num_patch_families; pf++)
    {
        if (verbose_flag == -6) {printf("\n--- Patch Family %d ---\n", zone[0].patch_families[pf][0].family_ID);}

        /*--------------------------------------------------------------*/
        /*	Patch family definitions & initializations                  */
        /*--------------------------------------------------------------*/
        p_ct = zone[0].patch_families[pf][0].num_patches_in_fam; // for simplicity

        /* Definitions */
        int incl_surf[p_ct];          // 0 = not included, 1 = include, 2 = gaining patch

        double dL[p_ct];         // loses of water from surface from patch, vol water
        double dG[p_ct];         // gains of water from surface from patch, vol water
        double wet_surf[p_ct];  // wetnes on surface above detention store size 

        /* Initializations */
        p_ct_incl_surf = 0;
        wet_mean_surf = 0;
        area_sum_loss = 0;
        area_sum_g = 0;
        dG_pot = 0;
        dL_pot = 0;
        surf_trans = 0; 


        /*--------------------------------------------------------------*/
        /*	Loop 1 - Get mean amount of water above det store capacity            */
        /*--------------------------------------------------------------*/
        if (verbose_flag == -6)
        {
            printf("|| Pre-transfer ||\n|   ID   | Include |  Area  |  Detention store | surface water above limit | pct |\n");
        }
        // start loop for number of patches in family, for the entire thing 
        for (i = 0; i < zone[0].patch_families[pf][0].num_patches_in_fam; i++)
        {
            /* Initializations */
            zone[0].patch_families[pf][0].patches[i][0].surface_transfer = 0;

            // check if detention store is >0
            if (zone[0].patch_families[pf][0].patches[i][0].detention_store > ZERO)
            {
                incl_surf[i] = 1;
            }
            else
            {
                incl_surf[i] = 0;
            }

            
            // if both sh coefficients are not 0, include patch; we :use sh_l, only g , don't need to check those 
            if (zone[0].patch_families[pf][0].num_patches_in_fam > 1)
            {
                /* include patches with detention store greater than 0 (incl_surf > 0) */
		/* RACHEL - the condition below will always be true based on line 120-126 */
                if (incl_surf[i] > ZERO )
                {
                    // incrament mean wetness based on storage (amount > detention_store_size) * area
                    if (zone[0].patch_families[pf][0].patches[i][0].detention_store > zone[0].patch_families[pf][0].patches[i][0].landuse_defaults[0][0].detention_store_size)
                    {
                        wet_surf[i] = (zone[0].patch_families[pf][0].patches[i][0].detention_store - 
                                    zone[0].patch_families[pf][0].patches[i][0].landuse_defaults[0][0].detention_store_size);
                    }
                    else {
                        wet_surf[i] = 0.0; 
                    }
                    
                    wet_mean_surf += wet_surf[i] * zone[0].patch_families[pf][0].patches[i][0].area;
                    p_ct_incl_surf += 1;
                }
            }

            // area sum for all patches able to lose 
            area_sum_loss += zone[0].patch_families[pf][0].patches[i][0].area;
            
            //if surf_g is 1, count in area for able to gain 
            if(zone[0].patch_families[pf][0].patches[i][0].landuse_defaults[0][0].surf_g==1) 
            {
                area_sum_g += zone[0].patch_families[pf][0].patches[i][0].area;
            }   

            if (verbose_flag == -6)
            {
                      //|   ID   | Include |  Area  |  Detention store | surface water above limit | pct |
                printf("| %4d |  %12d | %2.12f | %2.12f | %2.12f | %8d | \n",
                       zone[0].patch_families[pf][0].patches[i][0].ID,
                       incl_surf[i],
                       zone[0].patch_families[pf][0].patches[i][0].area,
                       zone[0].patch_families[pf][0].patches[i][0].detention_store,
                       wet_surf[i],
                       p_ct_incl_surf);
            }

        } // end loop 1

		/* RACHEL - if patch family only has one patch you still do below - but this needs to be flaged */
		/* because if there is only one patch you don't need to do any transfer */
		/* otherwise you will go through without seting wet_surf etc for a single patch patch family */

        if(zone[0].patch_families[pf][0].num_patches_in_fam>1){
    
            // Get mean wetness - vol water/(total patch family) area - units are meters depth
            if (area_sum_loss > ZERO) {
                wet_mean_surf /= area_sum_loss;
            }
            else {
                wet_mean_surf = 0.0;
            }

            if (verbose_flag == -6)
            {
                printf("| Mean surface water: %f |\n", wet_mean_surf);
                printf("| Total area: %f |\n", area_sum_loss);
                printf("\n|| Losing ( > mean) Patches ||\n");
                printf("|  ID        | Include |   Area   | Est.Loses (m3) | potential | \n");
            }

            /*--------------------------------------------------------------*/
            /*  loop 2, loop through losing (>mean) patches                 */
            /*--------------------------------------------------------------*/

            for (i = 0; i < zone[0].patch_families[pf][0].num_patches_in_fam; i++)
            {
                // if - included and surface is > mean (losers)
                if ( incl_surf[i] > 0 && (wet_surf[i] - wet_mean_surf) > zone[0].patch_families[pf][0].patches[i][0].landuse_defaults[0][0].surface_routing_threshold)
                {
                    dL[i] = (wet_surf[i] - wet_mean_surf) * zone[0].patch_families[pf][0].patches[i][0].area;
                }
                else if (incl_surf[i] >= 0 && -(wet_surf[i] - wet_mean_surf) > zone[0].patch_families[pf][0].patches[i][0].landuse_defaults[0][0].surface_routing_threshold && zone[0].patch_families[pf][0].patches[i][0].landuse_defaults[0][0].surf_g > 0)
                {
                // is a gaining patch
                    incl_surf[i] = 2;
                    dL[i] = 0;
                }
                else
                {
                    dL[i] = 0;
                } 
                
                dL_pot += dL[i];

                if (command_line[0].verbose_flag == -6)
                {
                    printf("| %4d |  %12d | %6.4f | %2.12f | %2.12f |\n",
                        zone[0].patch_families[pf][0].patches[i][0].ID,
                        incl_surf[i],
                        zone[0].patch_families[pf][0].patches[i][0].area,
                        dL[i],
                        dL_pot);
                }

            } // end loop 2

            if (verbose_flag == -6)
            {
                printf("| Total Estimated Losses | Actual: %f |\n", dL_pot);
                printf("|| Gaining ( < mean) Patches ||\n");
                printf("|  ID  | Include Surf |   Area   |  Gains (m3)  | potential  |\n");
            }

            /*--------------------------------------------------------------*/
            /*  loop 3, loop through gaining (<mean) patches              	*/
            /*--------------------------------------------------------------*/

            for (i = 0; i < zone[0].patch_families[pf][0].num_patches_in_fam; i++)
            {
                // if < mean wetness (gainers)
                if (incl_surf[i] == 2 && dL_pot > ZERO)
                {
                    if(zone[0].patch_families[pf][0].patches[i][0].landuse_defaults[0][0].surf_g==1){
                        dG[i] = (wet_mean_surf - wet_surf[i]) * zone[0].patch_families[pf][0].patches[i][0].area;
                    }else{   
                        dG[i] = 0;
                    }
                    // all of gaining water (dG) goes to, unless patch cannot gain
                    // surface detention store gain
                    if(zone[0].patch_families[pf][0].patches[i][0].landuse_defaults[0][0].surf_g==1){
                        
                        surf_trans = dG[i] / zone[0].patch_families[pf][0].patches[i][0].area;
                        zone[0].patch_families[pf][0].patches[i][0].surface_transfer = surf_trans;

                        dG_pot += zone[0].patch_families[pf][0].patches[i][0].surface_transfer*zone[0].patch_families[pf][0].patches[i][0].area;

                        area_sum_g += zone[0].patch_families[pf][0].patches[i][0].area;
                    }
                    
                    // area of all gaining patches 
                    
                }
                else
                {
                    dG[i] = 0;
                }

                if (command_line[0].verbose_flag == -6)
                {
                    printf("| %4d |  %12d | %6.4f | %2.12f | %2.12f |\n",
                        zone[0].patch_families[pf][0].patches[i][0].ID,
                        incl_surf[i],
                        zone[0].patch_families[pf][0].patches[i][0].area,
                        dG[i],
                        dG_pot);
                }
            }     // end loop 3

            
            

    /* RACHEL - I see an issue in that we are letting detention store water infiltation twice (e.g also in patch_daily_F etc.  */
    /* we could get around this by only infiltrating the additional water (eg. for gaining patches the surface_transfer */
    /* changed compute infiltration input to surface transfer instead of detention store*/

    printf("starting infiltration of det store after transfer\n"); 

            for (i = 0; i < zone[0].patch_families[pf][0].num_patches_in_fam; i++)
            {

            // take detention store losses from losing patch 
            if (incl_surf[i] == 1 & dL_pot>0) 
                {
                    // update loss values for gains less than potential, should have no impact if L actual = G actual
                    // if only two patches, and one has no gains (surf_g=0), then there should be no losses from the other (surface overflow)
                    // if three or more patches, and only one has no gains (surf_g=0), then need to figure out loops when not counting the no gains patch 

                    // distribute losses away from det stores
                    surf_trans = -1 * dL[i] / zone[0].patch_families[pf][0].patches[i][0].area;
                    
                    zone[0].patch_families[pf][0].patches[i][0].surface_transfer = surf_trans;
                    zone[0].patch_families[pf][0].patches[i][0].detention_store += zone[0].patch_families[pf][0].patches[i][0].surface_transfer;

                }

            
            // if gaining patch and surface trasnfer is greater than 0, infiltrate surface transfer 
            /* add infiltration - copied 673 - 788 from compute_subsurface_routing p*/
                if (incl_surf[i]==2 & zone[0].patch_families[pf][0].patches[i][0].surface_transfer > ZERO) 
                    if (zone[0].patch_families[pf][0].patches[i][0].rootzone.depth > ZERO) {
                        infiltration = compute_infiltration(verbose_flag,
                            zone[0].patch_families[pf][0].patches[i][0].sat_deficit_z, 
                            zone[0].patch_families[pf][0].patches[i][0].rootzone.S,
                            zone[0].patch_families[pf][0].patches[i][0].Ksat_vertical,
                            zone[0].patch_families[pf][0].patches[i][0].soil_defaults[0][0].Ksat_0_v,
                            zone[0].patch_families[pf][0].patches[i][0].soil_defaults[0][0].mz_v,
                            zone[0].patch_families[pf][0].patches[i][0].soil_defaults[0][0].porosity_0,
                            zone[0].patch_families[pf][0].patches[i][0].soil_defaults[0][0].porosity_decay,
                            (zone[0].patch_families[pf][0].patches[i][0].surface_transfer), time_int,
                            zone[0].patch_families[pf][0].patches[i][0].soil_defaults[0][0].psi_air_entry);
                        } else {
                        infiltration = compute_infiltration(verbose_flag,
                            zone[0].patch_families[pf][0].patches[i][0].sat_deficit_z, 
                            zone[0].patch_families[pf][0].patches[i][0].S,
                            zone[0].patch_families[pf][0].patches[i][0].Ksat_vertical,
                            zone[0].patch_families[pf][0].patches[i][0].soil_defaults[0][0].Ksat_0_v,
                            zone[0].patch_families[pf][0].patches[i][0].soil_defaults[0][0].mz_v,
                            zone[0].patch_families[pf][0].patches[i][0].soil_defaults[0][0].porosity_0,
                            zone[0].patch_families[pf][0].patches[i][0].soil_defaults[0][0].porosity_decay,
                            (zone[0].patch_families[pf][0].patches[i][0].surface_transfer), time_int,
                            zone[0].patch_families[pf][0].patches[i][0].soil_defaults[0][0].psi_air_entry);
                        }
                    else
                        infiltration = 0.0;
                    /*--------------------------------------------------------------*/
                    /* added an surface N flux to surface N pool	and		*/
                    /* allow infiltration of surface N				*/
                    /*--------------------------------------------------------------*/
                    if ((grow_flag > 0) && (infiltration > ZERO)) {
                        zone[0].patch_families[pf][0].patches[i][0].soil_ns.DON += ((infiltration
                                / zone[0].patch_families[pf][0].patches[i][0].detention_store) * zone[0].patch_families[pf][0].patches[i][0].surface_DON);
                        zone[0].patch_families[pf][0].patches[i][0].soil_cs.DOC += ((infiltration
                                / zone[0].patch_families[pf][0].patches[i][0].detention_store) * zone[0].patch_families[pf][0].patches[i][0].surface_DOC);
                        zone[0].patch_families[pf][0].patches[i][0].soil_ns.nitrate += ((infiltration
                                / zone[0].patch_families[pf][0].patches[i][0].detention_store) * zone[0].patch_families[pf][0].patches[i][0].surface_NO3);
                        zone[0].patch_families[pf][0].patches[i][0].surface_NO3 -= ((infiltration
                                / zone[0].patch_families[pf][0].patches[i][0].detention_store) * zone[0].patch_families[pf][0].patches[i][0].surface_NO3);
                        zone[0].patch_families[pf][0].patches[i][0].soil_ns.sminn += ((infiltration
                                / zone[0].patch_families[pf][0].patches[i][0].detention_store) * zone[0].patch_families[pf][0].patches[i][0].surface_NH4);
                        zone[0].patch_families[pf][0].patches[i][0].surface_NH4 -= ((infiltration
                                / zone[0].patch_families[pf][0].patches[i][0].detention_store) * zone[0].patch_families[pf][0].patches[i][0].surface_NH4);
                        zone[0].patch_families[pf][0].patches[i][0].surface_DOC -= ((infiltration
                                / zone[0].patch_families[pf][0].patches[i][0].detention_store) * zone[0].patch_families[pf][0].patches[i][0].surface_DOC);
                        zone[0].patch_families[pf][0].patches[i][0].surface_DON -= ((infiltration
                                / zone[0].patch_families[pf][0].patches[i][0].detention_store) * zone[0].patch_families[pf][0].patches[i][0].surface_DON);
                    }

                    /*--------------------------------------------------------------*/
                    /*	Determine if the infifltration will fill up the unsat	*/
                    /*	zone or not.						*/
                    /*	We use the strict assumption that sat deficit is the	*/
                    /*	amount of water needed to saturate the soil.		*/
                    /*--------------------------------------------------------------*/

                    if (infiltration > zone[0].patch_families[pf][0].patches[i][0].sat_deficit - zone[0].patch_families[pf][0].patches[i][0].unsat_storage
                                    - zone[0].patch_families[pf][0].patches[i][0].rz_storage) {
                        /*--------------------------------------------------------------*/
                        /*		Yes the unsat zone will be filled so we may	*/
                        /*		as well treat the unsat_storage and infiltration*/
                        /*		as water added to the water table.		*/
                        /*--------------------------------------------------------------*/
                        zone[0].patch_families[pf][0].patches[i][0].sat_deficit -= (infiltration
                                + zone[0].patch_families[pf][0].patches[i][0].unsat_storage + zone[0].patch_families[pf][0].patches[i][0].rz_storage);
                        /*--------------------------------------------------------------*/
                        /*		There is no unsat_storage left.			*/
                        /*--------------------------------------------------------------*/
                        zone[0].patch_families[pf][0].patches[i][0].unsat_storage = 0.0;
                        zone[0].patch_families[pf][0].patches[i][0].rz_storage = 0.0;
                        zone[0].patch_families[pf][0].patches[i][0].field_capacity = 0.0;
                        zone[0].patch_families[pf][0].patches[i][0].rootzone.field_capacity = 0.0;
                    } else if ((zone[0].patch_families[pf][0].patches[i][0].sat_deficit
                            > zone[0].patch_families[pf][0].patches[i][0].rootzone.potential_sat)
                            && (infiltration
                                    > zone[0].patch_families[pf][0].patches[i][0].rootzone.potential_sat
                                            - zone[0].patch_families[pf][0].patches[i][0].rz_storage)) {
                        /*------------------------------------------------------------------------------*/
                        /*		Just add the infiltration to the rz_storage and unsat_storage	*/
                        /*------------------------------------------------------------------------------*/
                        zone[0].patch_families[pf][0].patches[i][0].unsat_storage += infiltration
                                - (zone[0].patch_families[pf][0].patches[i][0].rootzone.potential_sat
                                        - zone[0].patch_families[pf][0].patches[i][0].rz_storage);
                        zone[0].patch_families[pf][0].patches[i][0].rz_storage = zone[0].patch_families[pf][0].patches[i][0].rootzone.potential_sat;
                    }
                    /* Only rootzone layer saturated - perched water table case */
                    else if ((zone[0].patch_families[pf][0].patches[i][0].sat_deficit > zone[0].patch_families[pf][0].patches[i][0].rootzone.potential_sat)
                            && (infiltration
                                    <= zone[0].patch_families[pf][0].patches[i][0].rootzone.potential_sat
                                            - zone[0].patch_families[pf][0].patches[i][0].rz_storage)) {
                        /*--------------------------------------------------------------*/
                        /*		Just add the infiltration to the rz_storage	*/
                        /*--------------------------------------------------------------*/
                        zone[0].patch_families[pf][0].patches[i][0].rz_storage += infiltration;
                    }

                    else if ((zone[0].patch_families[pf][0].patches[i][0].sat_deficit
                            <= zone[0].patch_families[pf][0].patches[i][0].rootzone.potential_sat)
                            && (infiltration
                                    <= zone[0].patch_families[pf][0].patches[i][0].sat_deficit - zone[0].patch_families[pf][0].patches[i][0].rz_storage
                                            - zone[0].patch_families[pf][0].patches[i][0].unsat_storage)) {
                        zone[0].patch_families[pf][0].patches[i][0].rz_storage += zone[0].patch_families[pf][0].patches[i][0].unsat_storage;
                        /* transfer left water in unsat storage to rootzone layer */
                        zone[0].patch_families[pf][0].patches[i][0].unsat_storage = 0;
                        zone[0].patch_families[pf][0].patches[i][0].rz_storage += infiltration;
                        zone[0].patch_families[pf][0].patches[i][0].field_capacity = 0;
                    }

                    if (zone[0].patch_families[pf][0].patches[i][0].sat_deficit < 0.0) {
                        zone[0].patch_families[pf][0].patches[i][0].detention_store -= (zone[0].patch_families[pf][0].patches[i][0].sat_deficit
                                - zone[0].patch_families[pf][0].patches[i][0].unsat_storage);
                        zone[0].patch_families[pf][0].patches[i][0].sat_deficit = 0.0;
                        zone[0].patch_families[pf][0].patches[i][0].unsat_storage = 0.0;
                    }

                    zone[0].patch_families[pf][0].patches[i][0].detention_store -= infiltration;
            }
            
        }

    } // end patch family loop

    return;
}
