/*--------------------------------------------------------------*/
/* 																*/
/*					execute_daily_output_event					*/
/*																*/
/*	execute_daily_output_event - outputs daily data			*/
/*																*/
/*	NAME														*/
/*	execute_daily_output_event - outputs daily data 	.		*/
/*																*/
/*	SYNOPSIS													*/
/*	void	execute_daily_output_event(						*/
/*					struct	world_object	*world,				*/
/*					struct	command_line_object *command_line,	*/
/*					struct	date	date,  						*/
/*					struct	world_output_file_object 	*outfile)*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	outputs spatial structure according to commandline			*/
/*	specifications to specific files, for daily info			*/
/*	a spatial structure is output only if its appropriate		*/
/*	option flag has been set and its ID matches a 				*/
/*	specified ID, or -999 which indicates all					*/
/*	units at that spatial scale are to be output				*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	We only permit one fileset per spatial modelling level.     */
/*	Each fileset has one file for each timestep.  				*/
/*																*/
/*	March 14, 1997	- 	RAF				*/
/*	Allowed output patch to also output the moss strata if	*/
/*		moss is present.				*/
/*--------------------------------------------------------------*/
#include <stdio.h>

#include "rhessys.h"
#include "output_filter/output_filter_output.h"

void execute_daily_output_event(struct world_object *world,
		struct command_line_object *command_line, struct date date,
		struct world_output_file_object *outfile) {
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
	void output_basin(int, struct basin_object*, struct date, FILE*);

	void output_hillslope(int, struct hillslope_object*, struct date, FILE*);

	void output_zone(int, int, struct zone_object*, struct date, FILE*);

	void output_patch(int, int, int, struct patch_object*, struct zone_object*,
			struct date, FILE*);

	void output_canopy_stratum(int, int, int, int, struct canopy_strata_object*,
			struct date, FILE*);

	void output_fire(int, int, int, int, struct canopy_strata_object*,
			struct date, FILE*);

	void output_shadow_strata(int, int, int, int, struct canopy_strata_object*,
			struct date, FILE*);

	void output_stream_routing(struct stream_network_object*, struct date,
			FILE*);
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int basinID, hillID, patchID, zoneID, stratumID, reachID;
	int b, h, p, z, c, s;
	/*--------------------------------------------------------------*/
	/*	check to see if there are any print options					*/
	/*--------------------------------------------------------------*/

	if (command_line[0].output_filter_flag) {
		/*----------------------------------------------------------------------*/
		/*  Handle output filter output.                                        */
		/*----------------------------------------------------------------------*/
		bool of_result = true;
		char *of_error = (char *)calloc(MAXSTR, sizeof(char));
		of_result = output_filter_output_daily(of_error, MAXSTR, command_line->verbose_flag,
				date, command_line->output_filter);
		if (!of_result) {
			fprintf(stderr, "output_filter_output_daily failed with error: %s\n", of_error);
			exit(EXIT_FAILURE);
		}
	} else if (command_line[0].legacy_output_flag) {
		/*----------------------------------------------------------------------*/
		/*  Handle legacy output: Make up the prefix for the output files.      */
		/*----------------------------------------------------------------------*/

		/*--------------------------------------------------------------*/
		/*	output stream_routing												*/
		/*--------------------------------------------------------------*/
		for (b = 0; b < world[0].num_basin_files; ++b) {
			for (s = 0; s < world[0].basins[b][0].stream_list.num_reaches;
					++s) {
				/*--------------------------------------------------------------*/
				/*	Construct the stream output files.							*/
				/*--------------------------------------------------------------*/
				if (command_line[0].stro != NULL) {
					reachID = command_line[0].stro->reachID;
					if ((world[0].basins[b][0].stream_list.stream_network[s].reach_ID
							== reachID) || (reachID == -999)) {
						output_stream_routing(
								&(world[0].basins[b]->stream_list.stream_network[s]),
								date, outfile->stream_routing->daily);
					}
				}
			}

		}
		if ((command_line[0].b != NULL) || (command_line[0].h != NULL)
				|| (command_line[0].z != NULL) || (command_line[0].p != NULL)
				|| (command_line[0].c != NULL)) {
			/*--------------------------------------------------------------*/
			/*	output basins												*/
			/*--------------------------------------------------------------*/
			for (b = 0; b < world[0].num_basin_files; ++b) {
				/*--------------------------------------------------------------*/
				/*	Construct the basin output files.							*/
				/*--------------------------------------------------------------*/
				if (command_line[0].b != NULL) {
					basinID = command_line[0].b->basinID;
					if ((world[0].basins[b][0].ID == basinID)
							|| (basinID == -999))
						output_basin(command_line[0].routing_flag,
								world[0].basins[b], date,
								outfile->basin->daily);
				}
				/*--------------------------------------------------------------*/
				/*	check to see if there are any lower print options			*/
				/*--------------------------------------------------------------*/
				if ((command_line[0].h != NULL) || (command_line[0].z != NULL)
						|| (command_line[0].p != NULL)
						|| (command_line[0].c != NULL)) {
					/*--------------------------------------------------------------*/
					/*	output hillslopes 											*/
					/*--------------------------------------------------------------*/
					for (h = 0; h < world[0].basins[b][0].num_hillslopes; ++h) {
						/*-----------------------------------------------------------*/
						/*	Construct the hillslope output files.						*/
						/*-----------------------------------------------------------*/
						if (command_line[0].h != NULL) {
							basinID = command_line[0].h->basinID;
							hillID = command_line[0].h->hillID;
							if ((world[0].basins[b][0].ID == basinID)
									|| (basinID == -999))
								if ((world[0].basins[b][0].hillslopes[h][0].ID
										== hillID) || (hillID == -999))
									output_hillslope(world[0].basins[b][0].ID,
											world[0].basins[b]->hillslopes[h],
											date, outfile->hillslope->daily);
						}
						/*-------------------------------------------------------------*/
						/*	check to see if there are any lower print options			*/
						/*-------------------------------------------------------------*/
						if ((command_line[0].z != NULL)
								|| (command_line[0].p != NULL)
								|| (command_line[0].c != NULL)) {
							/*---------------------------------------------------------*/
							/*	output zones												*/
							/*---------------------------------------------------------*/
							for (z = 0;
									z
											< world[0].basins[b][0].hillslopes[h][0].num_zones;
									++z) {
								/*------------------------------------------------------*/
								/*	Construct the zone output files.						  */
								/*-------------------------------------------------------*/
								if (command_line[0].z != NULL) {
									basinID = command_line[0].z->basinID;
									hillID = command_line[0].z->hillID;
									zoneID = command_line[0].z->zoneID;
									if ((world[0].basins[b][0].ID == basinID)
											|| (basinID == -999))
										if ((world[0].basins[b][0].hillslopes[h][0].ID
												== hillID) || (hillID == -999))
											if ((world[0].basins[b][0].hillslopes[h][0].zones[z][0].ID
													== zoneID)
													|| (zoneID == -999))
												output_zone(
														world[0].basins[b][0].ID,
														world[0].basins[b][0].hillslopes[h][0].ID,
														world[0].basins[b]->hillslopes[h]->zones[z],
														date,
														outfile->zone->daily);
								}
								/*-------------------------------------------------------*/
								/*	check to see if there are any lower print options		*/
								/*-------------------------------------------------------*/
								if ((command_line[0].p != NULL)
										|| (command_line[0].c != NULL)) {
									/*----------------------------------------------------*/
									/*	output patches 												*/
									/*---------------------------------------------------*/
									for (p = 0;
											p
													< world[0].basins[b][0].hillslopes[h][0].zones[z][0].num_patches;
											++p) {
										/*-------------------------------------------------*/
										/*	Construct the patch output files.					*/
										/*-------------------------------------------------*/
										if (command_line[0].p != NULL) {
											basinID =
													command_line[0].p->basinID;
											hillID = command_line[0].p->hillID;
											zoneID = command_line[0].p->zoneID;
											patchID =
													command_line[0].p->patchID;
											if ((world[0].basins[b][0].ID
													== basinID)
													|| (basinID == -999))
												if ((world[0].basins[b][0].hillslopes[h][0].ID
														== hillID)
														|| (hillID == -999))
													if ((world[0].basins[b][0].hillslopes[h][0].zones[z][0].ID
															== zoneID)
															|| (zoneID == -999))
														if ((world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].ID
																== patchID)
																|| (patchID
																		== -999)) {
															output_patch(
																	world[0].basins[b]->ID,
																	world[0].basins[b]->hillslopes[h]->ID,
																	world[0].basins[b]->hillslopes[h]->zones[z]->ID,
																	world[0].basins[b]->hillslopes[h]->zones[z]->patches[p],
																	world[0].basins[b]->hillslopes[h]->zones[z],
																	date,
																	outfile->patch->daily);
														}
										}
										/*------------------------------------------------*/
										/*	Construct the canopy_stratum output files		  */
										/*------------------------------------------------*/
										if (command_line[0].c != NULL) {
											/*----------------------------------------------*/
											/*	output canopy stratum 								*/
											/*----------------------------------------------*/
											for (c = 0;
													c
															< world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].num_canopy_strata;
													++c) {
												basinID =
														command_line[0].c->basinID;
												hillID =
														command_line[0].c->hillID;
												zoneID =
														command_line[0].c->zoneID;
												patchID =
														command_line[0].c->patchID;
												stratumID =
														command_line[0].c->stratumID;
												if ((world[0].basins[b][0].ID
														== basinID)
														|| (basinID == -999))
													if ((world[0].basins[b][0].hillslopes[h][0].ID
															== hillID)
															|| (hillID == -999))
														if ((world[0].basins[b][0].hillslopes[h][0].zones[z][0].ID
																== zoneID)
																|| (zoneID
																		== -999))
															if ((world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].ID
																	== patchID)
																	|| (patchID
																			== -999))
																if ((world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].canopy_strata[c][0].ID
																		== stratumID)
																		|| (stratumID
																				== -999)) {
																	output_canopy_stratum(
																			world[0].basins[b][0].ID,
																			world[0].basins[b][0].hillslopes[h][0].ID,
																			world[0].basins[b][0].hillslopes[h][0].zones[z][0].ID,
																			world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].ID,
																			world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->canopy_strata[c],
																			date,
																			outfile->canopy_stratum->daily);
																}
											} /* end stratum (c) for loop */
										} /* end if options */

										/*------------------------------------------------*/
										/*	Construct the fire output files		  */
										/*------------------------------------------------*/
										/*								if ( command_line[0].f != NULL ){
										 /*----------------------------------------------*/
										/*	output fire 								*/
										/*----------------------------------------------*/
										/*									for(c=0;
										 c < world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].num_canopy_strata;
										 ++c){
										 basinID = command_line[0].f->basinID;
										 hillID = command_line[0].f->hillID;
										 zoneID = command_line[0].f->zoneID;
										 patchID = command_line[0].f->patchID;
										 stratumID = command_line[0].f->stratumID;
										 if (( world[0].basins[b][0].ID == basinID)
										 || (basinID == -999))
										 if (( world[0].basins[b][0].hillslopes[h][0].ID == hillID)
										 || (hillID == -999))
										 if (( world[0].basins[b][0].hillslopes[h][0].zones[z][0].ID == zoneID)
										 || (zoneID == -999))
										 if (( world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].ID == patchID)
										 ||	(patchID == -999))
										 if (( world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].canopy_strata[c][0].ID == stratumID)
										 || (stratumID == -999)) {
										 output_fire(
										 world[0].basins[b][0].ID,
										 world[0].basins[b][0].hillslopes[h][0].ID,
										 world[0].basins[b][0].hillslopes[h][0].zones[z][0].ID,
										 world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].ID,
										 world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->canopy_strata[c],
										 date, outfile->fire->daily);
										 }
										 } /* end fire (f) for loop */
										//								} /* end if options */
									} /* end patch (p) for loop */
								} /* end if options */
							} /* end zone (z) for  loop*/
						} /* end if options */
					} /* end hillslope (h) for loop */
				} /* end if options */
			} /* end basin (b) for loop */
		} /* end if options */

	}
	return;
} /*end execute_daily_output_event*/
