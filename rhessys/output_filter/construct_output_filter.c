#include <string.h>

#include "rhessys.h"
#include "output_filter.h"
#include "index_struct_fields.h"
#include "output_filter/output_format_csv.h"


OutputFilter *parse(const char* input);
struct basin_object *find_basin(int basin_ID, struct world_object *world);
struct hillslope_object *find_hillslope_in_basin(int hillslope_ID, struct basin_object *basin);
struct zone_object *find_zone_in_hillslope(int zone_ID, struct hillslope_object *hillslope);
struct patch_object *find_patch(int patch_ID, int zone_ID, int hill_ID, struct basin_object *basin);


static bool returnWithError(char * const error, size_t error_len, char *error_mesg) {
	strncpy(error, error_mesg, error_len);
	return false;
}

static bool init_spatial_hierarchy_patch(OutputFilter *f,
		struct world_object * const w,
		bool verbose) {

	struct basin_object *b;

	if (f->patches == NULL) {
		fprintf(stderr, "init_spatial_hierarchy_patch: no patches defined but patch type was specified.");
		return false;
	}

	if (verbose) fprintf(stderr, "BEGIN init_spatial_hierarchy_patch\n");

	OutputFilterPatch *p = f->patches;
	while (p != NULL) {
		switch(p->output_patch_type) {
		case BASIN:
			if (verbose) {
				fprintf(stderr, "\tbasinID: %d\n", p->basinID);
			}
			p->basin = find_basin(p->basinID, w);
			if (p->basin == NULL) {
				fprintf(stderr, "init_spatial_hierarchy_patch: no basin with ID %d could be found.",
						p->basinID);
				return false;
			}
			break;
		case HILLSLOPE:
			if (verbose) {
				fprintf(stderr, "\tbasinID: %d, hillslopeID: %d\n",
						p->basinID, p->hillslopeID);
			}
			b = find_basin(p->basinID, w);
			if (b == NULL) {
				fprintf(stderr, "init_spatial_hierarchy_patch: basin %d could not be found, so could not locate hillslope %d.",
						p->basinID, p->hillslopeID);
				return false;
			}
			p->hill = find_hillslope_in_basin(p->hillslopeID, b);
			if (p->hill == NULL) {
				fprintf(stderr, "init_spatial_hierarchy_patch: hillslope %d could not be found in basin %d.",
						p->hillslopeID, p->basinID);
				return false;
			}
			break;
		case ZONE:
			if (verbose) {
				fprintf(stderr, "\tbasinID: %d, hillslopeID: %d, zoneID: %d\n",
						p->basinID, p->hillslopeID, p->zoneID);
			}
			b = find_basin(p->basinID, w);
			if (b == NULL) {
				fprintf(stderr, "init_spatial_hierarchy_patch: basin %d could not be found, so could not locate zone %d in hillslope %d.",
						p->basinID, p->zoneID, p->hillslopeID);
				return false;
			}
			struct hillslope_object *h = find_hillslope_in_basin(p->hillslopeID, b);
			if (h == NULL) {
				fprintf(stderr, "init_spatial_hierarchy_patch: hillslope %d could not be found in basin %d, so could not locate zone %d.",
						p->hillslopeID, p->basinID, p->zoneID);
				return false;
			}
			p->zone = find_zone_in_hillslope(p->zoneID, h);
			if (p->zone == NULL) {
				fprintf(stderr, "init_spatial_hierarchy_patch: zone %d could not be found in hillslope %d, basin %d.",
						p->zoneID, p->hillslopeID, p->basinID);
				return false;
			}
			break;
		case PATCH:
			if (verbose) {
				fprintf(stderr, "\tbasinID: %d, hillslopeID: %d, zoneID: %d, patchID: %d\n",
						p->basinID, p->hillslopeID, p->zoneID, p->patchID);
			}
			b = find_basin(p->basinID, w);
			if (b == NULL) {
				fprintf(stderr, "init_spatial_hierarchy_patch: basin %d could not be found, so could not locate patch %d in hillslope %d, patch %d.",
						p->basinID, p->patchID, p->hillslopeID, p->zoneID);
				return false;
			}
			p->patch = find_patch(p->patchID, p->zoneID, p->hillslopeID, b);
			if (p->patch == NULL) {
				fprintf(stderr, "init_spatial_hierarchy_patch: patch %d could not be found in zone %d, hillslope %d, basin %d.",
						p->patchID, p->zoneID, p->hillslopeID, p->basinID);
				return false;
			}
			break;
		}
		p = p->next;
	}

	if (verbose) fprintf(stderr, "END init_spatial_hierarchy_patch\n");

	return true;
}

static bool init_spatial_hierarchy(OutputFilter *f,
		struct world_object * const w,
		bool verbose) {
	switch(f->type) {
	case OUTPUT_FILTER_PATCH:
		return init_spatial_hierarchy_patch(f, w, verbose);
	case OUTPUT_FILTER_CANOPY_STRATA:
	default:
		fprintf(stderr, "init_spatial_hierarchy: output type %d is unknown or not yet implemented.", f->type);
		return false;
	}
	return true;
}

static bool init_output(OutputFilter *f) {
	switch(f->output->format) {
	case OUTPUT_TYPE_CSV:
		return output_format_csv_init(f);
	case OUTPUT_TYPE_NETCDF:
	default:
		fprintf(stderr, "init_output: output format type %d is unknown or not yet implemented.", f->output->format);
		return false;
	}
	return true;
}

static bool write_headers(OutputFilter *f) {
	switch(f->output->format) {
	case OUTPUT_TYPE_CSV:
		return output_format_csv_write_headers(f);
	case OUTPUT_TYPE_NETCDF:
	default:
		fprintf(stderr, "write_headers: output format type %d is unknown or not yet implemented.", f->output->format);
		return false;
	}
	return true;
}

bool construct_output_filter(char * const error, size_t error_len,
		struct command_line_object * const cmd,
		struct world_object * const world) {
	if (!cmd->output_filter_flag) {
		return returnWithError(error, error_len, "output_filter_flag is fall, not constructing output filter.");
	}

	// Parse output filter
	OutputFilter *filters = parse(cmd->output_filter_filename);
	if (filters == NULL) {
		return returnWithError(error, error_len, "unable to parse output filter.");
	}
	if (filters->parse_error) {
		return returnWithError(error, error_len, "output_filter_parser returned with an error.");
	}

	StructIndex_t *idx = index_struct_fields();

	bool status;
	for (OutputFilter *f = filters; f != NULL; f = f->next) {
		// Initialize spatial hierarchy objects referenced in filter
		status = init_spatial_hierarchy(f, world, cmd->verbose_flag);
		if (!status) {
			char *init_error = (char *)calloc(MAXSTR, sizeof(char));
			snprintf(init_error, MAXSTR, "unable to initialize output filter with path %s and filename %s.",
					f->output->path, f->output->filename);
			return returnWithError(error, error_len, init_error);
		}

		// TODO: Validate variables and write offsets to filter

		// Initialize output for filters
		if (f->output == NULL) {
			return returnWithError(error, error_len,
					"Output filter did not specify an output section.");
		}
		if (f->output->path == NULL) {
			return returnWithError(error, error_len,
					"Output filter did not specify output path.");
		}
		if (f->output->filename == NULL) {
			return returnWithError(error, error_len,
					"Output filter did not specify output filename.");
		}
		status = init_output(f);
		if (!status) {
			char *init_error = (char *)calloc(MAXSTR, sizeof(char));
			snprintf(init_error, MAXSTR, "unable to initialize output %s/%s for output filter.",
					f->output->path, f->output->filename);
			return returnWithError(error, error_len, init_error);
		}

		// Write header information for each output file
		status = write_headers(f);
		if (!status) {
			char *init_error = (char *)calloc(MAXSTR, sizeof(char));
			snprintf(init_error, MAXSTR, "unable to write headers for output %s/%s.",
					f->output->path, f->output->filename);
			return returnWithError(error, error_len, init_error);
		}
	}

	cmd->output_filter = filters;
	return true;
}
