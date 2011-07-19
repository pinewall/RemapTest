#include "data_func.h"
#include "data_bundle.h"
#include "name_cfg.h"
#include "models_cfg.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "grid_func.h"


data_bundle** data_bundles_atm_output_model;
data_bundle** data_bundles_atm_output_cpl;
data_bundle* data_bundles_atm_output_aggr;
data_bundle** data_bundles_atm_input_model;
data_bundle** data_bundles_atm_input_cpl;
data_bundle* data_bundles_atm_input_aggr;
data_bundle** data_bundles_lnd_output_model;
data_bundle** data_bundles_lnd_output_cpl;
data_bundle* data_bundles_lnd_output_aggr;
data_bundle** data_bundles_lnd_input_model;
data_bundle** data_bundles_lnd_input_cpl;
data_bundle* data_bundles_lnd_input_aggr;
data_bundle** data_bundles_ocn_output_model;
data_bundle** data_bundles_ocn_output_cpl;
data_bundle* data_bundles_ocn_output_aggr;
data_bundle** data_bundles_ocn_input_model;
data_bundle** data_bundles_ocn_input_cpl;
data_bundle* data_bundles_ocn_input_aggr;
data_bundle** data_bundles_sice_output_model;
data_bundle** data_bundles_sice_output_cpl;
data_bundle* data_bundles_sice_output_aggr;
data_bundle** data_bundles_sice_input_model;
data_bundle** data_bundles_sice_input_cpl;
data_bundle* data_bundles_sice_input_aggr;

void generate_all_data_bundles()
{
	int i;
	if (atm_model_set.num_models > 0) {
		data_bundles_atm_output_model = new data_bundle * [atm_model_set.num_models];
		data_bundles_atm_output_cpl = new data_bundle * [atm_model_set.num_models];
		data_bundles_atm_output_aggr = new data_bundle(&cfg_atm_output_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);		
		data_bundles_atm_input_model = new data_bundle * [atm_model_set.num_models];
		data_bundles_atm_input_cpl = new data_bundle * [atm_model_set.num_models];
		data_bundles_atm_input_aggr = new data_bundle(&cfg_atm_input_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);
		for (i = 0; i < atm_model_set.num_models; i ++) {
			data_bundles_atm_output_model[i] = new data_bundle(&cfg_atm_output_data_bundle,
			                                                         cpl_grids->get_grid(atm_model_set.models_cfg[i].grid_2D.name_grid)->get_num_points(),
			                                                         0);
			data_bundles_atm_input_model[i] = new data_bundle(&cfg_atm_input_data_bundle,
			                                                         cpl_grids->get_grid(atm_model_set.models_cfg[i].grid_2D.name_grid)->get_num_points(),
			                                                         0);
			data_bundles_atm_output_cpl[i] = new data_bundle(&cfg_atm_output_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);
			data_bundles_atm_input_cpl[i] = new data_bundle(&cfg_atm_input_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);
		}
	}

	if (lnd_model_set.num_models > 0) {
		data_bundles_lnd_output_model = new data_bundle * [lnd_model_set.num_models];
		data_bundles_lnd_output_cpl = new data_bundle * [lnd_model_set.num_models];
		data_bundles_lnd_output_aggr = new data_bundle(&cfg_lnd_output_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);		
		data_bundles_lnd_input_model = new data_bundle * [lnd_model_set.num_models];
		data_bundles_lnd_input_cpl = new data_bundle * [lnd_model_set.num_models];
		data_bundles_lnd_input_aggr = new data_bundle(&cfg_lnd_input_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);
		for (i = 0; i < lnd_model_set.num_models; i ++) {
			data_bundles_lnd_output_model[i] = new data_bundle(&cfg_lnd_output_data_bundle,
			                                                         cpl_grids->get_grid(lnd_model_set.models_cfg[i].grid_2D.name_grid)->get_num_points(),
			                                                         0);
			data_bundles_lnd_input_model[i] = new data_bundle(&cfg_lnd_input_data_bundle,
			                                                         cpl_grids->get_grid(lnd_model_set.models_cfg[i].grid_2D.name_grid)->get_num_points(),
			                                                         0);
			data_bundles_lnd_output_cpl[i] = new data_bundle(&cfg_lnd_output_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);
			data_bundles_lnd_input_cpl[i] = new data_bundle(&cfg_lnd_input_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);
		}
	}

	if (ocn_model_set.num_models > 0) {
		data_bundles_ocn_output_model = new data_bundle * [ocn_model_set.num_models];
		data_bundles_ocn_output_cpl = new data_bundle * [ocn_model_set.num_models];
		data_bundles_ocn_output_aggr = new data_bundle(&cfg_ocn_output_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);		
		data_bundles_ocn_input_model = new data_bundle * [ocn_model_set.num_models];
		data_bundles_ocn_input_cpl = new data_bundle * [ocn_model_set.num_models];
		data_bundles_ocn_input_aggr = new data_bundle(&cfg_ocn_input_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);
		for (i = 0; i < ocn_model_set.num_models; i ++) {
			data_bundles_ocn_output_model[i] = new data_bundle(&cfg_ocn_output_data_bundle,
			                                                         cpl_grids->get_grid(ocn_model_set.models_cfg[i].grid_2D.name_grid)->get_num_points(),
			                                                         0);
			data_bundles_ocn_input_model[i] = new data_bundle(&cfg_ocn_input_data_bundle,
			                                                         cpl_grids->get_grid(ocn_model_set.models_cfg[i].grid_2D.name_grid)->get_num_points(),
			                                                         0);
			data_bundles_ocn_output_cpl[i] = new data_bundle(&cfg_ocn_output_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);
			data_bundles_ocn_input_cpl[i] = new data_bundle(&cfg_ocn_input_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);
		}
	}

	if (sice_model_set.num_models > 0) {
		data_bundles_sice_output_model = new data_bundle * [sice_model_set.num_models];
		data_bundles_sice_output_cpl = new data_bundle * [sice_model_set.num_models];
		data_bundles_sice_output_aggr = new data_bundle(&cfg_sice_output_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);		
		data_bundles_sice_input_model = new data_bundle * [sice_model_set.num_models];
		data_bundles_sice_input_cpl = new data_bundle * [sice_model_set.num_models];
		data_bundles_sice_input_aggr = new data_bundle(&cfg_sice_input_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);
		for (i = 0; i < sice_model_set.num_models; i ++) {
			data_bundles_sice_output_model[i] = new data_bundle(&cfg_sice_output_data_bundle,
			                                                         cpl_grids->get_grid(sice_model_set.models_cfg[i].grid_2D.name_grid)->get_num_points(),
			                                                         0);
			data_bundles_sice_input_model[i] = new data_bundle(&cfg_sice_input_data_bundle,
			                                                         cpl_grids->get_grid(sice_model_set.models_cfg[i].grid_2D.name_grid)->get_num_points(),
			                                                         0);
			data_bundles_sice_output_cpl[i] = new data_bundle(&cfg_sice_output_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);
			data_bundles_sice_input_cpl[i] = new data_bundle(&cfg_sice_input_data_bundle,
			                                                         cpl_grids->get_grid(NAME_CPL_GRID_2D)->get_num_points(),
			                                                         0);
		}
	}
}

void finalize_all_data_bundles()
{
	int i;
	
	if (atm_model_set.num_models > 0) {
		for (i = 0; i < atm_model_set.num_models; i ++) {
			delete data_bundles_atm_output_model[i];
			delete data_bundles_atm_input_model[i];
			delete data_bundles_atm_output_cpl[i];
			delete data_bundles_atm_input_cpl[i];
		}
		delete [] data_bundles_atm_output_model;
		delete [] data_bundles_atm_output_cpl;
		delete data_bundles_atm_output_aggr;		
		delete [] data_bundles_atm_input_model;
		delete [] data_bundles_atm_input_cpl;
		delete data_bundles_atm_input_aggr;
	}

	if (lnd_model_set.num_models > 0) {
		for (i = 0; i < lnd_model_set.num_models; i ++) {
			delete data_bundles_lnd_output_model[i];
			delete data_bundles_lnd_input_model[i];
			delete data_bundles_lnd_output_cpl[i];
			delete data_bundles_lnd_input_cpl[i];
		}
		delete [] data_bundles_lnd_output_model;
		delete [] data_bundles_lnd_output_cpl;
		delete data_bundles_lnd_output_aggr;		
		delete [] data_bundles_lnd_input_model;
		delete [] data_bundles_lnd_input_cpl;
		delete data_bundles_lnd_input_aggr;
	}

	if (ocn_model_set.num_models > 0) {
		for (i = 0; i < ocn_model_set.num_models; i ++) {
			delete data_bundles_ocn_output_model[i];
			delete data_bundles_ocn_input_model[i];
			delete data_bundles_ocn_output_cpl[i];
			delete data_bundles_ocn_input_cpl[i];
		}
		delete [] data_bundles_ocn_output_model;
		delete [] data_bundles_ocn_output_cpl;
		delete data_bundles_ocn_output_aggr;		
		delete [] data_bundles_ocn_input_model;
		delete [] data_bundles_ocn_input_cpl;
		delete data_bundles_ocn_input_aggr;
	}

	if (sice_model_set.num_models > 0) {
		for (i = 0; i < sice_model_set.num_models; i ++) {
			delete data_bundles_sice_output_model[i];
			delete data_bundles_sice_input_model[i];
			delete data_bundles_sice_output_cpl[i];
			delete data_bundles_sice_input_cpl[i];
		}
		delete [] data_bundles_sice_output_model;
		delete [] data_bundles_sice_output_cpl;
		delete data_bundles_sice_output_aggr;		
		delete [] data_bundles_sice_input_model;
		delete [] data_bundles_sice_input_cpl;
		delete data_bundles_sice_input_aggr;
	}
}

int copy_in_atm_set_data(char *input_data)
{
	int i, copy_size;

	for (i = 0; i < atm_model_set.num_models; i ++)
		copy_size = data_bundles_atm_output_model[i]->copy_in_data(input_data);

	printf("copy size is %d\n", copy_size);
	return copy_size;
}

int copy_in_lnd_set_data(char *input_data)
{
	int i, copy_size;

	for (i = 0; i < lnd_model_set.num_models; i ++)
		copy_size = data_bundles_lnd_output_model[i]->copy_in_data(input_data);
	
	return copy_size;
}

int copy_in_ocn_set_data(char *input_data)
{
	int i, copy_size;

	for (i = 0; i < ocn_model_set.num_models; i ++)
		copy_size = data_bundles_ocn_output_model[i]->copy_in_data(input_data);
	
	return copy_size;
}

int copy_in_sice_set_data(char *input_data)
{
	int i, copy_size;

	for (i = 0; i < sice_model_set.num_models; i ++)
		copy_size = data_bundles_sice_output_model[i]->copy_in_data(input_data);
	
	return copy_size;
}

int copy_in_model_set_data(char *name_set, char *input_data)
{
	if (strcmp(name_set, NAME_ATM_SET) == 0)
		return copy_in_atm_set_data(input_data);
	else if (strcmp(name_set, NAME_LND_SET) == 0)
		return copy_in_lnd_set_data(input_data);
	else if (strcmp(name_set, NAME_OCN_SET) == 0)
		return copy_in_ocn_set_data(input_data);
	else if (strcmp(name_set, NAME_SICE_SET))
		return copy_in_sice_set_data(input_data);
	else {
		printf("error of input name set\n");
		return 0;
	}
}


int copy_out_atm_set_data(char *output_data)
{
	int i, copy_size;

	for (i = 0; i < atm_model_set.num_models; i ++)
		copy_size = data_bundles_atm_output_cpl[i]->copy_out_data(output_data);

	printf("copy size is %d\n", copy_size);
	return copy_size;
}

int copy_out_lnd_set_data(char *output_data)
{
	int i, copy_size;

	for (i = 0; i < lnd_model_set.num_models; i ++)
		copy_size = data_bundles_lnd_output_cpl[i]->copy_out_data(output_data);
	
	return copy_size;
}

int copy_out_ocn_set_data(char *output_data)
{
	int i, copy_size;

	for (i = 0; i < ocn_model_set.num_models; i ++)
		copy_size = data_bundles_ocn_output_cpl[i]->copy_out_data(output_data);
	
	return copy_size;
}

int copy_out_sice_set_data(char *output_data)
{
	int i, copy_size;

	for (i = 0; i < sice_model_set.num_models; i ++)
		copy_size = data_bundles_sice_output_cpl[i]->copy_out_data(output_data);
	
	return copy_size;
}

int copy_out_model_set_data(char *name_set, char *output_data)
{
	if (strcmp(name_set, NAME_ATM_SET) == 0)
		return copy_out_atm_set_data(output_data);
	else if (strcmp(name_set, NAME_LND_SET) == 0)
		return copy_out_lnd_set_data(output_data);
	else if (strcmp(name_set, NAME_OCN_SET) == 0)
		return copy_out_ocn_set_data(output_data);
	else if (strcmp(name_set, NAME_SICE_SET))
		return copy_out_sice_set_data(output_data);
	else {
		printf("error of output name set\n");
		return 0;
	}
}