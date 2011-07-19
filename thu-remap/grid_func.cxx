#include "grid.h"
#include "grid_func.h"
#include "remap_cfg.h"
#include "models_cfg.h"
#include "data_cfg.h"
#include <string.h>
#include <stdio.h>

grid_group *cpl_grids;

grid_group::grid_group()
{
	size = 16;
	count = 0;
	grids = new grid *[size];
}

grid_group::~grid_group()
{
	for (int i = 0; i < count; i ++)
		delete grids[i];
	delete [] grids;
}

void grid_group::realloc()
{
	grid **new_buf;
	int i;

	size *= 2;
	new_buf = new grid * [size];
	for (i = 0; i < count; i ++)
		new_buf[i] = grids[i];
	delete [] grids;
	grids = new_buf;
}

int grid_group::find_grid(char *name)
{
	for (int i = 0; i < count; i ++)
		if (strcmp(name, grids[i]->get_grid_name()) == 0)
			return i;
	return -1;
}

grid* grid_group::get_grid(char *name)
{
	int indx = find_grid(name);
	if (indx == -1) {
		printf("wrong name of grid %s\n", name);
		return NULL;
	}	
	return grids[indx];
}

void grid_group::generate_grid(char *name, 
                                          char *file, 
                                          int model_set_id,
                                          int file_type)
{
	printf("grid %s\n", name);
	if (find_grid(name) != -1) {
		printf("grid <%s> has been generated\n", name);
		return;
	}
	grids[count ++] = new grid(file, name, model_set_id, file_type);

	if (count >= size)
		realloc();
}

void generate_all_cpl_grids()
{
	int i;
	
	cpl_grids = new grid_group();

	if (cfg_atm_output_data_bundle.num_fields_2D > 0 || 
	    cfg_atm_input_data_bundle.num_fields_2D > 0)
		for (i = 0; i < atm_model_set.num_models; i ++) 
			cpl_grids->generate_grid(atm_model_set.models_cfg[i].grid_2D.name_grid,
			                                  atm_model_set.models_cfg[i].grid_2D.name_file,  
			                                  atm_set_id, 
			                                  atm_model_set.models_cfg[i].grid_2D.type_file);

	if (cfg_atm_output_data_bundle.num_fields_3D > 0 ||
	    cfg_atm_input_data_bundle.num_fields_3D > 0) 
		for (i = 0; i < atm_model_set.num_models; i ++) 
			cpl_grids->generate_grid(atm_model_set.models_cfg[i].grid_3D.name_grid,
			                                  atm_model_set.models_cfg[i].grid_3D.name_file,  
			                                  atm_set_id, 
			                                  atm_model_set.models_cfg[i].grid_3D.type_file);
	
	if (cfg_lnd_output_data_bundle.num_fields_2D > 0 ||
	    cfg_lnd_input_data_bundle.num_fields_2D > 0)
		for (i = 0; i < lnd_model_set.num_models; i ++) 
			cpl_grids->generate_grid(lnd_model_set.models_cfg[i].grid_2D.name_grid,
			                                  lnd_model_set.models_cfg[i].grid_2D.name_file,  
			                                  atm_set_id, 
			                                  lnd_model_set.models_cfg[i].grid_2D.type_file);
		
	if (cfg_lnd_output_data_bundle.num_fields_3D > 0 ||
	    cfg_lnd_input_data_bundle.num_fields_3D > 0)
		for (i = 0; i < lnd_model_set.num_models; i ++) 
			cpl_grids->generate_grid(lnd_model_set.models_cfg[i].grid_3D.name_grid,
			                                  lnd_model_set.models_cfg[i].grid_3D.name_file,  
			                                  atm_set_id, 
			                                  lnd_model_set.models_cfg[i].grid_3D.type_file);

	if (cfg_ocn_output_data_bundle.num_fields_2D > 0 ||
	    cfg_ocn_input_data_bundle.num_fields_2D > 0)
		for (i = 0; i < ocn_model_set.num_models; i ++) 
			cpl_grids->generate_grid(ocn_model_set.models_cfg[i].grid_2D.name_grid,
			                                  ocn_model_set.models_cfg[i].grid_2D.name_file,  
			                                  atm_set_id, 
			                                  ocn_model_set.models_cfg[i].grid_2D.type_file);
		
	if (cfg_ocn_output_data_bundle.num_fields_3D > 0 ||
	    cfg_ocn_input_data_bundle.num_fields_3D > 0)
		for (i = 0; i < ocn_model_set.num_models; i ++) 
			cpl_grids->generate_grid(ocn_model_set.models_cfg[i].grid_3D.name_grid,
			                                  ocn_model_set.models_cfg[i].grid_3D.name_file,  
			                                  atm_set_id, 
			                                  ocn_model_set.models_cfg[i].grid_3D.type_file);

	if (cfg_sice_output_data_bundle.num_fields_2D > 0 ||
	    cfg_sice_input_data_bundle.num_fields_2D > 0)
		for (i = 0; i < sice_model_set.num_models; i ++) 
			cpl_grids->generate_grid(sice_model_set.models_cfg[i].grid_2D.name_grid,
			                                  sice_model_set.models_cfg[i].grid_2D.name_file,  
			                                  atm_set_id, 
			                                  sice_model_set.models_cfg[i].grid_2D.type_file);
		
	if (cfg_sice_output_data_bundle.num_fields_3D > 0 ||
	    cfg_sice_input_data_bundle.num_fields_3D > 0)
		for (i = 0; i < sice_model_set.num_models; i ++)
			cpl_grids->generate_grid(sice_model_set.models_cfg[i].grid_3D.name_grid,
			                                  sice_model_set.models_cfg[i].grid_3D.name_file,  
			                                  atm_set_id, 
			                                  sice_model_set.models_cfg[i].grid_3D.type_file);
}

void finalize_all_cpl_grids()
{
	delete cpl_grids;
}



