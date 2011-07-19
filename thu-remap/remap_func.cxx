#include "dist_remap.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "remap_base.h"
#include "remap_cfg.h"
#include "grid.h"
#include "name_cfg.h"
#include "remap_func.h"
#include "data_bundle.h"
#include "models_cfg.h"
#include "grid_func.h"
#include "data_cfg.h"
#include "data_func.h"
#include "bilinear_remap.h"

//adding
#include "conserv_remap.h"


remap_operator_group *cpl_remap_operators;


remap_operator_group::remap_operator_group()
{
	count = 0;
	remap_operators = new Common_remap* [num_data_remaps];
	name_operators = new char * [num_data_remaps];

	for (int i = 0; i < num_data_remaps; i ++)
		name_operators[i] = new char [NAME_STR_SIZE];
}

remap_operator_group::~remap_operator_group()
{
	int i, j;
	
	for (i = 0; i < count; i ++) { 
		if (remap_operators[i] == NULL)
			continue;
		for (j = i + 1; j < count; j ++)
			if (remap_operators[j] == remap_operators[i])
				remap_operators[j] = NULL;
		delete remap_operators[i];
	}
	delete [] remap_operators;
}

int remap_operator_group::find_operator(char *grid1_name,
                                                       char *grid2_name, 
                                                       char *alg_name)
{
	
	for (int i = 0; i < count; i ++) 
		if (remap_operators[i]->search_remap_operator(grid1_name, grid2_name, alg_name) != NULL) {
		     	printf("find %s %s %s: %d\n", grid1_name, grid2_name, alg_name, count);
			return i;
		}
	return -1;
}

void remap_operator_group::generate_name(char *name_buf,
                                                            int *i_paras, 
                                                            double *d_paras)
{
	int i;
	char tmp_buf[1024];
	
	for (i = 0; i < REMAP_PARA_NUM; i ++) {
		sprintf(tmp_buf, "_i_%d", i_paras[i]);
		strcpy(name_buf + strlen(name_buf), tmp_buf);
	}
	for (i = 0; i < REMAP_PARA_NUM; i ++) {
		sprintf(tmp_buf, "_d_%f", d_paras[i]);
		strcpy(name_buf + strlen(name_buf), tmp_buf);
	}
}

void remap_operator_group::generate_operators()
{
	int i, indx;
	char name_operator[128];
	

	for (i = 0; i < num_data_remaps; i ++) {
		strcpy(name_operator, data_remaps_cfg[i].name_remap_alg);
		generate_name(name_operator, &(data_remaps_cfg[i].i_para1), &(data_remaps_cfg[i].d_para1));
		strcpy(name_operators[count], name_operator);

		if ((indx = find_operator(data_remaps_cfg[i].name_grid_src, data_remaps_cfg[i].name_grid_dst, name_operator)) != -1) {
			remap_operators[count ++] = remap_operators[indx];
			continue;
		}
		if (strcmp(NAME_DIST_REMAP_2D, data_remaps_cfg[i].name_remap_alg) == 0) {
			remap_operators[count] = new dist_remap_2D_operator(data_remaps_cfg[i].name_grid_src, 
			                                                                            data_remaps_cfg[i].name_grid_dst,
			                                                                            name_operator,
			                                                                            data_remaps_cfg[i].i_para1,
			                                                                            data_remaps_cfg[i].d_para1,
			                                                                            data_remaps_cfg[i].d_para2);                                                                                 
			remap_operators[count ++]->init_remap();
		}
        else if (strcmp(NAME_BILINEAR_REMAP_2D, data_remaps_cfg[i].name_remap_alg) == 0) {
			remap_operators[count] = new bilinear_remap_2D_operator(data_remaps_cfg[i].name_grid_src, 
			                                                                            data_remaps_cfg[i].name_grid_dst,
			                                                                            name_operator);                                                                                 
			remap_operators[count ++]->init_remap();
		}
		
        else if (strcmp(NAME_CONSERV_REMAP_2D, data_remaps_cfg[i].name_remap_alg) == 0) {
			remap_operators[count] = new conserv_remap_2D_operator(data_remaps_cfg[i].name_grid_src, 
			                                                                            data_remaps_cfg[i].name_grid_dst,
			                                                                            name_operator);                                                                            
			remap_operators[count ++]->init_remap();
		}
		
        else if (strcmp(NAME_BICUBIC_REMAP_2D, data_remaps_cfg[i].name_remap_alg) == 0) {

		}

        else {
            printf("Unknown Algorithm name!\n");
            printf("\t[bilinear_remap_2D]\n");
            printf("\t[distwgt_remap_2D]\n");
            printf("\t[conserv_remap_2D]\n");
            printf("\t[...]\n");
            return;
        }

		// ...
	}
}

void remap_operator_group::do_remap(char *name_data_field, 
                                                     char *name_grid_src, 
                                                     char *name_grid_dst, 
                                                     double *data_src, 
                                                     double *data_dst)
{
	int indx;

		printf("ok %s %s %s\n", name_data_field, name_grid_src, name_grid_dst);
	for (indx = 0; indx < num_data_remaps; indx ++)
		if (strcmp(data_remaps_cfg[indx].name_data_field, name_data_field) == 0 &&
		    strcmp(data_remaps_cfg[indx].name_grid_src, name_grid_src) == 0 &&
		    strcmp(data_remaps_cfg[indx].name_grid_dst, name_grid_dst) == 0)
			break;

	if (indx == num_data_remaps)
		printf("error: the config information of remap %s from %s to %s is needed\n", name_data_field, name_grid_src, name_grid_dst);
		
	remap_operators[indx]->cal_remap(data_src, data_dst);    
}

void generate_remap_operators()
{
	cpl_remap_operators = new remap_operator_group();
	cpl_remap_operators->generate_operators();
}

void do_remap(char *name_data_field, 
                     char *name_grid_src,
                     char *name_grid_dst, 
                     double *data_src, 
                     double *data_dst)
{
	cpl_remap_operators->do_remap(name_data_field,
	                                             name_grid_src, 
	                                             name_grid_dst, 
	                                             data_src, 
	                                             data_dst);
}

void finalize_remap_operators()
{
	delete cpl_remap_operators;
}

void remap_data_bundle(data_bundle *bundle_src,
                                 data_bundle *bundle_dst,
                                 int nfields_2d,
                                 int nfields_3d,
                                 char *name_grid2d_src,
                                 char *name_grid2d_dst,
                                 char *name_grid3d_src,
                                 char *name_grid3d_dst)
{
	int i, nfields;
	double *data_src, *data_dst;
	char *name_field;

	for (i = 0; i < nfields_2d; i ++) {
		data_src = (double *) bundle_src->get_field_data(i);
		data_dst = (double *) bundle_dst->get_field_data(i);
		name_field = bundle_src->get_field_name(i);
		cpl_remap_operators->do_remap(name_field, 
		                                                name_grid2d_src, 
		                                                name_grid2d_dst, 
		                                                data_src, 
		                                                data_dst);
	}
}

void remap_model_set_model2cpl(data_bundle **bundles_src,
                                             data_bundle **bundles_dst,
                                             model_set_cfg *cfg_model_set,
                                             data_bundle_cfg *cfg_data_bundle)
{
	int i;
	int nfields_2d, nfields_3d;

	nfields_2d = cfg_data_bundle->num_fields_2D;
	nfields_3d = cfg_data_bundle->num_fields_3D;

	for (i = 0; i < cfg_model_set->num_models; i ++) 
		remap_data_bundle(bundles_src[i], 
		                               bundles_dst[i], 
		                               nfields_2d, 
		                               nfields_3d, 
		                               cfg_model_set->models_cfg[i].grid_2D.name_grid,
		                               NAME_CPL_GRID_2D,
		                               cfg_model_set->models_cfg[i].grid_3D.name_grid,
		                               NAME_CPL_GRID_3D);
}

void remap_model_set_cpl2model(data_bundle **bundles_src,
                                             data_bundle **bundles_dst,
                                             model_set_cfg *cfg_model_set,
                                             data_bundle_cfg *cfg_data_bundle)
{
	int i;
	int nfields_2d, nfields_3d;

	nfields_2d = cfg_data_bundle->num_fields_2D;
	nfields_3d = cfg_data_bundle->num_fields_3D;

	for (i = 0; i < cfg_model_set->num_models; i ++) 
		remap_data_bundle(bundles_src[i], 
		                               bundles_dst[i], 
		                               nfields_2d, 
		                               nfields_3d, 
		                               NAME_CPL_GRID_2D,
		                               cfg_model_set->models_cfg[i].grid_2D.name_grid,
		                               NAME_CPL_GRID_3D,
		                               cfg_model_set->models_cfg[i].grid_3D.name_grid);
}

void remap_atm_set_atm2cpl()
{
	remap_model_set_model2cpl(data_bundles_atm_output_model, 
	                                             data_bundles_atm_output_cpl,
	                                             &atm_model_set,
	                                             &cfg_atm_output_data_bundle);
}

void remap_atm_set_cpl2atm()
{
	remap_model_set_cpl2model(data_bundles_atm_input_model, 
	                                             data_bundles_atm_input_cpl,
	                                             &atm_model_set,
	                                             &cfg_atm_input_data_bundle);
}

void remap_lnd_set_lnd2cpl()
{
	remap_model_set_model2cpl(data_bundles_lnd_output_model, 
	                                             data_bundles_lnd_output_cpl,
	                                             &lnd_model_set,
	                                             &cfg_lnd_output_data_bundle);
}

void remap_lnd_set_cpl2lnd()
{
	remap_model_set_cpl2model(data_bundles_lnd_input_model, 
	                                             data_bundles_lnd_input_cpl,
	                                             &lnd_model_set,
	                                             &cfg_lnd_input_data_bundle);
}

void remap_ocn_set_ocn2cpl()
{
	remap_model_set_model2cpl(data_bundles_ocn_output_model, 
	                                             data_bundles_ocn_output_cpl,
	                                             &ocn_model_set,
	                                             &cfg_ocn_output_data_bundle);
}

void remap_ocn_set_cpl2ocn()
{
	remap_model_set_cpl2model(data_bundles_ocn_input_model, 
	                                             data_bundles_ocn_input_cpl,
	                                             &ocn_model_set,
	                                             &cfg_ocn_input_data_bundle);
}

void remap_sice_set_sice2cpl()
{
	remap_model_set_model2cpl(data_bundles_sice_output_model, 
	                                             data_bundles_sice_output_cpl,
	                                             &sice_model_set,
	                                             &cfg_sice_output_data_bundle);
}

void remap_sice_set_cpl2sice()
{
	remap_model_set_cpl2model(data_bundles_sice_input_model, 
	                                             data_bundles_sice_input_cpl,
	                                             &sice_model_set,
	                                             &cfg_sice_input_data_bundle);
}

void remap_model_set(char *name_set, char *direction)
{
	if (strcmp(name_set, NAME_ATM_SET) == 0)
		if (strcmp(direction, "cpl2model") == 0)
			remap_atm_set_cpl2atm();
		else if (strcmp(direction, "model2cpl") == 0) 
			remap_atm_set_atm2cpl();
		else printf("wrong remap direction %s\n", direction);
	else if (strcmp(name_set, NAME_OCN_SET) == 0)
		if (strcmp(direction, "cpl2model") == 0)
			remap_ocn_set_cpl2ocn();
		else if (strcmp(direction, "model2cpl") == 0) 
			remap_ocn_set_ocn2cpl();
		else printf("wrong remap direction %s\n", direction);
	else if (strcmp(name_set, NAME_LND_SET) == 0)
		if (strcmp(direction, "cpl2model") == 0)
			remap_lnd_set_cpl2lnd();
		else if (strcmp(direction, "model2cpl") == 0) 
			remap_lnd_set_lnd2cpl();
		else printf("wrong remap direction %s\n", direction);
	else if (strcmp(name_set, NAME_SICE_SET) == 0)
		if (strcmp(direction, "cpl2model") == 0)
			remap_sice_set_cpl2sice();
		else if (strcmp(direction, "model2cpl") == 0) 
			remap_sice_set_sice2cpl();
		else printf("wrong remap direction %s\n", direction);
	else printf("wrong model name %s\n", name_set);
}

