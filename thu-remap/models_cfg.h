#ifndef MODELS_CFG_H
#define MODELS_CFG_H

#include "name_cfg.h"

#define MAX_NUM_MODELS	128
#define NAME_ATM_SET	"atm_models"
#define NAME_LND_SET	"lnd_models"
#define NAME_OCN_SET	"ocn_models"
#define NAME_SICE_SET	"sice_models"

struct grid_cfg
{
	char name_grid[NAME_STR_SIZE];
	char name_file[NAME_STR_SIZE];
	int type_file;							//0: binary file, 1: netcdf file
	int num_dimension;
};

struct model_cfg
{
	char name_model[NAME_STR_SIZE];
	grid_cfg grid_2D;
	grid_cfg grid_3D;
};

struct model_set_cfg
{
	char name_model_set[NAME_STR_SIZE];
	int num_models;
	model_cfg models_cfg[MAX_NUM_MODELS];
};


extern int num_model_set;
extern int atm_set_id;
extern int ocn_set_id;
extern int lnd_set_id;
extern int sice_set_id;
extern model_set_cfg atm_model_set;
extern model_set_cfg lnd_model_set;
extern model_set_cfg ocn_model_set;
extern model_set_cfg sice_model_set;

#endif
