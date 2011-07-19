#ifndef DATA_CFG_H

#include "name_cfg.h"

#define DATA_CFG_H

#define MAX_NUM_FIELDS	128

#define NAME_ATM_OUTPUT_FIELDS	"atm_output_data_bundle"
#define NAME_ATM_INPUT_FIELDS		"atm_input_data_bundle"
#define NAME_LND_OUTPUT_FIELDS	"lnd_output_data_bundle"
#define NAME_ROFF_OUTPUT_FIELDS	"roff_output_data_bundle"
#define NAME_LND_INPUT_FIELDS		"lnd_input_data_bundle"
#define NAME_OCN_OUTPUT_FIELDS	"ocn_output_data_bundle"
#define NAME_OCN_INPUT_FIELDS		"ocn_input_data_bundle"
#define NAME_SICE_OUTPUT_FIELDS	"sice_output_data_bundle"
#define NAME_SICE_INPUT_FIELDS		"sice_input_data_bundle"
#define NAME_ATM_OCN_FIELDS		"atm_ocn_data_bundle"



struct data_field_cfg 
{
	char name_field[NAME_STR_SIZE];
	int type_field;						// The data type of the field: 0 -> integer, 1 -> single precision, 2 -> double precision
};

struct data_bundle_cfg
{
	char name_fields[NAME_STR_SIZE];
	int num_fields_2D;
	int num_fields_3D;
	data_field_cfg all_fields[MAX_NUM_FIELDS];				// All fields of 3D are after all fields of 2D
};

extern data_bundle_cfg cfg_atm_output_data_bundle;
extern data_bundle_cfg cfg_atm_input_data_bundle;
extern data_bundle_cfg cfg_lnd_output_data_bundle;
extern data_bundle_cfg cfg_roff_output_data_bundle;
extern data_bundle_cfg cfg_lnd_input_data_bundle;
extern data_bundle_cfg cfg_ocn_output_data_bundle;
extern data_bundle_cfg cfg_ocn_input_data_bundle;
extern data_bundle_cfg cfg_sice_output_data_bundle;
extern data_bundle_cfg cfg_sice_input_data_bundle;
extern data_bundle_cfg cfg_atm_ocn_data_bundle;

#endif