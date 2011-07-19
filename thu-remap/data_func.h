#ifndef DATA_FUNC_H
#define DATA_FUNC_H

#include "data_bundle.h"

extern data_bundle** data_bundles_atm_output_model;
extern data_bundle** data_bundles_atm_output_cpl;
extern data_bundle* data_bundles_atm_output_aggr;
extern data_bundle** data_bundles_atm_input_model;
extern data_bundle** data_bundles_atm_input_cpl;
extern data_bundle* data_bundles_atm_input_aggr;
extern data_bundle** data_bundles_lnd_output_model;
extern data_bundle** data_bundles_lnd_output_cpl;
extern data_bundle* data_bundles_lnd_output_aggr;
extern data_bundle** data_bundles_lnd_input_model;
extern data_bundle** data_bundles_lnd_input_cpl;
extern data_bundle* data_bundles_lnd_input_aggr;
extern data_bundle** data_bundles_ocn_output_model;
extern data_bundle** data_bundles_ocn_output_cpl;
extern data_bundle* data_bundles_ocn_output_aggr;
extern data_bundle** data_bundles_ocn_input_model;
extern data_bundle** data_bundles_ocn_input_cpl;
extern data_bundle* data_bundles_ocn_input_aggr;
extern data_bundle** data_bundles_sice_output_model;
extern data_bundle** data_bundles_sice_output_cpl;
extern data_bundle* data_bundles_sice_output_aggr;
extern data_bundle** data_bundles_sice_input_model;
extern data_bundle** data_bundles_sice_input_cpl;
extern data_bundle* data_bundles_sice_input_aggr;

extern void generate_all_data_bundles();
extern void finalize_all_data_bundles();
extern int copy_in_model_set_data(char*, char*);
extern int copy_out_model_set_data(char*, char*);

#endif
