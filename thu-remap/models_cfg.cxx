#include "models_cfg.h"
#include "io.h"

int num_model_set = 4;
int atm_set_id = 0;
int ocn_set_id = 1;
int lnd_set_id = 2;
int sice_set_id = 3;

model_set_cfg atm_model_set = {NAME_ATM_SET,
                                                2,
                                                {{"gamil360x180",
                                                {"rect_grid_atm_360x180", "grids/T42_Gaussian_Grid.nc", NETCDF_FILE, 2},
                                                {"rect_grid_atm_360x180", "T42", NETCDF_FILE, 3}},
                                                {"gamil360x180",
                                                {"rect_grid_atm_360x180", "T42", 0, 2},
                                                {"rect_grid_atm_360x180", "T42", 0, 3}}}};

model_set_cfg lnd_model_set = {NAME_LND_SET,
                                                0,
                                                {}};

model_set_cfg ocn_model_set = {NAME_OCN_SET,
                                                1,
                                                {{"ocn360x180",
                                                {"rect_grid_ocn_360x196", "grids/T62_Gaussian_Grid.nc", NETCDF_FILE, 2},
                                                {"rect_grid_ocn_360x196", "T62", 0, 3}}}};

model_set_cfg sice_model_set = {NAME_SICE_SET,
                                                0,
                                                {}};

//grid_cfg atm_grids_cfg[] = {{"rect_grid_atm_360x180", "file_grid_atm_1x1", 0}};
//grid_cfg ocn_grids_cfg[] = {{"rect_grid_ocn_360x196", "file_grid_atm_0.5x0.5", 0}};
//grid_cfg ocn_grids_cfg[] = {{"rect_grid_ocn_360x196", "file_grid_atm_0.5x0.5_sparse4_lat-45to45_lon90to270", 0}};

//grid_cfg ocn_grids_cfg[] = {{"rect_grid_ocn_360x196", "file_grid_atm_0.5x0.5", 0}};

//grid_cfg atm_grids_cfg[] = {{"rect_grid_atm_360x180", "file_grid_atm_1x1", 0}};
//grid_cfg ocn_grids_cfg[] = {{"rect_grid_ocn_360x196", "file_grid_ocn_1x1", 0}};
//grid_cfg lnd_grids_cfg[] = {{"rect_grid_atm_360x180", "file_grid_atm_1x1", 0}};
//grid_cfg seaice_grids_cfg[] = {{"rect_grid_ocn_360x196", "file_grid_ocn_1x1", 0}};

