#include "remap_cfg.h"

int num_data_remaps = 1;




remap_cfg data_remaps_cfg[] = {{"Sa_z", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 4, 0, 120.0, 1.0}};
/*
                                                {"Sa_u", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0}, 
                                                {"Sa_v", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sa_tbot", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sa_ptem", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sa_shum", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sa_dens", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sa_pbot", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sa_pslv", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxa_lwdn", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxa_rainc", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxa_rainl", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxa_snowc", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxa_snowl", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxa_swndr", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxa_swvdr", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxa_swndf", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxa_swvdf", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxa_swnet", "rect_grid_atm_360x180", "rect_grid_ocn_360x196", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sx_tref", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sx_qref", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sx_avsdr", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sx_anidr", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sx_avsdf", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sx_anidf", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sx_t", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"So_t", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sx_snowh", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sx_ifrac", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Sx_ofrac", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxx_taux", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxx_tauy", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxx_lat", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxx_sen", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxx_lwup", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},
                                                {"Faxx_evap", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0},                                              
                                                {"Sa_z", "rect_grid_ocn_360x196", "rect_grid_atm_360x180", "dist_remap_2D", 6, 0, 120.0, 8.0}};

*/
