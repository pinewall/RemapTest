#include "dist_remap.h"
#include "grid.h"
#include "t_sort.h"
#include "grid_func.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "t_get_nbr.h"

dist_remap_2D_operator::dist_remap_2D_operator(char *grid1_name, 
                                                                      char *grid2_name, 
                                                                      char *alg_name,
                                                                      int n_n_pnts,
                                                                      double dist_t,
                                                                      double power_)
: Common_remap(grid1_name, grid2_name, alg_name)
{
	num_nearest_pnts = n_n_pnts;
	if (n_n_pnts <= 0)
		num_nearest_pnts = 4;
	dist_threshold = dist_t;
	power_p = power_;
	num_wgts = npts_dst * num_nearest_pnts;
	wgt_indx_src = new int [num_wgts];
	wgt_indx_dst = NULL;
	wgt_values = new double [num_wgts];
	count = 0;

	if (dist_threshold < 10.0)
		dist_threshold = 10.0;
	if (power_p < 1.0)
		power_p = 1.0;
}

dist_remap_2D_operator::~dist_remap_2D_operator()
{
	delete [] wgt_indx_src;
	delete [] wgt_values;
}


void dist_remap_2D_operator::search_neighbors_all_pnts_dst(double* lon_coords_src, 
                                                                                                    double* lat_coords_src,
                                                                                                    double* lon_coords_dst, 
                                                                                                    double* lat_coords_dst, 
                                                                                                    int* nlon_each_lat_src)
{
	int point_indx;
	double orig_dist_thrsh = dist_threshold;
	double sum_time = 0.0;
	bool flag;
	int pos_indx[18];
	double local_weights[18];
	int i, j;
	bool *mask_src, *mask_dst;
	
	mask_src = cpl_grids->get_grid(grid_name_src)->get_mymask();
	mask_dst = cpl_grids->get_grid(grid_name_dst)->get_mymask();


	for (point_indx = 0; point_indx < npts_dst; point_indx ++) 
		if (mask_dst[point_indx]) {
			get_nearest_nbrs(lat_coords_dst[point_indx], 
		                           lon_coords_dst[point_indx], 
		                           cpl_grids->get_grid(grid_name_src),
		                           dist_threshold, 
		                           num_nearest_pnts, 
		                           wgt_indx_src + num_nearest_pnts * point_indx,
		                           wgt_values + num_nearest_pnts * point_indx, 
		                           power_p);
		}
		else {
			for (i = num_nearest_pnts * point_indx; i < num_nearest_pnts * (point_indx + 1); i ++) {
				wgt_indx_src[i] = 0;
				wgt_values[i] = 0;
			}
		}
}

void dist_remap_2D_operator::init_remap() 
{
	int nlats_src;
	grid *grid_src, *grid_dst;
	double* lat_coords_src;
	double *lon_coords_src;
	double* lat_coords_dst;
	double *lon_coords_dst;

	grid_src = cpl_grids->get_grid(grid_name_src);
	grid_dst = cpl_grids->get_grid(grid_name_dst);
	nlats_src = grid_src->get_num_lats();
	lat_coords_src = grid_src->get_lat_coords();
	lon_coords_src = grid_src->get_lon_coords();
	lat_coords_dst = grid_dst->get_lat_coords();
	lon_coords_dst = grid_dst->get_lon_coords();
	
	
	search_neighbors_all_pnts_dst(lon_coords_src, 
	                                          lat_coords_src, 
	                                          lon_coords_dst, 
	                                          lat_coords_dst, 
	                                          grid_src->get_nlon_each_lat());
}

void dist_remap_2D_operator::cal_remap(double *data_src, 
                                                              double *data_dst)
{
	int i, j, lb,ub = 0;
	double value;
	int *local_point_indx_src = wgt_indx_src;
	double *local_weight = wgt_values;

	//#pragma omp parallel for private (i, j, lb, ub, value) schedule(dynamic, 1000)
	for (i = 0; i < npts_dst; i ++) {
		value = 0;
		lb = i * num_nearest_pnts;
		ub = lb + num_nearest_pnts;
		for (j = lb; j < ub; j ++) 
			value += data_src[local_point_indx_src[j]] * local_weight[j];
		data_dst[i] = value;
	}
}

