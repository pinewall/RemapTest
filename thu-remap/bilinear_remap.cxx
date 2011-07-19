#include "bilinear_remap.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "grid.h"
#include "grid_func.h"
#include "t_get_nbr.h"

bilinear_remap_2D_operator::bilinear_remap_2D_operator(char *grid1_name, 
												   char *grid2_name, 
												   char *alg_name)
: Common_remap(grid1_name, grid2_name, alg_name)
{
	grid *grid_src, *grid_dst;
	int i,j;

	grid_src = cpl_grids->get_grid(grid_name_src);
	grid_dst = cpl_grids->get_grid(grid_name_dst);

	nlon_each_lat_src=grid_src->get_nlon_each_lat();
	nlon_each_lat_dst=grid_dst->get_nlon_each_lat();

	lat_begindx_src = new int [num_lats_src];
	lat_begindx_dst = new int [num_lats_dst];
	num_wgts = npts_dst * 4;
	wgt_indx_src = new int [num_wgts];
	wgt_indx_dst = NULL;
	wgt_values = new double [num_wgts];
		
	lat_begindx_src[0] = 0;
	for(i = 1; i < num_lats_src; i ++) 
		lat_begindx_src[i] = lat_begindx_src[i - 1] + nlon_each_lat_src[i - 1];
		
	lat_begindx_dst[0] = 0;
	for(i = 1; i < num_lats_dst; i ++)
		lat_begindx_dst[i] = lat_begindx_dst[i - 1] + nlon_each_lat_dst[i - 1];
}

bilinear_remap_2D_operator::~bilinear_remap_2D_operator()
{
	delete [] lat_begindx_src;
	delete [] lat_begindx_dst;
	delete [] wgt_indx_src;
	delete [] wgt_values;
}


void bilinear_remap_2D_operator::wgt_bilinear_regular( double *lons_src,
                                                                         double *lats_src,
                                                                         double lon_dst,
                                                                         double lat_dst,
                                                                         double *wgts) 
{
	double t1,t2;
	double u;

	if (!((lon_dst >= lons_src[0]) &&
		(lons_src[1] >= lons_src[0]))) {
		if (lon_dst < lons_src[0])
			lon_dst += 360;
		if (lons_src[1] < lons_src[0])
			lons_src[1] += 360;
	}
	if (!((lon_dst >= lons_src[2]) &&
		(lons_src[3] >= lons_src[2]))) {
		if (lon_dst < lons_src[2])
			lon_dst += 360;
		if (lons_src[3] < lons_src[2])
			lons_src[3] += 360;
	}
	t1=(lon_dst-lons_src[0])/(lons_src[1]-lons_src[0]);
	t2=(lon_dst-lons_src[2])/(lons_src[3]-lons_src[2]);
	u=(lat_dst-lats_src[0])/(lats_src[1]-lats_src[0]);
	wgts[0] = (1 - u) * (1 - t1);
	wgts[1] = (1 - u) * t1;
	wgts[2] = u * (1 - t2);
	wgts[3] = u * t2;
}			 



bool bilinear_remap_2D_operator::wgt_bilinear_irregular( double *lons_src,
                                                                           double *lats_src,
                                                                           double lon_dst,
                                                                           double lat_dst,
                                                                           double *wgts) 
{
	double dth1, dth2, dth3, dthp;
	double dph1, dph2, dph3, dphp;
	double iguess, jguess;
	double mat1, mat2, mat3, mat4;
	int iter, max_iter = 100;
	double determinant, deli, delj;
	double converge = 0.0000000001;

	dth1 = lats_src[1] - lats_src[0];
	dth2 = lats_src[3] - lats_src[0];
	dth3 = lats_src[2] - lats_src[1] - dth2;
	dph1 = lons_src[1] - lons_src[0];
	dph2 = lons_src[3] - lons_src[0];
	dph3 = lons_src[2] - lons_src[1];
	if (dph1 >  540) 
		dph1 = dph1 - 360;
	if (dph2 >  540) 
		dph2 = dph2 - 360;
	if (dph3 >  540) 
		dph3 = dph3 - 360;
	if (dph1 < -540) 
		dph1 = dph1 + 360;
	if (dph2 < -540) 
		dph2 = dph2 + 360;
	if (dph3 < -540) 
		dph3 = dph3 + 360;
	dph3 = dph3 - dph2;
	iguess = 0.5;
	jguess = 0.5;

	for (iter = 0; iter < max_iter; iter ++) {
		dthp = lat_dst - lats_src[0] - dth1*iguess - dth2*jguess - dth3*iguess*jguess;
		dphp = lon_dst - lons_src[0];
		if (dphp >  540) dphp = dphp - 360;
		if (dphp < -540) dphp = dphp + 360;
		dphp = dphp - dph1*iguess - dph2*jguess - dph3*iguess*jguess;
		mat1 = dth1 + dth3*jguess;
		mat2 = dth2 + dth3*iguess;
		mat3 = dph1 + dph3*jguess;
		mat4 = dph2 + dph3*iguess;
		
		determinant = mat1*mat4 - mat2*mat3;
		deli = (dthp*mat4 - mat2*dphp) / determinant;
		delj = (mat1*dphp - dthp*mat3) / determinant;
            if (fabs(deli) < converge && fabs(delj) < converge)
			break;

            iguess = iguess + deli;
            jguess = jguess + delj;
	}

	if (iter < max_iter) {
		wgts[0] = (1.0-iguess) * (1.0-jguess);
		wgts[1] = iguess * (1.0-jguess); 
		wgts[2] = iguess * jguess; 
		wgts[3] = (1.0-iguess) * jguess;
		return true;
	}
	else return false;
}			 


void bilinear_remap_2D_operator::init_remap()
{
	int i, j, k, point_dst = 0;
	int bound_lat_indx1, bound_lat_indx2;
	int bound_lon_indx1, bound_lon_indx2;
	int left_nbr, right_nbr;
	double left_lon, right_lon;
	double left_lat, right_lat;
	double bilinear_lons[4];
	double bilinear_lats[2];
	grid *grid_src, *grid_dst;
	bool *mask_src, *mask_dst;
	int type_grid_src, type_grid_dst;
	double *max_lat_each_row_src, *min_lat_each_row_src;
	double last_lat = -10000;
	bool find_bilinear_point;
	int nlons_src;
	bool lat_in, lon_in;
	double lons_src_box[4], lats_src_box[4];
	int corner_indx[4];
	double max_lat, min_lat;
	double max_lon, min_lon, tmp_lon;
	double *bound_boxes_src;
	double *cornner_boxes_lat_src, *cornner_boxes_lon_src;

	grid_src = cpl_grids->get_grid(grid_name_src);
	grid_dst = cpl_grids->get_grid(grid_name_dst);
	mask_src = grid_src->get_mymask();
	mask_dst = grid_dst->get_mymask();
	type_grid_src = grid_src->get_type_grid();
	type_grid_dst = grid_dst->get_type_grid();
	max_lat_each_row_src = grid_src->get_max_lat_each_row();
	min_lat_each_row_src = grid_src->get_min_lat_each_row();
	nlons_src = grid_src->get_num_lons();

	if (type_grid_src == 0) {
		for (point_dst = 0; point_dst < npts_dst; point_dst ++) {
			if (last_lat != center_lats_dst[point_dst]) {
				last_lat = center_lats_dst[point_dst];
				bound_lat_indx1 = -1;
				for(i = 0; i < num_lats_src; i ++)
		        		if(last_lat - max_lat_each_row_src[i] < 0) {
						bound_lat_indx1 = i - 1;
						break;
					}
				bound_lat_indx2 = bound_lat_indx1 + 1;
			}
			if (!mask_dst[point_dst]) {
			  	wgt_indx_src[point_dst * 4 + 0] = 0;
				wgt_indx_src[point_dst * 4 + 1] = 0;
				wgt_indx_src[point_dst * 4 + 2] = 0;
				wgt_indx_src[point_dst * 4 + 3] = 0;
				wgt_values[point_dst * 4 + 0] = 0;
				wgt_values[point_dst * 4 + 1] = 0;
				wgt_values[point_dst * 4 + 2] = 0;
				wgt_values[point_dst * 4 + 3] = 0;
				continue;
			}
			if (bound_lat_indx1 == -1) {
				get_nearest_nbrs(center_lats_dst[point_dst], 
			                           center_lons_dst[point_dst],
			                           grid_src, 
			                           120.0, 
			                           4, 
			                           wgt_indx_src + point_dst * 4,
			                           wgt_values + point_dst * 4, 
			                           1.0);
				continue;
			}
			bound_lon_indx2 = 0;
			bound_lon_indx1 = 0;
			for(i = 0; i < nlons_src; i ++) {
				left_nbr = (i - 1 + nlons_src) % nlons_src;
				right_nbr = i + lat_begindx_src[bound_lat_indx1];
				left_nbr += lat_begindx_src[bound_lat_indx1];
				left_lon = center_lons_src[left_nbr];
				right_lon = center_lons_src[right_nbr];
				if (right_lon >= left_lon) {
					if (center_lons_dst[point_dst] - right_lon < 0 &&
					    center_lons_dst[point_dst] - left_lon >= 0) {
						bound_lon_indx1 = i - 1;
						break;
					}
				}
				else if (center_lons_dst[point_dst] - right_lon < 0) {
					bound_lon_indx1 = i - 1;
					break;
				}
			}
			bound_lon_indx2 = i;
			if (i == 0 || i == nlons_src) {
				bound_lon_indx1 = nlons_src - 1;
				bound_lon_indx2 = 0;
			} 
			if (!mask_src[bound_lon_indx1 + lat_begindx_src[bound_lat_indx1]] ||
			    !mask_src[bound_lon_indx2 + lat_begindx_src[bound_lat_indx1]] ||
			    !mask_src[bound_lon_indx1 + lat_begindx_src[bound_lat_indx2]] ||
			    !mask_src[bound_lon_indx2 + lat_begindx_src[bound_lat_indx2]]) {
				get_nearest_nbrs(center_lats_dst[point_dst], 
			                           center_lons_dst[point_dst],
			                           grid_src, 
			                           120.0, 
			                           4, 
			                           wgt_indx_src + point_dst * 4,
			                           wgt_values + point_dst * 4, 
			                           1.0);
			}
			else {
			  	wgt_indx_src[point_dst * 4 + 0] = bound_lon_indx1 + lat_begindx_src[bound_lat_indx1];
				wgt_indx_src[point_dst * 4 + 1] = bound_lon_indx2 + lat_begindx_src[bound_lat_indx1];
				wgt_indx_src[point_dst * 4 + 2] = bound_lon_indx1 + lat_begindx_src[bound_lat_indx2];
				wgt_indx_src[point_dst * 4 + 3] = bound_lon_indx2 + lat_begindx_src[bound_lat_indx2];
				bilinear_lons[0]= center_lons_src[bound_lon_indx1 + lat_begindx_src[bound_lat_indx1]];
				bilinear_lons[1] = center_lons_src[bound_lon_indx2 + lat_begindx_src[bound_lat_indx1]];
				bilinear_lons[2] = center_lons_src[bound_lon_indx1 + lat_begindx_src[bound_lat_indx2]];
				bilinear_lons[3] = center_lons_src[bound_lon_indx2 + lat_begindx_src[bound_lat_indx2]];
				bilinear_lats[0] = center_lats_src[lat_begindx_src[bound_lat_indx1]]; 
				bilinear_lats[1] = center_lats_src[lat_begindx_src[bound_lat_indx2]];
				wgt_bilinear_regular(bilinear_lons, bilinear_lats, center_lons_dst[point_dst], center_lats_dst[point_dst], wgt_values + 4*point_dst);
			}
		}
	}
	else if (type_grid_src == 1) {
		bound_boxes_src = new double [4*npts_src];
		cornner_boxes_lat_src = new double [4*npts_src];
		cornner_boxes_lon_src = new double [4*npts_src];
		for (i = 0; i < npts_src - nlons_src; i ++) {
			j = i % nlons_src;
			k = i / nlons_src;
			corner_indx[0] = (j - 1 + nlons_src) % nlons_src + k * nlons_src;
			corner_indx[1] = i;
			corner_indx[2] = corner_indx[1] + nlons_src;
			corner_indx[3] = corner_indx[0] + nlons_src;
			for (k = 0; k < 4; k ++) {
				lons_src_box[k] = center_lons_src[corner_indx[k]];
				lats_src_box[k] = center_lats_src[corner_indx[k]];
			}

			max_lat = lats_src_box[0];
			min_lat = lats_src_box[0];
			min_lon = lons_src_box[0];
			max_lon = lons_src_box[0];
			for (k = 1; k < 4; k ++) {
				if (max_lat < lats_src_box[k])
					max_lat = lats_src_box[k];
				if (min_lat > lats_src_box[k])
					min_lat = lats_src_box[k];
				tmp_lon = lons_src_box[k];
				if (tmp_lon - lons_src_box[0] >= 180)
					tmp_lon -= 360;
				else if (tmp_lon - lons_src_box[0] <= -180)
					tmp_lon += 360;
				if (min_lon > tmp_lon)
					min_lon = tmp_lon;
				if (max_lon < tmp_lon)
					max_lon = tmp_lon;
			}
			bound_boxes_src[i*4] = min_lat;
			bound_boxes_src[i*4+1] = max_lat;
			bound_boxes_src[i*4+2] = min_lon;
			bound_boxes_src[i*4+3] = max_lon;
			tmp_lon = center_lons_dst[point_dst];			
			for (k = 0; k < 4; k ++) {
				cornner_boxes_lat_src[i*4+k] = lats_src_box[k];
				cornner_boxes_lon_src[i*4+k] = lons_src_box[k];
			}
		}
		for (point_dst = 0; point_dst < npts_dst; point_dst ++) {
			if (last_lat != center_lats_dst[point_dst]) {
				last_lat = center_lats_dst[point_dst];
				bound_lat_indx1 = -1;
				for(i = 0; i < num_lats_src; i ++)
		        		if(min_lat_each_row_src[i] > last_lat) {
						bound_lat_indx1 = i - 1;
						break;
					}
				if (i == num_lats_src) {
					for (i --; i >= 0 && max_lat_each_row_src[i] < last_lat; i --);
					bound_lat_indx1 = i;
				}
				bound_lat_indx2 = bound_lat_indx1 + 1;
				for (; bound_lat_indx1 >= 0; bound_lat_indx1 --) 
					if (max_lat_each_row_src[bound_lat_indx1] <= last_lat) {
						break;
					}
				if (bound_lat_indx1 == -1)
					bound_lat_indx1 = 0;
			}
			if (!mask_dst[point_dst]) {
			  	wgt_indx_src[point_dst * 4 + 0] = 0;
				wgt_indx_src[point_dst * 4 + 1] = 0;
				wgt_indx_src[point_dst * 4 + 2] = 0;
				wgt_indx_src[point_dst * 4 + 3] = 0;
				wgt_values[point_dst * 4 + 0] = 0;
				wgt_values[point_dst * 4 + 1] = 0;
				wgt_values[point_dst * 4 + 2] = 0;
				wgt_values[point_dst * 4 + 3] = 0;
				continue;
			}
			find_bilinear_point = false;
			if (bound_lat_indx2 == num_lats_src)
				bound_lat_indx2 --;
			for (j = bound_lat_indx1; j < bound_lat_indx2; j ++) {
				for(i = 0; i < nlons_src; i ++) {
					corner_indx[1] = i + j * nlons_src;
					min_lat = bound_boxes_src[corner_indx[1] * 4];
					max_lat = bound_boxes_src[corner_indx[1] * 4+1];
					min_lon = bound_boxes_src[corner_indx[1] * 4+2];
					max_lon = bound_boxes_src[corner_indx[1] * 4+3];
					tmp_lon = center_lons_dst[point_dst];
					if (tmp_lon - cornner_boxes_lon_src[corner_indx[1]*4] >= 180)
						tmp_lon -= 360;
					else if (tmp_lon - cornner_boxes_lon_src[corner_indx[1]*4] <= -180)
						tmp_lon += 360;
					if (min_lat <= center_lats_dst[point_dst] && 
					    max_lat >= center_lats_dst[point_dst] && 
					    tmp_lon <= max_lon && 
					    tmp_lon >= min_lon)	{
						for (k = 0; k < 4; k ++) {
							lats_src_box[k] = cornner_boxes_lat_src[corner_indx[1]*4+k];
							lons_src_box[k] = cornner_boxes_lon_src[corner_indx[1]*4+k];
						}
						if (point_in_polygon(center_lats_dst[point_dst], 
						                            center_lons_dst[point_dst], 
						                            4, 
						                            lats_src_box, 
						                            lons_src_box)) {
							find_bilinear_point = true;
							corner_indx[0] = (i - 1 + nlons_src) % nlons_src + j * nlons_src;
							corner_indx[2] = corner_indx[1] + nlons_src;
							corner_indx[3] = corner_indx[0] + nlons_src;
							
							if (mask_src[corner_indx[0]] &&
							    mask_src[corner_indx[1]] &&
							    mask_src[corner_indx[2]] &&
							    mask_src[corner_indx[3]] &&
							    wgt_bilinear_irregular(lons_src_box, 
							                                 lats_src_box, 
							                                 center_lons_dst[point_dst], 
							                                 center_lats_dst[point_dst], 
							                                 wgt_values + 4*point_dst)) {
								wgt_indx_src[point_dst*4 + 0] = corner_indx[0];
								wgt_indx_src[point_dst*4 + 1] = corner_indx[1];
								wgt_indx_src[point_dst*4 + 2] = corner_indx[2];
								wgt_indx_src[point_dst*4 + 3] = corner_indx[3];
							}
							else {
								get_nearest_nbrs(center_lats_dst[point_dst], 
								                           center_lons_dst[point_dst],
								                           grid_src, 
								                           120.0, 
								                           4, 
								                           wgt_indx_src + point_dst * 4,
								                           wgt_values + point_dst * 4, 
								                           1.0);
							}
							break;
						}
					}
				}
				if (find_bilinear_point)
					break;
			}
			if (!find_bilinear_point) {
				get_nearest_nbrs(center_lats_dst[point_dst], 
				                           center_lons_dst[point_dst],
				                           grid_src, 
				                           120.0, 
				                           4, 
				                           wgt_indx_src + point_dst * 4,
				                           wgt_values + point_dst * 4, 
				                           1.0);
			}
		}
		delete [] bound_boxes_src;
		delete [] cornner_boxes_lat_src;
		delete [] cornner_boxes_lon_src;
	}
}

void bilinear_remap_2D_operator::cal_remap(double *data_src, double *data_dst)
{

	int i, j, lb,ub = 0;
	int *local_point_indx_src = wgt_indx_src;
	double *local_wgts = wgt_values;

	//#pragma omp parallel for private (i, j) schedule(dynamic, 1000)
	for (i = 0, j = 0; i < npts_dst; i ++) {
		data_dst[i] = data_src[local_point_indx_src[j]] * local_wgts[j] +
		                   data_src[local_point_indx_src[j + 1]] * local_wgts[j + 1] +
		                   data_src[local_point_indx_src[j + 2]] * local_wgts[j + 2] +
		                   data_src[local_point_indx_src[j + 3]] * local_wgts[j + 3];
		j += 4;
	}
}

bool bilinear_remap_2D_operator::point_in_polygon(double lat_dst, 
                                                                     double lon_dst, 
                                                                     int ncorner, 
                                                                     double *lats_src, 
                                                                     double*lons_src)
{
	double lon_diff1, lon_diff2;
	double lat_diff1, lat_diff2;
	int i, next_i;
	double cross_product, cross_product_last;
	
	
	lon_diff1 = lons_src[0] - lon_dst;
	if (lon_diff1 > 180)
		lons_src[0] = lons_src[0] - 360;
	else if (lon_diff1 < -180) 
		lons_src[0] = lons_src[0] + 360;
	for (i = 1; i < ncorner; i ++) {
		lon_diff1 = lons_src[i] - lons_src[0];
		if (lon_diff1 >  180)
			lons_src[i] = lons_src[i] - 360;
		else if (lon_diff1 < -180)
			lons_src[i] = lons_src[i] + 360;
	}
	for (i = 0; i < ncorner; i ++) {
		next_i = (i + 1) % ncorner;
		lat_diff1 = lats_src[next_i] - lats_src[i];
		lon_diff1 = lons_src[next_i] - lons_src[i];
		lat_diff2 = lat_dst - lats_src[i];
		lon_diff2 = lon_dst - lons_src[i];
		if (lon_diff1 >  540)
			lon_diff1 = lon_diff1 - 360;
		else if (lon_diff1 < -540)
			lon_diff1 = lon_diff1 + 360;
		if (lon_diff2 > 540)
			lon_diff2 = lon_diff2 - 360;
		else if (lon_diff2 < -540)
			lon_diff2 = lon_diff2 + 360;
		cross_product = lon_diff1*lat_diff2 - lon_diff2*lat_diff1;
		if (i == 0) 
			cross_product_last = cross_product;
		if (cross_product*cross_product_last < 0) 
			return false;
		cross_product_last = cross_product;
	}
	return true;
}

