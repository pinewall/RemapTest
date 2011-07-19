#include "grid.h"
#include <stdio.h>
#include <string.h>
#include "io.h"
#include "models_cfg.h"
#include <netcdf.h>
#include <stdlib.h>
#include <math.h>


#define NC_ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}


#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define ZERO  0.0
#define TWO 2
#define HALF 0.5
#define PI2   TWO*PI
#define PIH   HALF*PI
#define MAX_VALUE    9999
#define MIN_VALUE   -9999
#define LATITUDE 1
#define LATLON    2

grid::grid(char *filename, char *grid_name, int model_set_id, int file_type)
{
	double last_lat;
	double last_lon;
	double *tmp_lon_array1;
	double *tmp_lon_array2;
	double *tmp_double;
	int i, j, k, m, n, count, indx;
	int *tmp_nlat_each_lon1;
	int *tmp_nlat_each_lon2;
	int *tmp_int;
	int *tmp_current_indx_each_lon;
	FILE *fp_grid;
	int rcode;
	int fid;
	int did;
	int vid;
	unsigned int value;
	int grid_dims[3];
	unsigned long start[2];
	unsigned long size[2];
	int *temp_mask;
	double max_lat;
	double min_lat;
	
	my_model_set_id = model_set_id;

	/*input is binary file*/
	strcpy(my_name, grid_name);
	if (file_type == BINARY_FILE) {
		if ((fp_grid = fopen(filename, "rb")) == NULL)
			printf("cannot open the file of grid\n");
		fseek(fp_grid, 0, SEEK_SET);
		fread(&ndim, sizeof(int), 1, fp_grid);
		fread(&npnts, sizeof(int), 1, fp_grid);
		fread(&nlons, sizeof(int), 1, fp_grid);
		fread(&nlats, sizeof(int), 1, fp_grid);			
		if (ndim == 3)
			fread(&nhghts, sizeof(int), 1, fp_grid);
		else nhghts = -1;
		fread(&nvertexes, sizeof(int), 1, fp_grid);
		center_lon_coords = new double [npnts];
		center_lat_coords = new double [npnts];
		fread(center_lat_coords, sizeof(double), npnts, fp_grid);
		fread(center_lon_coords, sizeof(double), npnts, fp_grid);
		if (ndim == 3) {
			center_heights = new double [npnts];
			fread(center_heights, sizeof(double), npnts, fp_grid);
		}
		else center_heights = NULL;
		vertex_lat_coords = new double [nvertexes * npnts];
		vertex_lon_coords = new double [nvertexes * npnts];
		fread(vertex_lat_coords, sizeof(double), nvertexes * npnts, fp_grid);
		fread(vertex_lon_coords, sizeof(double), nvertexes * npnts, fp_grid);
		if (ndim == 3) {
			vertex_heights = new double [nvertexes * npnts];
			fread(vertex_heights, sizeof(double), nvertexes * npnts, fp_grid);
		}
		else vertex_heights = NULL;
		area_in = new double [npnts];
		areas = new double [npnts];
//		fread(area_in, sizeof(double), npnts, fp_grid);
		mask = new bool *[num_model_set];
		model_frac = new double * [num_model_set];
		for (i = 0; i < num_model_set; i ++) {
			mask[i] = NULL;
			model_frac[i] = NULL;
		}
		mask[model_set_id] = new bool [npnts];
		model_frac[model_set_id] = new double [npnts];
		fread(mask[model_set_id], sizeof(bool), npnts, fp_grid);
	}
	if (file_type == NETCDF_FILE) {
		if((rcode = nc_open(filename, NC_NOWRITE, &fid))) 
			NC_ERR(rcode);
		if((rcode = nc_inq_dimid(fid, "grid_rank", &did))) 
			NC_ERR(rcode);
		if((rcode = nc_inq_dimlen(fid, did, &value))) 
			NC_ERR(rcode);
		ndim = value;
		if((rcode = nc_inq_dimid(fid, "grid_size", &did))) 
			NC_ERR(rcode);
		if((rcode = nc_inq_dimlen(fid, did, &value))) 
			NC_ERR(rcode);		
		npnts = value;
		if((rcode = nc_inq_dimid(fid, "grid_corners", &did))) 
			NC_ERR(rcode); 
		if((rcode = nc_inq_dimlen(fid, did, &value))) 
			NC_ERR(rcode);
		nvertexes = value;
		if((rcode = nc_inq_varid(fid, "grid_dims", &vid))) 
			NC_ERR(rcode);
		value = ndim;
		if((rcode = nc_get_var_int(fid, vid, grid_dims))) 
			NC_ERR(rcode);
		nlons = grid_dims[0];
		nlats = grid_dims[1];
		nhghts = grid_dims[2];
		if (ndim < 3)
			nhghts = -1;
		center_lon_coords = new double [npnts];
		center_lat_coords = new double [npnts];
		vertex_lat_coords = new double [nvertexes * npnts];
		vertex_lon_coords = new double [nvertexes * npnts];
		area_in = new double [npnts];
		areas = new double [npnts];
		mask = new bool *[num_model_set];
		model_frac = new double * [num_model_set];
		for (i = 0; i < num_model_set; i ++) {
			mask[i] = NULL;
			model_frac[i] = NULL;
		}
		mask[model_set_id] = new bool [npnts];
		model_frac[model_set_id] = new double [npnts];
		if((rcode = nc_inq_varid(fid, "grid_center_lat", &vid))) 
			NC_ERR(rcode);
		if((rcode = nc_get_var_double(fid, vid, center_lat_coords)))
			NC_ERR(rcode);
		if((rcode = nc_inq_varid(fid,"grid_center_lon",&vid))) 
			NC_ERR(rcode);
		if((rcode = nc_get_var_double(fid, vid, center_lon_coords))) 
			NC_ERR(rcode);
		if((rcode = nc_inq_varid(fid, "grid_imask", &vid))) 
			NC_ERR(rcode);   
		temp_mask = new int [npnts];
		if((rcode = nc_get_var_int(fid, vid, temp_mask))) 
			NC_ERR(rcode);  
		for (i = 0; i < npnts; i ++)
			mask[model_set_id][i] = (temp_mask[i] == 1);
		delete [] temp_mask;
		if((rcode = nc_inq_varid(fid, "grid_corner_lat", &vid))) 
			NC_ERR(rcode);
		if((rcode = nc_get_var_double(fid, vid, vertex_lat_coords)))
			NC_ERR(rcode);
		if((rcode = nc_inq_varid(fid, "grid_corner_lon", &vid))) 
			NC_ERR(rcode);
		if((rcode = nc_get_var_double(fid, vid, vertex_lon_coords)))
			NC_ERR(rcode);
		center_heights = NULL;
		vertex_heights = NULL;
	}
#ifdef DEBUG_GRID
	check_grid_info(0);
#endif

	/*calculate the maximum and minimum latitudes of each row*/
	max_lat_each_row = new double [nlats];
	min_lat_each_row = new double [nlats];
	nlon_each_lat = new int [nlats];
	nlat_each_lon = new int [nlons];
	different_lons = new double [nlons];
	different_lats = new double [nlats];
	reverse_indx = new int [npnts];

	if (npnts == nlats * nlons) {		// logical rectangular grid
		for (i = 0; i < nlats; i ++) {
			nlon_each_lat[i] = nlons;
			max_lat = -100000;
			min_lat = 100000;
			for (j = 0; j < nlons; j ++) {
				if (center_lat_coords[i * nlons + j] > max_lat)
					max_lat = center_lat_coords[i * nlons + j];
				if (center_lat_coords[i * nlons + j] < min_lat)
					min_lat = center_lat_coords[i * nlons + j];
			}
			max_lat_each_row[i] = max_lat;
			min_lat_each_row[i] = min_lat;
		}
		for (i = 0; i < nlats; i ++)
			if (max_lat_each_row[i] != min_lat_each_row[i])
				break;
		if (i != nlats) {
			type_grid = 1;
			printf("irregular grid\n");
		}
		else {
			type_grid = 0;
			printf("regular grid\n");
		}
	}
	else {
		/*calculate the number of longitudes in each latitude*/
		for (i = 0; i < nlats; i ++)
			nlon_each_lat[i] = 0;
		last_lat = center_lat_coords[0];
		indx = 0;
		for (i = 0, count = 0; i < npnts; i ++) {
			if (last_lat == center_lat_coords[i])
				count ++;
			else {
				nlon_each_lat[indx ++] = count; 
				last_lat = center_lat_coords[i];
				count = 1;
			}
		}
		nlon_each_lat[indx ++] = count; 
		
		/*calculate the different latitudes in the grid*/
		for (i = 0, count = 0; i < nlats; i ++) {
			different_lats[i] = center_lat_coords[count];
			count += nlon_each_lat[i];
		}

		/*calculate the number of latitudes in each longitude and the different longitudes in the grid*/
		tmp_lon_array1 = new double [nlons];
		tmp_lon_array2 = new double [nlons];
		tmp_nlat_each_lon1 = new int [nlons];
		tmp_nlat_each_lon2 = new int [nlons];
		for (i = 0; i < nlons; i ++)
			tmp_nlat_each_lon1[i] = 0;
		for (i = 0; i < nlon_each_lat[0]; i ++) {
			tmp_nlat_each_lon1[i] = 1;
			tmp_lon_array1[i] = center_lon_coords[i];
		}
		count = nlon_each_lat[0];
		for (i = 1, j = nlon_each_lat[0]; i < nlats; i ++) {
			k = 0;
			m = 0;
			n = 0;
			while (k < nlon_each_lat[i] && m < count) {
				if (tmp_lon_array1[m] == center_lon_coords[j + k]) {
					tmp_nlat_each_lon2[n] = tmp_nlat_each_lon1[m] + 1;
					tmp_lon_array2[n ++] = tmp_lon_array1[m];
					m ++;
					k ++;
				}
				else if (tmp_lon_array1[m] < center_lon_coords[j + k]) {
					tmp_nlat_each_lon2[n] = tmp_nlat_each_lon1[m];
					tmp_lon_array2[n ++] = tmp_lon_array1[m];
					m ++;
				}
				else {
					tmp_nlat_each_lon2[n] = 1;
					tmp_lon_array2[n ++] = center_lon_coords[j + k];
					k ++;
				}
			}
			for (; k < nlon_each_lat[i]; k ++) {
				tmp_nlat_each_lon2[n] = 1;
				tmp_lon_array2[n ++] = center_lon_coords[j + k];
			}
			for (; m < count; m ++) {
				tmp_nlat_each_lon2[n] = tmp_nlat_each_lon1[m];
				tmp_lon_array2[n ++] = tmp_lon_array1[m];
			}
			count = n;
			tmp_double = tmp_lon_array1;
			tmp_lon_array1 = tmp_lon_array2;
			tmp_lon_array2 = tmp_double;
			tmp_int = tmp_nlat_each_lon1;
			tmp_nlat_each_lon1 = tmp_nlat_each_lon2;
			tmp_nlat_each_lon2 = tmp_int;
			j += nlon_each_lat[i];
		}
		for (i = 0; i < nlons; i ++) {
			nlat_each_lon[i] = tmp_nlat_each_lon1[i];
			different_lons[i] = tmp_lon_array1[i];
		}

		/*check whether the grid is retangular*/
		is_rect = true;
		for (i = 0; i < nlats; i ++)
			if (nlons != nlon_each_lat[i])
				is_rect = false;
		for (i = 0; i < nlons; i ++)
			if (nlats != nlat_each_lon[i])
				is_rect = false;
			
		/*calculate the reverse index of each point when transforming the lat-major matrix into lon-major matrix*/
		tmp_current_indx_each_lon = new int [nlons];
		tmp_current_indx_each_lon[0] = 0;
		for (i = 1; i < nlons; i ++)
			tmp_current_indx_each_lon[i] = tmp_current_indx_each_lon[i - 1] + nlat_each_lon[i];
		for (i = 0, j = 0; i < nlats; i ++) {
			k = 0;
			m = 0;
			n = 0;
			while (k < nlon_each_lat[i] && m < nlons) {
				if (tmp_lon_array1[m] == center_lon_coords[j + k]) {
					reverse_indx[j + k] = tmp_current_indx_each_lon[m];
					tmp_current_indx_each_lon[m] ++;
					m ++;
					k ++;
				}
				else if (tmp_lon_array1[m] < center_lon_coords[j + k]) 
					m ++;
				else {
					printf("never happen case when building reverse indexes for grid\n");
				}
			}
			j += nlon_each_lat[i];
		}	

		delete [] tmp_current_indx_each_lon;
		delete [] tmp_lon_array1;
		delete [] tmp_lon_array2;
		delete [] tmp_nlat_each_lon1;
		delete [] tmp_nlat_each_lon2;
	}

#ifdef DEBUG_GRID
	check_grid_info(1);
#endif
}

void grid::test_bug()
{
	delete [] areas;
}
grid::~grid()
{
	delete [] center_lon_coords;
	delete [] center_lat_coords;
	if (center_heights != NULL)
		delete [] center_heights;
	delete [] vertex_lat_coords;
	delete [] vertex_lon_coords;
	if (vertex_heights != NULL)
		delete [] vertex_heights;
	delete [] area_in;
	delete [] areas;
	for (int i = 0; i < num_model_set; i ++) 
		if (mask[i] != NULL) {
			delete [] mask[i];
			delete [] model_frac[i];
		}
	delete [] mask;
	delete [] model_frac;
	delete [] nlat_each_lon;
	delete [] nlon_each_lat;
	delete [] reverse_indx;
	delete [] different_lats;
	delete [] max_lat_each_row;
	delete [] min_lat_each_row;
	delete [] different_lons;
	
	//adding
  delete [] bound_box; 						
	delete [] dims;	                                  
  delete [] bin_addr;                             
  delete [] bin_lats;           
  delete [] bin_lons;
	
}

#ifdef DEBUG_GRID
void grid::check_grid_info(int check_type)
{
	int i, j;
	double *tmp_lat_coords;
	double *tmp_lon_coords;

	/*check whether the grid is lat-major and sorted after loading it*/
	if (check_type == 0) {
		for (i = 1; i < npnts; i ++) {
			if (center_lat_coords[i] < center_lat_coords[i - 1]) {
				printf("grid is not lat-major and sorted: case 1\n");
				break;
			}
			if (center_lat_coords[i] == center_lat_coords[i - 1] &&
			    center_lon_coords[i] < center_lon_coords[i - 1]) {
				printf("grid is not lat-major and sorted: case 2 %d\n", i);
				break;
			}
		}
	}

	/*check the correctness of nlats, nlons and reverse matrix*/
	if (check_type == 1) {
		if (!is_rect)
			printf("%s is not a rect grid\n", my_name);
		tmp_lat_coords = new double [npnts];
		tmp_lon_coords = new double [npnts];
		for (i = 0; i < nlats; i ++)
			if (nlon_each_lat[i] == 0) {
				printf("nlats of grid is wrong\n");
				break;
			}
		for (i = 0; i < nlons; i ++)
			if (nlat_each_lon[i] == 0) {
				printf("nlons of grid is wrong\n");
				break;
			}
		for (i = 0; i < npnts; i ++) {
			tmp_lat_coords[reverse_indx[i]] = tmp_lat_coords[i];
			tmp_lon_coords[reverse_indx[i]] = tmp_lon_coords[i];
		}
		/*check whether the reversed grid is lon-major and sorted*/
		for (i = 1; i < npnts; i ++) {
			if (tmp_lon_coords[i] < tmp_lon_coords[i - 1]) {
				printf("the reversed grid is not lon-major and sorted: case 1\n");
				break;
			}
			if (tmp_lon_coords[i] == tmp_lon_coords[i - 1] &&
			    tmp_lat_coords[i] < tmp_lat_coords[i - 1]) {
				printf("the reversed grid is not lon-major and sorted: case 2\n");
				break;
			}
		}
		delete [] tmp_lat_coords;
		delete [] tmp_lon_coords;
	}
}
#endif


