/*
??????????
*/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "dist_remap.h"
#include "grid.h"
#include "grid_func.h"
#include "data_func.h"
#include "io.h"
#include "remap_func.h"
#include "remap_cfg.h"
#include <sys/time.h>
#include <string.h>
#include "models_cfg.h"

// modified by ssq
#include <netcdf.h>
#include <math.h>

void wtime(double *t)
{
  static int sec = -1;
  struct timeval tv;
  gettimeofday(&tv, (__timezone_ptr_t)0);
  if (sec < 0) 
   sec = tv.tv_sec;
  *t = (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
}

void Print_Dst_To_File(char* str, int npts_dst_para, double* lon_coords_dst_para, double* lat_coords_dst_para, double* data_dst_para)
{
	FILE* out;
	if ((out = fopen(str, "a+")) == NULL)
	{
		printf("Cannot open file!\n");
	}
	else
	{
		int count = 0;
		while (count < npts_dst_para)
		{
			char out1[20] , out2[20] , out3[20];
			const char* out4 = ",";
			const char* out5 = ":";
			const char* out6 = "\n";
			sprintf(out1 , "res-%d: %lf" , count, lon_coords_dst_para[count]);
			sprintf(out2 , "%lf" , lat_coords_dst_para[count]);
			sprintf(out3 , "%lf" , data_dst_para[count]);

			fputs((const char*)out1, out);
			fputs(out4 , out);

			fputs((const char*)out2 , out);	
			fputs(out5 , out); 

			fputs((const char*)out3 , out);
			fputs(out6 , out);

			count++;
		}
	}
	fclose(out);
}


int main (int argc, char **argv)
{
	int i, j;
	int read_size;
	double * data_src_array;				// the src grid_2D data to be remapped
	int size_input_data_file;
	double * data_dst_array;				// the dst grid_2D data after remapping				
    double * data_dst_remap;
    double * data_dst_error;

	FILE *fp_data_src;
	FILE *fp_data_dst;
	grid *grid_src, *grid_dst;
	char *name_data_field = data_remaps_cfg[0].name_data_field;
	char *grid_name1 = data_remaps_cfg[0].name_grid_src;
	char *grid_name2 = data_remaps_cfg[0].name_grid_dst;
	double time1, time2;

	if (argc < 5) {
		printf("usage <grid1> <grid2> <test_case:[1-3]> <algorithm>\n");
		exit(0);
	}

	strcpy(atm_model_set.models_cfg[0].grid_2D.name_file, argv[1]);
	strcpy(ocn_model_set.models_cfg[0].grid_2D.name_file, argv[2]);
	strcpy(data_remaps_cfg[0].name_remap_alg, argv[4]);
	generate_all_cpl_grids();
	grid_src = cpl_grids->get_grid(grid_name1);
	grid_dst = cpl_grids->get_grid(grid_name2);

	wtime(&time1);	
	generate_remap_operators();
	wtime(&time2);
	
	printf("initialize_time %f\n", time2 - time1);

	data_src_array = new double [grid_src->get_num_points()];
	data_dst_array = new double [grid_dst->get_num_points()];	
	data_dst_remap = new double [grid_dst->get_num_points()];	
	data_dst_error = new double [grid_dst->get_num_points()];	
    memset(data_src_array, 0, sizeof(double) * grid_src->get_num_points());
    memset(data_dst_array, 0, sizeof(double) * grid_dst->get_num_points());
    memset(data_dst_remap, 0, sizeof(double) * grid_dst->get_num_points());
    memset(data_dst_error, 0, sizeof(double) * grid_dst->get_num_points());

    if (strcmp(argv[3], "1") == 0)
    {
        double length = 0.1 * 2 * 3.14159265359;

        double grid_tmp;

        for (int i = 0; i < grid_src->get_num_points(); i++)
        {
            data_src_array[i] = cos(grid_src->get_lat_coords()[i]) 
                              * cos(grid_src->get_lon_coords()[i]);
            grid_tmp = acos(-data_src_array[i]) / length;
            if (grid_tmp <= 1.0)
            {
                data_src_array[i] = 2.0 + cos(PI * grid_tmp);
            }
            else
            {
                data_src_array[i] = 1.0;
            }

            if (! grid_src->get_mymask()[i])
            {
                data_src_array[i] = 0.0;
            }
        }

        for (int i = 0; i < grid_dst->get_num_points(); i++)
        {
            data_dst_array[i] = cos(grid_dst->get_lat_coords()[i]) 
                              * cos(grid_dst->get_lon_coords()[i]);
            grid_tmp = acos(-data_dst_array[i]) / length;
            if (grid_tmp <= 1.0)
            {
                data_dst_array[i] = 2.0 + cos(PI * grid_tmp);
            }
            else
            {
                data_dst_array[i] = 1.0;
            }

            if (! grid_dst->get_mymask()[i])
            {
                data_dst_array[i] = 0.0;
            }
        }
    }
    else if (strcmp(argv[3], "2") == 0)
    {
        for (int i = 0; i < grid_src->get_num_points(); i++)
        {
            if (grid_src->get_mymask()[i])
            {
                data_src_array[i] = 2.0 + cos(grid_src->get_lat_coords()[i])* cos(grid_src->get_lat_coords()[i]) * cos(2.0 * grid_src->get_lon_coords()[i]);
            }
            else
            {
                data_src_array[i] = 0.0;
            }
        }

        for (int i = 0; i < grid_dst->get_num_points(); i++)
        {
            if (grid_dst->get_mymask()[i])
            {
                data_dst_array[i] = 2.0 + cos(grid_dst->get_lat_coords()[i])* cos(grid_dst->get_lat_coords()[i]) * cos(2.0 * grid_dst->get_lon_coords()[i]);
            }
            else
            {
                data_dst_array[i] = 0.0;
            }
        }
    }
    else if (strcmp(argv[3], "3") == 0)
    {
        for (int i = 0; i < grid_src->get_num_points(); i++)
        {
            if (grid_src->get_mymask()[i])
            {
                data_src_array[i] = 2.0 + pow(cos(grid_src->get_lat_coords()[i]),16) * cos(16 * grid_src->get_lon_coords()[i]);
            }
            else
            {
                data_src_array[i] = 0;
            }
        }

        for (int i = 0; i < grid_dst->get_num_points(); i++)
        {
            if (grid_dst->get_mymask()[i])
            {
                data_dst_array[i] = 2.0 + pow(cos(grid_dst->get_lat_coords()[i]),16) * cos(16 * grid_dst->get_lon_coords()[i]);
            }
            else
            {
                data_dst_array[i] = 0;
            }
        }
    }
    else
    {
        printf("Unknown Test Case\n");
        printf("\t[1]: cosine hill at lon=pi,lat=0\n");
        printf("\t[2]: pseudo-spherical harmonic l=2, m=2\n");
        printf("\t[3]: pseudo-spherical harmonic l=32, m=16\n");
        return -1;
    }

	for (i = 0; i < grid_src->get_num_points(); i += grid_src->get_num_points()) {
		wtime(&time1);
		do_remap(name_data_field, grid_name1, grid_name2, data_src_array + i, data_dst_remap);
		//do_remap(name_data_field, grid_name2, grid_name1, data_dst, data_src + i / sizeof(double));
		wtime(&time2);

		printf("remap_time %f\n", time2 - time1);
	
		//Print_Dst_To_File("data_dst1", grid_dst->get_num_points(), grid_dst->get_lon_coords(), grid_dst->get_lat_coords(), data_dst);
	}

    // modified by ssq
    int ncid;
    int status;

    // variable ID
    int nc_srcgrdcntrlat_id, nc_srcgrdcntrlon_id;
    int nc_dstgrdcntrlat_id, nc_dstgrdcntrlon_id;
    int nc_srcgrdimask_id, nc_dstgrdimask_id;
    int nc_srcgrdarea_id, nc_dstgrdarea_id;
    int nc_srcgrdfrac_id, nc_dstgrdfrac_id;
    int nc_srcarray_id, nc_dstarray_id;
    int nc_srcgradlat_id, nc_srcgradlon_id;
    int nc_dstremap1_id, nc_dstremap2_id;
    int nc_dsterror1_id, nc_dsterror2_id;

    status = nc_create ("THUoutput.nc", NC_CLOBBER, &ncid);

    char remap_name[1024];
    strcpy(remap_name, grid_src->get_grid_name());
    strcat(remap_name, "  -->  ");
    strcat(remap_name, grid_dst->get_grid_name());
    status = nc_put_att_text (ncid, NC_GLOBAL, "title", strlen(remap_name), remap_name);

    // define grid size dimensions
    int nc_grid1size_id[2];
    int nc_grid2size_id[2];
    status = nc_def_dim (ncid, "grid1_dim1", grid_src->get_num_lats(), &(nc_grid1size_id[0]));
    status = nc_def_dim (ncid, "grid1_dim2", grid_src->get_num_lons(), &(nc_grid1size_id[1]));
    status = nc_def_dim (ncid, "grid2_dim1", grid_dst->get_num_lats(), &(nc_grid2size_id[0]));
    status = nc_def_dim (ncid, "grid2_dim2", grid_dst->get_num_lons(), &(nc_grid2size_id[1]));

    // define grid center latitude array
    status = nc_def_var (ncid, "src_grid_center_lat", NC_DOUBLE, 2, nc_grid1size_id, &nc_srcgrdcntrlat_id);
    status = nc_put_att_text (ncid, nc_srcgrdcntrlat_id, "units", 7, "radians");
    status = nc_def_var (ncid, "dst_grid_center_lat", NC_DOUBLE, 2, nc_grid2size_id, &nc_dstgrdcntrlat_id);
    status = nc_put_att_text (ncid, nc_dstgrdcntrlat_id, "units", 7, "radians");

    // define grid center longitude array
    status = nc_def_var (ncid, "src_grid_center_lon", NC_DOUBLE, 2, nc_grid1size_id, &nc_srcgrdcntrlon_id);
    status = nc_put_att_text (ncid, nc_srcgrdcntrlon_id, "units", 7, "radians");
    status = nc_def_var (ncid, "dst_grid_center_lon", NC_DOUBLE, 2, nc_grid2size_id, &nc_dstgrdcntrlon_id);
    status = nc_put_att_text (ncid, nc_dstgrdcntrlon_id, "units", 7, "radians");

    // define grid mask
    status = nc_def_var (ncid, "src_grid_imask", NC_INT, 2, nc_grid1size_id, &nc_srcgrdimask_id);
    status = nc_put_att_text (ncid, nc_srcgrdimask_id, "units", 8, "unitless");
    status = nc_def_var (ncid, "dst_grid_imask", NC_INT, 2, nc_grid2size_id, &nc_dstgrdimask_id);
    status = nc_put_att_text (ncid, nc_dstgrdimask_id, "units", 8, "unitless");
    
    // define grid area
    status = nc_def_var (ncid, "src_grid_area", NC_DOUBLE, 2, nc_grid1size_id, &nc_srcgrdarea_id);
    status = nc_def_var (ncid, "dst_grid_area", NC_DOUBLE, 2, nc_grid2size_id, &nc_dstgrdarea_id);

    // define grid frac
    status = nc_def_var (ncid, "src_grid_frac", NC_DOUBLE, 2, nc_grid1size_id, &nc_srcgrdfrac_id);
    status = nc_def_var (ncid, "dst_grid_frac", NC_DOUBLE, 2, nc_grid2size_id, &nc_dstgrdfrac_id);

    // define source arrays
    status = nc_def_var (ncid, "src_array", NC_DOUBLE, 2, nc_grid1size_id, &nc_srcarray_id);
    status = nc_def_var (ncid, "src_grad_lat", NC_DOUBLE, 2, nc_grid1size_id, &nc_srcgradlat_id);
    status = nc_def_var (ncid, "src_grad_lon", NC_DOUBLE, 2, nc_grid1size_id, &nc_srcgradlon_id);

    // define destination arrays
    status = nc_def_var (ncid, "dst_array", NC_DOUBLE, 2, nc_grid2size_id, &nc_dstarray_id);
    status = nc_def_var (ncid, "dst_remap1", NC_DOUBLE, 2, nc_grid2size_id, &nc_dstremap1_id);
    status = nc_def_var (ncid, "dst_remap2", NC_DOUBLE, 2, nc_grid2size_id, &nc_dstremap2_id);

    // define error arrays
    status = nc_def_var (ncid, "dst_error1", NC_DOUBLE, 2, nc_grid2size_id, &nc_dsterror1_id);
    status = nc_def_var (ncid, "dst_error2", NC_DOUBLE, 2, nc_grid2size_id, &nc_dsterror2_id);

    // end definination stage
    status = nc_enddef(ncid);


    /* write some grid info */

    // write grid center latitude array
    status = nc_put_var_double (ncid, nc_srcgrdcntrlat_id, grid_src->get_lat_coords());
    status = nc_put_var_double (ncid, nc_dstgrdcntrlat_id, grid_dst->get_lat_coords());

    // write grid center longitude array
    status = nc_put_var_double (ncid, nc_srcgrdcntrlon_id, grid_src->get_lon_coords());
    status = nc_put_var_double (ncid, nc_dstgrdcntrlon_id, grid_src->get_lon_coords());

    // write grid mask
    int * grid1_imask = new int [grid_src->get_num_points()];
    int * grid2_imask = new int [grid_dst->get_num_points()];
    for (int i = 0; i < grid_src->get_num_points(); i++)
        grid1_imask[i] = (grid_src->get_mymask()[i]) ? 1 : 0;
    for (int i = 0; i < grid_dst->get_num_points(); i++)
        grid2_imask[i] = (grid_dst->get_mymask()[i]) ? 1 : 0;
    status = nc_put_var_int (ncid, nc_srcgrdimask_id, grid1_imask);
    status = nc_put_var_int (ncid, nc_dstgrdimask_id, grid2_imask);
    delete [] grid1_imask;
    delete [] grid2_imask;

    // write grid area
    status = nc_put_var_double (ncid, nc_srcgrdarea_id, grid_src->get_cell_area());
    status = nc_put_var_double (ncid, nc_dstgrdarea_id, grid_dst->get_cell_area());

    // write grid frac [added later]

    // write source and destination grids
    status = nc_put_var_double (ncid, nc_srcarray_id, data_src_array);
    status = nc_put_var_double (ncid, nc_dstarray_id, data_dst_array);
    status = nc_put_var_double (ncid, nc_dstremap1_id, data_dst_remap);
    status = nc_put_var_double (ncid, nc_dstremap2_id, data_dst_remap);

    // write errors
    for (int i = 0; i < grid_dst->get_num_points(); i++)
    {
        data_dst_error[i] = data_dst_remap[i] - data_dst_array[i];
    }
    status = nc_put_var_double (ncid, nc_dsterror1_id, data_dst_error);
    status = nc_put_var_double (ncid, nc_dsterror2_id, data_dst_error);

    // close file
    status = nc_close(ncid);

	finalize_all_cpl_grids();
	finalize_remap_operators();
	return 0;
}
