#include "remap_base.h"
#include "grid_func.h"
#include <string.h>
#include <stdio.h>


Common_remap::Common_remap(char *grid1_name,
                                                  char *grid2_name, 
                                                  char *alg_name)
{
    strcpy(grid_name_src, grid1_name);
    strcpy(grid_name_dst, grid2_name);
    strcpy(algorithm_name, alg_name);
    next = NULL;
    localize_grid_fields();
}

Common_remap *Common_remap::search_remap_operator(char *grid1_name,
                                                                                       char *grid2_name, 
                                                                                       char *alg_name)
{
    if (strcmp(grid1_name, grid_name_src) == 0 &&
        strcmp(grid2_name, grid_name_dst) == 0 &&
        strcmp(alg_name, algorithm_name) == 0)
        return this;

    if (next != NULL)
        return next->search_remap_operator(grid1_name, grid2_name, alg_name);

    return NULL;
}

void Common_remap::localize_grid_fields()
{
	grid *grid_src, *grid_dst;


    /* Localize the fields of source grid */
    grid_src = cpl_grids->get_grid(grid_name_src);
    npts_src = grid_src->get_num_points();
    num_vertexes_src = grid_src->get_num_vertexes();
    vertex_lats_src = grid_src->get_vertex_lat_coords();	   
    vertex_lons_src = grid_src->get_vertex_lon_coords();
    frac_src = grid_src->get_single_frac(0);
    mask_src = grid_src->get_single_mask(0);
    area_src = grid_src->get_cell_area();
    center_lats_src = grid_src->get_lat_coords();
    center_lons_src = grid_src->get_lon_coords();	
    num_lats_src = grid_src->get_num_lats();
    num_lons_src= grid_src->get_num_lons();

    /* Localize the fields of source grid */
    grid_dst = cpl_grids->get_grid(grid_name_dst);  
    npts_dst = grid_dst->get_num_points();
    num_vertexes_dst = grid_dst->get_num_vertexes();  
    vertex_lats_dst = grid_dst->get_vertex_lat_coords();
    vertex_lons_dst = grid_dst->get_vertex_lon_coords();
    frac_dst = grid_dst->get_single_frac(0);
    mask_dst = grid_dst->get_single_mask(0);
    area_dst = grid_dst->get_cell_area();
    center_lats_dst = grid_dst->get_lat_coords(); 
    center_lons_dst = grid_dst->get_lon_coords(); 
    num_lats_dst = grid_dst->get_num_lats();
    num_lons_dst = grid_dst->get_num_lons();
}
