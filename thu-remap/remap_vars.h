#define REMAP_VARS_H

class remap_vars 
{
	public:
		int  max_links_map1;    // current size of link arrays
              int  num_links_map1;    // actual number of links for remapping
              int  max_links_map2;    // current size of link arrays
              int  num_links_map2;    // actual number of links for remapping
              int  num_maps;              // num of remappings for this grid pair

		int  fld_choice;               // choice of field to be interpolated
              int  norm_opt;                // option for normalization (conserv only)
              int  resize_increment;   // default amount to increase array size
       
              int *grid1_add_map1;   // grid1 address for each link in mapping 1
              int *grid2_add_map1;   // grid2 address for each link in mapping 1
              int *grid1_add_map2;   // grid1 address for each link in mapping 2
              int *grid2_add_map2;   // grid2 address for each link in mapping 2

              double *wts_map1;        // map weights for each link (num_wts,max_links)
              double *wts_map2;        // map weights for each link (num_wts,max_links)
	
	   	int  num_wts;                 // num of weights used in remapping     

	       bool luse_grid_centers;// use centers for bounding boxes
              bool luse_grid1_area;   // use area from grid file
              bool luse_grid2_area;   // use area from grid file

	       
			  
	public:
		remap_vars(int, int);
		~remap_vars();

              // int get_num_wts() {return num_wts; }

       	void reset_remap_vars(int, int);		
 
};



