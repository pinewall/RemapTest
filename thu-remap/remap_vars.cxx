#include "remap_vars.h"

#include <stdio.h>
#include <stdlib.h>


#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))

remap_vars::remap_vars(int grid1_size, int grid2_size) 
{
 
  int i;

  // !!!!!!!
  num_maps = 2;                 
  num_wts    = 3;                              //map_type_conserv
  norm_opt   = 3;                              //norm_opt_frcarea
  fld_choice  = 2;
  luse_grid1_area = false;
  luse_grid2_area = false;
  // !!!!!!!
 
	num_links_map1 = 0;
       max_links_map1 = 4*grid2_size;    // can be changed later

	if (num_maps > 1) {
		num_links_map2 = 0;
		max_links_map1 = MAX(4*grid1_size, 4*grid2_size);
		max_links_map2 = max_links_map1;
		}

	resize_increment = 0.1*MAX(grid1_size, grid2_size);

       //allocate address and weight arrays for mapping 1
       grid1_add_map1 = new int[max_links_map1];
       grid2_add_map1 = new int[max_links_map1];
       wts_map1 = new double[num_wts*max_links_map1];

       for (i=0;i<max_links_map1;i++){
	   	grid1_add_map1[i] = -1;
		grid2_add_map1[i] = -1;
	   	}

	for (i=0;i<num_wts*max_links_map1;i++)
	   	wts_map1[i] = -1;
	   
       //allocate address and weight arrays for mapping 2 if necessary 
       if (num_maps > 1){
	  	grid1_add_map2 = new int[max_links_map2];
              grid2_add_map2 = new int[max_links_map2];
              wts_map2 = new double[num_wts*max_links_map2];
      	}

	for (i=0;i<max_links_map2;i++){
	   	grid1_add_map2[i] = -1;
		grid2_add_map2[i]  = -1;
	   	}

	for (i=0;i<num_wts*max_links_map2;i++)
	   	wts_map2[i] = -1;
 
}

remap_vars::~remap_vars() 
{
 	delete [] grid1_add_map1;
	delete [] grid2_add_map1;
	delete [] wts_map1;
	delete [] grid1_add_map2;
	delete [] grid2_add_map2;
	delete [] wts_map2; 
}

void remap_vars::reset_remap_vars(int nmap, int increment) {
      int mxlinks,i;            // size of link arrays
      int *add1_tmp;        // temp array for resizing address arrays
      int *add2_tmp;        // temp array for resizing address arrays
      double *wts_tmp;    // temp array for resizing weight arrays

      if (nmap == 1) { 
	  	// resize map 1 arrays if required
              add1_tmp = grid1_add_map1;                 // ????????????????????????
              add2_tmp = grid2_add_map1;
              wts_tmp  = wts_map1;

	       // deallocate originals and increment max_links then reallocate arrays at new size
             mxlinks = max_links_map1;
             max_links_map1 = mxlinks + increment;

             grid1_add_map1 = new int[max_links_map1];
             grid2_add_map1 = new int[max_links_map1];
             wts_map1 = new double[num_wts*max_links_map1];

             //restore original values from temp arrays and deallocate temps
             mxlinks = MIN(mxlinks, max_links_map1);
             for (i=0;i<mxlinks;i++) {
			grid1_add_map1[i] = add1_tmp[i];
                        grid2_add_map1[i] = add2_tmp[i];
		}
             for (i=0;i<num_wts*mxlinks;i++) 
		        wts_map1[i] = wts_tmp[i];

             delete [] add1_tmp;
	     delete [] add2_tmp;
	     delete [] wts_tmp;
	}

	else { // nmap==2

	      // allocate temporaries to hold original values
              add1_tmp = grid1_add_map2;                 // ????????????????????????
              add2_tmp = grid2_add_map2;
              wts_tmp  = wts_map2;
          
              //deallocate originals and increment max_links then reallocate arrays at new size
             mxlinks = max_links_map2;
             max_links_map2 = mxlinks + increment;

             grid1_add_map2 = new int[max_links_map2];
             grid2_add_map2 = new int[max_links_map2];
             wts_map2 = new double[num_wts*max_links_map2];

             //restore original values from temp arrays and deallocate temps
             mxlinks = MIN(mxlinks, max_links_map2);
             for (i=0;i<mxlinks;i++) {
			grid1_add_map2[i] = add1_tmp[i];
                     grid2_add_map2[i] = add2_tmp[i];
             }
             for (i=0;i<num_wts*mxlinks;i++) 
			 wts_map2[i] = wts_tmp[i];

             delete [] add1_tmp;
	      delete [] add2_tmp;
	      delete [] wts_tmp;
	}

}
