#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>

#include "conserv_remap.h"
#include "grid.h"
#include "grid_func.h"
#include "remap_base.h"


#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))

#define ZERO  0.0
#define ONE   1.0
#define TWO   2.0
#define THREE 3.0
#define HALF  0.5
#define QUART 0.25
#define TINY  (1.e-14)              //?????????????????????????
#define PI2   TWO*PI
#define PIH   HALF*PI

#define NORM_OPT_NONE       1
#define NORM_OPT_DSTAREA 2
#define NORM_OPT_FRCAREA 3

#define north_thresh  2.00     // threshold for coord transf.
#define south_thresh -2.00    // threshold for coord transf.


#define default_size 1000

#define MAX_VALUE    9999
#define MIN_VALUE   -9999
#define LATITUDE 1
#define LATLON    2


void conserv_remap_2D_operator::cal_remap(double *data_src, double *data_dst)
{
     
       int tmp_nw = rvars->num_wts;
       int *p,*q;  
       p = rvars->grid2_add_map1;
       q = rvars->grid1_add_map1;

       for (int i=0; i < npts_dst; i++)
		data_dst[i] = 0;
       for (int i=0; i < rvars->num_links_map1; i++) {
	   	data_dst[p[i]]+=data_src[q[i]]*wgt_values[i];
       }
}


void conserv_remap_2D_operator::init_remap() 
{
	grid *grid1, *grid2;
	double *tmp_wgt;
	int i,j,tmp_nw;
 
       first_call = true;
	wgt_size = 6;
	rvars = new remap_vars(npts_src,npts_dst);
       
	link_add_src = new int[2*npts_src];
	link_add_dst = new int[2*npts_dst];

	cal_grid_bound_info(grid_name_src);
	cal_grid_bound_info(grid_name_dst);
       conserv_remap();

       num_wgts = rvars->num_links_map1;
       wgt_indx_src = rvars->grid1_add_map1;
       wgt_indx_dst = rvars->grid2_add_map1;
       wgt_values = new double[rvars->num_links_map1];
	   
       tmp_wgt = rvars->wts_map1;
       tmp_nw = rvars->num_wts;

	for (i=0, j=0; i<rvars->num_links_map1; i++) {
           wgt_values[i] = tmp_wgt[j];
           j = j + tmp_nw;
        }
	   
}


conserv_remap_2D_operator::~conserv_remap_2D_operator() {
       delete	rvars;
	delete	[]link_add_src;
	delete	[]link_add_dst;
}

conserv_remap_2D_operator::conserv_remap_2D_operator(char *grid1_name,  char *grid2_name,  char *alg_name) 
: Common_remap(grid1_name, grid2_name, alg_name)
{
	lthresh = false;
	luse_last = false;
	avoid_pole_count = 0;
	avoid_pole_offset = TINY;
	first_call = true;
}

void conserv_remap_2D_operator:: conserv_remap()
{
    // local variables
    int max_subseg = 10000; // max number of subsegments per segment to prevent infinite loop

    int grid1_add;          // current linear address for grid1 cell
    int grid2_add;          // current linear address for grid2 cell
    int grid_add;

    int min_add;            // addresses for restricting search of
	int max_add;            // destination grid

	int i, l, n, tmp_n, m, nwgt;            // generic counters
    int corner;             // corner of cell that segment starts from
    int next_corn;          // corner of cell that segment ends on
    int num_subseg;         // number of subsegments 

    bool lcoinc;         // flag for coincident segments
    bool lrevers;        // flag for reversing direction of segment
    bool lbegin;         // flag for first integration of a segment

    // allocatable 
    bool * srch_mask;    // mask for restricting searches

    double  intrsct_lat, intrsct_lon;                 // lat/lon of next intersect
    double  beglat, endlat, beglon, endlon;    // endpoints of current seg.
    double  norm_factor;                                // factor for normalizing wts

    // allocatable 
    double * grid1_centroid_lat;   // centroid coords on each grid
	double * grid1_centroid_lon; 
    double * grid2_centroid_lat;
	double * grid2_centroid_lon;  
       
    double begseg[2];      // begin lat/lon for full segment

    double weights[6];     // local wgt array

	int x,y;

    //for test
    pole_count = 0;
    inter_count = 0;
    store_count1 = 0;
    store_count2 = 0;
    store_count3 = 0;
    store_count4 = 0;
    store_count_all = 0;
    store_count_return1 = 0;
    store_count_return2 = 0;
	   
    //initialize centroid arrays
    grid1_centroid_lat = new double[npts_src];
    grid1_centroid_lon = new double[npts_src];
    grid2_centroid_lat = new double[npts_dst];
    grid2_centroid_lon = new double[npts_dst];
	   
	
    for (n=0; n<npts_src; n++) {
		grid1_centroid_lat[n] = ZERO;
		grid1_centroid_lon[n] = ZERO;
		}

	for (n=0; n<npts_dst; n++) {
		grid2_centroid_lat[n] = ZERO;
		grid2_centroid_lon[n] = ZERO;
		}

    // ntegrate around each cell on grid1
	srch_mask = new bool[npts_dst];


   //srch_add = new int[default_size]; 
   //srch_corner_lat = new double[num_vertexes_dst*default_size];
   //srch_corner_lon = new double[num_vertexes_dst*default_size];
	   

	for(grid1_add = 0; grid1_add<npts_src; grid1_add++) {	

       //printf("grid1_add: %d\n", grid1_add);
		
		grid_add  = grid1_add*4;
		//restrict searches first using search bins
		min_add  = npts_dst-1;     //min_add  = npts_dst*4;
		max_add = 0;
	  
		for (n=0; n<num_srch_bins; n++) {
			m=2*n;
			if (grid1_add >= bin_addr_src[m]  &&  grid1_add <= bin_addr_src[m+1]) {
				min_add  = MIN(min_add, bin_addr_dst[m]);
				max_add = MAX(max_add, bin_addr_dst[m+1]);
				}
			}
              // if(grid1_add == 1)
		//printf("grid1_add   min_add : %d    max_add : %d",min_add,max_add);
		   
		// further restrict searches using bounding boxes
		num_srch_cells = 0;
		for(grid2_add = min_add; grid2_add<=max_add; grid2_add++) { 
			m = grid2_add*4;
			srch_mask[grid2_add]= 
				(bound_box_dst[m] <= bound_box_src[grid_add+1]) && 
				(bound_box_dst[m+1] >= bound_box_src[grid_add]) && 
				(bound_box_dst[m+2] <= bound_box_src[grid_add+3]) && 
				(bound_box_dst[m+3] >= bound_box_src[grid_add+2]);

			if (srch_mask[grid2_add])  num_srch_cells++;
			}
			  
		//create search arrays
		//printf("!!!  \n");
              //printf("grid1_add: %d\n",grid1_add);
		//printf("num_srch_cells : %d\n",num_srch_cells);

              //printf("num_srch_cells: %d\n", num_srch_cells);
		srch_add = new int[num_srch_cells]; 
		//printf("ttttttt\n");
		
		srch_corner_lat = new double[num_vertexes_dst*num_srch_cells];
		srch_corner_lon = new double[num_vertexes_dst*num_srch_cells];

		srch_corner_lat_dim1 = num_vertexes_dst;
		srch_corner_lat_dim2 = num_srch_cells;
			  
		n = 0;
		// gather1
		for(grid2_add=min_add; grid2_add<=max_add; grid2_add++) {
			if (srch_mask[grid2_add]) {
				srch_add[n] = grid2_add;
				m = n * num_vertexes_dst;
				l   = grid2_add * num_vertexes_dst;
				for ( i=0; i<num_vertexes_dst; i++ ) {
					srch_corner_lat[i+m]  = vertex_lats_dst[i+l];
					srch_corner_lon[i+m] =  vertex_lons_dst[i+l];
					}
				n++;
				}
			}
		
	      /*
	      if(grid1_add==3670) {
		  	printf("3670  srch_corner_lat:  \n");
			for(n=0; n<num_vertexes_dst*num_srch_cells; n++)
				printf("%f   \n",srch_corner_lat[n]);	
                     printf("srch_corner_lat is  over\n");
			}
             */
     
		// integrate around this cell
		for (corner=0; corner<num_vertexes_src; corner++) {
			next_corn = (corner + 1)%num_vertexes_src;   // 1 -- num_vertexes_src-1 -- 0

			// define endpoints of the current segment
			beglat = vertex_lats_src[corner+grid_add];
			endlat = vertex_lats_src[next_corn+grid_add];
			beglon = vertex_lons_src[corner+grid_add];
			endlon = vertex_lons_src[next_corn+grid_add];

			lrevers = false;

			// to ensure exact path taken during both sweeps, always integrate segments in the same direction (SW to NE).
			if ((endlat < beglat) || (endlat == beglat && endlon < beglon)) {
				beglat  = vertex_lats_src[next_corn+grid_add];
				beglon  = vertex_lons_src[next_corn+grid_add];
				endlat  = vertex_lats_src[corner+grid_add];
				endlon  = vertex_lons_src[corner+grid_add];
				lrevers = true;
				}

			begseg[0] = beglat;
			begseg[1] = beglon;
			lbegin = true;
			num_subseg = 0;

			//if this is a constant-longitude segment, skip the rest since the line integral contribution will be zero.
			if (endlon != beglon) { 
                             //printf("if (endlon != beglon)   corner = %d\n",corner);
				
				// integrate along this segment, detecting intersections and computing the line integral for each sub-segment
				while (beglat != endlat || beglon != endlon) {
					// prevent infinite loops if integration gets stuck near cell or threshold boundary
					num_subseg = num_subseg + 1;
					if (num_subseg > max_subseg) {
						//printf("integration stalled: num_subseg exceeded limit\n Stop!!!!!\n");
						exit(2);
						}

					
					//printf("grid2_add:    %d\n ", grid2_add);
					//printf("intrsct_lat:    %f\n ", intrsct_lat);
					//printf("intrsct_lon:    %f\n ", intrsct_lon);
					//printf("beglat:    %f\n ", beglat);
					//printf("beglon:    %f\n ", beglon);	  
					//printf("endlat:    %f\n ", endlat);
					//printf("endlon:    %f\n ", endlon);
					//printf("begseg:    %f   %f\n ", begseg[0],begseg[1]);
						
					// find next intersection of this segment with a grid line on grid 2.
			              lcoinc = intersection(&grid2_add,&intrsct_lat,&intrsct_lon, beglat, beglon, endlat, endlon, begseg, lbegin, lrevers);
					//printf("after grid2_add:    %d\n ", (grid2_add+1));
					//printf("after begseg:    %f   %f\n ", begseg[0],begseg[1]);
					//if(lcoinc)  printf("after locinc :  T\n ");
					//else  printf("after locinc : F\n ");
						
					
					lbegin = false;

					// compute line integral for this subsegment.
					if (grid2_add != -1) 
						line_integral(weights, rvars->num_wts, beglon, intrsct_lon, beglat, intrsct_lat, center_lats_src[grid1_add], center_lons_src[grid1_add], center_lats_dst[grid2_add], center_lons_dst[grid2_add]);
					else 
						line_integral(weights, rvars->num_wts, beglon, intrsct_lon, beglat, intrsct_lat, center_lats_src[grid1_add], center_lons_src[grid1_add], center_lats_src[grid1_add], center_lons_src[grid1_add]);
			
					// if integrating in reverse order, change sign of weights
					
					if (lrevers)  for ( i=0; i<6; i++) weights[i] = -weights[i];

					// store the appropriate addresses and weights. also add contributions to cell areas and centroids.

                                   if ( grid2_add != -1 && mask_src[grid1_add]) { 
						store_count1++;
						store_link_cnsrv(grid1_add, grid2_add, weights);
						frac_src[grid1_add] = frac_src[grid1_add] +  weights[0];
						frac_dst[grid2_add] = frac_dst[grid2_add] +  weights[rvars->num_wts];
						}

					area_src[grid1_add] = area_src[grid1_add] + weights[rvars->num_wts];
                                   grid1_centroid_lat[grid1_add] = grid1_centroid_lat[grid1_add] + weights[rvars->num_wts+1];
                                   grid1_centroid_lon[grid1_add] = grid1_centroid_lon[grid1_add] + weights[rvars->num_wts+2];

					// reset beglat and beglon for next subsegment. 
					beglat = intrsct_lat;
					beglon = intrsct_lon;
					}// end do
				} 

			// end of segment
			
		}//  end do
             
		//finished with this cell: deallocate search array and start on next cell


		delete [] srch_add;
		delete [] srch_corner_lat;
		delete [] srch_corner_lon;

	} //end do

       delete []  srch_mask;

	//printf("pole_count : %d\n ", pole_count);
  
       // integrate around each cell on grid2
	srch_mask = new bool[npts_src];

	for (grid2_add=0; grid2_add<npts_dst; grid2_add++) {

		grid_add  = grid2_add*4;

		// restrict searches first using search bins
		min_add = npts_src-1;
		max_add = 0;

		for( n=0; n<num_srch_bins;n++ ) {
			m = n*2;
			if (grid2_add >= bin_addr_dst[m] && grid2_add <= bin_addr_dst[m+1]) {
				min_add = MIN(min_add, bin_addr_src[m]);
				max_add = MAX(max_add, bin_addr_src[m+1]);
				}
			}

		//if(grid2_add == 0)
		//printf("grid2_add   min_add : %d    max_add : %d",min_add,max_add);

		// further restrict searches using bounding boxes
		num_srch_cells = 0;
		for( grid1_add = min_add; grid1_add<=max_add; grid1_add++ ) {
			m = grid1_add*4;
			srch_mask[grid1_add]= 
				(bound_box_src[m] <= bound_box_dst[grid_add+1]) && 
				(bound_box_src[m+1] >= bound_box_dst[grid_add]) && 
				(bound_box_src[m+2] <= bound_box_dst[grid_add+3]) && 
				(bound_box_src[m+3] >= bound_box_dst[grid_add+2]);

			if (srch_mask[grid1_add])  num_srch_cells++;
			}

		srch_add = new int[num_srch_cells]; 
		srch_corner_lat = new double[num_vertexes_src*num_srch_cells];
		srch_corner_lon = new double[num_vertexes_src*num_srch_cells];
		srch_corner_lat_dim1 = num_vertexes_src;
		srch_corner_lat_dim2 = num_srch_cells; 
		

		n = 0;
		// gather2
		for(grid1_add=min_add; grid1_add<=max_add; grid1_add++) {
			if (srch_mask[grid1_add]) {
				srch_add[n] = grid1_add;
				m = n * num_vertexes_src;
				l = grid1_add * num_vertexes_src;
				for ( i=0; i<num_vertexes_src; i++ ) {
					srch_corner_lat[i+m]  = vertex_lats_src[i+l];
					srch_corner_lon[i+m] =  vertex_lons_src[i+l];
					}
				n++;
				}
			}

		// integrate around this cell
		for (corner=0; corner<num_vertexes_dst; corner++) {
			next_corn = (corner + 1)%num_vertexes_dst;   // 1 -- num_vertexes_src-1 -- 0

			// define endpoints of the current segment
			beglat = vertex_lats_dst[corner+grid_add];
			endlat = vertex_lats_dst[next_corn+grid_add];
			beglon = vertex_lons_dst[corner+grid_add];
			endlon = vertex_lons_dst[next_corn+grid_add];

			lrevers = false;

			// to ensure exact path taken during both sweeps, always integrate in the same direction
			if ((endlat < beglat) || (endlat == beglat && endlon < beglon)) {
				beglat  = vertex_lats_dst[next_corn+grid_add];
				beglon  = vertex_lons_dst[next_corn+grid_add];
				endlat  = vertex_lats_dst[corner+grid_add];
				endlon  = vertex_lons_dst[corner+grid_add];
				lrevers = true;
 				}
                     //printf("corner:%d\n", (corner+1));
			//printf("beglat:%f\n", beglat);
			//printf("beglon:%f\n", beglon);
			//printf("endlat:%f\n", endlat);
			//printf("endlon:%f\n", endlon);
			

			begseg[0] = beglat;
			begseg[1] = beglon;
			lbegin = true;

			// if this is a constant-longitude segment, skip the rest since the line integral contribution will be zero.
			if (endlon != beglon) {
				num_subseg = 0;

				// integrate along this segment, detecting intersections and computing the line integral for each sub-segment
				while (beglat!=endlat || beglon!=endlon) {
					// prevent infinite loops if integration gets stuck near cell or threshold boundary 
					num_subseg = num_subseg + 1;

					if (num_subseg > max_subseg) {
						//printf("integration stalled: num_subseg exceeded limit\n Stop\n");
						exit(2);
						}

					// find next intersection of this segment with a line on grid 2.
			              lcoinc = intersection(&grid1_add,&intrsct_lat,&intrsct_lon, beglat, beglon, endlat, endlon, begseg, lbegin, lrevers);
				       lbegin = false;

			       	// compute line integral for this subsegment.

 				       if (grid1_add != -1)
					   	line_integral(weights, rvars->num_wts, beglon, intrsct_lon, beglat, intrsct_lat, center_lats_src[grid1_add], center_lons_src[grid1_add], center_lats_dst[grid2_add], center_lons_dst[grid2_add]);
			 	      else
					 	line_integral(weights, rvars->num_wts, beglon, intrsct_lon, beglat, intrsct_lat, center_lats_dst[grid2_add], center_lons_dst[grid2_add], center_lats_dst[grid2_add], center_lons_dst[grid2_add]);
		

				     if (lrevers) 
					    for ( i=0; i<6; i++) weights[i] = -weights[i];


				     if (!lcoinc && grid1_add != -1) {
					   if (mask_src[grid1_add]) {
					   	  store_count2++;
			              	  store_link_cnsrv(grid1_add, grid2_add, weights);
						  frac_src[grid1_add] = frac_src[grid1_add] + weights[0];
						  frac_dst[grid2_add] = frac_dst[grid2_add] + weights[rvars->num_wts];
						  }
					   }

				     area_dst[grid2_add] = area_dst[grid2_add] + weights[rvars->num_wts];
				     grid2_centroid_lat[grid2_add] = grid2_centroid_lat[grid2_add] + weights[rvars->num_wts+1];
				     grid2_centroid_lon[grid2_add] = grid2_centroid_lon[grid2_add] + weights[rvars->num_wts+2];

				     // reset beglat and beglon for next subsegment.
				     beglat = intrsct_lat;
				     beglon = intrsct_lon;
				     }
		        }
			// end of segment
		}
		// finished with this cell: deallocate search array and start on next cell
		delete [] srch_add;
	    delete [] srch_corner_lat;
	    delete [] srch_corner_lon;
	}  //end do

	delete [] srch_mask;

	//printf("pole_count : %d\n ", pole_count);
	// North Pole
       
	weights[0] =  PI2;
	weights[1] =  PI*PI;
	weights[2] =  ZERO;
	weights[3] =  PI2;
	weights[4] =  PI*PI;
	weights[5] =  ZERO;

	grid1_add = -1;

	//pole_loop1: do n=1,npts_src
	for ( n=0; n<npts_src; n++) { 
		if (area_src[n] < -THREE*PIH && center_lats_src[n] > ZERO) {
			grid1_add = n;
			break;
			}
		} // end do pole_loop1

	grid2_add = -1;

       //pole_loop2: do n=1,npts_dst
	for ( n=0; n<npts_dst; n++) {
		if (area_dst[n] < -THREE*PIH && center_lats_dst[n] > ZERO) {
			grid2_add = n;
			break;
			}
		} // end do pole_loop2

	if (grid1_add !=-1) {
		area_src[grid1_add] = area_src[grid1_add] + weights[0];
		grid1_centroid_lat[grid1_add] = grid1_centroid_lat[grid1_add] + weights[1];
		grid1_centroid_lon[grid1_add] = grid1_centroid_lon[grid1_add] + weights[2];
		}

	if (grid2_add !=-1) { 
		area_dst[grid2_add] = area_dst[grid2_add] + weights[rvars->num_wts+0];
		grid2_centroid_lat[grid2_add] = grid2_centroid_lat[grid2_add] + weights[rvars->num_wts+1];
		grid2_centroid_lon[grid2_add] = grid2_centroid_lon[grid2_add] + weights[rvars->num_wts+2];
		}

	if (grid1_add != -1 && grid2_add !=-1) {
		store_count3++;
              store_link_cnsrv(grid1_add, grid2_add, weights);
		frac_src[grid1_add] = frac_src[grid1_add] + weights[0];
		frac_dst[grid2_add] = frac_dst[grid2_add] + weights[rvars->num_wts];
		}

	// South Pole
	weights[0] =  PI2;
	weights[1] =  -PI*PI;
	weights[2] =  ZERO;
	weights[3] =  PI2;
	weights[4] =  -PI*PI;
	weights[5] =  ZERO;

	grid1_add = -1;
       //pole_loop3: do n=1,npts_src
	for ( n=0; n<npts_src; n++ ) {
		if (area_src[n] < -THREE*PIH && center_lats_src[n] > ZERO){
			grid1_add = n;
			break;
			}
		} //end do pole_loop3

	grid2_add = -1;

	//pole_loop4: do n=1,npts_dst
	for ( n=0; n<npts_dst; n++ ) {
		if (area_dst[n] < -THREE*PIH && center_lats_dst[n] > ZERO){
			grid2_add = n;
			break;
			}
		} // end do pole_loop4

	if (grid1_add !=-1) {
		area_src[grid1_add] = area_src[grid1_add] + weights[0];
		grid1_centroid_lat[grid1_add] = grid1_centroid_lat[grid1_add] + weights[1];
		grid1_centroid_lon[grid1_add] = grid1_centroid_lon[grid1_add] + weights[2];
		}

	if (grid2_add !=-1) {
		area_dst[grid2_add] = area_dst[grid2_add] + weights[rvars->num_wts];
		grid2_centroid_lat[grid2_add] = grid2_centroid_lat[grid2_add] + weights[rvars->num_wts+1];
		grid2_centroid_lon[grid2_add] = grid2_centroid_lon[grid2_add] + weights[rvars->num_wts+2];
		}

	if (grid1_add != -1 && grid2_add != -1) {
		store_count4++;
		store_link_cnsrv(grid1_add, grid2_add, weights);
		frac_src[grid1_add] = frac_src[grid1_add] +  weights[0];
		frac_dst[grid2_add] = frac_dst[grid2_add] +  weights[rvars->num_wts];
		}

	// finish centroid computation
	for(i=0; i<npts_src; i++)
		if(area_src[i]!=ZERO) {
			grid1_centroid_lat[i] = grid1_centroid_lat[i]/area_src[i];
			grid1_centroid_lon[i] = grid1_centroid_lon[i]/area_src[i];
			}

	for(i=0; i<npts_dst; i++)
		if(area_dst[i]!=ZERO) {
			grid2_centroid_lat[i] = grid2_centroid_lat[i]/area_dst[i];
			grid2_centroid_lon[i] = grid2_centroid_lon[i]/area_dst[i];
			}

	// include centroids in weights and normalize using destination area if requested
       for (n=0; n<rvars->num_links_map1; n++) {
	   	grid1_add = rvars->grid1_add_map1[n];
		grid2_add = rvars->grid2_add_map1[n];

	       tmp_n = n* rvars->num_wts;
		   
		for (nwgt=0;nwgt<rvars->num_wts;nwgt++) {
			
			weights[nwgt] = rvars->wts_map1[nwgt+tmp_n];
			if (rvars->num_maps >1) 
				weights[rvars->num_wts+nwgt] = rvars->wts_map2[nwgt+tmp_n];
			}

		if(rvars->norm_opt==NORM_OPT_DSTAREA) {
			if (area_dst[grid2_add]!= ZERO) {
				if (rvars->luse_grid2_area) 
					norm_factor = ONE/grid2_area_in[grid2_add];
				else
					norm_factor = ONE/area_dst[grid2_add];
				}
			else
				norm_factor = ZERO;
			}
		else if(rvars->norm_opt==NORM_OPT_FRCAREA){
			if (frac_dst[grid2_add] != ZERO) {
				if (rvars->luse_grid2_area)
					norm_factor = area_dst[grid2_add]/(frac_dst[grid2_add]*grid2_area_in[grid2_add]);
				else
					norm_factor = ONE/frac_dst[grid2_add];
				}
			else
				norm_factor = ZERO;
			}
		else
			if (rvars->norm_opt==NORM_OPT_NONE) norm_factor = ONE;

		rvars->wts_map1[0+tmp_n] =  weights[0] * norm_factor;
		rvars->wts_map1[1+tmp_n] = (weights[1] - weights[0]* grid1_centroid_lat[grid1_add])* norm_factor;
		rvars->wts_map1[2+tmp_n] = (weights[2] - weights[0]* grid1_centroid_lon[grid1_add])* norm_factor;

		if (rvars->num_maps > 1) {
			if(rvars->norm_opt==NORM_OPT_DSTAREA) {
				if (area_src[grid1_add]!= ZERO) {
					if (rvars->luse_grid1_area)
						norm_factor = ONE/grid1_area_in[grid1_add];
					else
						norm_factor = ONE/area_src[grid1_add];
					}
				else
					norm_factor = ZERO;
				}	
			else if(rvars->norm_opt==NORM_OPT_FRCAREA){
				if (frac_src[grid1_add] != ZERO){
					if (rvars->luse_grid1_area)
						norm_factor = area_src[grid1_add]/(frac_src[grid1_add]*grid1_area_in[grid1_add]);
					else
						norm_factor = ONE/frac_src[grid1_add];
					}
				else
					norm_factor = ZERO;
				}
			else
				if (rvars->norm_opt==NORM_OPT_NONE) norm_factor = ONE;

		       rvars->wts_map2[0+tmp_n] =  weights[rvars->num_wts]*norm_factor;
		       rvars->wts_map2[1+tmp_n] = (weights[rvars->num_wts+1]-weights[rvars->num_wts]*grid2_centroid_lat[grid2_add])*norm_factor;
		       rvars->wts_map2[2+tmp_n] = (weights[rvars->num_wts+2]-weights[rvars->num_wts]*grid2_centroid_lon[grid2_add])*norm_factor;
		}
	}
	   
	for(i=0;i<npts_src;i++)
		if(area_src[i]!=ZERO)
			frac_src[i] = frac_src[i]/area_src[i];

	for(i=0; i<npts_dst; i++)
		if(area_dst[i]!=ZERO)
			frac_dst[i] = frac_dst[i]/area_dst[i];

	// perform some error checking on final weights
	for(n=0; n<npts_src; n++){
		grid1_centroid_lat[n] = ZERO;
		grid1_centroid_lon[n] = ZERO;
		}
	
	for(i=0;i<npts_dst;i++){
		grid2_centroid_lat[i] = ZERO;
		grid2_centroid_lon[i] = ZERO;
		}

	for(n=0; n<rvars->num_links_map1; n++){
		grid1_add = rvars->grid1_add_map1[n];
		grid2_add = rvars->grid2_add_map1[n];

		tmp_n = n*rvars->num_wts;

		grid2_centroid_lat[grid2_add] = grid2_centroid_lat[grid2_add] + rvars->wts_map1[tmp_n];
		if (rvars->num_maps > 1)
			grid1_centroid_lat[grid1_add] = grid1_centroid_lat[grid1_add] + rvars->wts_map2[tmp_n];
		}

	for(n=0; n<npts_dst; n++){
		if(rvars->norm_opt==NORM_OPT_DSTAREA)
			norm_factor = frac_dst[grid2_add];
		else if(rvars->norm_opt==NORM_OPT_FRCAREA)
			norm_factor = ONE;
		else{ 
			if(rvars->norm_opt==NORM_OPT_NONE){
				if (rvars->luse_grid2_area)
					norm_factor = grid2_area_in[grid2_add];
				else
					norm_factor = area_dst[grid2_add];
				}
			}
		}

	if (rvars->num_maps > 1) {
		for(n=0; n<npts_src; n++) {
			if(rvars->norm_opt==NORM_OPT_DSTAREA)
				norm_factor = frac_src[grid1_add];
			else if(rvars->norm_opt==NORM_OPT_FRCAREA)
				norm_factor = ONE;
			else { 
				if(rvars->norm_opt==NORM_OPT_NONE) {
					if (rvars->luse_grid1_area)
						norm_factor = grid1_area_in[grid1_add];
					else
						norm_factor = area_src[grid1_add];
					}
				}
			}
		}

       //printf("store_count1 : %d\n",store_count1);
	//printf("store_count2 : %d\n",store_count2);
	//printf("store_count3 : %d\n",store_count3);
	//printf("store_count4 : %d\n",store_count4);
	//printf("store_count_return1 : %d\n",store_count_return1);
	//printf("store_count_return2 : %d\n",store_count_return2);

	//delete [] grid1_centroid_lat;
	//delete [] grid1_centroid_lon; 
       //delete [] grid2_centroid_lat;
	//delete [] grid2_centroid_lon; 
}


/*
  this routine finds the next intersection of a destination grid 
  line with the line segment given by beglon, endlon, etc.
  a coincidence flag is returned if the segment is entirely 
  coincident with an ocean grid line.  the cells in which to search
  for an intersection must have already been restricted in the
  calling routine.
*/

bool conserv_remap_2D_operator::intersection(int * location, double * intrsct_lat, double * intrsct_lon, double beglat, double beglon, double endlat, double endlon, double * begseg, bool lbegin, bool lrevers) {
      
      // local variables 
      bool tmp_lcoinc;                                     // flag segments which are entirely coincident with a grid line
      int n, next_n,tmp_n,cell, pole_loc;
      bool loutside;                                           // flags points outside grid
      double lon1, lon2;                                  // local longitude variables for segment
      double lat1, lat2;                                   // local latitude  variables for segment
      double grdlon1, grdlon2;                      // local longitude variables for grid cell
      double grdlat1, grdlat2;                       // local latitude  variables for grid cell
      double vec1_lat, vec1_lon;                  // vectors and cross products used
      double vec2_lat, vec2_lon;                  // during grid search
      double cross_product; 
      double eps, offset;                         // small offset away from intersect
      double s1, s2, determ;                      // variables used for linear solve to
      double mat1, mat2, mat3, mat4, rhs1, rhs2;  // find intersection

      bool srch_loop;

      //printf("enter intersection\n");
      inter_count++;
    
      // initialize defaults, flags, etc.
      *location    = -1;
      tmp_lcoinc  = false;
      srch_loop   = true;
      *intrsct_lat = endlat;
      *intrsct_lon = endlon;

      if (num_srch_cells == 0) {
	  	//printf("end intersection   num_srch_cells == 0\n");
	  	return tmp_lcoinc;
		}

      if (beglat > north_thresh || beglat < south_thresh) {

		if (lthresh) {
                      //printf("lthresh is true!      last_loc: %d\n ",last_loc);
			*location = last_loc;
		}
		
		//printf("bf location: %d      \n ",(*location+1));
		tmp_lcoinc = pole_intersection(location,intrsct_lat,intrsct_lon,beglat,beglon,endlat,endlon,begseg,lrevers);
      //printf("af location:                 %d\n ",(*location+1));
	       	 	  
		if (lthresh) {
			last_loc  = *location;
            		intrsct_lat_off = *intrsct_lat;
            		intrsct_lon_off = *intrsct_lon;
			}
		return tmp_lcoinc;
		}
      
      loutside = false;
      
      if (lbegin) {
	  	lat1 = beglat;
		lon1 = beglon;
		}
      else {
	  	lat1 = intrsct_lat_off;
		lon1 = intrsct_lon_off;
		}
      lat2 = endlat;
      lon2 = endlon;
      if ((lon2-lon1) > THREE*PIH)
	  	lon2 = lon2 - PI2;
      else 
	  	if ((lon2-lon1) < (-THREE*PIH))
			lon2 = lon2 + PI2;

      s1 = ZERO;

      //search for location of this segment in ocean grid using cross product method to determine whether a point is enclosed by a cell      
      do { // srch_loop: do 

	       // if last segment crossed threshold, use that location
	  	if (lthresh) {
			for(cell=0; cell<num_srch_cells; cell++) // do cell=1,num_srch_cells   
				if (srch_add[cell] == last_loc) {
					*location = last_loc;
					eps = TINY;
					srch_loop=false;
					break;
					}
				}
		if(!srch_loop) {
			  //printf("exit srch_loop_1\n");
			  break;
			}
		
		// otherwise normal search algorithm
		for (cell=0; cell<num_srch_cells; cell++) { // cell_loop: do cell=1,num_srch_cells
			for (n=0; n<srch_corner_lat_dim1; n++){       // corner_loop: do n=1,srch_corners      
				next_n = (n+1)%srch_corner_lat_dim1;    //??? MOD(n,srch_corners) + 1
                            //here we take the cross product of the vector making up each cell side with the vector formed by the vertex and search point.  if all the cross products are positive, the point is contained in the cell.
                            tmp_n = cell*srch_corner_lat_dim1;

				vec1_lat = srch_corner_lat[next_n+tmp_n] - srch_corner_lat[n+tmp_n];
                            vec1_lon = srch_corner_lon[next_n+tmp_n] - srch_corner_lon[n+tmp_n];
                            vec2_lat = lat1 - srch_corner_lat[n+tmp_n];
                            vec2_lon = lon1 - srch_corner_lon[n+tmp_n];

				// if endpoint coincident with vertex, offset the endpoint

				if (vec2_lat == 0 && vec2_lon == 0) {
					lat1 = lat1 + (1.e-10)*(lat2-lat1);                   
					lon1 = lon1 + (1.e-10)*(lon2-lon1);
					vec2_lat = lat1 - srch_corner_lat[n+tmp_n];
					vec2_lon = lon1 - srch_corner_lon[n+tmp_n];
					}

	                    // check for 0,2pi crossings
	                    if (vec1_lon > PI) 
					vec1_lon = vec1_lon - PI2;  
			      else 
				       if (vec1_lon < (-PI))
					   	vec1_lon = vec1_lon + PI2;

			      if (vec2_lon > PI)
				  	vec2_lon = vec2_lon - PI2;
			      else 
				 	if (vec2_lon < (-PI))
						vec2_lon = vec2_lon + PI2;

		             cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat;    
					 
 //if the cross product for a side is zero, the point lies exactly on the side or the side is degenerate(zero length).  if degenerate, set the cross product to a positive number.  otherwise perform another cross product between the side and the segment itself. 
 //if this cross product is also zero, the line is coincident with the cell boundary - perform the dot product and only choose the cell if the dot product is positive (parallel vs anti-parallel).
     
			      if (cross_product == ZERO) {
				  	if (vec1_lat != ZERO || vec1_lon != ZERO) {
						vec2_lat = lat2 - lat1;
						vec2_lon = lon2 - lon1;

						if (vec2_lon > PI)
							vec2_lon = vec2_lon - PI2;
						else 
							if (vec2_lon < (-PI))
								vec2_lon = vec2_lon + PI2;

						cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat;
						}
					else
						cross_product = ONE;

					if (cross_product == ZERO){
						tmp_lcoinc = true;
						cross_product = vec1_lon*vec2_lon + vec1_lat*vec2_lat;
						if (lrevers) 
							cross_product = -cross_product;
						}
					}
				  
                              // if cross product is less than zero, this cell doesn't work
				  if (cross_product < ZERO) {
				  	//printf("exit corner_loop_1\n");
					break;
					}

			} //end do corner_loop

                     // if cross products all positive, we found the location
                     if (n >= srch_corner_lat_dim1) {
				*location = srch_add[cell];
				// if the beginning of this segment was outside the grid, invert the segment so the intersection found will be the first intersection with the grid
				if (loutside) {
					lat2 = beglat;
					lon2 = beglon;
					*location = 0;
					eps  = -TINY;
					}
				else
					eps  = TINY;
                            srch_loop = false;
				break;
				}
	             // otherwise move on to next cell
	             } //end cell_loop
              if(!srch_loop)
		{
			  //printf("exit srch_loop_2\n");
			  break;
		}
          

              // if still no cell found, the point lies outside the grid. take some baby steps along the segment to see if any part of the segment lies inside the grid.  
             
		loutside = true;
		s1 = s1 + 0.001;
		lat1 = beglat + s1*(endlat - beglat);
		lon1 = beglon + s1*(lon2   - beglon);

		// reached the end of the segment and still outside the grid return no intersection

		if (s1 >= ONE) {
                     //printf("end intersection   s1 >= ONE\n");
			return tmp_lcoinc;
			}
	} while(srch_loop);  
	//end srch_loop
     
       //now that a cell is found, search for the next intersection. loop over sides of the cell to find intersection with side must check all sides for coincidences or intersections   
       for (n=0; n<srch_corner_lat_dim1; n++)  { // intrsct_loop: do n=1,srch_corners
	   	next_n = (n+ 1)%srch_corner_lat_dim1;
              tmp_n = cell*srch_corner_lat_dim1;
		grdlon1 = srch_corner_lon[n+tmp_n];
		grdlon2 = srch_corner_lon[next_n+tmp_n];
		grdlat1 = srch_corner_lat[n+tmp_n];
		grdlat2 = srch_corner_lat[next_n+tmp_n];

		// set up linear system to solve for intersection
        
              mat1 = lat2 - lat1;
              mat2 = grdlat1 - grdlat2;
              mat3 = lon2 - lon1;
              mat4 = grdlon1 - grdlon2;
              rhs1 = grdlat1 - lat1;
              rhs2 = grdlon1 - lon1;

		if (mat3 > PI)
			mat3 = mat3 - PI2;
		else 
			if (mat3 < (-PI))
				mat3 = mat3 + PI2;

		if (mat4 > PI)
			mat4 = mat4 - PI2;
		else
			if (mat4 < (-PI))
				mat4 = mat4 + PI2;

		if (rhs2 > PI) 
			rhs2 = rhs2 - PI2;
		else 
			if (rhs2 < (-PI)) 
				rhs2 = rhs2 + PI2;

		determ = mat1*mat4 - mat2*mat3;


// if the determinant is zero, the segments are either parallel or coincident.  coincidences were detected above so do nothing.
// if the determinant is non-zero, solve for the linear parameters s for the intersection point on each line segment.
// if 0<s1,s2<1 then the segment intersects with this side.return the point of intersection (adding a small number so the intersection is off the grid line).

              if (fabs(determ) > 1.e-30) {
			s1 = (rhs1*mat4 - mat2*rhs2)/determ;
			s2 = (mat1*rhs2 - rhs1*mat3)/determ;

			if (s2 >= ZERO && s2 <= ONE && s1 > ZERO && s1 <= ONE) {
				// recompute intersection based on full segment so intersections are consistent for both sweeps
                            if (!loutside) {
					mat1 = lat2 - begseg[0];
					mat3 = lon2 - begseg[1];
					rhs1 = grdlat1 - begseg[0];
					rhs2 = grdlon1 - begseg[1];
					}
				else {
					mat1 = begseg[0] - endlat;
					mat3 = begseg[1] - endlon;
					rhs1 = grdlat1 - endlat;
					rhs2 = grdlon1 - endlon;
				}

				if (mat3 > PI)
					mat3 = mat3 - PI2;
				else 
					if (mat3 < (-PI))
						mat3 = mat3 + PI2;

				if (rhs2 > PI)
					rhs2 = rhs2 - PI2;
				else
				       if (rhs2 < (-PI))
					   	rhs2 = rhs2 + PI2;

				determ = mat1*mat4 - mat2*mat3;
            // sometimes due to roundoff, the previous  determinant is non-zero, but the lines are actually coincident.  if this is the case, skip the rest.

				if (determ != ZERO) {
					s1 = (rhs1*mat4 - mat2*rhs2)/determ;
					s2 = (mat1*rhs2 - rhs1*mat3)/determ;

					offset = s1 + eps/determ;

					if (offset > ONE) offset = ONE;

					if (!loutside) {
						*intrsct_lat = begseg[0] + mat1*s1;
						*intrsct_lon = begseg[1] + mat3*s1;
						intrsct_lat_off = begseg[0] + mat1*offset;
						intrsct_lon_off = begseg[1] + mat3*offset;
						}
					else{  
						*intrsct_lat = endlat + mat1*s1;
						*intrsct_lon = endlon + mat3*s1;
						intrsct_lat_off = endlat + mat1*offset;
						intrsct_lon_off = endlon + mat3*offset;
						}
					break;  //exit intrsct_loop
					}
				}
			}
          // no intersection this side, move on to next side  
	   } //end do intrsct_loop

         //if the segment crosses a pole threshold, reset the intersection to be the threshold latitude.  only check if this was not a threshold segment since sometimes coordinate transform can end up on other side of threshold again.

         if (lthresh) {
		if (*intrsct_lat < north_thresh || *intrsct_lat > south_thresh)
        	lthresh = false;
		}
	 else {
	 	if (lat1 > ZERO && *intrsct_lat > north_thresh){
			*intrsct_lat = north_thresh + TINY;
			intrsct_lat_off = north_thresh + eps*mat1;
			s1 = (*intrsct_lat - begseg[0])/mat1;
			*intrsct_lon     = begseg[1] + s1*mat3;
			intrsct_lon_off = begseg[1] + (s1+eps)*mat3;
			last_loc = *location;
			lthresh = true;
			}
		else {
			if (lat1 < ZERO && *intrsct_lat < south_thresh){
				*intrsct_lat = south_thresh - TINY;
				intrsct_lat_off = south_thresh + eps*mat1;
				s1 = (*intrsct_lat - begseg[0])/mat1;
				*intrsct_lon     = begseg[1] + s1*mat3;
				intrsct_lon_off = begseg[1] + (s1+eps)*mat3;
				last_loc = *location;
				lthresh = true;
				}
			}
		} 

	  // if reached end of segment, do not use x,y intersect on next entry

	  if (*intrsct_lat == endlat && *intrsct_lon == endlon) luse_last = false;

	  //printf("end intersection!\n");

        //printf("intersection : \n");
        //printf("location : %d                  inter_count : %d\n ", (*location+1), inter_count);
        //printf("intrsct_lat :  %lf\n ",*intrsct_lat);
        //printf("intrsct_lon : %lf\n ",* intrsct_lon);
        
         return tmp_lcoinc; 
}


/*
  this routine is identical to the intersection routine except that a coordinate transformation (using a Lambert azimuthal
  equivalent projection) is performed to treat polar cells more accurately.
*/
                                                                             
bool conserv_remap_2D_operator::pole_intersection(int *location, double *intrsct_lat, double * intrsct_lon, double beglat, double beglon, double endlat, double endlon, double * begseg, bool lrevers)
{
      //local variables
      bool tmp_lcoinc;                                              // flag segment coincident with grid line
      int n,next_n, tmp_n, cell, pole_loc; 
      bool loutside;                                                   // flags points outside grid
      double pi4, rns;                                                // north/south conversion
      double x1, x2;                                                  // local x variables for segment
      double y1, y2;                                                  // local y variables for segment
      double begx, begy;                                          // beginning x,y variables for segment
      double endx, endy;                                          // beginning x,y variables for segment
      double begsegx, begsegy;                               // beginning x,y variables for segment
      double grdx1, grdx2;                                       // local x variables for grid cell
      double grdy1, grdy2;                                       // local y variables for grid cell
      double vec1_y, vec1_x;                                  // vectors and cross products used
      double vec2_y, vec2_x;                                  // during grid search
      double cross_product, eps;                            // eps=small offset away from intersect
      double s1, s2, determ;                                    // variables used for linear solve to
      double mat1, mat2, mat3, mat4, rhs1, rhs2;  // find intersection

      double *srch_corner_x;                                  // x of each corner of srch cells
      double *srch_corner_y;                                  // y of each corner of srch cells

      bool srch_loop = true;

       pole_count = pole_count+1;
      //printf("pole_enter intersection\n");
      //if(pole_count==964)
		//	printf("*location 0 :  %d\n", *location);

      // initialize defaults, flags, etc.
      if (!lthresh) 
	  	*location = -1; 
	  
      *intrsct_lat = endlat;
      *intrsct_lon = endlon;
      tmp_lcoinc = false;
      loutside = false;
      s1 = ZERO;

      

      // convert coordinates
      srch_corner_x = new double[srch_corner_lat_dim1*srch_corner_lat_dim2]; 
      srch_corner_y = new double[srch_corner_lat_dim1*srch_corner_lat_dim2]; 
     
      if (beglat > ZERO) {
	  	pi4 = QUART*PI;
		rns = ONE;
		}
      else {
	  	pi4 = -QUART*PI;
              rns = -ONE;
		}

      if (luse_last){
	  	x1 = intrsct_x;
		y1 = intrsct_y;
		}
      else {
	  	x1 = rns*TWO*sin(pi4 - HALF*beglat)*cos(beglon);
              y1 = TWO*sin(pi4 - HALF*beglat)*sin(beglon);
              luse_last = true;
		}
      x2 = rns*TWO*sin(pi4 - HALF*endlat)*cos(endlon);
      y2 =       TWO*sin(pi4 - HALF*endlat)*sin(endlon);

      for (int i=0; i<srch_corner_lat_dim1*srch_corner_lat_dim2; i++){                                            //optimization
              srch_corner_x[i]  = rns*TWO*sin(pi4 - HALF*srch_corner_lat[i])*cos(srch_corner_lon[i]);
		srch_corner_y[i]  =       TWO*sin(pi4 - HALF*srch_corner_lat[i])*sin(srch_corner_lon[i]);
	  	}
	  	
      begx = x1;
      begy = y1;
      endx = x2;
      endy = y2;
      begsegx = rns*TWO*sin(pi4 - HALF*begseg[0])*cos(begseg[1]);
      begsegy =       TWO*sin(pi4 - HALF*begseg[0])*sin(begseg[1]);
      intrsct_x = endx;
      intrsct_y = endy;

/*
  search for location of this segment in ocean grid using cross
  product method to determine whether a point is enclosed by a cell
*/
     //srch_loop: do
     do {
        // if last segment crossed threshold, use that location
        if (lthresh)
	 {
			for(cell=0; cell<num_srch_cells; cell++) {       
				if (srch_add[cell] == *location) {
					eps = TINY;
					srch_loop=false;
					break;
					}
			}
     	 }

        if(!srch_loop) break;

        // otherwise normal search algorithm
        for (cell=0; cell<num_srch_cells; cell++)  {  // cell_loop: do cell=1,num_srch_cells                           //optimization
		for (n=0; n<srch_corner_lat_dim1; n++)  {  // corner_loop: do n=1,srch_corners     
			       next_n = (n+1)%srch_corner_lat_dim1;        //??? MOD(n,srch_corners) + 1
                            tmp_n = cell*srch_corner_lat_dim1;             //??? check
		              // here we take the cross product of the vector making up each cell side with the vector formed by the vertex and search point.  if all the cross products are positive, the point is contained in the cell.
                            vec1_x = srch_corner_x[next_n+tmp_n] - srch_corner_x[n+tmp_n];
			       vec1_y = srch_corner_y[next_n+tmp_n] - srch_corner_y[n+tmp_n]; 
			       vec2_x = x1 - srch_corner_x[n+tmp_n];
			       vec2_y = y1 - srch_corner_y[n+tmp_n];

				// if endpoint coincident with vertex, offset the endpoint

				if (vec2_x == 0 && vec2_y == 0){
					x1 = x1 + (1.0e-10)*(x2-x1);                     //????????
					y1 = y1 + (1.0e-10)*(y2-y1);
					vec2_x = x1 - srch_corner_x[n+tmp_n];
					vec2_y = y1 - srch_corner_y[n+tmp_n];
					}

				cross_product = vec1_x*vec2_y - vec2_x*vec1_y;

                            if (cross_product == ZERO) {
					if (vec1_x != ZERO || vec1_y != ZERO) {
						vec2_x = x2 - x1;
						vec2_y = y2 - y1;
						cross_product = vec1_x*vec2_y - vec2_x*vec1_y;
						}
					else
						cross_product = ONE;

					if (cross_product == ZERO) {
						tmp_lcoinc = true;
						cross_product = vec1_x*vec2_x + vec1_y*vec2_y;
						if (lrevers)  cross_product = -cross_product;
						}
					}

				// if cross product is less than zero, this cell doesn't work

				if (cross_product < ZERO) break;
				} //end do corner_loop

		        //if(pole_count==964)
				//        printf("nnnnnnnnnnnnnnnnnnnnnnn :  %d\n", n);
                 

			// if cross products all positive, we found the location
			if (n >= srch_corner_lat_dim1) {
				*location = srch_add[cell];
                            	  		   
			       //if(pole_count==964)
				   //     printf("*location3 :  %d\n", *location);
				if (loutside) {
					x2 = begx;
					y2 = begy;
					*location = -1;
					 //if(pole_count==964)
					 //	printf("*location3 :  %d\n", *location);
					eps  = -TINY;
					}
				else
					eps  = TINY;
				
                            srch_loop = false;
			       break;
				}
		     // otherwise move on to next cell
	} //end do cell_loop

       if (!srch_loop) break;
        
        // if no cell found, the point lies outside the grid. take some baby steps along the segment to see if any part of the segment lies inside the grid.  
        loutside = true;
        s1 = s1 + 0.001;
        x1 = begx + s1*(x2 - begx);
        y1 = begy + s1*(y2 - begy);

        // reached the end of the segment and still outside the grid return no intersection
        if (s1 >= ONE) {
		    delete [] srch_corner_x;
        	delete [] srch_corner_y;
        	luse_last = false;

		//	if(pole_count==964)
		//		printf("111111111 location : %d\n",*location);
        	return tmp_lcoinc;
        }

     }while (srch_loop);//end  srch_loop

	 //if(pole_count==964)
	 //	printf("111111111 s1 : %f\n", s1);
     
     
     // now that a cell is found, search for the next intersection. loop over sides of the cell to find intersection with side must check all sides for coincidences or intersections

     for (n=0; n<srch_corner_lat_dim1; n++)    // intrsct_loop: do n=1,srch_corners
     {          
	 	next_n = (n+1)%srch_corner_lat_dim1;     //next_n = mod(n,srch_corners) + 1
              tmp_n = cell*srch_corner_lat_dim1;          //??????????????
		grdy1 = srch_corner_y[n+tmp_n];
              grdy2 = srch_corner_y[next_n+tmp_n];
              grdx1 = srch_corner_x[n+tmp_n];
              grdx2 = srch_corner_x[next_n+tmp_n];
			  
              // set up linear system to solve for intersection
              mat1 = x2 - x1;
              mat2 = grdx1 - grdx2;
		mat3 = y2 - y1;
		mat4 = grdy1 - grdy2;
		rhs1 = grdx1 - x1;
		rhs2 = grdy1 - y1;
		determ = mat1*mat4 - mat2*mat3;

       
              // if the determinant is zero, the segments are either parallel or coincident.  coincidences were detected above so do nothing.
              // if the determinant is non-zero, solve for the linear parameters s for the intersection point on each line segment.
              // if 0<s1,s2<1 then the segment intersects with this side. return the point of intersection (adding a small number so the intersection is off the grid line).

	      if (fabs(determ) > 1.e-30) {
		  	s1 = (rhs1*mat4 - mat2*rhs2)/determ;
			s2 = (mat1*rhs2 - rhs1*mat3)/determ;

			if (s2 >= ZERO && s2 <= ONE && s1 > ZERO && s1 <= ONE) {
				// recompute intersection based on full segment so intersections are consistent for both sweeps
				if (!loutside) {
					mat1 = x2 - begsegx;
					mat3 = y2 - begsegy;
					rhs1 = grdx1 - begsegx;
					rhs2 = grdy1 - begsegy;
					}
				else {
					mat1 = x2 - endx;
					mat3 = y2 - endy;
					rhs1 = grdx1 - endx;
					rhs2 = grdy1 - endy;
					}

				determ = mat1*mat4 - mat2*mat3;

                            //sometimes due to roundoff, the previous  determinant is non-zero, but the lines are actually coincident.  if this is the case, skip the rest.
                            if (determ != ZERO) 
				{
					s1 = (rhs1*mat4 - mat2*rhs2)/determ;
					s2 = (mat1*rhs2 - rhs1*mat3)/determ;

					if (!loutside) {
						intrsct_x = begsegx + s1*mat1;
						intrsct_y = begsegy + s1*mat3;
						}
					else {  
						intrsct_x = endx + s1*mat1;
						intrsct_y = endy + s1*mat3;
						}

		                     // convert back to lat/lon coordinates
		                     *intrsct_lon = rns*atan2(intrsct_y,intrsct_x);
					if (*intrsct_lon < ZERO) 
						*intrsct_lon = *intrsct_lon + PI2;

					if (fabs(intrsct_x) > (1.0e-10))
					       *intrsct_lat = (pi4 - asin(rns*HALF*intrsct_x/cos(*intrsct_lon)))*TWO;
					else if (fabs(intrsct_y) > (1.0e-10)) 
						*intrsct_lat = (pi4 - asin(HALF*intrsct_y/sin(*intrsct_lon)))*TWO;
					else
						*intrsct_lat = TWO*pi4;

					// add offset in transformed space for next pass.
					if (s1 - eps/determ < ONE) {
						intrsct_x = intrsct_x - mat1*(eps/determ);
						intrsct_y = intrsct_y - mat3*(eps/determ);
						}
					else {
						if (!loutside) {
							intrsct_x = endx;
							intrsct_y = endy;
							*intrsct_lat= endlat;
                                                 *intrsct_lon = endlon;
							}
						else {
							intrsct_x = begsegx;
							intrsct_y = begsegy;
				                     *intrsct_lat = begseg[0];
				                     *intrsct_lon = begseg[1];
							}
						}

					break;  //exit intrsct_loop
				}
			}
		  }
     } //end do intrsct_loop
     
     delete [] srch_corner_x;
     delete [] srch_corner_y;

     //if segment manages to cross over pole, shift the beginning endpoint in order to avoid hitting pole directly (it is ok for endpoint to be pole point)
     if (fabs(intrsct_x) < 1.e-10 && fabs(intrsct_y) < 1.e-10 && (endx!=ZERO && endy!=0)) {
	 	if (avoid_pole_count > 2) {
			avoid_pole_count = 0;
			avoid_pole_offset = 10*avoid_pole_offset;   //10.*avoid_pole_offset;   ????
			}
		cross_product = begsegx*(endy-begsegy) - begsegy*(endx-begsegx);
		*intrsct_lat = begseg[0];

		if (cross_product*(*intrsct_lat) > ZERO) {
			*intrsct_lon = beglon    + avoid_pole_offset;
			begseg[1]   = begseg[1] + avoid_pole_offset;
			}
		else {
			*intrsct_lon = beglon    - avoid_pole_offset;
			begseg[1]   = begseg[1] - avoid_pole_offset;
			}

		avoid_pole_count = avoid_pole_count + 1;
		luse_last = false;
	}
      else {
	  	avoid_pole_count = 0;
		avoid_pole_offset = TINY;
		}

      //  if the segment crosses a pole threshold, reset the intersection to be the threshold latitude.  only check if this was not a threshold segment since sometimes coordinate transform can end up on other side of threshold again.
      if (lthresh)
      {
        if (*intrsct_lat < north_thresh || *intrsct_lat > south_thresh)
        	lthresh = false;
       }
      else 
      {
      	if (beglat > ZERO && *intrsct_lat < north_thresh)
      	{
      	   mat4 = endlat - begseg[0];
          mat3 = endlon - begseg[1];
          if (mat3 > PI) 
          	mat3 = mat3 - PI2;
          else 
          	if (mat3 < (-PI)) 
          		mat3 = mat3 + PI2;
          *intrsct_lat = north_thresh - TINY;
          s1 = (north_thresh - begseg[0])/mat4;
          *intrsct_lon = begseg[1] + s1*mat3;
          luse_last = false;
          lthresh = true;
        }
        else 
        {
        	if (beglat < ZERO && *intrsct_lat > south_thresh)
        	{
        		mat4 = endlat - begseg[0];
        		mat3 = endlon - begseg[1];
            if (mat3 > PI) 
            	mat3 = mat3 - PI2;
            else 
            	if (mat3 < (-PI)) 
            		mat3 = mat3 + PI2;
            *intrsct_lat = south_thresh + TINY;
            s1 = (south_thresh - begseg[0])/mat4;
            *intrsct_lon = begseg[1] + s1*mat3;
            luse_last = false;
            lthresh = true;
          }
        }
      }
      
     // if reached end of segment, do not use x,y intersect on next entry
     if (*intrsct_lat == endlat && *intrsct_lon == endlon) luse_last = false;

     //printf("end pole_intersection : \n");

	
     /*
              printf("pole_intersection : \n");
              printf("location : %d                  pole_count : %d\n ", (*location+1), pole_count);
              printf("intrsct_lat : %f\n ",*intrsct_lat);
              printf("intrsct_lon : %f\n ",* intrsct_lon);
	 
             if(lthresh)  printf("lthresh :  T\n ");
              else          printf("lthresh :  F\n ");	 
           */
     //printf("between   location = %d          pole_count:      %d\n", (*location+1), pole_count);
     //if (*location==-1)  printf("-1  pole_count:  %d     location = %d\n", pole_count);
     return tmp_lcoinc; 
}



/*
  this routine computes the line integral of the flux function 
  that results in the interpolation weights. the line is defined
  by the input lat/lon of the endpoints.
*/
void conserv_remap_2D_operator::line_integral(double *weights, int num_wts, double in_phi1, double in_phi2, double theta1, double theta2, double grid1_lat, double grid1_lon, double grid2_lat, double grid2_lon ){     
     /* local variables */   
     double dphi, sinth1, sinth2, costh1, costh2, fac, phi1, phi2, phidiff1, phidiff2, sinint, f1, f2, fint;

     //printf("enter line_integral\n");

     /*
       weights for the general case based on a trapezoidal approx to
       the integrals.
     */
     sinth1 = sin(theta1);
     sinth2 = sin(theta2);
     costh1 = cos(theta1);
     costh2 = cos(theta2);
  
     dphi = in_phi1 - in_phi2;
     if (dphi>PI)
     	  dphi = dphi - PI2;
     else 
       if (dphi<(-PI))
          dphi = dphi + PI2;
      
     dphi = HALF * dphi;
     
     /*
       the first weight is the area overlap integral. the second and
       fourth are second-order latitude gradient weights.
     */ 
      weights[0] = dphi*(sinth1 + sinth2);   
      weights[num_wts] = dphi*(sinth1 + sinth2);
      weights[1] = dphi*(costh1 + costh2 + (theta1*sinth1 + theta2*sinth2));    
      weights[num_wts+1] = dphi*(costh1 + costh2 + (theta1*sinth1 + theta2*sinth2));

      
     /*
       the third and fifth weights are for the second-order phi gradient
       component.  must be careful of longitude range.
     */

      f1 = HALF*(costh1*sinth1 + theta1);
      f2 = HALF*(costh2*sinth2 + theta2);

      phi1 = in_phi1 - grid1_lon;
      if (phi1 > PI)
        phi1 = phi1 - PI2;
      else
	  if (phi1 < (-PI))
	  	phi1 = phi1 + PI2;

      phi2 = in_phi2 - grid1_lon;
      
      if (phi2 > PI)
        phi2 = phi2 - PI2;
      else 
	  if (phi2 < (-PI))
	  	phi2 = phi2 + PI2;
      
      if ((phi2-phi1) < PI && (phi2-phi1) > (-PI))
        weights[2] = dphi*(phi1*f1 + phi2*f2);          //weights(3)
      else {
	  	if (phi1 > ZERO)
			fac = PI;
		else
			fac = -PI;
              fint = f1 + (f2-f1)*(fac-phi1)/fabs(dphi);
              weights[2] = HALF*(phi1*(phi1-fac)*f1 - phi2*(phi2+fac)*f2 + fac*(phi1+phi2)*fint);   //weights(3)
              }

      phi1 = in_phi1 - grid2_lon;
      if (phi1 > PI)
	  	phi1 = phi1 - PI2;
      else
	  	if (phi1 < (-PI)) 
			phi1 = phi1 + PI2;
  
      phi2 = in_phi2 - grid2_lon;
      if (phi2 > PI)
	  	phi2 = phi2 - PI2;
      else 
	  	if (phi2 < (-PI))
			phi2 = phi2 + PI2;
      

      if ((phi2-phi1) <  PI && (phi2-phi1) > (-PI))
	  	weights[num_wts+2] = dphi*(phi1*f1 + phi2*f2);    //weights(num_wts+3)
      else {
	  	if (phi1 > ZERO)
			fac = PI;
              else
			fac = -PI;
              fint = f1 + (f2-f1)*(fac-phi1)/fabs(dphi);
              weights[num_wts+2] = HALF*(phi1*(phi1-fac)*f1 - phi2*(phi2+fac)*f2 + fac*(phi1+phi2)*fint);
	       }

      //printf("aaa\n");  
      // for(int c=0; c<6; c++)
      //	    printf("weights  :   %f\n",weights[c]);

	  //printf("end line_integral !\n");
}




/*
  this routine stores the address and weight for this link in
  the appropriate address and weight arrays and resizes those
  arrays if necessary.
*/    
  
void conserv_remap_2D_operator::store_link_cnsrv(int add1, int add2, double * weights)
{
	 
     int min_link, max_link, tmp_c;        // link index
     int num_wts = rvars ->num_wts;
     int i, h, j;

	 store_count_all++;


     //printf("store_link_cnsrv\n");
   
     // if all weights are zero, do not bother storing the link
     for(i = 0; i<wgt_size; i++)
        if (weights[i] != ZERO)
        	break;
     if (i == wgt_size){
	 	store_count_return1++;
        return;
     	}
	
    
     // restrict the range of links to search for existing links
     if (first_call) {
	 	for(i=0;i<2*npts_src;i++)
			 link_add_src[i]=-1;
		for(i=0;i<2*npts_dst;i++)
			 link_add_dst[i]=-1;

		first_call = false;
        min_link = 0;
        max_link = -1;
	    }
      else {
	  	min_link = MIN(link_add_src[2*add1], link_add_dst[2*add2]);
	    max_link = MAX(link_add_src[2*add1+1], link_add_dst[2*add2+1]);
		if (min_link == -1) {
			min_link = 0;
            max_link = -1;
			}
		}

      //if(store_count_all==9)
	  //
	  //{
	  //printf("link_add_src[2*add1]: %d     link_add_dst[2*add2]: %d\n",link_add_src[2*add1]+1, link_add_dst[2*add2]+1);
	  //printf("link_add_src[2*add1+1]: %d     link_add_dst[2*add2+1]: %d\n",link_add_src[2*add1+1]+1, link_add_dst[2*add2+1]+1);
	  //}
	  	
	  
	  //printf("num_links_map1: %d\n",rvars->num_links_map1);
      //printf("add1: %d\n",add1+1);
	  //printf("add2: %d\n",add2+1);
      //printf("min_link: %d\n",min_link+1);
	  //printf("max_link: %d\n",max_link+1);

      for (i=min_link; i<=max_link; i++){
	  	if (add1 == rvars ->grid1_add_map1[i] && add2 == rvars ->grid2_add_map1[i]) { 
				for (j=0; j<num_wts; j++)
					rvars ->wts_map1[num_wts*i+j] += weights[j];
				if (rvars ->num_maps == 2)
					for (h=0; h<num_wts; h++)
						rvars ->wts_map2[num_wts*i+h] += weights[num_wts+h];
                		store_count_return2++;
				return;
        }     
	 }

/*
  if the link does not yet exist, increment number of links and 
  check to see if remap arrays need to be increased to accomodate 
  the new link.  then store the link.
*/
      rvars->num_links_map1++;
      if (rvars->num_links_map1 > rvars ->max_links_map1) 
	  	rvars->reset_remap_vars(1, rvars->resize_increment);

      tmp_c = rvars->num_links_map1-1;
      rvars->grid1_add_map1[tmp_c] = add1;
      rvars->grid2_add_map1[tmp_c] = add2;

      //printf("tmp_c    %d \n", tmp_c+1);
      //printf("add2      %d \n", add2+1);

      tmp_c = tmp_c*num_wts;
      for (i=0; i<num_wts; i++)
	  	rvars->wts_map1[tmp_c+i] = weights[i];
	  	
	  
      if (rvars->num_maps > 1) {
	  	rvars->num_links_map2++;
		if (rvars->num_links_map2 > rvars->max_links_map2) 
			rvars->reset_remap_vars(2,rvars->resize_increment);
		tmp_c = rvars->num_links_map2-1;
		rvars->grid1_add_map2[tmp_c] = add1;
		rvars->grid2_add_map2[tmp_c] = add2;
		tmp_c = tmp_c*num_wts;
        for (i=0; i<num_wts; i++)
			rvars->wts_map2[tmp_c+i] = weights[num_wts+i];    //wts_map2(:,num_links_map2) = weights(num_wts+1:2*num_wts)  ???
      	}

      if (link_add_src[2*add1] == -1) 
	  	link_add_src[2*add1] = rvars->num_links_map1-1;
      if (link_add_dst[2*add2] == -1) 
	  	link_add_dst[2*add2] = rvars->num_links_map1-1;
      
      link_add_src[2*add1+1] = rvars->num_links_map1-1;
      link_add_dst[2*add2+1] = rvars->num_links_map1-1;

      //printf("!!!!!!!!!!\n");  
      //for(tmp_c=0;tmp_c<6;tmp_c++)
      //	    printf("weights  :   %f\n",weights[tmp_c]);
      //printf("end store_link_cnsrv\n");  
}


//adding
void conserv_remap_2D_operator::cal_grid_bound_info(char *grid_name)
{
	int npnts;
	double *center_lons, *center_lats;
	double *vertex_lons, *vertex_lats;
	int nvertexes;
	double *bound_box;
	int nlons, nlats;
	int *bin_addr;
	double *bin_lats, *bin_lons;
	double *tmp_center_lons, *tmp_center_lats;
	double *tmp_vertex_lons, *tmp_vertex_lats;
	int i, j, n, m;
	int  ip1, jp1, n_add, e_add, ne_add;
	double minlat, maxlat, minlon, maxlon;
	double minlat1, maxlat1, minlon1, maxlon1;
	double dlat, dlon;
	double tmp_iv, tmp_iv1;	  
  	double  tmp_lats[4]; 
	double  tmp_lons[4]; 	
       int count;
	int type_bin = LATITUDE;    
	num_srch_bins = 300;                  
       bool is_gcenter = false;

	/* Localize the fields of source or destination grid */
	if (strcmp(grid_name, grid_name_src) == 0) {
		npnts = npts_src;
		center_lons = center_lons_src;
		center_lats = center_lats_src;
		vertex_lons = vertex_lons_src;
		vertex_lats = vertex_lats_src;
		nvertexes = num_vertexes_src;
		nlons = num_lons_src;
		nlats = num_lats_src;
	}
	else {
		npnts = npts_dst;
		center_lons = center_lons_dst;
		center_lats = center_lats_dst;
		vertex_lons = vertex_lons_dst;
		vertex_lats = vertex_lats_dst;
		nvertexes = num_vertexes_dst;
		nlons = num_lons_dst;
		nlats = num_lats_dst;
	}

	tmp_center_lons = new double [npnts];
	tmp_center_lats = new double [npnts];
	tmp_vertex_lons = new double [npnts*nvertexes]; 
	tmp_vertex_lats = new double [npnts*nvertexes]; 

	/* Convert latitude and logitude units from degree to rad */
	if(true) {
		for(int i = 0; i < npnts; i ++) {
			tmp_center_lons[i] = RADIAN(center_lons[i]);
			tmp_center_lats[i] = RADIAN(center_lats[i]);
		}
		for(int i = 0; i < npnts * nvertexes; i++) {
			tmp_vertex_lons[i] = RADIAN(vertex_lons[i]);
			tmp_vertex_lats[i] = RADIAN(vertex_lats[i]);
		}
	}	

	/* Force longitudes between 0 and 2PI and latitudes between -PI/2 and PI/2 */
      	for(i = 0; i < npnts; i ++) {
  		if (tmp_center_lons[i] > PI2)     
			tmp_center_lons[i] = tmp_center_lons[i]  - PI2;
		if (tmp_center_lons[i] < ZERO)  
			tmp_center_lons[i] = tmp_center_lons[i] + PI2;
		if (tmp_center_lats[i] > PIH)     
			tmp_center_lats[i] = PIH;
		if (tmp_center_lats[i] < -PIH)  
			tmp_center_lats[i] = -PIH;
	}  
     	for(i = 0; i < npnts*nvertexes; i ++) {
		if (tmp_vertex_lons[i] > PI2)   
			tmp_vertex_lons[i] = tmp_vertex_lons[i]  - PI2;
		if (tmp_vertex_lons[i] < ZERO)  
			tmp_vertex_lons[i] = tmp_vertex_lons[i] + PI2;
		if (tmp_vertex_lats[i] > PIH)     
			tmp_vertex_lats[i] = PIH;
		if (tmp_vertex_lats[i] < -PIH)   
			tmp_vertex_lats[i] = -PIH;
	}
      

      	bound_box = new double[4*npnts];

	if (! is_gcenter) { 
		for (i = 0; i < npnts; i ++) {
			minlat = MAX_VALUE;
			minlon = MAX_VALUE;
                     maxlat = MIN_VALUE;
			maxlon = MIN_VALUE;
			//printf("%f\t%f\t%f\t%f\n",minlat,minlon,maxlat,maxlon);
			for(j = 0; j < nvertexes; j ++) {
				n = 4 * i + j;
				if (minlat > tmp_vertex_lats[n])  
					minlat = tmp_vertex_lats[n];
				if(maxlat < tmp_vertex_lats[n])  
					maxlat = tmp_vertex_lats[n];
				if(minlon > tmp_vertex_lons[n])  
					minlon = tmp_vertex_lons[n];
				if(maxlon < tmp_vertex_lons[n]) 
					maxlon = tmp_vertex_lons[n];
				}
			bound_box[4*i]   = minlat;
			bound_box[4*i+1] = maxlat;
			bound_box[4*i+2] = minlon;
			bound_box[4*i+3] = maxlon; 
		}
	}
       else {
              for( n=0; n<npnts; n++) {
			j = n/nlons +1;
                	i = (n+1) - (j-1)*nlons;
			if (i < nlons) 
				ip1 = i + 1;
                      else {
				//assume cyclic
                            ip1 = 1;
				//but if it is not, correct
                            e_add = (j - 1)*nlons + ip1;
				if (fabs(tmp_center_lats[e_add] -tmp_center_lats[n]) > PIH)
						ip1 = i;
				}

			if (j < nlats)
				jp1 = j+1;
			else  {
				//assume cyclic           ??????????????????????? 
				jp1 = 1;
                            //but if it is not, correct
                            n_add = (jp1 - 1)*nlons + i;
                            if (fabs(tmp_center_lats[n_add] -tmp_center_lats[n]) > PIH)
					       jp1 = j;
				}
				
			n_add = (jp1 - 1)*nlons + i;
                     e_add = (j - 1)*nlons + ip1;
                     ne_add = (jp1 - 1)*nlons + ip1;

                     //find N,S and NE lat/lon coords and check bounding box
                     tmp_lats[0] = tmp_center_lats[n];
                     tmp_lats[1] = tmp_center_lats[e_add];
                     tmp_lats[2] = tmp_center_lats[ne_add];
                     tmp_lats[3] = tmp_center_lats[n_add];
							
                     tmp_lons[0] = tmp_center_lons[n];
                     tmp_lons[1] = tmp_center_lons[e_add];
                     tmp_lons[2] = tmp_center_lons[ne_add];
                     tmp_lons[3] = tmp_center_lons[n_add];

                     if(tmp_lats[0]<tmp_lats[1]) {
				minlat  = tmp_lats[0];
				maxlat = tmp_lats[1];
			}
			else { 
				maxlat = tmp_lats[0];
				minlat  = tmp_lats[1];
			}

			if(tmp_lats[2]<tmp_lats[3]) {
				minlat1=tmp_lats[2];
				maxlat1=tmp_lats[3];
			}
			else {
				maxlat1=tmp_lats[2];
				minlat1=tmp_lats[3];
			}
			if(minlat > minlat1) 
				minlat = minlat1;
			if(maxlat < maxlat1)
			 	maxlat = maxlat1;
			if(tmp_lons[0]<tmp_lons[1]) {
				minlon  = tmp_lons[0];
				maxlon = tmp_lons[1];
			}
			else { 
				maxlon = tmp_lons[0];
				minlon  = tmp_lons[1];
			}

			if(tmp_lons[2]<tmp_lons[3]) {
				minlon1=tmp_lons[2];
				maxlon1=tmp_lons[3];
			}
			else {
				maxlon1=tmp_lons[2];
				minlon1=tmp_lons[3];
			}
			if(minlon > minlon1) 
			 	minlon = minlon1;
			if(maxlon < maxlon1)
			 	maxlon = maxlon1;
    
			bound_box[n*4]     = minlat;
                     bound_box[n*4+1] = maxlat;
                     bound_box[n*4+2] = minlon;
                     bound_box[n*4+3] = maxlon;
		}
       }
	for(i=0; i<npnts; i++)
	   	if(fabs(bound_box[4*i+3] - bound_box[4*i+2]) > PI) {
			bound_box[4*i+2]  = ZERO;
			bound_box[4*i+3] = PI2;
		}

	for(i = 0; i < npnts; i ++) {
              //try to check for cells that overlap poles
              if(tmp_center_lats[i] > bound_box[4*i+1])
			bound_box[4*i+1] = PIH;
			 
		if(tmp_center_lats[i] < bound_box[4*i])
		   	bound_box[4*i] = -PIH;
	}
	if (type_bin==LATITUDE) {  // Using latitude bins to restrict search

		bin_addr = new int[2*num_srch_bins];
		bin_lats = new double[2*num_srch_bins];
		bin_lons = new double[2*num_srch_bins]; 

        	dlat = PI/num_srch_bins;

	       //printf("dlat:  %lf\n",    dlat);
		//printf("num_srch_bins:  %d\n",num_srch_bins);		

              for(i=0; i<num_srch_bins; i++) {
                     j=i*2;                                        //???
			bin_lats[j] = i*dlat - PIH;
                     bin_lats[j+1] = (i+1)*dlat - PIH;
                     bin_lons[j] = ZERO;
                     bin_lons[j+1] = PI2;
			bin_addr[j] = npnts;       //bin_addr[j] = npnts *4;
                     bin_addr[j+1] = -1;
		}

              count = 0;
	       for(i = 0; i < npnts; i ++) {
		   	n = i * 4;
		   	for(j = 0; j < num_srch_bins; j ++) {
                            m = j * 2;
                            if (bound_box[n] <= bin_lats[m+1] && bound_box[n+1] >= bin_lats[m]) {
				       bin_addr[m]     =  MIN(i,  bin_addr[m]);
	                            bin_addr[m+1] =  MAX(i, bin_addr[m+1]);
				}
			}	  
		}
	}

	else if (type_bin==LATLON)                // test not cover this piece of code
	{  // Using lat/lon boxes to restrict search
		dlat = PI/num_srch_bins;
              dlon = PI2/num_srch_bins;

		bin_addr = new int[2*num_srch_bins*num_srch_bins];
		bin_lats  = new double[2*num_srch_bins*num_srch_bins];
		bin_lons = new double[2*num_srch_bins*num_srch_bins]; 

              n = 0;
		tmp_iv = 0;
		tmp_iv1 = 0;
              for (i=0; i<npnts; i++) {
			tmp_iv = i*dlat - PIH;
			tmp_iv1 = (i+1)*dlat - PIH;
			for (j=0; j<npnts; j++)
			{
				bin_lats[n] = tmp_iv;
                            bin_lats[n+1] = tmp_iv1;
                            bin_lons[n] = j*dlon;
                            bin_lons[n+1] = (j+1)*dlon;
                            bin_addr[j] = npnts;       //bin_addr[n] = npnts *4;
                            bin_addr[n+1] = -1;
				n=n+2;
			}
                 }
                 num_srch_bins = num_srch_bins*num_srch_bins;
					
              for(i=0; i<npnts; i++) 
		{
			n = i*4;
			for(j=0; j<num_srch_bins; j++) 
			{
				m = j*2;
				if ( bound_box[n] <= bin_lats[m+1] && bound_box[n+1] >= bin_lats[m] && bound_box[n+2]  <= bin_lons[m+1] && bound_box[n+3] >= bin_lons[m]) 
				{
					bin_addr[m]     =  MIN(i, bin_addr[m]);
	                            bin_addr[m+1] =  MAX(i,bin_addr[m+1]);
					//bin_addr[m]     =  MIN(n, bin_addr[m]);
	                            //bin_addr[m+1] =  MAX(n,bin_addr[m+1]);   
				}
			}
                }	  
	}
	else 
	{
	   	printf("unknown search restriction method\n");
		exit(1);
	}

	if (strcmp(grid_name, grid_name_src) == 0) {
		center_lons_src = tmp_center_lons;
		center_lats_src = tmp_center_lats;
		vertex_lons_src = tmp_vertex_lons;
		vertex_lats_src = tmp_vertex_lats;
		bin_addr_src = bin_addr;
		bound_box_src = bound_box;
	}
	else {
		center_lons_dst = tmp_center_lons;
		center_lats_dst = tmp_center_lats;
		vertex_lons_dst = tmp_vertex_lons;
		vertex_lats_dst = tmp_vertex_lats;
		bin_addr_dst = bin_addr;
		bound_box_dst = bound_box;
	}

	delete [] bin_lats;
	delete [] bin_lons;
} 


