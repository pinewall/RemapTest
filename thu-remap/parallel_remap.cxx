#include "parallel_remap.h"
#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "grid_func.h"

/* This function builds parallel remap weight information according to the 
    decomposition information and the sequential remap weight information */
void Parallel_remap::build_parallel_remap(int *decomp_indx_src, 
                                                                int *decomp_indx_dst,
                                                                int grid_size_src,
                                                                int grid_size_dst,
                                                                int num_wgts_seq,
                                                                int *wgt_indx_src_seq,
                                                                int *wgt_indx_dst_seq,
                                                                double *wgt_values_seq)
{
    int i, j;
    int nwgt_per_dst_cell;
    int *decomp_hash_table_src, *decomp_hash_table_dst;


    /* Initialize the hash tables for looking up the indexes of cells in 
         the decomposition of source and destination grids */
    decomp_hash_table_src = new int [grid_size_src];
    decomp_hash_table_dst = new int [grid_size_dst];
    for (i = 0; i < grid_size_src; i ++)
        decomp_hash_table_src[i] = -1;
    for (i = 0; i < grid_size_dst; i ++)
        decomp_hash_table_dst[i] = -1;
    for (i = 0; i < npts_src; i ++) 
        decomp_hash_table_src[decomp_indx_src[i]] = i;
    for (i = 0; i < npts_dst; i ++) 
        decomp_hash_table_dst[decomp_indx_dst[i]] = i;

    /* Compute the number of weights in the parallel remap */
    if (wgt_indx_dst_seq == NULL) {
        nwgt_per_dst_cell = num_wgts_seq / grid_size_dst;
        num_wgts = npts_dst * nwgt_per_dst_cell;
    }
    else {
        num_wgts = 0;
        for (i = 0; i < num_wgts_seq; i ++)
            if (decomp_hash_table_dst[wgt_indx_dst_seq[i]] != -1)
                num_wgts ++;
    }

    /* Allocate the memory space for parallel remap information */
    if (wgt_indx_dst_seq == NULL)
        wgt_indx_dst = NULL;
    else
        wgt_indx_dst = new int [num_wgts];
    wgt_indx_src = new int [num_wgts];
    wgt_values = new double [num_wgts];

    /* Compute the weight information arrays for parallel remap */
    if (wgt_indx_dst_seq == NULL) {
        for (i = 0; i < npts_dst; i ++) {
            for (j = 0; j < nwgt_per_dst_cell; j ++) {
                wgt_indx_src[i*nwgt_per_dst_cell+j] = wgt_indx_src_seq[decomp_indx_src[i]*nwgt_per_dst_cell+j];
                wgt_values[i*nwgt_per_dst_cell+j] = wgt_values_seq[decomp_indx_src[i]*nwgt_per_dst_cell+j];
            }
        }
    }
    else {
        num_wgts = 0;
        for (i = 0; i < num_wgts_seq; i ++)
            if (decomp_hash_table_dst[wgt_indx_dst_seq[i]] != -1) {
                wgt_indx_src[num_wgts] = wgt_indx_src_seq[i];
                wgt_indx_dst[num_wgts] = wgt_indx_dst_seq[i];
                wgt_values[num_wgts] = wgt_values_seq[i];
                num_wgts ++;
            }
    }

    delete [] decomp_hash_table_src;
    delete [] decomp_hash_table_dst;    
}



Parallel_remap::Parallel_remap(char *decomp1_name, 
                                              char *decomp2_name, 
                                              char *alg_name)
: Common_remap(decomp1_name, decomp2_name, alg_name)
{
    /* Get the index array of the decomposition for source and desitination grids */

    /* Get the weight information of the corresponding sequential remap algorithm */

    /* Call the function to build parallel remap information */
    printf("ok\n");
}
