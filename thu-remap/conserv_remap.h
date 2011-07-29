
#define CONSERV_REMAP_H

#include "remap_base.h"
#include "remap_vars.h"

class conserv_remap_2D_operator : public Common_remap
{
       //variables (save) in intersection
       //private:
    protected:

		int last_loc;                       // save location when crossing threshold
		bool lthresh;          		// flags segments crossing threshold bndy  !also in pole_intersection?е?intent(inout)
		double intrsct_lat_off;          // lat coords offset for next search
		double intrsct_lon_off;         // lon coords offset for next search

		//variables (save) in pole_intersection

         	// save last intersection to avoid roundoff during coord transformation
          	bool luse_last;
          	double intrsct_x;                 // x,y for intersection
          	double intrsct_y;

          	// variables necessary if segment manages to hit pole
          	int avoid_pole_count;         // count attempts to avoid pole
          	double avoid_pole_offset;  // endpoint offset to avoid pole

          	 // variables (save) in store_link_cnsrv
          	bool first_call;

          	// module variables
          	int num_srch_cells;   // num cells in restricted search arrays     ?? save
		int * srch_add;       // global address of cells in srch arrays        ?? save

		double *srch_corner_lat;    // lat of each corner of srch cells     ?? save
		int srch_corner_lat_dim1;
		int srch_corner_lat_dim2;	   
		double *srch_corner_lon;   // lon of each corner of srch cells     ?? save    


		//interface 
		int num_srch_bins;
		double *bound_box_src;
		double *bound_box_dst;
		int *bin_addr_src;
		int *bin_addr_dst;

		double *grid1_area_in;                        // ??????????? no initialization 
		double *grid2_area_in;                        // ??????????? no test case

		int wgt_size;                                         // the size of local weight array(size is 6 currently)
		int *link_add_src;                                     // min,max link add to restrict search
		int *link_add_dst;                                     // min,max link add to restrict search
         
		remap_vars *rvars; 



		//for test
		int pole_count;
		int inter_count;

		int store_count1;
		int store_count2;
		int store_count3;
		int store_count4;
		int store_count_return1;
		int store_count_return2;
		
		int store_count_all;

	public:
		conserv_remap_2D_operator(char *, char *,  char *);
		~conserv_remap_2D_operator();
		void cal_grid_bound_info(char *);

		void conserv_remap();
		
		bool intersection(int*, double*, double*, double, double, double, double, double*, bool, bool);
		void line_integral(double*, int, double, double, double, double, double, double, double, double);    
              bool pole_intersection(int*, double*, double*, double, double, double, double, double*, bool);
	       void store_link_cnsrv(int, int, double*);

		void init_remap(); 
		void cal_remap(double *, double *);
		
};

class conserv2_remap_2D_operator : public conserv_remap_2D_operator
{
    public:
        conserv2_remap_2D_operator(char *, char *, char *);
        ~conserv2_remap_2D_operator();
        void cal_remap(double *, double *);
};
