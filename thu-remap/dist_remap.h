#ifndef DIST_REMAP_H
#define DIST_REMAP_H

#include "remap_base.h"

class dist_remap_2D_operator : public Common_remap
{
	private:
		int count;					// Number of useful elements in two arrays
		double dist_threshold;
		double power_p;
		int num_nearest_pnts;		// Number of nearest points around a point. This value is specified by user
		
	public:
		dist_remap_2D_operator(char*, char*, char*, int, double, double);
		~dist_remap_2D_operator();
		void search_neighbors_all_pnts_dst(double*, double*, double*, double*, int*);
		void init_remap(); 
		void cal_remap(double*, double *);
};

#endif
