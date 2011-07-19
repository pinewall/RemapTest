#ifndef BILINEAR_REMAP_H
#define BILINEAR_REMAP_H

#include "remap_base.h"
#include "remap_func.h"
class bilinear_remap_2D_operator : public Common_remap
{	
	private:		
		int *nlon_each_lat_src;			//
		int *nlon_each_lat_dst;			//
		int *lat_begindx_src;				//
		int *lat_begindx_dst;				//

	public:
		bilinear_remap_2D_operator(char*, char*, char*);
		~bilinear_remap_2D_operator();
		void wgt_bilinear_regular(double*, double*, double, double, double*);
		bool wgt_bilinear_irregular(double*, double*, double, double, double*);
		void init_remap(); 
		void cal_remap(double *, double *);
		bool point_in_polygon(double, double, int, double*, double*);
};

#endif
