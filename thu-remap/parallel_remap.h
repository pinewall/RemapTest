#ifndef PARALLEL_REMAP
#define PARALLEL_REMAP

#include "remap_base.h"


class Parallel_remap : public Common_remap
{

		
	public:
		Parallel_remap(char *, char *, char *);
		~Parallel_remap(){}
		void build_parallel_remap(int*, int*, int, int, int, int*, int*, double*);
		void init_remap(){}
		void cal_remap(double *, double *){}
};


#endif