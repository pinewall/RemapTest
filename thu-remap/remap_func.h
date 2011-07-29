#ifndef REMAP_FUNC_H
#define REMAP_FUNC_H

#include "remap_base.h"

#define NAME_DIST_REMAP_2D		    "dist_remap_2D"
#define NAME_BILINEAR_REMAP_2D		"bilinear_remap_2D"
#define NAME_BICUBIC_REMAP_2D		"bicubic_remap_2D"
#define NAME_SPLINE_REMAP_2D		"spline_remap_2D"
#define NAME_POLYN_REMAP_2D		    "polyn_remap_2D"
#define NAME_CONSERV_REMAP_2D		"conserv_remap_2D"
#define NAME_CONSERV2_REMAP_2D		"conserv2_remap_2D"


class remap_operator_group
{
	private:
		int count;
		Common_remap **remap_operators;
		char **name_operators;

	public:
		remap_operator_group();
		~remap_operator_group();
		void generate_name(char*, int*, double*);
		void generate_operators();
		int find_operator(char *, char *, char *);
		void do_remap(char *, char *, char *, double *, double *);
};

extern void generate_remap_operators();
extern void do_remap(char*, char*, char*, double*, double*);
extern void finalize_remap_operators();
extern void remap_model_set(char*, char*);

#endif
