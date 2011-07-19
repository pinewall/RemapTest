#ifndef GRID_FUNC_H
#define GRID_FUNC_H

#include "grid.h"

#define NAME_CPL_GRID_2D	"rect_grid_ocn_360x196"			//modified to "cpl_grid_2d" in future
#define NAME_CPL_GRID_3D	"cpl_grid_3d"

class grid_group
{
	private: 
		int size;
		int count;
		grid **grids;

	public:
		grid_group();
		~grid_group();
		int find_grid(char*);
		void generate_grid(char*, char*, int, int);
		void realloc();
		grid *get_grid(char*);
};


extern grid_group *cpl_grids;
extern void generate_all_cpl_grids();
extern void finalize_all_cpl_grids();

#endif