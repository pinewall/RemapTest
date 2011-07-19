#ifndef REMAP_CFG_H

#include "name_cfg.h"

#define REMAP_CFG_H
#define REMAP_PARA_NUM	2


struct remap_cfg
{
	char name_data_field[NAME_STR_SIZE];
	char name_grid_src[NAME_STR_SIZE];
	char name_grid_dst[NAME_STR_SIZE];
	char name_remap_alg[NAME_STR_SIZE];
	int i_para1;
	int i_para2;
	double d_para1;
	double d_para2;
};

extern int num_data_remaps;
extern remap_cfg data_remaps_cfg[];
#endif