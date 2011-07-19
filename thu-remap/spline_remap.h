#ifndef SPLINE_REMAP_H
#define SPLINE_REMAP_H


class spline_remap_2D_operator : public Common_remap
{
	private:
		int remap_directions[2]; //0:  no remap; 1: remap at lat direction; 2: remap at lon direction
	public:
		spline_remap_2D_operator(char*, char*, char*, int);
		
};

#endif
