#ifndef data_bundle_H
#define data_bundle_H

#include "data_cfg.h"
#include "grid.h"

class data_bundle
{
	private:
		int num_fields;			// The number of fields in the object
		char **name_fields;		// The name of each field
		char *data_buf;			// The buffer of hold the data of all fields
		char **field_data;		// The start position of each field in the data buffer
		int *type_fields;			// The data type of each field
		int *length_fields;		// The number of data in each field
		int *counts;
	public:
		data_bundle(data_bundle_cfg *, int, int);
		char *get_field_data(int i) {return field_data[i]; }
		char *get_field_name(int i)	{return name_fields[i]; }
		~data_bundle();
		int copy_in_data(char *);
		int copy_out_data(char *);
};

#endif
