#include "data_bundle.h"
#include "name_cfg.h"
#include "models_cfg.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

data_bundle::data_bundle(data_bundle_cfg *data_cfg, 
                                int num_data_in_grid_2D, 
                                int num_data_in_grid_3D)
{
	int buf_size;
	int i;

	num_fields = data_cfg->num_fields_2D + data_cfg->num_fields_3D;
	name_fields = new char * [num_fields];
	type_fields = new int [num_fields];
	field_data = new char* [num_fields];
	length_fields = new int [num_fields];
	counts = new int [num_fields];
	
	for (i = 0; i < num_fields; i ++) {
		name_fields[i] = new char [NAME_STR_SIZE];
		strcpy(name_fields[i], data_cfg->all_fields[i].name_field);
		type_fields[i] = data_cfg->all_fields[i].type_field;
		counts[i] = 0;
	}

	for (i = 0; i < data_cfg->num_fields_2D; i ++)
		length_fields[i] = num_data_in_grid_2D;
	for (; i < num_fields; i ++)
		length_fields[i] = num_data_in_grid_3D;

	for (i = 0, buf_size = 0; i < num_fields; i ++)
		if (type_fields[i] == 0)
			buf_size += length_fields[i] * sizeof(int);
		else if (type_fields[i] == 1)
			buf_size += length_fields[i] * sizeof(float);
		else if (type_fields[i] == 2)
			buf_size += length_fields[i] * sizeof(double);
#ifdef DEBUG_data_bundle
		else printf("error in create data fields\n");
#endif

	data_buf = new char [buf_size];

	printf("buf size is %d\n", buf_size);

	for (i = 0, buf_size = 0; i < num_fields; i ++) {
		field_data[i] = data_buf + buf_size;
		if (type_fields[i] == 0)
			buf_size += length_fields[i] * sizeof(int);
		else if (type_fields[i] == 1)
			buf_size += length_fields[i] * sizeof(float);
		else if (type_fields[i] == 2)
			buf_size += length_fields[i] * sizeof(double);
	}
#ifdef DEBUG_data_bundle
		else printf("error in create data fields\n");
#endif
}

int data_bundle::copy_in_data(char *input_data)
{
	int i, whole_size = 0;

	for (i = 0; i < num_fields; i ++)
		whole_size += length_fields[i] * sizeof(double);
	for (i = 0; i < whole_size; i ++)
		data_buf[i] = input_data[i];

	return whole_size;
}

int data_bundle::copy_out_data(char *output_data)
{
	int i, whole_size = 0;

	for (i = 0; i < num_fields; i ++)
		whole_size += length_fields[i] * sizeof(double);
	for (i = 0; i < whole_size; i ++)
		output_data[i] = data_buf[i];

	return whole_size;
}

data_bundle::~data_bundle()
{
	for (int i = 0; i < num_fields; i ++)
		delete [] name_fields[i];
	
	delete [] name_fields;
	delete [] data_buf;
	delete [] field_data;
	delete [] type_fields;
	delete [] length_fields;
	delete [] counts;
}