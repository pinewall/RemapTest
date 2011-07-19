icpc -O3 -g -I/opt/netCDF/include *.cxx -c
icpc -o thu-remap *.o -L/opt/netCDF/lib/ -lnetcdf

