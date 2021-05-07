mex -c -compatibleArrayDims mod_main.f90
mex -compatibleArrayDims sus_iter.f90 mod_main.o
mex -compatibleArrayDims comm_iter.f90 mod_main.o