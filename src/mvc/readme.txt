Instructions to use OpenMP with Matlab (Linux specific)
1. Set the following variables BEFORE starting Matlab
	export LD_PRELOAD=/usr/lib64/libgomp.so.1.0.0
	export OMP_NUM_THREADS=12
2. Start matlab
3. Compile mex code
4. Copy the generated .mexa64/.mexmaci64/.mexw64 file to src folder (parent)
