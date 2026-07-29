nsys profile -t cuda,nvtx,openacc,osrt -f true -o my_report mpiexec -np 1 ./build/a.out
