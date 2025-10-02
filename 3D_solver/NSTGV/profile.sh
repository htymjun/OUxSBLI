nohup nsys profile -t cuda,nvtx,openacc,osrt -f true -o my_report mpiexec -np 2 ./a.out &

