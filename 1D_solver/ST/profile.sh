timestamp=$(date +"%Y%m%d_%H%M%S")

nsys profile \
    -t cuda,nvtx,openacc,osrt \
    -f true \
    -o "my_report_${timestamp}" \
    mpiexec -np 1 ./build/a.out &

ncu \
    --set full \
    --target-processes all \
    --import-source yes \
    -o "my_report_ncu_${timestamp}" \
    mpiexec -np 1 ./build/a.out &
