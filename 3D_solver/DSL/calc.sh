rm ./data/*.d
rm ./data/*.vtr
nohup mpiexec -n 2 ./a.out &
cp mod_globals.f90 ./data
cp set.f90 ./data

