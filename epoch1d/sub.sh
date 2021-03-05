#########################################################################
# File Name: sub.sh
# Author: Xieyu
# mail: xieyu123@pku.edu.cn
# Created Time: 2020年07月11日 星期六 13时01分05秒
#########################################################################
#!/bin/bash
echo Data/1D/be7/ | mpirun -np 16 ./bin/epoch1d 
echo Data/1D/be75/ | mpirun -np 16 ./bin/epoch1d 
echo Data/1D/be8/ | mpirun -np 16 ./bin/epoch1d 
echo Data/1D/be85/ | mpirun -np 16 ./bin/epoch1d 
echo Data/1D/be9/ | mpirun -np 16 ./bin/epoch1d 

