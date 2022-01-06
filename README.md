# LSOPs-Benchmark-2022

Peilan Xu, Wenjian Luo, Xin Lin, and Zeneng She, "A Large-Scale Continuous Optimization Benchmark Suite with Versatile Coupled Heterogeneous Modules"

Version: 1.0

Developers: Peilan Xu

## Core files:
benchmark_func2022.m: Main program to compute the fitness values for test functions.

base_functions.m: The base functions used in the benchmark suite. 

## Call

fitness = feval('benchmark_func2022', population, func_num);

population: the matrix with N rows and D columns, where N is the size of the population, and D is the dimension of the problem;

func_num: the index of the test function.

## Contact:
Questions and bug reports should be sent to:

Peilan Xu: xpl@mail.ustc.edu.cn

