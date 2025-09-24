#!/bin/bash
export CC=/usr/bin/gcc
gcccmd="$CC -std=gnu23 -Wp,-D_GNU_SOURCE,-DNODEBUG,-U_DEFAULT_SOURCE,-DVER=1.0 -Ofast -flto -funroll-loops -fno-asynchronous-unwind-tables -fomit-frame-pointer -fpeel-loops -funswitch-loops -fno-ident -fmerge-all-constants -falign-functions=64 -falign-loops=64 -falign-jumps=64 -fwhole-program -fipa-pta -fipa-cp-clone -fgcse-sm -fgcse-las -fgcse-after-reload -fivopts -ftree-loop-distribution -ftree-loop-vectorize -ftree-slp-vectorize -fvect-cost-model=unlimited -fassociative-math -freciprocal-math -fno-signed-zeros -fno-trapping-math -fno-math-errno -fno-stack-protector -fno-plt -fopenmp -march=native -mtune=native -Wl,-O3,--as-needed,--strip-all,--gc-sections,--relax,--sort-common picalc.c -lgmp -o picalc"
echo "${gcccmd}"
${gcccmd}
exit 0
