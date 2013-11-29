all : life xlife

life : life.c
	mpicc -o life life.c

xlife : xlife.c
	mpicc -o xlife xlife.c
