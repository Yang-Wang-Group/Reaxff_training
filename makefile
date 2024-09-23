CC=gcc

OBJS =  init.o reaxffio.o deviance.o mcmosa.o output.o main.o pso.o
OPT = -c -Wall -g -O 

autofit:$(OBJS) 
	$(CC) $(OBJS) -g -o autofit -lm 
main.o:main.c
	$(CC) main.c $(OPT)
init.o:init.c
	$(CC) init.c $(OPT)
reaxffio.o:reaxffio.c
	$(CC) reaxffio.c $(OPT)
deviance.o:deviance.c
	$(CC) deviance.c $(OPT)
mcmosa.o:mcmosa.c
	$(CC) mcmosa.c $(OPT)
output.o:output.c
	$(CC) output.c $(OPT)
pso.o:pso.c
	$(CC) pso.c $(OPT)
clean:
	rm *.o
