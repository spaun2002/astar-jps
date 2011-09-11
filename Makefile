testAStar: AStar.o IndexPriorityQueue.o TestAStar.o
	gcc -g -lm  IndexPriorityQueue.o AStar.o TestAStar.o -o testAStar

AStar.o: AStar.c AStar.h IndexPriorityQueue.h
	gcc -march=native -O2 -Wall -W -std=c99 AStar.c -c -o AStar.o

TestAStar.o: TestAStar.c AStar.h
	gcc -march=native -O2 -Wall -W -std=c99 TestAStar.c -c -o TestAStar.o

IndexPriorityQueue.o: IndexPriorityQueue.c IndexPriorityQueue.h
	gcc -march=native -O2 -Wall -W -std=c99 IndexPriorityQueue.c -c -o IndexPriorityQueue.o