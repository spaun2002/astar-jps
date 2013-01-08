#include "AStar.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <vector>

int main (int argc, char **argv)
{
	if (argc != 2) {
		fprintf (stderr, "testAStar <scenfile>\n");
		fprintf (stderr, "(where <scenfile> is of the format used in http://www.aiide.org/benchmarks/)\n");
		exit (1);
	}

	FILE *scenFile = fopen (argv[1], "r");
	if (!scenFile) {
		perror ("couldn't open given scenario file: ");
		exit (1);
	}

	fscanf (scenFile, "version 1.0\n");
	char mapFileBuf[255];
	int bucket, height, width, startX, startY, goalX, goalY, optimal;
	double something;
	fscanf (scenFile, "%i %s %i %i %i %i %i %i %i %lf\n", 
		&bucket, mapFileBuf, &width, &height, &startX, &startY,
		&goalX, &goalY, &optimal, &something);

	FILE *mapFile = fopen (mapFileBuf, "r");
	if (!mapFile) {
		fprintf (stderr, "couldn't open map file %s: %s\n", 
			 mapFileBuf, strerror(errno));
		exit (1);
	}

	fscanf (mapFile, "type octile\nheight %i\nwidth %i\nmap\n", &height, &width);

  std::vector<char> grid(width*height);
	memset (grid.data(), 0, sizeof (char) * width * height);

  std::vector<char> buf(width+1); //< space for a new line

	for (int i = 0; i < height; i++) {
		fread (buf.data(), 1, width + 1, mapFile);
		for (int j = 0; j < width; j++) {
			if (buf[j] == '.' || buf[j] == 'G')
				grid[width*i+j] = 0;
			else
				grid[width*i+j] = 1;
		}
	}
	int doContinue = 1;
	do {
		int solLen = 0;

    AStarJPS::SimpleMap<char> map(grid.data(), width, height);

    auto result = AStarJPS::AStarCompute(map, solLen, AStarJPS::Coordinate(startX, startY), AStarJPS::Coordinate(goalX, goalY));

		if (solLen > optimal) {
			fprintf (stderr, "validity error! In map %s, from (%i,%i) to (%i, %i), expected length %i, was length %i\n", mapFileBuf, startX, startY, goalX, goalY, optimal, solLen);
			exit (1);
		}
		doContinue = fscanf(scenFile,"%i %s %i %i %i %i %i %i %i %lf\n",
				     &bucket, mapFileBuf, &width, &height, 
				     &startX, &startY, &goalX,
				     &goalY, &optimal, &something);
	} while (doContinue > 0);
}

