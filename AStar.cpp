#include "AStar.hpp"
#include "IndexPriorityQueue.hpp"
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <queue>

// Distance metrics, you might want to change these to match your game mechanics

// Chebyshev distance metric for distance estimation by default
inline static double estimateDistance (Coordinate start, Coordinate end)
{
	return std::max(abs(start.x - end.x), 
                  abs(start.y - end.y));
}

// Since we only work on uniform-cost maps, this function only needs
// to see the coordinates, not the map itself.
// Euclidean geometry by default.
// Note that since we jump over points, we actually have to compute 
// the entire distance - despite the uniform cost we can't just collapse
// all costs to 1
static double preciseDistance (Coordinate start, Coordinate end)
{
  using namespace std;
  auto x = abs(start.x - end.x);
  auto y = abs(start.y - end.y);
  return sqrt(x*x + y*y);
}

// Below this point, not a lot that there should be much need to change!

//struct AStar {
//	const char *grid;
//	Coordinate bounds;
//	node start;
//	node goal;
//	queue *open;
//	char *closed;
//	double *gScores;
//	node *cameFrom;
//	int *solutionLength;
//};

// The order of directions is: 
// N, NE, E, SE, S, SW, W, NW 
//typedef unsigned char direction;
//#define NO_DIRECTION 8
//typedef unsigned char directionset;

enum class Direction
{
  N = 0,
  NE,
  E,
  SE,
  S, 
  SW, 
  W, 
  NW,
  NO_DIRECTION = 8
};

typedef std::set<Direction> DirectionSet;

typedef int node;

class AStar_jsp
{
public:
  const char* m_grid;
  Coordinate m_bounds;
  node m_start;
  node m_goal;
  queue* m_open;
  char *m_closed;
  double *m_gScores;
  node *m_cameFrom;
  int *m_solutionLength;

public:
  bool isEnterable (const Coordinate& coord);
  DirectionSet forcedNeighbours (const Coordinate& coord, Direction dir);
  void addToOpenSet (int node, int nodeFrom);
  bool HasForcedNeighbours(Coordinate coord, Direction dir);
  int jump(Direction dir, int start);
  int nextNodeInSolution (int *target, int node);
  int *recordSolution ();
  Direction directionWeCameFrom (int node, int nodeFrom);
};

// return and remove a direction from the set
// returns NO_DIRECTION if the set was empty
static Direction nextDirectionInSet (DirectionSet& dirs)
{
  if (dirs.empty())
    return Direction::NO_DIRECTION;

  auto result = *dirs.begin();
  dirs.erase(dirs.begin());

  return result;
}

static void addDirectionToSet (DirectionSet& dirs, Direction dir)
{
  dirs.insert(dir);
}

static void addDirectionsToSet(DirectionSet& dirs, const DirectionSet& newDirections)
{
  dirs.insert(begin(newDirections), end(newDirections));
}

/* Coordinates are represented either as pairs of an x-coordinate and
   y-coordinate, or map indexes, as appropriate. getIndex and getCoord
   convert between the representations. */
static int getIndex (Coordinate bounds, Coordinate c)
{
	return c.x + c.y * bounds.x;
}

int astar_getIndexByWidth (int width, int x, int y)
{
	return x + y * width;
}

static Coordinate getCoord (Coordinate bounds, int c)
{
	return Coordinate ( c % bounds.x, c / bounds.x );
}

void astar_getCoordByWidth (int width, int node, int *x, int *y)
{
	*x = node % width;
	*y = node / width;
}


// is this coordinate contained within the map bounds?
static int contained (Coordinate bounds, Coordinate c)
{
	return c.x >= 0 && c.y >= 0 && c.x < bounds.x && c.y < bounds.y;
}

// is this coordinate within the map bounds, and also walkable?
bool AStar_jsp::isEnterable (const Coordinate& coord)
{
	node node = getIndex (this->m_bounds, coord);
	return contained (this->m_bounds, coord) && 
		this->m_grid[node];
}

static int directionIsDiagonal (Direction dir)
{
  //return (dir == Direction::NE || dir == Direction::NW || dir == Direction::SE || dir == Direction::SW);
  return (static_cast<int>(dir) % 2) != 0;
}

// the coordinate one tile in the given direction
static Coordinate adjustInDirection (Coordinate c, Direction dir)
{
	// we want to implement "rotation" - that is, for instance, we can
	// subtract 2 from the direction "north" and get "east"
	// C's modulo operator doesn't quite behave the right way to do this,
	// but for our purposes this kluge should be good enough
  switch (dir)
  {
  case Direction::N:  return Coordinate(c.x,     c.y - 1);
  case Direction::NE: return Coordinate(c.x + 1, c.y - 1);
  case Direction::E:  return Coordinate(c.x + 1, c.y    );
  case Direction::SE: return Coordinate(c.x + 1, c.y + 1);
  case Direction::S:  return Coordinate(c.x,     c.y + 1);
  case Direction::SW: return Coordinate(c.x - 1, c.y + 1);
  case Direction::W:  return Coordinate(c.x - 1, c.y    );
  case Direction::NW: return Coordinate(c.x - 1, c.y - 1);
  default: 
    return Coordinate( -1, -1 );
  }
}

//
// rotationFactor - numbers from -7 to 7 including
//
static Direction RotateDirection(Direction dir, int rotationFactor)
{
  // 7 0 1
  // 6 x 2
  // 5 4 3
  int tmp = static_cast<int>(dir);
  tmp += rotationFactor + 65536;
  tmp %= static_cast<int>(Direction::NO_DIRECTION);
  return static_cast<Direction>(tmp);
}

// logical implication operator
static int implies (int a, int b)
{
	return a ? b : 1;	
}

/* Harabor's explanation of exactly how to determine when a cell has forced
   neighbours is a bit unclear IMO, but this is the best explanation I could
   figure out. I won't go through everything in the paper, just the extra
   insights above what I thought was immediately understandable that it took
   to actually implement this function.

   First, to introduce the problem, we're looking at the immedate neighbours
   of a cell on the grid, considering what tile we arrived from.

   ...  This is the basic situation we're looking at. Supposing the top left
   -X.  period is cell (0,0), we're coming in to cell (1, 1) from (0, 1).
   ...  

   ...  The other basic case, the diagonal case. All other cases are obviously
   .X.  derivable from these two cases by symmetry.
   /..

   The question is: Given that some tiles might have walls, *how many tiles
   are there that we can reach better by going through the center tile than
   any other way?* (for the horizontal case it's ok to only be able to reach
   them as well some other as through the center tile too)

   In case there are no obstructions, the answers are simple: In the horizontal
   or vertical case, the cell directly ahead; in the diagonal case, the three
   cells ahead.

   The paper is pretty abstract about what happens when there *are* 
   obstructions, but fortunately the abstraction seems to collapse into some
   fairly simple practical cases:

   123  Position 4 is a natural neighbour (according to the paper's terminology)
   -X4  so we don't need to count it. Positions 1, 2, 5 and 6 are accessible
   567  without going through the center tile. This leaves positions 3 and 7
   to be looked at.

   Considering position 3 (everything here also follows for 7 by symmetry):
   If 3 is obstructed, then it doesn't matter what's in position in 2.
   If 3 is free and 2 is obstructed, 3 is a forced neighbour.
   If 3 is free and 2 is free, 3 is pruned (not a forced neighbour)

   i.e. logically, 
   3 is not a forced neighbour iff (3 is obstructed) implies (2 is obstructed).

   Similar reasoning applies for the diagonal case, except with bigger angles.
   
 */
/*
static int hasForcedNeighbours (astar_t *astar, Coordinate coord, int dir)
{
#define ENTERABLE(n) isEnterable (astar, \
	                          adjustInDirection (coord, dir + (n)))
	if (directionIsDiagonal (dir))
		return !implies (ENTERABLE (-2), ENTERABLE (-3)) ||
		       !implies (ENTERABLE (2), ENTERABLE (3));
	else 
		return !implies (ENTERABLE (-1), ENTERABLE (-2)) ||
		       !implies (ENTERABLE (1), ENTERABLE (2));
#undef ENTERABLE
}
*/
DirectionSet AStar_jsp::forcedNeighbours (const Coordinate& coord, 
                                          Direction dir)
{
  DirectionSet dirs;
	if (dir == Direction::NO_DIRECTION)
		return dirs;

#define ENTERABLE(n) isEnterable (adjustInDirection (coord, RotateDirection(dir, n)))

	if (directionIsDiagonal (dir)) {

		if (!implies (ENTERABLE (6), ENTERABLE (5)))
			addDirectionToSet (dirs, RotateDirection(dir,6));
		if (!implies (ENTERABLE (2), ENTERABLE (3)))
			addDirectionToSet (dirs, RotateDirection(dir,2));
	
  }
	else {
	
    if (!implies (ENTERABLE (7), ENTERABLE (6)))
			addDirectionToSet (dirs, RotateDirection(dir, 7));
		if (!implies (ENTERABLE (1), ENTERABLE (2)))
			addDirectionToSet (dirs, RotateDirection(dir,1));
	
  }	
#undef ENTERABLE	
	return dirs;
}

static DirectionSet MakeFullDirectionSet()
{
  DirectionSet dirs;
  dirs.insert(Direction::N);
  dirs.insert(Direction::NE);
  dirs.insert(Direction::E);
  dirs.insert(Direction::SE);
  dirs.insert(Direction::S);
  dirs.insert(Direction::SW);
  dirs.insert(Direction::W);
  dirs.insert(Direction::NW);
  return dirs;
}

static DirectionSet FindStrightDirections(Direction dir)
{
  DirectionSet dirs;
  if (!directionIsDiagonal(dir))
    return dirs;

  dirs.insert(RotateDirection(dir, 1));
  dirs.insert(RotateDirection(dir, -1));
  return dirs;
}

static DirectionSet naturalNeighbours (Direction dir)
{
	if (dir == Direction::NO_DIRECTION)
		return MakeFullDirectionSet();

	DirectionSet dirs;
	addDirectionToSet (dirs, dir);
  addDirectionsToSet(dirs, FindStrightDirections(dir));
	return dirs;
}

void AStar_jsp::addToOpenSet (int node, 
                              int nodeFrom)
{
	Coordinate nodeCoord = getCoord (this->m_bounds, node);
	Coordinate nodeFromCoord = getCoord (this->m_bounds, nodeFrom);

	if (!exists (this->m_open, node)) {
		this->m_cameFrom[node] = nodeFrom;
		this->m_gScores[node] = this->m_gScores[nodeFrom] + preciseDistance(nodeFromCoord, nodeCoord);
		insert (this->m_open, node, this->m_gScores[node] + estimateDistance(nodeCoord, getCoord (this->m_bounds, this->m_goal)));
	}
	else if (this->m_gScores[node] > this->m_gScores[nodeFrom] + preciseDistance(nodeFromCoord, nodeCoord)) 
  {
		this->m_cameFrom[node] = nodeFrom;
		double oldGScore = this->m_gScores[node];
		this->m_gScores[node] = this->m_gScores[nodeFrom] + preciseDistance (nodeFromCoord, nodeCoord);
		double newPri = priorityOf (this->m_open, node) - oldGScore + this->m_gScores[node];
		changePriority (this->m_open, node, newPri);
	}	
}

bool AStar_jsp::HasForcedNeighbours(Coordinate coord, Direction dir)
{
  return !forcedNeighbours(coord, dir).empty();
}

// directly translated from "algorithm 2" in the paper
int AStar_jsp::jump(Direction dir, int start)
{
	Coordinate coord = adjustInDirection (getCoord (this->m_bounds, start), dir);
	int node = getIndex(this->m_bounds, coord);
	if (!isEnterable (coord))
		return -1;

	if (node == this->m_goal || HasForcedNeighbours(coord, dir)) 
  {
		return node;
	}

  auto strights = FindStrightDirections(dir);
  if (!strights.empty())
  {
    for (auto it = begin(strights); it != end(strights); ++it)
    {
      if (jump(*it, node) >= 0)
        return node;
    }
  }

	return jump (dir, node);
}

// path interpolation between jump points in here
int AStar_jsp::nextNodeInSolution (int *target,
                               int node)
{
	Coordinate c = getCoord (this->m_bounds, node);
	Coordinate cTarget = getCoord (this->m_bounds, *target);

	if (c.x < cTarget.x) 
		c.x++;
	else if (c.x > cTarget.x)
		c.x--;

	if (c.y < cTarget.y) 
		c.y++;
	else if (c.y > cTarget.y)
		c.y--;

	node = getIndex (this->m_bounds, c);

	if (node == *target)
		*target = this->m_cameFrom[*target];

	return node;
}

// a bit more complex than the usual A* solution-recording method,
// due to the need to interpolate path chunks
int* AStar_jsp::recordSolution ()
{
	int rvLen = 1;
	*this->m_solutionLength = 0;
	int target = this->m_goal;
	int *rv = static_cast<int*>(malloc (rvLen * sizeof (int)));
	int i = this->m_goal;

	for (;;) {
		i = nextNodeInSolution (&target, i);
		rv[*this->m_solutionLength] = i;
		(*this->m_solutionLength)++;
		if (*this->m_solutionLength >= rvLen) {
			rvLen *= 2;
			rv = static_cast<int*>(realloc (rv, rvLen * sizeof (int)));
			if (!rv)
				return NULL;
		}
		if (i == this->m_start)
			break;
	}

	(*this->m_solutionLength)--; // don't include the starting tile
	return rv;
}

static Direction MergeDirections(Direction d1, Direction d2)
{

  //  7 0 1
  //  6 x 2
  //  5 4 3
  //
  // We can merge only directions 0,2,4 and 6
  //

  if (d1 == Direction::NO_DIRECTION)
    return d2;
  if (d2 == Direction::NO_DIRECTION)
    return d1;

  if (directionIsDiagonal(d1) || directionIsDiagonal(d2))
    return Direction::NO_DIRECTION;

  int tmp1 = static_cast<int>(d1);
  int tmp2 = static_cast<int>(d2);

  if (std::abs(tmp1-tmp2) > 2)
    return Direction::NO_DIRECTION;

  return static_cast<Direction>((tmp1+tmp2)/2);
}

static Direction directionOfMove (Coordinate to, Coordinate from)
{

  Direction dir_x = Direction::NO_DIRECTION;
  Direction dir_y = Direction::NO_DIRECTION;

  if (from.x < to.x)
    dir_x = Direction::E;
  else if (from.x > to.x)
    dir_x = Direction::W;

  if (from.y < to.y)
    dir_y = Direction::S;
  else if (from.y > to.y)
    dir_y = Direction::N;

  return MergeDirections(dir_x, dir_y);
}

Direction AStar_jsp::directionWeCameFrom (int node, int nodeFrom)
{
	if (nodeFrom == -1)
		return Direction::NO_DIRECTION;

	return directionOfMove (getCoord (this->m_bounds, node), 
				getCoord (this->m_bounds, nodeFrom));
}

int *astar_compute (const char *grid, 
		    int *solLength, 
		    int boundX, 
		    int boundY, 
		    int start, 
		    int end)
{
	*solLength = -1;
	AStar_jsp astar;
	Coordinate bounds(boundX, boundY);

	int size = bounds.x * bounds.y;


	if (start >= size || start < 0 || end >= size || end < 0)
		return NULL;

	Coordinate startCoord = getCoord (bounds, start);
	Coordinate endCoord = getCoord (bounds, end);

	if (!contained (bounds, startCoord) || !contained (bounds, endCoord))
		return NULL;

	queue *open = createQueue();
  std::vector<char> closed(size);
  std::vector<double> gScores(size);
  std::vector<int> cameFrom(size);

	astar.m_solutionLength = solLength;
	astar.m_bounds = bounds;
	astar.m_start = start;
	astar.m_goal = end;
	astar.m_grid = grid;
	astar.m_open = open;
	astar.m_closed = closed.data();
	astar.m_gScores = gScores.data();
	astar.m_cameFrom = cameFrom.data();

	memset (closed.data(), 0, sizeof(closed));

	gScores[start] = 0;
	cameFrom[start] = -1;

	insert (open, start, estimateDistance (startCoord, endCoord));
	while (open->size) {
		int node = findMin (open)->value; 
		Coordinate nodeCoord = getCoord (bounds, node);
		if (nodeCoord.x == endCoord.x && nodeCoord.y == endCoord.y) {
			freeQueue (open);
			return astar.recordSolution ();
		}

		deleteMin (open);
		closed[node] = 1;

    Direction from = astar.directionWeCameFrom (node, cameFrom[node]);

    DirectionSet fneighbours = astar.forcedNeighbours(nodeCoord, from);
    DirectionSet nneighbours = naturalNeighbours(from);

		DirectionSet dirs(fneighbours);
    dirs.insert(nneighbours.begin(), nneighbours.end());

		for (Direction dir = nextDirectionInSet(dirs); dir != Direction::NO_DIRECTION; dir = nextDirectionInSet(dirs))
		{
			int newNode = astar.jump (dir, node);
			Coordinate newCoord = getCoord (bounds, newNode);

			// this'll also bail out if jump() returned -1
			if (!contained (bounds, newCoord))
				continue;

			if (closed[newNode])
				continue;
			
			astar.addToOpenSet (newNode, node);

		}
	}
	freeQueue (open);
	return NULL;
}



int *astar_unopt_compute (const char *grid, 
		    int *solLength, 
		    int boundX, 
		    int boundY, 
		    int start, 
		    int end)
{
	AStar_jsp astar;
	Coordinate bounds(boundX, boundY);

	int size = bounds.x * bounds.y;


	if (start >= size || start < 0 || end >= size || end < 0)
		return NULL;

	Coordinate startCoord = getCoord (bounds, start);
	Coordinate endCoord = getCoord (bounds, end);

	if (!contained (bounds, startCoord) || !contained (bounds, endCoord))
		return NULL;

	queue *open = createQueue();
  std::vector<char> closed(size);
  std::vector<double> gScores(size);
  std::vector<int> cameFrom(size);

	astar.m_solutionLength = solLength;
	*astar.m_solutionLength = -1;
	astar.m_bounds = bounds;
	astar.m_start = start;
	astar.m_goal = end;
	astar.m_grid = grid;
	astar.m_open = open;
	astar.m_closed = closed.data();
	astar.m_gScores = gScores.data();
	astar.m_cameFrom = cameFrom.data();

	memset (closed.data(), 0, sizeof(closed));

	gScores[start] = 0;
	cameFrom[start] = -1;

	insert (open, start, estimateDistance (startCoord, endCoord));

	while (open->size) {
		int node = findMin (open)->value; 
		Coordinate nodeCoord = getCoord (bounds, node);
		if (nodeCoord.x == endCoord.x && nodeCoord.y == endCoord.y) {
			freeQueue (open);
			return astar.recordSolution ();
		}

		deleteMin (open);
		closed[node] = 1;

    DirectionSet dirs = MakeFullDirectionSet();

		for (Direction dir = nextDirectionInSet(dirs); dir != Direction::NO_DIRECTION; dir = nextDirectionInSet(dirs))
		{
			Coordinate newCoord = adjustInDirection (nodeCoord, dir);
			int newNode = getIndex (bounds, newCoord);

			if (!contained (bounds, newCoord) || !grid[newNode])
				continue;

			if (closed[newNode])
				continue;
			
			astar.addToOpenSet (newNode, node);

		}
	}
	freeQueue (open);
	return NULL;
}

