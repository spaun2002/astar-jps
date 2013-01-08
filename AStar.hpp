#ifndef ASTAR_H_
#define ASTAR_H_

#include <vector>
#include <type_traits>
#include <queue>


namespace AStarJPS
{


struct Coordinate {
  union
  {
    struct {
      int x;
      int y;
    };
    struct {
      size_t width;
      size_t height;
    };
  };
  Coordinate(int _x, int _y) : x(_x), y(_y){}
  Coordinate(){};
};

class IMap
{
public:
  virtual ~IMap() = 0 {};
  virtual size_t GetIndexByCoord(const Coordinate& coord) const = 0;
  virtual Coordinate GetCoordByIndex(size_t index) const = 0;
  virtual bool IsEnterable(size_t index) const = 0;
  virtual Coordinate GetBounds() const = 0;
  virtual size_t GetMaxIndex() const = 0;
  virtual bool Contains(const Coordinate& coord) const = 0;
};

class MapSimpleCoordConvImpl : virtual public IMap
{
  size_t GetIndexByCoord(const Coordinate& coord) const
  {
    return coord.x + this->GetBounds().width * coord.y;
  }
  Coordinate GetCoordByIndex(size_t index) const
  {
    if (index >= 0 && index <= this->GetMaxIndex())
    {
      auto width = this->GetBounds().width;
      return Coordinate(index % width, index / width);
    }
    return Coordinate(-1,-1);
  }
};

template <class _Tx>
struct simple_isenterable
{
  bool operator()(const typename std::enable_if<std::is_convertible<_Tx, bool>::value, _Tx>::type & val) const
  {
    return static_cast<bool>(val);
  }
};

template <class _UnderlayingType, 
          class _IsEnterable = simple_isenterable<_UnderlayingType> >
class SimpleMap : public MapSimpleCoordConvImpl
{
public:
  typedef _UnderlayingType data_t;
  typedef SimpleMap<_UnderlayingType, _IsEnterable> type_t;

  SimpleMap(const data_t* data, size_t width, size_t height)
    : m_data(data)
    , m_width(width)
    , m_height(height)
    , m_maxIndex(width*height-1)
    , m_pr()
  {

  }

  ~SimpleMap()
  {

  }

  virtual bool IsEnterable( size_t index ) const
  {
    return m_pr(GetCellDataByIndex(index));
  }

  virtual Coordinate GetBounds() const
  {
    return Coordinate(m_width, m_height);
  }

  virtual size_t GetMaxIndex() const
  {
    return m_maxIndex;
  }

  data_t GetCellDataByIndex(size_t index) const
  {
    if (index > GetMaxIndex())
      throw std::out_of_range("Index is out of range");

    return m_data[index];
  }

  virtual bool Contains( const Coordinate& coord ) const
  {
    return coord.x >= 0 
      && coord.x < m_width 
      && coord.y >= 0 
      && coord.y < m_height;
  }

private:
  const data_t* m_data;
  size_t m_width;
  size_t m_height;
  size_t m_maxIndex;
  _IsEnterable m_pr;
};

/* Run A* pathfinding over uniform-cost 2d grids using jump point search.

   grid: 0 if obstructed, non-0 if non-obstructed (the value is ignored beyond that).
   solLength: out-parameter, will contain the solution length
   boundX: width of the grid
   boundY: height of the grid
   start: index of the starting node in the grid
   end: index of the ending node in the grid

   return value: Array of node indexes making up the solution, of length solLength, in reverse order
 */

std::vector<int> AStarCompute(const IMap& map, 
                              int& solutionLength, 
                              const Coordinate& start, 
                              const Coordinate& end);

std::vector<int> AStarUnoptCompute(const IMap& map, 
                                   int& solutionLength, 
                                   const Coordinate& start, 
                                   const Coordinate& end);

//std::vector<int> astar_unopt_compute (const char *grid, 
//		    int *solLength, 
//		    int boundX, 
//		    int boundY, 
//		    int start, 
//		    int end);


///* Compute cell indexes from cell coordinates and the grid width */
//int astar_getIndexByWidth (int width, int x, int y);
//
///* Compute coordinates from a cell index and the grid width */
//void astar_getCoordByWidth (int width, int node, int *x, int *y);

}//namespace AStarJPS

#endif
