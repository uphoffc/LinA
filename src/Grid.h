#ifndef GRID_H_
#define GRID_H_

#include <cstring>

template<typename T>
class Grid {
public:
  Grid(int X, int Y);
  ~Grid();
  
  /// Implements periodic boundary conditions
  inline int periodicIndex(int i, int N) {
    int pi;
    if (i >= 0) {
      pi = i % N;
    } else {
      pi = N-1 - (-i-1)%N;
    }
    return pi;
  }
  
  inline T& get(int x, int y) {
    return m_data[periodicIndex(y, m_Y) * m_X + periodicIndex(x, m_X)];
  }
  
  inline int X() const {
    return m_X;
  }
  
  inline int Y() const {
    return m_Y;
  }

private:
  int m_X;
  int m_Y;
  T* m_data;
};

template<typename T>
Grid<T>::Grid(int X, int Y)
  : m_X(X), m_Y(Y)
{
  m_data = new T[X*Y];
  memset(m_data, 0, X*Y*sizeof(T));
}

template<typename T>
Grid<T>::~Grid() {
  delete[] m_data;
}

#endif // GRID_H_
