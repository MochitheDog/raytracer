#include "sphere.h"
class Chessboard : public Model3D {
public:
  Chessboard() {
    isChessboard = true;
  }
  ~Chessboard();
  
};

Model3D *add_chessboard(Model3D *, Point, float, float [], float [], float [], float, float, int);
Vector chessboard_normal()
{
  Vector n = {0,1,0}; // up
  return n;
}