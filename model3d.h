#include <vector>
class Spheres;
class Chessboard;

class Model3D {
public:

  virtual float intersect_model(Point, Vector, Model3D*, Point*);
  virtual float intersect_model(Point, Vector, Spheres*, Point*);
  int index;               // identifies a sphere; must be greater than 0
  float mat_ambient[3];    // material property used in Phong model
  float mat_diffuse[3];
  float mat_specular[3];
  float mat_shineness;

  float reflectance;       // this number [0,1] determines how much 
                           // reflected light contributes to the color
                           // of a pixel
  Model3D* next;
  bool isChessboard = false; // TODO remove

};