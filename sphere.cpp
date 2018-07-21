#include "sphere.h"
#include <stdlib.h>
#include <math.h>
#include <cfloat> // FLT_MAX
#include <iostream>
#include <stdio.h>
extern Spheres *scene;
extern Spheres* chessboard;
/**********************************************************************
 * This function intersects a ray with a given sphere 'sph'. You should
 * use the parametric representation of a line and do the intersection.
 * The function should return the parameter value for the intersection, 
 * which will be compared with others to determine which intersection
 * is closest. The value -1.0 is returned if there is no intersection
 *
 * If there is an intersection, the point of intersection should be
 * stored in the "hit" variable
 **********************************************************************/
float intersect_sphere(Point o, Vector v, Spheres *sph, Point *hit) {
  // how to calculate intersection?
  // Parametric representation of a line: origin + t*Direction vector
  //https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
	//https://en.wikipedia.org/wiki/Line-sphere_intersection
  
  Vector L = get_vec(sph->center,o);
  float a = vec_dot(v,v);
  float b = 2*(vec_dot(v, L));
  float c = vec_dot(L, L) - (sph->radius*sph->radius);
  float discrim = (b*b) - (4*a*c); // b^2 - 4ac
  float x1, x2, t;
 
  if (discrim < 0)
  {
    return -1.0; // no intersection
  }
  else // one or two intersections
  {
    // Quadratic formula
    x1 = (-b + sqrt(discrim)) / (2*a);
    x2 = (-b - sqrt(discrim)) / (2*a);
    // If x is negative, intersection is behind eye
    // if both x negative, don't draw
    // else draw the closest point in front of the camera
    if (x1 < 0 && x2 < 0)
    {
      return -1.0;
    }
    else
    {
      if (x1 < x2)
      {
        if (x1 > 0)
        {
          t = x1;
        }
        else
        {
          t = x2;
        }
      }
      else // x1 > x2
      {
        if (x2 > 0)
        {
          t = x2;
        }
        else
        {
          t = x1;
        }
      }
      
      if (hit != nullptr)
      {
        Point intersectpoint = get_point(o, vec_scale(v, t));
        hit->x = intersectpoint.x;
        hit->y = intersectpoint.y;
        hit->z = intersectpoint.z;
      }
      
      return t;
    }
  }
}

/*********************************************************************
 * This function returns a pointer to the sphere object that the
 * ray intersects first; NULL if no intersection. You should decide
 * which arguments to use for the function. For exmaple, note that you
 * should return the point of intersection to the calling function.
 **********************************************************************/
Spheres* intersect_scene(Point eye, Vector ray, Point* hit, Spheres* ignoreMe ) {

  // Iterate over objects
  Spheres* curr_sphere = scene;
  Spheres* closest_sphere = nullptr;
  float closest_dist = FLT_MAX; // from cfloat
  
  while (curr_sphere != nullptr)
  {
    Point hitpoint;
    // For reflection etc- can't reflect off itself
    if (curr_sphere != nullptr && curr_sphere == ignoreMe) 
    {
      curr_sphere = curr_sphere->next;
      continue;
    }
    float dist = intersect_sphere(eye, ray, curr_sphere, hit);
    
    if (dist != -1.0 && dist < closest_dist)
    {
      closest_sphere = curr_sphere;
      closest_dist = dist;
    }
    curr_sphere = curr_sphere->next;
  }
	return closest_sphere;
}

/*****************************************************
 * This function adds a sphere into the sphere list
 *
 * You need not change this.
 *****************************************************/
Spheres *add_sphere(Spheres *slist, Point ctr, float rad, float amb[],
		    float dif[], float spe[], float shine, 
		    float refl, float refract, int sindex ) {
  Spheres *new_sphere;

  new_sphere = (Spheres *)malloc(sizeof(Spheres));
  new_sphere->index = sindex;
  new_sphere->center = ctr;
  new_sphere->radius = rad;
  (new_sphere->mat_ambient)[0] = amb[0];
  (new_sphere->mat_ambient)[1] = amb[1];
  (new_sphere->mat_ambient)[2] = amb[2];
  (new_sphere->mat_diffuse)[0] = dif[0];
  (new_sphere->mat_diffuse)[1] = dif[1];
  (new_sphere->mat_diffuse)[2] = dif[2];
  (new_sphere->mat_specular)[0] = spe[0];
  (new_sphere->mat_specular)[1] = spe[1];
  (new_sphere->mat_specular)[2] = spe[2];
  new_sphere->mat_shineness = shine;
  new_sphere->reflectance = refl;
  new_sphere->next = nullptr;
  Vector n = {0,0,0};
  new_sphere->normal = n; // not used
  new_sphere->refraction = refract;
  if (slist == nullptr) { // first object
    slist = new_sphere;
  } else { // insert at the beginning
    new_sphere->next = slist;
    slist = new_sphere;
  }

  return slist;
}

Spheres *add_chessboard(Spheres *slist, Point ctr, float rad, float amb[],
		    float dif[], float spe[], float shine, 
		    float refl, int sindex, Vector normal) {
  Spheres *new_sphere;

  new_sphere = (Spheres *)malloc(sizeof(Spheres));
  new_sphere->index = sindex;
  new_sphere->center = ctr;
  new_sphere->radius = rad;
  (new_sphere->mat_ambient)[0] = amb[0];
  (new_sphere->mat_ambient)[1] = amb[1];
  (new_sphere->mat_ambient)[2] = amb[2];
  (new_sphere->mat_diffuse)[0] = dif[0];
  (new_sphere->mat_diffuse)[1] = dif[1];
  (new_sphere->mat_diffuse)[2] = dif[2];
  (new_sphere->mat_specular)[0] = spe[0];
  (new_sphere->mat_specular)[1] = spe[1];
  (new_sphere->mat_specular)[2] = spe[2];
  new_sphere->mat_shineness = shine;
  new_sphere->reflectance = refl;
  new_sphere->next = nullptr;
  new_sphere->normal = normal;
  new_sphere->isChessboard = true;
  if (slist == nullptr) { // first object
    slist = new_sphere;
  } else { // insert at the beginning
    new_sphere->next = slist;
    slist = new_sphere;
  }
  
  return slist;
}

/******************************************
 * computes a sphere normal - done for you
 ******************************************/
Vector sphere_normal(Point q, Spheres *sph) {
  Vector rc;
  rc = get_vec(sph->center, q);
  normalize(&rc);
  return rc;
}

// Return point of intersection on chessboard in hit variable 
// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
bool intersect_board(Point origin, Vector vec, Point* hit)
{
  Vector norm = chessboard->normal;
  normalize(&norm);
  Vector v = vec;
  normalize(&v);
  float dotnl = vec_dot(norm, v);
  if (dotnl > 0)
  {
    Point planeCenter = chessboard->center;
    Vector originToPlane = get_vec(planeCenter, origin);
    float dist = vec_dot(originToPlane, norm) / dotnl;
    Point intersectpoint = get_point(origin, vec_scale(vec, dist));
    if (dist < 0) return false;
    hit->x = intersectpoint.x;
    hit->y = intersectpoint.y;
    hit->z = intersectpoint.z;
    // 8 by 8 board
    if(intersectpoint.x >= 5.0 || intersectpoint.x < -5.0)
      return false;
    if(intersectpoint.z >= 2 || intersectpoint.z < -6)
      return false;
    
    return true;
  }
  return false;
}