#include <stdio.h>
#include <GL/glut.h>
#include <math.h>
#include "global.h"
#include "sphere.h"
#include <iostream>
#include <algorithm>
#include <utility>
//
// Global variables
//
extern int win_width;
extern int win_height;

extern GLfloat frame[WIN_HEIGHT][WIN_WIDTH][3];  

extern float image_width;
extern float image_height;

extern Point eye_pos;
extern float image_plane;
extern RGB_float background_clr;
extern RGB_float null_clr;

extern Spheres *scene;
extern Spheres* chessboard;
// light 1 position and color
extern Point light1;
extern float light1_ambient[3];
extern float light1_diffuse[3];
extern float light1_specular[3];

// global ambient term
extern float global_ambient[3];

// light decay parameters
extern float decay_a;
extern float decay_b;
extern float decay_c;

extern bool shadow_on;
extern bool reflection_on;
extern bool chessboard_on;
extern bool supersample_on;
extern bool stochastic_on;
extern bool refraction_on;
extern int step_max;

/////////////////////////////////////////////////////////////////////

// Return true if blocked by another object, false if not
bool shadowCheck(Point p, Vector shadowRay, Spheres *spheres, Spheres* ignoreMe)
{
  // check if shadow ray intersects object on the way to light source
  Spheres* curr_sphere = spheres;
  float lightDist = vec_len(shadowRay);
  while (curr_sphere != nullptr)
  {
    if (curr_sphere == ignoreMe)
    {
      curr_sphere = curr_sphere->next;
      continue;
    }
    float dist = intersect_sphere(p, shadowRay, curr_sphere, nullptr);
    if (dist != -1.0)
    {
      return true;
    }
    curr_sphere = curr_sphere->next;
  }
  return false;
}

/*********************************************************************
 * Phong illumination 
 * p - pixel point
 * sph - this object being illuminated
 * l - light vector (normalized)
 * n - normal vector (normalized)
 * vv - eye vector (normalized)
 * r - reflection vector (normalized)
 *********************************************************************/
RGB_float phong(Point p, Vector vRay, Spheres *sph ) {
	RGB_float color;
  Vector l = get_vec(p,light1); // light vector
  Vector n;
  n = sphere_normal(p,sph); // normal vector
  
  Vector toEye = get_vec(p, eye_pos); // eye vector
  normalize(&l);
  normalize(&n);
  Vector r = vec_minus(vec_scale(n,2*vec_dot(n,l)),l); // reflection vector
  normalize(&r);
  normalize(&toEye);
  
  RGB_float globalAmbient;
  globalAmbient.r = global_ambient[0] *sph->mat_ambient[0]; // i hope this is k_ga
  globalAmbient.g = global_ambient[1] *sph->mat_ambient[1];
  globalAmbient.b = global_ambient[2] *sph->mat_ambient[2];

  // Shadow ray: from point to light source
  // if shadow ray hits light source, compute light contribution
  // if shadow ray meets object, don't.
  Vector shadowRay = l;
  if (shadow_on) // only if shadows option is on
  {
    if (shadowCheck(p, shadowRay, scene, sph)) // blocked
    {
      // if not lit by light, it's lit by... global ambient?
      return globalAmbient;
    }
  }
  RGB_float ambient, diffuse, specular;
  ambient.r = light1_ambient[0] * sph->mat_ambient[0] + globalAmbient.r;
  ambient.g = light1_ambient[1] * sph->mat_ambient[1] + globalAmbient.g;
  ambient.b = light1_ambient[2] * sph->mat_ambient[2] + globalAmbient.b;
  
  float n_dot_l = vec_dot(n,l);
  diffuse.r = light1_diffuse[0] * sph->mat_diffuse[0] * n_dot_l;
  diffuse.g = light1_diffuse[1] * sph->mat_diffuse[1] * n_dot_l;
  diffuse.b = light1_diffuse[2] * sph->mat_diffuse[2] * n_dot_l;
  // ray of reflection
  float r_dot_v = vec_dot(r, toEye);
  float r_dot_v_N = pow(r_dot_v, sph->mat_shineness);
  specular.r = light1_specular[0] * sph->mat_specular[0] * r_dot_v_N;
  specular.g = light1_specular[1] * sph->mat_specular[1] * r_dot_v_N;
  specular.b = light1_specular[2] * sph->mat_specular[2] * r_dot_v_N;

  float lightDist = vec_len(get_vec(light1, p)); 
  float abcd = 1.0 / // decay factor
          (decay_a + (decay_b*lightDist) + (decay_c*(lightDist*lightDist)));
  diffuse = clr_scale(diffuse,abcd);
  specular = clr_scale(specular,abcd);
  // check normal-light direction-
  // don't add specular if facing away from light source
  if (n_dot_l >= 0.1) 
  {
    color = clr_add(diffuse, specular);
  }
  else
  {
    color = diffuse;
  }
  color = clr_add(color, ambient);
	return color;
}

RGB_float chessboard_phong(Point p, Vector vRay, bool whiteSquare)
{
  // Returns black or white depending on p's x and z
  RGB_float color;
  Vector l = get_vec(p,light1); // light vector
  Vector n = chessboard->normal; // normal vector
  Vector toEye = get_vec(p, eye_pos); // eye vector
  normalize(&l);
  normalize(&n);
  Vector r = vec_minus(vec_scale(n,2*vec_dot(n,l)),l); // reflection vector
  normalize(&r);
  normalize(&toEye);
  RGB_float globalAmbient;
  globalAmbient.r = global_ambient[0] *chessboard->mat_ambient[0]; // i hope this is k_ga
  globalAmbient.g = global_ambient[1] *chessboard->mat_ambient[1];
  globalAmbient.b = global_ambient[2] *chessboard->mat_ambient[2];

  // Shadow ray: from point to light source
  // if shadow ray hits light source, compute light contribution
  // if shadow ray meets object, don't.
  Vector shadowRay = get_vec(p,light1);
  //std::cout << "phong check p " << p.x << "," << p.y << "," << p.z << std::endl;
  if (shadow_on) // only if shadows option is on
  {
    if (shadowCheck(p, shadowRay, scene, nullptr)) // blocked
    {
      // if not lit by light, it's lit by... global ambient?
      //std::cout << "blocked" << std::endl;
      return globalAmbient;
    }
  }
  float fcolor;
  if (whiteSquare)
    fcolor = -1.0; // yeah idk what's up
  else
    fcolor = 0.0;

  RGB_float ambient, diffuse, specular;
  ambient.r = light1_ambient[0] * fcolor + globalAmbient.r;
  ambient.g = light1_ambient[1] * fcolor + globalAmbient.g;
  ambient.b = light1_ambient[2] * fcolor + globalAmbient.b;
  //std::cout << "amb " << ambient.r << std::endl;
  
  float n_dot_l = vec_dot(n,l);
  diffuse.r = light1_diffuse[0] * fcolor * n_dot_l;
  diffuse.g = light1_diffuse[1] * fcolor * n_dot_l;
  diffuse.b = light1_diffuse[2] * fcolor * n_dot_l;
  //std::cout << "df " << diffuse.r << " " << diffuse.g << " " << diffuse.b << std::endl;
  //std::cout << "df ndl " << n_dot_l << std::endl;
  // ray of reflection
  float r_dot_v = vec_dot(r, toEye);
  float r_dot_v_N = pow(r_dot_v, chessboard->mat_shineness);
  specular.r = light1_specular[0] * fcolor * r_dot_v_N;
  specular.g = light1_specular[1] * fcolor * r_dot_v_N;
  specular.b = light1_specular[2] * fcolor * r_dot_v_N;

  float lightDist = vec_len(get_vec(light1, p)); 
  float abcd = 1.0 / // decay factor
          (decay_a + (decay_b*lightDist) + (decay_c*(lightDist*lightDist)));
  diffuse = clr_scale(diffuse,abcd);
  specular = clr_scale(specular,abcd);
  if (n_dot_l >= 0.1)
  {
    color = clr_add(diffuse, specular);
  }
  else
  {
    color = diffuse;
  }
  color = clr_add(color, ambient);
  //std::cout << "phong board: " << color.r << ", "<<color.g << ", " << color.b << std::endl;
  return color;
}

int random1()
{
  if (rand()%2 == 1)
    return 1;
  else
    return -1;
}

Vector randomVector(Vector normal) {
  // https://gamedev.stackexchange.com/questions/26789/random-vector-within-a-cone
  Vector v;
  float x, y, z;
  do {
    x = 2.0*((float)rand() / RAND_MAX)-1; //2*(number between 0 and 1)
    y = 2.0*((float)rand() / RAND_MAX)-1;
    z = 2.0*((float)rand() / RAND_MAX)-1;
    v = {x, y, z};
  } while ( /*x*x+y*y > 1.0 && */vec_dot(v,normal) < 0);
  normalize(&v);
  //std::cout << v.x << ", " << v.y << ", " << v.z << std::endl;
  return v;
}

// Vector getRefractionRay(float refractionIndex, Vector incident, Point hitpoint, Spheres* obj)
// {
//   // https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
//   float i = vec_dot(sphere_normal(hitpoint, obj),incident);
//   i = clamp(i, -1, 1);
//   float e = 1;
//   float ior = refractionIndex;
//   float n = sphere_normal(hitpoint, obj);
//   if (i < 0)
//   {
//     i = -i;
//   }
//   else
//   {
//     swap(i, e);
//     n = -n;
//   }
//   float div = e/ior;
//   float k = 1 - div*div*(1-i*i);
//   Vector zero = {0,0,0};
//   if (k < 0)
//     return zero;
//   else
//   {
//     return vec_add(vec_scale(incident,div), vec_dot(div*i-sqrtf(k), n)); 
//   }
// }
/************************************************************************
 * This is the recursive ray tracer - you need to implement this!
 * You should decide what arguments to use.
 ************************************************************************/
RGB_float recursive_ray_trace(Vector ray, Point e, int iteration, Spheres* ignoreMe) {

  Point* hit = new Point();
  // See if ray intersects any objects in scene
  Spheres* obj = intersect_scene(e, ray, hit, ignoreMe);
  RGB_float color;
  Vector n; // normal vector
  if (obj != nullptr)
  {
    color = phong(*hit, ray, obj);
    // Reflected light: color = phongcolor + reflectance * reflected light accumulated thru recursive rt
    if (reflection_on && iteration < step_max)
    {
      // get reflection ray
      Vector l = get_vec(*hit,e); // (light but not light) vec from previous point
      
      n = sphere_normal(*hit,(Spheres*)obj);
      normalize(&l);
      normalize(&n);
      Vector r = vec_minus(vec_scale(n,2*vec_dot(n,l)),l); // reflection vector
      normalize(&r);
      RGB_float ref_color = recursive_ray_trace(r, *hit, iteration+1, obj);
      if (ref_color.r != background_clr.r ||
          ref_color.g != background_clr.g ||
          ref_color.b != background_clr.b)
      {
        ref_color = clr_scale(ref_color, obj->reflectance);
        if (ref_color.r < 0)
          ref_color.r = 0;
        if (ref_color.g < 0)
          ref_color.g = 0;
        if (ref_color.b < 0)
          ref_color.b = 0;
        color = clr_add(color, ref_color);
      } 
    }
    // if (refraction_on && iteration < step_max)
    // {
    //   // https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
    //   Vector hitNormal;
    //   Point refractOut;
    //   Vector refraction_ray;

      
    //   //float cos = clamp(vec_dot(I,N), -1, 1);
    //   //Vector norm = sph->normal;
    //   //Vector refractionVector = {0,0,0};
    //   //RGB_float refraction_color = {0,0,0};
    //   //float ndotincident = vec_dot(norm,)
    //   //color = clr_add(color, refraction_color);
    // }
    if (stochastic_on && iteration <= step_max)
    {
      RGB_float sto_color = {0,0,0};
      for (int i = 0; i < STOCHASTIC_RAYS; i++)
      {
        RGB_float temp_color;
        // diffuse reflections... shoot rays in random directions from hit
        n = sphere_normal(*hit,(Spheres*)obj);
        Vector randomRay = randomVector(n);
        normalize(&randomRay);
        //std::cout << "baah" << std::endl;
        temp_color = recursive_ray_trace(randomRay, *hit, iteration+1, obj);
        if (temp_color.r == background_clr.r &&
          temp_color.g == background_clr.g &&
          temp_color.b == background_clr.b) continue;
        //std::cout << sto_color.r << ", " << sto_color.g << ", " << sto_color.b << std::endl;
        temp_color = clr_scale(temp_color, obj->reflectance);
        if (temp_color.r < 0)
          temp_color.r = 0;
        if (temp_color.g < 0)
          temp_color.g = 0;
        if (temp_color.b < 0)
          temp_color.b = 0;
        sto_color = clr_add(sto_color, temp_color);
      }
      sto_color.r = sto_color.r / STOCHASTIC_RAYS;
      sto_color.g = sto_color.g / STOCHASTIC_RAYS;
      sto_color.b = sto_color.b / STOCHASTIC_RAYS;
      color = clr_add(color, sto_color);
    }

    // no weird shadows when saved
    if (color.r < 0)
      color.r = 0;
    if (color.g < 0)
      color.g = 0;
    if (color.b < 0)
      color.b = 0;

    return color;
  }
  // didn't hit a sphere- check if hit the chessboard
  
  if (chessboard_on)
  {
    Point* boardhit = new Point();
    if (intersect_board(e, ray, boardhit))
    {
      bool square_white = false;
      int x = ((int)boardhit->x);
      int z = ((int)boardhit->z);
      if (x % 2 == 0 && z%2 == 0)
        square_white = false;
      else if (x % 2 == 0 && z%2 != 0)
        square_white = true;
      else if (x % 2 != 0 && z%2 == 0)
        square_white = true;
      else
        square_white = false;

      color = chessboard_phong(*boardhit, ray, square_white);
      if (reflection_on && iteration < step_max)
      {
        // get reflection ray
        Vector l = get_vec(*hit,e); // (light but not light) vec from previous point
        Vector n = chessboard->normal; 
        normalize(&l);
        normalize(&n);
        Vector r = vec_minus(vec_scale(n,2*vec_dot(n,l)),l); // reflection vector
        normalize(&r);
        RGB_float ref_color = recursive_ray_trace(r, *boardhit, iteration+1, nullptr);
        //std::cout << ref_color.r << ", " << ref_color.g << ", " << ref_color.b << std::endl;
        
        ref_color = clr_scale(ref_color, chessboard->reflectance);
        color = clr_add(color, ref_color);
      }
      return color;
    }
  }
  //std::cout << "board" << std::endl;
  return background_clr;

}

/*********************************************************************
 * This function traverses all the pixels and cast rays. It calls the
 * recursive ray tracer and assign return color to frame
 *
 * You should not need to change it except for the call to the recursive
 * ray tracer. Feel free to change other parts of the function however,
 * if you must.
 *********************************************************************/
void ray_trace() {
  //int i, j;
  float x_grid_size = image_width / float(win_width);
  float y_grid_size = image_height / float(win_height);
  float x_start = -0.5 * image_width;
  float y_start = -0.5 * image_height;
  RGB_float ret_color;
  Point cur_pixel_pos;
  Vector ray;

  // ray is cast through center of pixel
  cur_pixel_pos.x = x_start + 0.5 * x_grid_size;
  cur_pixel_pos.y = y_start + 0.5 * y_grid_size;
  cur_pixel_pos.z = image_plane;
  // Ray for each pixel
  for (int i=0; i<win_height; i++) {
    for (int j=0; j<win_width; j++) {
      ray = get_vec(eye_pos, cur_pixel_pos);
      normalize(&ray);
      ret_color = recursive_ray_trace(ray, eye_pos, 0, nullptr);
      // Parallel rays can be cast instead using below
      //
      // ray.x = ray.y = 0;
      // ray.z = -1.0;
      // ret_color = recursive_ray_trace(cur_pixel_pos, ray, 1);
      
      if (supersample_on)
      {
        // move in 1/4 pixel increments from middle of pixel
        // take 5 samples (4 corners+middle) and average the colours
        Point ss_pix = cur_pixel_pos;

        ss_pix.x = cur_pixel_pos.x - 0.25*x_grid_size;
        ss_pix.y = cur_pixel_pos.y - 0.25*y_grid_size;
        ray = get_vec(eye_pos, ss_pix);
        normalize(&ray);
        ret_color = clr_add(ret_color, recursive_ray_trace(ray, eye_pos, 0, nullptr));

        ss_pix.x = cur_pixel_pos.x + 0.25*x_grid_size;
        ss_pix.y = cur_pixel_pos.y + 0.25*y_grid_size;
        ray = get_vec(eye_pos, ss_pix);
        normalize(&ray);
        ret_color = clr_add(ret_color, recursive_ray_trace(ray, eye_pos, 0, nullptr));

        ss_pix.x = cur_pixel_pos.x + 0.25*x_grid_size;
        ss_pix.y = cur_pixel_pos.y - 0.25*y_grid_size;
        ray = get_vec(eye_pos, ss_pix);
        normalize(&ray);
        ret_color = clr_add(ret_color, recursive_ray_trace(ray, eye_pos, 0, nullptr));

        ss_pix.x = cur_pixel_pos.x - 0.25*x_grid_size;
        ss_pix.y = cur_pixel_pos.y + 0.25*y_grid_size;
        ray = get_vec(eye_pos, ss_pix);
        normalize(&ray);
        ret_color = clr_add(ret_color, recursive_ray_trace(ray, eye_pos, 0, nullptr));

      }

      frame[i][j][0] = GLfloat(ret_color.r);
      frame[i][j][1] = GLfloat(ret_color.g);
      frame[i][j][2] = GLfloat(ret_color.b);

      cur_pixel_pos.x += x_grid_size;
    }

    cur_pixel_pos.y += y_grid_size;
    cur_pixel_pos.x = x_start;
  }


}
