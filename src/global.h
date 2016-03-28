/*
 * global.h
 *
 *  Created on: Feb 16, 2016
 *      Author: Rui
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <vector>
#include "tools.h"
#include "camera.h"
#include "surface.h"
#include "sphere.h"
#include "light.h"

class Camera;
class Surface;

extern Camera cam;
extern std::vector<Surface * > mysurfaces;
extern Lights lights;

#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2

#endif /* GLOBAL_H_ */
