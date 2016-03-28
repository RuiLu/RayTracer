/*
 * camera.h
 *
 *  Created on: Feb 16, 2016
 *      Author: Rui
 */

#ifndef CAMERA_H_
#define CAMERA_H_

#include "global.h"
#include "tools.h"

class Camera {

public:

	void initCamera(Point eye, Vector dir, double lr, double tb, double fl, int width, int height);
	void createBasis(Vector &dir);

	Point eye;
	Vector u, v, w;

	double focal_length;
	double l, r, t, b;

	int width, height;

	static Vector UP;

private:
};

#endif /* CAMERA_H_ */
