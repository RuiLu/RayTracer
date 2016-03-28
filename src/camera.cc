/*
 * camera.cc
 *
 *  Created on: Feb 16, 2016
 *      Author: Rui
 */

#include "camera.h"

using namespace std;

Vector Camera::UP (0, 1, 0);

void Camera::initCamera(Point eye, Vector direction, double lr, double tb, double fl, int width, int height) {
    this->eye = eye;
    this->width = width;
    this->height = height;
    this->focal_length = fl;
    this->r = lr / 2.0;
    this->l = -1. * r;
    this->t = tb / 2.0;
    this->b = -1. * t;
    createBasis(direction);
}

void Camera::createBasis(Vector &dir) {
	dir.normalize();
	w = -dir;
	w.normalize();

	if (1 - fabs (w.getY()) < .0001) {
		u.setVector(1.0, 0.0, 0.0);
		v = u^w;
	} else {
		u = UP^w;
		v = u^w;
	}

//	u = dir^UP;
//	v = u^dir;
//	w = -dir;
	u.normalize();
	v.normalize();
//	w.normalize();
}
