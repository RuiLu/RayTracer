/*
 * sphere.h
 *
 *  Created on: Feb 16, 2016
 *      Author: Rui
 */

#ifndef SPHERE_H_
#define SPHERE_H_

#include "tools.h"
#include "surface.h"
#include "global.h"
#include <iostream>
#include <cmath>

class Sphere : public Surface {

public:
	Sphere(Point o, double radius);
	double getRadius();
	Point getPoint();
	void setIntersectionPoint(double x, double y, double z);

	virtual bool isIntersected(Ray ray,double t0, double t1, HitRecord *hitRecord);
	virtual Vector getSurfaceNormal();
	virtual Imf::Rgba getColor(Lights &lights, Ray primaryRay);
	virtual Point getIntersectionPoint();
	virtual string getType();
	virtual void setType(string type);

private:
	Point o;
	Point intersectionPoint;
	double radius;
	string type;

};

#endif /* SPHERE_H_ */
