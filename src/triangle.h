/*
 * triangle.h
 *
 *  Created on: Feb 22, 2016
 *      Author: Rui
 */

#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include <cmath>
#include <iostream>
#include "global.h"
#include "tools.h"
#include "surface.h"

class Triangle : public Surface{

public:
	Triangle(Point v0, Point v1, Point v2);
	std::vector<Point> vertices;
	virtual bool isIntersected(Ray ray, double t0, double t1, HitRecord *hitRecord);
	virtual Imf::Rgba getColor(Lights &lights, Ray primaryRay);
	virtual Vector getSurfaceNormal();
	void setIntersectionPoint(double x, double y, double z);
	Vector getNormal(Point v0, Point v1, Point v2);
	bool isInside(Point ip_p);
	virtual Point getIntersectionPoint();
	virtual string getType();
	virtual void setType(string type);

private:
	Point v0;
	Point v1;
	Point v2;
	Vector normal;
	Point intersectionPoint;
	string type;

};


#endif /* TRIANGLE_H_ */
