/*
 * sphere.cc
 *
 *  Created on: Feb 16, 2016
 *      Author: Rui
 */

#include "sphere.h"

using namespace std;

Sphere::Sphere(Point o, double radius) : Surface(o) {
	this->o = o;
	this->radius = radius;
	this->bbox.minP = Point(o.getX() - radius, o.getY() - radius, o.getZ() - radius);
	this->bbox.maxP = Point(o.getX() + radius, o.getY() + radius, o.getZ() + radius);
}

double Sphere::getRadius() {
	return radius;
}

Point Sphere::getPoint(){
	return o;
}

Point Sphere::getIntersectionPoint() {
	return this->intersectionPoint;
}

void Sphere::setIntersectionPoint(double x, double y, double z) {
	this->intersectionPoint.setPoint(x, y, z);
}

bool Sphere::isIntersected(Ray primaryRay, double t0, double t1, HitRecord *hitRecord) {
//
//	if (!bbox.hitBox(primaryRay, t0, t1, hitRecord)) {
//		return false;
//	}

	std::vector<double> hits;

    double x0 = primaryRay.origin.getX();
    double y0 = primaryRay.origin.getY();
    double z0 = primaryRay.origin.getZ();

    double xd = primaryRay.dir.getX();
    double yd = primaryRay.dir.getY();
    double zd = primaryRay.dir.getZ();

    double xc = this->o.getX();
    double yc = this->o.getY();
    double zc = this->o.getZ();

    double radius = this->radius;

    double A = primaryRay.dir * primaryRay.dir;
    double B = 2 * (xd * (x0 - xc) + yd * (y0 - yc) + zd * (z0 - zc));
    double C = (x0 - xc) * (x0 - xc) + (y0 - yc) * (y0 - yc) + (z0 - zc) * (z0 - zc) - radius * radius;

    double discriminant = B * B - 4 * C;

    if (discriminant < 0.0) {
    	return false;
    } else if (discriminant == 0.0) {
    	double t = (-1.0 * B) / (2 * A);
    	// The only intersection point is behind or on the camera
    	if (t <= 0.001) {
    		return false;
    	}

//    	std::cout<<t<<std::endl;

    	double ip_x = x0 + t * xd;
    	double ip_y = y0 + t * yd;
    	double ip_z = z0 + t * zd;

//    	Point ip = Point(ip_x, ip_y, ip_z);
//    	hitRecord->intersectionPoint = ip;

    	setIntersectionPoint(ip_x, ip_y, ip_z);
    	hits.push_back(t);
    	this->getHitRecord(hitRecord, hits, primaryRay);

    	return true;

    } else {
    	double t0 = (-1.0 * B - sqrt(discriminant)) / (2 * A);
    	double t1 = (-1.0 * B + sqrt(discriminant)) / (2 * A);

    	// Both intersection points are behind or on the camera
    	if (t1 <= 0.) return false;

    	// t0 is behind the camera and t1 is in front of the camera,
    	// let's return false first
    	if (t0 <= 0.001) return false;

    	//Both intersection points are in front of the camera
    	double ip_x = x0 + t0 * xd;
    	double ip_y = y0 + t0 * yd;
    	double ip_z = z0 + t0 * zd;

    	setIntersectionPoint(ip_x, ip_y, ip_z);
    	hits.push_back(t0);
    	this->getHitRecord(hitRecord, hits, primaryRay);
    	return true;
    }
}

Vector Sphere::getSurfaceNormal() {

	double ip_x = this->intersectionPoint.getX();
	double ip_y = this->intersectionPoint.getY();
	double ip_z = this->intersectionPoint.getZ();

	double xc = this->o.getX();
	double yc = this->o.getY();
	double zc = this->o.getZ();

	Vector surface_normal = Vector(ip_x - xc, ip_y - yc, ip_z - zc);

//	double res = (xc - ip_x) * (xc - ip_x) + (yc - ip_y) * (yc - ip_y) + (zc - ip_z) * (zc - ip_z);
//	res = res - this->radius * this->radius;
//	std::cout<<res<<std::endl;

	surface_normal.normalize();
	return surface_normal;
}


Imf::Rgba Sphere::getColor(Lights &lights, Ray primaryRay) {
	Imf::Rgba c = Imf::Rgba(0,0,0,1);
	return c;
}

string Sphere::getType() {
	return this->type;
}

void Sphere::setType(string type) {
	this->type = type;
}
