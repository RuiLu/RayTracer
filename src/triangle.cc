/*
 * triangle.cc
 *
 *  Created on: Feb 22, 2016
 *      Author: Rui
 */

#include "triangle.h"


Triangle::Triangle(Point v0, Point v1, Point v2) : Surface() {


	this->v0 = v0;
	this->v1 = v1;
	this->v2 = v2;
	this->normal = getNormal(v0, v1, v2);

	vertices.clear();
	vertices.push_back(v0);
	vertices.push_back(v1);
	vertices.push_back(v2);

	double centerX = (v0.getX() + v1.getX() + v2.getX()) / 3.0;
	double centerY = (v0.getY() + v1.getY() + v2.getY()) / 3.0;
	double centerZ = (v0.getZ() + v1.getZ() + v2.getZ()) / 3.0;
	this->center = Point(centerX, centerY, centerZ);

	double minX = min(vertices[0].getX(), min(vertices[1].getX(), vertices[2].getX()));
	double minY = min(vertices[0].getY(), min(vertices[1].getY(), vertices[2].getY()));
	double minZ = min(vertices[0].getZ(), min(vertices[1].getZ(), vertices[2].getZ()));
	double maxX = max(vertices[0].getX(), max(vertices[1].getX(), vertices[2].getX()));
	double maxY = max(vertices[0].getY(), max(vertices[1].getY(), vertices[2].getY()));
	double maxZ = max(vertices[0].getZ(), max(vertices[1].getZ(), vertices[2].getZ()));

	this->bbox.minP = Point(minX, minY, minZ);
	this->bbox.maxP = Point(maxX, maxY, maxZ);

}

Vector Triangle::getNormal(Point v0, Point v1, Point v2) {
	Vector A = v1 - v0;
	Vector B = v2 - v0;
	double cx = A.getY() * B.getZ() - A.getZ() * B.getY();
	double cy = A.getZ() * B.getX() - A.getX() * B.getZ();
	double cz = A.getX() * B.getY() - A.getY() * B.getX();

	Vector c = Vector(cx, cy, cz);
	c.normalize();
	return c;
}

Vector Triangle::getSurfaceNormal() {
	return this->normal;
}

Point Triangle::getIntersectionPoint() {
	return this->intersectionPoint;
}

bool Triangle::isIntersected(Ray primaryRay, double t0, double t1, HitRecord *hitRecord) {

	std::vector<double> hits;
	Vector diff = this->v0 - primaryRay.origin;
	double up = diff * this->normal;
	double down = primaryRay.dir * this->normal;
	double t = up / down;
	double intersected = this->normal * primaryRay.dir;

	if (intersected == 0.) {
//		std::cout<<"This ray is parallel with this triangle."<<std::endl;
		return false;
	}

	if (t <= 0.001) {
//		std::cout<<"this pixel of triangle is behind the camera."<<std::endl;
		return false;
	}

	double ip_x = primaryRay.origin.getX() + t * primaryRay.dir.getX();
	double ip_y = primaryRay.origin.getY() + t * primaryRay.dir.getY();
	double ip_z = primaryRay.origin.getZ() + t * primaryRay.dir.getZ();
	Point ip_p = Point(ip_x, ip_y, ip_z);

	if (!isInside(ip_p)) {
		return false;
	}

	setIntersectionPoint(ip_x, ip_y, ip_z);

	hits.push_back(t);
	this->getHitRecord(hitRecord, hits, primaryRay);

	return true;

}

bool Triangle::isInside(Point ip_p) {

	Vector tmp1, tmp2, tmp3;
	double tmp4;

	Point a = this->v0;
	Point b = this->v1;
	Point c = this->v2;

	tmp1 = b - a;
	tmp2 = ip_p - a;
	tmp3 = crossProduct(tmp1, tmp2);
	tmp4 = tmp3 * this->normal;
	if (tmp4 < 0) {
//		std::cout<<"Intersection point is not in the triangle"<<std::endl;
		return false;
	}

	tmp1 = c - b;
	tmp2 = ip_p - b;
	tmp3 = crossProduct(tmp1, tmp2);
	tmp4 = tmp3 * this->normal;
	if (tmp4 < 0) {
//		std::cout<<"Intersection point is not in the triangle"<<std::endl;
		return false;
	}

	tmp1 = a - c;
	tmp2 = ip_p - c;
	tmp3 = crossProduct(tmp1, tmp2);
	tmp4 = tmp3 * this->normal;
	if (tmp4 < 0) {
//		std::cout<<"Intersection point is not in the triangle"<<std::endl;
		return false;
	}

	return true;
}

void Triangle::setIntersectionPoint(double x, double y, double z) {
	this->intersectionPoint.setPoint(x, y, z);
}

Imf::Rgba Triangle::getColor(Lights &lights, Ray primaryRay) {

	Imf::Rgba c = Imf::Rgba(0,0,0,1);

	Vector n = this->normal;

	double lamr = 0.0, lamg = 0.0, lamb = 0.0;
	double bpr = 0.0, bpg = 0.0, bpb = 0.0;

	for (int y = 0; y < lights.pointLights.size(); y++) {

		PointLight *pl = lights.pointLights[y];

		Vector l = pl->getPoint() - this->intersectionPoint;
		l.normalize();

		// the normalized vector that is opposite to camera primary ray -> v
		Vector v = -primaryRay.dir;
		v.normalize();

		// diffuse shading
		double lDotNormal = n * l;

		lDotNormal = lDotNormal < 0.0 ? 0.0 : lDotNormal;

		lamr += this->getDiffuse().r * pl->getIntensity_p().r * lDotNormal;
		lamg += this->getDiffuse().g * pl->getIntensity_p().g * lDotNormal;
		lamb += this->getDiffuse().b * pl->getIntensity_p().b * lDotNormal;

		// Blinn-Phong shading
		Vector h = v + l;

		if (h.length() > 0.0) {

			h.normalize();

			double exponent_res = n * h;
			exponent_res = exponent_res < 0.0 ? 0.0 : exponent_res;
			exponent_res = pow (exponent_res, this->getPhongExp());

			bpr += this->getSpecular().r * pl->getIntensity_p().r * exponent_res;
			bpg += this->getSpecular().g * pl->getIntensity_p().g * exponent_res;
			bpb += this->getSpecular().b * pl->getIntensity_p().b * exponent_res;

		}
	}

	c.r = lamr + bpr;
	c.g = lamg + bpg;
	c.b = lamb + bpb;
	c.a = 1.0;

	return c;
}

string Triangle::getType() {
	return this->type;
}

void Triangle::setType(string type) {
	this->type = type;
}
