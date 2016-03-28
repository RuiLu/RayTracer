/*
 * surface.h
 *
 *  Created on: Feb 16, 2016
 *      Author: Rui
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include <string>
#include <iostream>
#include <limits>
#include "tools.h"
#include "light.h"

using namespace std;
/*
 * I found the idea of demo code in this part is neat and elegant,
 * So I learn this from the demo code.
 */

class Material;
class Surface;
class BoundingBox;

typedef struct hit_record_ {

	int hitNums;
	std::vector<double> hits;

	Point intersectionPoint;
	double phong_exp;
	Vector normal;
	Imf::Rgba ambient;
	Imf::Rgba diffuse;
	Imf::Rgba specular;
	Imf::Rgba ideal;


} HitRecord;

class BoundingBox {
public:
	Point minP;
	Point maxP;

	BoundingBox();
	BoundingBox(Point, Point);
	static BoundingBox combine(BoundingBox &, BoundingBox &);
	static BoundingBox getBoundingBox(std::vector<Surface *> &);
	bool hitBox(Ray, double, double, HitRecord *);
};

class Surface {

public:
	void setMaterial(Material material);
	Imf::Rgba getDiffuse();
	Imf::Rgba getSpecular();
	Imf::Rgba getIdealSpecular();
	double getPhongExp();

	BoundingBox bbox;
	void getHitRecord(HitRecord *, std::vector<double>, Ray );

	Point center;
	Surface();
	Surface(Point center);

	bool isLeaf;

	virtual bool isIntersected(Ray ,double, double, HitRecord *);
	virtual Vector getSurfaceNormal();
	virtual Vector getl(PointLight );
	virtual Imf::Rgba getColor(Lights &, Ray);
	virtual Point getIntersectionPoint();
	virtual string getType();
	virtual void setType(string type);

private:
	Imf::Rgba diffuse;
	Imf::Rgba specular;
	Imf::Rgba ideal_specular;
	double phong_exponent;
	Point intersectionPoint;
	string type;

};

class Material {

public:

	Material(Imf::Rgba diffuse,
			 Imf::Rgba specular,
			 Imf::Rgba ideal_specular,
			 double phong_exponent);

	Imf::Rgba diffuse;
	Imf::Rgba specular;
	Imf::Rgba ideal_specular;
	double phong_exponent;

};

class Plane : public Surface {

public:
	Plane(Vector normal, double scalar_value);
	void setIntersectionPoint(double x, double y, double z);
	double getScalarValue();

	virtual Vector getSurfaceNormal();
	virtual Imf::Rgba getColor(Lights &, Ray);
	virtual bool isIntersected(Ray , double, double, HitRecord *);
	virtual Point getIntersectionPoint();
	virtual string getType();
	virtual void setType(string type);

private:
	Vector normal;
	double scalar_value;
	Point intersectionPoint;
	string type;

};

class BvhNode : public Surface {
public:
	Surface *left;
	Surface *right;
	BvhNode(std::vector<Surface *> &surfaces, int axis);
	void create();
	void split();
	virtual bool isIntersected(Ray, double, double, HitRecord *);
	virtual Vector getSurfaceNormal();

private:
	std::vector< Surface * > surfaces;
	std::vector< Surface * > leftSurfaces;
	std::vector< Surface * > rightSurfaces;
	int axis;

};

class Sort {
public:
	static std::vector< Surface * > quickSort(std::vector< Surface * > &surfaces, int axis);
};

#endif /* SURFACE_H_ */
