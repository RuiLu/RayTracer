/*
 * surface.cc
 *
 *  Created on: Feb 16, 2016
 *      Author: Rui
 */

#include "surface.h"

Surface::Surface() {
	this->isLeaf = false;
	this->phong_exponent = 0.0;
}

Surface::Surface(Point center) {
	this->center = center;
	this->isLeaf = false;
	this->phong_exponent = 0.0;
}

string Surface::getType() {
	std::cout << "error: surface::getType() should not be called - it's virtual!"
			    << std::endl;
	return "";
}

void Surface::setType(string type) {
	std::cout << "error: surface::setType() should not be called - it's virtual!"
		    << std::endl;
}

bool Surface::isIntersected(Ray ray, double t0, double t1, HitRecord *hitRecord) {

	std::cout << "error: surface::intersect should not be called - it's virtual!"
	    << std::endl;

	return false;
}

Vector Surface::getSurfaceNormal() {
	std::cout << "error: surface::getSurfaceNormal should not be called - it's virtual!"
		    << std::endl;
	return Vector(0, 0, 0);
}

Vector Surface::getl(PointLight pl) {
	std::cout << "error: surface::getl should not be called - it's virtual!"
		    << std::endl;
	return Vector(0, 0, 0);
}

Imf::Rgba Surface::getColor(Lights &lights, Ray ray) {
	std::cout << "error: surface::getColor should not be called - it's virtual!"
			    << std::endl;
	return Imf::Rgba (0,0,0,1);
}

void Surface::setMaterial(Material material)
{
	this->diffuse = material.diffuse;
	this->specular = material.specular;
	this->ideal_specular = material.ideal_specular;
	this->phong_exponent = material.phong_exponent;
}

Point Surface::getIntersectionPoint() {
	std::cout << "error: surface::getIntersectionPoint() should not be called - it's virtual!"
				    << std::endl;
	return this->intersectionPoint;
}

Imf::Rgba Surface::getDiffuse() {
	return this->diffuse;
}

Imf::Rgba Surface::getSpecular() {
	return this->specular;
}

Imf::Rgba Surface::getIdealSpecular() {
	return this->ideal_specular;
}

double Surface::getPhongExp() {
	return this->phong_exponent;
}

void Surface::getHitRecord(HitRecord *rec, std::vector<double> hits, Ray ray) {
	if (hits.size() < 1) {
		rec->hitNums = 0;
		return;
	}

	for (int i = 0; i < hits.size(); i++) {
		rec->hits.push_back(hits[i]);
		rec->hitNums++;
	}

	rec->diffuse = this->getDiffuse();
	rec->specular = this->getSpecular();
	rec->ideal = this->getIdealSpecular();
	rec->normal = this->getSurfaceNormal();
	rec->phong_exp = this->getPhongExp();
}



/*
 * Material
 */

Material::Material(Imf::Rgba diffuse,
				   Imf::Rgba specular,
				   Imf::Rgba ideal_specular,
				   double phong_exponent)
{
	this->diffuse = diffuse;
	this->specular = specular;
	this->ideal_specular = ideal_specular;
	this->phong_exponent = phong_exponent;
}

/*
 * Plane
 */

Plane::Plane(Vector normal, double scalar_value) : Surface() {
	normal.normalize();
	this->normal = normal;
	this->scalar_value = scalar_value;

	double infinity = std::numeric_limits<double>::infinity();
	this->bbox.minP = Point(-infinity, -infinity, -infinity);
	this->bbox.maxP = Point(infinity, infinity, infinity);

}

Vector Plane::getSurfaceNormal() {
	return this->normal;
}

double Plane::getScalarValue() {
	return this->scalar_value;
}

bool Plane::isIntersected(Ray ray, double t0, double t1, HitRecord *hitRecord) {
//
//	if (!bbox.hitBox(ray, t0, t1, hitRecord)) {
//		return false;
//	}

	std::vector<double> hits;
	Point r0 = ray.origin;
	Vector rd = ray.dir;
	Vector pn = this->normal;
	double scalar_value = this->scalar_value;

	double vd = pn * rd;
	double v0 = -1.0 * (r0 * pn + scalar_value);
	double t = v0 /vd;

	if (t <= 0.001) return false;

	double ip_x = r0.getX() + t * rd.getX();
	double ip_y = r0.getY() + t * rd.getY();
	double ip_z = r0.getZ() + t * rd.getZ();

	setIntersectionPoint(ip_x, ip_y, ip_z);

	hits.push_back(t);
	this->getHitRecord(hitRecord, hits, ray);

	return true;

}

void Plane::setIntersectionPoint(double x, double y, double z) {
	this->intersectionPoint.setPoint(x, y, z);
}

Imf::Rgba Plane::getColor(Lights &lights, Ray ray)
{

		Imf::Rgba c = Imf::Rgba(0,0,0,1);

		Vector n = this->normal;

		for (int y = 0; y<lights.pointLights.size(); y++) {
			// the normalized vector from the closest intersection point
			// to the light -> l
			PointLight *pl = lights.pointLights[y];
;
			Vector l = pl->getPoint() - this->intersectionPoint;
			l.normalize();

			// the normalized vector that is opposite to camera primary ray -> v
			Vector v = -ray.dir;
			v.normalize();

			// diffuse shading
			double lDotNormal = n * l;
			lDotNormal = lDotNormal < 0.0 ? 0.0 : lDotNormal;

			c.r += this->getDiffuse().r * pl->getIntensity_p().r * lDotNormal;
			c.g += this->getDiffuse().g * pl->getIntensity_p().g * lDotNormal;
			c.b += this->getDiffuse().b * pl->getIntensity_p().b * lDotNormal;

			// Blinn-Phong shading
			Vector h = v + l;
			h = h / (h.length());

			double exponent_res = pow (n * h, this->getPhongExp());
			exponent_res = exponent_res < 0.0 ? 0.0 : exponent_res;

			c.r += this->getSpecular().r * pl->getIntensity_p().r * exponent_res;
			c.g += this->getSpecular().g * pl->getIntensity_p().g * exponent_res;
			c.b += this->getSpecular().b * pl->getIntensity_p().b * exponent_res;

		}

		return c;
}


Point Plane::getIntersectionPoint() {
	return this->intersectionPoint;
}

string Plane::getType() {
	return this->type;
}

void Plane::setType(string type) {
	this->type = type;
}

/*
 * Bounding box
 */
BoundingBox::BoundingBox() {
	this->minP = Point(0, 0, 0);
	this->maxP = Point(0, 0, 0);
}

BoundingBox::BoundingBox(Point minP, Point maxP) {
	this->minP = minP;
	this->maxP = maxP;
}

BoundingBox BoundingBox::combine(BoundingBox &left, BoundingBox &right) {
	double minX, minY, minZ;
	double maxX, maxY, maxZ;

	if (left.minP.getX() < right.minP.getX()) {
		minX = left.minP.getX();
	} else {
		minX = right.minP.getX();
	}
	if (left.minP.getY() < right.minP.getY()) {
		minY = left.minP.getY();
	} else {
		minY = right.minP.getY();
	}
	if (left.minP.getZ() < right.minP.getZ()) {
		minZ = left.minP.getZ();
	} else {
		minZ = right.minP.getZ();
	}
	if (left.maxP.getX() > right.maxP.getX()) {
		maxX = left.maxP.getX();
	} else {
		maxX = right.maxP.getX();
	}
	if (left.maxP.getY() > right.maxP.getY()) {
		maxY = left.maxP.getY();
	} else {
		maxY = right.maxP.getY();
	}
	if (left.maxP.getZ() > right.maxP.getZ()) {
		maxZ = left.maxP.getZ();
	} else {
		maxZ = right.maxP.getZ();
	}

	Point newMinP = Point(minX, minY, minZ);
	Point newMaxP = Point(maxX, maxY, maxZ);

	BoundingBox bbox(newMinP, newMaxP);

	return bbox;
}

BoundingBox BoundingBox::getBoundingBox(std::vector< Surface * > &surfaces) {

	double minX = 0.0, minY = 0.0, minZ = 0.0;
	double maxX = 0.0, maxY = 0.0, maxZ = 0.0;

	if (surfaces.size() > 0) {
		minX = surfaces[0]->bbox.minP.getX();
		minY = surfaces[0]->bbox.minP.getY();
		minZ = surfaces[0]->bbox.minP.getZ();
		maxX = surfaces[0]->bbox.maxP.getX();
		maxY = surfaces[0]->bbox.maxP.getY();
		maxZ = surfaces[0]->bbox.maxP.getZ();

		for (int i = 1; i < surfaces.size(); i++) {
			// get the minimum point
			if (surfaces[i]->bbox.minP.getX() < minX) minX = surfaces[i]->bbox.minP.getX();
			if (surfaces[i]->bbox.minP.getY() < minY) minY = surfaces[i]->bbox.minP.getY();
			if (surfaces[i]->bbox.minP.getZ() < minZ) minZ = surfaces[i]->bbox.minP.getZ();
			// get the maximal point
			if (surfaces[i]->bbox.maxP.getX() > maxX) maxX = surfaces[i]->bbox.maxP.getX();
			if (surfaces[i]->bbox.maxP.getY() > maxY) maxY = surfaces[i]->bbox.maxP.getY();
			if (surfaces[i]->bbox.maxP.getZ() > maxZ) maxZ = surfaces[i]->bbox.maxP.getZ();
		}
	}

	Point newMinP = Point(minX, minY, minZ);
	Point newMaxP = Point(maxX, maxY, maxZ);
	BoundingBox bbox(newMinP, newMaxP);
	return bbox;
}

bool BoundingBox::hitBox(Ray ray, double t0, double t1, HitRecord *hitRecord) {

	double a;	// 1 / direction
	double txmin, txmax, tymin, tymax, tzmin, tzmax;
	double t_close, t_far;

	a = 1 / ray.dir.getX();
	if (a >= 0.) {
		txmin = a * (this->minP.getX() - ray.origin.getX());
		txmax = a * (this->maxP.getX() - ray.origin.getX());
	} else {
		txmax = a * (this->minP.getX() - ray.origin.getX());
		txmin = a * (this->maxP.getX() - ray.origin.getX());
	}

	a = 1 / ray.dir.getY();
	if (a >= 0.) {
		tymin = a * (this->minP.getY() - ray.origin.getY());
		tymax = a * (this->maxP.getY() - ray.origin.getY());
	} else {
		tymax = a * (this->minP.getY() - ray.origin.getY());
		tymin = a * (this->maxP.getY() - ray.origin.getY());
	}

	a = 1 / ray.dir.getZ();
	if (a >= 0.) {
		tzmin = a * (this->minP.getZ() - ray.origin.getZ());
		tzmax = a * (this->maxP.getZ() - ray.origin.getZ());
	} else {
		tzmax = a * (this->minP.getZ() - ray.origin.getZ());
		tzmin = a * (this->maxP.getZ() - ray.origin.getZ());
	}

	t_close = std::max(tzmin, std::max(txmin, tymin));
	t_far = std::min(tzmax, std::min(txmax, tymax));

	if (t_close > t_far) return false;

	if ((t_close >= t0 && t_close <= t1) || (t_far >= t0 && t_far <= t1)) {
		hitRecord->hits.push_back(t_close);

		if (t_close == tzmin) {
			hitRecord->normal = Vector(0, 0, 1);
		} else if (t_close == tymin) {
			hitRecord->normal = Vector(0, 1, 0);
		} else if (t_close == txmin) {
			hitRecord->normal = Vector(1, 0, 0);
		}

		return true;
	}
	else if (t_close < t0 && t_far > t1) {
		hitRecord->hits.push_back(t_close);

		if (t_close == tzmin) {
			hitRecord->normal = Vector(0, 0, 1);
		} else if (t_close == tymin) {
			hitRecord->normal = Vector(0, 1, 0);
		} else if (t_close == txmin) {
			hitRecord->normal = Vector(1, 0, 0);
		}

		return true;
	}
	else {
		return false;
	}
}

/*
 * BVH node
 */

BvhNode::BvhNode(std::vector< Surface * > &surfaces, int axis) {
	this->surfaces = surfaces;
	this->axis = axis;
	this->left = NULL;
	this->right = NULL;
}

void BvhNode::create() {

	int N = surfaces.size();
//	std::cout<<N<<std::endl;
	if (N == 1) {
		left = surfaces[0];
		left->isLeaf = true;
		right = NULL;
		bbox = surfaces[0]->bbox;
	} else if (N == 2) {
		left = surfaces[0];
		left->isLeaf = true;
		right = surfaces[1];
		right->isLeaf = true;
		bbox = BoundingBox::combine(surfaces[0]->bbox, surfaces[1]->bbox);
	}  else {
		BvhNode *nodeToExpand;
		this->split();
		nodeToExpand = new BvhNode(leftSurfaces, (axis + 1) % 3);
		nodeToExpand->create();
		left = nodeToExpand;

		nodeToExpand = new BvhNode(rightSurfaces, (axis + 1) % 3);
		nodeToExpand->create();
		right = nodeToExpand;
		bbox = BoundingBox::combine(left->bbox, right->bbox);
	}
}

void BvhNode::split() {
//	BoundingBox bbox = BoundingBox::getBoundingBox(surfaces);
	surfaces = Sort::quickSort(surfaces, axis);

	int middle = surfaces.size() / 2;

	this->leftSurfaces.clear();
	for (int i = 0; i < middle; i++) {
		Surface *s = surfaces[i];
		this->leftSurfaces.push_back(s);
	}
	this->rightSurfaces.clear();
	for (int i = middle; i < surfaces.size(); i++) {
		Surface *s = surfaces[i];
		this->rightSurfaces.push_back(s);
	}

}

bool BvhNode::isIntersected(Ray ray, double t0, double t1, HitRecord *hitRecord) {
	if (this->bbox.hitBox(ray, t0, t1, hitRecord)) {
		HitRecord lrec, rrec;
		lrec.hitNums = 0;
		rrec.hitNums = 0;
		bool hitLeft = false;
		bool hitRight = false;
		if (this->left != NULL) hitLeft = this->left->isIntersected(ray, t0, t1, &lrec);
		if (this->right != NULL) hitRight = this->right->isIntersected(ray, t0, t1, &rrec);

		if (hitLeft && hitRight) {
			if (lrec.hits[0] <= rrec.hits[0]) {
				*hitRecord = lrec;
			} else {
				*hitRecord = rrec;
			}
			return true;
		} else if (hitLeft) {
			*hitRecord = lrec;
			return true;
		} else if (hitRight) {
			*hitRecord = rrec;
			return true;
		} else {
			return false;
		}
	} else {
		return false;
	}
}


/*
 * Sort
 */
std::vector< Surface * > Sort::quickSort(std::vector< Surface * > &surfaces, int axis) {

	std::vector< Surface * > left, right, res;
	if (surfaces.size() <= 1) return surfaces;

	int middle = surfaces.size() / 2;
	Surface *pivot = surfaces[middle];
	surfaces.erase(surfaces.begin() + middle);

	for (int i = 0; i < surfaces.size(); i++) {
		if (axis == 0) {
			if (surfaces[i]->center.getX() <= pivot->center.getX()) {
				left.push_back(surfaces[i]);
			} else {
				right.push_back(surfaces[i]);
			}
		} else if (axis == 1) {
			if (surfaces[i]->center.getY() <= pivot->center.getY()) {
				left.push_back(surfaces[i]);
			} else {
				right.push_back(surfaces[i]);
			}
		} else if (axis == 2) {
			if (surfaces[i]->center.getZ() <= pivot->center.getZ()) {
				left.push_back(surfaces[i]);
			} else {
				right.push_back(surfaces[i]);
			}
		}
	}
	std::vector< Surface * > sortedL = Sort::quickSort(left, axis);
	std::vector< Surface * > sortedR = Sort::quickSort(right, axis);
	for (int i = 0; i < sortedL.size(); i++) {
		Surface *s = sortedL[i];
		res.push_back(s);
	}
	res.push_back(pivot);

	for (int i = 0; i < sortedR.size(); i++) {
		Surface *s = sortedR[i];
		res.push_back(s);
	}

	return res;

}

Vector BvhNode::getSurfaceNormal() {
	return Vector(0, 0, 0);
}
