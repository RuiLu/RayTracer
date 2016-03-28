/*
 * light.h
 *
 *  Created on: Feb 17, 2016
 *      Author: Rui
 */

#ifndef LIGHT_H_
#define LIGHT_H_


#include <iostream>
#include <vector>

#include "tools.h"

class PointLight {

public:
	PointLight(Point pos, Imf::Rgba intensity_p);
	Point getPoint();
	Imf::Rgba getIntensity_p();

private:
	Point pos;
	Imf::Rgba intensity_p;

};

class DirectionalLight {

public:
	DirectionalLight(Vector dir, Imf::Rgba intensity_d);
	Vector getDir();
	Imf::Rgba getIntensity_d();

private:
	Vector dir;
	Imf::Rgba intensity_d;
};

class AmbientLight {

public:
	AmbientLight(Imf::Rgba intensity_a);
	Imf::Rgba getIntensity_a();

private:
	Imf::Rgba intensity_a;

};

class AreaLight {
public:
	AreaLight(Point center, Vector dir, Vector u_dir, double len, Imf::Rgba intensity_area);
	Imf::Rgba getIntensity_area();
	Point getCenter();
	Vector getDir();
	Vector getUDir();
	double getLen();


private:
	Imf::Rgba intensity_area;
	Point center;
	Vector dir;
	Vector u_dir;
	double len;

};

typedef struct lights {
	std::vector< PointLight * > pointLights;
	std::vector< DirectionalLight * > directionalLights;
	std::vector< AmbientLight * > ambientLights;
	std::vector< AreaLight * > areaLights;
} Lights;

#endif /* LIGHT_H_ */
