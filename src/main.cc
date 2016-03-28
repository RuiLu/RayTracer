/*
 * main.cc
 *
 *  Created on: Feb 16, 2016
 *      Author: Rui
 */

#include <vector>
#include <iostream>
#include <string>
#include "global.h"
#include "surface.h"
#include "parser.h"
#include "camera.h"
#include "panel.h"

Camera cam;
std::vector< Surface * > mysurfaces;
Lights lights;

int main(int argc, const char * argv[])
{

	if (argc < 3 || argc > 5) {
		std::cout << "usage: <scenefile> <outputfile>" << std::endl;
		return -1;
	}

	clock_t start = clock();

	int r_samples = atoi(argv[3]);
	int s_samples = atoi(argv[4]);
	Parser parser;
	parser.parser(argv[1], cam, mysurfaces, lights);
	Panel *p = new Panel(r_samples, s_samples);
	Render *render;
	render = p->render(cam, mysurfaces, lights);
	Output::renderToFile(argv[2], render);

	// delete pointers
	for (auto s : mysurfaces) {
		delete s;
	}
	for (auto pl : lights.pointLights) {
		delete pl;
	}
	for (auto dl : lights.directionalLights) {
		delete dl;
	}
	for (auto al : lights.ambientLights) {
		delete al;
	}
	for (auto areal : lights.areaLights) {
			delete areal;
		}
	delete p;
	delete render;

	double duration = (clock() - start) / (double) CLOCKS_PER_SEC;
	cout<< duration << "s" << endl;

	return 0;
}


