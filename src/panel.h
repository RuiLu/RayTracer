/*
 * panel.h
 *
 *  Created on: Feb 17, 2016
 *      Author: Rui
 */

#ifndef PANEL_H_
#define PANEL_H_

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <stdio.h>
#include <vector>
#include <limits>
#include <cmath>
#include "camera.h"
#include "surface.h"
#include "sphere.h"
#include "tools.h"
#include "light.h"

class Panel {

public:
	Panel(int, int);
    Render *render(Camera &cam,
    			   std::vector< Surface * > &mysurfaces,
				   Lights &lights);
    Vector getColorBvh(Ray ray, double min_t, double max_t, int recurse_limit);
    void deleteTree(BvhNode *root);

private:

    static Imf::Rgba black;
    BvhNode *bvhTree;
    int r_samples;
    int s_samples;

};



#endif /* PANEL_H_ */
