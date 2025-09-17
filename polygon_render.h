#ifndef POLYGON_RENDER_H
#define POLYGON_RENDER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "xmath.h"
#include "xmath3d.h"
#include "xqueue.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { // Polygon structure
    vec_t* points; // Dynamic array of the points of the polygon
    int count;     // Number of the points of the polygon
} Polygon;

typedef struct { // Infill parameters structure
    double scanRadius;     // The radius or width of the scan line
    double overlapRatio;   // Overlap rate between scan lines
    double infillAngle;    // Scan direction angle at which the polygon will be filled
    double offsetDist; // Offset distance to be applied from the polygon edge inwards
} InfillParams;

Polygon offset_polygon(Polygon poly, InfillParams params); //Offsets the polygon

vec_t* generate_horizontal_path(Polygon poly, InfillParams params, int* outputCount); //Creates the infill path horizantally

vec_t* generate_spiral_path(Polygon poly, InfillParams params, int* outputCount); //Creates the infill path spiral (from center to edges)

vec_t* generate_infill_path(Polygon poly, InfillParams params, int* outputCount, int choice); //Creates an infill path for the polygon by using the user's choice

void free_polygon(Polygon* poly); //Releases memory allocated in the polygon structure

#endif