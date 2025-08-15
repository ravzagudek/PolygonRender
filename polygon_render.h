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
    double scan_radius;     // The radius or width of the scan line
    double overlap_ratio;   // Overlap rate between scan lines
    double infill_angle;    // Scan direction angle at which the polygon will be filled
    double offset_distance; // Offset distance to be applied from the polygon edge inwards
} InfillParams;

vec_t intersect(vec_t A1, vec_t B1, vec_t A2, vec_t B2); //Finds the intersection point of two line segments

vec_t find_normal_unit_vector(vec_t p1, vec_t p2); //Finds the normal unit vector between two given points

Polygon offset_polygon(Polygon poly, InfillParams params); //Offsets the polygon

double calc_poly_area(Polygon poly); //Calculates the polygon area by using Shoelace Formula

void reverse_poly_points(Polygon* poly); // Reverses the points of the given polygon

vec_t* queue_to_array(xqueue_t* q, int* count); // Turns the point queue into an array

Polygon rotate_polygon(Polygon poly, double angle); //Rotates the polygon by the given angle

vec_t* generate_horizontal_path(Polygon poly, InfillParams params, int* output_count); //Creates the infill path horizantally

vec_t* generate_spiral_path(Polygon poly, InfillParams params, int* output_count); //Creates the infill path spiral (from center to edges)

vec_t* generate_infill_path(Polygon poly, InfillParams params, int* output_count, int choice); //Creates an infill path for the polygon by using the user's choice

void export_to_bmp(Polygon poly, vec_t* path, int path_count, const char* filename); //Exports the created polygon and fill path to a .bmp file

void free_polygon(Polygon* poly); //Releases memory allocated in the polygon structure

#endif