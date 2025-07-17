#ifndef POLYGON_RENDER_H
#define POLYGON_RENDER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <xmath.h>
#include <xmath3d.h>

typedef struct { // Polygon structure
    vec_t* points; // Dynamic array of the points of the polygon
    int count;     // Number of the points of the polygon
} Polygon;

typedef struct InfillPointNode { // linked list node structure, used to dynamically store the points of the infill path
    vec_t point;                // Point hidden in the node
    struct InfillPointNode* next; // Pointer to the next node in the linked list
} InfillPointNode;

typedef struct { // Infill parameters structure
    double scan_radius;     // The radius or width of the scan line
    double overlap_ratio;   // Overlap rate between scan lines
    double infill_angle;    // Scan direction angle at which the polygon will be filled
    double offset_distance; // Offset distance to be applied from the polygon edge inwards
} InfillParams;

// Rotates the polygon by the given angle
Polygon rotate_polygon(Polygon poly, double angle);

// Finds the intersection point of two line segments
vec_t intersect(vec_t A1, vec_t B1, vec_t A2, vec_t B2);

// Offsets the polygon
Polygon offset_polygon(Polygon poly, InfillParams params);

// Creates an infill path for the polygon
vec_t* generate_infill_path(Polygon poly, InfillParams params, int* output_count);

// Creates a new linked list node
InfillPointNode* create_node(vec_t p);

// Adds a new node to the list
void add_to_list(InfillPointNode** head, vec_t p);

// Frees the linked list
void free_list(InfillPointNode* head);

// Turns the linked list into a Point array
vec_t* list_to_array(InfillPointNode* head, int* count);

// Exports the created polygon and fill path to a .bmp file
void export_to_bmp(Polygon poly, vec_t* path, int path_count, const char* filename);

// Releases memory allocated in the polygon structure
void free_polygon(Polygon* poly);

#endif
