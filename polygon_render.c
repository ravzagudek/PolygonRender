#include "polygon_render.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

// Auxiliary Functions

vec_t intersect(vec_t A1, vec_t B1, vec_t A2, vec_t B2) {
    double x1 = A1.x; double y1 = A1.y; double x2 = B1.x; double y2 = B1.y;
    double x3 = A2.x; double y3 = A2.y; double x4 = B2.x; double y4 = B2.y;

    double den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

    if (fabs(den) < DBL_EPSILON) return (vec_t){NAN, NAN, NAN};

    double px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / den;
    double py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / den;

    vec_t p = {px, py, 0.0f};
    if (p.x >= fmin(A1.x, B1.x) - DBL_EPSILON && p.x <= fmax(A1.x, B1.x) + DBL_EPSILON &&
    p.x >= fmin(A2.x, B2.x) - DBL_EPSILON && p.x <= fmax(A2.x, B2.x) + DBL_EPSILON &&
    p.y >= fmin(A1.y, B1.y) - DBL_EPSILON && p.y <= fmax(A1.y, B1.y) + DBL_EPSILON &&
    p.y >= fmin(A2.y, B2.y) - DBL_EPSILON && p.y <= fmax(A2.y, B2.y) + DBL_EPSILON) {
        return p;
    } else return (vec_t) {NAN, NAN, NAN};
    
}

vec_t find_normal_unit_vector(vec_t p1, vec_t p2){
    vec_t vec = {p2.x - p1.x, p2.y - p1.y, 0.0f};
    double len = sqrt(vec.x * vec.x + vec.y * vec.y);
    if (len < 1e-12) return (vec_t){0.0, 0.0, 0.0};
    return (vec_t) {vec.y/len, -vec.x/len, 0.0f};
}

vec_t normalize(vec_t p){
    double len = sqrt(p.x*p.x + p.y*p.y);
    if (len < 1e-12) return (vec_t){0.0, 0.0, 0.0};
    return (vec_t) {p.x/len, p.y/len, 0.0f};
}

double calc_poly_area(Polygon poly){
    double area = 0.0;
    for (int i=0; i<poly.count; i++){
        vec_t p1 = poly.points[i];
        vec_t p2 = poly.points[(i+1)%poly.count];
        area += (p1.x * p2.y) - (p2.x * p1.y);
    }
    return area*0.5;
}

void reverse_poly_points(Polygon* poly){
    for (int i = 0; i < poly->count/2; i++){
        vec_t temp = poly->points[i];
        poly->points[i] = poly->points[poly->count - 1 - i];
        poly->points[poly->count - 1 - i] = temp;
    }
}

vec_t* queue_to_array(xqueue_t* q, int* count) {
    if (q == NULL || count == NULL) return NULL;

    *count = xqueueSize(q);

    vec_t* pointArray = (vec_t*)malloc(*count * sizeof(vec_t));
    if (pointArray == NULL)  return NULL;

    xqueue_t temp = xqueueNew(sizeof(vec_t)); // Temporary queue
    if(temp.array == NULL) {
        free(pointArray);
        return NULL;
    }
    for (int i = 0; i < *count; i++){
        vec_t currPoint;
        xqueueDequeue(q, &currPoint);
        pointArray[i] = currPoint;
        xqueueEnqueue(&temp, &currPoint);
    }
    for (int i = 0; i < *count; i++){ // Fill the original queue again
        vec_t tempPoint;
        xqueueDequeue(&temp, &tempPoint);
        xqueueEnqueue(q, &tempPoint);
    }
    xqueueFree(&temp);

    return pointArray;
}

void free_polygon(Polygon* poly) {
    if (poly != NULL && poly->points != NULL) {
        free(poly->points);
        poly->points = NULL;
        poly->count = 0;
    }
}

int cmp(const void *a, const void *b){
    double dx = ((vec_t*)a)->x - ((vec_t*)b)->x;
    if (dx<0) return -1;
    if (dx>0) return 1;
    return 0;
}

// Main Functions

Polygon offset_polygon(Polygon poly, InfillParams params) {
    Polygon newPoly = {NULL, 0};
    if (poly.points==NULL || poly.count<=0) return newPoly;
    
    newPoly.points = malloc(sizeof(vec_t) * poly.count);
    if (newPoly.points == NULL) {
        newPoly.count = 0; 
        return newPoly;
    }
    newPoly.count = poly.count;
    
    double area = calc_poly_area(poly);
    double direction = (area>=0) ? -1.0 : 1.0;
    double d = -params.offsetDist * direction;
    for (int i = 0; i < poly.count; i++) {
        vec_t P1 = poly.points[(i - 1 + poly.count) % poly.count]; // previous corner
        vec_t P2 = poly.points[i]; // current corner
        vec_t P3 = poly.points[(i + 1) % poly.count]; // next corner

        vec_t unit_n1 = find_normal_unit_vector(P1, P2);
        vec_t unit_n2 = find_normal_unit_vector(P2, P3);
       
        vec_t A1 = {P1.x + d * unit_n1.x, P1.y + d * unit_n1.y, 0.0f};
        vec_t B1 = {P2.x + d * unit_n1.x, P2.y + d * unit_n1.y, 0.0f};
        vec_t A2 = {P2.x + d * unit_n2.x, P2.y + d * unit_n2.y, 0.0f};
        vec_t B2 = {P3.x + d * unit_n2.x, P3.y + d * unit_n2.y, 0.0f};

        vec_t ip = intersect(A1, B1, A2, B2);

        if (isnan(ip.x) || isnan(ip.y)) newPoly.points[i] = P2;
        else newPoly.points[i] = ip;
    }

    return newPoly;
}

Polygon rotate_polygon(Polygon poly, double angle) {
    Polygon newPoly = {NULL, 0};
    if (poly.points == NULL || poly.count <= 0) return newPoly;
    
    newPoly.points = malloc(sizeof(vec_t) * poly.count);
	if (newPoly.points == NULL) {
        newPoly.count = 0;
        return newPoly;
    }
    newPoly.count = poly.count;

    double sum_x = 0, sum_y = 0;
    double rad = angle * M_PI / 180.0;

    for (int i = 0; i < poly.count; i++) {
        sum_x += poly.points[i].x;
        sum_y += poly.points[i].y;
    }
    double cx = sum_x / poly.count;
    double cy = sum_y / poly.count;

    for (int i = 0; i < poly.count; i++) {
        double x = poly.points[i].x;
        double y = poly.points[i].y;
        newPoly.points[i].x = cos(rad) * (x - cx) - sin(rad) * (y - cy) + cx;
        newPoly.points[i].y = sin(rad) * (x - cx) + cos(rad) * (y - cy) + cy;
        newPoly.points[i].z = 0.0f; // Z coordinate remains 0 for 2D operations
    }

    return newPoly;
}

vec_t* generate_horizontal_path(Polygon poly, InfillParams params, int* outputCount) {
    if (outputCount == NULL) return NULL;
    *outputCount = 0;
    
    Polygon rotPoly = rotate_polygon(poly, params.infillAngle);
    if (rotPoly.count == 0) {
        free_polygon(&rotPoly); return NULL;
    }
    Polygon offsetPoly = offset_polygon(rotPoly, params);
    free_polygon(&rotPoly);
    rotPoly = offsetPoly;

    double minY = rotPoly.points[0].y; double maxY = rotPoly.points[0].y;
    for (int i = 0; i < rotPoly.count; i++) {
        if (rotPoly.points[i].y < minY) minY = rotPoly.points[i].y;
        if (rotPoly.points[i].y > maxY) maxY = rotPoly.points[i].y;
    }

    double scanSpacing = params.scanRadius * (1.0 - params.overlapRatio);

    xqueue_t infillQueue = xqueueNew(sizeof(vec_t));
    if (infillQueue.array == NULL) {
        free_polygon(&rotPoly); return NULL;
    }

    double epsilon = DBL_EPSILON * 1000.0;
    vec_t prevLineEnd = {NAN, NAN, NAN};
    bool forwardScan = true;

    for (double y = minY + epsilon; y <= maxY - epsilon; y += scanSpacing) {
        vec_t* currLineIntersects = (vec_t*)malloc(sizeof(vec_t) * rotPoly.count * 2);
        if (currLineIntersects == NULL) { 
            xqueueFree(&infillQueue);
            free_polygon(&rotPoly); return NULL;
        }
        int numIntersections = 0; 

        vec_t lineStart = {-1e9, y, 0.0f};
        vec_t lineEnd = {1e9, y, 0.0f};

        for (int j = 0; j < rotPoly.count; j++) {
            vec_t p1 = rotPoly.points[j];
            vec_t p2 = rotPoly.points[(j + 1) % rotPoly.count];

            /* Special case: If the y-value of the intersection point is very close to the scan line
			and the segment is horizontal, we should not add this point twice and do special processing (the entire edge is on the line). */
            if (fabs(p1.y - y) < epsilon && fabs(p2.y - y) < epsilon) {
                if (p1.x < p2.x) {
                    bool is_duplicate1 = false;
                    for (int k = 0; k < numIntersections; k++) {
                        if (fabs(currLineIntersects[k].x - p1.x) < epsilon && fabs(currLineIntersects[k].y - p1.y) < epsilon) {
                            is_duplicate1 = true; break;
                        }
                    }
                    if (!is_duplicate1) currLineIntersects[numIntersections++] = p1;

                    bool is_duplicate2 = false;
                    for (int k = 0; k < numIntersections; k++) {
                        if (fabs(currLineIntersects[k].x - p2.x) < epsilon && fabs(currLineIntersects[k].y - p2.y) < epsilon) {
                            is_duplicate2 = true; break;
                        }
                    }
                    if (!is_duplicate2) currLineIntersects[numIntersections++] = p2;

                } else {
                    bool is_duplicate1 = false;
                    for (int k = 0; k < numIntersections; k++) {
                        if (fabs(currLineIntersects[k].x - p2.x) < epsilon && fabs(currLineIntersects[k].y - p2.y) < epsilon) {
                            is_duplicate1 = true; break;
                        }
                    }
                    if (!is_duplicate1) currLineIntersects[numIntersections++] = p2;

                    bool is_duplicate2 = false;
                    for (int k = 0; k < numIntersections; k++) {
                        if (fabs(currLineIntersects[k].x - p1.x) < epsilon && fabs(currLineIntersects[k].y - p1.y) < epsilon) {
                            is_duplicate2 = true; break;
                        }
                    }
                    if (!is_duplicate2) currLineIntersects[numIntersections++] = p1;
                }
                continue;
            }

            vec_t ip = intersect(lineStart, lineEnd, p1, p2);
            if (!isnan(ip.x) && !isnan(ip.y) &&
                (ip.x >= fmin(p1.x, p2.x) - epsilon && ip.x <= fmax(p1.x, p2.x) + epsilon) &&
                (ip.y >= fmin(p1.y, p2.y) - epsilon && ip.y <= fmax(p1.y, p2.y) + epsilon)) {

                bool is_duplicate = false;
                for (int k = 0; k < numIntersections; k++) {
                    if (fabs(currLineIntersects[k].x - ip.x) < epsilon && fabs(currLineIntersects[k].y - ip.y) < epsilon) {
                        is_duplicate = true;
                        break;
                    }
                }

                if (!is_duplicate) {
                    currLineIntersects[numIntersections++] = ip;
                }
            }
        }

        qsort(currLineIntersects, numIntersections, sizeof(vec_t), cmp);

        if (numIntersections >= 2 && numIntersections % 2 == 0) {
            if (!isnan(prevLineEnd.x)) {
                if (forwardScan) {
                    xqueueEnqueue(&infillQueue, &prevLineEnd);
                    xqueueEnqueue(&infillQueue, &currLineIntersects[0]);
                } else {
                    xqueueEnqueue(&infillQueue, &prevLineEnd);
                    xqueueEnqueue(&infillQueue, &currLineIntersects[numIntersections - 1]);
                }
            }

            if (forwardScan) {
                for (int i = 0; i < numIntersections; i += 2) {
                    xqueueEnqueue(&infillQueue, &currLineIntersects[i]);
                    xqueueEnqueue(&infillQueue, &currLineIntersects[i + 1]);
                }
                prevLineEnd = currLineIntersects[numIntersections - 1];
            } else {
                for (int i = numIntersections - 1; i >= 0; i -= 2) {
                    xqueueEnqueue(&infillQueue, &currLineIntersects[i]);
                    xqueueEnqueue(&infillQueue, &currLineIntersects[i - 1]);
                }
                prevLineEnd = currLineIntersects[0];
            }
        }
        free(currLineIntersects);
        forwardScan = !forwardScan;
    }

    vec_t* infillPath = queue_to_array(&infillQueue, outputCount);
    Polygon resultPoly = {infillPath, *outputCount};
    Polygon rotated = rotate_polygon(resultPoly, -params.infillAngle);
    
    free(infillPath);
    xqueueFree(&infillQueue);
    free_polygon(&rotPoly);
    
    vec_t* result = rotated.points;
    return result;
}

vec_t* generate_spiral_path(Polygon poly, InfillParams params, int* outputCount){
    xqueue_t infillQueue = xqueueNew(sizeof(vec_t)); 
    if (infillQueue.array == NULL) {
        *outputCount = 0;
        return NULL;
    }

    double scanSpacing = -(params.scanRadius*(1-params.overlapRatio));

    Polygon currPoly = { NULL, poly.count };
    currPoly.points = malloc(sizeof(vec_t) * poly.count);
    if (!currPoly.points) { xqueueFree(&infillQueue); return NULL; }
    memcpy(currPoly.points, poly.points, sizeof(vec_t)*poly.count);
    double lastArea = calc_poly_area(currPoly);
    int i=1;
    while (true) {
        params.offsetDist = scanSpacing*i;
        Polygon nextPoly = offset_polygon(currPoly, params);
        // When the polygon area shrinks sufficiently or begins to grow, end the loop
        double currArea = calc_poly_area(nextPoly);
        if (fabs(currArea) < DBL_EPSILON * 100000.0 || fabs(currArea)>=fabs(lastArea) || nextPoly.count<3){
            double sumx=0, sumy=0, cx=0, cy=0;
            for (int i = 0; i < currPoly.count; i++){
                sumx += currPoly.points[i].x; sumy += currPoly.points[i].y;
            }
            vec_t center = {sumx / currPoly.count, sumy / currPoly.count, 0.0f}; 
            xqueueEnqueue(&infillQueue, &center);
            free_polygon(&currPoly); free_polygon(&nextPoly); 
            break;
        }
        // Save all the points of the rings in the same direction
        if (calc_poly_area(currPoly)*calc_poly_area(nextPoly)<0) reverse_poly_points(&nextPoly);
        for (int i = 0; i < nextPoly.count; i++){
            xqueueEnqueue(&infillQueue, &nextPoly.points[i]);
        }
        free_polygon(&nextPoly); lastArea=currArea;
        i++;
    }

    *outputCount = xqueueSize(&infillQueue);
    vec_t* tempArr = queue_to_array(&infillQueue, outputCount);
    xqueueFree(&infillQueue);
    xqueue_t reversedInfillQueue = xqueueNew(sizeof(vec_t));
    for (int i = *outputCount-1; i >= 0; i--) {
        xqueueEnqueue(&reversedInfillQueue, &tempArr[i]);
    }

    vec_t* result = queue_to_array(&reversedInfillQueue, outputCount);
    xqueueFree(&reversedInfillQueue);
    free(tempArr);
    return result;
}

vec_t* generate_infill_path(Polygon poly, InfillParams params, int* outputCount, int choice) {
    if (outputCount == NULL) return NULL;
    *outputCount = 0;

    switch (choice){
    case 1:
        return generate_horizontal_path(poly, params, outputCount);
    case 2:
        return generate_spiral_path(poly, params, outputCount);
    default:
        return NULL;
    }
}
