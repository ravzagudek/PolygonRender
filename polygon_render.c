#include "polygon_render.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

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

vec_t intersect(vec_t A1, vec_t B1, vec_t A2, vec_t B2) {
    double x1 = A1.x; double y1 = A1.y; double x2 = B1.x; double y2 = B1.y;
    double x3 = A2.x; double y3 = A2.y; double x4 = B2.x; double y4 = B2.y;

    double den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

    if (fabs(den) < DBL_EPSILON) {
        return (vec_t){NAN, NAN, NAN};
    }

    double px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / den;
    double py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / den;

    vec_t p = {px, py, 0.0f};
    if (
        p.x >= fmin(A1.x, B1.x) && p.x <= fmax(A1.x, B1.x) &&
        p.x >= fmin(A2.x, B2.x) && p.x <= fmax(A2.x, B2.x) &&
        p.y >= fmin(A1.y, B1.y) && p.y <= fmax(A1.y, B1.y) && 
        p.y >= fmin(A2.y, B2.y) && p.y <= fmax(A2.y, B2.y) ){
        return p;
    } else {
        return (vec_t) {NAN, NAN, NAN};
    }
    
}

vec_t find_normal_unit_vector(vec_t p1, vec_t p2){
    double vec_dx = p2.x - p1.x;
    double vec_dy = p2.y - p1.y;
    double len = sqrt(vec_dx * vec_dx + vec_dy * vec_dy);
    double nx = 0, ny = 0;
    if (len > DBL_EPSILON) {
        nx = -vec_dy/len ; ny = vec_dx/len;
    }
    vec_t n = {nx, ny, 0.0f};
    return n;
}

Polygon offset_polygon(Polygon poly, InfillParams params) {
    Polygon newPoly = {NULL, 0};
    if (poly.points==NULL || poly.count<=0) return newPoly;
    

    double offset_dist = params.offset_distance;
    newPoly.points = malloc(sizeof(vec_t) * poly.count);
    if (newPoly.points == NULL) {
        newPoly.count = 0; 
        return newPoly;
    }
    newPoly.count = poly.count;

    for (int i = 0; i < poly.count; i++) {
        vec_t P1 = poly.points[(i - 1 + poly.count) % poly.count]; // previous corner
        vec_t P2 = poly.points[i]; // current corner
        vec_t P3 = poly.points[(i + 1) % poly.count]; // next corner

        vec_t unit_n1 = find_normal_unit_vector(P1, P2);
        vec_t unit_n2 = find_normal_unit_vector(P2, P3);
       
        vec_t A1 = {P1.x + offset_dist * unit_n1.x, P1.y + offset_dist * unit_n1.y, 0.0f};
        vec_t B1 = {P2.x + offset_dist * unit_n1.x, P2.y + offset_dist * unit_n1.y, 0.0f};
        vec_t A2 = {P2.x + offset_dist * unit_n2.x, P2.y + offset_dist * unit_n2.y, 0.0f};
        vec_t B2 = {P3.x + offset_dist * unit_n2.x, P3.y + offset_dist * unit_n2.y, 0.0f};

        vec_t ip = intersect(A1, B1, A2, B2);

        // Checks that the intersection point is valid 
        if (isnan(ip.x) || isnan(ip.y)) {
            newPoly.points[i] = P2; // uses the original point
        } else {
            newPoly.points[i].x = ip.x;
            newPoly.points[i].y = ip.y;
            newPoly.points[i].z = 0.0f;
        }
    }
    return newPoly;
}

vec_t* queue_to_array(xqueue_t* q, int* count){
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

vec_t* generate_horizontal_path(Polygon poly, InfillParams params, int* output_count) {
    if (output_count == NULL) return NULL;
    *output_count = 0;
    
    Polygon rotPoly = rotate_polygon(poly, params.infill_angle);
    if (rotPoly.count == 0) {
        free_polygon(&rotPoly); return NULL;
    }

    double min_y = rotPoly.points[0].y; double max_y = rotPoly.points[0].y;

    for (int i = 0; i < rotPoly.count; i++) {
        if (rotPoly.points[i].y < min_y) min_y = rotPoly.points[i].y;
        if (rotPoly.points[i].y > max_y) max_y = rotPoly.points[i].y;
    }

    double scan_spacing = params.scan_radius * (1.0 - params.overlap_ratio);

    xqueue_t infill_queue = xqueueNew(sizeof(vec_t));
    if (infill_queue.array == NULL) {
        free_polygon(&rotPoly); return NULL;
    }

    double epsilon = DBL_EPSILON * 1000.0;
    vec_t previous_line_end_point = {NAN, NAN, NAN};
    bool forward_scan = true;

    for (double y = min_y + epsilon; y <= max_y - epsilon; y += scan_spacing) {
        vec_t* current_line_intersections = (vec_t*)malloc(sizeof(vec_t) * rotPoly.count * 2);
        if (current_line_intersections == NULL) { 
            xqueueFree(&infill_queue);
            free_polygon(&rotPoly); return NULL;
        }

        int numIntersections = 0; 

        vec_t line_start = {-1e12, y, 0.0f};
        vec_t line_end = {1e12, y, 0.0f};

        for (int j = 0; j < rotPoly.count; j++) {
            vec_t p1 = rotPoly.points[j];
            vec_t p2 = rotPoly.points[(j + 1) % rotPoly.count];

            /* Special case: If the y-value of the intersection point is very close to the scan line
			and the segment is horizontal, we should not add this point twice and do special processing (the entire edge is on the line). */
            if (fabs(p1.y - y) < epsilon && fabs(p2.y - y) < epsilon) {
                if (p1.x < p2.x) {
                    bool is_duplicate1 = false;
                    for (int k = 0; k < numIntersections; k++) {
                        if (fabs(current_line_intersections[k].x - p1.x) < epsilon && fabs(current_line_intersections[k].y - p1.y) < epsilon) {
                            is_duplicate1 = true; break;
                        }
                    }
                    if (!is_duplicate1) current_line_intersections[numIntersections++] = p1;

                    bool is_duplicate2 = false;
                    for (int k = 0; k < numIntersections; k++) {
                        if (fabs(current_line_intersections[k].x - p2.x) < epsilon && fabs(current_line_intersections[k].y - p2.y) < epsilon) {
                            is_duplicate2 = true; break;
                        }
                    }
                    if (!is_duplicate2) current_line_intersections[numIntersections++] = p2;

                } else {
                    bool is_duplicate1 = false;
                    for (int k = 0; k < numIntersections; k++) {
                        if (fabs(current_line_intersections[k].x - p2.x) < epsilon && fabs(current_line_intersections[k].y - p2.y) < epsilon) {
                            is_duplicate1 = true; break;
                        }
                    }
                    if (!is_duplicate1) current_line_intersections[numIntersections++] = p2;

                    bool is_duplicate2 = false;
                    for (int k = 0; k < numIntersections; k++) {
                        if (fabs(current_line_intersections[k].x - p1.x) < epsilon && fabs(current_line_intersections[k].y - p1.y) < epsilon) {
                            is_duplicate2 = true; break;
                        }
                    }
                    if (!is_duplicate2) current_line_intersections[numIntersections++] = p1;
                }
                continue;
            }

            vec_t ip = intersect(line_start, line_end, p1, p2);

            if (!isnan(ip.x) && !isnan(ip.y) &&
                (ip.x >= fmin(p1.x, p2.x) - epsilon && ip.x <= fmax(p1.x, p2.x) + epsilon) &&
                (ip.y >= fmin(p1.y, p2.y) - epsilon && ip.y <= fmax(p1.y, p2.y) + epsilon)) {

                bool is_duplicate = false;
                for (int k = 0; k < numIntersections; k++) {
                    if (fabs(current_line_intersections[k].x - ip.x) < epsilon && fabs(current_line_intersections[k].y - ip.y) < epsilon) {
                        is_duplicate = true;
                        break;
                    }
                }

                if (!is_duplicate) {
                    current_line_intersections[numIntersections++] = ip;
                }
            }
        }

        for (int i = 0; i < numIntersections - 1; i++) {
            for (int j = i + 1; j < numIntersections; j++) {
                if (current_line_intersections[i].x > current_line_intersections[j].x) {
                    vec_t temp = current_line_intersections[j];
                    current_line_intersections[j] = current_line_intersections[i];
                    current_line_intersections[i] = temp;
                }
            }
        }

        if (numIntersections >= 2 && numIntersections % 2 == 0) {
            if (!isnan(previous_line_end_point.x)) {
                if (forward_scan) {
                    xqueueEnqueue(&infill_queue, &previous_line_end_point);
                    xqueueEnqueue(&infill_queue, &current_line_intersections[0]);
                } else {
                    xqueueEnqueue(&infill_queue, &previous_line_end_point);
                    xqueueEnqueue(&infill_queue, &current_line_intersections[numIntersections - 1]);
                }
            }

            if (forward_scan) {
                for (int i = 0; i < numIntersections; i += 2) {
                    xqueueEnqueue(&infill_queue, &current_line_intersections[i]);
                    xqueueEnqueue(&infill_queue, &current_line_intersections[i + 1]);
                }
                previous_line_end_point = current_line_intersections[numIntersections - 1];
            } else {
                for (int i = numIntersections - 1; i >= 0; i -= 2) {
                    xqueueEnqueue(&infill_queue, &current_line_intersections[i]);
                    xqueueEnqueue(&infill_queue, &current_line_intersections[i - 1]);
                }
                previous_line_end_point = current_line_intersections[0];
            }
        }
        free(current_line_intersections);
        forward_scan = !forward_scan;
    }

    vec_t* infill_path = queue_to_array(&infill_queue, output_count);
    Polygon resultPoly = {infill_path, *output_count};
    Polygon rotated = rotate_polygon(resultPoly, -params.infill_angle);
    
    free(infill_path);
    xqueueFree(&infill_queue);
    free_polygon(&rotPoly);
    
    vec_t* result = rotated.points;
    return result;
}

vec_t* generate_spiral_path(Polygon poly, InfillParams params, int* output_count){
    xqueue_t infill_queue = xqueueNew(sizeof(vec_t)); 
    if (infill_queue.array == NULL) {
        *output_count = 0;
        return NULL;
    }

    Polygon currPoly = { NULL, poly.count };
    currPoly.points = malloc(sizeof(vec_t) * poly.count);
    if (!currPoly.points) { xqueueFree(&infill_queue); return NULL; }
    memcpy(currPoly.points, poly.points, sizeof(vec_t)*poly.count);
    
    double scan_spacing = params.scan_radius * (1.0 - params.overlap_ratio);
    double epsilon = DBL_EPSILON * 1000.0;
    if (!(scan_spacing > epsilon)) {
        xqueueFree(&infill_queue);
        return NULL;
    }

    params.offset_distance = -(scan_spacing);
    double lastArea = calc_poly_area(currPoly);
    while (true) {
        Polygon nextPoly = offset_polygon(currPoly, params);
        // When the polygon area shrinks sufficiently or begins to grow, end the loop
        double currArea = calc_poly_area(nextPoly);
        if (fabs(currArea) < DBL_EPSILON * 10000.0 || fabs(currArea)>=fabs(lastArea)){
            // Find the least polygon's center, free polygons and break the loop
            double sumx=0, sumy=0, cx=0, cy=0;
            for (int i = 0; i < currPoly.count; i++){
                sumx += currPoly.points[i].x; sumy += currPoly.points[i].y;
            }
            cx = sumx / currPoly.count; cy = sumy / currPoly.count;
            vec_t center = {cx, cy, 0.0f}; xqueueEnqueue(&infill_queue, &center);
            free_polygon(&currPoly); free_polygon(&nextPoly); 
            break;
        }
        // Save all the points of the rings in the same direction
        if (calc_poly_area(currPoly)*calc_poly_area(nextPoly)<0) reverse_poly_points(&nextPoly);
        for (int i = 0; i < nextPoly.count; i++){
            xqueueEnqueue(&infill_queue, &nextPoly.points[i]);
        }
        free_polygon(&currPoly); 
        currPoly = nextPoly; 
        lastArea=currArea;
    }

    *output_count = xqueueSize(&infill_queue);
    vec_t* tempArr = queue_to_array(&infill_queue, output_count);
    xqueueFree(&infill_queue);
    xqueue_t reversedInfillQueue = xqueueNew(sizeof(vec_t));
    for (int i = *output_count-1; i >= 0; i--) {
        xqueueEnqueue(&reversedInfillQueue, &tempArr[i]);
    }

    vec_t* result = queue_to_array(&reversedInfillQueue, output_count);
    xqueueFree(&reversedInfillQueue);
    free(tempArr);
    return result;
}

vec_t* generate_infill_path(Polygon poly, InfillParams params, int* output_count, int choice) {
    if (output_count == NULL) return NULL;
    *output_count = 0;

    switch (choice){
    case 1:
        return generate_horizontal_path(poly, params, output_count);
    case 2:
        return generate_spiral_path(poly, params, output_count);
    default:
        return NULL;
    }
}

void set_pixel(unsigned char* image, int width, int height, int x, int y, unsigned char r, unsigned char g, unsigned char b) {
    if (x >= 0 && x < width && y >= 0 && y < height) {
        int index = (y * width + x) * 3;
        image[index + 0] = b;
        image[index + 1] = g;
        image[index + 2] = r;
    }
}

void draw_line(unsigned char* image, int width, int height,
                vec_t p1_geo, vec_t p2_geo,
                double scale_x, double scale_y,
                double min_x, double min_y,
                unsigned char r, unsigned char g, unsigned char b) {

    int x1 = (int)((p1_geo.x - min_x) * scale_x);
    int y1 = (int)((p1_geo.y - min_y) * scale_y);
    int x2 = (int)((p2_geo.x - min_x) * scale_x);
    int y2 = (int)((p2_geo.y - min_y) * scale_y);

    x1 = fmax(0, fmin(width - 1, x1));
    y1 = fmax(0, fmin(height - 1, y1));
    x2 = fmax(0, fmin(width - 1, x2));
    y2 = fmax(0, fmin(height - 1, y2));

    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int sx = (x1 < x2) ? 1 : -1;
    int sy = (y1 < y2) ? 1 : -1;
    int err = dx - dy;

    while (1) {
        set_pixel(image, width, height, x1, y1, r, g, b);
		if (x1 == x2 && y1 == y2) break;
        int e2 = 2 * err;
        if (e2 > -dy) { err -= dy; x1 += sx; }
        if (e2 < dx) { err += dx; y1 += sy; }
    }
}

void export_to_bmp(Polygon poly, vec_t* path, int path_count, const char* filename) {
    double min_x = DBL_MAX, max_x = -DBL_MAX;
    double min_y = DBL_MAX, max_y = -DBL_MAX;

    if (poly.count > 0) {
        min_x = poly.points[0].x; max_x = poly.points[0].x;
        min_y = poly.points[0].y; max_y = poly.points[0].y;
        for (int i = 0; i < poly.count; i++) {
            if (poly.points[i].x < min_x) min_x = poly.points[i].x;
            if (poly.points[i].x > max_x) max_x = poly.points[i].x;
            if (poly.points[i].y < min_y) min_y = poly.points[i].y;
            if (poly.points[i].y > max_y) max_y = poly.points[i].y;
        }
    }

    if (path != NULL) {
        for (int i = 0; i < path_count; i++) {
            if (path[i].x < min_x) min_x = path[i].x;
            if (path[i].x > max_x) max_x = path[i].x;
            if (path[i].y < min_y) min_y = path[i].y;
            if (path[i].y > max_y) max_y = path[i].y;
        }
    }

    if (fabs(max_x - min_x) < DBL_EPSILON) {
        min_x -= 50.0; max_x += 50.0;
    }
    if (fabs(max_y - min_y) < DBL_EPSILON) {
        min_y -= 50.0; max_y += 50.0;
    }

    double current_range_x = max_x - min_x;
    double current_range_y = max_y - min_y;

    double padding_factor = 0.2;
    double padding_x = fmax(current_range_x * padding_factor, 10.0);
    double padding_y = fmax(current_range_y * padding_factor, 10.0);

    min_x -= padding_x;
    max_x += padding_x;
    min_y -= padding_y;
    max_y += padding_y;

    double range_x = max_x - min_x;
    double range_y = max_y - min_y;

    int width = 1200;
    int height = 1200;

    if (range_x < DBL_EPSILON) range_x = 1.0;
    if (range_y < DBL_EPSILON) range_y = 1.0;
    double scale = fmin((double)(width-1) / range_x, (double)(height-1) / range_y);
    double scale_x = scale;
    double scale_y = scale;

    unsigned char* image = (unsigned char*)calloc(width * height * 3, 1);
    if (image == NULL) return;

    for (int i = 0; i < width * height * 3; ++i) {
        image[i] = 255;
    }

    for (int i = 0; i < poly.count; i++) { // Draws the polygon frame (black)
        vec_t p1 = poly.points[i];
        vec_t p2 = poly.points[(i + 1) % poly.count];
        draw_line(image, width, height, p1, p2, scale_x, scale_y, min_x, min_y, 0, 0, 0);
    }

    for (int i = 0; i < path_count - 1; i++) { // Draws the infill path (red)
        draw_line(image, width, height, path[i], path[i+1], scale_x, scale_y, min_x, min_y, 255, 0, 0);
    }

    FILE* fp = fopen(filename, "wb");
    if (fp == NULL) {
        free(image);
        return;
    }

    unsigned char bmp_file_header[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
    unsigned int file_size = width * height * 3 + 54;
    bmp_file_header[2] = (unsigned char)(file_size);
    bmp_file_header[3] = (unsigned char)(file_size >> 8);
    bmp_file_header[4] = (unsigned char)(file_size >> 16);
    bmp_file_header[5] = (unsigned char)(file_size >> 24);
    fwrite(bmp_file_header, 1, 14, fp);

    unsigned char bmp_info_header[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
    bmp_info_header[4] = (unsigned char)(width);
    bmp_info_header[5] = (unsigned char)(width >> 8);
    bmp_info_header[6] = (unsigned char)(width >> 16);
    bmp_info_header[7] = (unsigned char)(width >> 24);
    bmp_info_header[8] = (unsigned char)(height);
    bmp_info_header[9] = (unsigned char)(height >> 8);
    bmp_info_header[10] = (unsigned char)(height >> 16);
    bmp_info_header[11] = (unsigned char)(height >> 24);
    fwrite(bmp_info_header, 1, 40, fp);

    fwrite(image, 1, width * height * 3, fp);

    fclose(fp);
    free(image);
}

void free_polygon(Polygon* poly) {
    if (poly != NULL && poly->points != NULL) {
        free(poly->points);
        poly->points = NULL;
        poly->count = 0;
    }
}