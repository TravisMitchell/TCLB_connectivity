#ifndef PREPROC_H

#include "geom.h"
#include "KDTree.h"

class Preprocessor {
public:
    int nx, ny, nz, x0, y0, z0;
    float xmin, xmax, ymin, ymax, zmin, zmax;
    int d_;
    int Q_;

    vector_t nodes;
    vector_c cells;
    vector_p cellPoints;
    coord_t d3q27_x[27] = {0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1,1};
    coord_t d3q27_y[27] = {0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1};
    coord_t d3q27_z[27] = {0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

    int ntri;
    STL_Tri * tris;
    KDTree triTree;
    int outCount = 0;
    long avgNoNodes = 0;
    long collideSTLCalls = 0;
    int minNoNodes = 1000;

    //void init(int xmin_, int xmax_, int ymin_, int ymax_, int zmin_, int zmax_);
    int loadSTL(const char* filename);
    void updateBounds(float p[3]);
    int generateDomain();
    void generateKDTree();
    int getNodeType(Coords_c pos, float* xinBounds, float* xoutBounds, float* youtBounds, int nXInlets, int nXOutlets, int nYOutlets);
    Coords_c latticeToSTLTransform(Coords c);
    bool in_box(Coords & c);
    bool triRayCollision(Coords_c & c1, Coords_c & c2, STL_Tri * tri);
    bool collideSTL(Coords & c1, Coords & c2);
    void writeToConnectivity(char* filename);
    void generateCellData(CellNode * node, float dx, map_c & M_c);
};

#define PREPROC_H 1
#endif