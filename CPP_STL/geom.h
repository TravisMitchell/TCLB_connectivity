#include <queue>
#include <map>
#include <vector>
#include <list>

#ifndef GEOM_H

typedef int coord_t;
typedef float coord_c;
typedef size_t idx_t;
const idx_t NULL_ID = (size_t) -1;

struct Coords {
    coord_t x,y,z;
    bool operator<(const Coords& other) const {
        if (x < other.x) return true;
        if (x > other.x) return false;
        if (y < other.y) return true;
        if (y > other.y) return false;
        if (z < other.z) return true;
        return false;
    }
};

struct Coords_c {
    coord_c x, y, z;
    bool operator<(const Coords_c& other) const {
        if(x < other.x) return true;
        if(x > other.x) return false;
        if(y < other.y) return true;
        if(y > other.y) return false;
        if(z < other.z) return true;
        return false;
    }
};

struct Ray {
    Coords_c p, q;
};

struct Node {
    Coords coords;
    idx_t id;
    idx_t con[27];
    int nodeType;
    bool finish;
    bool inter = false;
    bool body = false;
};


struct CellNode {
    Coords_c centre;
    idx_t con[8];
};

struct STL_Tri {
    float norm[3];
    float p1[3];
    float p2[3];
    float p3[3];
    float cent[3];
    short int v;
};

static float dotProduct(float v1[3], float v2[3]) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

static void crossProduct(float v1[3], float v2[3], float * ret) {
    ret[0] = v1[1] * v2[2] - v1[2] * v2[1];
    ret[1] = v1[2] * v2[0] - v1[0] * v2[2];
    ret[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

typedef std::queue<idx_t> queue_t;
typedef std::map<Coords, idx_t> map_t;
typedef std::vector<Node> vector_t;
typedef std::map<Coords_c, idx_t> map_c;
typedef std::vector<CellNode> vector_c;
typedef std::vector<Coords_c> vector_p;
typedef std::list<STL_Tri*> tri_list;

#define GEOM_H 1
#endif