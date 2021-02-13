/**
 * Class for the bounding-volume k-d tree which stores domain triangles in a hierarchical format
 **/
#ifndef KDTREE_H
#include "geom.h"

struct TreeNode {
    int id;
    float xmin, xmax, ymin, ymax, zmin, zmax;
    float splitPt;
    int dir; // 0 = x
    int level; // 0 = root
    bool leaf;
    tri_list tris;

    TreeNode * left;
    TreeNode * right;

    TreeNode() {
        
    };
};

class KDTree {
    public:
    TreeNode root;
    std::vector<TreeNode*> nodes;
    int currentID = 0;
    double avgCellSize();
    double countCellSize(TreeNode * current);
    int countCells(TreeNode * current);
    std::list<int> findNodes(TreeNode * current, Ray ray);
    bool boxRayIntersect(TreeNode * node, Ray ray);
    void subdivide(TreeNode * current, int levels);
};

#define KDTREE_H 1
#endif