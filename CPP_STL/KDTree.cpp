#include "KDTree.h"

/**
 * Recursive function to find the leaf bounding box of k-d tree containing the search point.
 **/
std::list<int> KDTree::findNodes(TreeNode * current, Ray ray) {
    std::list<int> containingNodes;

    if(current->leaf) {
        containingNodes.push_back(current->id);
        return containingNodes;
    }

    // if ray intersects left box, traverse that
    if(boxRayIntersect(current->left, ray)) {
        std::list<int> leftNodes = findNodes(current->left, ray);
        std::list<int>::iterator it;
        // add all that come out to the master list
        for(it = leftNodes.begin(); it != leftNodes.end(); it++)
            containingNodes.push_back(*it);
    }
    // if ray intersects right box, traverse that too
    if(boxRayIntersect(current->right, ray)) {
        std::list<int> rightNodes = findNodes(current->right, ray);
        std::list<int>::iterator it;
        // add all that come out to the master list
        for(it = rightNodes.begin(); it != rightNodes.end(); it++)
            containingNodes.push_back(*it);
    }

    return containingNodes;
}

bool KDTree::boxRayIntersect(TreeNode * node, Ray ray) {
    // need to fix to be actually good
    /*if(ray.p.x >= node->xmin && ray.p.x < node->xmax && ray.p.y >= node->ymin && ray.p.y < node->ymax && ray.p.z >= node->zmin && ray.p.z < node->zmax)
        return true;
    if(ray.q.x >= node->xmin && ray.q.x < node->xmax && ray.q.y >= node->ymin && ray.q.y < node->ymax && ray.q.z >= node->zmin && ray.q.z < node->zmax)
        return true;
    return false;*/

    float dx = ray.q.x - ray.p.x;
    float dy = ray.q.y - ray.p.y;
    float dz = ray.q.z - ray.p.z;

    float txmin = (node->xmin - ray.p.x) / dx;
    float txmax = (node->xmax - ray.p.x) / dx;

    if(txmin > txmax) std::swap(txmin, txmax);

    float tymin = (node->ymin - ray.p.y) / dy;
    float tymax = (node->ymax - ray.p.y) / dy;

    if(tymin > tymax) std::swap(tymin, tymax);

    if((txmin > tymax) || (tymin > txmax))
        return false;
    
    if(tymin > txmin)
        txmin = tymin;
    
    if(tymax < txmax)
        txmax = tymax;
    
    float tzmin = (node->zmin - ray.p.z) / dz;
    float tzmax = (node->zmax - ray.p.z) / dz;

    if(tzmin > tzmax) std::swap(tzmin, tzmax);

    if((txmin > tzmax) || (tzmin > txmax))
        return false;
    
    if(tzmin > txmin)
        txmin = tzmin;
    
    if(tzmax < txmax)
        txmax = tzmax;
    
    return true;

}

double KDTree::avgCellSize() {
    return countCellSize(&root) / countCells(&root);
}

double KDTree::countCellSize(TreeNode * current) {
    if(current == nullptr)
        return 0;
    if(current->leaf) {
        return current->zmax - current->zmin;
    }
    
    return countCellSize(current->left) + countCellSize(current->right);
}

int KDTree::countCells(TreeNode * current) {
    if(current == nullptr)
        return 0;
    if(current->leaf) {
        return 1;
    }
    
    return countCells(current->left) + countCells(current->right);
}

/**
 * Function recursively called to partition the k-d tree and propagate traingles down into each node.
 * 
 **/
void KDTree::subdivide(TreeNode * current, int levels) {
    // once we get to the specified depth, finish
    if(levels == 0) {
        current->leaf = true;
        return;
    }
    
    // once a node and its volume only contains a few triangles, don't keep dividing
    if(current->tris.size() < 10) {
        current->leaf = true;
        return;
    }
    // initialise left and right nodes
    TreeNode* left = new TreeNode();
    TreeNode* right = new TreeNode();
    
    // otherwise, do the partition
    if(current->dir == 0) {
        current->splitPt = (current->xmin + current->xmax) / 2;
        
        left->xmin = current->xmin;
        left->xmax = current->splitPt;
        left->ymin = current->ymin;
        left->ymax = current->ymax;
        left->zmin = current->zmin;
        left->zmax = current->zmax;
        
        right->xmin = current->splitPt;
        right->xmax = current->xmax;
        right->ymin = current->ymin;
        right->ymax = current->ymax;
        right->zmin = current->zmin;
        right->zmax = current->zmax;

        left->dir = 1;
        right->dir = 1;

        // propagate triangles down to children
        tri_list::iterator it;
        for(it = current->tris.begin(); it != current->tris.end(); ++it) {
            if((*it)->cent[0] <= current->splitPt || (*it)->p1[0] <= current->splitPt || (*it)->p2[0] <= current->splitPt || (*it)->p3[0] <= current->splitPt) {
                left->tris.push_back((*it));
            }
            if((*it)->cent[0] > current->splitPt || (*it)->p1[0] > current->splitPt || (*it)->p2[0] > current->splitPt || (*it)->p3[0] > current->splitPt) {
                right->tris.push_back((*it));
            }
        }
        
    } else if(current->dir == 1) {
        current->splitPt = (current->ymin + current->ymax) / 2;

        left->xmin = current->xmin;
        left->xmax = current->xmax;
        left->ymin = current->ymin;
        left->ymax = current->splitPt;
        left->zmin = current->zmin;
        left->zmax = current->zmax;
        
        right->xmin = current->xmin;
        right->xmax = current->xmax;
        right->ymin = current->splitPt;
        right->ymax = current->ymax;
        right->zmin = current->zmin;
        right->zmax = current->zmax;

        left->dir = 2;
        right->dir = 2;

        // propagate triangles down to children
        tri_list::iterator it;
        for(it = current->tris.begin(); it != current->tris.end(); ++it) {
            if((*it)->cent[1] <= current->splitPt || (*it)->p1[1] <= current->splitPt || (*it)->p2[1] <= current->splitPt || (*it)->p3[1] <= current->splitPt) {
                left->tris.push_back((*it));
            } 
            if((*it)->cent[1] > current->splitPt || (*it)->p1[1] > current->splitPt || (*it)->p2[1] > current->splitPt || (*it)->p3[1] > current->splitPt) {
                right->tris.push_back((*it));
            }
        }
    } else if(current->dir == 2) {
        current->splitPt = (current->zmin + current->zmax) / 2;

        left->xmin = current->xmin;
        left->xmax = current->xmax;
        left->ymin = current->ymin;
        left->ymax = current->ymax;
        left->zmin = current->zmin;
        left->zmax = current->splitPt;
        
        right->xmin = current->xmin;
        right->xmax = current->xmax;
        right->ymin = current->ymin;
        right->ymax = current->ymax;
        right->zmin = current->splitPt;
        right->zmax = current->zmax;
        
        left->dir = 0;
        right->dir = 0;

        // propagate triangles down to children
        tri_list::iterator it;
        for(it = current->tris.begin(); it != current->tris.end(); ++it) {
            if((*it)->cent[2] <= current->splitPt || (*it)->p1[2] <= current->splitPt || (*it)->p2[2] <= current->splitPt || (*it)->p3[2] <= current->splitPt) {
                left->tris.push_back((*it));
            }
            if((*it)->cent[2] > current->splitPt || (*it)->p1[2] > current->splitPt || (*it)->p2[2] > current->splitPt || (*it)->p3[2] > current->splitPt) {
                right->tris.push_back((*it));
            }
        }
    }

    left->level = current->level + 1;
    right->level = current->level + 1;
    left->id = currentID++;
    right->id = currentID++;
    nodes.push_back(left);
    nodes.push_back(right);
    current->left = left;
    current->right = right;
    /*printf("X\n");
    printf("Xmin: %e\n", current->left->xmin);
    printf("Y\n");
    printf("Left size: %d, right size: %d\n", current->left->tris.size(), current->right->tris.size());*/
    subdivide(current->left, levels - 1);
    subdivide(current->right, levels - 1);
}