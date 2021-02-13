#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "Preprocessor.h"
#include "string.h"
#include <map>
#include <queue>

int Preprocessor::generateDomain() {
    // queue to store generated nodes to visit in BFS-type fashion
    queue_t Q;
    // map to store generated nodes/coordinates
    map_t M;
    // map to store generated cell points
    map_c M_c;

    Node n;
    CellNode cn;

    // starting point for BFS
    n.coords.x = x0;    
    n.coords.y = y0;    
    n.coords.z = z0;

    cn.centre.x = x0;
    cn.centre.y = y0;
    cn.centre.z = z0;

    // if you want to add custom ranges to define nodetypes, do that here
    float* xinBounds = {};
    float* xoutBounds = {};
    float* yinBounds = {};
    float* youtBounds = {};

    for(int i = 0; i < 8; i++)
        cn.con[i] = 0;

    n.id = nodes.size();
    n.nodeType = 0;
    nodes.push_back(n);
    generateCellData(&cn, 1, M_c);
    cells.push_back(cn);
    M[n.coords] = n.id;
    // initialise BFS
    Q.push(n.id);
    int nCount = 0;

    while (! Q.empty()) {
        idx_t id = Q.front();
        // expand on all 27 connectivity directions
        for (int i=0; i<27; i++) {
            Coords c;
            c.x = (nodes[id].coords.x + d3q27_x[i] + nx) % nx;
            c.y = (nodes[id].coords.y + d3q27_y[i] + ny) % ny;
            c.z = (nodes[id].coords.z + d3q27_z[i] + nz) % nz;
            map_t::iterator it = M.find(c);
            idx_t new_id = id;
            if (it == M.end()) {
                if ((! nodes[id].finish)) {
                    nCount++;
                    // count number of generated nodes so far
                    if(nCount % 10000 == 0)
                        printf("Generated %d nodes\n", nCount);
                    Node m;
                    CellNode cm;
                    new_id = m.id = nodes.size();
                    m.coords = c;
                    cm.centre.x = c.x;
                    cm.centre.y = c.y;
                    cm.centre.z = c.z;

                    generateCellData(&cm, 1, M_c);
                    cells.push_back(cm);
                    m.finish = true;
                    Coords_c tc = latticeToSTLTransform(c);
                    int nt = getNodeType(tc, xinBounds, xoutBounds, youtBounds, 0,0,0);
                    
                    if(nt == 0) {
                        m.nodeType = 0;
                    }

                    nodes.push_back(m);
                    M[m.coords] = m.id;
                    Q.push(m.id);
                }
            } else {
                new_id = it->second;
            }
            nodes[id].con[i] = -1;
            nodes[id].con[i] = new_id;
            if (new_id != id) {
                if(id != new_id && !nodes[id].finish) {
                    // this is specific to the granular pack case - a hardcoded channel around the domain
                    bool inChannel = false;
                    /*if(nodes[new_id].coords.y > 10 && nodes[new_id].coords.y < 224) {
                        if(nodes[new_id].coords.x == 0 || nodes[new_id].coords.x == 223 || nodes[new_id].coords.z == 0 || nodes[new_id].coords.z == 213) {
                            nodes[new_id].finish = true;
                            inChannel = true;
                        }
                    }*/
                    if(!inChannel) {
                        // check for collision against the STL
                        if(!collideSTL(nodes[id].coords, nodes[new_id].coords)) {
                            nodes[new_id].finish = false;
                        } else {
                            nodes[new_id].body = true;
                        }
                    }
                }
            }
        }
        Q.pop();
    }
    // print domain size
    printf("size: %ld\n", nodes.size());
    if (nodes.size() < 40) {
        for (size_t i=0; i<nodes.size(); i++) {
            printf("%ld (%3d,%3d,%3d)", i, nodes[i].coords.x, nodes[i].coords.y, nodes[i].coords.z);
            for (int j=0; j<27; j++) {
                printf(" %ld", nodes[i].con[j]);
            }
            printf("\n");
        }
    }
}

int Preprocessor::getNodeType(Coords_c pos, float* xinBounds, float* xoutBounds, float* youtBounds, int nXInlets, int nXOutlets, int nYOutlets) {
    // if you want multiple node types, e.g. inlet, outlets, you can do something like the following:
    // if x_inlet, return 1
    /*for(int i = 0; i < nXInlets; i++) {
        if(pos.x >= xinBounds[i*6] && pos.y >= xinBounds[i*6 + 1] && pos.z >= xinBounds[i*6 + 2] && pos.x < xinBounds[i*6 + 3] && pos.y < xinBounds[i*6 + 4] && pos.z < xinBounds[i*6 + 5])
            return 1;
    }

    // if x_outlet, return 2
    for(int i = 0; i < nXOutlets; i++) {
        if(pos.x >= xoutBounds[i*6] && pos.y >= xoutBounds[i*6 + 1] && pos.z >= xoutBounds[i*6 + 2] && pos.x < xoutBounds[i*6 + 3] && pos.y < xoutBounds[i*6 + 4] && pos.z < xoutBounds[i*6 + 5])
            return 2;
    }

    // if y_outlet, return 3
    for(int i = 0; i < nYOutlets; i++) {
        if(pos.x >= youtBounds[i*6] && pos.y >= youtBounds[i*6 + 1] && pos.z >= youtBounds[i*6 + 2] && pos.x < youtBounds[i*6 + 3] && pos.y < youtBounds[i*6 + 4] && pos.z < youtBounds[i*6 + 5])
            return 3;
    }*/

    // else return 0 - MRT
    return 0;


}

/**
 * Function to generate cell data for a single node. Iterates over 8 corners of a cell and tries to add a point.
 **/
void Preprocessor::generateCellData(CellNode * node, float dx, map_c & M_c) {
    int cell_dxs[24] = {-1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1};
    for(int i = 0; i < 8; i++) {
        Coords_c pt;
        pt.x = node->centre.x + (dx/2) * cell_dxs[i*3];
        pt.y = node->centre.y + (dx/2) * cell_dxs[i*3 + 1];
        pt.z = node->centre.z + (dx/2) * cell_dxs[i*3 + 2];
        map_c::iterator it = M_c.find(pt);
        idx_t pt_id = NULL_ID;
        // if we haven't already got this point, create it
        if(it == M_c.end()) {
            pt_id = cellPoints.size();
            cellPoints.push_back(pt);
            // add it to our map
            M_c[pt] = pt_id;
        } else {
            pt_id = it->second;
        }
        // set the connectivity info in the passed node to include this
        node->con[i] = pt_id;
    }
}

/**
 * Function to transform lattice coords (int indices to float coords)
 **/
Coords_c Preprocessor::latticeToSTLTransform(Coords c) {
    Coords_c tc;
    // temporary for this one case - shift the thing up 10 pixels
    int nx_ = 214;
    int ny_ = 214;
    int dy_ = 10;

    tc.x = xmin + (c.x * ((xmax - xmin) / nx_));
    tc.y = ymin + ((c.y - dy_) * ((ymax - ymin) / ny_));
    tc.z = zmin + (c.z * ((zmax - zmin) / nz));
    return tc;
}

bool Preprocessor::in_box(Coords & c) {
    if (c.x <=  0) return false;
    if (c.x >= nx) return false;
    if (c.y <=  0) return false;
    if (c.y >= ny) return false;
    if (c.z <=  0) return false;
    if (c.z >= nz) return false;
    return true;
}

bool Preprocessor::triRayCollision(Coords_c & c1, Coords_c & c2, STL_Tri * tri) {
    float kEps = 0.0000001;
    // get points in float array form so it works with our functions
    float p[3] = {c1.x, c1.y, c1.z};
    float q[3] = {c2.x, c2.y, c2.z};
    
    float h[3], q_[3];
    float a, f, u, v;
    float e1[3] = {tri->p2[0] - tri->p1[0], tri->p2[1] - tri->p1[1], tri->p2[2] - tri->p1[2]};
    float e2[3] = {tri->p3[0] - tri->p1[0], tri->p3[1] - tri->p1[1], tri->p3[2] - tri->p1[2]};
    float d[3] = {q[0] - p[0], q[1] - p[1], q[2] - p[2]};
    crossProduct(d, e2, h);
    a = dotProduct(e1, h);
    if(a > -kEps && a < kEps) {
        return false;
    }
        
    f = 1/a;
    float s[3] = {p[0] - tri->p1[0], p[1] - tri->p1[1], p[2] - tri->p1[2]};
    u = f * dotProduct(s, h);
    if(u < 0.0 || u > 1.0) {
        return false;
    }
        
    crossProduct(s, e1, q_);
    v = f * dotProduct(d, q_);
    if(v < 0.0 || u + v > 1.0) {
        return false;
    }
    
    float t = f * dotProduct(e2, q_);
    if(t > kEps && t <= 1.0) {
        return true;
    }

    return false;

}

bool Preprocessor::collideSTL(Coords & c1, Coords & c2) {
    //if(!in_box(c2))
    //    return true;

    Coords_c tc1 = latticeToSTLTransform(c1);
    Coords_c tc2 = latticeToSTLTransform(c2);

    Ray r;
    r.p = tc1;
    r.q = tc2;
    std::list<int> nodeIDs = triTree.findNodes(&triTree.root, r);
    avgNoNodes += nodeIDs.size();
    collideSTLCalls += 1;
    if(nodeIDs.size() < minNoNodes)
        minNoNodes = nodeIDs.size();
    //if(nodeIDs.size() != 35)
        //printf("node list size: %d\n", nodeIDs.size());
    std::list<int>::iterator it;
    for(it = nodeIDs.begin(); it != nodeIDs.end(); it++) {
        TreeNode * vol = triTree.nodes[(*it)];
        tri_list::iterator tri_it;
        for(tri_it = vol->tris.begin(); tri_it != vol->tris.end(); tri_it++) {
            if(triRayCollision(tc1, tc2, (*tri_it))) {
                return true;
            }
        }
    }

    return false;
}

/**
 * Function to load triangular data into a linear datastructure.
 * Currently assumes plaintext stl file.
 **/
int Preprocessor::loadSTL(const char * filename) {
    FILE* file = fopen(filename, "r");
    int ret;
    
    char junk[80];
    ntri = 0;
    if(file == NULL) {
        printf("ERROR: can't open file\n");
        return -1;
    }
    // read header
    ret = fscanf(file, "%s %s (%d triangles", junk, junk, &ntri, junk);
    fgets(junk, 80, file);
    printf("Reading %d triangles\n", ntri);
    tris = (STL_Tri*)malloc(ntri * sizeof(STL_Tri));
    printf("A\n");
    // read triangle data
    for(int i = 0; i < ntri; i++) {
        // read normal
        ret = fscanf(file, "%s %s %f %f %f\n", junk, junk, &tris[i].norm[0], &tris[i].norm[1], &tris[i].norm[2]);
        fgets(junk, 80, file);
        // read vertices
        fscanf(file, "%s %f %f %f\n", junk, &tris[i].p1[0], &tris[i].p1[1], &tris[i].p1[2]);
        fscanf(file, "%s %f %f %f\n", junk, &tris[i].p2[0], &tris[i].p2[1], &tris[i].p2[2]);
        fscanf(file, "%s %f %f %f\n", junk, &tris[i].p3[0], &tris[i].p3[1], &tris[i].p3[2]);
        
        // calculate centre
        tris[i].cent[0] = (tris[i].p1[0] + tris[i].p2[0] + tris[i].p3[0]) / 3;
        tris[i].cent[1] = (tris[i].p1[1] + tris[i].p2[1] + tris[i].p3[1]) / 3;
        tris[i].cent[2] = (tris[i].p1[2] + tris[i].p2[2] + tris[i].p3[2]) / 3;
        
        // initialise bounds
        if(i == 0) {
            xmin = tris[i].p1[0];
            xmax = tris[i].p1[0];
            ymin = tris[i].p1[1];
            ymax = tris[i].p1[1];
            zmin = tris[i].p1[2];
            zmax = tris[i].p1[2];
        }
        // update min, max bounds
        updateBounds(tris[i].p1);
        updateBounds(tris[i].p2);
        updateBounds(tris[i].p3);
        
        // read other stuff at the end
        fgets(junk, 80, file);
        fgets(junk, 80, file);
    }
    printf("STL Bounds: %e, %e : %e, %e : %e, %e\n", xmin, xmax, ymin, ymax, zmin, zmax);
    return 0;
}

void Preprocessor::updateBounds(float p[3]) {
    xmin = std::min(xmin, p[0]);
    xmax = std::max(xmax, p[0]);
    ymin = std::min(ymin, p[1]);
    ymax = std::max(ymax, p[1]);
    zmin = std::min(zmin, p[2]);
    zmax = std::max(zmax, p[2]);
}

void Preprocessor::generateKDTree() {
    // initialise the root node
    TreeNode root;
    root.leaf = false;
    root.level = 0;
    root.id = triTree.currentID++;
    root.xmin = xmin;
    root.xmax = xmax;
    root.ymin = ymin;
    root.ymax = ymax;
    root.zmin = zmin;
    root.zmax = zmax;
    root.dir = 0;

    // load all triangles into root node
    for(int i = 0; i < ntri; i++)
        root.tris.push_back(&tris[i]);

    // call the function to recursively subdivide the tree
    triTree.root = root;
    triTree.nodes.push_back(&root);
    triTree.subdivide(&triTree.root, 18);
}

/**
 * Function to output connectivity information in correct format to a file.
 * Assumes nodes are correctly generated.
 **/
void Preprocessor::writeToConnectivity(char* filename) {
    // create a filestream and open the file
    std::ofstream confile;
    std::ofstream cellfile;
    char* confn = (char*) malloc(strlen(filename) + 4 + 1);
    strcpy(confn, filename);
    strcat(confn, ".con");

    char* cellfn = (char*) malloc(strlen(filename) + 5 + 1);
    strcpy(cellfn, filename);
    strcat(cellfn, ".cell");

    confile.open(confn);

    // write header information
    confile << "LATTICESIZE " << nodes.size() << "\n";
    confile << "BASE_LATTICE_DIM " << nx << " " << ny << " " << nz << "\n";
    confile << "d " << d_ << "\n";
    confile << "Q " << Q_ << "\n";
    confile << "OFFSET_DIRECTIONS\n";
    
    // write the offset directions to be read by TCLB
    for(int q = 0; q < Q_; q++) {
        confile << "[" << d3q27_x[q] << "," << d3q27_y[q] << "," << d3q27_z[q] << "]";
        if(q < Q_ - 1)
            confile << ",";
        else
            confile << "\n";
    }

    confile << "NODES\n";
    int intCount = 0;
    // write nodal coordinates and connectivity information
    for(int i = 0; i < nodes.size(); i++) {
        confile << i << " " << nodes[i].coords.x << " " << nodes[i].coords.y << " " << nodes[i].coords.z << " ";
        // write connectivity
        for(int q = 0; q < Q_; q++) {
            if(nodes[i].con[q] == (size_t)-1) {
                confile << -1 << " ";
            } else {
                confile << nodes[i].con[q] << " ";
            }
            
        }
        // write number of group labels
        confile << "1 ";
        // write group label
        if(nodes[i].coords.x == 0 && nodes[i].coords.y == -1 && nodes[i].coords.z == 20) {
            printf("Test node is: %d\n", nodes[i].finish);
        }
        if(nodes[i].finish)
            if(nodes[i].body)
                confile << "Body\n";
            else
                confile << "Wall\n";
        else {
            if(nodes[i].nodeType == 0)
                confile << "Collide\n";
            else
                printf("Error - unknown nodetype\n");
        }
            
        
        if(nodes[i].inter && in_box(nodes[i].coords))
            intCount++;
    }
    printf("Intersect nodes count: %d\n", intCount);

    confile.close();

    // write cell data
    cellfile.open(cellfn);
    cellfile << "N_POINTS " << cellPoints.size() << "\n";
    cellfile << "N_CELLS " << cells.size() << "\n";
    cellfile << "POINTS\n";
    
    for(int p = 0; p < cellPoints.size(); p++) {
        cellfile << cellPoints[p].x << " " << cellPoints[p].y << " " << cellPoints[p].z << "\n";
    }

    cellfile << "CELLS\n";

    for(int c = 0; c < cells.size(); c++) {
        for(int d = 0; d < 8; d++) {
            cellfile << cells[c].con[d];
            if(d < 7)
                cellfile << " ";
            else
                cellfile << "\n";
        }
    }
    
    cellfile.close();
}