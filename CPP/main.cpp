#include <stdio.h>
#include "lodepng.h"
#include <vector>
#include <string>
#include <queue>
#include <map>
#include <math.h>

struct point {
    long int x,y,z;
    bool operator<(const point& p) const {
        if (x<p.x) return true;
        if (x>p.x) return false;
        if (y<p.y) return true;
        if (y>p.y) return false;
        if (z<p.z) return true;
        if (z>p.z) return false;
        return false;
    }
};

typedef std::map<point, size_t> id_map_t;

id_map_t id_map;
id_map_t point_map;

struct element {
    point p;
    bool interior;
    bool vtu_export;
    std::array<size_t, 27> con;
    std::array<size_t, 8> cell;
    size_t component;
};

std::vector<element> lattice;
std::vector<point> points;

const int pb_len=80;
int pb_now = 0;

void pb_tick(const size_t i,const size_t n) {
    int now = floor(1.0 *pb_len*i/n);
    if (now != pb_now) {
        printf("[");
        int k=0;
        for (; k<now; k++) printf("=");
        for (; k<pb_len; k++) printf(" ");
        printf("]\r");
        fflush(stdout);
        pb_now = now;
    }
    if (i == n) printf("\n");
}


int main() {

    // Selection
    unsigned int dx=40, dy=40, dz=40, nx=400, ny=400, nz=400; // Export only a cropped window
    bool exportInteriorOnly = true; // Export only the interior cells
    bool writeText = false; // Write text file with zero and ones
    unsigned int Lx=512, Ly=512, Lz=512; // Total size
    size_t comp_sel = 1; // Number of components (from first) to export. Set to 0 for ALL    
    
    unsigned int error;
    char filename[1024];
    size_t count = 0;
    FILE* f = NULL;
    
    printf("Generating interior:\n");
    if (writeText) f = fopen("frac1.txt","w");
    for (int z = 0; z<Lz; z++) {
    //int z=0; {
        //printf("File: %d\n",z);
        sprintf(filename,"cut_%04d.png",z+1);
        unsigned int width, height;
        std::vector<unsigned char> image;
        error = lodepng::decode(image, width, height, filename);
        if(error) {
            fprintf(stderr, "LODEPNG: %s: %s (error %u)\n", filename, lodepng_error_text(error), error);
            return -1;
        } else if ((width != Lx) || (height != Ly)) {
            fprintf(stderr, "Wrong dimensions of %s\n", filename);
            return -1;
        }
        for (int y=0; y<Ly; y++) {
            for (int x=0; x<Lx; x++) {
                if (image[1+4*(x+width*y)] < 128) {
                    if (writeText) fprintf(f,"0 ");
                    point p;
                    p.x=x;
                    p.y=y;
                    p.z=z;
                    element el;
                    el.interior = true;
                    el.component = 0;
                    el.p = p;
                    lattice.push_back(el);
                    count++;
                } else {
                    if (writeText) fprintf(f,"1 ");
                }                
            }
            if (writeText) fprintf(f,"\n");
        }
        pb_tick(z+1,Lz);
    }

    printf("Generating map:\n");
    for (size_t i=0; i<lattice.size(); i++) {
        point p = lattice[i].p;
        if (id_map.find(p) != id_map.end()) {
            fprintf(stderr, "Element in the map\n");
            return -1;
        }
        id_map[p] = i;
        pb_tick(i+1,lattice.size());
    }
    
    size_t component = 0;
    std::queue< size_t > Q;
    std::vector< size_t > comp_size;
    printf("Generating components:\n");
    size_t marked = 0;
    for (size_t j=0; j<lattice.size(); j++) if (lattice[j].component == 0) {
        component++;
        size_t n = 0;
        lattice[j].component = component;
        marked++; n++;
        Q.push(j);
        while (! Q.empty()) {
            size_t i = Q.front();
            point p = lattice[i].p;
            int k=0;
            for (int z=-1; z<=1; z++) {
                for (int y=-1; y<=1; y++) {
                    for (int x=-1; x<=1; x++) {
                        point np;
                        np.x = (p.x + x + Lx) % Lx;
                        np.y = (p.y + y + Ly) % Ly;
                        np.z = (p.z + z + Lz) % Lz;
                        id_map_t::iterator it = id_map.find(np);
                        if (it != id_map.end()) {
                            size_t k = it->second;
                            size_t ncom = lattice[k].component;
                            if (ncom == 0) {
                                lattice[k].component = component;
                                marked++; n++;
                                Q.push(k);                                
                            } else if (ncom != component) {
                                printf("Error. this should not happen.\n");
                                exit(-1);
                            }
                        }
                    }
                }
            }
            Q.pop();
            pb_tick(marked,lattice.size());
        }
        comp_size.push_back(n);
    }
    printf("Connected components: %ld\n", component);
    printf("Sizes of first <30 components:");
    for (int i=0; i<30; i++) if (i<comp_size.size()) printf(" %ld", comp_size[i]);
    printf("\n");

    if (comp_sel > 0) {
        printf("Selecting subset:\n");
        printf("Only components: %d-%d\n",1,comp_sel);
        std::vector<size_t> reindex;
        size_t shift = 0;
        for (size_t j=0; j<lattice.size(); j++) {
            if (lattice[j].component <= comp_sel) {
                if (shift > 0) {
                    lattice[j - shift] = lattice[j];
                }
            } else {
                shift++;
            }
        }
        printf("Max shift: %ld\n", shift);
        lattice.resize(lattice.size() - shift);
        printf("Regenerating map:\n");
        id_map.clear();
        for (size_t i=0; i<lattice.size(); i++) {
            point p = lattice[i].p;
            if (id_map.find(p) != id_map.end()) {
                fprintf(stderr, "Element in the map\n");
                return -1;
            }
            id_map[p] = i;
            pb_tick(i+1,lattice.size());
        }
    }
    
    if (writeText) fclose(f);
    printf("Generating connections:\n");
    for (size_t i=0; i<lattice.size(); i++) {
        point p = lattice[i].p;
        int k=0;
        for (int z=-1; z<=1; z++) {
            for (int y=-1; y<=1; y++) {
                for (int x=-1; x<=1; x++) {
                    point np;
                    size_t id;
                    np.x = (p.x + x + Lx) % Lx;
                    np.y = (p.y + y + Ly) % Ly;
                    np.z = (p.z + z + Lz) % Lz;
                    id_map_t::iterator it = id_map.find(np);
                    if (it != id_map.end()) {
                        id = it->second;
                    } else if (lattice[i].interior) {
                        element el;
                        el.p = np;
                        el.interior = false;
                        id = lattice.size();
                        id_map[el.p] = id;
                        lattice.push_back(el);
                    } else {
                        id = i;
                    }
                    lattice[i].con[k] = id;
                    k++;
                }
            }
        }      
        lattice[i].vtu_export =
            (dx <= p.x) && (dx+nx > p.x) &&
            (dy <= p.y) && (dy+ny > p.y) &&
            (dz <= p.z) && (dz+nz > p.z) &&
            (lattice[i].interior || (!exportInteriorOnly));
        k=0;
        if (lattice[i].vtu_export) {
            for (int z=0; z<=1; z++) {
                for (int y=0; y<=1; y++) {
                    for (int x=0; x<=1; x++) {
                        point np;
                        size_t id;
                        np.x = p.x + x;
                        np.y = p.y + y;
                        np.z = p.z + z;
                        id_map_t::iterator it = point_map.find(np);
                        if (it != point_map.end()) {
                            id = it->second;
                        } else {
                            id = points.size();
                            point_map[np] = id;
                            points.push_back(np);
                        }
                        lattice[i].cell[k] = id;
                        k++;
                    }
                }
            }
        }
        pb_tick(i+1,lattice.size());
    }
    printf("Writing connectivity:\n");
    f = fopen("frac1.cxn","w");
    fprintf(f,"LATTICESIZE %lu\n",lattice.size());
    fprintf(f,"BASE_LATTICE_DIM %d %d %d\n",20,20,20); // this is a mockup
    fprintf(f,"d 3\n");
    fprintf(f,"Q 27\n");
    fprintf(f,"OFFSET_DIRECTIONS\n");
    for (int z=-1; z<=1; z++) {
        for (int y=-1; y<=1; y++) {
            for (int x=-1; x<=1; x++) {
                fprintf(f,"[%d,%d,%d]",x,y,z);
                if ((x!=1) || (y!=1) || (z!=1))
                    fprintf(f,",");
                else
                    fprintf(f,"\n");
            }
        }
    }
    fprintf(f,"NODES\n");
    for (size_t i=0; i<lattice.size(); i++) {
        fprintf(f,"%lu %lg %lg %lg", i, 0.5+lattice[i].p.x, 0.5+lattice[i].p.y, 0.5+lattice[i].p.z);
        for (int k=0; k<27; k++) fprintf(f," %lu", lattice[i].con[k]);
        char* label = "NA";
        if (lattice[i].interior) {
            label = "Interior";
        } else {
            label = "Wall";
        }
        if (lattice[i].vtu_export) {
            fprintf(f," 1 %s\n", label);
        } else {
            fprintf(f," 2 %s %s\n", label, "HIDE");
        }
        pb_tick(i+1,lattice.size());
    }
    fclose(f);
    printf("Writing points:\n");
    f = fopen("frac1.cell","w");
    fprintf(f,"N_POINTS %lu\n",points.size());
    size_t cells = 0;
    for (size_t i=0; i<lattice.size(); i++) if (lattice[i].vtu_export) cells++;
    
    fprintf(f,"N_CELLS %lu\n",cells);
    fprintf(f,"POINTS\n");
    for (size_t i=0; i<points.size(); i++) {
        fprintf(f,"%lg %lg %lg\n", (double) points[i].x, (double) points[i].y, (double) points[i].z);
        pb_tick(i+1,points.size());
    }
    fprintf(f,"CELLS\n");
    printf("Writing cells:\n");
    size_t k=0;
    for (size_t i=0; i<lattice.size(); i++) if (lattice[i].vtu_export) {
        fprintf(f,"%lu %lu %lu %lu %lu %lu %lu %lu\n",
            lattice[i].cell[0],
            lattice[i].cell[1],
            lattice[i].cell[3],
            lattice[i].cell[2],
            lattice[i].cell[4],
            lattice[i].cell[5],
            lattice[i].cell[7],
            lattice[i].cell[6]);
        pb_tick(k+1,cells);
        k++;
    }
    fclose(f);
    return 0;
}