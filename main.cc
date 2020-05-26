#include <vector>
#include <chrono> 

#include "doubledescription-alt.h"
#include<triangulation/dim3.h>
#include<triangulation/example3.h>
#include<surfaces/normalsurfaces.h>
#include<surfaces/normalcoords.h>
#include<surfaces/nsvectorquad.h>

#include <maths/ray.h>
#include "outputIterator.h"
//Extra
#include<enumerate/doubledescription.h>
#include<enumerate/enumconstraints.h>

//#define DISPLAYANS

using namespace regina;

int main() {
    std::vector<Triangulation<3>*> triangulations;
    
    triangulations.push_back(Example<3>::s2xs1());
    triangulations.push_back(Example<3>::rp2xs1());
    triangulations.push_back(Example<3>::poincareHomologySphere());
    triangulations.push_back(Example<3>::weeks());
    //triangulations.push_back(Example<3>::weberSeifert());
    int i = 0;
    for(auto tri : triangulations) {
        MatrixInt* subspace = makeMatchingEquations(tri, NS_QUAD);
        //vector of sets 
        EnumConstraints* enumConstraints = makeEmbeddedConstraints(tri, NS_QUAD);
        auto it = NSVectorQuadOutputIterator();

        cout << "Method 1: " << i << endl;
        auto start = chrono::high_resolution_clock::now();
        DoubleDescription::enumerateExtremalRays<NSVectorQuad>(it, *subspace, enumConstraints);
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        cout << duration.count() << endl;

        start = chrono::high_resolution_clock::now();
        cout << "Method 2: " << i << endl;
        DoubleDescriptionAlt::enumerateExtremalRaysAlt<NSVectorQuad>(*subspace, enumConstraints);
        stop = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        cout << duration.count() << endl;
        i++;
    }
    return 0;
}



