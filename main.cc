#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <type_traits>
#include <chrono> 

#include "setTrie.h"
#include "doubledescription-alt.h"
#include<triangulation/dim3.h>
#include<triangulation/example3.h>
#include<triangulation/detail/triangulation.h>
#include<surfaces/normalsurfaces.h>
#include<surfaces/normalcoords.h>
#include<surfaces/nsvectorquad.h>
#include<link/link.h>

#include <maths/ray.h>
#include "outputIterator.h"
//Extra
#include<enumerate/doubledescription.h>
#include<enumerate/enumconstraints.h>


using namespace regina;


int main(int argc, char *argv[]) {
    std::string inFile = std::string(argv[1]);
    std::string outFile = std::string(argv[2]);
    std::ifstream in (inFile, std::ifstream::in);
    std::ofstream out;
    out.open(outFile);
    int number;
    in >> number;
    for(int i = 0; i < number; i++) {
        std::string name;
        in >> name;
        int method; //Method >= 7 tells the program that a link signature is used instead
        in >> method;
        Triangulation<3>* triangulation;
        DoubleDescriptionAlt::RunOptions options;
        if (method <= 6) {
            options.algorithm = (DoubleDescriptionAlt::Algorithm) (method - 1);
            if (name == "weberSeifert") {
                triangulation = Example<3>::weberSeifert();
            } else if (name == "weeks") {
                triangulation = Example<3>::weeks();
            } else {
                triangulation = Triangulation<3>::fromIsoSig(name);
            }
        } else {
            options.algorithm = (DoubleDescriptionAlt::Algorithm) (method - 8);
            Link* link = Link::fromKnotSig(name);
            triangulation = link->complement();
            if (triangulation->isIdeal()) {
                triangulation->idealToFinite();
                if (!triangulation->intelligentSimplify()) {
                    triangulation->simplifyExhaustive();
                }
            }
        }
        MatrixInt* subspace = makeMatchingEquations(triangulation, NS_QUAD);
        EnumConstraints* enumConstraints = makeEmbeddedConstraints(triangulation, NS_QUAD);
        auto it = NSVectorQuadOutputIterator();
        //Execute algorithms
        auto start = chrono::high_resolution_clock::now();
        if (method == 0 || method == 7) { 
            DoubleDescription::enumerateExtremalRays<NSVectorQuad>(it, *subspace, enumConstraints);
        } else {
            DoubleDescriptionAlt::enumerateExtremalRaysAlt<NSVectorQuad>(*subspace, options);
        }
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        out << name << ": " << duration.count() << " Method: " << method << std::endl;
    }
    out.close();
    return 0;
}



