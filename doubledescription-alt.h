#ifndef __DOUBLEDESCRIPTION_ALT_H
#define __DOUBLEDESCRIPTION_ALT_H

#include <bitset>
#include <regina-core.h>
#include <enumerate/enumconstraints.h>
#include <enumerate/ordering.h>
#include <maths/ray.h>
#include <maths/matrix.h>
#include <iterator>
#include <vector>
#include <maths/rational.h>

using namespace std;
#define setbits                     std::bitset<128>

namespace regina {

class Ray;
class ProgressTracker;

class REGINA_API DoubleDescriptionAlt {
    public:
        enum Algorithm {USE_SIMPLE, USE_TRIE, USE_GRAPH, USE_MATRIX, USE_COMBINED, TEST_ALL}; 

        struct RunOptions {
            Algorithm algorithm;
        };

        template <typename RayClass>
        static void enumerateExtremalRaysAlt(const MatrixInt& subspace, RunOptions options);
    private:
        class RayAlt : public Ray {
            public:
                unsigned long timeAlive;
                setbits zeroSet; 
                vector<RayAlt*> neighbours;

                RayAlt (int unitIndex, const MatrixInt& subspace, vector<unsigned long>& ordering);
                RayAlt (RayAlt* ray1, RayAlt* ray2, int hyperPlane, const MatrixInt& subspace);

                template <typename RayClass>
                void recover(RayClass* dest, const MatrixInt& subspace);
        };

        DoubleDescriptionAlt();

        static bool isCompatible(RayAlt* ray1, RayAlt* ray2);

        static bool isAdjacent(const MatrixInt& subspace, SetTrie& setTrie, vector<RayAlt*>& src, RayAlt* ray1, RayAlt* ray2, 
            int currentHyperplane, vector<unsigned long>& hyperplaneOrdering, RunOptions options);

        static bool isAdjacentStandard(vector<RayAlt*>& src, RayAlt* ray1, RayAlt* ray2);

        static bool isAdjacentTrie(const MatrixInt& subspace, SetTrie& setTrie, RayAlt* ray1, RayAlt* ray2);

        static bool isAdjacentGraph(vector<RayAlt*>& src, RayAlt* ray1, RayAlt* ray2);

        static bool isAdjacentAlgebraic(const MatrixInt& constraints, RayAlt* ray1, RayAlt* ray2, int currentHyperplane, 
            vector<unsigned long>& hyperplaneOrdering);

        static bool intersectHyperplaneAlt(int currentHyperplane, vector<RayAlt*>& src, 
            vector<RayAlt*>& dest, const MatrixInt& subspace, vector<unsigned long>& hyperplaneOrdering, RunOptions options);
};

}
#include "doubledescription-alt-impl.h"
#endif