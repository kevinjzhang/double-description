#ifndef __DOUBLEDESCRIPTION_ALT_H
#define __DOUBLEDESCRIPTION_ALT_H

#include <regina-core.h>
#include <enumerate/enumconstraints.h>
#include <enumerate/ordering.h>
#include <maths/ray.h>
#include <maths/matrix.h>
#include <iterator>
#include <vector>
#include <maths/rational.h>

using namespace std;

namespace regina {

class Ray;
class ProgressTracker;

class REGINA_API DoubleDescriptionAlt {
    public:
        template <typename RayClass>
        static void enumerateExtremalRaysAlt(const MatrixInt& subspace, 
            const EnumConstraints* constraints);

    private:
        struct RayAlt {
            unsigned long timeAlive;
            vector<uint32_t> zeroSet; 
            vector<Rational> innerProductVector;

            RayAlt (int unitIndex, const MatrixInt& subspace, vector<unsigned long>& ordering);

            RayAlt (vector<uint32_t>& zeroSet, vector<Rational>& innerProductVector);

            int sign();

            template <typename RayClass>
            void recover(RayClass* dest, const MatrixInt& subspace);
        };

        DoubleDescriptionAlt();

        static bool isCompatible(RayAlt* ray1, RayAlt* ray2);

        static bool isAdjacent(vector<RayAlt*>& src, RayAlt* ray1, RayAlt* ray2);

        static RayAlt* constructRay(RayAlt* ray1, RayAlt* ray2, int hyperPlane, const MatrixInt& subspace);

        static bool intersectHyperplaneAlt(int currentHyperplane, vector<RayAlt*>& src, 
            vector<RayAlt*>& dest, const MatrixInt& subspace);
};

}
#include "doubledescription-alt-impl.h"
#endif