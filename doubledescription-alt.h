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
#define setbits                     bitset<128>

namespace regina {

class Ray;
class ProgressTracker;

class REGINA_API DoubleDescriptionAlt {
    public:
        template <typename RayClass>
        static void enumerateExtremalRaysAlt(const MatrixInt& subspace, 
            const EnumConstraints* constraints);

    private:
        class RayAlt : public Ray {
            public:
                unsigned long timeAlive;
                setbits zeroSet; 

                RayAlt (int unitIndex, const MatrixInt& subspace, vector<unsigned long>& ordering);
                RayAlt (RayAlt* ray1, RayAlt* ray2, int hyperPlane, const MatrixInt& subspace);

                template <typename RayClass>
                void recover(RayClass* dest, const MatrixInt& subspace);
        };

        DoubleDescriptionAlt();

        static vector<Ray> reduce(const MatrixInt& subspace);

        static bool isCompatible(RayAlt* ray1, RayAlt* ray2);

        static bool isAdjacent(vector<RayAlt*>& src, RayAlt* ray1, RayAlt* ray2);

        static bool isAdjacentAlgebraic(vector<Ray>& constraints, RayAlt* ray1, RayAlt* ray2);

        static bool intersectHyperplaneAlt(int currentHyperplane, vector<RayAlt*>& src, 
            vector<RayAlt*>& dest, const MatrixInt& subspace);
};

}
#include "doubledescription-alt-impl.h"
#endif