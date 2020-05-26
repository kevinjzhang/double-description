#ifndef __NS_VECTOR_QUAD_OUTPUT_H
#define __NS_VECTOR_QUAD_OUTPUT_H

#include<surfaces/nsvectorquad.h>

using namespace regina;
class NSVectorQuadOutputIterator {

    public:
        int index = 0;
        vector<NSVectorQuad*> ans;

        NSVectorQuadOutputIterator operator=(NSVectorQuad* vec) {
        #ifdef DISPLAYANS
            Ray ray = vec->coords();
            for (int i = 0; i < ray.size(); i++) {
                std::cout << ray[i] << " ";
            }
            std::cout << endl;
            if(index == ans.size()) {
                ans.push_back(vec);
            } else {
                ans[index] = vec;
            }
        #endif
            return *this;
        }

        NSVectorQuadOutputIterator operator++(int) {
        #ifdef DISPLAYANS
            if(index < ans.size()) {
                index++;
            }
        #endif
            return *this;
        }

        NSVectorQuadOutputIterator operator*() {
            return *this;
        }

};

#endif