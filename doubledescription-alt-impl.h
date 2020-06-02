#ifndef __DOUBLEDESCRIPTION_ALT_IMPL_H
#define __DOUBLEDESCRIPTION_ALT_IMPL_H

//Bit size is 32 rounded down to multiple of 3
//#define LOG_FILENAME                "log.txt"
//#define MAKEGRAPH                   1
#define DISPLAYANS
#define BIT_SIZE                    63
#define FIRST_BIT                   0x1249249249249249 //0|001001....
#define SECOND_BIT                  0x2492492492492492 //0|010010....
#define THIRD_BIT                   0x4924924924924924 //0|100100....

#define setbits                     uint64_t

// #define TIMING                      1
// #define DEBUG                       1


#include <algorithm>
#include <iterator>
#include <vector>
#include <list>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <omp.h>
#include <bitset>
#include <regina-core.h>
#include <regina-config.h>
#include <maths/integer.h>

using namespace std;

namespace regina {

void printVec(vector<Rational> innerProductVector) {
    for(auto val : innerProductVector) {
        cout << val << " ";
    }
    cout << endl;    
}

#ifdef MAKEGRAPH
//Static graph in optimised version
//Visited must have space preallocated
void getConnectedComponent(int node, unordered_map<int, vector<int>>& graph, 
    vector<bool>& visited, vector<int>& nodes) {
    visited[node] = true;
    nodes.push_back(node);
    for(int adj : graph[node]) {
        if(!visited[adj]) {
            getConnectedComponent(adj, graph, visited, nodes);
        }
    }
}
#endif

struct LexicographicalOrder {      
    const MatrixInt* subspace;
    LexicographicalOrder (const MatrixInt& subspace) {
        this->subspace = &subspace;
    }

    bool operator() (int r1, int r2) {
        unsigned long cols = subspace->columns();
        for (int i = 0; i < cols; i++) {
            bool comp1 = subspace->entry(r1, i).isZero();
            bool comp2 = subspace->entry(r2, i).isZero();
            if (comp1 && !comp2) {
                return true;
            } else if (comp2 && !comp1) {
                return false;
            } else if (subspace->entry(r1, i) != subspace->entry(r2, i)) {
                return subspace->entry(r1, i) < subspace->entry(r2, i);
            }
        }
        //Return true if equal
        return true;
    }
};

LargeInteger gcd(LargeInteger a, LargeInteger b) {
    if(a > b) {
        auto temp = a;
        a = b;
        b = temp;
    }
    while(a != 0) {
        auto rem = b % a;
        b = a;
        a = rem;
    }
    return b;
}

//Constructor based on unit vector
DoubleDescriptionAlt::RayAlt::RayAlt (int unitIndex, const MatrixInt& subspace, vector<unsigned long>& ordering) {
    unsigned long eqns = subspace.rows();
    unsigned long dim = subspace.columns();
    zeroSet = UINT64_MAX;
    innerProductVector = vector<Rational>(eqns, 0);
    zeroSet ^= (1 << unitIndex);
    for(int i = 0; i < eqns; i++) {
        innerProductVector[i] = subspace.entry(ordering[i], unitIndex);
    }
    timeAlive = 0;
}

//Constructor with given parameters
DoubleDescriptionAlt::RayAlt::RayAlt (setbits& zeroSet, vector<Rational>& innerProductVector) {
    this->zeroSet = zeroSet;
    this->innerProductVector = innerProductVector;
    timeAlive = 0;
}

int DoubleDescriptionAlt::RayAlt::sign() {
    if(innerProductVector[timeAlive] > 0) {
        return 1;
    } else if(innerProductVector[timeAlive] < 0) {
        return -1;
    } else {
        return 0;
    }
}

vector<Ray> DoubleDescriptionAlt::reduce(const MatrixInt& subspace) {
    unsigned long rows = subspace.rows();
    unsigned long cols = subspace.columns();
    vector<Ray> constraints;
    for (int i = 0; i < rows; i++) {
        Ray row(cols);
        for(int j = 0; j < cols; j++) {
            row.setElement(j, (LargeInteger)subspace.entry(i,j));
        }
        constraints.push_back(row);
    }
    //Reduction to diagonal
    vector<int> rowOrder;
    unordered_set<int> rowSet;
    for (int i = 0; i < cols; i++) {
        int pivot = -1;
        for (int j = 0; j < rows; j++) {
            if (!rowSet.count(j) && constraints[j][i] != 0) {
                if (pivot == -1) {
                    pivot = j;
                    rowOrder.push_back(j);
                    rowSet.insert(j);
                } else {
                    Ray tempRay = Ray(constraints[pivot]);
                    tempRay *= constraints[j][i];
                    constraints[j] *= constraints[pivot][i];
                    constraints[j] -= tempRay;
                }
            }
        }              
    } 
    vector<Ray> orderedConstraints;
    for(auto index : rowOrder) {
        orderedConstraints.push_back(constraints[index]);
    }
    return orderedConstraints;   
}

bool DoubleDescriptionAlt::isAdjacentAlgebraic(vector<Ray>& constraints, RayAlt* ray1, RayAlt* ray2) {
    unsigned long rows = constraints.size();
    unsigned long cols = constraints[0].size();
    setbits zeroSet = ray1->zeroSet & ray2->zeroSet;
    vector<uint32_t> indices;
    unordered_set<uint32_t> indexSet;
    for (int j = 0; j < cols; j++) {
        if ((1 << j) & ~zeroSet) { 
            indices.push_back(j);
            indexSet.insert(j);
        } 
    }

    vector<Ray> subspace;
    for (int i = 0; i < rows; i++) {
        Ray row(indices.size());
        for(int j = 0; j < indices.size(); j++) {
            row.setElement(j, constraints[i][indices[j]]);
        }
        subspace.push_back(row);
    }    
    return true;
}

template <typename RayClass>
void DoubleDescriptionAlt::RayAlt::recover(RayClass* dest, const MatrixInt& subspace) {
    //Additional constraints -> zeroset + bounding hyperplane
    unsigned long rows = subspace.rows();
    //Count the size of the zero set
    unsigned long rem = subspace.columns(); //Number of elements in final bit set
    unsigned long zeroSetSize = 0;
    //List of indices to iterate over
    vector<uint32_t> indices;
    //Set of indices to check memebership
    unordered_set<uint32_t> indexSet;
    for (int j = 0; j < rem; j++) {
        if ((1 << j) & ~zeroSet) { 
            indices.push_back(j);
            indexSet.insert(j);
        } 
    }
    //Create the set of constraints to to solve for
    vector<Ray> constraints;
    for (int i = 0; i < rows; i++) {
        Ray row(indices.size());
        for(int j = 0; j < indices.size(); j++) {
            row.setElement(j, (LargeInteger)subspace.entry(i, indices[j]));

        }
        constraints.push_back(row);
    }
    rows = constraints.size();
    //Linear algebra -> Solve equations to rows == (col - 1)
    vector<int> rowOrder;
    unordered_set<int> rowSet;
    for (int i = 0; i < indices.size() - 1; i++) {
        int pivot = -1;
        for (int j = 0; j < rows; j++) {
            if (!rowSet.count(j) && constraints[j][i] != 0) {
                if (pivot == -1) {
                    pivot = j;
                    rowOrder.push_back(j);
                    rowSet.insert(j);
                } else {
                    Ray tempRay = Ray(constraints[pivot]);
                    tempRay *= constraints[j][i];
                    constraints[j] *= constraints[pivot][i];
                    constraints[j] -= tempRay;
                }
            }
        }              
    }
    
    unordered_map<int, Rational> solutions;
    //Set last col to 1
    solutions[indices.size() - 1] = 1;
    //Iterate over to find all solutions
    unsigned long i = 2;
    for(int k = rowOrder.size() - 1; k >= 0; k--) {
        Rational res = 0;
        for(int j = 1; j < i; j++) {
            //Multiply
            res = res + solutions[indices.size() - j] * constraints[rowOrder[k]][indices.size() - j];
        }
        //Divide
        solutions[indices.size() - i] = res / -constraints[rowOrder[k]][indices.size() - i];
        i++;
    }

    //Get LCM of all divisors and multiply through
    Vector<LargeInteger> integerSolutions(subspace.columns());
    LargeInteger lcm = 1;
    for(auto& it : solutions) {
        auto denominator = (LargeInteger)it.second.denominator();
        LargeInteger factor = gcd(denominator, lcm);
        lcm = lcm * denominator / factor;
    }
    //Actual soltion values
    unsigned long k = 0;
    for(int j = 0; j < subspace.columns(); j++) {
        if(indexSet.count(j) != 0) {
            auto numerator = (LargeInteger)solutions[k].numerator();
            auto denominator = (LargeInteger)solutions[k].denominator();
            integerSolutions.setElement(j, lcm * numerator / denominator);
            k++;
        }
    #ifdef DISPLAYANS
        cout << integerSolutions[j] << " ";
    #endif
    }
    #ifdef DISPLAYANS
    cout << endl;
    #endif
    dest = new RayClass(integerSolutions);
}

//Uses the selected hyperplane to construct the inner product vector
DoubleDescriptionAlt::RayAlt* DoubleDescriptionAlt::constructRay(RayAlt* ray1, RayAlt* ray2, int hyperPlane, 
    const MatrixInt& subspace) {
    setbits zeroSet = ray1->zeroSet & ray2->zeroSet;    
    vector<Rational> innerProductVector;
    unsigned long eqns = subspace.rows();
    //Starting indices of each ray
    unsigned long rayIndex1 = subspace.rows() - ray1->innerProductVector.size();
    unsigned long rayIndex2 = subspace.rows() - ray2->innerProductVector.size();
    for(int i = hyperPlane + 1; i < eqns; i++) {
        //ray1 = w, ray2 = u, i - rayIndex = j(in formula), hyperplane - rayIndex = i(in formular)
        Rational val = (ray1->innerProductVector[i - rayIndex1] * ray2->innerProductVector[hyperPlane - rayIndex2]
            - ray2->innerProductVector[i - rayIndex2] * ray1->innerProductVector[hyperPlane - rayIndex1])
            / (ray2->innerProductVector[hyperPlane - rayIndex2] - ray1->innerProductVector[hyperPlane - rayIndex1]);
        innerProductVector.push_back(val);
    }
    return new RayAlt(zeroSet, innerProductVector);
}

template <typename RayClass>
void DoubleDescriptionAlt::enumerateExtremalRaysAlt(const MatrixInt& subspace, 
    const EnumConstraints* constraints) {
    unsigned long eqns = subspace.rows();
    unsigned long dim = subspace.columns(); //Represents number of dimensions times 3

#ifdef TIMING
    auto start = chrono::high_resolution_clock::now();
#endif
    vector<unsigned long> ordering(eqns);
    for(int i = 0; i < eqns; i++) {
        ordering[i] = i;
    }
    //Gives hyperplane ordering
#ifdef DEBUG
    cout << "Sorting hyperplanes" << endl;
#endif
    sort(ordering.begin(), ordering.end(), LexicographicalOrder(subspace));
#ifdef TIMING
    auto sorting = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(sorting - start);
    cout << "Ordering: " << duration.count() << endl;
#endif
    //Initialise vertex set
#ifdef DEBUG
    cout << "Initialising Vertex Set" << endl;
#endif    
    vector<RayAlt*> vertexSets[2];
    for(int i = 0; i < dim; i++) {
        vertexSets[0].push_back(new RayAlt(i, subspace, ordering));
    }
#ifdef TIMING
    auto initialisation = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(initialisation - sorting);
    cout << "Initialisation: " << duration.count() << endl;
#endif

    //Intersect hyperplanes
#ifdef DEBUG
    cout << "Intersecting hyperplanes" << endl;
#endif    

#ifdef MAKEGRAPH
    //Clear log file
    ofstream myFile;
    myFile.open(LOG_FILENAME, ios::trunc);
    myFile.close();
#endif

    int currentSet = 0;
    for(int i = 0; i < eqns; i++) {
#ifdef DEBUG
    cout << "Iteration:" << i << endl;
#endif
    #ifdef TIMING
        auto itstart = chrono::high_resolution_clock::now();
    #endif          
        intersectHyperplaneAlt(i, vertexSets[currentSet], vertexSets[1 - currentSet], subspace);
        currentSet = 1 - currentSet;
    #ifdef TIMING
        auto itend = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::microseconds>(itend - itstart);
        cout << "Iteration " << i << " duration: " << duration.count() << endl;        
    #endif    
    }
#ifdef TIMING
    auto intersectingHyperplanes = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(intersectingHyperplanes - initialisation);
    cout << "Intersecting Hyperplanes: " << duration.count() << endl;
#endif

#ifdef DEBUG
    cout << "Results: " << vertexSets[currentSet].size() << endl;
#endif  
    //Output results
    for(auto ray : vertexSets[currentSet]) {
        RayClass* result;
        ray->recover(result, subspace);
        delete ray;
    }
#ifdef TIMING
    auto recovery = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(recovery - intersectingHyperplanes);
    cout << "Recovery: " << duration.count() << endl;
#endif
}

bool DoubleDescriptionAlt::isAdjacent(vector<RayAlt*>& src, RayAlt* ray1, RayAlt* ray2) {
    setbits pattern = ray1->zeroSet & ray2->zeroSet;
    for (RayAlt* ray : src) {
        if(ray == ray1 || ray == ray2) {
            continue;
        }
        if ((pattern & ray->zeroSet) == pattern) {
            return false;
        }
    }
    return true;
}

bool DoubleDescriptionAlt::isCompatible(RayAlt* ray1, RayAlt* ray2) {
    //If less bits are used, zeroSet is 11111..., negation is all zeros
    setbits pattern = ~(ray1->zeroSet & ray2->zeroSet);
    //011, 110 then 101 pattern check. Check is done on the set of ones
    if ((((pattern & FIRST_BIT) << 1) & pattern) ||
        (((pattern & SECOND_BIT) << 1) & pattern) ||
        (((pattern & FIRST_BIT) << 2) & pattern)) { 
        return false;
    } else {
        return true;
    }
}

bool DoubleDescriptionAlt::intersectHyperplaneAlt(
    int currentHyperplane, vector<RayAlt*>& src, 
    vector<RayAlt*>& dest, const MatrixInt& subspace) {
    vector<RayAlt*> poset;
    vector<RayAlt*> negset;

#ifdef MAKEGRAPH
    unordered_map<int, vector<int>> compatibilityGraph;
    vector<bool> visited(src.size(), false);
    vector<vector<int>> components;
    for (int i = 0; i < src.size(); i++) {
        for (int j = i + 1; j < src.size(); j++) {
            if (isCompatible(src[i], src[j])) {
                compatibilityGraph[i].push_back(j);
                compatibilityGraph[j].push_back(i);
            }
        }
    }

    for (int i = 0; i < src.size(); i++) {
        if (!visited[i]) {
            vector<int> component;
            getConnectedComponent(i, compatibilityGraph, visited, component);
            components.push_back(component);
        }
    }
    ofstream myFile;
    myFile.open(LOG_FILENAME, ios::app);
    myFile << "Iteration: " << currentHyperplane << endl;
    myFile << "Components: " << components.size() << endl;
    for (auto component : components) {
        for(auto num : component) {
            myFile << num << " ";
        }
        myFile << endl;
    }
    myFile << "Degree: " << endl;
    for (int i = 0; i < src.size(); i++) {
        myFile << i << " - " << compatibilityGraph[i].size() << endl;
    }
    myFile.close();
#endif

    for (auto ray : src) {
        if(ray->innerProductVector[ray->timeAlive] > 0) {
            poset.push_back(ray);
        } else if(ray->innerProductVector[ray->timeAlive] < 0) {
            negset.push_back(ray);
        } else {
            ray->timeAlive++;
            dest.push_back(ray);
        }
    } 

#ifdef DEBUG
    cout << "Filtering - Poset:" << poset.size() << " Negset:" << negset.size() << " Dest:" << dest.size() << endl;
#endif 
    //Filtering and adjacency steps
    for(auto ray1 : poset) {
        for(auto ray2 : negset) {
            if(isCompatible(ray1, ray2) && isAdjacent(src, ray1, ray2)) {
                dest.push_back(constructRay(ray1, ray2, currentHyperplane, subspace));
            }
        }
    }   
    for (auto ray : poset) {
        delete ray;
    }
    for (auto ray : negset) {
        delete ray;
    }
    src.clear();
    return true;
}

} // namespace regina

#endif