#ifndef __DOUBLEDESCRIPTION_ALT_IMPL_H
#define __DOUBLEDESCRIPTION_ALT_IMPL_H


#include <algorithm>
#include <iterator>
#include <set>
#include <vector>
#include <list>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <regina-core.h>
#include <regina-config.h>
#include <maths/integer.h>
#include <maths/ray.h>
#include <mutex>

#include "setTrie.h"
#include "rayTrie.h"

using namespace std;
//Bit size is 32 rounded down to multiple of 3
// #define LOG_FILENAME                "log.txt"
// #define MAKEGRAPH                   1
// #define DISPLAYANS
#define USE_RAYTRIE

#define HYPERPLANE_ANALYSIS
#define BIT_SIZE                    126
#define setbits                     std::bitset<128>

// #define DEBUG                       1

//Constants
const setbits FIRST_BIT (string("0001001001001001001001001001001001001001001001001001001001001001001001001001001001001001001001001001001001001001001001001001001"));
const setbits SECOND_BIT(string("0010010010010010010010010010010010010010010010010010010010010010010010010010010010010010010010010010010010010010010010010010010"));
const setbits THIRD_BIT (string("0100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100100"));

namespace regina {

#ifdef MAKEGRAPH
//Static graph in optimised version
//Visited must have space preallocated
template <typename T>
void getConnectedComponent(int node, unordered_map<T, bool>& visited, vector<T>& nodes) {
    visited[node] = true;
    nodes.push_back(node);
    for(auto adj : node->neighbours) {
        if(!visited[adj]) {
            getConnectedComponent(adj, visited, nodes);
        }
    }
}
#endif

template <typename T>
void getConnectedComponent2(int node, unordered_map<T, vector<T>>& graph, vector<bool>& visited, vector<T>& nodes, vector<bool>& combinations) {
    visited[node] = true;
    nodes.push_back(node);
    for(auto adj : graph[node]) {
        if(!visited[adj] && !combinations[adj]) {
            getConnectedComponent2(adj, graph, visited, nodes, combinations);
        }
    }
}

template <class Number>
Number gcd(Number a, Number b) {
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

template <typename T>
void subtractRow(vector<T>& pivot, vector<T>& row, int startingIndex) {
    auto gcdVal = gcd(row[startingIndex], pivot[startingIndex]);
    auto pivotMult = row[startingIndex] / gcdVal;
    auto rowMult = pivot[startingIndex] / gcdVal;
    for (int i = startingIndex; i < row.size(); i++) {
        row[i] = row[i] * rowMult - pivot[i] * pivotMult;
    }
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
        if (!zeroSet.test(j)) { 
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

//Constructor based on unit vector
DoubleDescriptionAlt::RayAlt::RayAlt (int unitIndex, const MatrixInt& subspace, vector<unsigned long>& ordering) : Ray(subspace.rows()) {
    zeroSet.set();
    zeroSet.reset(unitIndex);
    for(int i = 0; i < subspace.rows(); i++) {
        setElement(i, (LargeInteger)subspace.entry(ordering[i], unitIndex));
    }
    timeAlive = 0;
}

//Uses the selected hyperplane to construct the inner product vector
DoubleDescriptionAlt::RayAlt::RayAlt (RayAlt* ray1, RayAlt* ray2, int hyperPlane, 
    const MatrixInt& subspace) : Ray(subspace.rows()) {
    zeroSet = ray1->zeroSet & ray2->zeroSet; 
    timeAlive = hyperPlane + 1;
    for(int i = hyperPlane + 1; i < subspace.rows(); i++) {
        elements[i] = ray1->elements[i] * ray2->elements[hyperPlane]
            - ray2->elements[i] * ray1->elements[hyperPlane];
    }
    scaleDown();
    if (ray2->elements[hyperPlane] < zero) {
        negate();
    }
}

template <typename RayClass>
void DoubleDescriptionAlt::enumerateExtremalRaysAlt(const MatrixInt& subspace, 
    RunOptions options) {
    unsigned long eqns = subspace.rows();
    unsigned long dim = subspace.columns();
#ifdef HYPERPLANE_ANALYSIS
    vector<setbits> data(eqns, 0);
    for (int i = 0; i < eqns; i++) {
        for (int j = 0; j < dim; j++) {
            if(!subspace.entry(i, j).isZero()) {
                data[i].set(j);
            }
        }
    }
    unordered_map<int, vector<int>> graph;
    for (int i = 0; i < data.size(); i++) {
        for (int j = i + 1; j < data.size(); j++) {
            if ((data[i] & data[j]).any()) {
                graph[i].push_back(j);
                graph[j].push_back(i);
            }
        }    
    }
    cout << "Dimension: " << dim << endl;
    for(int k = 1; k < 6; k++) {
        vector<bool> combinations(eqns);
        for(int i = combinations.size() - k - 1; i < combinations.size(); i++) {
            combinations[i] = true;
        }
        do {
            vector<bool> visited(eqns);
            vector<vector<int>> components;
            for(int x = 0; x < data.size(); x++) {
                if (!combinations[x] && !visited[x]) {
                    vector<int> component;
                    getConnectedComponent2(x, graph, visited, component, combinations);
                }
            }
            if (components.size() >= 2) {
                cout << k << " " << components[0].size() << " " << components[1].size() << endl;
            }
        } while(std::next_permutation(combinations.begin(), combinations.end()));
    }
    return;
#endif
    vector<unsigned long> ordering(eqns);
    for(int i = 0; i < eqns; i++) {
        ordering[i] = i;
    }
    //Gives hyperplane ordering
#ifdef DEBUG
    cout << "Dimension: " << dim << endl;
    cout << "Sorting hyperplanes" << endl;
#endif
    sort(ordering.begin(), ordering.end(), LexicographicalOrder(subspace));
    //Initialise vertex set
#ifdef DEBUG
    cout << "Initialising Vertex Set" << endl;
#endif    
    vector<RayAlt*> vertexSets[2];
    for(int i = 0; i < dim; i++) {
        vertexSets[0].push_back(new RayAlt(i, subspace, ordering));
    }
    //Intersect hyperplanes
#ifdef DEBUG
    cout << "Intersecting hyperplanes" << endl;
#endif    
    int currentSet = 0;
    for(int i = 0; i < eqns; i++) {
#ifdef DEBUG
    cout << "Iteration:" << i << endl;
#endif        
        intersectHyperplaneAlt(i, vertexSets[currentSet], vertexSets[1 - currentSet], subspace, ordering, options);
        currentSet = 1 - currentSet;
    }


#ifdef DEBUG
    cout << "Results: " << vertexSets[currentSet].size() << endl;
#endif  
    //Output results
    for(auto ray : vertexSets[currentSet]) {
        RayClass* result;
        ray->recover(result, subspace);
        delete ray;
    }
}

bool DoubleDescriptionAlt::isAdjacent(const MatrixInt& subspace, SetTrie& setTrie, vector<RayAlt*>& src, RayAlt* ray1, RayAlt* ray2, 
        int currentHyperplane, vector<unsigned long>& hyperplaneOrdering, RunOptions options) {
    //Quick check for algebraic adjacency (optional)
    setbits zeroSet = ray1->zeroSet & ray2->zeroSet;
    if (currentHyperplane + zeroSet.count() + 2 < subspace.columns()) {
        return false;
    }
    switch(options.algorithm) {
        case USE_SIMPLE:
            return isAdjacentStandard(src, ray1, ray2);
        case USE_TRIE:
            return isAdjacentTrie(subspace, setTrie, ray1, ray2);
        case USE_GRAPH:
            return isAdjacentGraph(src, ray1, ray2);
        case USE_MATRIX:
            return isAdjacentAlgebraic(subspace, ray1, ray2, currentHyperplane, hyperplaneOrdering);
        default :
            return false;          
    }    
}


bool DoubleDescriptionAlt::isAdjacentStandard(vector<RayAlt*>& src, RayAlt* ray1, RayAlt* ray2) {
    setbits pattern = ray1->zeroSet & ray2->zeroSet;
    for (RayAlt* ray : src) {
        if ((pattern & ray->zeroSet) == pattern) {
            if(ray == ray1 || ray == ray2) {
                continue;
            }
            return false;
        }
    }
    return true;
}

bool DoubleDescriptionAlt::isAdjacentTrie(const MatrixInt& subspace, SetTrie& setTrie, RayAlt* ray1, RayAlt* ray2) {
    setbits zeroSet = ray1->zeroSet & ray2->zeroSet;
    vector<int> set;
    for (int j = 0; j < subspace.columns(); j++) {
        if (zeroSet.test(j)) { 
            set.push_back(j);
        }
    }
    return !setTrie.isSubsetDFS(set, reinterpret_cast<intptr_t>(ray1), reinterpret_cast<intptr_t>(ray2));
}

bool DoubleDescriptionAlt::isAdjacentGraph(vector<RayAlt*>& src, RayAlt* ray1, RayAlt* ray2) {
    setbits pattern = ray1->zeroSet & ray2->zeroSet;
    RayAlt* ray3 = (ray1->neighbours.size() < ray2->neighbours.size()) ? ray1 : ray2;
    for (auto ray : ray3->neighbours) {
        if ((pattern & ray->zeroSet) == pattern) {
            if(ray == ray1 || ray == ray2) {
                continue;
            }
            return false;
        }
    }
    return true;
}

bool DoubleDescriptionAlt::isAdjacentAlgebraic(const MatrixInt& constraints, RayAlt* ray1, RayAlt* ray2, 
        int currentHyperplane, vector<unsigned long>& hyperplaneOrdering) {
    setbits zeroSet = ray1->zeroSet & ray2->zeroSet;
    vector<uint32_t> indices;
    unordered_set<uint32_t> indexSet;
    for (int j = 0; j < constraints.columns(); j++) {
        if (!zeroSet.test(j)) { 
            indices.push_back(j);
            indexSet.insert(j);
        } 
    }
    //Make a copy of the submatrix
    vector<vector<LargeInteger>> subspace(currentHyperplane, vector<LargeInteger>(indices.size()));
    for (int i = 0; i < currentHyperplane; i++) {
        for (int j = 0; j < indices.size(); j++) {
            subspace[i][j] = (LargeInteger)constraints.entry(hyperplaneOrdering[i], indices[j]);
        }      
    }
    //Creates a map from first occupied column to rows
    unordered_map<int, vector<int>> sortBuckets;
    for (int i = 0; i < currentHyperplane; i++) {
        for (int j = 0; j < indices.size(); j++) {
            if (subspace[i][j] != 0) {
                sortBuckets[j].push_back(i);
                break;
            }
        }
    }
    //Creates an ordering for matrix (ordering[0] is the first row etc.)
    list<pair<int, int>> ordering;
    for (int i = 0; i < indices.size(); i++) {
        for (auto val : sortBuckets[i]) {
            ordering.push_back({i, val});
        }
    }

    auto currentRow = ordering.begin();
    for (int i = 0; i < indices.size() && currentRow != ordering.end(); i++) {
        //Move to next column if nothing in the column
        if (subspace[currentRow->second][i] == 0) {
            continue;
        }
        //End is index of the element 1 after
        auto next = currentRow;
        next++;
        vector<pair<int, int>> batch;
        while (next != ordering.end()) {
            if (subspace[next->second][i] != 0) {
                subtractRow(subspace[currentRow->second], subspace[next->second], i);
                for (int k = i + 1; k < indices.size(); k++) {
                    if (subspace[next->second][k] != 0) {
                        next->first = k;
                        break;
                    }
                }
                //Row is kept
                if (next->first != i) {
                    batch.push_back(*next);
                } 
                next = ordering.erase(next);
            } else {
                break;
            }
        }
        sort(batch.begin(), batch.end());
        //Fix ordering
        int batchIndex = 0;
        auto it = currentRow;
        while (it != ordering.end() && batchIndex < batch.size()) {
            if (batch[batchIndex].first < it->first) { //Prepend
                ordering.insert(it, batch[batchIndex]);
                batchIndex++;
            } else {
                it++;
            }
        }
        //Remaining elements
        while (batchIndex < batch.size()) {
            ordering.push_back(batch[batchIndex]);
            batchIndex++;
        }
        currentRow++;
    }
    return (ordering.size() == indices.size() - 2);
}

bool DoubleDescriptionAlt::isCompatible(RayAlt* ray1, RayAlt* ray2) {
    //If less bits are used, zeroSet is 11111..., negation is all zeros
    setbits pattern = (ray1->zeroSet & ray2->zeroSet).flip();
    //011, 110 then 101 pattern check. Check is done on the set of ones
    if ((((pattern & FIRST_BIT) << 1) & pattern).any() ||
        (((pattern & SECOND_BIT) << 1) & pattern).any() ||
        (((pattern & FIRST_BIT) << 2) & pattern).any()) { 
        return false;
    } else {
        return true;
    }
}

bool DoubleDescriptionAlt::intersectHyperplaneAlt(
    int currentHyperplane, vector<RayAlt*>& src, 
    vector<RayAlt*>& dest, const MatrixInt& subspace, vector<unsigned long>& hyperplaneOrdering, RunOptions options) {
    vector<RayAlt*> poset;
    vector<RayAlt*> negset;
    for (auto ray : src) {
        if((*ray)[ray->timeAlive] > 0) {
            poset.push_back(ray);
        } else if((*ray)[ray->timeAlive] < 0) {
            negset.push_back(ray);
        } else {
            ray->timeAlive++;
            dest.push_back(ray);
        }
    } 

#ifdef DEBUG
    cout << "Filtering - Poset:" << poset.size() << " Negset:" << negset.size() << " Dest:" << dest.size() << endl;
#endif 
    //Set up setTrie if needed
    SetTrie setTrie = SetTrie();
    if (options.algorithm == USE_TRIE || options.algorithm == TEST_ALL) {
        auto start = chrono::high_resolution_clock::now();
        for(auto ray : src) {
            vector<int> set;
            for (int j = 0; j < subspace.columns(); j++) {
                if (ray->zeroSet.test(j)) { 
                    set.push_back(j);
                }   
            }
            setTrie.insert(set, reinterpret_cast<intptr_t>(ray));        
        }
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        if (options.algorithm == TEST_ALL) {
            cout << "Trie: " << duration.count() << endl;
        }
    }
    //Set up graph if needed
    if (options.algorithm == USE_GRAPH || options.algorithm == TEST_ALL) {
        auto start = chrono::high_resolution_clock::now();
#ifdef USE_RAYTRIE
        //RayTrie method
        RayTrie rayTrie = RayTrie(subspace.columns() / 3);
        for (int i = 0; i < src.size(); i++) {
            rayTrie.insert(src[i]);
        }

        #pragma omp parallel for
        for (int i = 0; i < src.size(); i++) {
            rayTrie.findAll(src[i]);
        }
#else
        #pragma omp parallel for
        for (int i = 0; i < src.size(); i++) {
            for (int j = i + 1; j < src.size(); j++) {
                if (isCompatible(src[i], src[j])) {
                    src[i]->neighbours.push_back(src[j]);
                    src[j]->neighbours.push_back(src[i]);
                }
            }
        }
#endif
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
        if (options.algorithm == TEST_ALL) {
            cout << "Graph: " << duration.count() << endl;
        }
    }
#ifdef MAKEGRAPH
            unordered_map<RayAlt*, int> rayIndex;
            unordered_map<RayAlt*, bool> visited;
            vector<vector<RayAlt*>> components;
            for (int i = 0; i < src.size(); i++) {
                rayIndex[src[i]] = i;
            }
            for (int i = 0; i < src.size(); i++) {
                if (!visited[i]) {
                    vector<RayAlt*> component;
                    getConnectedComponent(i, visited, component);
                    components.push_back(component);
                }
            }
            ofstream myFile;
            myFile.open(LOG_FILENAME, ios::app);
            myFile << "Iteration: " << currentHyperplane << endl;
            myFile << "Components: " << components.size() << endl;
            for (auto component : components) {
                for(auto ray : component) {
                    myFile << rayIndex[ray] << " ";
                }
                myFile << endl;
            }
            myFile << "Degree: " << endl;
            for (int i = 0; i < src.size(); i++) {
                myFile << rayIndex[src[i]] << " - " << src[i]->neighbours.size() << endl;
            }
            myFile.close();
#endif    
    mutex destMutex;    
    if (options.algorithm != TEST_ALL) {
        auto filter = [&destMutex, &options, &setTrie, &hyperplaneOrdering, &negset, &currentHyperplane, &subspace, &src, &dest](RayAlt* ray1) {
            for(auto ray2 : negset) {
                if(isCompatible(ray1, ray2) && isAdjacent(subspace, setTrie, src, ray1, ray2, currentHyperplane, hyperplaneOrdering, options)) {
                    destMutex.lock();
                    dest.push_back(new RayAlt(ray1, ray2, currentHyperplane, subspace));
                    destMutex.unlock();
                }
            }
        };
        #pragma omp parallel for
        for (int k = 0; k < poset.size(); k++) {
            filter(poset[k]);
        }
    } else {
        cout << poset.size() << " " << negset.size() << " " << src.size() << endl;
        for (int i = 0; i < 4; i++) {
            options.algorithm = (Algorithm) i;
            vector<RayAlt*> temp;
            auto start = chrono::high_resolution_clock::now();
            #pragma omp parallel for
            for (int k = 0; k < poset.size(); k++) {
                auto ray1 = poset[k];
                for(auto ray2 : negset) {
                    if(isCompatible(ray1, ray2) && isAdjacent(subspace, setTrie, src, ray1, ray2, currentHyperplane, hyperplaneOrdering, options)) {
                        destMutex.lock();
                        if (options.algorithm == USE_TRIE) {
                            dest.push_back(new RayAlt(ray1, ray2, currentHyperplane, subspace));
                        } else {
                            temp.push_back(new RayAlt(ray1, ray2, currentHyperplane, subspace));
                        }
                        destMutex.unlock();
                    }
                }        
            }
            auto stop = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
            cout << i << " " << duration.count() << endl;
        }
    }

    if (options.algorithm == USE_GRAPH || options.algorithm == TEST_ALL) {
        for (auto ray : dest) {
            ray->neighbours.clear();
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