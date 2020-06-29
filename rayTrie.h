#ifndef RAYTRIE_H
#define RAYTRIE_H

#include <vector>
#include <iostream>
#include <unordered_map>
#include <climits>
#include <set>
#include <bitset>
#include <algorithm>
#include <utility>
#include "doubledescription-alt.h"

#define setbits                     std::bitset<128>
#define ROOT                        -1

namespace regina {

class RayTrie {
    public:
        RayTrie(int size) {
            this->size = size;
            this->root = new Node();
        }

        void findAll(setbits zeroSet, DoubleDescriptionAlt::RayAlt* self) {
            std::vector<Node*> stack;
            std::unordered_map<int, int> levelToType;
            for (int i = 0; i < size; i++) { 
                for (int j = 0; j < 3; j++) { 
                    if (!zeroSet.test(3 * i + j)) {
                        levelToType[i] = j + 1;
                    }
                }
            }
            stack.push_back(root);
            while(!stack.empty()) {
                Node* node = stack.back();
                stack.pop_back();
                if (node->uid != nullptr && node->uid != self) {
                    self->neighbours.push_back(node->uid);
                }
                for (auto neighbour : node->neighbours) {
                    if (levelToType[neighbour->level] == 0 || levelToType[neighbour->level] == neighbour->type) {
                        stack.push_back(neighbour);
                    }
                }
            }
        }

        void insert(setbits zeroSet, DoubleDescriptionAlt::RayAlt* uid) {
            Node* current = root;
            for (int i = 0; i < size; i++) { //i = level
                for (int j = 0; j < 3; j++) { //j = type
                    if (!zeroSet.test(3 * i + j)) {
                        auto marker = std::find_if(current->neighbours.cbegin(), 
                                                current->neighbours.cend(),
                                                [=](const Node* n) {
                                                    return (i == n->level) && (j + 1 == n->type);
                                                });
                        //Proceed to next node if found
                        if (marker != current->neighbours.cend()) {
                            current = *marker;
                            continue;
                        }
                        //Construct and go to next node
                        Node* newNode = new Node();
                        newNode->level = i;
                        newNode->type = j + 1;
                        current->neighbours.push_back(newNode);
                        current = newNode;
                    }
                }
            }
            current->uid = uid;
        }

    private:
        struct Node {
            int level; //The corresponding level 
            int type; //1, 2 or 3
            DoubleDescriptionAlt::RayAlt* uid; //uid = 0 if not end
            std::vector<Node*> neighbours;
            Node () {
                this->level = ROOT;
                uid = nullptr;
            }
        };
        Node* root;
        int size;
};

}
#endif