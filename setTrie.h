#ifndef SETTRIE_H
#define SETTRIE_H

#include <vector>
#include <iostream>
#include <unordered_map>
#include <climits>
#include <set>
#include <queue>
#include <algorithm>
#include <utility>

#define ROOT                -1

class SetTrie {
    public:
        SetTrie() {
            root = new Node();
        }

        //Assumes list is ordered
        bool isSubset(std::vector<int> sequence, intptr_t uid1, intptr_t uid2) {
            std::queue<std::pair<Node*, int>> q;
            int index = 0;
            int last = sequence.back();
            sequence.push_back(INT_MAX);
            q.push({root, 0});
            while (!q.empty()) {
                auto node = q.front().first;
                int i = q.front().second;
                q.pop();
                if (node->val >= last && node->isEnd) {
                    if (node->uid == uid1 || node->uid == uid2) {
                        continue;
                    } else {
                        return true;
                    }
                }

                for (auto neighbour : node->neighbours) {
                    if (neighbour->val <= sequence[i]) {
                        if(neighbour->val == sequence[i]) {
                            q.push({neighbour, i + 1});
                        } else {
                            q.push({neighbour, i});
                        }
                    }
                }
            }
            return false;
        }

        bool isSubsetDFS(std::vector<int> sequence, intptr_t uid1, intptr_t uid2) {
            std::vector<std::pair<Node *, int>> stack;
            stack.emplace_back(root, 0);

            // search entire trie
            while (!stack.empty()) {
                // unpack last element on stack
                Node *node;
                int i;
                std::tie(node, i) = stack.back();
                stack.pop_back();

                // check if match has been found
                if (i >= sequence.size()) {
                    if (node->uid != uid1 && node->uid != uid2 && node->isEnd) {
                        return true;
                    }

                    continue;
                }

                // iterate through neighbours
                for (Node *next : node->neighbours) {
                    // check if next value is less
                    if (next->val < sequence[i]) {
                        // add to stack with same index
                        stack.emplace_back(next, i);
                        continue;
                    }

                    // check if next value is equal
                    if (next->val == sequence[i]) {
                        // add to stack with incremented index
                        stack.emplace_back(next, i + 1);
                        continue;
                    }
                }
            }

            return false;
        }

        void insert(std::vector<int> set, intptr_t uid) {
            // insert elements of set into Node neighbours
            Node* node = root;
            // iterate through set of elements
            for (int elem : set) {
                // search for element in neighbours
                auto marker = std::find_if(node->neighbours.cbegin(),
                                           node->neighbours.cend(),
                                           [=](const Node* n) {
                                               return elem == n->val;
                                           });

                // if elem is found, continue search
                if (marker != node->neighbours.cend()) {
                    node = *marker;
                    continue;
                }

                // if elem is not found, add new node
                Node* newNode = new Node();
                newNode->val = elem;

                node->neighbours.push_back(newNode);
                node = newNode;
            }
            // set isEnd of last Node
            node->isEnd = true;
            node->uid = uid;
        }

    private:
        struct Node {
            int val;
            intptr_t uid;
            std::vector<Node*> neighbours;
            bool isEnd;
            Node () {
                this->val = ROOT;
                isEnd = false;
                uid = 0;
            }
        };
        Node* root;
};

#endif