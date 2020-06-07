#ifndef SETTRIE_H
#define SETTRIE_H

#include <vector>
#include <iostream>
#include <unordered_map>
#include <set>
#include <queue>
#include <algorithm>
#include <utility>

#define ROOT                -1
#define BFS

class SetTrie {
    public:
        SetTrie() {
            root = new Node();
        }

        //Assumes list is ordered
        bool isSubset(std::set<int> set, intptr_t uid1, intptr_t uid2) {
            std::queue<Node*> q;
            int index = 0;
            int last = *set.crbegin();
            q.push(root);
            while (!q.empty()) {
                auto node = q.front();
                q.pop();
                if (node->val >= last && node->isEnd) {
                    if (node->uid == uid1 || node->uid == uid2) {
                        continue;
                    } else {
                        return true;
                    }
                }
                auto it = set.upper_bound(node->val); //Points to next value in set
                for (auto neighbour : node->neighbours) {
                    if (it == set.end() || neighbour->val <= *it) { //Valid neighbour to explore
                        q.push(neighbour);
                    }
                }
            }
            return false;
        }

        bool isSubsetDFS(std::set<int> set, intptr_t uid1, intptr_t uid2) {
            // check if set is included in trie
            const std::vector<int> sequence(set.cbegin(), set.cend());
            const size_t last_index = sequence.size() - 1;

            std::vector<std::pair<Node *, int>> stack;
            stack.emplace_back(root, 0);

            // search entire trie
            while (!stack.empty()) {
                // unpack last element on stack
                Node *node;
                int i;
                std::tie(node, i) = stack.back();
                stack.pop_back();

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
                        // check if match has been found
                        if (i == last_index && next->uid != uid1 && next->uid != uid2) {
                            return true;
                        }

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