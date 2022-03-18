/**
 * @file ITree.cpp
 * @author Robin Dowell 
 * @brief Interval trees on gIntervals
 * @version 0.1
 * @date 2022-03-17
 * 
 */
#include "ITree.h"

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

/**
 * @brief Construct a new Inode:: Inode object
 */
Inode::Inode() {
  left = NULL;
  right = NULL;
}

Inode::Inode(double v_center, std::vector<gInterval *> current) {
   center = v_center;

   // MidLeft: copy of {overlap intervals} sorted by left endpoint
   MidLeft = current;
   // std::sort(MidLeft.begin(), MidLeft.end(), compareStart());

   // MidRight: copy of {overlap.intervals} sorted by right endpoint
   MidRight = current;
   // std::sort(MidRight.begin(), MidRight.end(), compareEnd());

   left = right = NULL;
}
  /**
 * @brief Constructor from a set of gInterval 
 * Given a set of n intervals on the number line, we want to construct 
 * a data structure so that we can efficiently retrieve all intervals overlapping
 * another interval or point.
 * Assumes:  gInterval is a sorted list of intervals (first has smallest start; last largest stop)
 * @param segments
 */
Inode *constructTree(std::vector<gInterval *>segments){
  if (segments.empty()) return NULL;

  // center = median of interval endpoints
  double center = (double(segments[0]->stop)  + double(segments[segments.size()-1]->stop)) / 2.;

  std::vector<gInterval *> Left;
  std::vector<gInterval *> Right;
  std::vector<gInterval *> Current; 
  for (int i = 0; i < segments.size(); i++) {
      // compute iLeft: set of all intervals where right endpoint < center
      if (segments[i]->stop < center) {
          Left.push_back(segments[i]);
          // compute iRight: set of all intervals where left endpoint > center
      } else if (segments[i]->start > center) {
          Right.push_back(segments[i]);
          // compute {overlap intervals}: left endpoint <= center && right endpoint >=center
      } else {
          Current.push_back(segments[i]);
      }
  }

  // Store center and {overlap intervals} at current node
  Inode *root = new Inode(center, Current);

  // leftSubTree = constructTree(iLeft)
  if (!Left.empty()) {
    root->left 	= constructTree(Left);
  }
  // rightSubTree = constructTree(iRight)
  if (!Right.empty()) {
    root->right 	= constructTree(Right);
  }
 
  // return current node (tree)
  return root;
}

std::vector<gInterval *>Inode::searchPoint(Inode *root, double point) {
    // if query->high < center search left 
    // if query->low > center search right

    // Check all current intervals for overlap -- use left/right sorts for speedup
}

std::vector<gInterval *>Inode::overlapSearch(Inode *root, gInterval *) {
    // must look for all points in interval and concatenate (without 
    // repeats) the resulting intervals
}

void Inode::inorder(Inode *root) {
    if (root == NULL) { return;} 
    
    inorder(root->left);
    root->write_currentIntervals();
    inorder(root->right);
}

bool Inode::compareStart(gInterval *i1, gInterval *i2) {
   return (i1->start < i2->start);
}
bool Inode::compareEnd(gInterval *i1, gInterval *i2) {
   return (i1->stop < i2->stop);
}
bool Inode::Overlap(gInterval *i1, gInterval *i2) {
   return ((i1->start <= i2->stop) && (i2->start <= i1->stop));
}
std::string Inode::write_currentIntervals() {
   std::string intervals;
   intervals = "\nIntervals overlapping " + std::to_string(center);
   // iterate over contents of current node's MidLeft (i.e. intervals)
   std::vector<gInterval *>::iterator it;
   for (it = MidLeft.begin(); it != MidLeft.end(); it++) {
      intervals += "\n" + (*it)->write_out();
    }
    return intervals;
}