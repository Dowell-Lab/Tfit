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

#include "Intervals.h"

/**
 * @brief Construct a new Inode:: Inode object
 */
Inode::Inode() {
  left = NULL;
  right = NULL;
}

Inode::Inode(double v_center, std::vector<gInterval *> current) {
   center = v_center;
   OverlapCenter = current;

   left = right = NULL;
}

/**
 * @brief Output function for THIS node.
 * 
 * @return std::string 
 */
std::string Inode::write_currentIntervals() {
   std::string intervals;
   intervals = "\nIntervals overlapping " + std::to_string(center);
   // iterate over contents of current node's intervals

   std::vector<gInterval *>::iterator it;
   for (it = OverlapCenter.begin(); it != OverlapCenter.end(); it++) {
      intervals += "\n" + (*it)->write_out();
    }
    return intervals;
}

/******************************************************/

CITree::CITree() {
  root = NULL;
}

CITree::CITree(Inode *v_root) {
  root = v_root;
}

CITree::CITree(std::vector<gInterval *> setIntervals) {
   if (setIntervals.empty()) { root = NULL; return; }

  // Sort setIntervals by end value of intervals
  std::sort(setIntervals.begin(),setIntervals.end(), 
      [](gInterval *const &l, gInterval *const &r) { return l->stop < r->stop;});

   // Build the tree!
   root = constructTree(setIntervals);
}

std::vector<gInterval *>CITree::searchPoint(double point) {
  std::vector<gInterval *> hits;  // all overlapping intervals

  if (root == NULL) { return hits; } // should return empty vector?

  // Search those that overlap this center 
  std::vector<gInterval *>::iterator it;
  for (it = root->OverlapCenter.begin(); it != root->OverlapCenter.end(); it++) {
    if ((*it)->Contains(point)) { hits.push_back(*it);}
  }

  std::vector<gInterval *> rset;
  // We must also search half the children
  if (point < root->center) {
    CITree Lftleaf(root->left);
    rset = Lftleaf.searchPoint(point);
  } else {
    CITree Rtleaf(root->right);
    rset = Rtleaf.searchPoint(point);
  }
  hits.insert(hits.end(), rset.begin(), rset.end());

  return hits;
}

/*
std::vector<gInterval *>CITree::overlapSearch(gInterval *) {
    // must look for all points in interval and concatenate (without 
    // repeats) the resulting intervals
}
*/

/**
 * @brief writes out the entire subtree at root
 * 
 * @param root  Centered Interval tree of interest
 * @return std::string 
 */
std::string CITree::write_Full_Tree() {
  if (root == NULL) { return "";} 

  std::string tree_output;

  // The tmp is just to not output the "L:" for empty nodes
  CITree tmp(root->left);
  std::string tstring = tmp.write_Full_Tree();
  if (tstring != "") { tree_output += "\nL:" + tstring; }

  tree_output += root->write_currentIntervals();

  // The tmp is just to not output the "R:" for empty nodes
  CITree tmpR(root->right);
  tstring = tmpR.write_Full_Tree();
  if (tstring != "") { tree_output += "\nR:" + tstring; }

  return (tree_output);
}

/**
 * @brief Constructor from a set of gInterval 
 * Given a set of n intervals on the number line, we want to construct 
 * a data structure so that we can efficiently retrieve all intervals overlapping
 * another interval or point.
 * 
 * Assumes:  gInterval is a sorted list of intervals (first has smallest stop; last largest stop)
 * @param segments
 */
Inode *CITree::constructTree(std::vector<gInterval *>segments){
  // center = median of interval endpoints
  int median = segments.size()/2;  // if even finds median; if odd takes upper bound on median
  // std::cout << median << std::endl;
  double center = segments[median]->stop;

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

