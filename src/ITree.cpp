/**
 * @file ITree.cpp
 * @author Robin Dowell 
 * @brief Interval trees on bed intervals 
 * @version 0.1
 * @date 2022-03-17
 * 
 */
#include "ITree.h"

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include "Bed.h"

/**
 * @brief Construct a new Inode:: Inode object
 */
Inode::Inode() {
  left = NULL;
  right = NULL;
}

Inode::Inode(double v_center, std::vector<bed12 *> current) {
   center = v_center;
   OverlapCenter = current;

   left = right = NULL;
}

Inode::~Inode() {
  if (left != NULL) delete left;
  if (right != NULL) delete right;
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
   std::vector<bed12 *>::iterator it;
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

CITree::~CITree() {
  // delete(root);
}

CITree::CITree(std::vector<bed12 *> setIntervals) {
  constructTree(setIntervals);
}

std::vector<bed12 *>CITree::searchPoint(double point) {
  std::vector<bed12 *> hits;  // all overlapping intervals

  if (root == NULL) { return hits; } // should return empty vector?

  // Search those that overlap this center 
  std::vector<bed12 *>::iterator it;
  for (it = root->OverlapCenter.begin(); it != root->OverlapCenter.end(); it++) {
    if ((*it)->Contains(point)) { hits.push_back(*it);}
  }

  std::vector<bed12 *> rset;
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

std::vector<bed12 *>CITree::overlapSearch(bed12 *query) {
  std::vector<bed12 *> hits;  // all overlapping intervals

  bool debug = 0;   // this is manual programming switch for convenience

  if (debug) {
    std::cout << "Tree root: " << write_Root() << std::endl;
    std::cout << "overlapSearch query: " << query->write_out() << std::endl;
  }

  if (root == NULL) { return hits; } // should return empty vector?
  
  // Must consider all intervals at the current root.
  std::vector<bed12 *>::iterator it;
  for (it = root->OverlapCenter.begin(); it != root->OverlapCenter.end(); it++) {
    if ((*it)->Overlap(query)) { 
      if (debug) { std::cout << "Hit: " << (*it)->write_out() << std::endl; }
      hits.push_back(*it);
    }
  }

  // if l < c then we have to check left
  std::vector<bed12 *> rset;
  if (query->start < root->center) {
    CITree Lftleaf(root->left);
    rset = Lftleaf.overlapSearch(query);
    hits.insert(hits.end(), rset.begin(), rset.end());
  }

  // if h > c then we have to check right 
  if (query->stop > root->center) {
    CITree Rtleaf(root->right);
    rset = Rtleaf.overlapSearch(query);
    hits.insert(hits.end(), rset.begin(), rset.end());
  }

  return hits;
}
/**
 * @brief writes out the entire subtree at root
 * 
 * @param root  Centered Interval tree of interest
 * @return std::string 
 */
std::string CITree::write_Full_Tree() {
  if (root == NULL) { return "";} 

  std::string tree_output;
  std::string tstring;  // Temp variable

  // The tmp is just to not output the "L:" for empty nodes
  if (root->left != NULL) {
    CITree tmp(root->left);
    tstring = tmp.write_Full_Tree();
    if (tstring != "") { tree_output += "\nL:" + tstring; }
  }
  tree_output += root->write_currentIntervals();

  // The tmp is just to not output the "R:" for empty nodes
  if (root->right != NULL) {
    CITree tmpR(root->right);
    tstring = tmpR.write_Full_Tree();
    if (tstring != "") {
      tree_output += "\nR:" + tstring;
    }
  }

  return (tree_output);
}

std::string CITree::write_Root() {
   if (root == NULL) { return "";} 
   return root->write_currentIntervals(); 
}

/**
 * @brief Constructor from a set of bed4 
 * Given a set of n intervals on the number line, we want to construct 
 * a data structure so that we can efficiently retrieve all intervals overlapping
 * another interval or point.
 * 
 * Assumes:  bed4 is a sorted list of intervals (first has smallest stop; last largest stop)
 * @param segments
 */
Inode *CITree::constructTreeSortedSet(std::vector<bed12 *>segments){
  // center = median of interval endpoints
  int median = segments.size()/2;  // if even finds median; if odd takes upper bound on median
  double center = segments[median]->stop;

  std::vector<bed12 *> Left;
  std::vector<bed12 *> Right;
  std::vector<bed12 *> Current; 
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
    root->left 	= constructTreeSortedSet(Left);
  }
  // rightSubTree = constructTree(iRight)
  if (!Right.empty()) {
    root->right 	= constructTreeSortedSet(Right);
  }
 
  // return current node (tree)
  return root;
}

void CITree::constructTree(std::vector<bed12 *> setIntervals) {
  if (setIntervals.empty()) { root = NULL; return; }

  // Shouldn't this be responsible for all intervals originating
  // from the same chromosome??

  // Sort setIntervals by end value of intervals
  std::sort(setIntervals.begin(),setIntervals.end(), 
      [](bed12 *const &l, bed12 *const &r) { return l->stop < r->stop;});

   // Build the tree!
   root = constructTreeSortedSet(setIntervals);
}


void CITree::destroyTree() {
    delete(root);
    root = NULL;
}
