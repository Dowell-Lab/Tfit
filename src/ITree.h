/**
 * @file ITree.h
 * @author Robin Dowell
 * @brief Header file for ITree classes
 * @version 0.1
 */
#ifndef ITree_H
#define ITree_H

#include <string>
#include <vector>

#include "Intervals.h"

/**
 * @brief A node within the centered interval tree.
 * Centered interval trees use binary search trees where each 
 * subtree is also a centered interval tree.   This approach has
 * faster query times O(log n + k) than the augmented red/black 
 * based interval trees (which query in O(log n * k)). 
 */
class Inode{
public:
	double center;
	Inode * left;  //!< All intervals fully to left of center
	Inode * right; //!< All intervals fully to the right of center

	std::vector<gInterval *> OverlapCenter;	//<! All intervals overlapping center, sorted by end 

	// Constructors
	Inode();	// empty constructor
    Inode(double, std::vector<gInterval *>);

	/* FUNCTIONS: */
    std::string write_currentIntervals(); // print intervals on current node.

};

/**
 * @brief The actual tree object only contains a pointer to the root.
 * 
 */
class CITree {
    Inode *root;

public:
	// Constructors
    CITree();
    CITree(Inode *);
    CITree(std::vector<gInterval *>);  // sorts intervals then builds tree.

    // Functions:
    std::string write_Full_Tree();
    // Search tree for all intervals that overlap a given point 
    std::vector<gInterval *>searchPoint(double);
    // Search tree for all intervals that overlap a given interval 
    std::vector<gInterval *>overlapSearch(gInterval *);

    // Could create insert/delete functions for tree?

private: 
    // Recursively builds tree, assumes sorted vector of intervals
    Inode *constructTree(std::vector<gInterval *>);
};

#endif