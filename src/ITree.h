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
    // Destructor
    ~Inode();

	/* FUNCTIONS: */
    std::string write_currentIntervals(); // print intervals on current node.

};

/**
 * @brief The actual tree object only contains a pointer to the root.
 * 
 */
class CITree {
public:
    Inode *root;

	// Constructors
    CITree();
    CITree(Inode *);
    CITree(std::vector<gInterval *>);  // sorts intervals then builds tree.
    //Destructor
    ~CITree();

    // Functions:
    // Debugging functions:
    std::string write_Full_Tree();
    std::string write_Root();

    // Search tree for all intervals that overlap a given point 
    std::vector<gInterval *>searchPoint(double);
    // Search tree for all intervals that overlap a given interval 
    std::vector<gInterval *>overlapSearch(gInterval *);

    // Recursively builds tree, assumes sorted vector of intervals
    Inode *constructTree(std::vector<gInterval *>);
    void destroyTree();
};

#endif
