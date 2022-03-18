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

	std::vector<gInterval *> MidLeft;	//<! All intervals overlapping center, sorted by left
	std::vector<gInterval *> MidRight;	//<! All intervals overlapping center, sorted by right

	// Constructors
	Inode();	// empty constructor
    Inode(double, std::vector<gInterval *>);

	/* FUNCTIONS: */
    // Construct Interval Tree
	Inode *constructTree(std::vector<gInterval *>);
    // Search for all intervals that overlap a given point 
    std::vector<gInterval *>searchPoint(Inode *, double);
    // Search for all intervals that overlap a given interval 
    std::vector<gInterval *>overlapSearch(Inode *, gInterval *);

    // Delete a node from the tree
    // Insert a new interval into the tree

    // Print a Tree
    void inorder(Inode *); // In order traversal of tree (output function)

    // Helper functions:
    bool Overlap(gInterval *, gInterval *);
    bool compareStart(gInterval *, gInterval *);
    bool compareEnd(gInterval *, gInterval *);
    std::string write_currentIntervals(); // print intervals on current node.

};



#endif
