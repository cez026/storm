/*
 * OverlapGraphConstructor.h
 *
 *  Created on: Mar 5, 2015
 *      Author: qy2
 */

#ifndef OVERLAPGRAPHCONSTRUCTOR_H_
#define OVERLAPGRAPHCONSTRUCTOR_H_

#include "Config.h"
#include "SingleKeyHashTable.h"
#include "DoubleKeyHashTable.h"
#include "QueryDataset.h"
#include "HashTableMethod.h"
#include "SubjectDataset.h"
#include "Alignment.h"


class OverlapGraphConstructor {

	HashTableMethod * hashTableMethod;

	UINT64 totaledgenumber;
	ofstream filePointer;
	UINT16 minimumOverlapLength;
	bool isContainedAlignment(Alignment * subjectAlignment);
public:
	OverlapGraphConstructor(HashTableMethod * hashTableMethod);
	~OverlapGraphConstructor();
	bool start();
	bool start(bool transitiveEdgeRemoval);
	bool printEdgesToFile(string outFileName);
	bool printEdgesToFileWithoutTransitiveEdge(string outFileName);
	bool printEdgesToFileWithoutDuplicate(string outFileName);

};

#endif /* OVERLAPGRAPHCONSTRUCTOR_H_ */
