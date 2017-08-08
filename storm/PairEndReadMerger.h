/*
 * PairEndReadMerger.h
 *
 *  Created on: Jun 29, 2015
 *      Author: qy2
 */

#ifndef PAIRENDREADMERGER_H_
#define PAIRENDREADMERGER_H_
#include "Config.h"
#include "HashTableMethod.h"
#include "SubjectDataset.h"

class PairEndReadMerger {

	HashTableMethod * hashTableMethod;

	int maxMismatch;
	double maxErrorRate;


	UINT16 left_minimumOverlapLength;
	UINT16 right_minimumOverlapLength;
	UINT16 overall_minimumOverlapLength;

	bool allowmismatch;

	int minimumdepth;
	double recallpercent;

	UINT16 checkKeyLength;
public:
	PairEndReadMerger(HashTableMethod * hashTableMethod);
	~PairEndReadMerger();
	bool start();
	bool createPairedEndAlignment(map<UINT64, vector<int>*>* hashkeypositionmap,Alignment* alignment, QueryRead* rightqueryread,UINT16 keylength, UINT16 rightminimumoverlap, UINT16 overallminimumoverlap);
	bool checkRightAlignment(string subjectSequence,string rightQuerySequence,int alignIndex, bool direction, int tolerateError);
	bool printMergedToFile(string notMergedFileName, string mergedFileName);
};

#endif /* PAIRENDREADMERGER_H_ */
