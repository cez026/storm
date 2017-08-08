/*
 * PairedEndReadsMerger.h
 *
 *  Created on: Mar 29, 2015
 *      Author: qy2
 */

#ifndef PAIREDENDREADSMERGER_H_
#define PAIREDENDREADSMERGER_H_
#include "Config.h"
#include "HashTableMethod.h"
#include "SubjectDataset.h"

class PairedEndReadsMerger {

	HashTableMethod * hashTableMethod;

	int maxMismatch;
	double maxErrorRate;


	UINT16 left_minimumOverlapLength;
	UINT16 right_minimumOverlapLength;
	UINT16 overall_minimumOverlapLength;

	bool start_errornum();
	bool start_errorrate();
	bool start_errornum_seqfree();//an alternative implementation to avoid storing the subject sequence
	bool start_errorrate_seqfree();//an alternative implementation to avoid storing the subject sequence

public:
	PairedEndReadsMerger(HashTableMethod * hashTableMethod);
	~PairedEndReadsMerger();
	int checkAlignment(string subjectString, string queryString, UINT16 rightminimumoverlap, UINT16 overallminimumoverlap,int tolerateErrors, bool direction);
	int checkAlignment(string subjectString, string queryString, UINT16 rightminimumoverlap, UINT16 overallminimumoverlap,double tolerateErrorRate, bool direction);

	bool start();
	bool searchHashTable(SubjectEdge * subjectEdge);
	bool printMergedToFile(string notMergedFileName, string mergedFileName);
};

#endif /* PAIREDENDREADSMERGER_H_ */
