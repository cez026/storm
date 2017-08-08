/*
 * SingleKeyHashTable.h
 *
 *  Created on: Jan 6, 2015
 *      Author: qy2
 */

#ifndef SINGLEKEYHASHTABLE_H_
#define SINGLEKEYHASHTABLE_H_

#include "Config.h"
#include "HashTableMethod.h"
#include "HashTable.h"
#include "QueryDataset.h"
#include "Alignment.h"
//#include "AlignmentPairedEnd.h"

class SingleKeyHashTable : public HashTableMethod {

	// 00 = 0 means prefix of the forward string.
	// 01 = 1 means suffix of the forward string.
	// 10 = 2 means prefix of the reverse string.
	// 11 = 3 means suffix of the reverse string.
	UINT8 numberOfMode;
	UINT8 numberOfMSB;
	UINT8 numberOfLSB;

	UINT16 minimumOverlapLength;
	UINT16 hashKeyLength;
	UINT8 maxMismatch;
	double maxErrorRate;




	string getReadSubstring(UINT64 readID, UINT8 mode);// mode ={0: forwardprefix, 1: forwardsuffix, 2: reverseprefix, 3: reversesuffix}
//	bool doAlignment(Alignment* align, string mode, int subjectStart);
//	bool checkForContainedAlignment(Alignment* align, string mode, int subjectStart);
//	bool subjectWindowRange(int& startpoint, int& stoppoint, string mode, string& subjectRead);

public:
	SingleKeyHashTable();
	SingleKeyHashTable(QueryDataset * qDataset);
	~SingleKeyHashTable();
	bool createHashTables();
	bool insertQueryDataset(QueryDataset* d);
	bool insertQueryDataset_leftReadForMerge(QueryDataset* d);
	bool insertQueryDataset_rightAlignOnly(QueryDataset* d);
	bool insertQueryRead(QueryRead *read, string subString, UINT8 mode);
	vector<UINT64> * getListOfReads(string subString);
	vector<UINT64> * getListOfReads(string subString, UINT64& hashValue);
	bool searchHashTable(SubjectEdge * subjectEdge);
	bool searchHashTableFullScan(SubjectEdge * subjectEdge);
	bool searchHashTable(SubjectEdge * subjectEdge, bool fullScan);
	bool searchHashTableFullScanAndMarkQueryRead(SubjectEdge * subjectEdge);
	bool searchHashTableFullScanAndRecordKeys(SubjectEdge * subjectEdge);
	bool createAlignment(Alignment* subjectAlignment, UINT8 queryMode, UINT64 subjectKeyStart);
	bool doAlignment(Alignment* subjectAlignment,string& queryString, string& subjectString, int remainStart, int remainEnd);
//	bool singleKeySearch(edge & Edge);
//	bool singleKeySearch(SubjectAlignment & subjectAlign);
//	bool singleKeySearch(SubjectAlignmentPairedEnd & subjectAlignment);
};

#endif /* SINGLEKEYHASHTABLE_H_ */
