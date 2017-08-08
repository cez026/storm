/*
 * DoubleKeyHashTable.h
 *
 *  Created on: Feb 15, 2015
 *      Author: qy2
 */

#ifndef DOUBLEKEYHASHTABLE_H_
#define DOUBLEKEYHASHTABLE_H_

#include "Config.h"
#include "HashTableMethod.h"
#include "HashTable.h"
#include "QueryDataset.h"
#include "Alignment.h"
//#include "AlignmentPairedEnd.h"

class DoubleKeyHashTable : public HashTableMethod {

	//left key
	// 000 = 0 means prefix of the forward string.
	// 001 = 1 means suffix of the forward string.
	// 010 = 2 means prefix of the reverse string.
	// 011 = 3 means suffix of the reverse string.
	//right key
	// 100 = 4 means prefix of the forward string.
	// 101 = 5 means suffix of the forward string.
	// 110 = 6 means prefix of the reverse string.
	// 111 = 7 means suffix of the reverse string.

//	0000444444*******1111155555
//	2222666666*******3333377777
	UINT8 numberOfMode;
	UINT8 numberOfMSB;
	UINT8 numberOfLSB;

	UINT16 minimumOverlapLength;
	UINT16 hashKeyLength_left;
	UINT16 hashKeyLength_right;
	UINT8 maxMismatch;
	double maxErrorRate;



	bool wennDiagramTwoLists(vector<UINT64>* list1, vector<UINT64>* list2, vector<UINT64>* list1only, vector<UINT64>* list2only, vector<UINT64>* list12);

	string getReadSubstring(UINT64 readID, UINT8 mode);// mode ={0: forwardprefix, 1: forwardsuffix, 2: reverseprefix, 3: reversesuffix}
//	bool doAlignment(Alignment* align, string mode, int subjectStart);
//	bool checkForContainedAlignment(Alignment* align, string mode, int subjectStart);
//	bool subjectWindowRange(int& startpoint, int& stoppoint, string mode, string& subjectRead);

public:
	DoubleKeyHashTable();
	DoubleKeyHashTable(QueryDataset * qDataset);
	~DoubleKeyHashTable();
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
	bool createAlignment(Alignment* subjectAlignment, UINT8 querymode , UINT8 keymatchmode, UINT64 subjectKeyStart);
	bool doAlignmentWithSeed(Alignment* subjectAlignment,string& queryString, string& subjectString, int start, int end, int seedStart, int seedEnd);
//	bool doubleKeySearch(edge & Edge);
//	bool doubleKeySearch(SubjectAlignment & subjectAlign);
//	bool doubleKeySearch(SubjectAlignmentPairedEnd & subjectAlignment);
};

#endif /* DOUBLEKEYHASHTABLE_H_ */
