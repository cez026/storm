/*
 * PairedEndRead.h
 *
 *  Created on: Mar 27, 2015
 *      Author: qy2
 */
#include "Config.h"
#include "QueryRead.h"
#include "SubjectRead.h"


#ifndef PAIREDENDREAD_H_
#define PAIREDENDREAD_H_

class PairedEndQueryAlignment;
class CountMatrix;
class QueryRead;

class PairedEndRead {
public:
	PairedEndRead();
	~PairedEndRead();

	omp_lock_t lock;

	QueryRead * leftRead; //r1
	QueryRead * rightRead; //r2
	// Two bits are used to represent the orientation of the reads in a matepair.
	// 0 = 00 means the reverse of r1 and the reverse of r2 are matepairs.
	// 1 = 01 means the reverse of r1 and the forward of r2 are matepairs.
	// 2 = 10 means the forward of r1 and the reverse of r2 are matepairs.
	// 3 = 11 means the forward of r1 and the forward of r2 are matepairs.
	UINT8 orientation;
	string mergedSequence;

	vector<PairedEndQueryAlignment*>* pairEndAlignmentList;

	//+ATCG : the gap seqence is ATCG, no overlap between paired end reads; -ATCG: the two paired end reads are overlapped,
	// only + or - means the gap length is 0.

	map<string, int>* gapstringfrequencymap;

	map<int, CountMatrix*>* sizeToMatrix;
	map<int, int>* sizeToCount;

	void clearCountMatrix();
	void clearPairEndAlignmentList();
	void cleargapstringfrequencymap();
	bool addToPairAlignmentList(PairedEndQueryAlignment* pairedAlign);

	void initializeCountMatrix();
	bool insertMatrixToMap(int readDepthCutoff);
	bool insertgapstringfrequencymap();
	bool insertgapstringfrequencymap(string gapstring);

	bool reconstructSequence(int readDepthCutoff, double recallpercentage);
	bool reconstructSequencefromGapSeq(int readDepthCutoff, double recallpercentage);


};

//in PairedEndQueryAlignment, the paired end query reads are always forward, from left read to the right read
//subject read therefore can be either forward or reverse
class PairedEndQueryAlignment
{
public:
	PairedEndQueryAlignment(PairedEndRead* pairEndRead,SubjectRead* subjectRead);
	~PairedEndQueryAlignment();
	int subjectStart;
	int left_queryEnd;
	int subjectEnd;
	int right_queryStart;
	int right_queryEnd;
	PairedEndRead* pairEndRead;
	SubjectRead* subjectRead;
	bool subjectReadOrientation;//true:forward, false: reverse
};

class CountMatrix
{
public:
	CountMatrix(int rowSize, int ColumnSize);
	~CountMatrix();

	int matrixRowSize;
	int matrixColumnSize;
	vector<vector<UINT16>*>* matrixdata;

	bool addACTGCount(string sequence);
	bool addACTGCount(string sequence, int frequency);
	bool addACTGCount(PairedEndQueryAlignment* pairedEndQueryAlignment);
};

#endif /* PAIREDENDREAD_H_ */
