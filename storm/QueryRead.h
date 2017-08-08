/*
 * QueryRead.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef QUERYREAD_H_
#define QUERYREAD_H_

#include "Config.h"
#include "Alignment.h"
#include "PairedEndRead.h"
class Alignment;
class QueryRead {
	omp_lock_t lock;

	vector<Alignment*>* queryAlignmentList;
	string readSequence;
	string readRevSequence;
	string readName;
	UINT64 readID; 		// Unique Identification of the read. start from one. zero means a new read.This is used in hashtable.
	UINT32 frequency; 	// Frequency of the read. Number of times this read is present in the dataset. Used in some statistical analysis.
	string correctedRead;
	int newalignmentcounter;

	//paired end read information
	PairedEndRead * pairEndRead;
	UINT8 pairPosition; //0: no paired end read; 1: this is the left part; 2: this is the right part;
	void sortRightAlignments(vector<Alignment*>* alignList);
	void sortLeftAlignments(vector<Alignment*>* alignList);
	string alignmentToOmegaString(Alignment* align);
	bool checkRightAlignmentOverlap(Alignment* targetAlign,Alignment* interiorAlign);
	bool checkLeftAlignmentOverlap(Alignment* targetAlign,Alignment* interiorAlign);



public:
	bool flag4Removal; // this is a contained or duplicate read with smaller name (only keep the largest name for duplicates)
	QueryRead();
	QueryRead(string & sequence, string & name);
	~QueryRead();
	bool addAlignment(Alignment* subjectAlignment);
	bool addAlignmentWithLock(Alignment* subjectAlignment);
	bool correctErrors();
	void setRemovalTag(bool tag);
	bool needRemoval();
	void setName(string & name);
	void setSequence(string & sequence);
	void setSequence(string & sequence, string & reverseSequence);
	void setFrequency(UINT32 number);
	void setIdentifier(UINT64 id);
	string getName();
	string getSequence();
	UINT64 getIdentifier();
	UINT32 getFrequency();
	UINT32 getReadLength();
	UINT64 getAlignmentNumber();

	static string reverseComplement(const string & read);
	string reverseComplement();
	bool printAlignmentToFile(ofstream & filePointer);
	string getStringForPrintAlignmentToFile();
	string getStringForPrintNonTransitiveAlignmentToFile();
	string getStringForPrintAlignmentToFileWithoutDuplicate();

	PairedEndRead* getPairEndRead();
	void setPairEndRead(PairedEndRead* pairedEndRead);
	UINT8 getPairPosition();
	void setPairPosition(UINT8 pos);
	bool addRightMatePair(QueryRead *rightRead, UINT8 orientation);

	bool markTransitiveEdge();
	bool markTransitiveEdgeIncremental();
	void deleteTransitiveRemovableRedundantEdge();
//	bool markAndDeleteTransitiveEdge();

};

#endif /* QUERYREAD_H_ */
