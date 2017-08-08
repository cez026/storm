/*
 * Alignment.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include "Config.h"
#include "QueryRead.h"
#include "SubjectRead.h"

class QueryRead;




class Alignment {



public:
	Alignment();
	Alignment(SubjectRead* sRead);
	Alignment(QueryRead* qRead);
	Alignment(SubjectRead* sRead, QueryRead* qRead);
	~Alignment();

	/*
	// Don't need these two variables for error correction to save memory.
	string subjectReadSequence;
	string subjectReadName;
	*/
	SubjectRead * subjectRead;
	//alignment is subject to the subject reads, which is on top of the orientation graph

	QueryRead * queryRead;

	// orientation of the query read
	bool queryOrientation;		// true for forward, false for reverse
	bool subjectOrientation;

	// coordinates of the overlap alignment, which is defined by the query reads
	//  The following case will have a positive subject start position
	//  query:        XXXXXXMMMMMMMM
	//	subject:            MMMMMMMMXXXXXXX
	//  The following case will have a negative subject start position
	//  query:        MMMMMMXXXXXXXX
	//	subject: XXXXXMMMMMM
	int subjectStart;
	int queryEnd;
	int subjectEnd;

	UINT8 containedTag; //0-not contained alignment, 1-subject is contained in query 2-query is contained in subject
//128-transitive edge


	int orientationTranslate();
	bool isContainedAlignment();
	bool isLeftOrRightAlignment(bool& isRight, int& startPoint);
	int getEditDistance(); //return the size of the editInfor
	bool insertSubstitution(int position, char base);
	char reverseComplement(char base);
	string editInfoToOmegaString();
	void reverseAlignment();
	int reverseCoordinate(int oldcoordinate);
	string getSubjectSubSequence(int alignstart, int alignend);//alignstart, alignend are coordinates from our alignment system
	int distanceBetweenEditInfo(Alignment* align);

private:
	// coordinates are defined by the query reads
	// all insertion and deletion is on the subject read
	// upper case char means substitution
	// lower case char means insertion
	// 'D' means deletion
	//
	map<int, char>* editInfor;
	int calculateWenn(map<int, char>* editInfo1, map<int, char>* editInfo2);
};

class Edge{
public:
	QueryRead* source;
	QueryRead* destination;
	Edge(){source = NULL; destination = NULL;};
	~Edge(){};
};

class ContainedAlignment {



public:

	ContainedAlignment(SubjectRead* sRead, QueryRead* qRead);
	~ContainedAlignment();
	char reverseComplement(char base);
	/*
	// Don't need these two variables for error correction to save memory.
	string subjectReadSequence;
	string subjectReadName;
	*/
	SubjectRead * subjectRead;

	QueryRead * queryRead;
};



class SubjectEdge{
public:
	// alignment input
	SubjectRead* subjectRead;

	//hash key storage
	map<UINT64, vector<int>*>* hashkeypositionmap; //key in single hash table method, or left key in double hash table method

	// alignment results
	vector<Alignment*>* alignmentList;
	vector<ContainedAlignment*>* contained_alignmentList;

	// a list of query reads that are duplicated with the subject read
	// when a query is duplicate with a subject read, we compare their names and add the query read to the duplicate read list if its name is greater than the subject read name.
	vector<QueryRead *>* DuplicateReadList;
	SubjectEdge(SubjectRead* sRead);
	~SubjectEdge();
	bool addAlignment(Alignment* subjectAlignment);
	bool addContainedAlignment(ContainedAlignment* subjectAlignment);
	bool addDuplicateList(QueryRead* queryRead);
	bool addKeyToMap(UINT64 keyvalue, int position);
	void clearAlignmentList();


};

#endif /* ALIGNMENT_H_ */
