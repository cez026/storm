/*
 * QueryRead.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef SUBJECTREAD_H_
#define SUBJECTREAD_H_

#include "Config.h"
#include "Alignment.h"
class Alignment;
class SubjectRead {
//	omp_lock_t lock;
	string readSequence;
	string readName;
	INT64 linknum; //the number of alignment that is using this subject read

public:
	bool flag4Removal; // this is a contained or duplicate read
	SubjectRead();
	SubjectRead(string & sequence, string & name);
	~SubjectRead();
	bool needRemoval();
	void setName(string & name);
	void setSequence(string & sequence);
	string getName();
	string getSequence();
	UINT32 getReadLength();
	void setLinkNum(INT64 linknum);
	static string reverseComplement(const string & read);
	string reverseComplement();
	char reverseBaseAt(const string & read, int position);
	void increaselinknumberWithLock();
	void decreaselinknumberWithLock();
	bool isNoLinkage();
};

#endif /* SUBJECTREAD_H_ */
