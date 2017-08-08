/*
 * SubjectRead.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "QueryRead.h"

SubjectRead::SubjectRead() {

	readSequence = "";
	readName = "";
	linknum = 0;
	flag4Removal = false;
//	omp_init_lock(&lock);

}

SubjectRead::SubjectRead(string& sequence, string& name)
{
	readSequence = sequence;
	readName = name;
	linknum = 0;
	flag4Removal = false;
//	omp_init_lock(&lock);

}

SubjectRead::~SubjectRead() {
	// TODO Auto-generated destructor stub
//	omp_destroy_lock(&lock);
}

bool SubjectRead::needRemoval()
{
	return flag4Removal;
}
void SubjectRead::setName(string & name)
{
	readName = name;
}

void SubjectRead::setLinkNum(INT64 linknum)
{
	this->linknum = linknum;
}
void SubjectRead::setSequence(string & sequence)
{
	readSequence = sequence;
}

string SubjectRead::getName()
{
	return readName;
}
string SubjectRead::getSequence()
{
	return readSequence;
}


UINT32 SubjectRead::getReadLength()
{
	return this->readSequence.length();
}


string SubjectRead::reverseComplement(const string & read)
{
	UINT64 stringLength = read.length();
	string reverse(stringLength,'0');
	for(UINT64 i = 0;i < stringLength; i++)						// Then complement the string. Change A to T, C to G, G to C and T to A.
	{
		if( read[i] & 0X02 ) // C <==> G
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X04;
		else // A <==> T
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X15;
	}
	return reverse; // return the reverse complement as a string
}

char SubjectRead::reverseBaseAt(const string & read, int position)
{
	if(position<0 || position>=(int)read.length()) return 0;
	char reversebase = 0;
	if( read[position] & 0X02 ) // C <==> G
		reversebase = read[position] ^ 0X04;
	else // A <==> T
		reversebase = read[position] ^ 0X15;
	return reversebase;

}

string SubjectRead::reverseComplement()
{
	return QueryRead::reverseComplement(this->readSequence);
}

void SubjectRead::increaselinknumberWithLock()
{
//	omp_set_lock(&lock);
//	this->linknum=this->linknum+1;
//	omp_unset_lock(&lock);
#pragma omp atomic
	this->linknum++;
}
void SubjectRead::decreaselinknumberWithLock()
{
//	omp_set_lock(&lock);
//	if(this->linknum>0)this->linknum--;
//	omp_unset_lock(&lock);
#pragma omp atomic
	this->linknum--;
}
bool SubjectRead::isNoLinkage()
{
	return this->linknum<=0;
}

