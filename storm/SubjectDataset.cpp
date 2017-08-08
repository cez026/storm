/*
 * SubjectDataset.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "SubjectDataset.h"
#include "QueryDataset.h"

SubjectDataset::SubjectDataset() {
	// TODO Auto-generated constructor stub
	currentFileIndex = 0;
	currentFileType = UNDEFINED;
	chunkSize = Config::getStreamChunkSize();
	outOfDataFlag = false;
	newFileFlag = true;
	currentChunkSize = 0;
	subjectFilenameList.clear();
	this->subjectReadList=NULL;
	readFilterSize = Config::minimumOverlapLength;
	this->totalProcessedSize = 0;

}

SubjectDataset::~SubjectDataset() {
	// TODO Auto-generated destructor stub

	subjectFilenameList.clear();
	if(subjectReadList!=NULL)
	{
	for(UINT64 i = 0; i < subjectReadList->size(); i++)
	{
//		cout<<"deleting: "<<i<<"/"<<subjectReadList->size()<<" "<<subjectReadList->at(i)->getName()<<endl;
		if(subjectReadList->at(i)!=NULL)
		delete subjectReadList->at(i);
	}

	subjectReadList->clear();
	delete subjectReadList;
	}
}

bool SubjectDataset::clearSubjectReadList()
{
	if(subjectReadList!=NULL)
	{
	for(UINT64 i = 0; i < subjectReadList->size(); i++)
	{
		if(subjectReadList->at(i)!=NULL)
		delete subjectReadList->at(i);
	}

	subjectReadList->clear();
	delete subjectReadList;
	subjectReadList = NULL;
	}
	return true;
}

void SubjectDataset::clearNoLinkageReads()
{
	if(subjectReadList!=NULL)
	{

		omp_set_dynamic(0);
		omp_set_num_threads(Config::numberOfThreads);
		#pragma omp parallel
			{
		#pragma omp for schedule(dynamic)
			for(UINT64 i = 0; i < subjectReadList->size(); i++)
			{
				if(subjectReadList->at(i)!=NULL&&subjectReadList->at(i)->isNoLinkage())
				{
				delete subjectReadList->at(i);
				subjectReadList->at(i) = NULL;
				}
			}
			}

	}
}

void SubjectDataset::clearNoLinkageReadsAndDelete()
{

	if(subjectReadList!=NULL&&subjectReadList->size()>0)
	{
		INT64 topElement = 0;
		INT64 bottomElement = (int)this->subjectReadList->size()-1;
		UINT64 nonEmptySpace = this->subjectReadList->size();
		//bottomElement needs to point to the first good element from the bottom
		while(bottomElement>=0&&this->subjectReadList->at(bottomElement)->isNoLinkage())
		{
			delete this->subjectReadList->at(bottomElement);
			nonEmptySpace = nonEmptySpace -1 ;
			bottomElement--;
		}
		while(bottomElement>topElement)
		{
			//check top element
			if(this->subjectReadList->at(topElement)->isNoLinkage())
			{
				delete this->subjectReadList->at(topElement);
				nonEmptySpace = nonEmptySpace -1 ;
				this->subjectReadList->at(topElement) = this->subjectReadList->at(bottomElement);
				bottomElement--;
				topElement++;
				while(bottomElement>=topElement&&this->subjectReadList->at(bottomElement)->isNoLinkage())
				{
					delete this->subjectReadList->at(bottomElement);
					bottomElement--;
					nonEmptySpace = nonEmptySpace -1 ;
				}
			}
			else
			{
				topElement++;
			}
		}

		this->subjectReadList->resize(nonEmptySpace);

	}
}

bool SubjectDataset::setFilenameList(const vector<string> & sFilenames){

	subjectFilenameList = sFilenames;
	currentFileIndex = 0;
	return true;

}

bool SubjectDataset::setReadFilterSize(UINT16 size)
{
	this->readFilterSize = size;
	return true;
}

bool SubjectDataset::addSubjectRead(SubjectRead* subjectRead)
{
	if(this->subjectReadList==NULL)
		this->subjectReadList = new vector<SubjectRead*>;
	this->subjectReadList->push_back(subjectRead);
	return true;
}

bool SubjectDataset::loadNextChunk(vector<SubjectRead*> * chunkData, QueryDataset * querydataset)
{
	if(currentFileIndex>querydataset->getNumberOfUniqueReads())return false;
	for(UINT64 k=0;k < chunkSize;k++)
	{
		currentFileIndex++;
		if(currentFileIndex>querydataset->getNumberOfUniqueReads())break;
		SubjectRead * subjectRead= new SubjectRead();
		string readname = querydataset->getReadFromID(currentFileIndex)->getName();
		string readseq = querydataset->getReadFromID(currentFileIndex)->getSequence();
		subjectRead->setName(readname);
		subjectRead->setSequence(readseq);
		chunkData->push_back(subjectRead);

	}
	return true;
}

bool SubjectDataset::loadNextChunk(vector<SubjectRead*> * chunkData)
{
	if(Config::numberOfThreads == 1)
		return loadNextChunkSingle(chunkData);
	else
		return this->loadNextChunkParallel(chunkData, Config::numberOfThreads);
}

bool SubjectDataset::loadNextChunkParallel(vector<SubjectRead*> * chunkData, UINT16 numberOfThreads)
{
	if(this->outOfDataFlag) return false;

	if(this->newFileFlag)
	{
		string currentFileName = subjectFilenameList.at(currentFileIndex);

		currentFileStreamer.open (currentFileName.c_str());
		if(currentFileStreamer == NULL)
			cout<<"Unable to open file: " << currentFileStreamer <<endl;

		if(currentFileType == UNDEFINED)
		{
			string text;
			getline (currentFileStreamer,text);
			if(text[0] == '>')
				currentFileType = FASTA;
			else if(text[0] == '@')
				currentFileType = FASTQ;
			else
				cout<< "Unknown input file format."<<endl;
			currentFileStreamer.seekg(0, ios::beg);
		}
		newFileFlag = false;
	}

	UINT64 lineOffset;
	if(currentFileType == FASTA)
		lineOffset = 2;
	else if(currentFileType == FASTQ)
		lineOffset = 4;

	vector<string>* lineData = new vector<string>();
	UINT64 thiscurrentChunkSize = 0;
	UINT64 chunkCounter;
	for(chunkCounter=0;chunkCounter < chunkSize;chunkCounter++)
	{
		// if this is end of the file, break out this while loop
		if(currentFileStreamer.eof())
		{
			currentFileIndex++;
			if(currentFileIndex>=subjectFilenameList.size())
				outOfDataFlag = true;
			else newFileFlag = true;
			currentFileStreamer.close();
			break;
		}

		//vector<string> line;
		//string line1,line0, text, line1ReverseComplement,line2ReverseComplement;
		string text="";
		if(currentFileType == FASTA) 			// Fasta file
		{
			getline(currentFileStreamer,text);
			lineData->push_back(text);
			getline(currentFileStreamer,text,'>');
			text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
			lineData->push_back(text);

		}
		else if(currentFileType == FASTQ) 					// Fastq file.
		{
			for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
			{
				getline (currentFileStreamer,text);
				lineData->push_back(text);
			}

		}
	}


	omp_set_dynamic(0);
	omp_set_num_threads(numberOfThreads);
	#pragma omp parallel shared(chunkData, thiscurrentChunkSize)
		{
	#pragma omp for schedule(dynamic)
			for(UINT64 i = 0; i < chunkCounter; i++)
			{
				string line0, line1;
				line0 = lineData->at(i*lineOffset);
				line1 = lineData->at(i*lineOffset+1);
				//remove the > if it appears in the title or name
				string readname="";
				if(line0[0] == '>'|| line0[0] == '@')
					readname = line0.substr(1);
				else readname = line0;

				for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
				    *p = toupper(*p);
				if(line1.length() > Config::getminimumOverlapLength() && QueryDataset::qualityFilter(line1) ) // Test the read is of good quality.
				{

					SubjectRead * subjectRead= new SubjectRead();
					subjectRead->setName(readname);
					string line1ReverseComplement;
					line1ReverseComplement = subjectRead->reverseComplement(line1);
					if(line1.compare(line1ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
						subjectRead->setSequence(line1);
					else
						subjectRead->setSequence(line1ReverseComplement);

				#pragma omp critical(readcount)
					{
						chunkData->push_back(subjectRead);
						thiscurrentChunkSize++;
					}
				}

			}
		}

		currentChunkSize = thiscurrentChunkSize;

		lineData->clear();
		delete lineData;

	return true;
}

bool SubjectDataset::loadNextChunkSingle(vector<SubjectRead*> * chunkData)
{
	if(this->outOfDataFlag) return false;

	if(this->newFileFlag)
	{
		string currentFileName = subjectFilenameList.at(currentFileIndex);

		currentFileStreamer.open (currentFileName.c_str());
		if(currentFileStreamer == NULL)
			cout<<"Unable to open file: " << currentFileStreamer <<endl;

		if(currentFileType == UNDEFINED)
		{
			string text;
			getline (currentFileStreamer,text);
			if(text[0] == '>')
				currentFileType = FASTA;
			else if(text[0] == '@')
				currentFileType = FASTQ;
			else
				cout<< "Unknown input file format."<<endl;
			currentFileStreamer.seekg(0, ios::beg);
		}
		newFileFlag = false;
	}

	currentChunkSize = 0;
	for(UINT64 k=0;k < chunkSize;k++)
	{
		// if this is end of the file, break out this while loop
		if(currentFileStreamer.eof())
		{
			currentFileIndex++;
			if(currentFileIndex>=subjectFilenameList.size())
				outOfDataFlag = true;
			else newFileFlag = true;
			currentFileStreamer.close();
			break;
		}

		vector<string> line;
		string line1,line0, text, line1ReverseComplement,line2ReverseComplement;
		if(currentFileType == FASTA) 			// Fasta file
		{
			getline(currentFileStreamer,text);
			line.push_back(text);
			getline(currentFileStreamer,text,'>');
			line.push_back(text);

			line.at(1).erase(std::remove(line.at(1).begin(), line.at(1).end(), '\n'), line.at(1).end());
			line1 = line.at(1);								// The first string is in the 2nd line.

			line0 = line.at(0);                              //title line is the very first line
		}
		else if(currentFileType == FASTQ) 					// Fastq file.
		{
			for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
			{
				getline (currentFileStreamer,text);
				line.push_back(text);
			}
			line1 = line.at(1); 			// The first string is in the 2nd line.
			line0 = line.at(0);             //title line is the very first line
		}

		//remove the > if it appears in the title or name
		string readname="";
		if(line0[0] == '>'|| line0[0] == '@')
			readname = line0.substr(1);
		else readname = line0;

		for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
		    *p = toupper(*p);
		if(line1.length() > Config::getminimumOverlapLength() && QueryDataset::qualityFilter(line1) ) // Test the read is of good quality.
		{

			SubjectRead * subjectRead= new SubjectRead();
			subjectRead->setName(readname);
			line1ReverseComplement = subjectRead->reverseComplement(line1);
			if(line1.compare(line1ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
				subjectRead->setSequence(line1);
			else
				subjectRead->setSequence(line1ReverseComplement);
			chunkData->push_back(subjectRead);
			currentChunkSize++;
		}

		line.clear();

	}

	return true;
}


bool SubjectDataset::loadNextChunk(vector<SubjectEdge*>* chunkData)
{
	return this->loadNextChunkParallel(chunkData, Config::numberOfThreads);
}

bool SubjectDataset::loadNextChunkParallel(vector<SubjectEdge*>* chunkData, UINT16 numberOfThreads)
{
	if(this->outOfDataFlag) return false;

//	chunkData->reserve(this->chunkSize);

	if(this->newFileFlag)
	{
		string currentFileName = subjectFilenameList.at(currentFileIndex);

		currentFileStreamer.open (currentFileName.c_str());
		if(currentFileStreamer == NULL)
			cout<<"Unable to open file: " << currentFileStreamer <<endl;

		if(currentFileType == UNDEFINED)
		{
			string text;
			getline (currentFileStreamer,text);
			if(text[0] == '>')
				currentFileType = FASTA;
			else if(text[0] == '@')
				currentFileType = FASTQ;
			else
				cout<< "Unknown input file format."<<endl;
			currentFileStreamer.seekg(0, ios::beg);
		}
		newFileFlag = false;
	}

	UINT64 lineOffset;
	if(currentFileType == FASTA)
		lineOffset = 2;
	else if(currentFileType == FASTQ)
		lineOffset = 4;

	vector<string>* lineData = new vector<string>();
	UINT64 thiscurrentChunkSize = 0;
	UINT64 chunkCounter;
	for(chunkCounter=0;chunkCounter < chunkSize;chunkCounter++)
	{
		// if this is end of the file, break out this while loop
		if(currentFileStreamer.eof())
		{
			currentFileIndex++;
			if(currentFileIndex>=subjectFilenameList.size())
				outOfDataFlag = true;
			else newFileFlag = true;
			currentFileStreamer.close();
			break;
		}

		//vector<string> line;
		//string line1,line0, text, line1ReverseComplement,line2ReverseComplement;
		string text="";
		if(currentFileType == FASTA) 			// Fasta file
		{
			getline(currentFileStreamer,text);
			lineData->push_back(text);
			getline(currentFileStreamer,text,'>');
			text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
			lineData->push_back(text);

		}
		else if(currentFileType == FASTQ) 					// Fastq file.
		{
			for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
			{
				getline (currentFileStreamer,text);
				lineData->push_back(text);
			}

		}
	}


	omp_set_dynamic(0);
	omp_set_num_threads(numberOfThreads);
	#pragma omp parallel shared(chunkData, thiscurrentChunkSize)
		{
	#pragma omp for schedule(dynamic)
			for(UINT64 i = 0; i < chunkCounter; i++)
			{
				string line0, line1;
				line0 = lineData->at(i*lineOffset);
				line1 = lineData->at(i*lineOffset+1);
				//remove the > if it appears in the title or name
				string readname="";
				if(line0[0] == '>'|| line0[0] == '@')
					readname = line0.substr(1);
				else readname = line0;

				for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
				    *p = toupper(*p);
				if(line1.length() >= this->readFilterSize && QueryDataset::qualityFilter(line1) ) // Test the read is of good quality.
				{

					SubjectRead * subjectRead= new SubjectRead();
					subjectRead->setName(readname);
					string line1ReverseComplement;
					line1ReverseComplement = subjectRead->reverseComplement(line1);
					if(line1.compare(line1ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
						subjectRead->setSequence(line1);
					else
						subjectRead->setSequence(line1ReverseComplement);

					SubjectEdge * subjectEdge = new SubjectEdge(subjectRead);


				#pragma omp critical(readcount)
					{
						chunkData->push_back(subjectEdge);
						thiscurrentChunkSize++;
					}
				}

			}
		}

		currentChunkSize = thiscurrentChunkSize;

		lineData->clear();
		delete lineData;

		UINT64 oldsize = this->totalProcessedSize;
		this->totalProcessedSize = this->totalProcessedSize+currentChunkSize;
		int oldvalue = floor(oldsize/5000000);
		int newvalue = floor(this->totalProcessedSize/5000000);
		if(newvalue - oldvalue >=1)
		{
		CURMEM;
		if(currentFileIndex<subjectFilenameList.size())
			cout << "Now reading File: " << this->subjectFilenameList.at(currentFileIndex) << endl;

		cout << setw(10) << this->totalProcessedSize << " good reads have been streamed from the data file(s)."  << endl;
		CURTIME;
		}

	return true;
}

/*
bool SubjectDataset::loadNextChunk(vector<SubjectAlignment*> & chunkData)
{
	if(outOfDataFlag) return false;

	if(newFileFlag)
	{
		string currentFileName = subjectFilenameList.at(currentFileIndex);

		currentFileStreamer.open (currentFileName.c_str());
		if(currentFileStreamer == NULL)
			cout<<"Unable to open file: " << currentFileStreamer <<endl;

		if(currentFileType == UNDEFINED)
		{
			string text;
			getline (currentFileStreamer,text);
			if(text[0] == '>')
				currentFileType = FASTA;
			else if(text[0] == '@')
				currentFileType = FASTQ;
			else
				cout<< "Unknown input file format."<<endl;
			currentFileStreamer.seekg(0, ios::beg);
		}
		newFileFlag = false;
	}

	currentChunkSize = 0;
	for(UINT64 k=0;k < chunkSize;k++)
	{
		// if this is end of the file, break out this while loop
		if(currentFileStreamer.eof())
		{
			currentFileIndex++;
			if(currentFileIndex>=subjectFilenameList.size())
				outOfDataFlag = true;
			else newFileFlag = true;
			currentFileStreamer.close();
			break;
		}

		vector<string> line;
		string line1,line0, text, line1ReverseComplement,line2ReverseComplement;
		if(currentFileType == FASTA) 			// Fasta file
		{
			getline(currentFileStreamer,text);
			line.push_back(text);
			getline(currentFileStreamer,text,'>');
			line.push_back(text);

			line.at(1).erase(std::remove(line.at(1).begin(), line.at(1).end(), '\n'), line.at(1).end());
			line1 = line.at(1);								// The first string is in the 2nd line.

			line0 = line.at(0);                              //title line is the very first line
		}
		else if(currentFileType == FASTQ) 					// Fastq file.
		{
			for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
			{
				getline (currentFileStreamer,text);
				line.push_back(text);
			}
			line1 = line.at(1); 			// The first string is in the 2nd line.
			line0 = line.at(0);             //title line is the very first line
		}

		//remove the > if it appears in the title or name
		string readname="";
		if(line0[0] == '>')
			readname = line0.substr(1);
		else readname = line0;

		for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
		    *p = toupper(*p);
		if(line1.length() > Config::getminimumOverlapLength() && QueryDataset::qualityFilter(line1) ) // Test the read is of good quality.
		{


			SubjectAlignment* subAlign = new SubjectAlignment();
			subAlign->subjectReadSequence = line1;
			subAlign->subjectReadName = readname;
			chunkData.push_back(subAlign);
			currentChunkSize++;
		}



		line.clear();

	}

	return true;
}

bool SubjectDataset::loadNextChunk(vector<edge*> & chunkData)
{
	if(outOfDataFlag) return false;

	if(newFileFlag)
	{
		string currentFileName = subjectFilenameList.at(currentFileIndex);

		currentFileStreamer.open (currentFileName.c_str());
		if(currentFileStreamer == NULL)
			cout<<"Unable to open file: " << currentFileStreamer <<endl;

		if(currentFileType == UNDEFINED)
		{
			string text;
			getline (currentFileStreamer,text);
			if(text[0] == '>')
				currentFileType = FASTA;
			else if(text[0] == '@')
				currentFileType = FASTQ;
			else
				cout<< "Unknown input file format."<<endl;
			currentFileStreamer.seekg(0, ios::beg);
		}
		newFileFlag = false;
	}

	currentChunkSize = 0;
	for(UINT64 k=0;k < chunkSize;k++)
	{
		// if this is end of the file, break out this while loop
		if(currentFileStreamer.eof())
		{
			currentFileIndex++;
			if(currentFileIndex>=subjectFilenameList.size())
				outOfDataFlag = true;
			else newFileFlag = true;
			currentFileStreamer.close();
			break;
		}

		vector<string> line;
		string line1,line0, text, line1ReverseComplement,line2ReverseComplement;
		if(currentFileType == FASTA) 			// Fasta file
		{
			getline(currentFileStreamer,text);
			line.push_back(text);
			getline(currentFileStreamer,text,'>');
			line.push_back(text);

			line.at(1).erase(std::remove(line.at(1).begin(), line.at(1).end(), '\n'), line.at(1).end());
			line1 = line.at(1);								// The first string is in the 2nd line.

			line0 = line.at(0);                              //title line is the very first line
		}
		else if(currentFileType == FASTQ) 					// Fastq file.
		{
			for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
			{
				getline (currentFileStreamer,text);
				line.push_back(text);
			}
			line1 = line.at(1); 			// The first string is in the 2nd line.
			line0 = line.at(0);             //title line is the very first line
		}

		//remove the > if it appears in the title or name
		string readname="";
		if(line0[0] == '>')
			readname = line0.substr(1);
		else readname = line0;

		for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
		    *p = toupper(*p);
		if(line1.length() > Config::getminimumOverlapLength() && QueryDataset::qualityFilter(line1) ) // Test the read is of good quality.
		{


			edge* Edge = new edge();
			Edge->subjectReadSequence = line1;
			Edge->subjectReadName = readname;
			chunkData.push_back(Edge);
			currentChunkSize++;
		}



		line.clear();

	}

	return true;
}

bool SubjectDataset::loadNextChunk(vector<SubjectAlignmentPairedEnd*> & chunkData)
{
	if(outOfDataFlag) return false;

	if(newFileFlag)
	{
		string currentFileName = subjectFilenameList.at(currentFileIndex);

		currentFileStreamer.open (currentFileName.c_str());
		if(currentFileStreamer == NULL)
			cout<<"Unable to open file: " << currentFileStreamer <<endl;

		if(currentFileType == UNDEFINED)
		{
			string text;
			getline (currentFileStreamer,text);
			if(text[0] == '>')
				currentFileType = FASTA;
			else if(text[0] == '@')
				currentFileType = FASTQ;
			else
				cout<< "Unknown input file format."<<endl;
			currentFileStreamer.seekg(0, ios::beg);
		}
		newFileFlag = false;
	}

	currentChunkSize = 0;
	for(UINT64 k=0;k < chunkSize;k++)
	{
		// if this is end of the file, break out this while loop
		if(currentFileStreamer.eof())
		{
			currentFileIndex++;
			if(currentFileIndex>=subjectFilenameList.size())
				outOfDataFlag = true;
			else newFileFlag = true;
			currentFileStreamer.close();
			break;
		}

		vector<string> line;
		string line1,line0, text, line1ReverseComplement,line2ReverseComplement;
		if(currentFileType == FASTA) 			// Fasta file
		{
			getline(currentFileStreamer,text);
			line.push_back(text);
			getline(currentFileStreamer,text,'>');
			line.push_back(text);

			line.at(1).erase(std::remove(line.at(1).begin(), line.at(1).end(), '\n'), line.at(1).end());
			line1 = line.at(1);								// The first string is in the 2nd line.

			line0 = line.at(0);                              //title line is the very first line
		}
		else if(currentFileType == FASTQ) 					// Fastq file.
		{
			for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
			{
				getline (currentFileStreamer,text);
				line.push_back(text);
			}
			line1 = line.at(1); 			// The first string is in the 2nd line.
			line0 = line.at(0);             //title line is the very first line
		}

		//remove the > if it appears in the title or name
		string readname="";
		if(line0[0] == '>')
			readname = line0.substr(1);
		else readname = line0;

		for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
		    *p = toupper(*p);
		if(line1.length() > Config::getminimumOverlapLength() && QueryDataset::qualityFilter(line1) ) // Test the read is of good quality.
		{


			SubjectAlignmentPairedEnd* subAlign = new SubjectAlignmentPairedEnd();
			subAlign->subjectReadSequence = line1;
			subAlign->subjectReadName = readname;
			chunkData.push_back(subAlign);
			currentChunkSize++;
		}



		line.clear();

	}

	return true;
}
*/
