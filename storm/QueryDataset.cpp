/*
 * QueryDataset.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "QueryDataset.h"

QueryDataset::QueryDataset()
{
	// TODO Auto-generated constructor stub
//	cout << "constructor " << queryReadList.size() << endl;
	numberOfReads = 0;								// Number of total reads present in the dataset.
	numberOfUniqueReads = 0; 						// number of unique reads in the dataset.

	shortestReadLength = 0xFFFFFFFFFFFFFFFF;
	longestReadLength = 0;
	queryReadList=NULL;
	pairedEndReadList=NULL;

	dataset_minimumoverlaplength = Config::getminimumOverlapLength();
	dataset_QueryFilename = Config::getQueryDatasetFilename();
	dataset_leftminimumoverlaplength = Config::getleftminimumOverlapLength();
	dataset_rightminimumoverlaplength = Config::getrightminimumOverlapLength();
}

QueryDataset::QueryDataset(const string & QueryFilename)
{
	// TODO Auto-generated constructor stub
//	cout << "constructor " << queryReadList.size() << endl;
	numberOfReads = 0;								// Number of total reads present in the dataset.
	numberOfUniqueReads = 0; 						// number of unique reads in the dataset.

	shortestReadLength = 0xFFFFFFFFFFFFFFFF;
	longestReadLength = 0;
	queryReadList=NULL;
	pairedEndReadList=NULL;

	dataset_minimumoverlaplength = Config::getminimumOverlapLength();
	dataset_QueryFilename = QueryFilename;
	dataset_leftminimumoverlaplength = Config::getleftminimumOverlapLength();
	dataset_rightminimumoverlaplength = Config::getrightminimumOverlapLength();
}

QueryDataset::~QueryDataset()
{
	if(queryReadList!=NULL)
	{
	for(UINT64 i = 0; i < queryReadList->size(); i++)
	{
		delete queryReadList->at(i);
	}
	queryReadList->clear();
	delete queryReadList;
	}
	if(pairedEndReadList!=NULL)
	{
		for(UINT64 i = 0; i < pairedEndReadList->size(); i++)
		{
			delete pairedEndReadList->at(i);
		}
		pairedEndReadList->clear();
		delete pairedEndReadList;
	}
}

UINT64 QueryDataset::getNumberOfReads()
{
	return this->numberOfReads;
}
UINT64 QueryDataset::getNumberOfUniqueReads()
{
	return this->numberOfUniqueReads;
}

QueryRead * QueryDataset::getReadFromID(UINT64 ID)
{
	if(queryReadList!=NULL)
	{
		if(ID < 1 || ID > numberOfUniqueReads)	// ID outside the range.
		{

			cout << "ID " << ID << " out of bound."<<endl;
			return NULL;

		}
		else return queryReadList->at(ID - 1);
	}
	else return NULL;
}

/*
QueryRead * QueryDataset::getReadFromString(const string & read)
{
	if(queryReadList!=NULL)
	{
		UINT64 min = 0, max = getNumberOfUniqueReads()-1;
		string readReverse = QueryRead::reverseComplement(read);
		int comparator;
		if(read.compare(readReverse) < 0)
		{
			while (max >= min) 															// At first search for the forward string.
			{
				UINT64 mid = (min + max) / 2; 	// Determine which subarray to search.
				comparator = queryReadList->at(mid)->getSequence().compare(read.c_str());
				if(comparator == 0)
					return queryReadList->at(mid);
				else if (comparator < 0) 	// Change min index to search upper subarray.
					min = mid + 1;
				else if (comparator > 0) 	// Change max index to search lower subarray.
					max = mid - 1;
			}
		}
		else
		{
			while (max >= min) 																	// If forward string is not found then search for the reverse string
			{
				UINT64 mid = (min+max) / 2; 													// Determine which subarray to search
				comparator = queryReadList->at(mid)->getSequence().compare(readReverse.c_str());
				if( comparator == 0)
					return queryReadList->at(mid);
				else if (comparator < 0) 	// Change min index to search upper subarray.
					min = mid + 1;
				else if (comparator > 0) 	// Change max index to search lower subarray.
					max = mid - 1;
			}
		}
		cout<<"String not found "<<endl;
		return NULL;
	}
	else
	{
		cout<<"no reads in the dataset";
		return NULL;
	}
}
*/

bool QueryDataset::qualityFilter(string & sequence)
{

		//Returns true if the read contains only {A,C,G,T} and does not contain more than 80% of the same nucleotide
		UINT64 cnt[4] = {0,0,0,0};
		UINT64 readLength = sequence.length();
		for(UINT64 i = 0; i < readLength; i++) // Count the number of A's, C's , G's and T's in the string.
		{
			if(sequence[i]!= 'A' && sequence[i] != 'C' && sequence[i] != 'G' && sequence[i] != 'T')
				return false;
			cnt[(sequence[i] >> 1) & 0X03]++; // Least significant 2nd and 3rd bits of ASCII value used here
		}
		UINT64 threshold = sequence.length()*.8;	// 80% of the length.
		if(cnt[0] >= threshold || cnt[1] >= threshold || cnt[2] >= threshold || cnt[3] >= threshold)
			return false;	// If 80% bases are the same base.
		return true;

}
/*
bool QueryDataset::qualityFilter(string & sequence)
{
	UINT64 readLength = sequence.length();
		for(UINT64 i = 0; i < readLength; i++) // Count the number of A's, C's , G's and T's in the string.
		{
			if(sequence[i]!= 'A' && sequence[i] != 'C' && sequence[i] != 'G' && sequence[i] != 'T')
				return false;
		}
		return true;

}
*/
bool QueryDataset::duplicateFilter()
{
	if(queryReadList!=NULL)
	{
		UINT64 j = 0;
		QueryRead *temp;
		for(UINT64 i = 0; i < queryReadList->size(); i++)		// Move the unique reads in the top of the sorted list. Store the frequencey of the duplicated reads.
		{
			if(queryReadList->at(j)->getSequence()!= queryReadList->at(i)->getSequence())
			{
				j++;
				temp = queryReadList->at(j);
				queryReadList->at(j) = queryReadList->at(i);
				queryReadList->at(i) = temp;
			}
			else if(i!=j)
			{
				queryReadList->at(j)->setFrequency(queryReadList->at(j)->getFrequency() + 1);
				//set the name as the alphabetical larger name
				string newname = queryReadList->at(i)->getName();
				if(queryReadList->at(j)->getName()< newname)
					queryReadList->at(j)->setName(newname);
			}
		}
		numberOfUniqueReads = j+1;
		cout <<"Number of unique reads: " << numberOfUniqueReads << endl;
		for(UINT64 i = 0 ; i < queryReadList->size(); i++) 		// Assign ID's to the reads.
		{
			if(i < getNumberOfUniqueReads())
				queryReadList->at(i)->setIdentifier(i + 1);
			else
				delete queryReadList->at(i); 					// Free the unused reads.
		}
		queryReadList->resize(numberOfUniqueReads);				//Resize the list.

		return true;
	}
	else
	{
		cout<<"no reads in the dataset";
		return false;
	}
}



bool compareReads (QueryRead *read1, QueryRead *read2)
{
	return read1->getSequence() < read2->getSequence();
}
void QueryDataset::sortReads()
{
	if(queryReadList!=NULL)
	std::sort(queryReadList->begin(),queryReadList->end(), compareReads);	// Sort the reads lexicographically.
	else
	cout<<"no reads in the dataset";

}

bool QueryDataset::buildDataset()
{if(Config::numberOfThreads == 1)
	return this->buildDataset(this->dataset_QueryFilename);
else
	return this->buildDatasetParallel(this->dataset_QueryFilename, Config::numberOfThreads);
}



bool QueryDataset::buildDatasetParallel(const string & QueryFilename, UINT16 numberOfThreads)
{
	if(queryReadList==NULL)
		queryReadList = new vector<QueryRead *>();
	cout << "initial " << queryReadList->size() << endl;

	UINT64 goodReads = 0, badReads = 0;
	UINT64 totalNumberOfReads = this->numberOfReads;
	vector<QueryRead *>*& thisqueryReadList = this->queryReadList;
	string fileName = QueryFilename;
	cout << "Reading dataset: "  << " from file: " << fileName << endl;
	ifstream myFile;
	myFile.open (fileName.c_str());
	if(myFile == NULL)
		{cout<<"Unable to open file: " << fileName <<endl;return false;}
	enum FileType {FASTA, FASTQ, UNDEFINED};
	FileType fileType = UNDEFINED;
	if(fileType == UNDEFINED)
	{
		string text;
		getline (myFile,text);
		if(text[0] == '>')
			fileType = FASTA;
		else if(text[0] == '@')
			fileType = FASTQ;
		else
			cout<< "Unknown input file format."<<endl;
		myFile.seekg(0, ios::beg);
	}

	UINT64 lineOffset;
	if(fileType == FASTA)
		lineOffset = 2;
	else if(fileType == FASTQ)
		lineOffset = 4;


	while(!myFile.eof())
	{
		UINT64 chunkCounter;
		vector<string>* lineData = new vector<string>();
		for(chunkCounter = 0; !myFile.eof() && chunkCounter<Config::queryChunkSize; chunkCounter++ )
		{
			string text = "";
			if(fileType == FASTA) 			// Fasta file
			{

				getline (myFile,text);
				lineData->push_back(text);
				getline (myFile,text,'>');
				text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
				lineData->push_back(text);

			}
			else if(fileType == FASTQ) 					// Fastq file.
			{
				for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
				{
					getline (myFile,text);
					lineData->push_back(text);
				}
			}
		}

		omp_set_dynamic(0);
		omp_set_num_threads(numberOfThreads);
		#pragma omp parallel shared(goodReads, badReads, totalNumberOfReads, thisqueryReadList)
			{
		#pragma omp for schedule(dynamic)
				for(UINT64 i = 0; i < chunkCounter; i++)
				{
					string line0, line1;
					line0 = lineData->at(i*lineOffset);
					line1 = lineData->at(i*lineOffset+1);
					string readname;
					if(line0[0] == '>' || line0[0] == '@')
						readname = line0.substr(1);
					else readname = line0;

					for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
					    *p = toupper(*p);
					if(line1.length() > dataset_minimumoverlaplength && qualityFilter(line1) ) // Test the read is of good quality.
					{

						QueryRead *r1=new QueryRead();
						r1->setName(readname);
						string line1ReverseComplement;
						line1ReverseComplement = r1->reverseComplement(line1);
						if(line1.compare(line1ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
							r1->setSequence(line1,line1ReverseComplement);
						else
							r1->setSequence(line1ReverseComplement,line1);
						r1->setFrequency(1);



					#pragma omp critical(goodreadcount)
						{
						thisqueryReadList->push_back(r1);						// Store the first string in the dataset.
						totalNumberOfReads++;							// Counter of the total number of reads.
						goodReads++;
						}
					}
					else
					{
					#pragma omp critical(badreadcount)
						{
						badReads++;
						}
					}
				}

			}

		lineData->clear();
		delete lineData;

	}

	myFile.close();

	this->numberOfReads = totalNumberOfReads;
    cout << "File name: " << fileName << endl;
	cout << setw(10) << goodReads << " good reads in current dataset."  << endl;
	cout << setw(10) << badReads << " bad reads in current dataset." << endl;
	cout << setw(10) << goodReads + badReads << " total reads in current dataset." << endl;
	cout << setw(10) << numberOfReads << " good reads in all datasets." << endl << endl;


	//sort the reads by its alphabetical order; and remove the duplicate reads so then the readID can be assigned uniquely
	sortReads();
	duplicateFilter();


	return true;
}


bool QueryDataset::buildDataset(const string & QueryFilename)
{
	if(queryReadList==NULL)
		queryReadList = new vector<QueryRead *>();
	cout << "initial " << queryReadList->size() << endl;
	string fileName = QueryFilename;
	cout << "Reading dataset: "  << " from file: " << fileName << endl;
	ifstream myFile;
	myFile.open (fileName.c_str());
	if(myFile == NULL)
		{cout<<"Unable to open file: " << fileName <<endl;return false;}
	UINT64 goodReads = 0, badReads = 0;
	vector<string> line;
	string line1,line0, text, line1ReverseComplement,line2ReverseComplement;
	enum FileType {FASTA, FASTQ, UNDEFINED};
	FileType fileType = UNDEFINED;

	while(!myFile.eof())
	{
		if( (goodReads + badReads ) != 0 && (goodReads + badReads)%1000000 == 0)
			cout<< setw(10) << goodReads + badReads << " reads processed in dataset " << setw(2) << setw(10) << goodReads << " good reads." << setw(10) << badReads << " bad reads." << endl;
		if(fileType == UNDEFINED)
		{
			getline (myFile,text);
			if(text[0] == '>')
				fileType = FASTA;
			else if(text[0] == '@')
				fileType = FASTQ;
			else
				cout<< "Unknown input file format."<<endl;
			myFile.seekg(0, ios::beg);
		}
		line.clear();
		if(fileType == FASTA) 			// Fasta file
		{

			getline (myFile,text);
			line.push_back(text);
			getline (myFile,text,'>');
			line.push_back(text);

			line.at(1).erase(std::remove(line.at(1).begin(), line.at(1).end(), '\n'), line.at(1).end());
			line1 = line.at(1);								// The first string is in the 2nd line.

			line0 = line.at(0);                              //title line is the very first line


		}
		else if(fileType == FASTQ) 					// Fastq file.
		{
			for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
			{
				getline (myFile,text);
				line.push_back(text);
			}
			line1 = line.at(1); 			// The first string is in the 2nd line.
			line0 = line.at(0);             //title line is the very first line
		}

		//remove the > if it appears in the title or name
		string readname="";
		if(line0[0] == '>' || line0[0] == '@')
			readname = line0.substr(1);
		else readname = line0;

		for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
		    *p = toupper(*p);
		if(line1.length() > dataset_minimumoverlaplength && qualityFilter(line1) ) // Test the read is of good quality.
		{

			QueryRead *r1=new QueryRead();
			r1->setName(readname);
			line1ReverseComplement = r1->reverseComplement(line1);
			if(line1.compare(line1ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
				r1->setSequence(line1,line1ReverseComplement);
			else
				r1->setSequence(line1ReverseComplement,line1);
			r1->setFrequency(1);

			UINT64 len = r1->getReadLength();

			if(len > longestReadLength)
				longestReadLength = len;
			if(len < shortestReadLength)
				shortestReadLength = len;
//			cout << queryReadList.size() << endl;

//			cout << r1->getSequence() << endl;
			queryReadList->push_back(r1);						// Store the first string in the dataset.

			numberOfReads++;							// Counter of the total number of reads.
			goodReads++;
		}
		else
			badReads++;
	}

	myFile.close();

    cout << "File name: " << fileName << endl;
	cout << setw(10) << goodReads << " good reads in current dataset."  << endl;
	cout << setw(10) << badReads << " bad reads in current dataset." << endl;
	cout << setw(10) << goodReads + badReads << " total reads in current dataset." << endl;
	cout << setw(10) << numberOfReads << " good reads in all datasets." << endl << endl;


	//sort the reads by its alphabetical order; and remove the duplicate reads so then the readID can be assigned uniquely
	sortReads();
	duplicateFilter();


	return true;
}

bool QueryDataset::buildDatasetNoSortingDuplicateRemoving()
{
	return this->buildDatasetNoSortingDuplicateRemoving(this->dataset_QueryFilename, Config::numberOfThreads);
}

bool QueryDataset::buildDatasetNoSortingDuplicateRemoving(const string & QueryFilename, UINT16 numberOfThreads)
{
	if(queryReadList==NULL)
		queryReadList = new vector<QueryRead *>();
	cout << "initial " << queryReadList->size() << endl;

	UINT64 goodReads = 0, badReads = 0;
	UINT64 totalNumberOfReads = this->numberOfReads;
	vector<QueryRead *>*& thisqueryReadList = this->queryReadList;
	string fileName = QueryFilename;
	cout << "Reading dataset: "  << " from file: " << fileName << endl;
	ifstream myFile;
	myFile.open (fileName.c_str());
	if(myFile == NULL)
		{cout<<"Unable to open file: " << fileName <<endl;return false;}
	enum FileType {FASTA, FASTQ, UNDEFINED};
	FileType fileType = UNDEFINED;
	if(fileType == UNDEFINED)
	{
		string text;
		getline (myFile,text);
		if(text[0] == '>')
			fileType = FASTA;
		else if(text[0] == '@')
			fileType = FASTQ;
		else
			cout<< "Unknown input file format."<<endl;
		myFile.seekg(0, ios::beg);
	}

	UINT64 lineOffset;
	if(fileType == FASTA)
		lineOffset = 2;
	else if(fileType == FASTQ)
		lineOffset = 4;


	while(!myFile.eof())
	{
		UINT64 chunkCounter;
		vector<string>* lineData = new vector<string>();
		for(chunkCounter = 0; !myFile.eof() && chunkCounter<Config::queryChunkSize; chunkCounter++ )
		{
			string text = "";
			if(fileType == FASTA) 			// Fasta file
			{

				getline (myFile,text);
				lineData->push_back(text);
				getline (myFile,text,'>');
				text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
				lineData->push_back(text);

			}
			else if(fileType == FASTQ) 					// Fastq file.
			{
				for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
				{
					getline (myFile,text);
					lineData->push_back(text);
				}
			}
		}

		omp_set_dynamic(0);
		omp_set_num_threads(numberOfThreads);
		#pragma omp parallel shared(goodReads, badReads, totalNumberOfReads, thisqueryReadList)
			{
		#pragma omp for schedule(dynamic)
				for(UINT64 i = 0; i < chunkCounter; i++)
				{
					string line0, line1;
					line0 = lineData->at(i*lineOffset);
					line1 = lineData->at(i*lineOffset+1);
					string readname;
					if(line0[0] == '>' || line0[0] == '@')
						readname = line0.substr(1);
					else readname = line0;

					for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
					    *p = toupper(*p);
					if(line1.length() > dataset_minimumoverlaplength && qualityFilter(line1) ) // Test the read is of good quality.
					{

						QueryRead *r1=new QueryRead();
						r1->setName(readname);
						string line1ReverseComplement;
						line1ReverseComplement = r1->reverseComplement(line1);
						if(line1.compare(line1ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
							r1->setSequence(line1,line1ReverseComplement);
						else
							r1->setSequence(line1ReverseComplement,line1);

						r1->setFrequency(1);



					#pragma omp critical(goodreadcount)
						{
						thisqueryReadList->push_back(r1);						// Store the first string in the dataset.
						totalNumberOfReads++;							// Counter of the total number of reads.
						goodReads++;
						}
					}
					else
					{
					#pragma omp critical(badreadcount)
						{
						badReads++;
						}
					}
				}

			}

		lineData->clear();
		delete lineData;

	}

	myFile.close();

	this->numberOfReads = totalNumberOfReads;
	this->numberOfUniqueReads = this->numberOfReads;

    cout << "File name: " << fileName << endl;
	cout << setw(10) << goodReads << " good reads in current dataset."  << endl;
	cout << setw(10) << badReads << " bad reads in current dataset." << endl;
	cout << setw(10) << goodReads + badReads << " total reads in current dataset." << endl;
	cout << setw(10) << numberOfReads << " good reads in all datasets." << endl << endl;


	//assign IDs to the query reads in dataset
	for(UINT64 i = 0 ; i < this->queryReadList->size(); i++) 		// Assign ID's to the reads.
	{
		if(i < getNumberOfUniqueReads())
			queryReadList->at(i)->setIdentifier(i + 1);
		else
			delete queryReadList->at(i); 					// Free the unused reads.
	}
	queryReadList->resize(numberOfUniqueReads);				//Resize the list.


	return true;
}

bool QueryDataset::buildDatasetFromMatePairFile(const string& QueryFilename, UINT8 orient)
{
	switch(orient)
	{
	case 0:
		return this->buildDatasetFromMatePairFileRR(QueryFilename);
		break;
	case 1:
		return this->buildDatasetFromMatePairFileRF(QueryFilename);
		break;
	case 2:
		return this->buildDatasetFromMatePairFileFR(QueryFilename);
		break;
	case 3:
		return this->buildDatasetFromMatePairFileFF(QueryFilename);
		break;
	default:
		return false;
	}
}

bool QueryDataset::buildDatasetFromMatePairFileFF(const string& QueryFilename)
{
	if(this->pairedEndReadList==NULL)this->pairedEndReadList = new vector<PairedEndRead*>();
	string fileName = QueryFilename;
	cout << "Store paired-end information of dataset: " << " from file: " << fileName << endl;

	if(queryReadList==NULL)
		queryReadList = new vector<QueryRead *>();

	UINT64 goodMatePairs = 0, badMatePairs = 0;
	UINT64 totalNumberOfReads = this->numberOfReads;
	ifstream myFile;
	myFile.open (fileName.c_str());
	if(myFile == NULL)
	{
		cout<<"Unable to open file: "<<fileName;
		return false;
	}
	enum FileType {FASTA, FASTQ, UNDEFINED};
	FileType fileType = UNDEFINED;
	if(fileType == UNDEFINED)
	{
		string text;
		getline (myFile,text);
		if(text[0] == '>')
			fileType = FASTA;
		else if(text[0] == '@')
			fileType = FASTQ;
		else
			cout<< "Unknown input file format."<<endl;
		myFile.seekg(0, ios::beg);
	}

	UINT64 lineOffset;
	if(fileType == FASTA)
		lineOffset = 2;
	else if(fileType == FASTQ)
		lineOffset = 4;


	while(!myFile.eof())
	{
		UINT64 chunkCounter;
		vector<string>* lineData = new vector<string>();
		for(chunkCounter = 0; !myFile.eof() && chunkCounter<Config::queryChunkSize; chunkCounter++ )
		{
			string text = "";
			if(fileType == FASTA)	 				// Fasta file
			{
				getline (myFile,text);
				lineData->push_back(text);
				getline (myFile,text,'>');
				text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
				lineData->push_back(text);

				getline (myFile,text);
				lineData->push_back(text);
				getline (myFile,text,'>');
				text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
				lineData->push_back(text);

			}
			else if(fileType == FASTQ)  			// Fastq file.
			{
				for(UINT64 i = 0; i < 8; i++)		// Read the remaining 7 lines. Total of 8 lines represent two sequence in a fastq file.
				{
					getline (myFile,text);
					lineData->push_back(text);
				}
			}
		}


		omp_set_dynamic(0);
		omp_set_num_threads(Config::numberOfThreads);
		#pragma omp parallel shared(goodMatePairs, badMatePairs, totalNumberOfReads)
			{
		#pragma omp for schedule(dynamic)
				for(UINT64 i = 0; i < chunkCounter; i++)
				{
					string line0, line1, line2, line3;
					line0 = lineData->at(2*i*lineOffset);
					line1 = lineData->at(2*i*lineOffset+1);
					line2 = lineData->at((2*i+1)*lineOffset);
					line3 = lineData->at((2*i+1)*lineOffset+1);
					string readname1="";
					if(line0[0] == '>' || line0[0] == '@')
						readname1 = line0.substr(1);
					else readname1 = line0;

					string readname2="";
					if(line2[0] == '>' || line2[0] == '@')
						readname2 = line2.substr(1);
					else readname2 = line2;

					for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
					    *p = toupper(*p);
					for (std::string::iterator p = line3.begin(); line3.end() != p; ++p) // Change the case
					    *p = toupper(*p);
					if(line1.length() > dataset_leftminimumoverlaplength && qualityFilter(line1)
					&& line3.length() > dataset_rightminimumoverlaplength && qualityFilter(line3)) // Test the read is of good quality.
					{

//						UINT8 r1forward, r2forward;
						QueryRead *r1=new QueryRead();
						r1->setName(readname1);
						string line1ReverseComplement;
						line1ReverseComplement = r1->reverseComplement(line1);
//						if(line1.compare(line1ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
//							{r1->setSequence(line1,line1ReverseComplement);r1forward = 1;}
//						else
//							{r1->setSequence(line1ReverseComplement,line1);r1forward = 0;}

						r1->setSequence(line1,line1ReverseComplement);//r1forward = 1;
						r1->setFrequency(1);


						QueryRead *r2=new QueryRead();
						r2->setName(readname2);
						string line2ReverseComplement;
						line2ReverseComplement = r2->reverseComplement(line3);
//						if(line3.compare(line2ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
//							{r2->setSequence(line3,line2ReverseComplement);r2forward = 1;}
//						else
//							{r2->setSequence(line2ReverseComplement,line3);r2forward = 0;}
						r2->setSequence(line3,line2ReverseComplement);//r2forward = 1;
						r2->setFrequency(1);

//						UINT8 orientation = r1forward*2+r2forward;
//						r1->addRightMatePair(r2,orientation);
						r1->addRightMatePair(r2,3);

					#pragma omp critical(goodreadcount)
						{
						r1->setIdentifier(totalNumberOfReads+1);
						r2->setIdentifier(totalNumberOfReads+2);
						queryReadList->push_back(r1);						// Store the first string in the dataset.

						queryReadList->push_back(r2);						// Store the first string in the dataset.
						this->pairedEndReadList->push_back(r1->getPairEndRead());
						totalNumberOfReads += 2;	// Counter of the total number of reads.
						goodMatePairs += 2;
						}
					}
					else
					{
					#pragma omp critical(badreadcount)
						{
							badMatePairs +=2;
						}
					}
				}

			}

		lineData->clear();
		delete lineData;
	}
	myFile.close();

	this->numberOfReads = totalNumberOfReads;

	cout << "File name: " << fileName << endl;
	cout << setw(10) << goodMatePairs << " reads in " << setw(10) << goodMatePairs/2 << " mate-pairs are good." << endl;
	cout << setw(10) << badMatePairs << " reads in " << setw(10) << badMatePairs/2 << " mate-pairs are discarded." << endl << endl;

	this->numberOfUniqueReads = this->numberOfReads;
	return true;
}

bool QueryDataset::buildDatasetFromMatePairFileFR(const string& QueryFilename)
{
	if(this->pairedEndReadList==NULL)this->pairedEndReadList = new vector<PairedEndRead*>();
	string fileName = QueryFilename;
	cout << "Store paired-end information of dataset: " << " from file: " << fileName << endl;

	if(queryReadList==NULL)
		queryReadList = new vector<QueryRead *>();

	UINT64 goodMatePairs = 0, badMatePairs = 0;
	UINT64 totalNumberOfReads = this->numberOfReads;
	ifstream myFile;
	myFile.open (fileName.c_str());
	if(myFile == NULL)
	{
		cout<<"Unable to open file: "<<fileName;
		return false;
	}
	enum FileType {FASTA, FASTQ, UNDEFINED};
	FileType fileType = UNDEFINED;
	if(fileType == UNDEFINED)
	{
		string text;
		getline (myFile,text);
		if(text[0] == '>')
			fileType = FASTA;
		else if(text[0] == '@')
			fileType = FASTQ;
		else
			cout<< "Unknown input file format."<<endl;
		myFile.seekg(0, ios::beg);
	}

	UINT64 lineOffset;
	if(fileType == FASTA)
		lineOffset = 2;
	else if(fileType == FASTQ)
		lineOffset = 4;


	while(!myFile.eof())
	{
		UINT64 chunkCounter;
		vector<string>* lineData = new vector<string>();
		for(chunkCounter = 0; !myFile.eof() && chunkCounter<Config::queryChunkSize; chunkCounter++ )
		{
			string text = "";
			if(fileType == FASTA)	 				// Fasta file
			{
				getline (myFile,text);
				lineData->push_back(text);
				getline (myFile,text,'>');
				text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
				lineData->push_back(text);

				getline (myFile,text);
				lineData->push_back(text);
				getline (myFile,text,'>');
				text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
				lineData->push_back(text);

			}
			else if(fileType == FASTQ)  			// Fastq file.
			{
				for(UINT64 i = 0; i < 8; i++)		// Read the remaining 7 lines. Total of 8 lines represent two sequence in a fastq file.
				{
					getline (myFile,text);
					lineData->push_back(text);
				}
			}
		}


		omp_set_dynamic(0);
		omp_set_num_threads(Config::numberOfThreads);
		#pragma omp parallel shared(goodMatePairs, badMatePairs, totalNumberOfReads)
			{
		#pragma omp for schedule(dynamic)
				for(UINT64 i = 0; i < chunkCounter; i++)
				{
					string line0, line1, line2, line3;
					line0 = lineData->at(2*i*lineOffset);
					line1 = lineData->at(2*i*lineOffset+1);
					line2 = lineData->at((2*i+1)*lineOffset);
					line3 = lineData->at((2*i+1)*lineOffset+1);
					string readname1="";
					if(line0[0] == '>' || line0[0] == '@')
						readname1 = line0.substr(1);
					else readname1 = line0;

					string readname2="";
					if(line2[0] == '>' || line2[0] == '@')
						readname2 = line2.substr(1);
					else readname2 = line2;

					for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
					    *p = toupper(*p);
					for (std::string::iterator p = line3.begin(); line3.end() != p; ++p) // Change the case
					    *p = toupper(*p);
					if(line1.length() > dataset_leftminimumoverlaplength && qualityFilter(line1)
					&& line3.length() > dataset_rightminimumoverlaplength && qualityFilter(line3)) // Test the read is of good quality.
					{

//						UINT8 r1forward, r2forward;
						QueryRead *r1=new QueryRead();
						r1->setName(readname1);
						string line1ReverseComplement;
						line1ReverseComplement = r1->reverseComplement(line1);
//						if(line1.compare(line1ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
//							{r1->setSequence(line1,line1ReverseComplement);r1forward = 1;}
//						else
//							{r1->setSequence(line1ReverseComplement,line1);r1forward = 0;}

						r1->setSequence(line1,line1ReverseComplement);//r1forward = 1;
						r1->setFrequency(1);


						QueryRead *r2=new QueryRead();
						r2->setName(readname2);
						string line2ReverseComplement;
						line2ReverseComplement = r2->reverseComplement(line3);
//						if(line3.compare(line2ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
//							{r2->setSequence(line3,line2ReverseComplement);r2forward = 1;}
//						else
//							{r2->setSequence(line2ReverseComplement,line3);r2forward = 0;}
						r2->setSequence(line2ReverseComplement,line3);//r2forward = 0;
						r2->setFrequency(1);

//						UINT8 orientation = r1forward*2+r2forward;
//						r1->addRightMatePair(r2,orientation);
						r1->addRightMatePair(r2,2);

					#pragma omp critical(goodreadcount)
						{
						r1->setIdentifier(totalNumberOfReads+1);
						r2->setIdentifier(totalNumberOfReads+2);
						queryReadList->push_back(r1);						// Store the first string in the dataset.

						queryReadList->push_back(r2);						// Store the first string in the dataset.
						this->pairedEndReadList->push_back(r1->getPairEndRead());
						totalNumberOfReads += 2;	// Counter of the total number of reads.
						goodMatePairs += 2;
						}
					}
					else
					{
					#pragma omp critical(badreadcount)
						{
							badMatePairs +=2;
						}
					}
				}

			}

		lineData->clear();
		delete lineData;
	}
	myFile.close();

	this->numberOfReads = totalNumberOfReads;

	cout << "File name: " << fileName << endl;
	cout << setw(10) << goodMatePairs << " reads in " << setw(10) << goodMatePairs/2 << " mate-pairs are good." << endl;
	cout << setw(10) << badMatePairs << " reads in " << setw(10) << badMatePairs/2 << " mate-pairs are discarded." << endl << endl;

	this->numberOfUniqueReads = this->numberOfReads;
	return true;
}

bool QueryDataset::buildDatasetFromMatePairFileRF(const string& QueryFilename)
{
	if(this->pairedEndReadList==NULL)this->pairedEndReadList = new vector<PairedEndRead*>();
	string fileName = QueryFilename;
	cout << "Store paired-end information of dataset: " << " from file: " << fileName << endl;

	if(queryReadList==NULL)
		queryReadList = new vector<QueryRead *>();

	UINT64 goodMatePairs = 0, badMatePairs = 0;
	UINT64 totalNumberOfReads = this->numberOfReads;
	ifstream myFile;
	myFile.open (fileName.c_str());
	if(myFile == NULL)
	{
		cout<<"Unable to open file: "<<fileName;
		return false;
	}
	enum FileType {FASTA, FASTQ, UNDEFINED};
	FileType fileType = UNDEFINED;
	if(fileType == UNDEFINED)
	{
		string text;
		getline (myFile,text);
		if(text[0] == '>')
			fileType = FASTA;
		else if(text[0] == '@')
			fileType = FASTQ;
		else
			cout<< "Unknown input file format."<<endl;
		myFile.seekg(0, ios::beg);
	}

	UINT64 lineOffset;
	if(fileType == FASTA)
		lineOffset = 2;
	else if(fileType == FASTQ)
		lineOffset = 4;


	while(!myFile.eof())
	{
		UINT64 chunkCounter;
		vector<string>* lineData = new vector<string>();
		for(chunkCounter = 0; !myFile.eof() && chunkCounter<Config::queryChunkSize; chunkCounter++ )
		{
			string text = "";
			if(fileType == FASTA)	 				// Fasta file
			{
				getline (myFile,text);
				lineData->push_back(text);
				getline (myFile,text,'>');
				text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
				lineData->push_back(text);

				getline (myFile,text);
				lineData->push_back(text);
				getline (myFile,text,'>');
				text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
				lineData->push_back(text);

			}
			else if(fileType == FASTQ)  			// Fastq file.
			{
				for(UINT64 i = 0; i < 8; i++)		// Read the remaining 7 lines. Total of 8 lines represent two sequence in a fastq file.
				{
					getline (myFile,text);
					lineData->push_back(text);
				}
			}
		}


		omp_set_dynamic(0);
		omp_set_num_threads(Config::numberOfThreads);
		#pragma omp parallel shared(goodMatePairs, badMatePairs, totalNumberOfReads)
			{
		#pragma omp for schedule(dynamic)
				for(UINT64 i = 0; i < chunkCounter; i++)
				{
					string line0, line1, line2, line3;
					line0 = lineData->at(2*i*lineOffset);
					line1 = lineData->at(2*i*lineOffset+1);
					line2 = lineData->at((2*i+1)*lineOffset);
					line3 = lineData->at((2*i+1)*lineOffset+1);
					string readname1="";
					if(line0[0] == '>' || line0[0] == '@')
						readname1 = line0.substr(1);
					else readname1 = line0;

					string readname2="";
					if(line2[0] == '>' || line2[0] == '@')
						readname2 = line2.substr(1);
					else readname2 = line2;

					for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
					    *p = toupper(*p);
					for (std::string::iterator p = line3.begin(); line3.end() != p; ++p) // Change the case
					    *p = toupper(*p);
					if(line1.length() > dataset_leftminimumoverlaplength && qualityFilter(line1)
					&& line3.length() > dataset_rightminimumoverlaplength && qualityFilter(line3)) // Test the read is of good quality.
					{

//						UINT8 r1forward, r2forward;
						QueryRead *r1=new QueryRead();
						r1->setName(readname1);
						string line1ReverseComplement;
						line1ReverseComplement = r1->reverseComplement(line1);
//						if(line1.compare(line1ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
//							{r1->setSequence(line1,line1ReverseComplement);r1forward = 1;}
//						else
//							{r1->setSequence(line1ReverseComplement,line1);r1forward = 0;}

						r1->setSequence(line1ReverseComplement,line1);//r1forward = 0;
						r1->setFrequency(1);


						QueryRead *r2=new QueryRead();
						r2->setName(readname2);
						string line2ReverseComplement;
						line2ReverseComplement = r2->reverseComplement(line3);
//						if(line3.compare(line2ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
//							{r2->setSequence(line3,line2ReverseComplement);r2forward = 1;}
//						else
//							{r2->setSequence(line2ReverseComplement,line3);r2forward = 0;}
						r2->setSequence(line3,line2ReverseComplement);//r2forward = 1;
						r2->setFrequency(1);

//						UINT8 orientation = r1forward*2+r2forward;
//						r1->addRightMatePair(r2,orientation);
						r1->addRightMatePair(r2,1);

					#pragma omp critical(goodreadcount)
						{
						r1->setIdentifier(totalNumberOfReads+1);
						r2->setIdentifier(totalNumberOfReads+2);
						queryReadList->push_back(r1);						// Store the first string in the dataset.

						queryReadList->push_back(r2);						// Store the first string in the dataset.
						this->pairedEndReadList->push_back(r1->getPairEndRead());
						totalNumberOfReads += 2;	// Counter of the total number of reads.
						goodMatePairs += 2;
						}
					}
					else
					{
					#pragma omp critical(badreadcount)
						{
							badMatePairs +=2;
						}
					}
				}

			}

		lineData->clear();
		delete lineData;
	}
	myFile.close();

	this->numberOfReads = totalNumberOfReads;

	cout << "File name: " << fileName << endl;
	cout << setw(10) << goodMatePairs << " reads in " << setw(10) << goodMatePairs/2 << " mate-pairs are good." << endl;
	cout << setw(10) << badMatePairs << " reads in " << setw(10) << badMatePairs/2 << " mate-pairs are discarded." << endl << endl;

	this->numberOfUniqueReads = this->numberOfReads;
	return true;
}

bool QueryDataset::buildDatasetFromMatePairFileRR(const string& QueryFilename)
{
	if(this->pairedEndReadList==NULL)this->pairedEndReadList = new vector<PairedEndRead*>();
	string fileName = QueryFilename;
	cout << "Store paired-end information of dataset: " << " from file: " << fileName << endl;

	if(queryReadList==NULL)
		queryReadList = new vector<QueryRead *>();

	UINT64 goodMatePairs = 0, badMatePairs = 0;
	UINT64 totalNumberOfReads = this->numberOfReads;
	ifstream myFile;
	myFile.open (fileName.c_str());
	if(myFile == NULL)
	{
		cout<<"Unable to open file: "<<fileName;
		return false;
	}
	enum FileType {FASTA, FASTQ, UNDEFINED};
	FileType fileType = UNDEFINED;
	if(fileType == UNDEFINED)
	{
		string text;
		getline (myFile,text);
		if(text[0] == '>')
			fileType = FASTA;
		else if(text[0] == '@')
			fileType = FASTQ;
		else
			cout<< "Unknown input file format."<<endl;
		myFile.seekg(0, ios::beg);
	}

	UINT64 lineOffset;
	if(fileType == FASTA)
		lineOffset = 2;
	else if(fileType == FASTQ)
		lineOffset = 4;


	while(!myFile.eof())
	{
		UINT64 chunkCounter;
		vector<string>* lineData = new vector<string>();
		for(chunkCounter = 0; !myFile.eof() && chunkCounter<Config::queryChunkSize; chunkCounter++ )
		{
			string text = "";
			if(fileType == FASTA)	 				// Fasta file
			{
				getline (myFile,text);
				lineData->push_back(text);
				getline (myFile,text,'>');
				text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
				lineData->push_back(text);

				getline (myFile,text);
				lineData->push_back(text);
				getline (myFile,text,'>');
				text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
				lineData->push_back(text);

			}
			else if(fileType == FASTQ)  			// Fastq file.
			{
				for(UINT64 i = 0; i < 8; i++)		// Read the remaining 7 lines. Total of 8 lines represent two sequence in a fastq file.
				{
					getline (myFile,text);
					lineData->push_back(text);
				}
			}
		}


		omp_set_dynamic(0);
		omp_set_num_threads(Config::numberOfThreads);
		#pragma omp parallel shared(goodMatePairs, badMatePairs, totalNumberOfReads)
			{
		#pragma omp for schedule(dynamic)
				for(UINT64 i = 0; i < chunkCounter; i++)
				{
					string line0, line1, line2, line3;
					line0 = lineData->at(2*i*lineOffset);
					line1 = lineData->at(2*i*lineOffset+1);
					line2 = lineData->at((2*i+1)*lineOffset);
					line3 = lineData->at((2*i+1)*lineOffset+1);
					string readname1="";
					if(line0[0] == '>' || line0[0] == '@')
						readname1 = line0.substr(1);
					else readname1 = line0;

					string readname2="";
					if(line2[0] == '>' || line2[0] == '@')
						readname2 = line2.substr(1);
					else readname2 = line2;

					for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
					    *p = toupper(*p);
					for (std::string::iterator p = line3.begin(); line3.end() != p; ++p) // Change the case
					    *p = toupper(*p);
					if(line1.length() > dataset_leftminimumoverlaplength && qualityFilter(line1)
					&& line3.length() > dataset_rightminimumoverlaplength && qualityFilter(line3)) // Test the read is of good quality.
					{

//						UINT8 r1forward, r2forward;
						QueryRead *r1=new QueryRead();
						r1->setName(readname1);
						string line1ReverseComplement;
						line1ReverseComplement = r1->reverseComplement(line1);
//						if(line1.compare(line1ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
//							{r1->setSequence(line1,line1ReverseComplement);r1forward = 1;}
//						else
//							{r1->setSequence(line1ReverseComplement,line1);r1forward = 0;}

						r1->setSequence(line1ReverseComplement,line1);//r1forward = 0;
						r1->setFrequency(1);


						QueryRead *r2=new QueryRead();
						r2->setName(readname2);
						string line2ReverseComplement;
						line2ReverseComplement = r2->reverseComplement(line3);
//						if(line3.compare(line2ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
//							{r2->setSequence(line3,line2ReverseComplement);r2forward = 1;}
//						else
//							{r2->setSequence(line2ReverseComplement,line3);r2forward = 0;}
						r2->setSequence(line2ReverseComplement,line3);//r2forward = 0;
						r2->setFrequency(1);

//						UINT8 orientation = r1forward*2+r2forward;
//						r1->addRightMatePair(r2,orientation);
						r1->addRightMatePair(r2,0);

					#pragma omp critical(goodreadcount)
						{
						r1->setIdentifier(totalNumberOfReads+1);
						r2->setIdentifier(totalNumberOfReads+2);
						queryReadList->push_back(r1);						// Store the first string in the dataset.

						queryReadList->push_back(r2);						// Store the first string in the dataset.
						this->pairedEndReadList->push_back(r1->getPairEndRead());
						totalNumberOfReads += 2;	// Counter of the total number of reads.
						goodMatePairs += 2;
						}
					}
					else
					{
					#pragma omp critical(badreadcount)
						{
							badMatePairs +=2;
						}
					}
				}

			}

		lineData->clear();
		delete lineData;
	}
	myFile.close();

	this->numberOfReads = totalNumberOfReads;

	cout << "File name: " << fileName << endl;
	cout << setw(10) << goodMatePairs << " reads in " << setw(10) << goodMatePairs/2 << " mate-pairs are good." << endl;
	cout << setw(10) << badMatePairs << " reads in " << setw(10) << badMatePairs/2 << " mate-pairs are discarded." << endl << endl;

	this->numberOfUniqueReads = this->numberOfReads;
	return true;
}

vector<PairedEndRead *>* QueryDataset::getPairedEndReadList()
{
	return this->pairedEndReadList;
}
