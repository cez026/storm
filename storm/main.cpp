#include "Config.h"
#include "QueryDataset.h"
#include "HashTableMethod.h"
#include "SingleKeyHashTable.h"
#include "DoubleKeyHashTable.h"
#include "OmegaHashTable.h"
#include "OmegaGraphConstructor.h"
#include "OverlapGraphConstructor.h"
#include "QueryDatasetFilter.h"
//#include "PairedEndReadsMerger.h"
#include "PairEndReadMerger.h"
#include "GeneralQueryDatasetFilter.h"
//#include "ErrorCorrector.h"


int main(int argc, char **argv)
{


	CLOCKSTART;


	// parse command line options:
	if(!Config::setConfig(argc, argv)){
		cout << "Please follow the above help information." << endl;
		return false;
	}


//--------generate data sets-----------

	QueryDataset* queryDataset = new QueryDataset(Config::getQueryDatasetFilename());
	if(Config::getOperation()=="ConstructOverlapGraph")
	{
		{
		CLOCKSTART;
		MEMORYSTART;

		if (!queryDataset->buildDatasetNoSortingDuplicateRemoving()){
			cout << "Error: cannot build query dataset" << endl;
			return false;
		}
		else if(queryDataset->getNumberOfUniqueReads()==0){
			cout << "Unfortunately all the query data are filtered out. Further processing is aborted. Please double check with the data or command settings." <<endl;
			return false;
		}
		MEMORYSTOP;
		CLOCKSTOP;

		}
	}
	else if(Config::getOperation()=="ConstructGraph" || Config::getOperation()=="RemoveContainedReads")
	{

		{
		CLOCKSTART;
		MEMORYSTART;

		if (!queryDataset->buildDataset()){
			cout << "Error: cannot build query dataset" << endl;
			return false;
		}
		else if(queryDataset->getNumberOfUniqueReads()==0){
			cout << "Unfortunately all the query data are filtered out. Further processing is aborted. Please double check with the data or command settings." <<endl;
			return false;
		}
		MEMORYSTOP;
		CLOCKSTOP;

		}
	}
	else if(Config::getOperation()=="MergePairedEndReads")
	{
		//set the minimumOverlaplength to the same as left, since only the left end reads will be inserted into hashtable.
		Config::minimumOverlapLength = Config::left_minimumOverlapLength;

		{
		CLOCKSTART;
		MEMORYSTART;

		if (!queryDataset->buildDatasetFromMatePairFile(Config::getQueryDatasetFilename(),Config::pairedEndReadsInputMode)){
			cout << "Error: cannot build query dataset" << endl;
			return false;
		}
		else if(queryDataset->getNumberOfReads()==0){
			cout << "Unfortunately all the query data are filtered out. Further processing is aborted. Please double check with the data or command settings." <<endl;
			return false;
		}
		MEMORYSTOP;
		CLOCKSTOP;

		}
	}
	else
	{
		cout<<"wrong operation code"<<endl;
		return false;
	}

//-------------------------------------------------------------------------------
//---------Construct Overlap Graph (Can deal with transitive edge removal)-------
//--------single key hash table method------
if(Config::getOperation()=="ConstructOverlapGraph" && Config::getHashTableType() == "single")
{
	HashTableMethod* singleKeyHashTable = new SingleKeyHashTable(queryDataset);
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------singlekeyhashtable------" << endl;
	singleKeyHashTable->createHashTables();
//	singleKeyHashTable->insertQueryDataset_rightAlignOnly(queryDataset);
	singleKeyHashTable->insertQueryDataset(queryDataset);
	MEMORYSTOP;
	CLOCKSTOP;

		{
		CLOCKSTART;
		MEMORYSTART;
		cout << "-------overlapgraph constructing------" << endl;
		OverlapGraphConstructor * overlapgraphConstructor = new OverlapGraphConstructor(singleKeyHashTable);
		overlapgraphConstructor->start(Config::withTransitiveEdgeReduction);
		MEMORYSTOP;
		CLOCKSTOP;
		delete overlapgraphConstructor;
		}
	}
	delete singleKeyHashTable;
}
//---------double key hash table method-------
if(Config::getOperation()=="ConstructOverlapGraph" && Config::getHashTableType() == "double")
{
	HashTableMethod* doubleKeyHashTable = new DoubleKeyHashTable(queryDataset);
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------doublekeyhashtable------" << endl;
	doubleKeyHashTable->createHashTables();
//	doubleKeyHashTable->insertQueryDataset_rightAlignOnly(queryDataset);
	doubleKeyHashTable->insertQueryDataset(queryDataset);
	MEMORYSTOP;
	CLOCKSTOP;
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------overlapgraph constructing------" << endl;
	OverlapGraphConstructor * overlapgraphConstructor = new OverlapGraphConstructor(doubleKeyHashTable);
	overlapgraphConstructor->start(Config::withTransitiveEdgeReduction);
	MEMORYSTOP;
	CLOCKSTOP;
	delete overlapgraphConstructor;
	}
	}
	delete doubleKeyHashTable;
}

//---------omega hash table method------
if(Config::getOperation()=="ConstructOverlapGraph" && Config::getHashTableType() == "omega")
{

	OmegaHashTable *omegaHashTable=new OmegaHashTable();
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------omegahashtable------" << endl;
//	omegaHashTable->insertDataset_rightAlignOnly(queryDataset, Config::minimumOverlapLength);
	omegaHashTable->insertDataset(queryDataset, Config::minimumOverlapLength);
	MEMORYSTOP;
	CLOCKSTOP;
	}
	OmegaGraphConstructor * omegaGraph = new OmegaGraphConstructor(omegaHashTable);
	omegaGraph->start(Config::withTransitiveEdgeReduction);
	delete omegaGraph;
	delete omegaHashTable;
}



//-----------------------------------------------------
//---------Construct Overlap Graph (Old Approach)-------
//--------single key hash table method------
if(Config::getOperation()=="ConstructGraph" && Config::getHashTableType() == "single")
{
	HashTableMethod* singleKeyHashTable = new SingleKeyHashTable(queryDataset);
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------singlekeyhashtable------" << endl;
	singleKeyHashTable->createHashTables();
	singleKeyHashTable->insertQueryDataset(queryDataset);
	MEMORYSTOP;
	CLOCKSTOP;

		{
		CLOCKSTART;
		MEMORYSTART;
		cout << "-------overlapgraph constructing------" << endl;
		OverlapGraphConstructor * overlapgraphConstructor = new OverlapGraphConstructor(singleKeyHashTable);
		overlapgraphConstructor->start();
		MEMORYSTOP;
		CLOCKSTOP;
		delete overlapgraphConstructor;
		}
	}
	delete singleKeyHashTable;
}
//---------double key hash table method-------
if(Config::getOperation()=="ConstructGraph" && Config::getHashTableType() == "double")
{
	HashTableMethod* doubleKeyHashTable = new DoubleKeyHashTable(queryDataset);
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------doublekeyhashtable------" << endl;
	doubleKeyHashTable->createHashTables();
	doubleKeyHashTable->insertQueryDataset(queryDataset);
	MEMORYSTOP;
	CLOCKSTOP;
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------overlapgraph constructing------" << endl;
	OverlapGraphConstructor * overlapgraphConstructor = new OverlapGraphConstructor(doubleKeyHashTable);
	overlapgraphConstructor->start();
	MEMORYSTOP;
	CLOCKSTOP;
	delete overlapgraphConstructor;
	}
	}
	delete doubleKeyHashTable;
}

//---------omega hash table method------
if(Config::getOperation()=="ConstructGraph" && Config::getHashTableType() == "omega")
{
	OmegaHashTable *omegaHashTable=new OmegaHashTable();
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------omegahashtable------" << endl;
	omegaHashTable->insertDataset(queryDataset, Config::minimumOverlapLength);
	MEMORYSTOP;
	CLOCKSTOP;
	}
	OmegaGraphConstructor * omegaGraph = new OmegaGraphConstructor(omegaHashTable);
	omegaGraph->start();
	delete omegaGraph;
	delete omegaHashTable;
}

//--------------------------------------------------
//--------Merge Paired End reads--------------------
//--------Merge Paired End reads using single key hash table method------
if(Config::getOperation()=="MergePairedEndReads" && Config::getHashTableType() == "single")
{
	HashTableMethod* singleKeyHashTable = new SingleKeyHashTable(queryDataset);
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------singlekeyhashtable------" << endl;
	singleKeyHashTable->createHashTables();
	singleKeyHashTable->insertQueryDataset_leftReadForMerge(queryDataset);
	MEMORYSTOP;
	CLOCKSTOP;

		{
		CLOCKSTART;
		MEMORYSTART;
		cout << "-------merging paired end reads------" << endl;
//		PairedEndReadsMerger * readmerger = new PairedEndReadsMerger(singleKeyHashTable);
		PairEndReadMerger * readmerger = new PairEndReadMerger(singleKeyHashTable);
		readmerger->start();
		MEMORYSTOP;
		CLOCKSTOP;
		delete readmerger;
		}
	}
	delete singleKeyHashTable;
}

//--------Merge Paired End reads using Double key hash table method------
if(Config::getOperation()=="MergePairedEndReads" && Config::getHashTableType() == "double")
{
	HashTableMethod* doubleKeyHashTable = new DoubleKeyHashTable(queryDataset);
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------doublekeyhashtable------" << endl;
	doubleKeyHashTable->createHashTables();
	doubleKeyHashTable->insertQueryDataset_leftReadForMerge(queryDataset);
	MEMORYSTOP;
	CLOCKSTOP;

		{
		CLOCKSTART;
		MEMORYSTART;
		cout << "-------merging paired end reads------" << endl;
//		PairedEndReadsMerger * readmerger = new PairedEndReadsMerger(doubleKeyHashTable);
		PairEndReadMerger * readmerger = new PairEndReadMerger(doubleKeyHashTable);
		readmerger->start();
		MEMORYSTOP;
		CLOCKSTOP;
		delete readmerger;
		}
	}
	delete doubleKeyHashTable;
}

//----------------------------------------------------------------------------
//---------Remove Contained Reads---------------------------------------------
//---------using omega hash table to filter the identical and contained reads---
if(Config::getOperation()=="RemoveContainedReads" && Config::getHashTableType() == "omega")
{
	OmegaHashTable *omegaHashTable=new OmegaHashTable();
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------remove contained reads omegahashtable------" << endl;
	omegaHashTable->insertDataset_half(queryDataset, Config::minimumOverlapLength);
	MEMORYSTOP;
	CLOCKSTOP;
	}
	QueryDatasetFilter * filter = new QueryDatasetFilter(omegaHashTable);
	filter->start();
	delete filter;
	delete omegaHashTable;
}
//---------using single hash table to filter the identical and contained reads---
if(Config::getOperation()=="RemoveContainedReads" && Config::getHashTableType() == "single")
{
	HashTableMethod* singleKeyHashTable = new SingleKeyHashTable(queryDataset);
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------singlekeyhashtable------" << endl;
	singleKeyHashTable->createHashTables();
	singleKeyHashTable->insertQueryDataset(queryDataset);
	MEMORYSTOP;
	CLOCKSTOP;

		{
		CLOCKSTART;
		MEMORYSTART;
		cout << "-------remove the identical and contained reads------" << endl;
		GeneralQueryDatasetFilter * filter = new GeneralQueryDatasetFilter(singleKeyHashTable);
		filter->start();
		MEMORYSTOP;
		CLOCKSTOP;
		delete filter;

		}
	}
	delete singleKeyHashTable;
}
//---------using double hash table to filter the identical and contained reads---
if(Config::getOperation()=="RemoveContainedReads" && Config::getHashTableType() == "double")
{
	HashTableMethod* doubleKeyHashTable = new DoubleKeyHashTable(queryDataset);
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------doublekeyhashtable------" << endl;
	doubleKeyHashTable->createHashTables();
	doubleKeyHashTable->insertQueryDataset(queryDataset);
	MEMORYSTOP;
	CLOCKSTOP;

		{
		CLOCKSTART;
		MEMORYSTART;
		cout << "-------remove the identical and contained reads------" << endl;
		GeneralQueryDatasetFilter * filter = new GeneralQueryDatasetFilter(doubleKeyHashTable);
		filter->start();
		MEMORYSTOP;
		CLOCKSTOP;
		delete filter;

		}
	}
	delete doubleKeyHashTable;
}

	delete queryDataset;


	CLOCKSTOP;


	return true;
}



