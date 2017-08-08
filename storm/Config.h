 /*
 * Config.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef CONFIG_H_
#define CONFIG_H_

// Common headers:

//multi-thread library OPENMP
#include <omp.h>

// C headers:
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

// C++ headers:
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <utility>
#include <limits>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <map>
#include <exception>


using namespace std;

    //define CLOCKSTART clock_t begin = clock(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
    //define CLOCKSTOP clock_t end = clock(); cout << "Function " << __FUNCTION__ << "() finished in " << double(end - begin) / CLOCKS_PER_SEC<< " Seconds." << endl << endl;
    #define CLOCKSTART double begin = omp_get_wtime(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
    #define CLOCKSTOP double end = omp_get_wtime(); cout << "Function " << __FUNCTION__ << "() finished in " << double(end - begin) << " Seconds." << endl << endl;

#define CURTIME time_t curr=time(0);cout << "current time is: " << ctime(&curr) <<endl;
#define CURMEM INT64 mem = checkMemoryUsage();cout<<"Memory used: " << mem <<  " MB."<< endl;


// To keep time information of functions.
#define MEMORYSTART INT64 mem_start = checkMemoryUsage(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
#define MEMORYSTOP INT64 mem_end = checkMemoryUsage(); cout << "Function " << __FUNCTION__ << "() finished . " << "Memory used: " << mem_end << " - " <<  mem_start << " = "<< mem_end - mem_start << " MB."<< endl;
// Get the memory usage with a Linux kernel.
inline unsigned int checkMemoryUsage()
{
    // get KB memory into count
    unsigned int count=0;

   #if defined(__linux__)
    ifstream f("/proc/self/status"); // read the linux file
    while(!f.eof()){
        string key;
        f>>key;
        if(key=="VmData:"){     // size of data
            f>>count;
        break;
        }

    }
    f.close();
   #endif

    // return MBs memory (size of data)
    return (count/1024);
};

#define MAX_VECTOR 500000
typedef unsigned char UINT8;
typedef unsigned short UINT16;
typedef short INT16;
typedef unsigned long UINT32;
typedef long INT32;
typedef unsigned long long UINT64;
typedef long long INT64;

class Config {

public:
	static string operationCode;//"RemoveContainedReads","ConstructGraph","MergePairedEndReads","CorrectErrors"
	static string queryFilename;
	static vector<string> subjectFilenameList;
	static string outputfilename;
	static string hashtabletype;

	static UINT16 minimumOverlapLength; //it's left side minimumoverlaplength when used in the paired end merging process
	static UINT16 left_minimumOverlapLength; // equal to the minimumoverlaplength
	static UINT16 right_minimumOverlapLength;
	static UINT16 overall_minimumOverlapLength;
	static int minimumdepth;
	static double recallpercent;
	static UINT8 pairedEndReadsInputMode;
	static UINT16 hashKeyLength;
	static UINT16 hashKeyLength_left;
	static UINT16 hashKeyLength_right;
	static UINT64 streamChunkSize;
	static UINT64 queryChunkSize; //used for parallel I/O
	static bool isSingleKey;

	static bool isSingleEnd; //subject is always treated as single end, while the query is only paired end when using "MergePairedEndReads"
	static bool isFilter;
	static bool useID;
	static UINT64 startID;
	static bool isPrintIDmap;
	static string IDmapfilename;
	static bool withTransitiveEdgeReduction;

	static bool perfectMatch;
	static UINT16 maxMismatch; //only valid when perfectMatch is set to false
	static UINT16 maxIndel; //only valid when perfectMatch is set to false
	static bool useErrorRate; //defaultly false;
	static UINT16 maxErrorRate; //only valid when useErrorRate is true
	static UINT16 numberOfThreads;

	static void printHelp();

	static inline unsigned int checkMemoryUsage();

	Config();
	~Config();
	static bool setConfig(int argc, char **argv);
	static string getOperation();
	static vector<string> getSubjectDatasetFilenames(); // subject is dataset B streamed from hard drive
	static string getQueryDatasetFilename();	// query is database A loaded to memory
	static UINT16 getminimumOverlapLength();
	static UINT16 getleftminimumOverlapLength();
	static UINT16 getrightminimumOverlapLength();
	static UINT16 getHashKeyLength();
	static UINT16 getHashLeftKeyLength();
	static UINT16 getHashRightKeyLength();
	static UINT64 getStreamChunkSize();

	static bool isSingleKeyHashTable();
	static bool allowMismatch();
	static bool isQuerySingleEndReads();
	static string getHashTableType();

	static UINT16 getMaxMismatch();
	static UINT16 getMaxIndel();
	static bool isFilterOrAlign();



};





#endif /* CONFIG_H_ */
