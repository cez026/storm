/*
 * Config.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "Config.h"

Config::Config() {
	// TODO Auto-generated constructor stub

}

Config::~Config() {
	// TODO Auto-generated destructor stub
}

//Here is to create and initialize the static members
string Config::operationCode="";//"ConstructGraph","MergePairedEndReads","CorrectErrors"
string Config::queryFilename = "";
string Config::hashtabletype = "";

vector<string> Config::subjectFilenameList;
string Config::outputfilename = "out.txt";


UINT16 Config::minimumOverlapLength=40;
UINT16 Config::left_minimumOverlapLength=20;
UINT16 Config::right_minimumOverlapLength=20;
UINT16 Config::overall_minimumOverlapLength=40;
int Config::minimumdepth = 5;
double Config::recallpercent = 0.9;
// 0 = 00 means the reverse of r1 and the reverse of r2 are matepairs.
// 1 = 01 means the reverse of r1 and the forward of r2 are matepairs.
// 2 = 10 means the forward of r1 and the reverse of r2 are matepairs.
// 3 = 11 means the forward of r1 and the forward of r2 are matepairs.
UINT8 Config::pairedEndReadsInputMode=8;//0-reverse,reverse, 1-reverse,forward, 2-forward,reverse, 3-forward,forward
UINT16 Config::hashKeyLength = 39;//key needs to be smaller than minimumoverlap, if double key, 2*key<minoverlap
UINT16 Config::hashKeyLength_left = 19;
UINT16 Config::hashKeyLength_right = 20;
UINT64 Config::streamChunkSize = 400;
UINT64 Config::queryChunkSize = 800;
bool Config::isSingleKey = true;
bool Config::useID = false;
UINT64 Config::startID = 0;
bool Config::isPrintIDmap = false;
string Config::IDmapfilename = "idmap.txt";
bool Config::useErrorRate = false;
bool Config::withTransitiveEdgeReduction = true;//if the instance is ConstructOverlapGraph

bool Config::isSingleEnd = true; //subject is always treated as single end, while the query is only paired end when using "MergePairedEndReads"
bool Config::isFilter = false;

bool Config::perfectMatch = false;
UINT16 Config::maxMismatch = 1; //only valid when perfectMatch is set to false
UINT16 Config::maxErrorRate = 1; //only valid when useErrorRate is set to true, in percentage
UINT16 Config::maxIndel = 0; //only valid when perfectMatch is set to false

UINT16 Config::numberOfThreads = 1;




void Config::printHelp()
{
	string command = Config::operationCode;

	if(command == "RemoveContainedReads")
	{
	    std::cout << std::endl
	              << "  Usage:" << std::endl
	              << "    storm -i/--instance RemoveContainedReads [OPTION]...<PARAM>..." << std::endl
	              << std::endl
	              << "  <PARAM>" << std::endl
                  << "    -q/--query\t query file name (usually it's a small piece of large subject file(s))" << std::endl  // file in fasta/fastq format.
	              << "    -s/--subject\t subject file name(s) (comma separated)" << std::endl  // Single-end file name list
	              << "    -ht/--hashtable single/double/omega\t single hash table or double hash table method or omega hash table method" << std::endl
	              << "    in RemoveContainedReads omega is default" << std::endl

	              << "      single : please set up the key length -k" << std::endl
	              << "      double : please set up the left key length -lk and right key length -rk" << std::endl
	              << "      omega : key length will be minimum_overlap-1, will ignore -m setting because no mismatch is allowed." << std::endl

				  << "    -l\t minimum overlap length (default:40)" << std::endl
				  << "    -k\t single hash table key length (default:39, on need in omega hash table)" << std::endl
	              << "    -lk\t double hash table left key length (default:19)" << std::endl
	              << "    -rk\t double hash table right key length (default:20)" << std::endl
				  << "    -m\t maximum allowed mismatch (default:1)" << std::endl
				  << "    -mr\t maximum allowed mismatch rate in percentage [1 means 1% of the overlap length] " << std::endl
				  << "    -t\t number of threads (default:1 [single thread])" << std::endl
				  << "    -z\t stream chunk size of subject read file (default:400)" << std::endl
				  << "    -o/--out\t output file name (default:out.txt)" << std::endl

	              << std::endl
	              << "  [OPTION]" << std::endl
	              << "    -h/--help\t only print out the help contents" << std::endl
				  << "    -id/--ID\t use and print IDs in the fasta file rather than the names. It should be followed by a specific number,the ID numbering will start from number+1" << std::endl
				  << "    --IDmap\t print the ID mapping file, followed by the file name"<< std::endl
	//			  << "    -f/--filter" << std::endl
	//			  << "    -a/--align" << std::endl
				  << std::endl

		<< "Example: ./storm -i RemoveContainedReads --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out cleanreads.fasta -l 40 -t 4" << std::endl
		<< "Example: ./storm -i RemoveContainedReads --ID 0 --IDmap mapfile.txt --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out cleanreads.fasta -l 40 -t 4" << std::endl
		<< "Example: ./storm -i RemoveContainedReads -ht omega --ID 100000 --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out cleanreads.fasta -l 40 -t 4" << std::endl
		<< "Example: ./storm -i RemoveContainedReads -ht single --ID 100000 --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out cleanreads.fasta -l 40 -m 1 -k 32 -t 4" << std::endl
		<< "Example: ./storm -i RemoveContainedReads -ht double --ID 100000 --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out cleanreads.fasta -l 40 -mr 1 -lk 16 -rk 16 -t 4" << std::endl;


	}
	else if(command == "ConstructOverlapGraph")
	{
	    std::cout << std::endl
	              << "  Usage:" << std::endl
	              << "    storm -i/--instance ConstructOverlapGraph [OPTION]...<PARAM>..." << std::endl
	              << std::endl
	              << "  <PARAM>" << std::endl
				  << "    --TransitiveReduction (default) or --noTransitiveReduction if the running instance is 'ConstructOverlapGraph'" <<std::endl
	              << "    -q/--query\t query file name (usually it's a small piece of large subject file(s))" << std::endl  // file in fasta/fastq format.
	              << "    -s/--subject\t subject file name(s) (comma separated)" << std::endl  // Single-end file name list
	              << "    -ht/--hashtable single/double/omega\t single hash table or double hash table method or omega hash table method" << std::endl
	              << "    single hash table method is default setting" << std::endl

	              << "      single : please set up the key length -k" << std::endl
	              << "      double : please set up the left key length -lk and right key length -rk" << std::endl
	              << "      omega : key length will be minimum_overlap-1, will ignore -m setting because no mismatch is allowed." << std::endl

				  << "    -l\t minimum overlap length (default:40)" << std::endl
	              << "    -k\t single hash table key length (default:39, on need in omega hash table)" << std::endl
	              << "    -lk\t double hash table left key length (default:19)" << std::endl
	              << "    -rk\t double hash table right key length (default:20)" << std::endl
				  << "    -m\t maximum allowed mismatch (default:1)" << std::endl
				  << "    -mr\t maximum allowed mismatch rate in percentage [1 means 1% of the overlap length] " << std::endl
				  << "    -t\t number of threads (default:1 [single thread])" << std::endl
				  << "    -z\t stream chunk size of subject read file (default:400)" << std::endl
				  << "    -o/--out\t output file name (default:out.txt)" << std::endl

	              << std::endl
	              << "  [OPTION]" << std::endl
	              << "    -h/--help\t only print out the help contents" << std::endl
	//			  << "    -f/--filter" << std::endl
	//			  << "    -a/--align" << std::endl
				  << std::endl

		<< "Example: ./storm -i ConstructOverlapGraph -ht omega --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out outgraph.txt -l 40 -t 4" << std::endl
		<< "Example: ./storm -i ConstructOverlapGraph -ht single --TransitiveReduction --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out outgraph.txt -l 40 -m 1 -k 32 -t 4 -z 1000" << std::endl
		<< "Example: ./storm -i ConstructOverlapGraph -ht double --noTransitiveReduction --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out outgraph.txt -l 40 -mr 1 -lk 16 -rk 16 -t 4" << std::endl;


	}
	else if(command == "ConstructGraph")
	{
	    std::cout << std::endl
	              << "  Usage:" << std::endl
	              << "    storm -i/--instance ConstructGraph [OPTION]...<PARAM>..." << std::endl
	              << std::endl
	              << "  <PARAM>" << std::endl
	              << "    -q/--query\t query file name (usually it's a small piece of large subject file(s))" << std::endl  // file in fasta/fastq format.
	              << "    -s/--subject\t subject file name(s) (comma separated)" << std::endl  // Single-end file name list
	              << "    -ht/--hashtable single/double/omega\t single hash table or double hash table method or omega hash table method" << std::endl
	              << "    single hash table method is default setting" << std::endl

	              << "      single : please set up the key length -k" << std::endl
	              << "      double : please set up the left key length -lk and right key length -rk" << std::endl
	              << "      omega : key length will be minimum_overlap-1, will ignore -m setting because no mismatch is allowed." << std::endl

	              << "    -l\t minimum overlap length (default:40)" << std::endl
				  << "    -k\t single hash table key length (default:39, on need in omega hash table)" << std::endl
	              << "    -lk\t double hash table left key length (default:19)" << std::endl
	              << "    -rk\t double hash table right key length (default:20)" << std::endl
				  << "    -m\t maximum allowed mismatch (default:1)" << std::endl
				  << "    -mr\t maximum allowed mismatch rate in percentage [1 means 1% of the overlap length] " << std::endl
				  << "    -t\t number of threads (default:1 [single thread])" << std::endl
				  << "    -z\t stream chunk size of subject read file (default:400)" << std::endl
				  << "    -o/--out\t output file name (default:out.txt)" << std::endl

	              << std::endl
	              << "  [OPTION]" << std::endl
	              << "    -h/--help\t only print out the help contents" << std::endl
	//			  << "    -f/--filter" << std::endl
	//			  << "    -a/--align" << std::endl
				  << std::endl

		<< "Example: ./storm -i ConstructGraph -ht omega --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out outgraph.txt -l 40 -t 4" << std::endl
	    << "Example: ./storm -i ConstructGraph -ht single --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out outgraph.txt -l 40 -m 1 -k 32 -t 4 -z 1000" << std::endl
		<< "Example: ./storm -i ConstructGraph -ht double --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out outgraph.txt -l 40 -mr 1 -lk 16 -rk 16 -t 4" << std::endl;


	}
	else if(command == "MergePairedEndReads")
	{
	    std::cout << std::endl
	              << "  Usage:" << std::endl
	              << "    storm -i/--instance MergePairedEndReads [OPTION]...<PARAM>..." << std::endl
	              << std::endl
	              << "  <PARAM>" << std::endl
	              << "    -q/--query\t query file name (usually it's a small piece of large subject file(s))" << std::endl  // file in fasta/fastq format.
				  << "	  -or/--orient\t orientation for the paired end reads in query file" << std::endl
				  << "	  0-reverse,reverse; 1-reverse,forward; 2-forward,reverse (usually for Illumina reads); 3-forward,forward"<<std::endl
	              << "    -s/--subject\t subject file name(s) (comma separated)" << std::endl  // Single-end file name list
	              << "    -ht/--hashtable single/double/omega\t single hash table or double hash table method or omega hash table method" << std::endl
	              << "    single hash table method is default setting" << std::endl
				  << "	  omega is not allowed in MergePairedEndReads, only single or double hash table is allowed" <<std::endl

	              << "      single : please set up the key length -k" << std::endl
	              << "      double : please set up the left key length -lk and right key length -rk" << std::endl


				  << "    -ll\t minimum overlap length (default:20, only needed in MergePairedEndReads)" << std::endl
				  << "    -rl\t minimum overlap length (default:20, only needed in MergePairedEndReads)" << std::endl
				  << "    -ol\t minimum overlap length (default:40, only needed in MergePairedEndReads)" << std::endl
	              << "    -k\t single hash table key length (default:39)" << std::endl
	              << "    -lk\t double hash table left key length (default:19)" << std::endl
	              << "    -rk\t double hash table right key length (default:20)" << std::endl
				  << "    -m\t maximum allowed mismatch (default:1)" << std::endl
				  << "    -mr\t maximum allowed mismatch rate in percentage [1 means 1% of the overlap length] " << std::endl
				  << "    -c\t minum coverage depth required to merge the reads (default: 5) " << std::endl
				  << "    -b\t nucleotide recall percentage (default: 0.9) " << std::endl
				  << "    -t\t number of threads (default:1 [single thread])" << std::endl
				  << "    -z\t stream chunk size of subject read file (default:400)" << std::endl
				  << "    -o/--out\t output file name (default:out.txt)" << std::endl

	              << std::endl
	              << "  [OPTION]" << std::endl
	              << "    -h/--help\t only print out the help contents" << std::endl
	//			  << "    -f/--filter" << std::endl
	//			  << "    -a/--align" << std::endl
				  << std::endl

					<< "Example: ./storm -i MergePairedEndReads -ht single --query qreads.fasta --orient 2 --subject sreads1.fasta,sreads2.fasta --out outreads -m 0 -ll 30 -rl 30 -ol 60 -k 29 -t 4" << std::endl
					<< "Example: ./storm -i MergePairedEndReads -ht double --query qreads.fasta --orient 2 --subject sreads1.fasta,sreads2.fasta --out outreads -m 1 -ll 10 -rl 10 -ol 20 -lk 4 -rk 4 -t 4"<< std::endl
					<< "Example: ./storm -i MergePairedEndReads -ht single --query qreads.fasta --orient 2 --subject sreads1.fasta,sreads2.fasta --out outreads -m 0 -c 3 -b 0.8 -ll 30 -rl 30 -ol 60 -k 29 -t 4" << std::endl
					<< "Example: ./storm -i MergePairedEndReads -ht double --query qreads.fasta --orient 2 --subject sreads1.fasta,sreads2.fasta --out outreads -m 1 -c 5 -b 0.9 -ll 10 -rl 10 -ol 20 -lk 4 -rk 4 -t 4"<< std::endl;


	}
	else
	{
    std::cout << std::endl
              << "  Usage:" << std::endl
              << "    storm [OPTION]...<PARAM>..." << std::endl
              << std::endl
              << "  <PARAM>" << std::endl
			  << "    -i/--instance\t RemoveContainedReads,ConstructGraph,ConstructOverlapGraph,MergePairedEndReads [or CorrectErrors(NOT YET)]" << std::endl
			  << "    --TransitiveReduction (default) or --noTransitiveReduction if the running instance is 'ConstructOverlapGraph'" <<std::endl
			  << "	  Therefore, 'ConstructOverlapGraph' uses different framework than 'ConstructGraph' and it only reports out-stretched edges rather than all types of edges" << std::endl
              << "    -q/--query\t query file name (usually it's a small piece of large subject file(s))" << std::endl  // file in fasta/fastq format.
              << "    -s/--subject\t subject file name(s) (comma separated)" << std::endl  // Single-end file name list
              << "    -ht/--hashtable single/double/omega\t single hash table or double hash table method or omega hash table method" << std::endl
              << "    single hash table method is default setting except that in RemoveContainedReads omega is default" << std::endl
			  << "	  omega is not allowed in MergePairedEndReads, only single or double hash table is allowed" <<std::endl

              << "      single : please set up the key length -k" << std::endl
              << "      double : please set up the left key length -lk and right key length -rk" << std::endl
              << "      omega : key length will be minimum_overlap-1, will ignore -m setting because no mismatch is allowed." << std::endl

			  << "	  -or/--orient\t orientation for the paired end reads in query file(required in MergePairedEndReads)" << std::endl
			  << "	  0-reverse,reverse; 1-reverse,forward; 2-forward,reverse (usually for Illumina reads); 3-forward,forward "<<std::endl

              << "    -l\t minimum overlap length (default:40, no need in MergePairedEndReads)" << std::endl
			  << "    -ll\t minimum overlap length (default:20, only needed in MergePairedEndReads)" << std::endl
			  << "    -rl\t minimum overlap length (default:20, only needed in MergePairedEndReads)" << std::endl
			  << "    -ol\t minimum overlap length (default:40, only needed in MergePairedEndReads)" << std::endl
              << "    -k\t single hash table key length (default:39)" << std::endl
              << "    -lk\t double hash table left key length (default:19)" << std::endl
              << "    -rk\t double hash table right key length (default:20)" << std::endl
			  << "    -m\t maximum allowed mismatch (default:1)" << std::endl
			  << "    -mr\t maximum allowed mismatch rate in percentage [1 means 1% of the overlap length] " << std::endl
			  << "    -t\t number of threads (default:1 [single thread])" << std::endl
			  << "    -z\t stream chunk size of subject read file (default:400)" << std::endl
			  << "    -o/--out\t output file name (default:out.txt)" << std::endl

              << std::endl
              << "  [OPTION]" << std::endl
              << "    -h/--help\t only print out the help contents" << std::endl
			  << "    -id/--ID\t use and print IDs in the fasta file rather than the names. It should be followed by a specific number,the ID numbering will start from number+1" << std::endl
//			  << "    -f/--filter" << std::endl
//			  << "    -a/--align" << std::endl
			  << std::endl

				<< "Example: ./storm -i RemoveContainedReads --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out cleanreads.fasta -l 40 -t 4" << std::endl
				<< "Example: ./storm -i RemoveContainedReads --ID 0 --IDmap mapfile.txt --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out cleanreads.fasta -l 40 -t 4" << std::endl
				<< "Example: ./storm -i RemoveContainedReads -ht omega --ID 100000 --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out cleanreads.fasta -l 40 -t 4" << std::endl
				<< "Example: ./storm -i RemoveContainedReads -ht single --ID 100000 --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out cleanreads.fasta -l 40 -m 1 -k 32 -t 4" << std::endl
				<< "Example: ./storm -i RemoveContainedReads -ht double --ID 100000 --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out cleanreads.fasta -l 40 -mr 1 -lk 16 -rk 16 -t 4" << std::endl
				<< "Example: ./storm -i ConstructOverlapGraph -ht omega --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out outgraph.txt -l 40 -t 4" << std::endl
				<< "Example: ./storm -i ConstructOverlapGraph -ht single --TransitiveReduction --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out outgraph.txt -l 40 -m 1 -k 32 -t 4 -z 1000" << std::endl
				<< "Example: ./storm -i ConstructOverlapGraph -ht double --noTransitiveReduction --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out outgraph.txt -l 40 -mr 1 -lk 16 -rk 16 -t 4" << std::endl
				<< "Example: ./storm -i ConstructGraph -ht omega --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out outgraph.txt -l 40 -t 4" << std::endl
			    << "Example: ./storm -i ConstructGraph -ht single --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out outgraph.txt -l 40 -m 1 -k 32 -t 4 -z 1000" << std::endl
				<< "Example: ./storm -i ConstructGraph -ht double --query qreads.fasta --subject sreads1.fasta,sreads2.fasta --out outgraph.txt -l 40 -mr 1 -lk 16 -rk 16 -t 4" << std::endl
				<< "Example: ./storm -i MergePairedEndReads -ht single --query qreads.fasta --orient 2 --subject sreads1.fasta,sreads2.fasta --out outreads -m 0 -ll 30 -rl 30 -ol 60 -k 29 -t 4" << std::endl
				<< "Example: ./storm -i MergePairedEndReads -ht double --query qreads.fasta --orient 2 --subject sreads1.fasta,sreads2.fasta --out outreads -m 1 -ll 10 -rl 10 -ol 20 -lk 4 -rk 4 -t 4"<< std::endl
				<< "Example: ./storm -i MergePairedEndReads -ht single --query qreads.fasta --orient 2 --subject sreads1.fasta,sreads2.fasta --out outreads -m 0 -c 3 -b 0.8 -ll 30 -rl 30 -ol 60 -k 29 -t 4" << std::endl
				<< "Example: ./storm -i MergePairedEndReads -ht double --query qreads.fasta --orient 2 --subject sreads1.fasta,sreads2.fasta --out outreads -m 1 -c 5 -b 0.9 -ll 10 -rl 10 -ol 20 -lk 4 -rk 4 -t 4"<< std::endl;

	}
}

bool Config::setConfig(int argc, char **argv)
{
	Config::subjectFilenameList.clear();

	vector<string> argumentsList;
	cout << "PRINTING ARGUMENTS" << endl;
	for(int i = 0; i < argc; i++)
	{
		cout << argv[i] << ' ';
	}
	cout << endl;
	while(argc--)
			argumentsList.push_back(*argv++);

	if(argumentsList.size() == 1)
	{
        Config::printHelp();
		return false;
	}

	for(UINT64 i = 1; i <= argumentsList.size()-1; i++)
	{
		if (argumentsList[i] == "-h" || argumentsList[i] == "--help")
	    {
	        Config::printHelp();
			return false;
	    }
/*		if (argumentsList[i] == "-f" || argumentsList[i] == "--filter")
	    {
	        Config::isFilter = true;

	    }
		if (argumentsList[i] == "-a" || argumentsList[i] == "--align")
	    {
	        Config::isFilter = false;

	    }
*/

		if (argumentsList[i] == "-id" || argumentsList[i] == "--ID")
	    {
	        Config::useID = true;
	        Config::startID = atol(argumentsList[++i].c_str());
	    }
		else if (argumentsList[i] == "--IDmap")
	    {
	        Config::isPrintIDmap = true;
	        Config::IDmapfilename = argumentsList[++i];
	    }
		else if(argumentsList[i] == "--TransitiveReduction" )
		{
			Config::withTransitiveEdgeReduction = true;
		}
		else if(argumentsList[i] == "--noTransitiveReduction" )
		{
			Config::withTransitiveEdgeReduction = false;
		}
		else if (argumentsList[i] == "-i" || argumentsList[i] == "--instance")
	    {
	        Config::operationCode = argumentsList[++i];

	    }
		else if (argumentsList[i] == "-or" || argumentsList[i] == "--orient")
	    {
			Config::pairedEndReadsInputMode = atoi(argumentsList[++i].c_str());

	    }
		else if(argumentsList[i] == "-s" || argumentsList[i] == "--subject")
		{
			string inputFilenames=argumentsList[++i];
            stringstream ss(inputFilenames);
            string item;

			while (getline(ss, item, ','))
			{
				Config::subjectFilenameList.push_back(item);
			}
		}
		else if(argumentsList[i] == "-q" || argumentsList[i] == "--query")
		{
			string inputFilename=argumentsList[++i];
			Config::queryFilename = inputFilename;

		}
		else if(argumentsList[i] == "-o" || argumentsList[i] == "--out")
		{
			string outputFilename=argumentsList[++i];
			Config::outputfilename = outputFilename;

		}
		else if(argumentsList[i] == "-ht" || argumentsList[i] == "--hashtable")
		{
			string type=argumentsList[++i];
			Config::hashtabletype = type;

		}

		else if (argumentsList[i] == "-l")
			Config::minimumOverlapLength = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-ll")
		{
			Config::left_minimumOverlapLength = atoi(argumentsList[++i].c_str());
			Config::minimumOverlapLength = Config::left_minimumOverlapLength;
		}
		else if (argumentsList[i] == "-rl")
			Config::right_minimumOverlapLength = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-ol")
			Config::overall_minimumOverlapLength = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-k")
			Config::hashKeyLength = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-lk")
			Config::hashKeyLength_left = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-rk")
			Config::hashKeyLength_right = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-m")
					Config::maxMismatch = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-mr")
		{
			Config::useErrorRate = true;
			Config::maxErrorRate = atoi(argumentsList[++i].c_str());
		}
		else if (argumentsList[i] == "-c")
				Config::minimumdepth = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-b")
				Config::recallpercent = atof(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-t")
					Config::numberOfThreads = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-z")
					Config::streamChunkSize = atoi(argumentsList[++i].c_str());
		else
		{
          Config::printHelp();
          return false;
		}
	}



    if( (Config::subjectFilenameList.size() == 0) || Config::queryFilename == "")
    {
        cout << "missed -subject or -query input files!" << std::endl;
        printHelp();
		return false;
    }
    if(Config::operationCode == "")
    {
    	cout<< "operation instance code is missing" << std::endl;
    	printHelp();
    	return false;
    }

    if(Config::operationCode=="MergePairedEndReads"&&Config::pairedEndReadsInputMode==8)
    {
    	cout<< "please specify --orient/-or for paired end reads orientation in your input query file " << std::endl;
    	printHelp();
    	return false;

    }

    //default setting for the hash table type
    if(Config::hashtabletype=="")
    {
    	if(Config::getOperation()=="RemoveContainedReads")
    		Config::hashtabletype="omega";
    	else
    		Config::hashtabletype="single";
    }

    if(Config::hashtabletype=="omega")
    {
    	Config::maxMismatch=0;Config::maxErrorRate=0;

    }



    return true;
}

string Config::getOperation()
{
	return Config::operationCode;
}
vector<string> Config::getSubjectDatasetFilenames()
{
	return Config::subjectFilenameList;
}
string Config::getQueryDatasetFilename()
{
	return Config::queryFilename;
}
UINT16 Config::getminimumOverlapLength()
{
	return Config::minimumOverlapLength;
}
UINT16 Config::getleftminimumOverlapLength()
{
	return Config::left_minimumOverlapLength;
}
UINT16 Config::getrightminimumOverlapLength()
{
	return Config::right_minimumOverlapLength;
}
UINT16 Config::getHashKeyLength()
{
	return Config::hashKeyLength;
}
UINT16 Config::getHashLeftKeyLength()
{
	return Config::hashKeyLength_left;
}
UINT16 Config::getHashRightKeyLength()
{
	return Config::hashKeyLength_right;
}
UINT64 Config::getStreamChunkSize()
{
	return Config::streamChunkSize;
}

bool Config::isSingleKeyHashTable()
{
	return Config::isSingleKey;
}
string Config::getHashTableType()
{
	return Config::hashtabletype;
}
bool Config::allowMismatch()
{
	return !Config::perfectMatch;
}
bool Config::isQuerySingleEndReads()
{
	return Config::isSingleEnd;
}

UINT16 Config::getMaxMismatch()
{
	return Config::maxMismatch;
}
UINT16 Config::getMaxIndel()
{
	return Config::maxIndel;
}
bool Config::isFilterOrAlign()
{
	return Config::isFilter;
}
