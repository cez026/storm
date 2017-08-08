# storm
STORM (A Scalable Tool for Merging Paired-End Reads with Variable Insert Sizes)
This is STORM version 1.5

* STORM is a tool for merging paired-end reads with variable insert sizes (i.e. for both overlapping and non-overlapping pairs)
* STORM results in longer reads with higher accuracy compared to current popular merging tools
* STORM is also able to remove contained reads or identical reads
* STORM can generate overlap graphs as input of Omega2 

https://github.com/qiumingyao/omega2

https://bitbucket.org/omicsbio/omega2


### Set up / Install Instructions ###

* Download the package
* Go to the Release folder, and "make"
* run "./storm -h"

### Parameters ###
./storm 

  Usage:
    storm [OPTION]...[PARAM]...


  [PARAM]

    -i/--instance	 MergePairedEndReads, RemoveContainedReads, ConstructOverlapGraph

    -q/--query	 query file name (usually it's a small piece of large subject file(s))

    -s/--subject	 subject file name(s) (comma separated)

    -ht/--hashtable single/double	 single hash table or double hash table method

       single hash table method is default setting

      single : please set up the key length -k

      -k	 single hash table key length (default:39)

      double : please set up the left key length -lk and right key length -rk

      -lk	 double hash table left key length (default:19)

      -rk	 double hash table right key length (default:20)

    -or/--orient	 orientation for the paired end reads in query file(required in MergePairedEndReads)
	  0-reverse,reverse; 1-reverse,forward; 2-forward,reverse (usually for Illumina reads); 3-forward,forward 

    -l	 minimum overlap length (default:40, no need in MergePairedEndReads)

    -ll	 minimum overlap length (default:20, only needed in MergePairedEndReads)

    -rl	 minimum overlap length (default:20, only needed in MergePairedEndReads)

    -ol	 minimum overlap length (default:40, only needed in MergePairedEndReads)

    -m	 maximum allowed mismatch (default:1)

    -mr	 maximum allowed mismatch rate in percentage [1 means 1% of the overlap length] 

    -t	 number of threads (default:1 [single thread])

    -z	 stream chunk size of subject read file (default:400)

    -o/--out	 output file name (default:out.txt)

  [OPTION]

    -h/--help	 only print out the help contents

### Examples ###

Example1: ./storm -i MergePairedEndReads -ht single --query qreads.fasta --orient 2 --subject sreads1.fasta,sreads2.fasta --out outreads -m 0 -ll 30 -rl 30 -ol 60 -k 29 -t 4

Example2: ./storm -i MergePairedEndReads -ht double --query qreads.fasta --orient 2 --subject sreads1.fasta,sreads2.fasta --out outreads -m 1 -ll 10 -rl 10 -ol 20 -lk 4 -rk 4 -t 4

Example3: ./storm -i MergePairedEndReads -ht single --query qreads.fasta --orient 2 --subject sreads1.fasta,sreads2.fasta --out outreads -m 0 -c 3 -b 0.8 -ll 30 -rl 30 -ol 60 -k 29 -t 4

Example4: ./storm -i MergePairedEndReads -ht double --query qreads.fasta --orient 2 --subject sreads1.fasta,sreads2.fasta --out outreads -m 1 -c 5 -b 0.9 -ll 10 -rl 10 -ol 20 -lk 4 -rk 4 -t 4

### Contacts ###

Qiuming Yao: yao.ornl@gmail.com

Chongle Pan: panc@ornl.gov
