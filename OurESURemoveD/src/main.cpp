//latest version
// modified on Jun 26, 2016
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <math.h>
#include "graph.h"
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <list>

#include "MyTree.h"
#include <ext/hash_map>

#define EPOC 10
//#define Debug

int frequency_thr = -1;
unsigned long long subgraph_THR = 0xffffffffffff;
int* subgraphDegree;
int* graphDegs;
hash_map<char,int> AsciiToInt;
hash_map<int,char> IntToAscii;
hash_map<std::string,long long int> graphInt;
hash_map<std::string,long long int> reordered_graphInt;
hash_map<std::string,long long int> treeInt;
using namespace std;

int subgraphSize = -1, num_random_graphs = 0;
float subgraphDensity = 2;
unsigned long long callNautyCount = 0;

//float AdjStrTime = 0, densityCheckTime = 0;

//g stores the input graph
Graph *g;

//isRand determines whether the enumeration of random networks is commenced or not
//directed indicates whether the input network is directed or not
bool isRand, directed;

extern unsigned long enumerated_class;

int* treeChildrenSize;
string* BFSVec;

/****************************************************************
****************************************************************/

void print_usage (FILE * stream, int exit_code)
{
	fprintf (stream, "Usage: OurESU options[inputfile...] \n ");
    fprintf (stream,
		 "\t-h	--help\t\t\tDisplay this usage information. \n"
		 "\t-i	--input filename\tInput filename.\n"
		 "\t-o	--output path\t\tOutput directory.\n"
		 "\t-r 	--random number \tNumber of random graphs (default = 0).\n"
		 "\t-s 	--size motifsize \tMotif size.\n"
		 "\t-u 	--undirected \tUndirected input network\n"
		 "\t-l	--limit \tlimitation of Enumerated subgraphs\n"
		 "\t-d	--density \tdensity of subgraph\n"
		 "\t-f	--threshold \tthreshold of frequency\n");
	     
    exit (exit_code);
}

/****************************************************************
****************************************************************/

bool ReadData(const char *path) {
	
	

	register int i, j;
	int graphSize;
	FILE * inFile = fopen(path, "r");
	
	if (inFile == NULL) {
		printf("Error opening %s file.\n", path);
		return false;
	}
	
	if (!feof(inFile))
		fscanf(inFile, "%d\n", &graphSize);	

	printf("graphSize = %d\n", graphSize);
	subgraphDegree = new int[graphSize + 1];
	//graphDegs = new int[graphSize +1];
	
	g = new Graph(graphSize, subgraphSize);
	
	while (!feof(inFile)) {
		fscanf(inFile, "%d %d\n", &i, &j);
		if(i == j) continue;
		g->addEdge(i, j);
		
		if(!directed)
		{
			g->addEdge(j, i);
		}
	}



	g->Finalize(directed);
	fclose(inFile);

	return true;
}

/***********************************************************************************
 * This function enumerates the subgraphs related to each vertex of inpur network. *
***********************************************************************************/
void Enumerate() {
	register int v;
	//g->preComC = 0;

	Subgraph* sub = new Subgraph(1, subgraphSize, g->Size());// creat 3 arrays: visited[154],vertices[4], children[154]
	for (v = 1; v <= g->Size(); v++)
	{
#ifdef Debug
		printf("+ Exploring Node %d ...\n", v);
#endif Debug

		// if(g->subgraphCounter >= 100)
		// 	break;//useful

		sub->subgraphSize = 1;
		sub->lastVertex = sub->vertices[0] = v;

		//sub->verticesEdges[0] = 0;//Tao add
		sub->verticesAdj[0] = "";

		
		sub->visited[v] = true;
		//cout << "Before Nexts, everything is OK"<<endl;
		//g->Nexts(sub, subgraphSize, 0, g->quaT->root); 
		g->Nexts(sub, subgraphSize, 0);

		sub->visited[v] = false;
		
	}

	delete sub;
}
/****************************************************************
****************************************************************/
std::string graphDegreeSequence(string adj, int subgraphSize)
{

    register int i, j, index = 0;
    //int totalLength = adj.size();
    vector<int> degreeVec;
    for(i = 0; i<subgraphSize; i++)
        degreeVec.push_back(0);



    for(i = 0; i < subgraphSize; i++)
    {
        for(j = i+1; j<subgraphSize; j++)
        {
        	index = j*(j-1)/2+i;
            if(adj[index] == '1')
            {
                degreeVec[i] += 1;
                degreeVec[j] += 1;
            }
            //index ++;
        }
    }
    // for(i = 0; i < subgraphSize; i++)
    // {
    //     for(j = i+1; j<subgraphSize; j++)
    //     {
    //         if(adj[index] == '1')
    //         {
    //             degreeVec[i] += 1;
    //             degreeVec[j] += 1;// undirected!
    //         }
    //         index ++;
    //     }
    // }
   
    std::sort(degreeVec.begin(), degreeVec.end());
   
    string ds = "";
    ds.reserve(subgraphSize);
    for(vector<int>::iterator it = degreeVec.begin(); it != degreeVec.end(); it++)
    {
        ds.push_back(IntToAscii.find(*it-1)->second); 
    }

   
    return ds;
}
bool isTree(string degreeSeq, int subgraphSize)
{
	
	int degree_sum = 0;
	for (string::iterator it = degreeSeq.begin(); it != degreeSeq.end(); it++)
	{
		degree_sum = degree_sum + AsciiToInt.find(*it+1)->second;
	}
	
	if (degree_sum == 2*subgraphSize -2)
		return TRUE;
	else
		return FALSE;
	
}
/****************************************************************
****************************************************************/
string calculateCam(string graphIn, int subgraphSize)
{
    //cout<<"graphIn"<<graphIn<<endl;
    string tgraph(subgraphSize*subgraphSize,'0');
   

    register int i, j, index;
    for(i = 0; i < subgraphSize; i++)
    {
        for (j = i+1; j<subgraphSize; j++)
        {
            index = j*(j-1)/2 +i;
            if(graphIn[index] == '2')
            {
                tgraph[i*subgraphSize + j] = '1';
                tgraph[j*subgraphSize + i] = '1';//undirected
            }
            else if(graphIn[index] == '1')
            {
            	tgraph[i*subgraphSize + j] = '1';
            }
            else if(graphIn[index] == '3')
            {
            	tgraph[j*subgraphSize + i] = '1';
            }
        }
    }
    // index = 0;
    // for(i = 0; i < subgraphSize; i++)
    // {
    //     for (j = i+1; j<subgraphSize; j++)
    //     {
    //         if(graphIn[index] == '1')
    //         {
    //             tgraph[i*subgraphSize + j] = '1';
    //             tgraph[j*subgraphSize + i] = '1';//undirected
    //         }
    //         index++;
    //     }
    // }

    int  n, m, v, k;
    unsigned nCanCode;
    char sCanCode[subgraphSize*subgraphSize];
    set *gv;
    graph g[subgraphSize*subgraphSize];
    graph canong[subgraphSize*subgraphSize];
    nvector lab[subgraphSize],ptn[subgraphSize],orbits[subgraphSize];
    static DEFAULTOPTIONS(options);
    setword workspace[160*subgraphSize];
    /*init for nauty*/
    options.writeautoms = FALSE;
    options.writemarkers = FALSE;
    options.getcanon = TRUE;
    options.defaultptn = TRUE;//wangtao
    options.digraph = TRUE;// FALSE->True,modified by wangtao 4/20th/2015
    statsblk(stats);
    

    n = subgraphSize;
    m = (n+WORDSIZE-1)/WORDSIZE;

    for(v=0; v<n; v++)
    {
        gv=GRAPHROW(g,v,m);
        EMPTYSET(gv,m);
        for(i=0; i<n; i++)
            if(tgraph[v*subgraphSize+i] != '0')//here make sure that the graph is undirected!!!no,depends on tgraph.
                ADDELEMENT(gv,i);
    }//initialize matrix g

   
  
    nauty(g,lab,ptn,NILSET,orbits,&options,&stats,workspace,160*subgraphSize,m,n,canong);
    nCanCode = 0;
    k=0;
    for(i=0; i<n; i++)
    {
        gv=GRAPHROW(canong,i,m);
        for(j=0; j<n; j++)
        {
            nCanCode = nCanCode<<1;
            nCanCode+=ISELEMENT(gv,j);
            sCanCode[k++] = (char)(ISELEMENT(gv,j)+48);//ascii(0) == 48
        }
    }//output cam string
    sCanCode[k]='\0';
    string re = sCanCode;
    
    return re;
}
/****************************************************************
****************************************************************/
bool isConnected(int i, int j, string adjStr)
{
	int a, b, index;
	if(i < j)
	{
		a = i;
		b = j;
	}
	else
	{
		a = j;
		b = i;
	}
	index = (subgraphSize-1)*a - a*(a+1)/2 + b -1;//(2*(subgraphSize-1)-a +1)*a/2 +b-a-1;
	if(adjStr[index] == '1')
		return TRUE;
	else
		return FALSE;
}
bool isConnected2(int i, int j, string adjStr)
{
	int a, b, index;
	if(i < j)
	{
		a = i;
		b = j;
	}
	else
	{
		a = j;
		b = i;
	}
	index = b*(b-1)/2 +a;
	if(adjStr[index] != '0')
		return TRUE;
	else
		return FALSE;
}
string subGetTreeString(int* subVertices, int root, string adjStr )
{
		
		MyTree tree(subgraphSize);
	
		bool used[subgraphSize];
		for(int i = 0; i< subgraphSize; i++)
		{
			used[i] = false;
		}
		list<int> queue;
  		queue.push_back(1);//always transform the root as 1
  		int lookup[subgraphSize+1];
  		lookup[0] = 0;//no sense
  		lookup[1] = root;
  		used[root] = true;// mark root as used, root's index is root1Ind
  		int ind = 2; // other nodes coding from 2

  		while(!queue.empty())
  		{
  			int k = queue.front();
  			queue.pop_front();
  			for (int i = 0; i < subgraphSize; ++i)
  			{
  				if (isConnected2(lookup[k], subVertices[i], adjStr) && used[i]==false)
  				{
  					tree.addEdge(k, ind);
  					queue.push_back(ind);
  					used[i] = true;
  					lookup[ind] = subVertices[i];
  					ind++;
  					
  				}
  			}

  		}

  		treeChildrenSize = new int[subgraphSize+1];
  		for(int i = 1; i< subgraphSize +1; i++)
		{
		   treeChildrenSize[i] = tree.ChildrenNum(i);
		}

		BFSVec = new string[subgraphSize+1];//Here!!!
		
		string outStr = tree.TreeCamGen(subgraphSize); //tree1.TreeCamGen(int** subgraphTree, int level);
		
		delete [] BFSVec;
		delete [] treeChildrenSize;
		return outStr;
}


string GetTreeCanStr(string adjStr, int subgraphSize)
{
	register int i,j, ind = 0;
	int* subgraphDegree = new int[subgraphSize];
	int* tempSubgraph = new int[subgraphSize];
	for (int i = 0; i < subgraphSize; ++i)
	{
		subgraphDegree[i] = 0;
		tempSubgraph[i] = i;
	}
	// for(i = 0; i < subgraphSize; i++)
	// {
	// 	for (j = i+1; j < subgraphSize; j++)
	// 	{
	// 		if(adjStr[ind] == '1')
	// 		{
	// 			subgraphDegree[i] += 1;
	// 			subgraphDegree[j] += 1; //undirected!!
	// 		}
	// 		ind++;
	// 	}
	// }
	//cout<<"start"<<endl;	
	for(i = 0; i < subgraphSize; i++)
    {
        for (j = i+1; j<subgraphSize; j++)
        {
            ind = j*(j-1)/2 +i;
            if(adjStr[ind] != '0')
            {
                subgraphDegree[i] += 1;
                subgraphDegree[j] += 1;//undirected
            }
        }
    }

    // for (int i = 0; i < subgraphSize; ++i)
    // {
    // 	cout<<subgraphDegree[i]<<endl;
    // }


	int* subVerticesDegree =  new int[subgraphSize];
	for (i = 0; i < subgraphSize; ++i)
	{
		subVerticesDegree[i] = subgraphDegree[i];
		//cout<<subVerticesDegree[i]<<" ";
	}
	//cout<<endl;
	int remainNum = subgraphSize;
	int toDelNode;
	while(remainNum>2)
	{
		std::vector<int> deleted_nodes;
		for (int i = 0; i < subgraphSize; ++i)
		{
			if (subgraphDegree[tempSubgraph[i]] == 1 && tempSubgraph[i]!=-1 ) // leaf node with degree ==1
			{
				// remove tempSubgraph[i]
				toDelNode = tempSubgraph[i];
				tempSubgraph[i] = -1;
				remainNum -= 1;
				deleted_nodes.push_back(toDelNode);
				subgraphDegree[toDelNode] = 0;
			}
		
		}
		//update the degree of remaining nodes
		std::vector<int>::iterator it;
		for (it = deleted_nodes.begin(); it < deleted_nodes.end(); ++it)
		{
			int deleted_node = *it;
			//cout<<"d"<<deleted_node;
			for (int i = 0; i < subgraphSize; ++i)
			{
				if (tempSubgraph[i]!=-1 && isConnected2(tempSubgraph[i], deleted_node, adjStr))
				{
					subgraphDegree[tempSubgraph[i]] -= 1;
				}
			}

		}
		// cout<<endl;
		// for (int i = 1; i < subgraphSize+1; ++i)
		// {
		// 	cout<<subgraphDegree[i]<<"s ";
		// 	/* code */
		// }
		// cout<<endl;
	}

	int root1 = -1, root2 = -1;
	int root1Ind = -1, root2Ind = -1;
	//cout<<"remainNum: "<<remainNum<<endl;
	if (remainNum ==1 )
	{
		for (int i = 0; i < subgraphSize; ++i)
		{
			if (tempSubgraph[i]!=-1)
			{
				root1 = tempSubgraph[i];
				root1Ind = i;
			}
		}
	}
	else if (remainNum ==2)
	{
		for (int i = 0; i < subgraphSize; ++i)
		{
			if (tempSubgraph[i]!=-1 && root1 == -1)
			{
				root1 = tempSubgraph[i];
				root1Ind = i;
			}
			else if (tempSubgraph[i]!= -1 && root1 != -1)
			{
				root2 = tempSubgraph[i];
				root2Ind = i;
			}
		}
		if (subVerticesDegree[root1Ind] > subVerticesDegree[root2Ind]) // select root which has more children
		{
			root2 = -1;
		}
		else if (subVerticesDegree[root2Ind] > subVerticesDegree[root1Ind])
		{
			root1 = -1;
		}
		// else // subVerticesDegree[root1Ind] == subVerticesDegree[root2Ind]
		// {
		// 	root2 = -1;
		// }

	}


	 //cout<<"root1: " <<root1<<" root2: "<<root2<<endl;
	 //cout<<"root1Ind: " <<root1Ind<<" root2Ind: "<<root2Ind<<endl;

	string outStr1 = ">", outStr2 = ">";
	//Call Tree part, input a layered matrix, o as the root
	if (root1 == -1 && root2 == -1)
	{
		cout<<"ERROR: root1 == -1 and root2 == -1"<<endl;
		exit(1);
	}
	int* subVertices = new int[subgraphSize];
	for (int i = 0; i < subgraphSize; ++i)
	{
		subVertices[i] = i;
	}
	if (root1 != -1)
	{
		outStr1 = subGetTreeString(subVertices, root1, adjStr);
		
	}
	if (root2 != -1)
	{
		outStr2 = subGetTreeString(subVertices, root2, adjStr);
	}

	

	if (outStr1 <= outStr2)
	{
		return outStr1;
	}
	else
	{
		return outStr2;
	}
	delete [] subgraphDegree;
	delete [] tempSubgraph;
	delete [] subVerticesDegree;

}

struct ordering {
    bool operator ()(pair<int, int> const& a, 
                     pair<int, int> const& b) {
        return a.second > b.second;
    }
};


/****************************************************************
****************************************************************/
int main(int argc, char *argv[]) {

	for(int ini = -1; ini < 93; ini++)
	{
		char temp = '!' + ini;
		AsciiToInt[temp] = ini;
		IntToAscii[ini] = temp;
	}



	double total_random_time = 0 , main_time;
	clock_t start_random_time, end_random_time;
	directed = true;
	register int i, j;
	long long unsigned subgraphCounterMain;

	int next_option;
	const char *const short_options = "h:i:o:r:s:d:f:l:u";
	const struct option long_options[] = {
		{"help",   0, NULL, 'h'},
		{"input",  1, NULL, 'i'},
		{"output", 1, NULL, 'o'},
		{"random", 1, NULL, 'r'},
		{"size",   1, NULL, 's'},
		{"density", 1, NULL, 'd'},
		{"frequency", 1, NULL, 'f'},
		{"limit", 1, NULL, 'l'},
		{"undirected",   0, NULL, 'u'},
		{NULL,     0, NULL,  0 }
	};
	
	char *program_name;
    char input_filename[256], output_directory[256];

    int verbose = 0;
    strcpy(output_directory, "result");

    program_name = argv[0];
    do {
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);
	
		switch (next_option) {
			case 'h':
	    		print_usage (stdout, 0);

			case 'i':
				strcpy(input_filename, optarg);
	    		break;
	    	
			case 'o':
				strcpy(output_directory, optarg);
	    		break;
			
			case 'r':
				 num_random_graphs = atoi(optarg);
	    		break;

			case 's':
				subgraphSize = atoi(optarg);
	    		break;
	    	case 'd':
	    		subgraphDensity = atof(optarg);
	    		break;
	    	case 'f':
	    		frequency_thr = atoi(optarg);
	    		break;
	    	case 'l':
	    		subgraph_THR = atof(optarg);
	    		break;
			case 'u':
				directed = false;
	    		break;
			
			case '?':
	    		print_usage (stdout, 1);
				
			case -1:		/* Done with options. */
			    break;
			
			default:		/* Something else: unexpected. */
                print_usage (stderr, -1);
		}
    } while (next_option != -1);
    cout<<"subgraph_THR = " << subgraph_THR<<endl;
    //directed = false;//Tao 2016 Jun 25, only work on undirected graph now.

	if (input_filename == NULL) {
		fprintf(stderr, "Input Argument Error: Please specify an input filename using \"-i <FILE NAME>\".\n");
        print_usage (stderr, -1);
	}
	
	if (subgraphSize == -1) {
		fprintf(stderr, "Input Argument Error: Please specify a motif size using \"-s <MOTIF SIZE>\".\n");
        print_usage (stderr, -1);
	}
	
	printf("Motif Size: %d\n", subgraphSize);
	printf("Input Graph: %s\n", input_filename);
	
	if (!ReadData(input_filename))
		return 0;

	printf("ReadData Finished!!!\n");
	cout<<"Nodes:"<<g->nV<<endl;
	cout<<"Edges:"<<g->Edges()<<endl;
	//g->setPath(output_directory);

	clock_t startTime = clock();
	//for main graph
	isRand = false;
	g->subgraphCounter = 0;
	g->notClassSubCounter = 0;
	clock_t mainStartTime = clock();
	printf("before enumerate, everything is OK!\n");
	Enumerate();

	printf("after enumerate, everything is OK!!\n");
	clock_t mainEndTime = clock();

	hash_map<string, pair< vector<const string*>, long long int> > degreeSeqPair;
	hash_map<string, pair< vector<const string*>, long long int> >::iterator itor;
	cout<<"graphIntSize = "<<graphInt.size()<<endl;
	hash_map<string, string> adjStr2degSeq;
	hash_map<string, string>::iterator aditer; 
	// string graphIntOut_filename = "graphInt_size_"+std::to_string(subgraphSize)+".txt";
	// ofstream graphIntOut;
	// graphIntOut.open(graphIntOut_filename);
	// for(hash_map<std::string, long long int>::iterator it = graphInt.begin(); it!= graphInt.end(); it++)
	// {
	// 	graphIntOut<<it->second<<endl;
	// }
	// graphIntOut.close();
	for(hash_map<std::string, long long int>::iterator it = graphInt.begin(); it!= graphInt.end(); it++)
	{
		int index = 0;
		string adj = it->first;
		//cout<<"adj: "<<adj<<endl;
	    std::vector< pair<int, int> > nodeDegreeVec;
	    for (i = 0; i < subgraphSize; ++i)
	    {
	    	nodeDegreeVec.push_back(make_pair(i, 0));
	    }
	    for(i = 0; i < subgraphSize; i++)
	    {
	        for(j = i+1; j<subgraphSize; j++)
	        {
	        	index = j*(j-1)/2+i;
	            if(adj[index] == '2') //Tao modified on Feb. 19 2017
	            {
	                nodeDegreeVec[i].second += 1;
	                nodeDegreeVec[j].second += 1;
	            }
	            else if(adj[index] == '1')// i point to j, but j doesn't point to i
	            {
	            	nodeDegreeVec[i].second += 1; //add 1 to outdegree of i
	            }
	            else if(adj[index] == '3') // i doesn't point to j, but j point to i
	            {
	            	nodeDegreeVec[j].second += 1;
	            }
	        }
	    }
	    string reorderedAdj = "";

	  
   
	    sort(nodeDegreeVec.begin(), nodeDegreeVec.end(), ordering());
	
		

	    if(it->second > 0)
	    {
		    reorderedAdj.reserve(subgraphSize*subgraphSize/2);

		    for (int a = 1; a < subgraphSize; ++a)
		    {
		    	for(int b = 0; b < a; ++b)
		    	{
		    		i = nodeDegreeVec[a].first;
		    		j = nodeDegreeVec[b].first;
		    		if(i>j)
		    		{
		    			int temp = i;
		    			i = j;
		    			j = temp;
		    			index = j*(j-1)/2+i;
		    			if(adj[index] == '0') // Tao modified on Feb 19 2017
			    		{
			    			reorderedAdj.push_back('0'); 
			    		}
			    		else if(adj[index] == '2')
			    		{
			    			reorderedAdj.push_back('2');
			    		}
			    		else if(adj[index] == '1')
			    		{
			    			reorderedAdj.push_back('1');
			    		}
			    		else if(adj[index] == '3')
			    		{
			    			reorderedAdj.push_back('3');
			    		}
			    		else
		    				cout<<"NOT possible! "<< adj[index]<<endl;

		    		}
		    		else
		    		{
		    			index = j*(j-1)/2+i;
		    			if(adj[index] == '0') // Tao modified on Feb 19 2017
			    		{
			    			reorderedAdj.push_back('0'); 
			    		}
			    		else if(adj[index] == '2')
			    		{
			    			reorderedAdj.push_back('2');
			    		}
			    		else if(adj[index] == '1')
			    		{
			    			reorderedAdj.push_back('3');
			    		}
			    		else if(adj[index] == '3')
			    		{
			    			reorderedAdj.push_back('1');
			    		}
			    		else
		    				cout<<"NOT possible! "<< adj[index]<<endl;

		    		}
		    		
		    		//cout<<"a:"<<a<<" b:"<<b<<" i:"<<i<<" j:"<<j<<" index:"<<index<<" adj[index]:"<<adj[index]<<endl;
		    	}
		    }
		}// end if
		else
		{
			cout<<"no way"<<endl;
			//reorderedAdj = adj;
		}
		
		//reorderedAdj = adj;
		  

		
	    //cout<<"reorderedAdj: "<< reorderedAdj<<endl;
	    hash_map<std::string, long long int>::iterator reitor = reordered_graphInt.find(reorderedAdj);
	    string graphDegSeq = "";

	    long long int tempCount = it->second;
	 
	    if (reitor == reordered_graphInt.end())
	    {

	    	//cout<<reorderedAdj<<endl;

	    	reordered_graphInt[reorderedAdj] = tempCount;

	    	// reitor = reordered_graphInt.find(reorderedAdj);
	    	// const std::string *graphAdjMatStr = &(reitor->first);
	    	
	    	//graphDegSeq.reserve(subgraphSize);
	    
		    for(std::vector<pair<int, int> >::iterator it3 = nodeDegreeVec.begin(); it3 != nodeDegreeVec.end(); it3++)
		    {
		    	int deg = it3->second;
		    	
		        graphDegSeq.push_back(IntToAscii.find(deg-1)->second); 
		    }
		    adjStr2degSeq[reorderedAdj] = graphDegSeq;

	    }
	    else
	    {
	    	reordered_graphInt[reorderedAdj] += tempCount;
	    }

	}
	clock_t graphReorderTime = clock();
	cout<<"reordered_graphIntSize = "<<reordered_graphInt.size()<<endl;
	int tmp_frequency_thr = 0;
	for(hash_map<std::string, long long int>::iterator it = reordered_graphInt.begin(); it!= reordered_graphInt.end(); it++)
	{
		
		// add by Tao 2016-10-23, set the frequency_thr to be the size of biggest AdjString bin
		if(tmp_frequency_thr < it->second)
		{
			tmp_frequency_thr = it->second;
		}
		const std::string *graphAdjMatStr = &(it->first);
		std::string graphDegSeq = adjStr2degSeq[*graphAdjMatStr];
		long long int tempCount = it->second;
		itor = degreeSeqPair.find(graphDegSeq);
		if( itor == degreeSeqPair.end())
		{
			vector<const string*> tmpstringvec;
			tmpstringvec.push_back(graphAdjMatStr);
			pair< vector<const string*>, long long int > tmpPair = make_pair(tmpstringvec, tempCount);

			degreeSeqPair[graphDegSeq] = tmpPair;
		}
		else
		{
			(itor->second).first.push_back(graphAdjMatStr);
			(itor->second).second += tempCount;
		}
	}

	clock_t CalDegreeSeqTime = clock();
	//frequency_thr = 0;
	if(frequency_thr == -1)
	{
		frequency_thr = tmp_frequency_thr;
	}
	cout<<"frequency_thr = "<< frequency_thr<<endl;
	cout<<"degreeSeqPair.size = "<< degreeSeqPair.size()<<endl;
	double TreeCam_time = 0, GraphCam_time = 0;
	hash_map<string, long long int> finalGraph;
	for(itor = degreeSeqPair.begin(); itor != degreeSeqPair.end(); itor++)
	{
		if((itor->second).second < frequency_thr )
		{
			continue;
		}
		//cout << itor->first<<endl;
		//add Tree part here!!! 2016-10-24
		if(isTree(itor->first, subgraphSize) && !directed)
		{
			clock_t tmpstart = clock();
			//cout<< "Tree: " << itor->first<<"  "<< (itor->second).second<<endl;
			for(int iv = 0; iv < (itor->second).first.size(); iv++ )
			{
				string adjStr = *(((itor->second).first)[iv]);
				//cout<<"before GetTreeCanStr OK"<<endl;
				//cout<<adjStr<<endl;
				string treeStr = GetTreeCanStr(adjStr, subgraphSize);//!!!!Here !!!! Need the new function
				//cout<<"treeStr"<<treeStr<<endl;	
				hash_map<std::string, long long int>::iterator iter = treeInt.find(treeStr);
				if(iter == treeInt.end()) //BUG!! graphInt.end()->treeInt.end()
				{
					treeInt[treeStr] = reordered_graphInt.find(adjStr)->second;
				}
				else
				{
					(iter->second) += reordered_graphInt.find(adjStr)->second;
				}
			}
			clock_t tmpend = clock();
			TreeCam_time += (difftime(tmpend, tmpstart)/(double)CLOCKS_PER_SEC);
		}
		else
		{
			clock_t tmpstart = clock();
			for(int iv = 0; iv < (itor->second).first.size(); iv++ )
			{
				//cout<<"here"<<endl;
				string tempCam = *(((itor->second).first)[iv]);
				string cam = calculateCam(tempCam, subgraphSize);
				callNautyCount += 1;

				hash_map<std::string, long long int>::iterator it2 =  finalGraph.find(cam);
	            if (it2 == finalGraph.end())
	            {

	                finalGraph[cam] = reordered_graphInt.find(tempCam)->second;
	            }
	            else
	            {
	                (it2->second) += reordered_graphInt.find(tempCam)->second;
	            }
			}
			clock_t tmpend = clock();
			GraphCam_time += (difftime(tmpend, tmpstart)/(double)CLOCKS_PER_SEC);
		}
	}

	ofstream frequentCams("frequentCams.txt");
	//count graph part
	int graph_class = 0, tree_class = 0;
	for(hash_map<std::string, long long int>::iterator it = finalGraph.begin(); it!= finalGraph.end(); it++)
	{
		if (it->second >= frequency_thr)
		{
			graph_class += 1;
			frequentCams<<it->first<<"\t"<<it->second<<endl;
		}
	}	

	//count tree part
	hash_map<string, long long int> finalTree;
	for(hash_map<std::string, long long int>::iterator it = treeInt.begin(); it!= treeInt.end(); it++)
	{
		if (it->second >= frequency_thr)
		{
			tree_class += 1;
			//cout<< it->first<<"  "<< it->second<<endl;
		}
	}


	subgraphCounterMain = g->subgraphCounter;
	enumerated_class = graph_class + tree_class;

	frequentCams.close();
	delete g;
	//delete IsD;
	delete subgraphDegree;
	delete graphDegs;
	//cout<<"1"<<endl;
	clock_t endTime = clock();
	main_time = difftime(mainEndTime, mainStartTime)/(double)CLOCKS_PER_SEC;
	double total_time = difftime(endTime, startTime)/(double)CLOCKS_PER_SEC;
	double graphReorder_time = difftime(graphReorderTime, mainEndTime)/(double)CLOCKS_PER_SEC;
	double degreeSeq_time = difftime(CalDegreeSeqTime, graphReorderTime)/(double)CLOCKS_PER_SEC;


	printf("\n===========RESULTS===========\n");
	
	printf("\nMotif Size: %d\n", subgraphSize);
	printf("\nTotal number of subgraphs in original network: %llu\n", subgraphCounterMain);
	printf("Number of frequent non-isomorphic classes: %lu\n", enumerated_class);
	printf("Number of frequent non-isomorphic graph classes: %lu\n", graph_class);
	printf("Number of frequent non-isomorphic tree classes: %lu\n", tree_class);
	printf("\nTime Used for Enumerate:      %f\n", main_time); 
	
	printf("\nTotal Time Used: %f\n", total_time); 
	printf("\ngraphReorder_time: %f\n", graphReorder_time);
	printf("\ndegreeSeq_time: %f\n", degreeSeq_time);
	printf("\nTreeCam_time: %f\n", TreeCam_time);
	printf("\nGraphCam_time: %f\n", GraphCam_time);
	printf("\n=============================\n");
	printf("Call Nauty %lld\n", callNautyCount );
	
	//printf("Time for density check: %f\n", densityCheckTime);
	//printf("Time for adjStr: %f\n", AdjStrTime);
	return 0;
}
