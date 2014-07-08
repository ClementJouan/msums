#include <cstdlib>

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "spmain.h"
#include "spinout.h"
#include "sputil.h"
#include "sidna.h"
#include "options.h"

using namespace std;


template<typename T>
void splitStr(const string & str, int pos, T & x1, T & x2)
	{
	x1 = lexical_cast<T>(str.substr(0, pos));
	x2 = lexical_cast<T>(str.substr(pos+1, str.size()-pos-1));
	}

void selectStats(const vector<vector<string> > & sel, 
	vector<bool> & sstats, vector<bool> & pstats, StatDict & dict);
void selectPops(const vector<vector<string> > & sel,
	vector<bool> & pops, vector<vector<bool> > & pairs);

void dependencies(vector<bool> & compSingles, vector<vector<bool> > & pairs,
	vector<bool> & sstats, vector<bool> & pstats);

int main(int argc,char *argv[])
	{
	string inputfilename = "spinput.txt";
	string statfilename = "ABCstat.txt";
	bool printPerLocus = false;
	bool printHelp = false;
	bool listStats = false;
	
	typedef OpList::Collector OLIter;
	OpList popList, statList;
	OLIter keepS(statList, "+"), dropS(statList, "-"), 
		   keepP(popList, "+"), dropP(popList, "-");

	OptionParser p;
	UnaryOption<string> o_ifn(p, 'i', "init", inputfilename);
	UnaryOption<string> o_sfn(p, 'o', "output", statfilename);
	SwitchOption o_ppl(p, 'l', "print_per_locus", printPerLocus);
	SwitchOption o_hlp(p, 'h', "help", printHelp);
	SwitchOption o_list(p, 'a', "list_stats", listStats);
	NaryOption<string, OLIter> o_ds(p, 's', "dropStats", dropS);
	NaryOption<string, OLIter> o_ks(p, 'S', "keepStats", keepS);
	NaryOption<string, OLIter> o_dp(p, 'p', "dropPops", dropP);
	NaryOption<string, OLIter> o_kp(p, 'P', "keepPops", keepP);
	
	try {
	p.parse(argv+1, argc-1);
	} catch (ArgException & e)
		{
		error(string("Error in main: unknown command line argument: ") +
		   	e.what());
		}

	if (printHelp)
		{
		cout <<
			"analms [-i FILE] [-o FILE] [-h] [-l] [[-s|-S STATLIST]...] "
			"[[-p|-P POPLIST]...]\n\n";
		cout <<
"-i, --init FILE       Use this initialization file. [spinput.txt]\n"
"-o, --output FILE     Use this output file. [ABCstat.txt]\n"
"-h, --help            Print the help text and exit.\n"
"-a, --list_stats      Print list of available stats and exit.\n"
"-l, --print_per_locus Print stats per locus, don't print aggregate values.\n"
"                      [false]\n"
"-S, --keepStats LIST  Print stats in LIST. LIST is a space-separated list of\n"
"                      stat names. 'all' selects all stats.\n"
"-S, --dropStats LIST  Do *not* print stats in LIST (see above).\n"
"-P, --keepPops LIST   Show populations or pairs in LIST. Populations can be\n"
"                      specified as single number (e.g. '2'), or range \n"
"                      (e.g. '1-5'). Pairs are specifies as p1xp2 (e.g. '1x2').\n"
"                      'all' selects all populations, 'allxall' all pairs.\n"
"-p, --dropPops LIST   Do *not* show populations/pairs in LIST (see above).\n";
		exit(0);
		}

	if (listStats)
		{
		cout << "sorry, not implemented yet :-(\n";

		exit(0);
		}


// **** starts initializations by getting initial conditions

	string datafilename;

	/** Number of sequences per population x locus.*/
	vector<vector<int> > n_sequences;
	/** Number of sites per locus. Assumes all sequences per locus are the
	    same length. */
	vector<int> n_sites;
	int n_datasets;

	try {
	get_initial_conditions(inputfilename, 
		n_sequences, n_sites, n_datasets, datafilename);
	// dump initial conditions 
	print_initial_conditions(cout, 
		inputfilename, n_sequences, n_sites, n_datasets, datafilename); 
	} catch (exception & e) 
		{error(e.what());}

	const int n_pops = n_sequences.size();
	const int n_loci = n_sites.size();


// **** now we know #pops we can process command line arguments

	vector<bool> showPops(n_pops, true);
	vector<vector<bool> > showPairs;
	showPairs.reserve(n_pops);
	for (int p=0; p<n_pops; p++)
		showPairs.push_back(vector<bool>(n_pops, true));
	// process command line selection of pops/pairs
	selectPops(popList.getList(), showPops, showPairs);

	StatDict statDict;
	vector<bool> showSStats(_sstats_last, true), showPStats(_pstats_last, true);

	// process command line selection of stats
	// and convert to simple list
	selectStats(statList.getList(), showSStats, showPStats, statDict);

	vector<bool> compSingles = showPops;
	// make sure only stats that are required are computed
	dependencies(compSingles, showPairs, showSStats, showPStats);


// **** prepare header of output file

	ofstream abcstatfile(statfilename.c_str());
	if (!abcstatfile.good())
		error("Error in main: cannot open file " + statfilename);

	try {
	initialize_write_ABCstat(abcstatfile, 
		showPops, showPairs, showSStats, showPStats, statDict, !printPerLocus);
	} catch (exception & e)
		{error(e.what());}

// **** prepare data containers

	/** Sequence by species x locus x sample. */
	vector<vector<Sample> > sequences;
	// pre-allocate for efficiency
	sequences.reserve(n_pops);

	int maxNSeq = 0;

	// building array for sequence data
	for (int p=0; p<n_pops; p++)
		{
		sequences.push_back(vector<Sample>());
		sequences.reserve(n_loci);
		for (int l=0; l<n_loci; l++)
			{
			sequences[p].push_back(Sample());
			sequences[p][l].reserve(n_sequences[p][l]);
			for (int s=0; s<n_sequences[p][l]; s++)
				// fill with empty strings
				sequences[p][l].push_back("");

			maxNSeq = max(maxNSeq, n_sequences[p][l]);
			}
		}

	/** Outgroup data. */
	Sample outgroup(n_loci);
	/** Vector with number of polymorphic sites per locus from 0 to n_loci-1. */
	vector<int> n_polym_sites(n_loci);
	
	ifstream input(datafilename.c_str());
	if (!input.good())
		error("Error in reading data set file: cannot open file " +
			datafilename);

// **** prepare results objects

	vector<SingleStats> pop_results;
	pop_results.reserve(n_pops);
	for (int p=0; p<n_pops; p++)
		pop_results.push_back(SingleStats(n_loci));

	PairStats pair_results(n_loci);


// **** start reading datasets and calculating
// **** dataset = replicate, i.e. 1 dataset = loci x pops x samples

	for(int dataset=0; dataset<n_datasets; ++dataset) 
		{ 
		const int imain = n_datasets/100;
		if (imain>0 && (dataset % imain)==0) 
			{
			string t_str;
			get_time_short(t_str);
			cerr << "\nReplicate: " << dataset << " Time: " << t_str;
			}

		try {
		get_dataset(input, dataset, 
			n_sequences, n_sites, sequences, outgroup, n_polym_sites);
		} catch (exception & e)
			{error(e.what());}

		// new data set
		abcstatfile << "\n"; 

		if (printPerLocus)
			analyse(abcstatfile, sequences, stats_per_locus);
		else
			analyse(abcstatfile, sequences, stats_aggregate);
		}

	return 0;
	}    /*end Main*/

void analyse(ostream & abcstatfile, const vector<vector<Sample> > & sequences,
	vector<AnalysisBase *> single_stats)
	{
	const size_t n_loci = sequences[0].size();
	const size_t n_pops = sequences.size();

	for(size_t l=0; l<n_loci; l++)	
		{
		// new locus
		abcstatfile << "\n" << dataset << "\t" << l;

		for (size_t p=0; p<n_pops; p++)
			{
			if (!compSingles[p])
				continue;

			for (int stat=0; stat<single_stats.size(); stat++)
				abcstatfile << '\t' <<
					single_stats[stat]->analyse(sequences[p][l]);
			}
		
		for (size_t pop1=0; pop1<n_pops-1; pop1++)
			for (size_t pop2=pop1+1; pop2<n_pops; pop2++)
				{
				if (!showPairs[pop1][pop2])
					continue;

				const Sample p1 = sequences[pop1][l];
				const Sample p2 = sequences[pop2][l];
				for (size_t stat=0; stat<pair_stats.size(); stat++)
					{
					pair_stats[stat]->analyse(p1, p2);

					abcstatfile << '\t' << pair_stats[stat].mean()
						<< '\t' << pair_stats.var();
					}
				}
		}
	}


void analyse_aggr(ostream & abcstatfile, 
	const vector<vector<Sample> > & sequences,
	vector<AnalysisBase *> single_stats)
	{
	const size_t n_pops = sequences.size();

	abcstatfile << dataset;

	for (size_t p=0; p<n_pops; p++)
		{
		for (size_t stat=0; stat<single_stats.size(); stat++)
			{
			Aggregate aggregate;
			for (size_t l=0; l<sequences[p].size(); l++)
				agreggate(single_stats[stat]->analyse(sequences[p]));

			aggregate.analyse();
			abcstatfile << '\t' << aggregate.mean() << '\t'<< aggregate.var();
			}
		}

	vector<Aggregate> aggregates();

	for (size_t pop1=0; pop1<n_pops-1; pop1++)
		for (size_t pop2=pop1+1; pop2<n_pops; pop2++)
			{
			if (!showPairs[pop1][pop2])
				continue;

			const vector<Sample> p1 = sequences[pop1];
			const vector<Sample> p2 = sequences[pop2];

			aggregates.clear();
			aggregates.resize(pair_stats.size());

			for (size_t l=0; l<p1.size(); l++)
				{
				// PairSample caches analysis results
				PairSample ps(p1[l], p2[l]);
				for (size_t stat=0; stat<pair_stats.size(); stat++)
					aggregate[stat](pair_stats[stat]->analyse(ps));
				}

			for (size_t stat=0; stat<pair_stats.size(); stat++)
				{
				aggregate[stat].analyse();
				abcstatfile << '\t' << aggregate[stat].mean() <<
				   	'\t' << aggregate[stat].var();
				}
			}
	}

void selectStats(
	const vector<vector<string> > & sel, 
	const vector<string> & all, set<string> & use)
	{
	for (size_t s=0; s<sel.size(); s++)
		{
		cerr << "sS: " << sel[s][0];
		const bool to = sel[s][0][0] == '+';
		for (size_t t=1; t<sel[s].size(); t++)
			{
			cerr << ", " << sel[s][t];
			if (sel[s][t] == "all")
				{
				for (size_t s=0; s<all.size(); s++)
					if (to)
						use.insert(all[s]);
					else
						use.erase(all[s]);
				continue;
				}

				error(
					"Error in parameters: stat '" + sel[s][t] + "' unknown!");

			if (to)
				use.insert(sel[s][t]);
			else
				use.erase(sel[s][t]);
			}
		cerr << endl;
		}
	}

void checkPopIndex(size_t i, size_t n)
	{
	if (i >= n)
		error(string("Error in parameter: population index ") + 
			lexical_cast<string>(i) + " out of range");
	}

void selectPops(const vector<vector<string> > & sel,
	vector<bool> & pops, vector<vector<bool> > & pairs)
	{
	// go through all selection sequences
	for (size_t s=0; s<sel.size(); s++)
		{
		// prefix determines whether this adds or removes pops/pairs
		const bool to = sel[s][0][0] == '+';
		
		// go through list of operations (skipping prefix)
		for (size_t t=1; t<sel[s].size(); t++)
			{
//			cerr << "-pP " << sel[s][t] << endl;
			const string & str = sel[s][t];
			// all pops
			if (str == "all")
				{
				for (size_t p=0; p<pops.size(); p++)
					pops[p] = to;
				continue;
				}
			// all pairs
			if (str == "all" || str == "allxall")
				{
				for (size_t p=0; p<pairs.size()-1; p++)
					for (size_t p2=p+1; p2<pairs[p].size(); p2++)
						pairs[p][p2] = to;
				continue;
				}

			size_t c;
			for (c=0; c<str.size(); c++)
				{
				if (str[c] == '-')	// pop range
					{
					size_t p1, p2;

					try {
					splitStr(str, c, p1, p2);
					} catch(exception & e)
						{error(string("Error in parameter: '") + str + 
							"' is not a valid population range");}

					checkPopIndex(p1, pops.size());
					checkPopIndex(p2, pops.size());

					for (size_t p=p1; p<=p2; p++)
						{
						pops[p] = to;
						for (size_t pp=0; pp<pops.size(); pp++)
							pairs[min(p, pp)][max(p, pp)] = pops[p] && pops[pp];
						}
					break;
					}
				if (str[c] == 'x') // pair
					{
					int p1, p2;

					try {
					splitStr(str, c, p1, p2);
					} catch(exception & e)
						{error(string("Error in parameter: '") + str + 
							"' is not a valid pair");}

					checkPopIndex(p1, pops.size());
					checkPopIndex(p2, pops.size());

					pairs[min(p1,p2)][max(p1,p2)] = to;
					break;
					}
				}

			if (c == str.size()) // no '-' or 'x' found => single pop
				{
				size_t p;
				try {
				p = lexical_cast<size_t>(str);
				} catch(exception & e)
					{error(string("Error in parameter: '") + str + 
						"' is not a valid population or pair");}

				checkPopIndex(p, pops.size());

				pops[p] = to;
				for (size_t pp=0; pp<pops.size(); pp++)
					pairs[min(p, pp)][max(p, pp)] = pops[p] && pops[pp];
				}
			}
		}
	}



void dependencies(vector<bool> & compSingles, vector<vector<bool> > & pairs,
	vector<bool> & sstats, vector<bool> & pstats)
	{
	// check whether any pops are supposed to be shown at all
	bool popsShown = false;
	for (size_t i=0; i<compSingles.size(); i++)
		if ((popsShown = popsShown || compSingles[i]))
			break;
	// if not, single stats can be dropped
	if (! popsShown)
		for (size_t s=0; s<sstats.size(); s++)
			sstats[s] = false;


	// check whether any pairs are supposed to be shown at all
	bool pairsShown = false;
	for (size_t p=0; p<pairs.size()-1; p++)
		for (size_t pp=p+1; pp<pairs.size(); pp++)
			if ((pairsShown = pairsShown || pairs[p][pp]))
				break;
	// if not, pair stats can be dropped
	if (!pairsShown)
		for (size_t s=0; s<pstats.size(); s++)
			pstats[s] = false;

	bool pairsNeedSingles = false;
	bool singlesNeeded = false;
	// check whether any of the pair stats requires comp.
	// of single stats
	for (size_t i=0; i<pstats.size(); i++)
			// show this stat
		if (pstats[i] && 
			// needs single stat results
			requiresSinglePopResults(PStatsI(i)))
			{
			pairsNeedSingles = true;
			break;
			}

	// or whether single stats will be shown
	for (size_t i=0; i<sstats.size(); i++)
		// there's a single stat that needs to be shown
		if (sstats[i])
			{
			singlesNeeded = true;
			break;
			}
	
	// no single stats to be shown, thus none to be computed
	if (!singlesNeeded)
		for (size_t p=0; p<compSingles.size(); p++)
			compSingles[p] = false;

	if (!pairsNeedSingles)
		// all fine, single pops independent of pairs
		return;

	// here all pairs imply compute of single pops
	for (size_t p1=0; p1<pairs.size()-1; p1++)
		for (size_t p2=p1+1; p2<pairs.size(); p2++)
			if (pairs[p2][p1])
				{
				compSingles[p1] = true;
				compSingles[p2] = true;
				}
	}
