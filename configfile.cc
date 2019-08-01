#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "sputil.h"
#include "configfile.h"

using namespace std;

typedef SPIOException SPIOE;

void ConfigFile::read(istream & inp)
 	{
	size_t n_loci, npops;

	if (!get_value(inp, n_loci))
		throw SPIOE(ERR_LOC
			" Error in reading initial conditions: cannot read nloci");
	if (!get_value(inp, npops))
	   throw SPIOE(ERR_LOC
			" Error in reading initial conditions: cannot read npops");

	for(size_t p=0; p<npops; p++)
		{
		_n_sequences.push_back(vector<size_t>());
		_n_sequences.back().resize(n_loci);
		}

	_n_sites.resize(n_loci);

	for (size_t l=0; l<n_loci; l++)
		{
		for (size_t  p=0; p<npops; p++)
			{
			int v;
			if (!get_value<int>(inp, v))
				{
				throw SPIOE(ERR_LOC
					" Error in reading initial conditions: cannot read nSeq at locus "
						+ lexical_cast<string>(l));
				}
			_n_sequences[p][l] = v;
			}

		if (!get_value(inp, _n_sites[l]))
			{
			throw SPIOE(ERR_LOC
				" Error in reading initial conditions: cannot read nSites at locus "
					+ lexical_cast<string>(l));
			}
		}

	if (!get_value(inp, _n_datasets))
		throw SPIOE(ERR_LOC
			" Error in reading inital conditions: cannot read ndatasets");
	if (!get_value(inp, _datafilename))
	   throw SPIOE(ERR_LOC	
			" Error in reading initial conditions: cannot read data file name");
	} 


/* print initial values of parameters read from input file as a check*/
void ConfigFile::dump(ostream & out)
	{
	const size_t nLoci = _n_sites.size();

	out << "\n\tnumber of populations: " << _n_sequences.size();
	out << "\n\tnumber of loci: " << nLoci;
	for(size_t l=0;l<nLoci;l++)	/*loop over loci*/
		{
		for (size_t p=0; p<_n_sequences.size(); p++){
			if (_n_sequences[p][l] == 0){
				throw SPIOE(ERR_LOC
					"Error, number of sequence cannot be 0 at locus " + lexical_cast<string>(l));
			}

			out << "\n\tnumber of sequences in population " << p << 
				" at locus " << l << ":" << _n_sequences[p][l];
		}

		if (_n_sites[l] == 0){
			throw SPIOE(ERR_LOC
				"Error, number of sites cannot be 0 at locus " + lexical_cast<string>(l));
		out << "\n\t\tnumber of sites at locus " << l << ": " << _n_sites[l];
		}
		//cout << _n_sites[l] << "\n";
	}

	out << "\n\tnumber of replicate datasets: " << _n_datasets; 
	out << "\n\tname of the dataset file: " << _datafilename;
	out << "\n\n";
	}	  /*end of print_initial_conditions*/

