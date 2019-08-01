#ifndef RNDmin_h
#define RNDmin_h

#include "anafunctors.h"
#include "stats_multi.h"
#include "pairsample.h"

template<class SAMPLEV>
struct RNDmin : public AnalysisBase<SAMPLEV>
	{
	double analyse(const SAMPLEV & samples) const
		{
		VERIFY_MSG(samples.size() >= 3, 
			"Error: RNDmin requires at least 3 populations");


		int n_size_0 = samples[0]->size(); 
		int n_size_1 = samples[1]->size(); 
		int n_size_2 = samples[2]->size(); 
		int n_sites_0 = samples[0]->n_sites();
		int n_sites_1 = samples[1]->n_sites();
		int n_sites_2 = samples[2]->n_sites();

		double Dmin = (sum_pair_mini(
			samples[0]->alleles().begin(), samples[0]->alleles().begin() + std::max(n_sites_0,n_sites_1),
			samples[1]->alleles().begin(), samples[1]->alleles().begin() + std::max(n_sites_0,n_sites_1)))/(n_size_0*n_size_1);
		
		double dXO = (sum_pair_diff(
			samples[0]->alleles().begin(), samples[0]->alleles().begin() + std::max(n_sites_0,n_sites_2),
			samples[2]->alleles().begin(), samples[2]->alleles().begin() + std::max(n_sites_0,n_sites_2)))/(n_size_0*n_size_2);

		double dYO = (sum_pair_diff(
			samples[1]->alleles().begin(), samples[1]->alleles().begin() + std::max(n_sites_1,n_sites_2),
			samples[2]->alleles().begin(), samples[2]->alleles().begin() + std::max(n_sites_2,n_sites_2)))/(n_size_1*n_size_2);

		double RND_min = Dmin/ ((dXO + dYO)/2);
		if (Dmin == 0 & dYO == 0 & dXO == 0){
			RND_min = 0;
			}
		if (dYO == 0 & dXO == 0){
			RND_min = 0;
			}
		//cout << "Dmin: " << Dmin << " " << "dYO: " << dYO << " " << "dXO: " << dXO <<  " " << "RNDmin: " << RND_min << " " << "n_sites: " << n_sites_1 << "\n";	
		return RND_min;

		}
	};

template<class SAMPLEV>
struct RNDminO : public AnalysisBase<SAMPLEV>
	{
	double analyse(const SAMPLEV & samples) const
		{
		VERIFY_MSG(samples.size() >= 3, 
			"Error: RNDmin requires at least 3 populations");


		int n_size_0 = samples[0]->size(); 
		int n_size_1 = samples[1]->size(); 
		int n_size_2 = 1; 
		int n_sites_0 = samples[0]->n_sites();
		int n_sites_1 = samples[1]->n_sites();
		int n_sites_2 = 0;

		double Dmin = (sum_pair_mini(
			samples[0]->alleles().begin(), samples[0]->alleles().begin() + std::max(n_sites_0,n_sites_1),
			samples[1]->alleles().begin(), samples[1]->alleles().begin() + std::max(n_sites_0,n_sites_1)))/(n_size_0*n_size_1);
		
		double dXO = (sum_pair_diff(
			samples[0]->alleles().begin(), samples[0]->alleles().begin() + std::max(n_sites_0,n_sites_2),
			samples[2]->alleles().begin(), samples[2]->alleles().begin() + std::max(n_sites_0,n_sites_2)))/(n_size_0*n_size_2);

		double dYO = (sum_pair_diff(
			samples[1]->alleles().begin(), samples[1]->alleles().begin() + std::max(n_sites_1,n_sites_2),
			samples[2]->alleles().begin(), samples[2]->alleles().begin() + std::max(n_sites_2,n_sites_2)))/(n_size_1*n_size_2);

		double RND_minO = Dmin/ ((dXO + dYO)/2);
		if (Dmin == 0 & dYO == 0 & dXO == 0){
			RND_minO = 0;
			}
		if (dYO == 0 & dXO == 0){
			RND_minO = 0;
			}
		//cout << "Dmin: " << Dmin << " " << "dYO: " << dYO << " " << "dXO: " << dXO <<  " " << "RNDmin: " << RND_minO << " " << "n_sites: " << n_sites_1 << " " << "n_size1: " << n_size_1 << "\n";	
		return RND_minO;

		}
	};









#endif
