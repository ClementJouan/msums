#ifndef maskeFILE_H
#define maskeFILE_H

#include "sputil.h"

void read_maske(istream & inp, vector<vector<bool> > & maske)
	{
	string str;
	
	while(true)
		{
		skip_space(inp, str);
		
		if (! inp.good())
			break;

		maske.push_back(vector<bool>());
		maske.back().resize(str.size());

		for (size_t i=0; i<str.size(); i++)
			maske.back()[i] = str[i] == '1';
		}
	}


#endif	// maskeFILE_H
