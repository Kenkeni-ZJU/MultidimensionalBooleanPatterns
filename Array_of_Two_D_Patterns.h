/*
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH
THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#ifndef _ARRAY_OF_TWO_D_PATTERNS_H_
#define _ARRAY_OF_TWO_D_PATTERNS_H_


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <random>
#include <time.h>

#include "Two_D_Pattern.h"

using namespace std;

class Array_of_Two_D_Patterns
{

public: // ------------------- For now...

	unsigned int max_number_of_patterns;	
	unsigned int number_of_patterns;

	Two_D_Pattern ** pattern;

	Array_of_Two_D_Patterns(unsigned int _max_number_of_patterns=10000, ostream & err_msg_out = cout);	// 10,000
	~Array_of_Two_D_Patterns();

// NOTE: comparison is based only on OTU_Index_1 and OTU_Index_1 (order is specific)
	bool is_Pattern_Present(Two_D_Pattern * TD_p, ostream & err_msg_out = cout);	


// NOTE: No Pattern copy will be created! The pointer to the pattern will be added to the array.
	bool add_Pattern(Two_D_Pattern * TD_p, ostream & err_msg_out = cout);

	void show_Statistics();
	void show_Statistics(ofstream & out_f);

	void show_All();
	void show_All(ofstream & out_f);



};
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
void Array_of_Two_D_Patterns::show_Statistics()
{
	cout << "\n============== Array_of_Two_D_Patterns::show_Statistics ========== (START)\n";

	cout<<"number_of_patterns="<<number_of_patterns<<endl;

	unsigned int i;
	for (i = 0; i < number_of_patterns; i++)
	{
		cout << "OTU_index1=" << pattern[i]->OTU_Index_1 << "\t";
		cout << "OTU_index2=" << pattern[i]->OTU_Index_2 << "\t";
		cout << "score=" << pattern[i]->score << "\t";
		cout << "population threshold=" << pattern[i]->population_threshold << "\t";
		cout << "threshold_x1=" << pattern[i]->threshold_x1 << "\t";
		cout << "threshold_x2=" << pattern[i]->threshold_x2 << endl;
	}

	cout << "\n============== Array_of_Two_D_Patterns::show_Statistics ========== (END)\n";

}
//-------------------------------------------------------------------------------------------------------------

void Array_of_Two_D_Patterns::show_Statistics(ofstream & out_f)
{


	out_f << "OTU 1 ID\tOTU 1 average abundance\tOTU 1 taxonomy\tOTU 2 ID\tOTU 2 average abundance\tOTU 2 taxonomy\tobserved_population_threshold\tscore\tthr_x1\tthr_x2\tp_00\tp_01\tp_10\tp_11\tpattern type" << endl;

	unsigned int i;
	for (i = 0; i < number_of_patterns; i++)
	{
	
		out_f << pattern[i]->OTU_Index_1 <<"\t" << pattern[i]->OTU_p_1->OTUs_average << pattern[i]->OTU_p_1->taxonomy_text << "\t"
			  << pattern[i]->OTU_Index_2 <<"\t" << pattern[i]->OTU_p_2->OTUs_average << pattern[i]->OTU_p_2->taxonomy_text << "\t"
			<< pattern[i]->population_threshold << "\t" << pattern[i]->score << "\t" << pattern[i]->threshold_x1 << "\t" << pattern[i]->threshold_x2 << "\t"
			<< pattern[i]->p_00 << "\t" << pattern[i]->p_01 << "\t" << pattern[i]->p_10 << "\t" << pattern[i]->p_11 << "\t";

		if (pattern[i]->pattern_type == 0) out_f << "unidentified pattern" << endl;
		if (pattern[i]->pattern_type == 1) out_f << "co-presence" << endl;
		if (pattern[i]->pattern_type == 2) out_f << "co-exclusion" << endl;
		if (pattern[i]->pattern_type == 3) out_f << "one way relations" << endl;
		if (pattern[i]->pattern_type >= 4) out_f << "unidentified pattern 2" << endl;


	}
}
//-------------------------------------------------------------------------------------------------------------

void Array_of_Two_D_Patterns::show_All(ofstream & out_f)
{
	out_f << "number_of_patterns=" << number_of_patterns << endl;

	unsigned int i;
	unsigned int ii;
	for (i = 0; i < number_of_patterns; i++)
	{
		out_f << "\n==============================================================\n";

		out_f << "OTU 1\t taxonomy\tOTU 2\ttaxonomy\tobserved_population_threshold\tscore\tthr_x1\tthr_x2\tp_00\tp_01\tp_10\tp_11" << endl;

		out_f << pattern[i]->OTU_Index_1 << pattern[i]->OTU_p_1->taxonomy_text << "\t" << pattern[i]->OTU_Index_2 << pattern[i]->OTU_p_2->taxonomy_text
			<< "\t" << pattern[i]->population_threshold << "\t" << pattern[i]->score << "\t" << pattern[i]->threshold_x1 << "\t" << pattern[i]->threshold_x2 << "\t"
			<< pattern[i]->p_00 << "\t" << pattern[i]->p_01 << "\t" << pattern[i]->p_10 << "\t" << pattern[i]->p_11 << endl;




		out_f << "OTU 1\t";
		for (ii = 0; ii < pattern[i]->OTU_p_1->number_of_samples; ii++) out_f << pattern[i]->OTU_p_1->abundance_value[ii] << "\t";
		out_f << endl;
		out_f << "OTU 2\t";
		for (ii = 0; ii < pattern[i]->OTU_p_2->number_of_samples; ii++) out_f << pattern[i]->OTU_p_2->abundance_value[ii] << "\t";
		out_f << endl << endl;

	}
}
//-------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------

Array_of_Two_D_Patterns::Array_of_Two_D_Patterns(unsigned int _max_number_of_patterns, ostream & err_msg_out)
{
	number_of_patterns=0;
	
	max_number_of_patterns=_max_number_of_patterns;

	pattern = new Two_D_Pattern *[max_number_of_patterns];

	unsigned int i;
	for (i = 0; i < max_number_of_patterns; i++) pattern[i] = NULL;
}
//-------------------------------------------------------------------------------------------------------------

Array_of_Two_D_Patterns::~Array_of_Two_D_Patterns()
{
	unsigned int i;
	for (i = 0; i < max_number_of_patterns; i++) if(pattern[i] != NULL) delete pattern[i];
	delete [] pattern;
}
//-------------------------------------------------------------------------------------------------------------

bool Array_of_Two_D_Patterns::is_Pattern_Present(Two_D_Pattern * TD_p, ostream & err_msg_out)
{
	unsigned int i;
	for (i = 0; i < number_of_patterns; i++)
	{
		if(pattern[i]->OTU_Index_1 != TD_p->OTU_Index_1) continue;
		if(pattern[i]->OTU_Index_2 != TD_p->OTU_Index_2) continue;
		if (pattern[i]->population_threshold != TD_p->population_threshold) continue;
		if (pattern[i]->score != TD_p->score) continue;
		return true;
	}
		
	return false;
}
//-------------------------------------------------------------------------------------------------------------

bool Array_of_Two_D_Patterns::add_Pattern(Two_D_Pattern * TD_p, ostream & err_msg_out)
{
	if (number_of_patterns >= max_number_of_patterns) return false;
	
	pattern[number_of_patterns] = TD_p;
	number_of_patterns++;

	return true;
}
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
#endif //Array_of_Two_D_Patterns
