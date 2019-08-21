/*
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH
THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#ifndef _THREE_D_PATTERN_H_
#define _THREE_D_PATTERN_H_


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <random>
#include <time.h>

#include "OTU_Profile.h"

using namespace std;

class Three_D_Pattern
{

public: 

	unsigned int OTU_Index_1;	// OTU index in the OTU_Table object
	unsigned int OTU_Index_2;
	unsigned int OTU_Index_3;

	OTU_Profile * OTU_p_1;
	OTU_Profile * OTU_p_2;
	OTU_Profile * OTU_p_3;

	double score;					// real (observed) value
	double population_threshold;	// real (observed) value

	double threshold_x1;
	double threshold_x2;
	double threshold_x3;

	double p_000;
	double p_001;
	double p_010;
	double p_011;
	double p_100;
	double p_101;
	double p_110;
	double p_111;


// -------------------------------

	Three_D_Pattern();	// creates EMPTY pattern. All zeros
	~Three_D_Pattern();
};
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------

Three_D_Pattern::Three_D_Pattern()
{
	OTU_Index_1 = 0;
	OTU_Index_2 = 0;
	OTU_Index_3 = 0;

	OTU_p_1 = NULL;
	OTU_p_2 = NULL;
	OTU_p_3 = NULL;

	score=0.;
	population_threshold=0.;

	threshold_x1 = 0.;
	threshold_x2 = 0.;
	threshold_x3 = 0.;
}
//-------------------------------------------------------------------------------------------------------------

Three_D_Pattern::~Three_D_Pattern()
{
	if (OTU_p_1 != NULL) delete OTU_p_1;
	if (OTU_p_2 != NULL) delete OTU_p_2;
	if (OTU_p_3 != NULL) delete OTU_p_3;
}
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
#endif //Three_D_Pattern
