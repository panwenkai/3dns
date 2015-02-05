//  _________________________________________________________
// |
// |   Time.cpp   Heat flow related routines
// |
// |
// |   (C) 1998-99  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________

#include "Clock.h"
#include "Compat.h"

// ______________________________________________________________
//  CLOCK::CLOCK()
//		default constructor
// ______________________________________________________________
CLOCK::CLOCK()
{
	memset(this, 0, sizeof(CLOCK));
 	
}; //endmethod


// ______________________________________________________________
//  CLOCK::~CLOCK()
//		default destructor
// ______________________________________________________________
CLOCK::~CLOCK()
{
	delete [] eraBound;			//coordinates of era boundaries
	delete [] stepMin;			//minimum step
	delete [] stepMax;			//maximum step
 	
}; //endmethod



//___________________________________________________
// CLOCK::CheckAlarm
//	  Compare positions of 2 clocks
//___________________________________________________
bool CLOCK::CheckAlarm(const CLOCK &clk)
{
//	if (!this->Valid())
//		return(false);

	return (clk.tCurrent >= this->tCurrent );
};


//___________________________________________________
// CLOCK::Decrement
//	  Decrement clock
//___________________________________________________
bool CLOCK::Decrement()
{
	if (tPrev == tCurrent)	//already decremented or at beginning
		return(false);

	tCurrent = tPrev;		//previous time

	if (tCurrent < eraBound[eCurrent]) //backed out of era
		eCurrent--;

	curTime = (double)(tCurrent * FS_TO_S);
	curDTime = 0;

	return(true);
}; //endif
	

//___________________________________________________
// CLOCK::Increment
//	  Increment clock
//___________________________________________________
bool CLOCK::Increment(__int64 alarm)
{
	__int64 stepNext;

	if (!Valid())			//also checks for end of clock condition
		return(false);
	
	stepNext = stepCurrent;

	if ((alarm > 0) && ((tCurrent + stepNext) > alarm)) //hit an alarm
		stepNext = (alarm - tCurrent);

	if ((tCurrent + stepNext) >= eraBound[eCurrent+1]) //beyond end of era
	{
		eCurrent++;		//begin new era
		stepNext = (eraBound[eCurrent] - tCurrent);
		stepCurrent = stepMin[eCurrent];
	};

	if (stepNext == 0)		//clock MUST advance
		return(false);	

	tPrev = tCurrent;		//save old time
	tCurrent += stepNext;	//move ahead
	stepPrev = stepNext;		//remember size of previous step
	
	curTime = (double)(tCurrent * FS_TO_S);
	curDTime = (double)(stepPrev * FS_TO_S);
	
	if (eCurrent >= eMax)
		return(false);

	return(true);

}; //endmethod

bool CLOCK::Increment(CLOCK &alarmClock)
{
	return( this->Increment(alarmClock.tCurrent));
};


//___________________________________________________
// CLOCK::Reset
//	  Initialize clock parameters and set to zero
//___________________________________________________
bool CLOCK::Reset(unsigned int numEras, double *eLength, double *stepLo, double *stepHi)
{
	unsigned int i;

	if (numEras >= MAX_ERAS)
		return(false);

	eMax = numEras;
	
	eraBound = new __int64[numEras+1];		//coordinates of era boundaries

	stepMin =		new __int64[numEras+1];		//minimum step size
	stepMax =		new __int64[numEras+1];		//maximum step size	
	
	eraBound[0] = 0;

	for (i = 0; i<numEras; i++)
	{
		eraBound[i+1] = eraBound[i] + (__int64)(eLength[i] * S_TO_FS);

		stepMin[i] = (__int64)((double)((stepLo[i] < stepHi[i]) ? stepLo[i] : stepHi[i]) * S_TO_FS);
		stepMax[i] = (__int64)((double)((stepLo[i] > stepHi[i]) ? stepLo[i] : stepHi[i]) * S_TO_FS);

		if ((stepMin[i] <= 0) || (stepMax[i] <= 0))
			return(false);			//must be positive
	};

	this->Reset();

	return(true);
}; //endmethod

bool CLOCK::Reset()		//reset times only
{
	eCurrent = 0;		//current era
	
	tPrev = 0;			//previous time
	stepPrev = 0;		//previous timestep

	tCurrent = 0;					//current time
	stepCurrent = stepMin[0];		//begin at minimum (default) clock

	curTime = 0.0;
	curDTime = 0.0;

	return(true);
}; //endmethod


//___________________________________________________
// CLOCK::StepAdjust
//	Adjust time step if possible
//___________________________________________________
bool CLOCK::StepAdjust(double factor)
{
	__int64 newStep;

	if (!Valid() || (factor < 0.0))
		return(false);		//only positive factors

	else if (factor == 0.0)		//reset to default
		newStep = stepMin[eCurrent];

	else 			//calculate new step	
		newStep = ( (__int64)((stepPrev > 0) ? stepPrev : stepCurrent) * 
				(__int64)(factor * 10000)) / ((__int64)10000);

	if (newStep > stepMax[eCurrent])   //upper bound
		newStep = stepMax[eCurrent];

	else if (newStep < (__int64)(MIN_CLOCK * S_TO_FS))	//lower bound
		newStep = (__int64)(MIN_CLOCK * S_TO_FS);

		stepCurrent = newStep;
											//update to new step if larger
	return(true);
}; //endfunc


//___________________________________________________
// CLOCK::Valid
//	  Check if clock is valid
//___________________________________________________
bool CLOCK::Valid()
{
	return ((eMax > 0) && (tCurrent < eraBound[eMax]));
}; //endmethod


//_________________________________________
// CLOCK::operator==
//	test for equality of strings
//_________________________________________
bool CLOCK::operator==(const CLOCK &rhs)
{
	return (this->tCurrent == rhs.tCurrent);
};


