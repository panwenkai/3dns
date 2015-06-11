//  _________________________________________________________
// |
// |   Clock.h    Time related functions
// |
// |   (C) 2000  Columbia University, MSME
// |   (C) 2000-01	Harvard University, DEAS
// |__________________________________________________________
#define S_TO_NS      ((double)(1E9))
#define S_TO_FS		 ((double)(1E15))
#define S_TO_PS		 ((double)(1E12))
#define FS_TO_S		 ((double)(1E-15))
#define PS_TO_S		 ((double)(1E-12))

#define MIN_CLOCK	 ((double)1E-14)  //minimum time step
#define MAX_CLOCK	 ((double)1E-5)	  //maximum time step
#define MAX_ERAS	 ((int) 32)	//maximum number of eras in one session
#define MAX_TIME	 ((double)1E-3)	//maximum simulation time

//_________________________________________
// Class definitions
//_________________________________________
class CLOCK
{
	public:
		CLOCK(const CLOCK&);				//copy constructor
		CLOCK();					//default constructor
		~CLOCK();							//destructor
		CLOCK& operator=(const CLOCK&);	//assignment operator
		
		bool operator==(const CLOCK &rhs);		//equality operator for CLOCK
		bool CheckAlarm(const CLOCK&);  //checks for alarm condition relative to another clock
		bool Decrement();    //decrements clock one step
		bool Increment(__int64 alarm=0);	//increment clock one step
		bool Increment(CLOCK &alarmClock);
		bool Reset(unsigned int numEras, double *eLength, double *stepLo, double *stepHi);
		bool Reset();		//reset time to 0 only
		bool StepAdjust(double factor=0.0);
	
		unsigned int eCurrent;

		double curTime;		//current clock
		double curDTime;	//size of most recent time step taken
		
	private:
		//fundamental clock units are FEMTOSECONDS 
		bool Valid();

		unsigned int eMax;			//total number of eras

		__int64 *eraBound;			//coordinates of era boundaries
		__int64 *stepMin;			//minimum step
		__int64 *stepMax;			//maximum step

		__int64 stepPrev;			//previous step taken
		__int64 stepCurrent;		//current step size
		
		__int64 tPrev;			//previous time if must back-up
		__int64 tCurrent;		//current time

}; //endclass
