#ifndef _CPUTIME_H_
#define _CPUTIME_H_
#include <iostream.h>
#include <iomanip.h>

double mytime();

class init_time
    {
public:
    double dummy;
    init_time() { dummy = mytime(); }
    friend class cpu_time;
    };

class cpu_time
    {
public:
    double time;		// in seconds
    friend ostream & operator << (ostream & s, const cpu_time & t);
    cpu_time()
	{ time = mytime(); }
    void mark() 
	{ time = mytime(); }
    cpu_time sincemark();

    static init_time init;
    };

#ifdef THIS_IS_MAIN
init_time cpu_time::init;
#endif

#endif //_CPUTIME_H_
