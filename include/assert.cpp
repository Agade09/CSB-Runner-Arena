#pragma once
#include <bits/stdc++.h>
using namespace std;

#ifdef TESTS
#define ASSERT(cond,error_message) \
	if(!(cond)){ \
		cerr << error_message << endl; \
		throw(0); \
	}
#else
#define ASSERT(cond,error_message)
#endif

/*inline void ASSERT(const bool check,const string error_message){
	if(tests && !check){
		cerr << error_message << endl;
		throw(0);
	}
}*/