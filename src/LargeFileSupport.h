#pragma once
#include <iostream>
#include <fstream>
#include <string>
#ifdef WIN32
#define ftell64(a)     _ftelli64(a)
#define fseek64(a,b,c) _fseeki64(a,b,c)
typedef __int64 off_type;
#elif __APPLE__
//#error mac osx GNU C compiler
typedef off_t off64_t;
#define fopen64(a,b)    fopen(a,b)
#define ftello64(a)     ftello(a)
#define fseeko64(a,b,c) fseeko(a,b,c)
#define ftell64(a)      ftello(a)
#define fseek64(a,b,c)  fseeko(a,b,c)
typedef off_t off_type;
#else
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef __off64_t off_type;
#endif
