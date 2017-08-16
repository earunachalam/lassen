#ifndef DEBUGTOOLS_H
#define DEBUGTOOLS_H

#include <csignal>
#include <cxxabi.h>
#include <iostream>

// print name and value of variable
#define P2S(a) std::cout << __FILE__ << " L" << __LINE__ << ": " << #a << ": " << (a) << std::endl
// print name and class of variable
#define GC(A) {int status; char * demangled = abi::__cxa_demangle(typeid(A).name(),0,0,&status); std::cout << __LINE__ << ": " #A << "\t" << demangled <<"\n";}


#endif
