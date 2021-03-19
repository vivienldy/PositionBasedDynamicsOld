#pragma once

#ifndef __TIMER_H
#define __TIMER_H

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <map>
#include <fstream>
#include <string>
#include <ostream>

using namespace std;

struct Timer 
{
	std::chrono::time_point<std::chrono::steady_clock> then;
	chrono::duration <double, milli> last, total, diff;
    std::map<string, chrono::duration <double, milli>> func2time;
    Timer();
    void Tick(), Tock();
	chrono::duration <double, milli> Diff();
	void PrintFuncTime();
	double GetFuncTime();
	void GetFuncTime(string funcName);
    void Write2Map(string funcName);
	void SaveMapOut(string filePath);
	char* GetLocalTimeStamp();
};

#endif
