#include "timer.h"
using namespace std;

Timer::Timer (): last(0), total(0), diff(0)
{
    Tick();
}

void Timer::Tick () 
{
    then = std::chrono::steady_clock::now();
}

void Timer::Tock () 
{
	auto now = std::chrono::steady_clock::now();
    last = now - then;
    total += last;
    then = now;
}

chrono::duration <double, milli> Timer::Diff()
{
	return last;
}

void Timer::PrintFuncTime()
{
	cout << last.count() << " ms" << endl;
	//cout << func2time[funcName].count() << " ms" << endl;
}

double Timer::GetFuncTime()
{
	 return last.count();
}

void Timer::GetFuncTime(string funcName)
{
	//cout << last.count() << " ms" << endl;
	cout << func2time[funcName].count() << " ms" << endl;
}

void Timer::Write2Map(string funcName)
{
	cout << funcName << "insert" << endl;
	func2time.insert(map<string, chrono::duration <double, milli>>::value_type(funcName, last));
}

void Timer::SaveMapOut(string filePath)
{
	ofstream ofs(filePath);
	if (!ofs.is_open())
		return;
	for (map<string, chrono::duration <double, milli>>::iterator iter = func2time.begin(); iter != func2time.end(); ++iter)
	{
		ofs << iter->first << " ------ " << iter->second.count() << " ms\n";
	}

	ofs.flush();
	ofs.close();
}

char* Timer::GetLocalTimeStamp()
{
	char str[50];
	time_t now = time(NULL);
	strftime(str, 50, "%x %X", localtime(&now));
	cout << str << endl;
	return str;
}
