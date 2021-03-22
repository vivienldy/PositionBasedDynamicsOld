#ifndef __BUFFER_DEBUG_MODULE_H__
#define __BUFFER_DEBUG_MODULE_H__


//AB测试？回归测试

//A B 结果要一样

//callstack



#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "Utility.h"

#define GREEN_START					"\033[40;32m" 
#define RED_START						"\033[40;31m" 
#define YELLOW_START				"\033[40;33m" 
#define BLUE_START						"\033[40;36m" 
#define COLOR_END						"\033[0m" 
#define GREEN_BOLD_START		"\033[1m\033[32m"      /* Bold Green */

#define ERROR						"\033[40;31m  [ERROR]" 


class BufferDebugModule
{
private:
	BufferDebugModule()
	{
		m_sStoragePath = "D:/0319CCDTest/BufferArray/";
	}
	~BufferDebugModule()
	{

	}

	static BufferDebugModule* m_Instance;
	std::string m_sStoragePath;

	// __FUNCTION__LINE_BUFFER_NAME__
	//
 	std::map<int, std::unordered_map<std::string, int*>>    m_IntMap;
	std::map<int, std::unordered_map<std::string, float*>> m_FloatMap;


	std::string _toRefString(std::string str) { return std::string("\"") + str + std::string("\""); };

public:

	int IntBufferCompare(int* a, int* b, unsigned int size)
	{
		int c = 0;
		for (int i = 0; i < size; i++)
			c += abs(a[i] - b[i]);
		return c;
	}

	float FloatBufferCompare(float* a, float* b, unsigned int size)
	{
		float c = 0;
		for (int i = 0; i < size; i++){
			float vv = a[i] - b[i];
			c += vv * vv;
		}
		return sqrt(c);
	}


	static BufferDebugModule* GetInstance() 
	{
		if (m_Instance == nullptr)
			m_Instance = new BufferDebugModule();
		return m_Instance;
	};

	
	template<typename T,typename S>
	bool CompareStack(int cookTimes, std::string key, Buffer<T>& buffer ,S tolerance)
	{
		if (typeid(T).name() == std::string("int"))
		{
			if (m_IntMap.find(cookTimes) == m_IntMap.end())
			{
				std::cout << ERROR << "Non-CookTimes : " << cookTimes << COLOR_END << std::endl;
				return false;
			}

			if (m_IntMap[cookTimes].find(key) == m_IntMap[cookTimes].end())
			{
				std::cout << ERROR << "Non-Key : " << _toRefString(key)<<"@CookTime:"<< cookTimes << COLOR_END << std::endl;
				return false;
			}

			int tol =  IntBufferCompare(m_IntMap[cookTimes][key], (int*)buffer.m_Data.data(), buffer.GetSize());
			if (tol > tolerance)
			{
				std::cout << ERROR << "Too biggg .. referance : " << tolerance << " current : " << tol << COLOR_END << std::endl;
				return false;
			}

			std::cout << GREEN_START << "Compare @: " << cookTimes
				<< " by key: " << _toRefString(key)<< " with " << buffer.GetName().c_str()
				<< " tolerance : " << tol << " succ !!!" << COLOR_END << std::endl;
			return true;
		}
		if (typeid(T).name() == std::string("struct float3"))
		{
			if (m_FloatMap.find(cookTimes) == m_FloatMap.end())
			{
				std::cout << ERROR << "Non-CookTimes : " << cookTimes << COLOR_END << std::endl;
				return false;
			}

			if (m_FloatMap[cookTimes].find(key) == m_FloatMap[cookTimes].end())
			{
				std::cout << ERROR << "Non-Key : " << _toRefString(key)<< "@CookTime:" << cookTimes << COLOR_END << std::endl;
				return false;
			}

			float tol = FloatBufferCompare(m_FloatMap[cookTimes][key], (float*)buffer.m_Data.data(), buffer.GetSize());
			if (tol > tolerance)
			{
				std::cout << ERROR << "Too biggg .. referance : " << tolerance << " current : " << tol << COLOR_END << std::endl;
				return false;
			}

			std::cout << GREEN_START << "Compare @: " << cookTimes 
				<< " by key: " << _toRefString(key)<< " with " << buffer.GetName().c_str() 
				<< " tolerance : " << tol << " succ !!!" << COLOR_END << std::endl;
			return true;
		}

		std::cout << ERROR << " TypeId : "<< typeid(T).name() << COLOR_END << std::endl;
		return false;
	}

	template<typename T>
	void PushStack(int cookTimes, std::string key, Buffer<T>& buffer)
	{
		if (typeid(T).name() == std::string("int"))
		{
			if (m_IntMap[cookTimes].find(key) != m_IntMap[cookTimes].end())
			{
				std::cout << ERROR << "Same key in stack ? " << COLOR_END << std::endl;
				return;
			}
			int* data = new int[buffer.GetItems()];
			std::memcpy(data, buffer.m_Data.data(), buffer.GetCopySize());
			m_IntMap[cookTimes][key] = data;
			std::cout << YELLOW_START << "PushStack " << typeid(T).name() << " data buffer succ @" << cookTimes << "cookTimes with key : " << _toRefString(key) << COLOR_END << std::endl;
		}
		if (typeid(T).name() == std::string("struct float3"))
		{
			if (m_FloatMap[cookTimes].find(key) != m_FloatMap[cookTimes].end())
			{
				std::cout << ERROR << "Same key in stack ? " << COLOR_END << std::endl;
				return;
			}
			float* data = new float[buffer.GetItems()];
			std::memcpy(data, buffer.m_Data.data(), buffer.GetCopySize());
			m_FloatMap[cookTimes][key] = data;
			std::cout << YELLOW_START << "PushStack " << typeid(T).name() << " data buffer succ @" << cookTimes << "cookTimes with key : " << _toRefString(key) << COLOR_END << std::endl;
		}
 	}

	void Load()
	{
		std::cout<< YELLOW_START << "Succ ENTRY" <<COLOR_END << std::endl;
	}
	
};

BufferDebugModule*  BufferDebugModule::m_Instance = nullptr;

#if 1
	#define PUSH_STACK(cooktime,key,buffer) BufferDebugModule::GetInstance()->PushStack(cooktime,key,buffer)
#else
	#define PUSH_STACK(cooktime,key,buffer) 
#endif

#endif 