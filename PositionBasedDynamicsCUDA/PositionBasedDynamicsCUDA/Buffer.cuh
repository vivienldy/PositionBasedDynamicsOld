#pragma once
#ifndef BUFFER_H
#define BUFFER_H

#include "cuda_runtime.h"
#include "helper_math.h"
#include "Utility.h"
#include <unordered_map>


using namespace std;
static std::unordered_map<std::string, int> gTypeMapId;


template<class T>
class Buffer
{
public:
	Buffer() : m_BufferSize(0)
	{
		m_sName = "";
		m_DevicePtr = nullptr;
	}

	typedef T ThisType;
	vector<T> m_Data;

	template<class S>
	S* GetDevicePtr() { return  (S*)m_DevicePtr; }
	void* GetDevicePtr() { return  m_DevicePtr; }

	void** GetDevicePtrRaw() { return (void**)&m_DevicePtr; }

	int GetSize() { return m_Data.size(); }
	int GetCopySize() { return m_Data.size() * sizeof(T); }

	int GetItems()
	{
		return m_Data.size() *this->EvalTupeSize();
	}

	int GetCapacity() { return m_Data.capacity(); }

	uint2 EvalBlockSize(int nThread) { return make_uint2(ceil((float)GetSize() / nThread), nThread); }

	void Save(std::ofstream& os)
	{
		if (m_sName.size() == 0)
			return;
		if (!os)
			return;
		os << m_sName << "|" << typeid(T).name() << "|";
		for (int i = 0; i < m_Data.size(); i++)
			os << m_Data[i] << ";";
		os << std::endl;
	}

	void Save(std::string path)
	{
		std::ofstream ofs(path);
		Save(ofs);
		ofs.close();
	}

	int EvalTupeSize()
	{
		if (typeid(T).name() == std::string("int"))
			return 1;
		if (typeid(T).name() == std::string("float"))
			return 1;
		if (typeid(T).name() == std::string("struct float3"))
			return 3;
		if (typeid(T).name() == std::string("struct int2"))
			return 2;
		return 0;
	}

	void Read(std::string path)
	{
		/*
 		std::ifstream ifs(path);
		if (!os)
		{
			std::cout << "ERROR Path @" << path << std::endl;
			return;
		}
		const std::string text((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>()));

		auto ret = UT::splitString(text,"|");
		m_sName = ret[0];
 		int tupeSize = Buffer<T>::EvalTupeSize(ret[1]);
		*/
	}


	void SetName(std::string n)
	{
		m_sName = n;
	}

	inline bool LoadToHost()
	{

		if (!m_DevicePtr)
		{
			printf("[ %s ] Device Not Registered!\n", m_sName.c_str());
			return false;
		}
		cudaError_t cudaStatus = cudaMemcpy(m_Data.data(), this->GetDevicePtr(), this->GetCopySize(), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			printf("[ %s ]  cudaMemcpy failed\n", m_sName.c_str());
			return false;
		}
		return true;
	}

	inline bool LoadToDevice(bool printSuccInfo = false)
	{
		if (!m_DevicePtr)
		{
			printf("[ %s ] Device Not Registered!\n", m_sName.c_str());
			return false;
		}
		if (m_iDeviceBufferSize < m_Data.size())
		{
			cudaFree(m_DevicePtr);
			m_DevicePtr = nullptr;
			DeviceMalloc();
		}
		cudaError_t cudaStatus = cudaMemcpy(this->GetDevicePtr(), m_Data.data(), this->GetCopySize(), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			printf(" [%s] cudaMemcpy failed whit size [%d]\n", m_sName.c_str(), this->GetCopySize());
			return false;
		}

		m_iDeviceBufferSize = m_Data.size();
		if (printSuccInfo)
			printf("Load [%s] To Device Succeeded\n", m_sName.c_str());
		return true;
	}

	inline bool DeviceMalloc()
	{
		if (this->GetCopySize() == 0)
		{
			printf("[%s] Buffer Size = 0 !!! \n", m_sName.c_str());
			return false;
		}
		if (m_DevicePtr != nullptr)
		{
			printf("[%s] DeviceRegistered!\n", m_sName.c_str());
			return false;
		}

		cudaError_t cudaStatus = cudaMalloc(this->GetDevicePtrRaw(), this->GetCopySize());
		if (cudaStatus != cudaSuccess)
		{
			printf("[%s] CudaMalloc Failed\n", m_sName.c_str());
			return false;
		}
		/*else
		{
			printf("Malloc device memory succ use | %12.3f MB | [%3s ]\n",
				(float)GetCopySize() / 1024.0f / 1024.0f, m_sName.c_str());
		}*/
		return true;
	}

	inline bool MallocAndLoadToDevice()
	{
		if (!this->m_DevicePtr)
			this->DeviceMalloc();
		if (!this->LoadToDevice())
			return false;
		return true;
	}

	std::string GetName()
	{
		return m_sName;
	}

private:
	std::string m_sName;
	void* m_DevicePtr;
	long int m_iDeviceBufferSize;
	long int m_BufferSize;
};

typedef Buffer<int> BufferInt;
typedef Buffer<unsigned int> BufferUInt;
typedef Buffer<float> BufferFloat;
typedef Buffer<float3> BufferVector3f;
typedef Buffer<int2> BufferInt2;

#endif