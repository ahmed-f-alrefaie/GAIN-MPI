
#pragma once
#include <mpi.h>
#include <cstdio>
#include <cstdarg>


//This is inherited by every class in order to manage data between different processes
class BaseProcess{
protected:
	int m_process_id;
	int m_num_processes;
	BaseProcess();
public:
	void Log(const char* format,...);
	void LogError(const char* format,...);
	void LogErrorAndAbort(const char* format,...);
	void VerboseLog(const char* format,...);
	//void BlockAll();
};
