#include "BaseProcess.h"


BaseProcess::BaseProcess()
{
	MPI_Comm_size(MPI_COMM_WORLD,&m_num_processes);
	MPI_Comm_rank(MPI_COMM_WORLD,&m_process_id);

};
void BaseProcess::Log(const char* format,...){
    if(m_process_id == 0){
    	va_list args;
    	va_start( args, format );
    	printf("[%d]: ",m_process_id);
    	vprintf(  format, args );
    	va_end( args );	
    }


}
void BaseProcess::VerboseLog(const char* format,...){
    	
    va_list args;
    va_start( args, format );
    printf("[%d]: ",m_process_id);
    vprintf( format, args );
    va_end( args );	

}


void BaseProcess::BaseProcess::LogError(const char* format,...){
    	
    va_list args;
    va_start( args, format );
    fprintf(stderr,"[%d]: ",m_process_id);
    vfprintf(stderr, format, args );
    va_end( args );	

}

void BaseProcess::LogErrorAndAbort(const char* format,...){

    va_list args;
    va_start( args, format );
    fprintf(stderr,"[%d]: ",m_process_id);
    vfprintf(stderr, format, args );
    va_end( args );

    MPI_Abort(MPI_COMM_WORLD,0);

}
