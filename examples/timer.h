#ifndef TIMER_H
#define TIMER_H 1

#include <iostream>
#include <string>
#include <sys/time.h>

class Timer{

private:
    timeval _start;
    timeval _stop;
    
public:
    
    void
    start(){
        gettimeofday(&_start, 0);
    }
    
    void
    stop(){
        gettimeofday(&_stop, 0);
    }
    
    double
    elapsed(){
        long long int usecs_start = _start.tv_sec * 1000000. + _start.tv_usec;
        long long int usecs_stop =  _stop.tv_sec * 1000000. + _stop.tv_usec;
        return double(usecs_stop - usecs_start) / 1000000.;
    }
    
    void
    out(int proc, std::string point, bool timings) {
    	
    	if (timings) {
	    	stop();
	    	
    		MPI::COMM_WORLD.Barrier();
    		if (point == "load") {
    	
    			if (proc == 0) {
    				std::cout << "Times for initialisation/data loading:" << std::endl;
    			}
    		
    		}
	    	else if (point == "mesh") {
    		
    			if (proc == 0) {
    				std::cout << "\nTimes for mesh generation and refinement:" << std::endl;
	    		}
	
    		}
    		else if (point == "assembly") {
    	
    			if (proc == 0) {
    				std::cout << "\nTimes for system assembly:" << std::endl;
    			}
    		
    		}	
	    	else if (point == "solving") {
    	
    			if (proc == 0) {
    				std::cout << "\nTimes for solution solving:" << std::endl;
    			}
    		
    		}
    		else {
    			return;
    		}
    	
	    	MPI::COMM_WORLD.Barrier();
    		
    		std::cout << "      Node " << proc << ": " << elapsed() << " seconds." << std::endl;
    		
    		MPI::COMM_WORLD.Barrier();
    		
	    	start();	
	    	
	    }
    }	
    
    void
    outNoMPI(std::string point, bool timings) {
    	
    	if (timings) {
	    	stop();
	    	
    		if (point == "load") {
    	
    				std::cout << "Times for initialisation/data loading:" << std::endl;
    		
    		}
	    	else if (point == "mesh") {
    	
    				std::cout << "\nTimes for mesh generation and refinement:" << std::endl;
	
    		}
    		else if (point == "assembly") {
    	
    				std::cout << "\nTimes for system assembly:" << std::endl;
    				    		
    		}	
	    	else if (point == "solving") {
    	
    				std::cout << "\nTimes for solution solving:" << std::endl;
    		
    		}
    		else {
    			return;
    		}
    		
    		std::cout << "      Serial" << proc << ": " << elapsed() << " seconds." << std::endl;

	    	start();	
	    	
	    }
    }			
};

#endif // TIMER_H