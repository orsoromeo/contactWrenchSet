#include <iostream>
#include <thread>
#include "setoper.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <Eigen/Dense>
#include <realtime_tools/realtime_publisher.h>

int main() {

	typedef Eigen::Vector3d ContactWrenchSet;
//	realtime_tools::RealtimeBuffer<ContactWrenchSet> contact_wrench_buffer_;

	int cws_thread_test();
	//  void distancePointHyperPlane(const Eigen::VectorXd);

	clock_t start = clock();
	int i;
	for (i=0; i<10000;i++){
		cws_thread_test();
	}
	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf ("seconds: %f  \n", seconds);
	return 0;
}

//This function will be called from a thread

void call_from_thread() {
	std::cout << "Hello, World" << std::endl;
}

int cws_thread_test() {
	//Launch a thread
	std::thread t1(call_from_thread);

	//Join the thread with the main thread
	t1.join();

	return 0;
}

