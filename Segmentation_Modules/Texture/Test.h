/*
 * Test.h
 *
 *  Created on: Apr 6, 2016
 *      Author: pmacias
 */

#ifndef TEST_H_
#define TEST_H_

class Test {
public:
	Test();
	virtual ~Test();
	void hello(std::string password = "Hola Test");
};
#ifndef ITK_MANUAL_INSTANTIATION
#include "Test.txx"
#endif

#endif /* TEST_H_ */
