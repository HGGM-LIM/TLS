/*
 * log.cxx
 *
 *  Created on: May 2, 2011
 *      Author: mario
 */

#include "log.h"
#include <sstream>
#include <stdlib.h>

/*
 * Default values
 */
int m_priority=log4cpp::Priority::INFO;
int m_output =  LOG_TO_CONSOLE;
log4cpp::Category* m_category = &log4cpp::Category::getRoot();

void set_priority(int newp){
	m_priority = newp;
}

void init_logging(const std::string name){
	std::string env_var = "EXTRACT_AIR_LOG_LEVEL";

	/* See  log4cpp::Priority for valid values:
	 *    EMERG  = 0,
	 *    FATAL  = 0,
	 *	  ALERT  = 100,
	 *	  CRIT   = 200,
	 *	  ERROR  = 300,
	 *	  WARN   = 400,
	 *	  NOTICE = 500,
	 * 	  INFO   = 600,
	 *	  DEBUG  = 700,
	 * 	  NOTSET = 800
	 */

	char * pl = getenv(env_var.c_str());
	if (pl !=NULL) {
	  m_priority=atoi(pl);
	  std::cout << "Log level set from environment var " << env_var << " is " << m_priority << std::endl;
	}


	m_category->setPriority(m_priority);

	if ((m_output & LOG_TO_FILE) != 0) {

		std::string filename = name + ".log";

		log4cpp::Layout *layout = new log4cpp::BasicLayout();

		log4cpp::FileAppender *appender = new log4cpp::FileAppender("FileAppender", filename);
			appender->setLayout(layout);

		m_category->addAppender(appender);
	}

	if ((m_output & LOG_TO_CONSOLE) != 0) {
		log4cpp::Layout *layout = new log4cpp::BasicLayout();


		log4cpp::OstreamAppender *appender = new log4cpp::OstreamAppender("OstreamAppender", &std::cout);
			appender->setLayout(layout);

		m_category->addAppender(appender);
	}

}

log4cpp::Category* get_default_logger() {

	return m_category;
}

log4cpp::Category* get_logger(const std::string name) {


	log4cpp::Category* cat = &log4cpp::Category::getInstance(name);
	cat->setPriority(m_priority);

	return cat;
}
