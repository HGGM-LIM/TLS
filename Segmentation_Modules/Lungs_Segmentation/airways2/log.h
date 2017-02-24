#ifndef LOG_H
#define LOG_H

#include	<log4cpp/Category.hh>
#include	<log4cpp/OstreamAppender.hh>
#include	<log4cpp/FileAppender.hh>
#include 	<log4cpp/SimpleLayout.hh>

#define LOG_TO_FILE					0b1000
#define LOG_TO_CONSOLE				0b0100


typedef log4cpp::Category* pLogType;

log4cpp::Category* get_logger(const std::string );
log4cpp::Category* get_default_logger(void);
void init_logging(const std::string);
void set_priority(int);

#endif // LOG_H





