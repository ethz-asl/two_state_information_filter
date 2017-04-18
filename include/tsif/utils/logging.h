#ifndef TSIF_LOGGING_HPP_
#define TSIF_LOGGING_HPP_

#include <iostream>

#if TSIF_VERBOSE > 0
#define TSIF_LOG(msg) std::cout << msg << std::endl
#define TSIF_LOGIF(con,msg) if(con) TSIF_LOG(msg)
#define TSIF_LOGW(msg) std::cout << "\033[33m" << __FILE__ << "(" << __LINE__ << "): " << msg << "\033[0m" << std::endl
#define TSIF_LOGWIF(con,msg) if(con) TSIF_LOGW(msg)
#define TSIF_LOGE(msg) std::cout << "\033[31m" << __FILE__ << "(" << __LINE__ << "): " << msg << "\033[0m" << std::endl
#define TSIF_LOGEIF(con,msg) if(con) TSIF_LOGE(msg)
#else
#define TSIF_LOG(msg)
#define TSIF_LOGIF(con,msg)
#define TSIF_LOGW(msg)
#define TSIF_LOGWIF(con,msg)
#define TSIF_LOGE(msg)
#define TSIF_LOGEIF(con,msg)
#endif

#endif /* TSIF_LOGGING_HPP_ */
