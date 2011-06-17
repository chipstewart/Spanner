#pragma once

#ifndef WIN32
#ifndef SAFEFUNCTIONS_H
#define SAFEFUNCTIONS_H

#include <stdarg.h>

typedef int errno_t;

inline int sprintf_s(char *buffer, size_t sizeOfBuffer, const char *format, ...) {
  va_list argp;
  va_start(argp, format);
  unsigned int err = vsprintf(buffer, format, argp);
  va_end(argp);
  return err;
}

inline errno_t strncpy_s(char *strDest, size_t sizeInBytes, const char *strSource, size_t count) {
	strncpy(strDest, strSource, count);
	return 0;
}

inline errno_t strcat_s(char *strDestination, size_t sizeInBytes, const char *strSource) {
	strcat(strDestination, strSource);
	return 0;
}

inline errno_t fopen_s(FILE** pFile, const char *filename, const char *mode) {
	*pFile = fopen(filename, mode);
	return 0;
}

inline errno_t localtime_s(struct tm* _tm, const time_t *time) {
	memcpy(_tm, localtime(time), sizeof(*_tm));
	return 0;
}

inline errno_t gmtime_s(struct tm* _tm, const time_t *time) {
	memcpy(_tm, gmtime(time), sizeof(*_tm));
	return 0;
}

inline errno_t asctime_s(char* buffer, size_t sizeInBytes, const struct tm *_tm) {
	char* time = asctime(_tm);
	unsigned int timeLen = strlen(time);
	strncpy(buffer, time, timeLen);
	buffer[timeLen] = 0;
	return 0;
}

inline char* strtok_s(char* strToken, const char *strDelimit, char **context) {
	return strtok(strToken, strDelimit);
}

#endif
#endif
