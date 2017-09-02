int mm_verbose = 3;
int mm_dbg_flag = 0;
double mm_realtime0;

#include "minimap.h"

#ifdef WIN32
#include "gettimeofday.h"

// taken from https://stackoverflow.com/questions/5272470/c-get-cpu-usage-on-linux-and-windows
#include <windows.h>
double cputime() {
  HANDLE hProcess = GetCurrentProcess();
  FILETIME ftCreation, ftExit, ftKernel, ftUser;
  SYSTEMTIME stKernel;
  SYSTEMTIME stUser;

  GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser);
  FileTimeToSystemTime(&ftKernel, &stKernel);
  FileTimeToSystemTime(&ftUser, &stUser);

  double kernelModeTime = ((stKernel.wHour * 60.) + stKernel.wMinute * 60.) + stKernel.wSecond * 1. + stKernel.wMilliseconds / 1000.;
  double userModeTime = ((stUser.wHour * 60.) + stUser.wMinute * 60.) + stUser.wSecond * 1. + stUser.wMilliseconds / 1000.;

  return kernelModeTime + userModeTime;
}
#else
#include <sys/resource.h>
#include <sys/time.h>

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}
#endif /* WIN32 */

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

#include "ksort.h"

#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8) 

#define sort_key_64(x) (x)
KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 8)

KSORT_INIT_GENERIC(uint32_t)
