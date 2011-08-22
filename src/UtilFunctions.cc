#include "R.h"
#include "UtilFunctions.h"

void
printVector(int len, int* const vec)
{
	Rprintf("c(");
	for (int i = 0; i < len; i++)
	{
		Rprintf("%d ", vec[i]);
	}
	Rprintf(")\n");
	R_FlushConsole();
}
