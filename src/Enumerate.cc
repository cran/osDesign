#include "R.h" // R memory io
#include "Rmath.h" // R math functions
#include "Rinternals.h" // R internal functions
#include <list> // List Library
#include <algorithm> // for max and min

#include "UtilFunctions.h"

/*
enumerate.nn = function(MM,NN){
  if( length(MM) == 1){
    if(sum(NN) == sum(MM)){return(c())}
    else{return()}
  }
  else{if(length(MM)>1){
    result = c()
    for(i in max(0,NN[2]-sum(MM[2:length(MM)])):min(MM[1],NN[2])){
      res = enumerate.nn(MM[2:length(MM)], NN - c(MM[1]-i,i))
      result = cbind(result, rbind(i, res))
    }
  }}
  return(result)
}
*/


double
EnumerateCount(int *MM, int MMlen, int sumMM, int NN[2]){
	if( MMlen == 1 ){
		if(NN[0]+NN[1] == MM[0]){
			return 1;
		}
		return 0;
	}
	else{
		const int lower = std::max(0, NN[1] - (sumMM - MM[0]));
		const int upper = std::min(MM[0], NN[1]);

		double result = 0;
		for(int i = lower; i <= upper; i++)
		{
			NN[0] -= MM[0]-i;
			NN[1] -= i; 
			result += EnumerateCount(MM+1, MMlen-1, sumMM-MM[0], NN);
			NN[0] += MM[0]-i;
			NN[1] += i;
		}
		return result;
	}
}

/*
double
LogLikelihood(int *MM, int *Nyx, int *beta, int size){
	double *temp = alloca(size*sizeof(double));
	for(int i = 0; i < size; i++){
		temp[i] = lchoose(MM[i], Nyx[i]);
	}
}

void
EnumerateCall(int *MM, int MMlen, int MMlenTot, int sumMM, int NN[2], int *partial, double &result, double (*Func)(int*,int*,int*)){
	if( MMlen == 1 ){
		if(NN[0]+NN[1] == MM[0]){
			partial[MMlenTot-1] = NN[1];
			result += Func(partial);
		}
	}
	else{
		const int lower = std::max(0, NN[1] - (sumMM - MM[0]));
		const int upper = std::min(MM[0], NN[1]);
		for(int i = lower; i <= upper; i++)
		{
			NN[0] -= MM[0]-i;
			NN[1] -= i;
			partial[MMlenTot - MMlen] = i; 
			EnumerateCall(MM+1, MMlen-1, MMlenTot, sumMM-MM[0], NN, partial, result, Func);
			NN[0] += MM[0]-i;
			NN[1] += i;
		}
	}
}
*/



void
EnumerateWindow_into(int *MM, int MMlen, int MMlenTot, int sumMM, int NN[2], int *partial, int* &result, double &start, double &window){
	if( MMlen == 1 ){
		if(NN[0]+NN[1] == MM[0]){
#ifdef DEBUGGING
				partial[MMlenTot-1] = NN[1];
				printVector(MMlenTot, partial);
#endif
			if( start <= 0 && window > 0){
				memcpy(result, partial, (MMlenTot-1)*sizeof(int));
				result[MMlenTot-1] = NN[1];
				result += MMlenTot;
				window--;
			}
			else{ start--; }
		}
	}
	else{
		const int lower = std::max(0, NN[1] - (sumMM - MM[0]));
		const int upper = std::min(MM[0], NN[1]);
		for(int i = lower; i <= upper && window > 0; i++)
		{
			NN[0] -= MM[0]-i;
			NN[1] -= i;
			partial[MMlenTot - MMlen] = i; 
			EnumerateWindow_into(MM+1, MMlen-1, MMlenTot, sumMM-MM[0], NN, partial, result, start, window);
			NN[0] += MM[0]-i;
			NN[1] += i;
		}
	}
}

SEXP EnumerateWindowWrap_into(int *MM, const int MMlen, int NN[2], double start, double window){
	int* partial=(int*)malloc(MMlen*sizeof(int));
	const int MMsum = sumVector<int>(MMlen, MM);

	int total = window;

	SEXP r_result;
	PROTECT(r_result = allocMatrix(INTSXP, MMlen, total));
	int* r_result_buf = INTEGER(r_result);

	EnumerateWindow_into(MM, MMlen, MMlen, MMsum, NN, partial, r_result_buf, start, window);
	
	free(partial);
	return r_result;
}


void
EnumerateStep_into(int *MM, int MMlen, int MMlenTot, int sumMM, int NN[2], int *partial, int* &result){
	if( MMlen == 1 ){
		if(NN[0]+NN[1] == MM[0]){
#ifdef DEBUGGING
				partial[MMlenTot-1] = NN[1];
				printVector(MMlenTot, partial);
#endif
			memcpy(result, partial, (MMlenTot-1)*sizeof(int));
			result[MMlenTot-1] = NN[1];
			result += MMlenTot;
		}
	}
	else{
		const int lower = std::max(0, NN[1] - (sumMM - MM[0]));
		const int upper = std::min(MM[0], NN[1]);
		for(int i = lower; i <= upper; i++)
		{
			NN[0] -= MM[0]-i;
			NN[1] -= i;
			partial[MMlenTot - MMlen] = i; 
			EnumerateStep_into(MM+1, MMlen-1, MMlenTot, sumMM-MM[0], NN, partial, result);
			NN[0] += MM[0]-i;
			NN[1] += i;
		}
	}
}

SEXP EnumerateWrap_into(int *MM, const int MMlen, int NN[2]){
	int* partial=(int*)malloc(MMlen*sizeof(int));
	const int MMsum = sumVector<int>(MMlen, MM);

	int total = EnumerateCount(MM, MMlen, MMsum, NN);

	SEXP r_result;
	PROTECT(r_result = allocMatrix(INTSXP, MMlen, total));
	int* r_result_buf = INTEGER(r_result);

	EnumerateStep_into(MM, MMlen, MMlen, MMsum, NN, partial, r_result_buf);
	
	free(partial);
	return r_result;
}

double
EnumerateCountWrap_into(int *MM, const int MMlen, int NN[2]){
	const int MMsum = sumVector<int>(MMlen, MM);

	return EnumerateCount(MM, MMlen, MMsum, NN);
}


void EnumerateStep(int *MM, int MMlen, int MMlenTot, int sumMM, int NN[2], int *partial, std::list<int*>& result){
	if( MMlen == 1 ){
		if(NN[0]+NN[1] == MM[0]){
#ifdef DEBUGGING
				partial[MMlenTot-1] = NN[1];
				printVector(MMlenTot, partial);
#endif
			int* tempResult = (int*)malloc(MMlenTot * sizeof(int));
			memcpy(tempResult, partial, (MMlenTot-1)*sizeof(int));
			tempResult[MMlenTot-1] = NN[1];
			result.push_back(tempResult);
		}
	}
	else{
		const int lower = std::max(0, NN[1] - (sumMM - MM[0]));
		const int upper = std::min(MM[0], NN[1]);
		for(int i = lower; i <= upper; i++)
		{
			NN[0] -= MM[0]-i;
			NN[1] -= i;
			partial[MMlenTot - MMlen] = i; 
			EnumerateStep(MM+1, MMlen-1, MMlenTot, sumMM-MM[0], NN, partial, result);
			NN[0] += MM[0]-i;
			NN[1] += i;
		}
	}
}

SEXP EnumerateWrap(int *MM, const int MMlen, int NN[2]){
	std::list<int*> result;
	int* partial=(int*)malloc(MMlen*sizeof(int));

	EnumerateStep(MM, MMlen, MMlen, sumVector(MMlen, MM), NN, partial, result);
	
	SEXP r_result;
	PROTECT(r_result = allocMatrix(INTSXP, MMlen, result.size()));
	
	int at = 0;
	for (std::list<int*>::iterator i = result.begin(), e = result.end(); i != e; i++) {
		memcpy(&(INTEGER(r_result)[at*MMlen]), *i, MMlen*sizeof(int));
		at++;
		free(*i);
	}

	free(partial);
	return r_result;
}


extern "C" { 

	/* Input: Vector of margin totals, MM
	 *        Vector of outcome totals (length = 2), NN 
	 * Output: Matrix with all possible cell values for N1.  Dimensions (MM-1)*(unknown)
	 */
	SEXP Enumerate(const SEXP _MM, const SEXP _NN, const SEXP _MMlen) {
		SEXP MM = coerceVector(_MM, INTSXP);
		SEXP NN = coerceVector(_NN, INTSXP);

		int MMlen = INTEGER(coerceVector(_MMlen, INTSXP))[0];

#ifdef DEBUGGING
		Rprintf("%d %d\n", INTEGER(MM)[0], INTEGER(MM)[1]);
		Rprintf("%d\n", MMlen);
#endif

		SEXP result = EnumerateWrap_into(INTEGER(MM), MMlen, INTEGER(NN));

		// R_FlushConsole();
		R_ProcessEvents();

		UNPROTECT_PTR(result);
		return result; // Return Nothing.
	}

	SEXP		
	EnumerateCount(const SEXP _MM, const SEXP _NN, const SEXP _MMlen) {
		SEXP MM = coerceVector(_MM, INTSXP);
		SEXP NN = coerceVector(_NN, INTSXP);

		int MMlen = INTEGER(coerceVector(_MMlen, INTSXP))[0];

#ifdef DEBUGGING
		Rprintf("%d %d\n", INTEGER(MM)[0], INTEGER(MM)[1]);
		Rprintf("%d\n", MMlen);
#endif

		SEXP result;
		PROTECT (result = allocMatrix(REALSXP, 1, 1));

		REAL(result)[0] = EnumerateCountWrap_into(INTEGER(MM), MMlen, INTEGER(NN));

		R_ProcessEvents();

		UNPROTECT_PTR(result);
		return result; // Return Nothing.
	}


	
	SEXP EnumerateWindow(const SEXP _MM, const SEXP _NN, const SEXP _MMlen, const SEXP _start, const SEXP _window) {
		SEXP MM = coerceVector(_MM, INTSXP);
		SEXP NN = coerceVector(_NN, INTSXP);

		int MMlen = INTEGER(coerceVector(_MMlen, INTSXP))[0];
		double start = REAL(coerceVector(_start, REALSXP))[0];
		double window = REAL(coerceVector(_window, REALSXP))[0];

#ifdef DEBUGGING
		Rprintf("%d %d\n", INTEGER(MM)[0], INTEGER(MM)[1]);
		Rprintf("%d\n", MMlen);
#endif

		SEXP result = EnumerateWindowWrap_into(INTEGER(MM), MMlen, INTEGER(NN), start, window);

		// R_FlushConsole();
		R_ProcessEvents();

		UNPROTECT_PTR(result);
		return result; // Return Nothing.
	}

}
