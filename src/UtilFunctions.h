#pragma once

template<typename T>
T
sumVector(int len, T const * vec)
{
	T result = 0;
	for (int i = 0; i < len; i++)
	{
		result += vec[i];
	}
	return result;
}

//int
//sumVector(int len, int* const vec);

void
printVector(int len, int* const vec);
