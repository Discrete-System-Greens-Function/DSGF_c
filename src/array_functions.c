#include <array_functions.h>

void double_linspace(double beginning, double end, int num_elements, double filler_array[]){
	double step = (end-beginning) / (num_elements - 1);

	for(int i = 0; i < num_elements; i++){
		filler_array[i] = beginning + i * step;
	}
}
