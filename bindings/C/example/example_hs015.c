#include <stdio.h>
#include "Uno_C_API.h"

int main() {
	int32_t uno_major, uno_minor, uno_patch;
	uno_get_version(&uno_major, &uno_minor, &uno_patch);
	printf("Uno v%d.%d.%d\n", uno_major, uno_minor, uno_patch);
	return 0;
}