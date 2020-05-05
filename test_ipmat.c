#include <math.h>
#include <stdio.h>
#include "ip_lib.h"

int main () {
    print_ipmat(ip_mat_create(5, 6, 3, 0));
    return 0;
}

