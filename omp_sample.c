#include <stdio.h>
#include <omp.h>

int main(void) {
  #pragma omp parallel
  {
  int i;
  #pragma omp for
  /* insert directive here */
    for (i=0; i<20; i++) {
      printf("myid=%d : i=%d\n", omp_get_thread_num(), i);
    }
  }
  return 0;
}