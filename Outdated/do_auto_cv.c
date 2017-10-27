/* do_auto_cv.c
 *
 * by: Mark D. Butala
 *
 */


#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>


const char *suffix_list[] = {
  "oct12_25_mk",
  "oct19_nov1_mk",
  "nov2_15_mk",
  "nov9_22_mk",
  "",
};


int main(int argc, char **argv) {
  int i;
  char buffer[256];
  
  for (i=0; strcmp(suffix_list[i], "") != 0; i++) {
    sprintf(buffer, "./auto_cv -s %s &> out_%s",
            suffix_list[i], suffix_list[i]);

    printf("%s\n", buffer);
    fflush(stdout);
    assert(system(buffer) == 0);
  }
  
  return 0;
}
