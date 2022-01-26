#include <stdio.h>
#include <argp.h>
#include <argz.h>
#include <stdlib.h> // argz vectors are malloc'd so we need the free function.

#include "parse_args.h"

int main(int argc, char **argv) {
  struct arguments arguments;
  if (argp_parse(&argp, argc, argv, 0, 0, &arguments) == 0) {
    const char *prev = NULL;
    char *word;
    while ((word = argz_next(arguments.argz, arguments.argz_len, prev))) {
      printf(" %s", word);
      prev = word;
    }
    printf("\n");
    free(arguments.argz);
  }
  return 0;
}

