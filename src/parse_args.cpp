#include "parse_args.h"

#include <stdio.h>
#include <argp.h>
#include <argz.h>
#include <stdlib.h> // argz vectors are malloc'd so we need the free function.

const char *argp_program_version = "Tfit v1.1 (in dev)";
const char *argp_program_bug_address = "<dowellde@colorado.edu>";

static char doc[] = "Transcription Fit (Tfit): Nascent transcription to RNA polymerase behavior\n";

// Needs to have options from current Tfit.  Currently morse code example.
struct argp_option options[] = {
    {0,0,0,0, "Morse Code Options:", 7},
    {"dot", 'd', "NUM", OPTION_ARG_OPTIONAL, "Show a dot on the screen."},
    {"dash", 888, 0, 0, "Show a dash on the screen."},
    {0,0,0,0, "Informational Options:", -1},
    {"SOS", 999, 0, 0, "Give some help in morse code"},
    {0}
};

// This is a brief usage statement with pre options and post options portions.
struct argp argp = {options, parse_opt, "WORD\nWORD WORD", 
    "Show some dots and dashes on the screen.\v"
      "A final newline is also shown regardless of whether any options were given."};

// This is the main parsing function.  Ultimately will need to check
// option quality/identity, parse them into objects, and error handle.  
int parse_opt(int key, char *arg, struct argp_state * state) {
  struct arguments *a = (struct arguments *)state->input;
  switch(key) {
    case 'd': {
      unsigned int i;
      unsigned int dots = 1;
      if (arg != NULL) {
        dots = atoi(arg);
      }
      for (i = 0; i < dots; i++) printf(".");
      break;
    }
    case 888:
      printf("-");
      break;
    case 777:
      return parse_opt('d',(char *)"3", state);
    case ARGP_KEY_ARG: {
      argz_add(&a->argz, &a->argz_len, arg);
      break;
    }
    case 999: {
      parse_opt('d', (char *)"3", state);
      printf(" ");
      parse_opt(888, NULL, state);
      parse_opt(888, NULL, state);
      parse_opt(888, NULL, state);
      printf(" ");
      parse_opt('d', (char *)"3", state);
      printf("\n");
      exit(0);
      break;
    }
    case ARGP_KEY_INIT: {
      a->argz=0;
      a->argz_len = 0;
      break;
    }
    case ARGP_KEY_END: {
      size_t count = argz_count(a->argz, a->argz_len);
      if (count > 2) {
        argp_failure(state,1,0,"too many arguments");
      } else if (count < 1) {
        argp_failure(state,1,0,"too few arguments");
      }
      break;
    }
  }
  return 0;
}
