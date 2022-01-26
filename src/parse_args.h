/**
 * @file parse_args.h
 * @author Robin Dowell
 * @brief Beginnings of an argp interface
 * @version 0.1
 * @date Mon Jan 24 12:42:22 PM MST 2022
 */

#ifndef parse_args_H 
#define parse_args_H 

#include <stdio.h>
#include <argp.h>
#include <argz.h>
#include <stdlib.h> // argz vectors are malloc'd so we need the free function.

struct arguments {
  char *argz;
  size_t argz_len;
};

extern struct argp_option options[];
extern struct argp argp;

int parse_opt(int key, char *arg, struct argp_state *state);

#endif
