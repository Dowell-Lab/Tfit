#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  function_tests.py
#  
#  Copyright 2018 Jack Dempsey <jackdempsey@rgnt2-102-74-dhcp.int.colorado.edu>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#
#  
#This script is made for running different tests on the Tfit/pytho3_src codebase
#each funciton should have a description of the test it runs and the proper files/functions it should be used with
#
#
#

def test1(number):
	
	print('Test: ',number)


def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
