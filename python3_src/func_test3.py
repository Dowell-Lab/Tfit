#!/Library/Frameworks/Python.framework/Versions/3.6/bin/python3
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
import simulate as sm
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import math as m

def test1(number):

	print('Test: ',number)
	return 0

def callSimulate3():
	####################Python3 Version####################
	#print("Testing Simulate3 Function")
	#run simulations 50- 100 times, get mean and variance for mean and variance

	X1 	= sm.runOne(mu=0, s=20, l=5, lr=100, ll=-100, we=0.5,wl=0.25, wr=0.25, pie=0.5, pil=0.1, pir=0.9, N=10000, SHOW=False , bins=200, noise=False )
	mu1 	= 0.5*(sm.weighted_mean(X1[:,0], X1[:,1])+sm.weighted_mean(X1[:,0], X1[:,2]))
	l1 	= 0.5*(sm.weighted_mean(X1[:,0], X1[:,1])-sm.weighted_mean(X1[:,0], X1[:,2]))
	var1 = m.sqrt(sm.weird_variance(X1[:,0], X1[:,1], mu1, 1))
	var2 = m.sqrt(sm.weird_variance(X1[:,0], X1[:,2], mu1, -1))

	#print ("Mu1: ",mu1,"l1: ", l1)
	#print ("Some Var1: ",m.sqrt(sm.weird_variance(X1[:,0], X1[:,1], mu1, 1)))
	#print ("Some Var2: ",m.sqrt(sm.weird_variance(X1[:,0], X1[:,2], mu1, -1)))


	return mu1, l1, var1, var2

def testSimulate3():
		for i in range(1,50):
			curr = callSimulate3()
			if i == 1:
				mu_sum = curr[0]
				l1_sum = curr[1]
				var1_sum = curr[2]
				var2_sum = curr[3]

			mu_sum = mu_sum + curr[0]
			l1_sum = l1_sum + curr[1]
			var1_sum = var1_sum + curr[2]
			var2_sum = var2_sum + curr[3]

		#calculating means, dont want to declare new variable
		mu_sum = mu_sum/50
		l1_sum = l1_sum/50
		var1_sum = var1_sum/50
		var2_sum = var2_sum/50

		print (mu_sum, l1_sum, var1_sum, var2_sum)

		'''#If I was to write results to a file:
		results = open("Python2_Results","w")

		for i in range(1,50):
			print "run" , i
			results.write("run: ")
			results.write(repr(i))
			results.write( ": ")
			results.write(repr(testSimulate2()))

		results.close()
		'''


def main(args):
    return 0

if __name__ == '__main__':

	testSimulate3()
