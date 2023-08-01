#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 14:08:31 2023

@author: hannahgermaine
"""

def single_option_loop(question, options):
	"""This function asks the user for input on a single selection from a list"""
	ask_loop = 1
	while ask_loop == 1:
		print(question)
		[print('\t' + str(i) + ': ' + str(options[i]) + '\n') for i in range(len(options))]
		answer = input("Enter index of selection: ")
		try:
			int_answer = int(answer)
			if 0 <= int_answer <= len(options):
				print("Selected option: " + str(options[int_answer]))
				ask_loop = 0
			else:
				print("Please try again. Incorrect entry.")
		except:
			print("Please try again. Incorrect entry.")
	
	return int_answer

def list_loop(question, options):
	"""This function asks the user for a list of inputs"""
	ask_loop = 1
	while ask_loop == 1:
		print(question)
		[print('\t' + str(i) + ': ' +  str(options[i]) + '\n') for i in range(len(options))]
		answer = input("Enter indices from above list, comma separated: ")
		try:
			split_answer = answer.split(',')
			individual_selections = [int(o) for o in split_answer]
			if len(individual_selections) > 0:
				ask_loop = 0
			else:
				print("Please try again. Incorrect entry.")
		except:
			print("Please try again. Incorrect entry.")
	
	return individual_selections

def single_int_loop(question):
	"""This function asks the user for an input of a float"""
	ask_loop = 1
	while ask_loop == 1:
		print(question)
		answer = input("Enter integer value: ")
		try:
			int_answer = int(answer)
			print("Given value: " + str(answer))
			ask_loop = 0
		except:
			print("Please try again. Incorrect entry.")
	
	return int_answer

def multi_int_loop(question):
	"""This function asks the user for a list of integers"""
	ask_loop = 1
	while ask_loop == 1:
		print(question)
		answer = input("Enter integer values, comma separated: ")
		try:
			split_answer = answer.split(',')
			individual_selections = [int(o) for o in split_answer]
			if len(individual_selections) > 0:
				ask_loop = 0
			else:
				print("Please try again. Incorrect entry.")
		except:
			print("Please try again. Incorrect entry.")
	
	return individual_selections

def single_float_loop(question):
	"""This function asks the user for an input of a float"""
	ask_loop = 1
	while ask_loop == 1:
		print(question)
		answer = input("Enter float value: ")
		try:
			float_answer = float(answer)
			print("Given value: " + answer)
			ask_loop = 0
		except:
			print("Please try again. Incorrect entry.")
	
	return float_answer

def multi_float_loop(question):
	"""This function asks the user for a list of floats"""
	ask_loop = 1
	while ask_loop == 1:
		print(question)
		answer = input("Enter float values, comma separated: ")
		try:
			individual_selections = [float(o) for o in answer.split(',')]
			if len(individual_selections) > 0:
				ask_loop = 0
			else:
				print("Please try again. Incorrect entry.")
		except:
			print("Please try again. Incorrect entry.")
	
	return individual_selections

def more_stuff_loop():
	"""This function asks the user if they want to keep running things or end 
	the session"""
	ask_loop = 1
	while ask_loop == 1:
		repeat_input = input("Would you like to perform any other actions [y/n]? ")
		if repeat_input != 'y' and repeat_input != 'n':
			print("Incorrect entry, please try again.")
		elif repeat_input == 'y':
			answer = 1
			ask_loop = 0
		elif repeat_input == 'n':
			answer = 0
			ask_loop = 0
			
	return answer