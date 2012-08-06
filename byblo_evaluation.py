#!/usr/bin/env python -c
#
# Byblo Thesaurus Comparison Tool
#
# Author: Joanne Robert <joanne254@gmail.com>
#         Hamish Morgan <hamish.morgan@sussex.ac.uk>
# Version: 2.0 (2nd August 2012)
#

'''
Utility for comparing thesauri produced by Byblo, in terms of the degree of
overlap over their respective neighbours lists.
'''

__author__='Joanne Robert, Hamish Morgan'
__version__='1.1.0'
__copyright__ = ''
__license__ = ''

import time
import string
import argparse
import math
import logging
import copy
import random
import sys

#################################################################

def main():	
	# 
	# print "Thesaurus Comparison Tool v" + __version__

	parser = argparse.ArgumentParser(
		description='Uility for comparing thesauri produced by Byblo, in terms of the degree of overlap between their respective neighbours lists.',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	# verbose option
	verbosity = parser.add_mutually_exclusive_group()
	verbosity.add_argument('-v', '--verbose', dest='loglevel', 
		action='store_const', const=logging.INFO, default=logging.WARN, 
		help='print extra information to console')
	verbosity.add_argument('-vv', '--very-verbose', dest='loglevel', 
		action='store_const', const=logging.DEBUG, default=logging.WARN, 
		help='print LOTS of debugging information to console')
		
	# file names
	parser.add_argument('files', metavar='file', 
		type=argparse.FileType('r'), nargs=2, 
		action='store', help='thesauri files')
	
	# comparison method
	
	parser.add_argument('-w', dest='weighting', 
		choices=__WEIGHTINGS__.keys(), default='score',
		help='weighting scheme used to encode the importance of elements in each neighbours list')
	parser.add_argument('-k', '--max-rank', 
		type=int, dest='k', action='store', default=None, 
		help='Maximum expected rank (for \'rank\' weighting scheme')

	parser.add_argument('-m', dest='measure', 
		choices=__MEASURES__.keys(), 
		default='cosine',
		help='measure used to compute neighbour list similarity')
		
	parser.add_argument('--test', action='store_true', dest='runtest')	
		
		
	parser.
	args = parser.parse_args()
	
	# Configure logging infrastructure
	logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s')
	log = logging.getLogger()
	log.setLevel(args.loglevel)

	start_time = time.time()
	
	if args.runtest:
		
		test(args.files[0], log=log)
		
	else:
		sim_score2 = thesaurus_file_similarity(
			args.files[0], args.files[1], 
			measure=__MEASURES__[args.measure], 
			weighting=__WEIGHTINGS__[args.weighting], 
			log=log, maxrank=args.k)	
		print sim_score2
	
	end_time = time.time()
	
	if log.isEnabledFor(logging.INFO):
		log.info("All done in %f seconds." % (end_time - start_time))
	

#################################################################



# ======================================
# Neighbours List Similarity Measures
# ======================================


def cosine(list1, list2):
	'''
	Compute the cosine of the angle between the vectors defined by the 
	sparse list of tuples. This is the measure origionally proposed by Lin
	and Weeds.
	'''
	dotproduct = 0
	sqnorm1 , sqnorm2 = 0 , 0	
	i , j = 0 , 0
	while i < len(list1) and j < len(list2):				
		if list1[i][0] < list2[j][0]:
			sqnorm1 += list1[i][1] ** 2
			i += 1
		elif list1[i][0] > list2[j][0]:
			sqnorm2 += list2[j][1] ** 2
			j += 1			
		else: # list1[i][0] == list2[j][0]:
			dotproduct += list1[i][1] * list2[j][1]
			sqnorm1 += list1[i][1] ** 2
			sqnorm2 += list2[j][1] ** 2
			i += 1
			j += 1
	while i < len(list1):
		sqnorm1 += list1[i][1] ** 2
		i += 1
	while j < len(list2):
		sqnorm2 += list2[j][1] ** 2
		j += 1
	if sqnorm1 <= 0 or sqnorm2 <= 0:
		return 0
	return dotproduct / math.sqrt(sqnorm1 * sqnorm2)


def jaccard(list1, list2):
	'''
	Weight Jaccard similarity measure, calculates the similarity between 
	lists as the intersection over the union of multisets.
	'''
	shared = 0
	union = 0	
	i , j = 0 , 0
	while i < len(list1) and j < len(list2):				
		if list1[i][0] < list2[j][0]:
			union += list1[i][1]
			i += 1
		elif list1[i][0] > list2[j][0]:
			union += list2[j][1]
			j += 1			
		else: # list1[i][0] == list2[j][0]:
			shared += min(list1[i][1], list2[j][1])
			union += max(list1[i][1], list2[j][1])
			i += 1
			j += 1
	while i < len(list1):
		union += list1[i][1]
		i += 1
	while j < len(list2):
		union += list2[j][1]
		j += 1
	if union == 0:
		return 0
	return 1.0 * shared / union


def average_precision(list1, list2):
	'''
	Compute the similarity between list1 and list2 as the precision averaged
	across recall deltas as k-cutoff increases. 
	'''
	olist1 = sorted(list1, key=lambda x: -x[1])
	E = {e for e,_ in list2}
	total = 0
	for k in xrange(0,len(olist1)):
		if len(E.intersection(olist1[k])) == 1:
			sublist1 = entry_order(olist1[:k+1])
			total += precision(sublist1, list2)
	return total  / len(list2)


def precision(list1, list2):
	'''
	Precision measures, computes similarity as the proportion of entries in
	list1 that are also present in list2, where presence is weighted by 
	neighbours score.
	'''
	shared = 0
	list1sum = 0	
	i , j = 0 , 0
	while i < len(list1) and j < len(list2):				
		if list1[i][0] < list2[j][0]:
			list1sum += list1[i][1]
			i += 1
		elif list1[i][0] > list2[j][0]:
			j += 1			
		else: # list1[i][0] == list2[j][0]:
			shared += list1[i][1]
			list1sum += list1[i][1]
			i += 1
			j += 1
	while i < len(list1):
		list1sum += list1[i][1]
		i += 1
	if list1sum == 0:
		return 0
	return 1.0 * shared / list1sum

def recall(list1, list2):
	'''
	Recall measures, computes similarity as the proportion of entries in
	list2 that are also present in list1, where presence is weighted by 
	neighbours score.
	'''
	return precision(list2, list1)


def fscore(list1, list2, beta=1):
	'''
	F-Beta score is the weighted hardmonic mean average of precision and 
	recall. Harmonic rather than arithmetic since precision and recall are 
	obstensibly rates, and we want the average rate over a fixed period.
	'''
	p = precision(list1, list2)
	r = recall(list1, list2)
	if p + r == 0:
		return 0
	beta_sq = (beta ** 2)
	return (beta_sq + 1) * p * r / (beta_sq * p + r)


# ======================================
# Neighbours List Re-Weighting Schemes
# ======================================


def score(neighs):
	'''
	No-op weighting that simply uses the similarity score
	'''
	return neighs


def rank(neighs):
	'''
	Replace score with rank ascending. Not generally useful on it\'s own, but
	a precurser to other weightings schemes.
	'''
	return [(e, r + 1) for r, (e, s) in 
		enumerate(sorted(neighs, key=lambda x: -x[1]))]


def inv_rank(neighs):
	'''
	Replace score with inverse-rank (1 over rank).
	'''
	return [(e, 1.0 / r) for e,r in rank(neighs)]


def k_minus_rank(neighs, k=None):	
	'''
	Replace score with rank descending from max size k.
	'''
	if k == None: k = len(neighs)
	assert k >= len(neighs)
	return [(e, (k + 1) - r) for e,r in rank(neighs)]


def binary(neighs):
	'''
	Replace score with 1; equally weight all positions in the neighburs list
	'''
	return [(e, 1 if s > 0 else 0) for e,s in neighs]


# ======================================
# Main thesaurus handling code
# ======================================


def score_order(lst):
	return sorted(lst, key=lambda x: -x[1])

def entry_order(lst):
 	return sorted(lst, key=lambda x: x[0])

def mean(X):
	'''Arithmetic mean of list X'''
	return 1.0 * sum(X) / len(X) \
		if len(X) > 0 else float('NaN')
	
def stddev(X, mu=None, correction=1):
	'''Standard deviation of list X'''
	if mu == None: 
		mu = mean(X)
	return math.sqrt(1.0 * sum([(x - mu) ** 2  for x in X]) / (len(X) - correction)) \
		if len(X) - correction > 0 else float("inf")


def thesaurus_file_similarity(filep1, filep2, 
		measure=cosine, weighting=k_minus_rank, maxrank=None,
		log=logging.getLogger()):
	'''
	'''
	log = log.getChild('thesaurus_similarity')
	log.info('-----------------------------------------')
	
	# Load the thesauri files into neighbours lists
	neighs1 = extract_terms(filep1, log=log)
	neighs2 = extract_terms(filep2, log=log)
	
	log.info('-----------------------------------------')
	
	# Sort in ascending order of base entry (to allow merging)
	log.info('Sorting first neighbours list.')
	neighs1.sort()
	log.info('Sorting second neighbours list.')
	neighs2.sort()
		
	log.info('-----------------------------------------')
	
	# Calculate the similarities
	sims = thesaurus_similarity(neighs1, neighs2, 
		measure=measure, weighting=weighting, maxrank=maxrank, log=log)

	log.info('-----------------------------------------')

	# thesaurus_similarity only returns scores where matchin base entries 
	# where found so need to add misses
	hits = [s for _,s in sims]
	misses = [0 for _ in xrange(1, len(neighs1) + len(neighs2) - (2 * len(hits)))]	
	results = hits + misses

	print hits
	# Print a whole bunch of stats if verbose output is enabled
	if log.isEnabledFor(logging.INFO):
		log.info('Matched  %d / %d base entries' % (len(hits), len(results)))		
		log.info('Similarity Range = [%f, %f]' % (
			min(hits) if len(misses) == 0 else 0,
			max(hits) if len(hits) > 0 else 0))

		log.info('Hits Mean Similarity = %f +/- %f (@95%%)' % (
			mean(hits), stddev(hits) * 1.959963984540054))

		log.info('Mean Similarity = %f +/- %f (@95%%)' % (
			mean(results), stddev(results) * 1.959963984540054))
			
		log.info('Best matches:\n\t' + 
			string.join(['%s => %f' % (e,s) 
				for e,s in score_order(sims)[0:20]], '\n\t'))
	
	return mean(results)



def thesaurus_similarity(neighs1, neighs2, 
		measure=cosine, weighting=k_minus_rank, maxrank=0,
		log=logging.getLogger()):
	'''
	Produce a list of similarity tuples, where each tuple is a base entry with
	the neighbours list similarity between each thesauri. If a base entry 
	doesn't exist in one or other of the thesauri no score is produced.
	
	Neighbours lists are assumed to have been previous sorted 
	lexicographically in ascending order of entry.
	'''
	log = log.getChild('compute')
	if log.isEnabledFor(logging.INFO):
		log.info('Calculating %s similarities between %s weighted neighbours lists.' % (
			measure.func_name, weighting.func_name))
	sims = []
	i,j = 0,0
	
	while i < len(neighs1) and j < len(neighs2) :
		if neighs1[i][0] == neighs2[j][0]:
			if log.isEnabledFor(logging.DEBUG):	
				log.debug('entry: %s' % neighs1[i][0])			
			
			sim = neighbours_list_similarity(neighs1[i][1], neighs2[j][1],
				measure=measure, weighting=weighting, maxrank=maxrank)
			
			sims.append( (neighs1[i][0], sim) )			
			i += 1
			j += 1
		elif neighs1[i][0] < neighs2[j][0]:
			i += 1
		else: #if th1[i][0] > th2[j][0]:
			j += 1
		if log.isEnabledFor(logging.INFO) and max(i,j) % 1000 == 0: 
			log.info('Calculated %d similarities. (%.1f%% complete)' % (
				len(sims), 100.0 * (i+j) / (len(neighs1)+len(neighs2))))
	if log.isEnabledFor(logging.INFO):
		log.info('Calculated %d similarities. (%.1f%% complete)' % (
			len(sims), 100.0))
	return sims


def neighbours_list_similarity(list1, list2, 
		measure=cosine, weighting=k_minus_rank, maxrank=0,
		log=logging.getLogger()):
			
	kwargs = dict()
	if weighting == k_minus_rank:
		kwargs['k'] = max(maxrank, len(list1), len(list2))

	wlist1 = weighting(sorted(list1), **kwargs)
	wlist2 = weighting(sorted(list2), **kwargs)	
	sim = measure(wlist1, wlist2)
	if log.isEnabledFor(logging.DEBUG):	
		log.debug('reweighting 1 (entry:before=>after): %s' % 
			['%s:%.3f=>%.3f' % (e,x,y) for (e,x),(_,y) in zip(sorted(list1), wlist1)])
		log.debug('reweighting 2 (entry:before=>after): %s' % 
			['%s:%.3f=>%.3f' % (e,x,y) for (e,x),(_,y) in zip(sorted(list2), wlist2)])
	return sim


def extract_terms(filep, log=logging.getLogger()):
	'''
	Read all entry neighbours lists from the given file, returning a list 
	containing one element for each base entry. Each element consists of pair; 
	the base_entry string, and another list of neighbour/score tuples.	
	'''
	log = log.getChild('load')
	if log.isEnabledFor(logging.INFO): 
		log.info('Reading thesaurus \'%s\'.' % filep.name)
	terms = []
	base_entry_count = 0
	neighbour_count = 0
	for line in filep:
		f = string.split(line, "\t")
		base_entry = f[0]
		neighbours = [(f[i], float(f[i+1])) for i in xrange(1, len(f), 2)]
		terms.append( (base_entry, neighbours) )
		base_entry_count += 1
		neighbour_count += len(neighbours)
		if log.isEnabledFor(logging.INFO) and base_entry_count % 20000 == 0:
			log.info('Read %d entries with %d neighbours.' % (
				base_entry_count, neighbour_count))
	if log.isEnabledFor(logging.INFO):
		log.info('Completed: Read %d entries with %d neighbours.' % (
			base_entry_count, neighbour_count))
	filep.close()
	return terms


# ==============================================
# Test various measure and weighting combinations
# ==============================================

def degrade(neighs, N=0, repeats=1):
	for i,(b,nl) in enumerate(neighs):
		neighs[i] = (b,score_order(nl)) 
	for _,lst in neighs:
		for _ in xrange(1,repeats):
			i = random.randint(0, max(len(lst)-1, N))
			j = i + 1
			if j < len(lst):
				(lst[i] , lst[j]) = ((lst[i][0], lst[j][1]), (lst[j][0], lst[i][1]))
			elif j == len(lst):
				lst[i] = (hex(random.randint(1,2**32))[2:], lst[i][1] * 0.9)
	for i,(b,nl) in enumerate(neighs):
		neighs[i] = (b,entry_order(nl)) 



#
#
#	
def test(filep, N=0, repeats=5000, interval=1000, log=logging.getLogger()): 
	'''
	Perform a test to evaluate how a particular measure performs as the 
	thesaurus degrades away from the one provided.
	'''
	# Laod the neighbours lists and sort by score descending
	neighs1 = extract_terms(filep, log=log)
	neighs1 = [ (b,sorted(nl,key=lambda x: -x[1])) for b,nl in neighs1]
	neighs1.sort()

	# 2. Take an exact copy of the neighbours list
	neighs2 = copy.deepcopy(neighs1)
	
	# Print the table headers
	sys.stdout.write('%7s' % '')
	for w in __WEIGHTINGS__.keys():
		for _ in __MEASURES__.keys():	
			sys.stdout.write('%8s' % w[:7])
	print
	sys.stdout.write('%7s' % 'repeat')
	for _ in __WEIGHTINGS__.keys():
		for m in __MEASURES__.keys():	
			sys.stdout.write('%8s' % m[:7])
	print
			
	# 3. For each iteration degrade the nieghbours list slightly by making
	#    random swaps. A swap involves reversing their position and updating 
	#	 the score to reflect the change. 	
	r = 0
	while r < repeats:
		sys.stdout.write('%7d' % (r))
		for w in __WEIGHTINGS__.keys():
			for m in __MEASURES__.keys():	
				sims = thesaurus_similarity(neighs1, neighs2, log=log,
					weighting=__WEIGHTINGS__[w], measure=__MEASURES__[m])
				sys.stdout.write('%8.3f' % mean([s for _,s in sims])),
		print
		degrade(neighs2, N=100, repeats=interval)
		r += interval
	#



#################################################################

__MEASURES__ = {'cosine':cosine, 'jaccard':jaccard, 
		'precision':precision, 'recall':recall, 'fscore':fscore,
		'ap':average_precision}
__WEIGHTINGS__ = {'score':score, 'rank':k_minus_rank, 
		'invrank':inv_rank, 'binary':binary}
	
__DEFAULT_MEASURE__ = cosine
__DEFAULT_WEIGHTING__ = score



if __name__=='__main__':
	main()


#################################################################


	    