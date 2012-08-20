#!/usr/bin/env python
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
__copyright__ = 'TODO'
__license__ = 'TODO'

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
		description=
'''Uility for comparing thesauri produced by Byblo, in terms of the degree 
of overlap between their respective neighbours lists.''',
		formatter_class=argparse.RawDescriptionHelpFormatter)

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
		help='Maximum expected rank (for \'rank\' weighting scheme)')

	parser.add_argument('-m', dest='measure', 
		choices=__MEASURES__.keys(), 
		default='cosine',
		help='measure used to compute neighbour list similarity')

	parser.add_argument('--version', action='version', 
		version=__version__)
	
	parser.epilog = 'Measures: \n\n  '
	parser.epilog += '\n  '.join(['%s %s' % (m,f.__doc__) for m,f in __MEASURES__.items()])
	parser.epilog += '\nWeightings: \n\n  '
	parser.epilog += '\n  '.join(['%s %s' % (m,f.__doc__) for m,f in __WEIGHTINGS__.items()])
		
		
	args = parser.parse_args()
	
	# Configure logging infrastructure
	logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s')
	log = logging.getLogger()
	log.setLevel(args.loglevel)

	start_time = time.time()

	(mu, sigma) = thesaurus_file_similarity(
		args.files[0], args.files[1], 
		measure=__MEASURES__[args.measure], 
		weighting=__WEIGHTINGS__[args.weighting], 
		log=log, maxrank=args.k)	
	print 'Similarity = %f +/- %f (@95%%)' % (mu, sigma * 1.959963984540054)

	end_time = time.time()
	
	if log.isEnabledFor(logging.INFO):
		log.info("All done in %f seconds." % (end_time - start_time))


#################################################################



# ======================================
# Neighbours List Similarity Measures
# ======================================


def cosine(list1, list2):
	'''
	Compute the cosine of the angle between the vectors. This is the measure 
	origionally proposed by Lin and Weeds.
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
	normaliser =  math.sqrt(sqnorm1 * sqnorm2) \
		if sqnorm1 > 0 and sqnorm2 > 0 else 0
	
	return (dotproduct, normaliser)


def jaccard(list1, list2):
	'''
	Weighted-Jaccard similarity measure, calculates the similarity between 
	lists as the intersection over the union of multisets. For traditional 
	Jaccard use binary weighting.
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
	return (shared, union)


def average_precision(list1, list2):
	'''
	Average Precision computes the similarity between two neighbours lists 
	as the precision averaged across all recall values, for all head sublist 
	of the first nieghbours. Can be thought of as the area under the 
	precision / recall curve for all ranked cut-offs.
	'''
	list1_ranked = score_order(list1)
	list2_entries = {e for e,_ in list2}
	total = 0
	for k,_ in enumerate(list1_ranked):
		if list1_ranked[k][0] in list2_entries:
			total += precision(entry_order(list1_ranked[:k+1]), list2)
	return (total, len(list2))


def precision(list1, list2):
	'''
	Precision measures, computes similarity as the (weighted) proportion of
 	entries in first neighbours list that are present in second. 
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
	return (shared, list1sum)

def recall(list1, list2):
	'''
	Recall measures, computes similarity as the (weighted) proportion of 
	entries in the second neighbours list that are also present in first.
	'''
	return precision(list2, list1)


def fscore(list1, list2, beta=1):
	'''
	F-Score is the hardmonic mean average of precision and recall. It 
	provides a conservative over-all performance from thoses statistics.
	'''
	p = precision(list1, list2)
	r = recall(list1, list2)

	p = 1.0 * p[0] / p[1] if p[1] else 0
	r = 1.0 * r[0] / r[1] if r[1] else 0
	
	beta_sq = (beta ** 2)
	return ((beta_sq + 1) * p * r / (beta_sq * p + r), 1) \
		if p + r else (0,1)



def random_measure(list1, list2):
	'''
	A measure that simply returns a uniform random number in the range [0,1]
	'''
	return (random.random(),1)
	


# ======================================
# Neighbours List Re-Weighting Schemes
# ======================================


def score(neighs):
	'''
	Uses the similarity score from the neighbours list as the weighting.
	(The method origionally proposed by Dekan Lin.)
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
	Replace score with the inverse-rank (1 over rank); i.e the largest score 
	element will have weight 1/1, the next 1/2, then 1/3, 1/4 etc...
	'''
	return [(e, 1.0 / r) for e,r in rank(neighs)]


def random_weighting(neighs):
	'''
	Replace score with a random number sampled uniformly from the range [0,1]
	'''
	return [(e,random.random()) for e,_ in neighs]


def k_minus_rank(neighs, k=None):	
	'''
	Replace score with rank descending from max size (k); i.e the i^th larges 
	score neighbour, in a list of size k, will be given weighting of k - i.	
	'''
	if k == None: k = len(neighs)
	assert k >= len(neighs)
	return [(e, (k + 1) - r) for e,r in rank(neighs)]


def binary(neighs):
	'''
	Replace score with 1 and implicitely everything else takes score 0; has
	the effect of equally weighting all positions in the neighburs list.
	'''
	return [(e, 1 if s > 0 else 0) for e,s in neighs]



# ======================================
# Aggregation methods
# ======================================

from numpy import prod
from numpy import exp
from numpy import log


def ratio_mean(X, p=1, micro=False):
	'''
	Mean that supports arbitrary power generalisation for both micro and macro
	averaging a list of ratios.
	
	Means of power <= 0 are undefined for list containing zero elements. As a
	work around zero elements are ignored
	'''
	if p == float('-inf'):
		return min(1.0 * n / d for n,d in X if d)
	elif p == float('inf'):
		return max(1.0 * n / d for n,d in X if d)
	elif p == 0:  # lim p->0 (geometric mean)
		if micro: # micro-average
			denom = sum(d for _,d in X if d>0)
			if denom == 0: return 0
			nom =  sum(d * log(1.0 * n / d) for n,d in X if d>0 and n>0)
			mu = exp( nom  / denom ) 
			corr = float(sum(1 for n,d in X if d>0 and n>0)) / len(X)
			return (corr) * mu
		else: # macro-average
			N = len(X)
			N1 = sum(1 for n,d in X if d>0 and n>0)
			if N1 == 0:
				return 0
			mu = prod([1.0 * n / d for n,d in X if d>0 and n>0]) ** (1.0 / N1)
			corr = float(N1) / N			
			return corr * mu
	elif p < 0: # (including the harmonic mean at p=-1)
		if micro: # micro-average
			denom = float(sum(d for n,d in X if d>0 and n>0))
			if denom == 0: return 0
			nom = sum(float(d) * ((float(n) / float(d)) ** p ) for n,d in X if d>0 and n>0)
			if nom == 0: return 0
			corr = float(sum(1 for n,d in X if d>0 and n>0)) / sum(1 for n,d in X if d>0)
			return (corr) * (float( nom / denom ) ** (1.0 / p))
		else: # macro-average
			denom = sum(1 for n,d in X if d>0 and n>0)
			if denom == 0: return 0
			nom = sum((1.0 * n / d) ** p for n, d in X if d>0 and n>0)
			if nom == 0: return 0
			corr = float(denom) / len(X) 
		 	return corr * (( nom / denom ) ** (1.0 / p))
	else: # for other p: 0 > p > inf 
		if micro: # micro-average
			denom = float(sum(d for n,d in X if d))
			if denom == 0: return 0
			nom = sum(float(d) * ((float(n) / float(d)) ** p ) for n,d in X if d and n)
			if nom == 0: return 0
			mu = float( nom / denom ) ** float(1.0 / float(p))
			corr = float(sum(1 for n,d in X if d>0 and n>0)) / len(X)
			return corr * mu 
		else: # macro-average
			denom = len(X)
			if denom == 0: return 0
			nom = sum((1.0 * n / d) ** p for n, d in X if d and n)
			if nom == 0: return 0
		 	return ( nom / denom ) ** (1.0 / p)


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
	(mu, sigma) = thesaurus_similarity(neighs1, neighs2, 
		measure=measure, weighting=weighting, maxrank=maxrank, log=log)
	
	return (mu, sigma)


def thesaurus_similarity(neighs1, neighs2, 
		measure=cosine, weighting=k_minus_rank, maxrank=None,
		log=logging.getLogger(), power=1, micro=True):
	'''
	'''
	# Calculate the similarities
	sims = thesaurus_similarities(neighs1, neighs2, 
		measure=measure, weighting=weighting, maxrank=maxrank, log=log)
	
	return (ratio_mean([x for e,x in sims], p=power, micro=micro), 0)
	
	# 	
	# # Print a whole bunch of stats if verbose output is enabled
	# if log.isEnabledFor(logging.INFO):
	# 	log.info('Matched  %d / %d base entries' % (len(hits), len(results)))		
	# 	log.info('Similarity Range = [%f, %f]' % (
	# 		min(hits) if len(misses) == 0 else 0,
	# 		max(hits) if len(hits) > 0 else 0))
	# 
	# 	log.info('Hits Mean Similarity = %f +/- %f (@95%%)' % (
	# 		mean(hits), stddev(hits) * 1.959963984540054))
	# 
	# 	log.info('Mean Similarity = %f +/- %f (@95%%)' % (
	# 		mean(results), stddev(results) * 1.959963984540054))
	# 		
	# 	log.info('Best matches:\n\t' + 
	# 		string.join(['%s => %f' % (e,s) 
	# 			for e,s in score_order(sims)[0:20]], '\n\t'))
	# 
	# return (mean(results), stddev(results))




def thesaurus_similarities(neighs1, neighs2, 
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
			sim = neighbours_list_similarity(neighs1[i][1], [],
				measure=measure, weighting=weighting, maxrank=maxrank)
			sims.append( (neighs1[i][0], sim) )	
			i += 1
		else: #if th1[i][0] > th2[j][0]:
			sim = neighbours_list_similarity([], neighs2[j][1],
				measure=measure, weighting=weighting, maxrank=maxrank)			
			sims.append( (neighs2[j][0], sim) )			
			j += 1
		if log.isEnabledFor(logging.INFO) and max(i,j) % 1000 == 0: 
			log.info('Calculated %d similarities. (%.1f%% complete)' % (
				len(sims), 100.0 * (i+j) / (len(neighs1)+len(neighs2))))
	while i < len(neighs1):
		sim = neighbours_list_similarity(neighs1[i][1], [],
			measure=measure, weighting=weighting, maxrank=maxrank)
		sims.append( (neighs1[i][0], sim) )	
		i += 1
	while j < len(neighs2):
		sim = neighbours_list_similarity([], neighs2[j][1],
			measure=measure, weighting=weighting, maxrank=maxrank)			
		sims.append( (neighs2[j][0], sim) )			
		j += 1
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

	# Re-weighting the vectors
	if list1 == list2:
		wlist1 = weighting(sorted(list1), **kwargs)	
		wlist2 = wlist1
	else:
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



#################################################################

__MEASURES__ = {
		'cosine':cosine, 
		'jaccard':jaccard, 
		'precision':precision, 
		'recall':recall, 
		'fscore':fscore,
		'ap':average_precision
		}
__WEIGHTINGS__ = {
		'score':score, 
		'rnk':k_minus_rank,
		'invrank':inv_rank, 
		'binary':binary,
		# 'rnd':random_weighting
		}
	
__DEFAULT_MEASURE__ = cosine
__DEFAULT_WEIGHTING__ = score



if __name__=='__main__':
	main()


#################################################################


	    