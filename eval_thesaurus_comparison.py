#!/usr/bin/env python
#
# Evaluation script for Byblo Thesaurus Comparison Tool
#
# Author: Hamish Morgan <hamish.morgan@sussex.ac.uk>
# Version: 2.0 (2nd August 2012)
#

'''
'''

__author__='Hamish Morgan'
__version__='1.0.0'
__copyright__ = 'TODO'
__license__ = 'TODO'

import thesaurus_comparison as tc
import time
import string
import argparse
import math
import logging
import copy
from numpy import random
import sys
import numpy as np

#################################################################

def main():	
	# 
	# print "Thesaurus Comparison Tool v" + __version__

	parser = argparse.ArgumentParser(
		description=
'''Evaluate the various thesaurus comparison measures and weighting schemes.''',
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
		type=argparse.FileType('r'), nargs=1, 
		action='store', help='thesauri files')
	
	# comparison method
	

	parser.add_argument('-k', '--max-rank', 
		type=int, dest='k', action='store', default=None, 
		help='Maximum expected rank (for \'rank\' weighting scheme)')
		
	parser.add_argument('--version', action='version', 
		version=__version__)		
		
	args = parser.parse_args()
	
	# Configure logging infrastructure
	logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s')
	log = logging.getLogger()
	log.setLevel(args.loglevel)

	start_time = time.time()
	
	repeats=70
	interval=1
	
	for filep in args.files:
		print 'Running test on \'%s\':' % filep.name
		run_eval(filep, repeats=repeats, interval=interval, log=log, 
			maxrank=args.k)
	
	end_time = time.time()
	
	if log.isEnabledFor(logging.INFO):
		log.info("All done in %f seconds." % (end_time - start_time))
	

#################################################################





#################################################################


POWERS = {
	'p-2': -2,
	'hm': -1, 
	'gm':0,
	'am':1,
	'rms':2,
	'rmq':3
}
AVGTYPE = {
	 'mic':True,
	'mac':False
}

MEASURES = {
		'cosine':tc.cosine, 
		'jaccard':tc.jaccard, 
		'precision':tc.precision, 
		'recall':tc.recall, 
		'fscore':tc.fscore,
		# 'ap':tc.average_precision
		}
WEIGHTINGS = {
		'score':tc.score, 
		'rnk':tc.k_minus_rank,
		'invrank':tc.inv_rank, 
		'binary':tc.binary,
		# 'rnd':tc.random_weighting
		}
			
#
#
#	
def run_eval(filep, repeats=5000, interval=1000, log=logging.getLogger(), maxrank=None): 
	'''
	Perform a test to evaluate how a particular measure performs as the 
	thesaurus degrades away from the one provided.
	'''
	# Laod the neighbours lists and sort by score descending
	neighs1 = tc.extract_terms(filep, log=log)
	log.info('Sorting thesaurus base entries and neighbours.')
	neighs1 = [ (b,sorted(nl,key=lambda x: -x[1])) for b,nl in neighs1]
	neighs1.sort()

	log.info('Coping thesaurus.')
	# 2. Take an exact copy of the neighbours list
	neighs2 = copy.deepcopy(neighs1)

	stats = dict()
	# 3. For each iteration degrade the nieghbours list slightly by making
	#    random swaps. A swap involves reversing their position and updating 
	#	 the score to reflect the change. 	
	r = 0
	while r < repeats:
		stats[r] = dict()
		log.info('Computing similarities for all weighting/measure combinations.')
		for w in WEIGHTINGS.keys():
			stats[r][w] = dict()
			for m in MEASURES.keys():	
				stats[r][w][m] = dict()
				for p in POWERS.keys():
					stats[r][w][m][p] = dict()
					for mi in AVGTYPE.keys():
						(mu,sigma) = tc.thesaurus_similarity(
							neighs1, neighs2, log=log,
							weighting=WEIGHTINGS[w], 
							measure=MEASURES[m],
							maxrank=maxrank, 
							power=POWERS[p],
							micro=AVGTYPE[mi])
						stats[r][w][m][p][mi] =  (mu,sigma)
		log.info('Degrating thesaurus copy for %d repeats' % interval)
		degrade(neighs2, N=100, repeats=interval)
		print 'degredation iterations: %d / %d' % (r, repeats)
		r += interval
		
	wgts = sorted(WEIGHTINGS.keys())
	msrs = sorted(MEASURES.keys())
	reps = sorted([r for r,_ in stats.items()])
	pwrs = sorted(POWERS.keys())
	mics = sorted(AVGTYPE.keys())
	
	#
	# 
	# table0 = [[(r, mu) \
	# 	for m,(mu,st) in sorted(wstat.items()) \
	# 		for w,wstat in sorted(rstat.items())] \
	# 			for r,rstat in sorted(stats.items())]

		# 
		# table0 = [[ wstat for (m, mstat) in sorted(rstat.items()) \
		# 					for (w, wstat) in sorted(mstat.items()) ] \
		# 	for r,rstat in sorted(stats.items())]
		# 	
	
	print '\nComplete results table:'
	
	table0 = [[stats[r][w][m][p][mi] \
		for mi in mics for p in pwrs for m in msrs for w in wgts] for r in reps]
	table0_hrow = reps
	table0_hcol = [ '%s/%s/%s/%s' % (mi[:4],p[:4],m[:4],w[:4]) \
		for mi in mics for p in pwrs  for m in msrs for w in wgts]	
	print_stats_table(table0, hrow=table0_hrow, hcol=table0_hcol)
	plot_stats(table0,  hrow=table0_hrow, hcol=table0_hcol,
		filename=filep.name+'-thesaurus_comparison-results')

	foo = lambda X: (tc.mean(X), tc.stddev(X))
	
	print '\nAverage measure performance table:'
	table1 = [[ foo([ stats[r][w][m][p][mi][0] for mi in mics for p in pwrs for w in wgts]) for m in msrs] for r in reps ]
	table1_hrow = table0_hrow
	table1_hcol = msrs
	print_stats_table(table1, hrow=table1_hrow, hcol=table1_hcol)
	plot_stats(table1, hrow=table1_hrow, hcol=table1_hcol,
		filename=filep.name+'-thesaurus_comparison-results-avg-measure')


	print '\nAverage weighting performance table:'
	table2 = [[ foo([stats[r][w][m][p][mi][0] for mi in mics for p in pwrs for m in msrs]) for w in wgts] for r in reps ]
	table2_hrow = table0_hrow
	table2_hcol = wgts
	# 
	print_stats_table(table2, hrow=table2_hrow, hcol=table2_hcol)
	plot_stats(table2, hrow=table2_hrow, hcol=table2_hcol,
		filename=filep.name+'-thesaurus_comparison-results-avg-weighting')



	print '\nAverage power performance table:'
	table2 = [[ foo([stats[r][w][m][p][mi][0] for mi in mics for m in msrs for w in wgts]) for p in pwrs ] for r in reps ]
	table2_hrow = table0_hrow
	table2_hcol = pwrs
	# 
	print_stats_table(table2, hrow=table2_hrow, hcol=table2_hcol)
	plot_stats(table2, hrow=table2_hrow, hcol=table2_hcol,
		filename=filep.name+'-thesaurus_comparison-results-avg-power')

	


	print '\nAverage avg-type performance table:'
	table2 = [[ foo([stats[r][w][m][p][mi][0] for p in pwrs for m in msrs for w in wgts]) for mi in mics ] for r in reps ]
	table2_hrow = table0_hrow
	table2_hcol = mics
	# 
	print_stats_table(table2, hrow=table2_hrow, hcol=table2_hcol)
	plot_stats(table2, hrow=table2_hrow, hcol=table2_hcol,
		filename=filep.name+'-thesaurus_comparison-results-avg-avgtype')

	

def print_stats_table(table, hrow=None, hcol=None):
	if hrow == None: hrow = [''] * len(table)
	if hcol == None: hcol = [''] * max([len(r) for r in table])
	
	
	# strtab = [['%.2f' % (mu) \
	strtab = [['%.2f~%-.2f ' % (mu,st*1.96) \
		for (mu,st) in row] for row in table]
	hrow = ['%s ' % h for h in hrow]
	hcol = ['%s ' % str(h) for h in hcol]
	
	strtab = [hcol]  + strtab
	
	hrow = [''] + hrow
	
	strtab = [[hrow[i]] + row for (i,row) in enumerate(strtab)]
	
	widths = [[len(c) for c in r] for r in strtab ]
	widths = [max(col) for col in zip(*widths)]
	fmts = ['%'+str(w+1)+'s' for w in widths]
	
	out = ''
	for i,row in enumerate(strtab):
		for i,col in enumerate(row):
			out += fmts[i] % col
		out += '\n'
	print out


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cbook as cbook
import matplotlib.ticker as ticker


def plot_stats(table, hrow=None, hcol=None, filename=None):
	data = [[x for (x,y) in r] for r in table]
	err = [[y for (x,y) in r] for r in table]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	lines = list()
	
	clrs = ['b','g','r','c','m','y','k']
	lns = ['-','--',':', '-.']
	styles = [c+l for c in clrs for l in lns]
	for i,(r,e) in enumerate(zip(zip(*data), zip(*err))):
		lines.append(ax.errorbar(hrow, r,# yerr=e,linewidth=2, elinewidth=0.5, 
			linestyle=lns[(i/len(clrs)) % len(lns)],
			color=clrs[i % len(clrs)])) 
		#, 
	if len(lines) <= 30:
		ax.legend( lines, hcol, ncol=3 )
	ax.set_xlabel('Degredation repeats')
	ax.set_ylabel('Similarity')
	if filename==None:
		plt.show()
	else:
		fig.set_size_inches(8, 5)
		fig.savefig(filename + '.eps', format='eps')
	
	
# 
# def degrade(neighs, N=0, repeats=1):
# 	for i,(b,nl) in enumerate(neighs):
# 		neighs[i] = (b,tc.score_order(nl)) 
# 	for _,lst in neighs:
# 		for _ in xrange(1,repeats):
# 			i = random.randint(0, max(len(lst)-1, N))
# 			j = i + 1
# 			if j < len(lst):
# 				(lst[i] , lst[j]) = ((lst[i][0], lst[j][1]), (lst[j][0], lst[i][1]))
# 			elif j == len(lst):
# 				lst[i] = (hex(random.randint(1,2**32))[2:], lst[i][1] * 0.9)
# 	# Add a new base entry
# 	base = hex(random.randint(1,2**32))[2:]
# 	nei = hex(random.randint(1,2**32))[2:]
# 	sco = random.random()
# 	neighs.append((base, [(nei,sco)]))
# 	for i,(b,nl) in enumerate(neighs):
# 		neighs[i] = (b,tc.entry_order(nl)) 
# 
# 

RAND = random.RandomState(seed=2)

def degrade(neighs, N=0, repeats=1):
	# 
	# for i,(b,nl) in enumerate(neighs):
	# 	neighs[i] = (b,tc.score_order(nl)) 
	
	for _ in xrange(0,repeats):
		for _ in xrange(1,len(neighs)):
		# for _ in xrange(0,len(neighs)):
			lst_idx = RAND.randint(-1, len(neighs)+1)
			if lst_idx < 0:
				# delete a base entry
				i = RAND.randint(0, len(neighs))
				neighs[i:] = neighs[i+1:]
			elif lst_idx >= len(neighs): 
				x = [1]# Add a new base entry			
				base = hex(RAND.randint(1,2**32))[2:]
				nei = hex(RAND.randint(1,2**32))[2:]
				sco = RAND.rand()
				neighs.append((base, [(nei,sco)]))
			else: 				
				lst = tc.score_order(neighs[lst_idx][1])
				i = RAND.randint(0, max(len(lst)-1, N))
				j = i + 1
				if j < len(lst):
					(lst[i] , lst[j]) = ((lst[i][0], lst[j][1]), (lst[j][0], lst[i][1]))
				elif j == len(lst):
					lst[i] = (hex(RAND.randint(1,2**32))[2:], lst[i][1] * 0.9)
				neighs[lst_idx] = (neighs[lst_idx][0], tc.entry_order(lst))
	# 
	# for i,(b,nl) in enumerate(neighs):
	# 	neighs[i] = (b,tc.entry_order(nl)) 




#################################################################


if __name__=='__main__':
	main()


#################################################################
