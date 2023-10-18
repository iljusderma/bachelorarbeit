#%%
"TASEP"
import functools
import numpy as np
import matplotlib.pyplot as plt
import time

class Chain:
	def __init__(self, L, rates):
		self.L = L
		self.alpha = rates[0]
		self.beta = rates[1]
		self.p = rates[2]
		self.q = rates[3]
		self.state = None
		self.hop_counter = 0

	def reset_hop_counter(self):
		self.hop_counter = 0

	def initialize_state(self, ratio=0.5, mode='random'):
		one_num = int(self.L*ratio)
		self.state = np.zeros(self.L)
		self.state[:one_num] = np.ones(one_num)
		if mode == 'random':
			np.random.shuffle(self.state)
		elif mode =='lbulk':
			pass
		elif mode =='rbulk':
			self.state = np.flip(self.state, axis=0)
		else:
			print(mode, 'is unknown. Use ["random", "lbulk", "rbulk"] instead.')
	
	def eject(self):
		if self.state[-1] == 1:
				self.state[-1] = np.random.choice(2, 1, p=[self.beta, 1 - self.beta])[0]
	
	def hop(self):
		# sublattice-parallel update
		# first: (1,2)(3,4)(5,6)...
		# then: (2,3)(4,5)
		# for higher efficiency I tried to loop through occupied sites not all sites
		# generating the random values each is very inefficient
		occupied_sites = np.where(self.state==1)[0]
		random_val = np.random.choice(2, 2*len(occupied_sites), p=[1 - self.p, self.p])
		for index, site in enumerate(occupied_sites):
			if site%2 == 0 and self.state[site + 1] == 0 and random_val[index] == 1:
				# perform hop
				self.state[site] -= 1
				self.state[site + 1] += 1
				if site == self.L/2:
					self.hop_counter += 1
		# last site excluded bc i+1 out of bounds
		occupied_sites = np.where(self.state[:-1]==1)[0]
		random_val = np.random.choice(2, len(occupied_sites), p=[1 - self.p, self.p])
		for index, site in enumerate(occupied_sites):
			if site%2 == 1 and self.state[site + 1] == 0 and random_val[index] == 1:
				# perform hop
				self.state[site] -= 1
				self.state[site + 1] += 1
				if site == self.L/2:
					self.hop_counter += 1

	def inject(self):
		if self.state[0] == 0:
				self.state[0] = np.random.choice(2, 1, p=[1-self.alpha, self.alpha])[0]

def iterate(iterations, chain, current_steps=100):
		density = np.zeros(iterations)
		current = np.zeros(iterations//current_steps) 
		for i in range(iterations):
			chain.eject()
			chain.hop()
			chain.inject()
			# calculate density
			density[i] = np.sum(chain.state)/chain.L
			if (i+1)%current_steps == 0:
				current[i//current_steps] = chain.hop_counter/current_steps
				chain.reset_hop_counter()
		return density, current

def plot_density(density, probs):
	iterations = len(density)
	fig = plt.figure()
	ax = fig.add_subplot(111, title='Evolution of density', xlabel='iterations', ylabel='density')
	ax.plot(np.arange(iterations), density, marker='.', lw=0)
	if probs[0] > probs[1] and probs[1] < 0.5:
		ax.plot(np.full((iterations,), 1-probs[1]), color='red', label='HD')
	if probs[0] < probs[1] and probs[0] < 0.5:
		ax.plot(np.full((iterations,), probs[0]), color='red', label='LD')
	elif probs[0] > 0.5 and probs[1] > 0.5:
		ax.plot(np.full((iterations,), 0.5), color='red', label='MC')
	else:
		pass
	ax.legend()
	plt.show()

def plot_current(current):
	iterations = len(current)*10
	fig = plt.figure()
	ax = fig.add_subplot(111, title='Evolution of current', xlabel='iterations', ylabel='current')
	ax.plot(np.arange(iterations)[::10], current, marker='.', lw=0)
	plt.show()

def main():
	print(__doc__)
	L = 1000
	iterations = 10*1000
	# [alpha, beta, p, q]
	rates = [0.3, 0.6, 0.7, 0]
	chain = Chain(L, rates)
	chain.initialize_state()
	density, current = iterate(iterations, chain)
	#plot_density(density, rates)
	plot_current(current)

if __name__=='__main__':
	main()