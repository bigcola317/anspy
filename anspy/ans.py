import cffi
from ctypes.util import find_library
import numpy as np
from bitstring import BitArray
import sys
import os
from pathlib import Path
from abc import ABC, abstractmethod
import array


# Base class for any ANS system. Provides base functionality
# from a defined symbol spread to encoding/decoding.
# Extensions provide their own custom symbol spread via the 
# abstract `_getSymbolSpread()` method. Symbols are assumed to be coded
# by their enumeration i.e. mapped to the natural numbers [0, n_symbols).
# Symbol spread length `L`, which determines the precision of the probability
# quantization, is set through `lnL` and given by self.L = 2**lnL
class BaseANS(ABC):

	# Inherited classes should call super().__init__() at the
	# end of their own __init__()
	def __init__(self, lnL):
		self.ffi = cffi.FFI()
		self.ffi.cdef("""
				int* getCustomSymbolSpread(float* prob, int m, int lnL);
				int* getBinarySymbolSpread(float prob, int bits, int lnL);
			""")
		# Find libans shared library
		site_packages = next(p for p in sys.path if 'site-packages' in p)
		lib = next(f for f in os.listdir(site_packages) if 'libans' in f)
		libpath = Path(site_packages) / Path(lib)
		# Dynamically link shared library
		self.C = self.ffi.dlopen(str(libpath.resolve()))
		# Initialize the class
		self.lnL = lnL
		self.L = 2**lnL
		self.__configure()


	# Override in inherited classes, must return the symbol spread
	# as a list of symbols.
	@abstractmethod
	def _getSymbolSpread(self):
		pass


	def __createTables(self, symbol_spread):
		# Get Ls for each symbol s, defined as number of occurrences
		# in symbol spread function
		symbols, Ls = np.unique(np.array(symbol_spread), return_counts=True)
		# Create decoding table: list of pairs (symbol, state)
		# where state is the x-th appearance of symbol s, plus Ls
		appearance_dict = {symbols[i]:Ls[i] for i in range(len(symbols))}
		self.decoding_table = [{'symbol': symbol, 'state': 0} for symbol in symbol_spread]
		for i in range(len(self.decoding_table)):
			symbol = self.decoding_table[i]['symbol']
			self.decoding_table[i]['state'] = appearance_dict[symbol]
			appearance_dict[symbol] += 1
		# print('Decoding table:', self.decoding_table)
		# Create encoding table: 
		# dictionary {symbol: [(state, next_state), ...], ...}
		self.encoding_table = {symbol:None for symbol in symbols}
		for s, symbol in enumerate(symbols):
			transition_pairs = []
			occurrences = [i for i,x in enumerate(symbol_spread) if x==symbol]
			for i in range(Ls[s]):
				transition_pairs.append(
					(Ls[s]+i, self.__denormalize_state(occurrences.pop(0))) )
			self.encoding_table[symbol] = transition_pairs
		# print('Encoding table:', self.encoding_table)

	
	def __configure(self):
		symbol_spread = self._getSymbolSpread()
		# print('spread:', symbol_spread)
		self.__createTables(symbol_spread)


	# Sequence expected as a list of symbols.
	# Encodes the sequence in backward direction so that the decoding
	# can be done in forward direction.
	# Returns the state in normalized form and the bitstream.
	def encode(self, sequence):
		X = 0
		x = self.__denormalize_state(X)
		bits = BitArray()
		for symbol in reversed(sequence):
			x, step_bits = self.__index_encodingTable(x, symbol)
			bits += step_bits
		X = self.__normalize_state(x)
		return X, bits


	# Expects the state in normalized form and the
	# bitstream (as output by `encode`).
	# Returns the decoded sequence as a list of symbols.
	def decode(self, state, bits):
		x = self.__denormalize_state(state)
		sequence = []
		while True:
			symbol, x, done = self.__index_decodingTable(x, bits)
			if done:
				break
			sequence.append(symbol)
		return sequence


	def __normalize_state(self, state):
		return state - self.L


	def __denormalize_state(self, state):
		return state + self.L


	# Returns new state and bits to write.
	# Expects non-normalized state as input.
	def __index_encodingTable(self, state, symbol):
		maxXs = self.encoding_table[symbol][-1][0]
		write_bits = BitArray()
		while state > maxXs:
			write_bits += BitArray(bin=str(state%2))
			state = state // 2 # state >> 1
		for pstate, nstate in self.encoding_table[symbol]:
			if state == pstate:
				return nstate, write_bits


	# Expects non-normalized state as input.
	# `bits` is the bitstream from where new bits are
	# taken to form the actual state.
	# Returns the decoded symbol and next state as from
	# the decoding table. Also returns a flag which indicates
	# if there actually was a symbol to read or if we reached the end.
	def __index_decodingTable(self, state, bits):
		state2 = state
		while state < self.L:
			bit = bits[-1]
			del bits[-1] # Pop last item
			state = state*2 + bit
		done = (state==16 and bits.len==0)
		X = self.__normalize_state(state)
		return self.decoding_table[X]['symbol'],\
				self.decoding_table[X]['state'], done



# Creates an ANS for binary strings of length `bits`
# with probability `prob` of an asserted bit. The precision
# in the probability quantization is given by `lnL`.
class BinaryANS(BaseANS):

	def __init__(self, prob, bits, lnL):
		self.prob = prob
		self.bits = bits
		super().__init__(lnL)


	def _getSymbolSpread(self):
		return self.ffi.unpack(
			self.C.getBinarySymbolSpread(self.prob, self.bits, self.lnL), self.L)



# Wrapper around BaseANS which provides base functionality
# from a defined alphabet (and associated symbol probabilities `prob`)
# to encoding/decoding. `prob` specifies the probabilities of the
# symbols in the order of their enumeration.
class CustomANS(BaseANS):

	def __init__(self, prob, lnL):
		self.prob = prob
		super().__init__(lnL)


	def _getSymbolSpread(self):
		# Convert `prob` to a cdata pointer
		prob = array.array('d', self.prob)
		prob = self.ffi.from_buffer('float[]', prob)
		# Invoke C library function
		return self.ffi.unpack(
			self.C.getCustomSymbolSpread(prob, len(self.prob), self.lnL), self.L)		
