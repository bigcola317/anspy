#include "ANStoolkit.cpp"

extern "C"
{

	// Returns the symbol spread for a length L=2^lnL table, approximating
	// the probability distribution `prob` of an alphabet of size `m`. 
	extern int* getCustomSymbolSpread(prec* prob, int m, int lnL) {
		ANS test=init_from_distribution(prob, m);
		test.quantize_prec(lnL);
		test.spread_fast();
		int* symbol_spread = new int[1<<lnL];
		for (int i=0; i<1<<lnL; i++) symbol_spread[i] = test.s[i];
		return symbol_spread;
	}


	// Returns the symbol spread for a length L=2^lnL table, approximating
	// the probability distribution of a blocked binary alphabet of block
	// size `bits` with probability `prob` of a bit set to 1.
	extern int* getBinarySymbolSpread(float prob, int bits, int lnL) {
		ANS test=init_binary(prob, bits);
		test.quantize_prec(lnL);
		test.spread_fast();
		int* symbol_spread = new int[1<<lnL];
		for (int i=0; i<1<<lnL; i++) symbol_spread[i] = test.s[i];
		return symbol_spread;
	}

}
