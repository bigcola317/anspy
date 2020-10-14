#include "ANStoolkit.cpp"

extern "C"
{
	extern int* getBinarySymbolSpread(float prob, int bits, int lnL) {
		ANS test=init_binary(prob, bits);
		test.quantize_prec(lnL);
		test.spread_fast();
		int* symbol_spread = new int[1<<lnL];
		for (int i=0; i<1<<lnL; i++) symbol_spread[i] = test.s[i];
		return symbol_spread;
	}

    extern int add(int a, int b) {
		return a+b;
	}
}
