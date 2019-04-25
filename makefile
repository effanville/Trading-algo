hellomake: main.o  StockDecision.o BuySell.o dayencoding.o quicksort.o
	g++ -o	ver1	main.o  StockDecision.o BuySell.o dayencoding.o quicksort.o
	
main.o: main.cc
	g++ -O -c -std=c++11 main.cc 

quicksort.o: quicksort.cc quicksort.h
	g++ -O -c -std=c++11 quicksort.cc

StockDecision.o: StockDecision.cc StockDecision.h
	g++ -O -c -std=c++11 StockDecision.cc
	
BuySell.o: BuySell.cc BuySell.h
	g++ -O -c -std=c++11 BuySell.cc
	
dayencoding.o: dayencoding.cc dayencoding.h
	g++ -O -c -std=c++11 dayencoding.cc
	

