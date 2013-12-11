GHC = ghc

all: FlowgramAli_4frame

FlowgramAli_4frame: FlowgramAli_4frame.lhs
	ghc --make -O2 FlowgramAli_4frame
	cp FlowgramAli_4frame ../../../../scripts/
# On ghc7+ consider using -with-rtsopts="-H100M" for increased speed

clean:
	rm -f *.o
	rm -f *.hi
	rm -f FlowgramAli_4frame
	rm -f ../../../../scripts/FlowgramAli_4frame
