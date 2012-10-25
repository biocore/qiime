GHC = ghc

all: FlowgramAli_4frame

FlowgramAli_4frame: FlowgramAli_4frame.lhs
	ghc --make -O2 FlowgramAli_4frame
	cp FlowgramAli_4frame ../../../../scripts/
# On ghc7+ consider using -with-rtsopts="-H100M" for increased speed

clean:
	rm *.o
	rm *.hi
	rm FlowgramAli_4frame
	rm ../../../../scripts/FlowgramAli_4frame
