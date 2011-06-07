GHC = ghc

all: FlowgramAli_4frame

FlowgramAli_4frame: FlowgramAli_4frame.lhs
	ghc --make -O2 FlowgramAli_4frame
# On ghc7+ consider using -with-rtsopts="-H100M" for increased speed

install:
	cp FlowgramAli_4frame ../bin/
clean:
	rm *.o
	rm *.hi
	rm FlowgramAli_4frame
