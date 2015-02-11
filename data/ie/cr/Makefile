
all: syncin crtrD6samp.dat crtrD3samp.dat

figdir = ../../../doc/iedoc/fig

crtrD6samp.dat: crD6n4.dat crD6n5.dat crD6n6.dat crD6n7.dat syncin
	./ozcrtrgsl -M8192 -D6 --crtr=$@ --yr
	make -C $(figdir) -B cmpcorl.eps
	make -C $(figdir) -B cmpcorl.pdf

crtrD3samp.dat: crD3n4.dat crD3n5.dat crD3n6.dat crD3n7.dat crD3n8.dat
	./ozcrtrfftw_q -D3 --crtr=$@ --yr

syncin:
	rsync -avz cz1@he.host.utmb.edu:lwork/vir1/ie/cr/crD6n*.dat .
	rsync -avz cz1@he.host.utmb.edu:lwork/vir1/ie/cr/yrtrD6n*.dat .
	rsync -avz cz1@he.host.utmb.edu:lwork/vir1/ie/cr/crD3n*.dat .
	rsync -avz cz1@he.host.utmb.edu:lwork/vir1/ie/cr/yrtrD3n*.dat .
	rsync -avz cz1@he.host.utmb.edu:lwork/vir1/ie/cr/crtr*.dat .

