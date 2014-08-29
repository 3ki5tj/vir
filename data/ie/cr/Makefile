
all: syncin crtrD6samp.dat

crtrD6samp.dat: crD6n4.dat crD6n5.dat crD6n6.dat crD6n7.dat
	./ozcrtrgsl -M8192 -D6 --crtr=$@

syncin:
	rsync -avz cz1@he.host.utmb.edu:lwork/vir/ie/cr/crD6n*.dat .
	rsync -avz cz1@he.host.utmb.edu:lwork/vir/ie/cr/yrtrD6n*.dat .
	rsync -avz cz1@he.host.utmb.edu:lwork/vir/ie/cr/crtr*.dat .

