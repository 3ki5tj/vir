host=129.109.88.204

refresh: sync
	./ievirextra.py -D2 --lamc -w
	./ievirextra.py -D3 --lamc -w
	./ievirextra.py -D4 --lamc -w
	./ievirextra.py -D5 --lamc -w
	./ievirextra.py -D6 --lamc -w
	./ievirextra.py -D7 --lamc -w
	./ievirextra.py -D8 --lamc -w

sync:
	rsync -avz --exclude="DHT*.dat" --exclude="snapshot*" \
	  cz1@$(host):lwork/vir1/ie/lamc/*PYl*.dat .
