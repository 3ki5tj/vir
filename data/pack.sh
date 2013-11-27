#!/bin/bash

# pack data

bkdir="cz1@129.109.88.204:/Bossman/cz1/vir"

bk_data="$bkdir/data"
bk_web="$bkdir/web"
bk_doc="$bkdir/doc"
bk_prog="$bkdir/prog"

rsync -avz ../web/hsvirial/index.html $bk_web/
rsync -avz ../web/hsvirial/css $bk_web/
rsync -avz ../web/hsvirial/js $bk_web/
rsync -avzL ../web/hsvirial/virsampjava.png $bk_web/
make -C ../web/hsvirial/prog/java jar
rsync -avzL ../web/hsvirial/VirSampApp.jar $bk_web/
make -C ../web/hsvirial/prog/c zip
rsync -avzL ../web/hsvirial/vircode.zip $bk_web/

make -C ../doc zip
rsync -avz ../doc/virdoc.zip $bk_doc/

target_mc="mcdata"
target_py="pydata"

fns="Bring.dat Z.dat "
fns+="stampede/D2 "
fns+="stampede/D3 "
fns+="stampede/D4 "
fns+="stampede/D5 "
fns+="stampede/D6 "
fns+="stampede/D7 "
fns+="stampede/D8 "
fns+="stampede/benchD20n64 stampede/benchD20r1n64 "
fns+="stampede/D9n8 stampede/D9r1n20 stampede/D9r1n64Z "
fns+="stampede/D10n8 "
fns+="stampede/D11n8 stampede/D11r1n32 stampede/D11r1n64Z "
fns+="stampede/D12n8 "
fns+="stampede/D13n8 stampede/D13r1n64 stampede/D11r1n64Z "
fns+="stampede/D14n8 "
fns+="stampede/D15r1n8 stampede/D15r1n64 "
fns+="stampede/D16r1n8 "
fns+="stampede/D17r1n8 stampede/D17r1n64 "
fns+="stampede/D18r1n8 "
fns+="stampede/D19r1n8 stampede/D19r1n64 "
fns+="stampede/D20r1n8 stampede/D20r1n64 "
fns+="stampede/D25r2n8 stampede/D25r1n64 stampede/D25r2n64 "
fns+="stampede/D30r2n8 stampede/D30r2n64 "
fns+="stampede/D35r2n8 stampede/D35r2n64 "
fns+="stampede/D40r3n8 stampede/D40r3n64 "
fns+="stampede/D45r3n8 stampede/D45r3n64 "
fns+="stampede/D50r3n8 stampede/D50r3n64 "
fns+="stampede/D60r3n8 stampede/D60r3n64 "
fns+="stampede/D70r4n8 stampede/D70r4n64 "
fns+="stampede/D80r4n8 stampede/D80r4n64 "
fns+="stampede/D90r4n8 stampede/D90r4n64 "
fns+="stampede/D100r4n8 stampede/D100r4n64 "

fns+="kraken/D2r1n20Z "
fns+="kraken/D3r1n16Z "
fns+="kraken/D4r1n12Z "
fns+="kraken/D5r1n16Z "
fns+="kraken/D6r1n20Z "
fns+="kraken/D7r1n12Z "
fns+="kraken/D7r1n24Z "
fns+="kraken/D8r1n28Z "
fns+="kraken/D9r1n32Z "

for fn in $fns; do
  rsync -avz --exclude="*D*n*.o*" --exclude="fb*" --exclude="MTSEED" \
    --exclude="a.out" --exclude="a.mic" --exclude="*.bak1" \
    --exclude="ZrD*.dat" --exclude="ZrD*.dat[0-9]+" --exclude="run[0-9]*" \
    --exclude="*.tmp*" --exclude="fbnr.dat" --exclude="bad" \
    --exclude="_*" --exclude="checkn*" --exclude="bak*" --exclude="*.bdb" \
    --delete-after \
    $fn $target_mc/
done
zip -r mcdata.zip $target_mc
rsync -avz mcdata.zip $bk_data/

fns="pyhs/tconv*.txt pyhs/vir*.txt pyhs/Xt*.txt pyhs/Zt*.txt"
for fn in $fns; do
  rsync -avz $fn $target_py/
done
zip -r pydata.zip $target_py
rsync -avz pydata.zip $bk_data/

