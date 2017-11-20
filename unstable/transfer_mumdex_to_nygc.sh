#! /bin/bash

ssh -A andmac ssh andrewsp-488@mumdel.nygenome.org '"if [ -e mumdex ] ; then mv mumdex oldmumdex/'$(date +%s)'; fi"'

(cd ~/mumdex ; make zip)

scp ~/mumdex.zip andmac:

ssh -A andmac scp mumdex.zip andrewsp-488@mumdel.nygenome.org:

ssh -Att andmac ssh -tt andrewsp-488@mumdel.nygenome.org '"unzip mumdex.zip"'

exit 0

ssh -Att andmac ssh -tt andrewsp-488@mumdel.nygenome.org '"unzip mumdex.zip && cd mumdex && export SGE_ROOT=/opt/sge && make -j 3"'

echo done transfer and compile

(cd ~/mumdex ; make -j > /dev/null)

