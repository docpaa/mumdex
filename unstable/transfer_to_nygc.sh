#! /bin/bash

scp -pr $1 andmac:transfer/$1

ssh -A andmac scp -r transfer/$1 andrewsp-488@mumdel.nygenome.org:$1


