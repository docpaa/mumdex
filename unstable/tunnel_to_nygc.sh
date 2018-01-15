#! /bin/bash

echo do locally: ssh -p 2229 -R 2222:localhost:22 andrewsp-488@127.0.0.1

ssh -L 2229:mumdel.nygenome.org:22 andmac

