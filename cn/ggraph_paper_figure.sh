#! /bin/bash

#
# ggraph_paper_figure.sh
#
# Allows a user to reproduce the specific figure used in the ggraph paper
# 
# Assumes that mumdex code is compiled and available in ~/mumdex,
# the hg19 genome is located at ~/hg19/chrAll.fa,
# that gene and cytoband definition files are in ~/hg19/chrAll.fa.bin/ 
# and that the sample files {m,f,d,s}.txt are in the current directory.
# The gene definition, cytoband definition and sample files are available at 
# http://mumdex.com/ggraph/
#

abspos_left=2983775000
(cd ~/mumdex ; make -j > /dev/null) &&
(
echo 'Note: to exactly reproduce the figure from the paper:'
echo 'you need to run on a mac, to get the exact fonts used in the paper'
echo 'you need to manually turn on cytoband display (near top left)'
echo 'you need to increase point display size by two levels (near bottom right)'
echo 'have the mouse focus remain in the window, in the central graphing region'
echo 'take a window-screenshot: command-shift-4 space, then click on window'
) &&
~/mumdex/ggraph --geometry 1800x800+0+0 --initial-view $abspos_left $((abspos_left+1000000)) 0 2.5 cn ~/analysis/mums/hg19/chrAll.fa abspos,ratio,seg {m,f,d,s}.txt



