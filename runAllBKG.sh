#! /bin/bash


SECONDS=0
# For Tracks tree:
echo "photon+laser background samples"
root -l -b -q process_track_tree_draw_v6.C'("fileListBKG.txt")'
duration=$SECONDS
echo "Total time taken for this process ---- $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
