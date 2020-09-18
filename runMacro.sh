#! /bin/bash

# For Tracks tree:
root -l -b -q process_track_tree_draw_v2.C'("f_root_t0_list.txt")'

# For Hits tree:
root -l -b -q process_hits_tree_draw_v2.C'("f_root_t0_list.txt")'
