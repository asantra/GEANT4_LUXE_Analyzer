#! /bin/bash


SECONDS=0
### For Tracks tree:
echo "w0_3000nm JETI40"
root -l -b -q process_track_tree_draw_v3.C'("/nfs/dust/ilc/user/oborysov/hics_list/list_root_hics_165gev_w0_3000nm.txt")'
duration=$SECONDS
echo "Total time taken for this process ---- $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

#### hits tree
SECONDS=0
echo "w0_3000nm JETI40 Hits"
root -l -b -q process_hits_tree_draw_v2.C'("/nfs/dust/ilc/user/oborysov/hics_list/list_root_hics_165gev_w0_20000nm.txt")'
duration=$SECONDS
echo "Total time taken for this process ---- $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

exit 1

SECONDS=0
echo "w0_8000nm JETI40"
root -l -b -q process_track_tree_draw_v3.C'("/nfs/dust/ilc/user/oborysov/hics_list/list_root_hics_165gev_w0_8000nm.txt")'
# do some work
duration=$SECONDS
echo "Total time taken for this process ---- $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."



SECONDS=0
echo "w0_8000nm Phase2"
root -l -b -q process_track_tree_draw_v3.C'("/nfs/dust/ilc/user/oborysov/hics_list/list_root_hics_phase2_165gev_w0_8000nm.txt")'
duration=$SECONDS
echo "Total time taken for this process ---- $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
