#python 03_prepare_reference.py data/4i3q.cif




#obabel 4i3q.cif -O 4i3q_raw.pdb

#python ../scripts/zenodo_data_2/scripts/trans_rot_4i3q.py \
#    4i3q_raw.pdb \
#    4i3q_kuvek_ref.pdb




#python3 scripts/triangular_lattice_sphere.py -r 20 -s 6 --hemisphere -o sphere.pdb
#The Check: Verify sphere.pdb contains exactly 261 lines (260 ATOM records + 1 END record). This confirms you have the distal-only hemisphere needed for CYP active sites.

