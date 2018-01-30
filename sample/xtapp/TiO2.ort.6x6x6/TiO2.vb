# file map data
 &filemap
 basename = 'TiO2',
 number_PP_file = 2/
 ps-Ti ps-Ti.ichr
 ps-O  ps-O.ichr

# symmetry data
 &SYMMETRY
 SYMMETRY_FORMAT = 'reciprocal',
 denom_trans =   2, 
 NUMBER_SYM_OP = 8/
  1  0  0  0  1  0  0  0  1  0  0  0
 -1  0  0  0 -1  0  0  0 -1  0  0  0
 -1  0  0  0 -1  0  0  0  1  0  0  0
  1  0  0  0  1  0  0  0 -1  0  0  0
  1  0  0  0 -1  0  0  0 -1  1  1  1
 -1  0  0  0  1  0  0  0  1  1  1  1
 -1  0  0  0  1  0  0  0 -1  1  1  1
  1  0  0  0 -1  0  0  0  1  1  1  1

# atom data
  4 22 #Ti 
  6 8  #O 
  1 0.000000   0.000000   0.000000
  1 0.500000   0.500000   0.500000
  2 0.304269   0.304269   0.000000
  2 0.695731   0.695730   0.000000
  2 0.804275   0.195727   0.500000
  2 0.195725   0.804274   0.500000

# k-points data
&smpl_kpt
dos_mode = 'METHFESSEL_PAXTON',  
bz_mesh = 12, 
bz_number_tile = 1/ 
  6 6 6  
  2 2 2  

# main data
&tappinput
lattice_factor = 1.0d0,! 1.889726878d0, 
LATTICE_LIST = 8.679500,  0.d0,        0.d0,        0.d0,        8.679500,  0.d0,       0.d0,        0.d0,       5.591700,
cutoff_wave_function = 8.0, 
number_element = 2,
number_atom = 6,
number_band = 50,
store_wfn = 1,
initial_lpt=1, 
control_uptime = 3600,
SCF_CONVERGE = 1.0E-015,
xc_type = 'PBE',
elec_kbt=0.01/

# struct_opt data
&struct_opt
number_cycle = 0/

# str_opt_constr data
 1
 0
 
# trace band data
&trace_band
number_band = 50,
number_trace_block = 4/
  'M'   'Z'   'G'   'X'   'Y'  # nkpt(nbk+1)
 0.5d0 0.0d0 0.0d0 0.5d0 0.0d0 # ak(1,nbk+1)
 0.0d0 0.0d0 0.0d0 0.0d0 0.5d0 # ak(2,nbk+1)
 0.5d0 0.5d0 0.0d0 0.0d0 0.0d0 # ak(3,nbk+1)
      40    40    40    40    
