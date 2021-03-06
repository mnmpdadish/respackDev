# file map data
 &filemap
 basename = 'Si',
 number_PP_file = 1/
 ps-Si ps-Si.ichr

# symmetry data
 &SYMMETRY
 SYMMETRY_FORMAT = 'reciprocal',
 NUMBER_SYM_OP = 24/
    1  0  0    0  1  0    0  0  1     0  0  0 ! rg(3,3), pg(3)
    0  1  0    0  0  1    1  0  0     0  0  0
    0  0  1    1  0  0    0  1  0     0  0  0
    1  0  0    0  0  1    0  1  0     0  0  0
    0  0  1    0  1  0    1  0  0     0  0  0
    0  1  0    1  0  0    0  0  1     0  0  0
   -1 -1 -1    0  1  0    0  0  1     0  0  0
   -1 -1 -1    0  0  1    0  1  0     0  0  0
   -1 -1 -1    1  0  0    0  0  1     0  0  0
   -1 -1 -1    0  0  1    1  0  0     0  0  0
   -1 -1 -1    1  0  0    0  1  0     0  0  0
   -1 -1 -1    0  1  0    1  0  0     0  0  0
    0  1  0   -1 -1 -1    0  0  1     0  0  0
    0  0  1   -1 -1 -1    0  1  0     0  0  0
    1  0  0   -1 -1 -1    0  0  1     0  0  0
    0  0  1   -1 -1 -1    1  0  0     0  0  0
    1  0  0   -1 -1 -1    0  1  0     0  0  0
    0  1  0   -1 -1 -1    1  0  0     0  0  0
    0  1  0    0  0  1   -1 -1 -1     0  0  0
    0  0  1    0  1  0   -1 -1 -1     0  0  0
    1  0  0    0  0  1   -1 -1 -1     0  0  0
    0  0  1    1  0  0   -1 -1 -1     0  0  0
    1  0  0    0  1  0   -1 -1 -1     0  0  0
    0  1  0    1  0  0   -1 -1 -1     0  0  0

# atom data
 4 14 
 1  0.000  0.000  0.000  
 1  0.250  0.250  0.250  

# k-points data
&smpl_kpt
dos_mode = 'COS',
dos_mesh = 15, 15, 15,
bz_mesh = 12, 
bz_number_tile = 1/ 
 6 6 6 
 2 2 2  

# main data
&tappinput
lattice_factor = 10.318 
LATTICE_LIST = 0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 0.5, 0.0, 0.5,
cutoff_wave_function = 6.0,
number_element = 1,
number_atom = 2,
number_band = 50,
store_wfn = 1,
initial_lpt=1, 
control_uptime = 3600,
SCF_CONVERGE = 1.0E-015,
xc_type = 'PBE'/ 

# struct_opt data
&struct_opt
number_cycle = 0/

# str_opt_constr data
 1
 0
 
# trace band data
 &trace_band
 number_band = 50,
 number_trace_block = 4,   
 /
  'L'   'G'   'X'   'W'   'L'    
 0.500 0.000 0.500 0.500 0.500 
 0.500 0.000 0.000 0.250 0.500 
 0.500 0.000 0.500 0.750 0.500 
      30    30    30    30    
