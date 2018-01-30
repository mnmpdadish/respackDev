# file map data
 &filemap
 basename = 'La2CuO4',
 number_PP_file = 3/
 ps-Cu ps-Cu.ichr
 ps-La ps-La.ichr
 ps-O  ps-O.ichr

# symmetry data
 &SYMMETRY
 SYMMETRY_FORMAT = 'reciprocal',
 NUMBER_SYM_OP = 16/
  1  0  0  0  1  0  0  0  1  0  0  0
 -1  0  0  0 -1  0  0  0 -1  0  0  0
  1  0 -1  1  0  0  1 -1  0  0  0  0
 -1  0  1 -1  0  0 -1  1  0  0  0  0
  0  1 -1  1  0 -1  0  0 -1  0  0  0
  0 -1  1 -1  0  1  0  0  1  0  0  0
  0  1  0  0  1 -1 -1  1  0  0  0  0
  0 -1  0  0 -1  1  1 -1  0  0  0  0
  0 -1  1  0 -1  0  1 -1  0  0  0  0
  0  1 -1  0  1  0 -1  1  0  0  0  0
 -1  0  1  0 -1  1  0  0  1  0  0  0
  1  0 -1  0  1 -1  0  0 -1  0  0  0
 -1  0  0 -1  0  1 -1  1  0  0  0  0
  1  0  0  1  0 -1  1 -1  0  0  0  0
  0 -1  0 -1  0  0  0  0 -1  0  0  0
  0  1  0  1  0  0  0  0  1  0  0  0

# atom data
11 29 
 9 57 
 6  8
 2      0.361006958   0.361006958   0.000000000
 2      0.638993042   0.638993042   0.000000000
 1      0.000000000   0.000000000   0.000000000
 3      0.500000000   0.000000000   0.500000000
 3      0.000000000   0.500000000   0.500000000
 3      0.185459862   0.185459862   0.000000000
 3      0.814540138   0.814540138   0.000000000

# k-points data
&smpl_kpt
dos_mode = 'METHFESSEL_PAXTON', 
bz_mesh = 12,
bz_number_tile = 1/
  6  6  6    
  2  2  2  

# main data
&tappinput
lattice_factor=1.889726878,  
lattice_list=-1.89085,1.89085,6.62435,1.89085,-1.89085,6.62435,1.89085,1.89085,-6.62435, 
cutoff_wave_function = 9.0,
number_element = 3,
number_atom = 7,
number_band = 60,
store_wfn = 1,
initial_lpt = 1,
xc_type = 'PBE',
control_uptime = 72000.0, 
elec_kbt=0.01, 
SCF_CONVERGE = 1.0E-015/

# struct_opt data
&struct_opt
number_cycle = 0/

# str_opt_constr data
 1
 0

# trace band data
 &trace_band
 number_band = 60,
 number_trace_block = 5/
  'G'    'X'    'P'     'N'    'G'    'Z'  # nkpt(nbk+1)
 0.0d0  0.0d0  0.25d0  0.0d0  0.0d0  0.5d0 # ak(1,nbk+1)
 0.0d0  0.0d0  0.25d0  0.5d0  0.0d0  0.5d0 # ak(2,nbk+1)
 0.0d0  0.5d0  0.25d0  0.0d0  0.0d0 -0.5d0 # ak(3,nbk+1)
      20     20      20     20     20      # nkfi(nbk)
