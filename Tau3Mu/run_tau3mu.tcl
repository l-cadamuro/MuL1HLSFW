puts "@@@ Starting"
puts "@@@ Opening project"


########### done by hand
# open_project -reset tau3mu_no_inline_closer_elem_class_2
open_project -reset tau3mu_unroll
set_top tau_3mu

open_solution -reset "solution_1"

add_files src/mini_sorter.cpp
add_files src/miniSorter.cpp
# add_files src/tau3mu.cpp

###################################################

add_files -tb tau3mu_tb.cpp

set_part {xc7vx690tffg1927-2}
create_clock -period 4.17 -name default

puts "@@@ C SIM"
csim_design

puts "@@@ C SYNTH"
csynth_design

puts "@@@ C/RTL COSYM"
cosim_design -trace_level all
# export_design -format ip_catalog  -vendor "cern-cms"

exit
