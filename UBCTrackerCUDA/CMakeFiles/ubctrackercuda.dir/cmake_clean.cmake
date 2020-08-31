file(REMOVE_RECURSE
  "libubctrackercuda.pdb"
  "libubctrackercuda.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/ubctrackercuda.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
