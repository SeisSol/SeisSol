file(REMOVE_RECURSE
  "CMakeFiles/SeisSol-codegen"
  "src/generated_code/gpulike_subroutine.cpp"
  "src/generated_code/init.cpp"
  "src/generated_code/init.h"
  "src/generated_code/kernel.cpp"
  "src/generated_code/kernel.h"
  "src/generated_code/subroutine.cpp"
  "src/generated_code/subroutine.h"
  "src/generated_code/tensor.cpp"
  "src/generated_code/tensor.h"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/SeisSol-codegen.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
