cmake_minimum_required (VERSION 2.6)
project (seqan_sandbox_inverse_mapper_demos)

seqan_setup_demos(sandbox_inverse_mapper_demos)  # Add meta target and make executables depend on it.

seqan_add_all_executables(sandbox_inverse_mapper_demos)
