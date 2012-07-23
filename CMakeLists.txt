cmake_minimum_required (VERSION 2.6)
project (seqan_sandbox_inverse_mapper)

# We have to register our own include directory here, though.
seqan_setup_includes(include SeqAnSandboxInverse_mapper)
# SeqAnCore is added include in ../core/CMakeLists.txt.  We only have to
# register that we want it here, too.  Same for SeqAnExtras.
seqan_make_seqan_available(SeqAnCore)
seqan_make_seqan_available(SeqAnExtras)
seqan_make_seqan_available(SeqAnSandboxEz_seq_io)

seqan_add_all_subdirectories()  # include_directory(*)
