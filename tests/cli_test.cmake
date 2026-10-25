# Command line interface test for gcrank-build and gcrank-query.
#
# Invoked by CTest as:
#   cmake -DCLI=<gcrank-build> -DQUERY=<gcrank-query> -DDATA=<tests/data> -DTMP=<scratch> -P cli_test.cmake
#
foreach(required CLI QUERY DATA TMP)
  if(NOT DEFINED ${required})
    message(FATAL_ERROR "cli_test.cmake: -D${required}=... is required")
  endif()
endforeach()

file(REMOVE_RECURSE "${TMP}")
file(MAKE_DIRECTORY "${TMP}")
file(MAKE_DIRECTORY "${TMP}/out")

function(run_cli tag)
  execute_process(
    COMMAND "${CLI}" ${ARGN}
    RESULT_VARIABLE code
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err
    WORKING_DIRECTORY "${TMP}")
  set("${tag}_code" "${code}" PARENT_SCOPE)
  set("${tag}_out" "${out}" PARENT_SCOPE)
  set("${tag}_err" "${err}" PARENT_SCOPE)
endfunction()

function(run_query tag)
  execute_process(
    COMMAND "${QUERY}" ${ARGN}
    RESULT_VARIABLE code
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err
    WORKING_DIRECTORY "${TMP}")
  set("${tag}_code" "${code}" PARENT_SCOPE)
  set("${tag}_out" "${out}" PARENT_SCOPE)
  set("${tag}_err" "${err}" PARENT_SCOPE)
endfunction()

function(expect_code tag want)
  if(NOT "${${tag}_code}" STREQUAL "${want}")
    message(FATAL_ERROR "expected exit code ${want} from '${tag}' but got '${${tag}_code}'
stdout: ${${tag}_out}
stderr: ${${tag}_err}")
  endif()
endfunction()

function(expect_matches tag needle)
  if(NOT "${${tag}_out}${${tag}_err}" MATCHES "${needle}")
    message(FATAL_ERROR "output of '${tag}' does not match '${needle}':
stdout: ${${tag}_out}
stderr: ${${tag}_err}")
  endif()
endfunction()

# --- usage and version ------------------------------------------------------

run_cli(help --help)
expect_code(help 0)
expect_matches(help "Usage:")
expect_matches(help "--overwrite")

run_cli(version --version)
expect_code(version 0)
expect_matches(version "gcrank-build 0.2.0")

run_cli(noargs)
expect_code(noargs 1)
expect_matches(noargs "--fasta is required")

run_cli(unknown --frobnicate)
expect_code(unknown 1)

run_cli(badthreads --fasta x.fa --out out --threads 0)
expect_code(badthreads 1)

run_cli(badnumber --fasta x.fa --out out --chunk-size 12x)
expect_code(badnumber 1)

# --- input errors -----------------------------------------------------------

run_cli(missing --fasta "${TMP}/nope.fa" --out out)
expect_code(missing 2)

# --- a successful build -----------------------------------------------------

run_cli(build --fasta "${DATA}/test_small.fa" --out out)
expect_code(build 0)
set(index "${TMP}/out/test_small.gcrank")
foreach(file manifest.json contigs.tsv genome.bv genome.rb valid.bv valid.rb
             genome.bv.crc32 genome.rb.crc32 valid.bv.crc32 valid.rb.crc32 contigs.tsv.crc32)
  if(NOT EXISTS "${index}/${file}")
    message(FATAL_ERROR "expected ${index}/${file} to exist")
  endif()
endforeach()

# --- refusing to overwrite, then overwriting --------------------------------

run_cli(again --fasta "${DATA}/test_small.fa" --out out)
expect_code(again 3)
expect_matches(again "already exists")

run_cli(overwrite --fasta "${DATA}/test_small.fa" --out out --overwrite)
expect_code(overwrite 0)

# --- an explicit name and verbose output ------------------------------------

run_cli(named --fasta "${DATA}/test_small.fa" --out out --name custom --verbose --threads 2)
expect_code(named 0)
if(NOT EXISTS "${TMP}/out/custom.gcrank/manifest.json")
  message(FATAL_ERROR "--name custom did not produce custom.gcrank")
endif()
expect_matches(named "custom.gcrank")
if(NOT EXISTS "${TMP}/out/test_small.gcrank/manifest.json")
  message(FATAL_ERROR "the earlier index disappeared")
endif()

# --- gzip input -------------------------------------------------------------

run_cli(gzip --fasta "${DATA}/test_small.fa.gz" --out out --name from_gzip)
expect_code(gzip 0)

# --- stdin (non-seekable input) ---------------------------------------------

if(UNIX)
  execute_process(
    COMMAND sh -c "cat '${DATA}/test_small.fa' | '${CLI}' --fasta - --out '${TMP}/out' --name from_stdin"
    RESULT_VARIABLE stdin_code OUTPUT_VARIABLE stdin_out ERROR_VARIABLE stdin_err)
  if(NOT stdin_code EQUAL 0)
    message(FATAL_ERROR "stdin build failed: ${stdin_err}")
  endif()
  foreach(file contigs.tsv genome.bv genome.rb valid.bv valid.rb)
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E compare_files
              "${TMP}/out/from_stdin.gcrank/${file}" "${DATA}/test_small.gcrank/${file}"
      RESULT_VARIABLE same)
    if(NOT same EQUAL 0)
      message(FATAL_ERROR "the stdin build differs from the golden index in ${file}")
    endif()
  endforeach()

  # A gzip stream on stdin is auto-detected by zlib.
  execute_process(
    COMMAND sh -c "cat '${DATA}/test_small.fa.gz' | '${CLI}' --fasta - --out '${TMP}/out' --name from_stdin_gz"
    RESULT_VARIABLE stdingz_code OUTPUT_VARIABLE stdingz_out ERROR_VARIABLE stdingz_err)
  if(NOT stdingz_code EQUAL 0)
    message(FATAL_ERROR "stdin gzip build failed: ${stdingz_err}")
  endif()
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files
            "${TMP}/out/from_stdin_gz.gcrank/genome.bv" "${DATA}/test_small.gcrank/genome.bv"
    RESULT_VARIABLE stdingz_same)
  if(NOT stdingz_same EQUAL 0)
    message(FATAL_ERROR "the stdin gzip build differs from the golden index")
  endif()
endif()

# --- empty (but valid) FASTA with --verbose --------------------------------

file(WRITE "${TMP}/empty.fa" "")
run_cli(empty --fasta "${TMP}/empty.fa" --out out --name empty --verbose)
expect_code(empty 0)
expect_matches(empty "n/a")
if(NOT EXISTS "${TMP}/out/empty.gcrank/manifest.json")
  message(FATAL_ERROR "the empty index was not written")
endif()

# --- contig names must be non-empty and unique ------------------------------

file(WRITE "${TMP}/dup.fa" ">chr1\nACGT\n>chr1\nACGT\n")
run_cli(dup --fasta "${TMP}/dup.fa" --out out --name dup)
expect_code(dup 2)
expect_matches(dup "duplicate contig name")

file(WRITE "${TMP}/noname.fa" ">\nACGT\n")
run_cli(noname --fasta "${TMP}/noname.fa" --out out --name noname)
expect_code(noname 2)
expect_matches(noname "empty contig name")

# --- a truncated gzip is an input error, not a shorter genome ---------------

if(UNIX)
  file(SIZE "${DATA}/test_small.fa.gz" gz_size)
  math(EXPR gz_half "${gz_size} / 2")
  execute_process(
    COMMAND dd if=${DATA}/test_small.fa.gz of=${TMP}/truncated.fa.gz bs=1 count=${gz_half}
    RESULT_VARIABLE dd_code OUTPUT_QUIET ERROR_QUIET)
  if(NOT dd_code EQUAL 0)
    message(FATAL_ERROR "dd failed to create the truncated gzip fixture")
  endif()
  run_cli(truncated --fasta "${TMP}/truncated.fa.gz" --out out --name truncated)
  expect_code(truncated 2)
  expect_matches(truncated "unexpected end of file")
endif()

# --- the CLI reproduces the committed golden index --------------------------

file(MAKE_DIRECTORY "${TMP}/golden")
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E env "SOURCE_DATE_EPOCH=1700000000"
          "${CLI}" --fasta test_small.fa --out "${TMP}/golden" --name test_small
  RESULT_VARIABLE golden_code
  OUTPUT_VARIABLE golden_out
  ERROR_VARIABLE golden_err
  WORKING_DIRECTORY "${DATA}")
if(NOT golden_code EQUAL 0)
  message(FATAL_ERROR "golden build failed: ${golden_err}")
endif()
foreach(file manifest.json contigs.tsv genome.bv genome.rb valid.bv valid.rb)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E compare_files
            "${TMP}/golden/test_small.gcrank/${file}" "${DATA}/test_small.gcrank/${file}"
    RESULT_VARIABLE same)
  if(NOT same EQUAL 0)
    message(FATAL_ERROR "the CLI build differs from the golden index in ${file}")
  endif()
endforeach()

file(REMOVE_RECURSE "${TMP}")
message(STATUS "gcrank-build CLI test passed")

# --- gcrank-query -----------------------------------------------------------

file(MAKE_DIRECTORY "${TMP}/q")
file(WRITE "${TMP}/q/tiny.fa" ">chr1\nACGTACGTNN\n>chr2\nGGCCGGCC\n")
run_cli(qbuild --fasta "${TMP}/q/tiny.fa" --out "${TMP}/q" --name tiny)
expect_code(qbuild 0)

file(WRITE "${TMP}/q/good.bed"
  "chr1\t0\t8\nchr1\t2\t6\nchr1\t8\t10\nchr2\t0\t8\n# comment\ntrack name=x\nbrowser position chr1\n")

run_query(bed --index "${TMP}/q/tiny.gcrank" --bed "${TMP}/q/good.bed")
expect_code(bed 0)
expect_matches(bed "#chrom\tstart\tend\tpct_gc\tnum_gc\tnum_valid\tseq_len")
expect_matches(bed "chr1\t0\t8\t0.500000\t4\t8\t8")
expect_matches(bed "chr1\t2\t6\t0.500000\t2\t4\t4")
expect_matches(bed "chr2\t0\t8\t1.000000\t8\t8\t8")
expect_matches(bed "chr1\t8\t10\t0.000000\t0\t0\t2")

run_query(amb --index "${TMP}/q/tiny.gcrank" --bed "${TMP}/q/good.bed" --remove-ambiguous)
expect_code(amb 0)
expect_matches(amb "chr1\t8\t10\tnan\t0\t0\t2")
expect_matches(amb "chr1\t0\t8\t0.500000\t4\t8\t8")

run_query(tile --index "${TMP}/q/tiny.gcrank" --window 4 --step 4)
expect_code(tile 0)
expect_matches(tile "chr1\t0\t4\t0.500000\t2\t4\t4")
expect_matches(tile "chr2\t4\t8\t1.000000\t4\t4\t4")

file(WRITE "${TMP}/q/bad.bed" "chr1\t0\t12\n")
run_query(badlen --index "${TMP}/q/tiny.gcrank" --bed "${TMP}/q/bad.bed")
expect_code(badlen 2)
expect_matches(badlen "leaves contig")

file(WRITE "${TMP}/q/nope.bed" "nope\t0\t1\n")
run_query(nope --index "${TMP}/q/tiny.gcrank" --bed "${TMP}/q/nope.bed")
expect_code(nope 2)
expect_matches(nope "unknown contig")

run_query(noidx --bed "${TMP}/q/good.bed")
expect_code(noidx 1)
expect_matches(noidx "--index is required")

run_query(both --index "${TMP}/q/tiny.gcrank" --bed "${TMP}/q/good.bed" --window 4)
expect_code(both 1)
expect_matches(both "mutually exclusive")

run_query(missing --index "${TMP}/q/nope.gcrank" --bed "${TMP}/q/good.bed")
expect_code(missing 2)

run_query(qversion --version)
expect_code(qversion 0)
expect_matches(qversion "gcrank-query 0.2.0")

file(REMOVE_RECURSE "${TMP}")
message(STATUS "gcrank-query CLI test passed")
