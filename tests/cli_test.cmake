# Command line interface test for gcrank-build.
#
# Invoked by CTest as:
#   cmake -DCLI=<path to gcrank-build> -DDATA=<tests/data> -DTMP=<scratch> -P cli_test.cmake
#
foreach(required CLI DATA TMP)
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
expect_matches(version "gcrank-build 1.0.0")

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
