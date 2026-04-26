# Instructions
- Modify existing or add new unit tests under cpp-tests when changing or adding new functionality
- After making code changes always run all unit tests under cpp-tests and fix issues until those tests pass
- Do not relax test tolerances to make tests pass.
- Ask clarifying questions rather than making assumption if something is unclear or ambiguous or a choice needs to be made between multiple options
- Fail fast and loud. Dont try to recover from unexpected states by falling back to some defaults unless the user explicitly asked for it. 
- When failing use the RUNTIME_ERROR macro which auto-fills in the location of the error. 

## Documentation
- Update the README file and any other documentation so that they accurately describe the current state after any changes. Do not refer to old state or previous implementations.
- Add Doxygen usage documentation to header files using /** ... */ and including tags like @param and @return. Implementation focused longer internal comments belong in the .cpp files.
- Technical implementation documentation can be found in the doc-impl subdirectory organised by topic in separate markdown files. Keep this up to date.

## Testing
- Document tests you write using comments so its easy to understand what is being tested. 
  - use either ARRANGE/ACT/ASSERT or GIVEN/WHEN/THEN style comment-sections in the tests
- Run the tests from the build/tests directory so that the resource-relative paths resolve correctly. If you run the tests from the project root, the resources will not be found and many tests will fail.
- The suite is very verbose, so run it with output captured to a log file so you can verify the exit code and the final test summary cleanly. 
    - for running tests, use something like `cd /home/eero/projects/qlc3d/build/tests && ./cpp-test > /tmp/cpp-test.log 2>&1; status=$?; echo EXIT:$status; tail -n 40 /tmp/cpp-test.log`
    - update the above path if it changes
  