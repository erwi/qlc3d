# Instructions
- Modify existing or add new unit tests under cpp-tests when changing or adding new functionality
- After making code changes always run all unit tests under cpp-tests and fix issues until those tests pass
- Do not relax test tolerances to make tests pass.
- Ask clarifying questions rather than making assumption if something is unclear or ambiguous or a choice needs to be made between multiple options
- Code should fail fast and loud. Dont try to recover from unexpected states, missing arguments, null pointers etc. by falling back to some defaults unless the user explicitly asked for it. 
- When failing use the RUNTIME_ERROR macro which auto-fills in the location of the error. 

## Documentation
- Glossary of terminology in this project can be found in GLOSSARY.md. 
  - Use these consistently in code and documentation.
- Implementation focused documentation is located in the doc-impl subdirectory. This describes the code from a developers point of view.
  - First read the INDEX.md file to find relevant files only
  - Read the relevant files first before looking at code. It will help identify the correct code files. 
- A user-facing documentation or user manual is found in the qlcd/doc/README.md file
- Update all relevant documentation files as part of any other changes.
  - All documentation should accurately reflect the current state. Delete old documentation related to past implementations. Dont keep a history of changes. Only current state matters. 
- Add Doxygen usage documentation to header files using /** ... */ and including tags like @param and @return. Implementation focused longer internal comments belong in the .cpp files.
  - Add Doxygen comments to all new code 
  - Existing function or class changes: check that the Doxygen comments exist and accurately represent the current state


## Testing
- Document tests you write using comments so its easy to understand what is being tested. 
  - use either ARRANGE/ACT/ASSERT or GIVEN/WHEN/THEN style comment-sections in the tests
- Run the tests from the build/tests directory so that the resource-relative paths resolve correctly. If you run the tests from the project root, the resources will not be found and many tests will fail.
- The suite is very verbose, so run it with output captured to a log file so you can verify the exit code and the final test summary cleanly. 
    - for running tests, use something like `cd /home/eero/projects/qlc3d/build/tests && ./cpp-test > /tmp/cpp-test.log 2>&1; status=$?; echo EXIT:$status; tail -n 40 /tmp/cpp-test.log`
    - update the above path if it changes
  