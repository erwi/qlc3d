# Instructions
- Modify existing or add new unit tests under cpp-tests when changing or adding new functionality
- After making code changes 
  - always run all unit tests under cpp-tests and fix issues until those tests pass
  - use the code review skill to inspect all changed code and fix any issues found
  - finally run a full build with current build configuration from a clean so that all files are built and linked correctly.
- Ask clarifying questions rather than making assumption if something is unclear or ambiguous or a choice needs to be made between multiple options

## Documentation
- Glossary of terminology in this project can be found in GLOSSARY.md. 
  - Use these consistently in code and documentation.
- Implementation focused documentation is located in the doc-impl subdirectory. This describes the code from a developers point of view.
  - First read the INDEX.md file to find relevant files only
  - Read the relevant files first before looking at code. It will help identify the correct code files and saves tokens. 
- A user-facing documentation or user manual is found in the qlcd/doc/README.md file
  - this should be kept up to date with the current state of the code. If you make changes to the code, update this file as well.
  - especially this relates to settings file keys and their default values.
- Update all relevant documentation files as part of any other changes.
  - All documentation should accurately reflect the current state. Delete old documentation related to past implementations. Don't keep a history of changes. Only current state matters.
  - When you encounter a bug, add a note to known-bugs.md with a file:line citation.
    - Include date of discovery as well as a short severity description.
    - If you fix a bug, remove the note from known-bugs.md.
- Add Doxygen usage documentation to header files using /** ... */ and including tags like @param and @return. Implementation focused longer internal comments belong in the .cpp files.
  - Add Doxygen comments to all new code 
  - Existing function or class changes: check that the Doxygen comments exist and accurately represent the current state
  - Add this at least to enums, classes, functions, structs


## Testing
- Run the tests from the build/tests directory so that the resource-relative paths resolve correctly. If you run the tests from the project root, the resources will not be found and many tests will fail.
- The suite is very verbose, so run it with output captured to a log file so you can verify the exit code and the final test summary cleanly. 
    - for running tests, use something like `cd /home/eero/projects/qlc3d/build/tests && ./cpp-test > /tmp/cpp-test.log 2>&1; status=$?; echo EXIT:$status; tail -n 40 /tmp/cpp-test.log`
    - update the above path in this file if it changes
  