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