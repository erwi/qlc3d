# Coding Standards
This file contains coding standards for AI assisted development.

## Guidelines
- Follow widely known best practices for software development.
  - SOLID principle
  - DRY principle
  - KISS principle

## Tests
- All new code must have unit tests.
- do not relax test tolerances to make tests pass without permission from human operator. 
- Write tests in BDD style using GIVEN/WHEN/THEN or ARRANGE/ACT/ASSERT style comments.
- Tests are part of documentation. They should be readable and understandable by humans.

## Error Handling
- The general philosophy is to fail fast and loud. 
  - Don't try to recover from unexpected states, missing arguments, null pointers etc. by falling back to some defaults unless the user explicitly asked for it.
- Use the RUNTIME_ERROR macro for error handling. It automatically fills in the location of the error.

## Comments and Documentation
- Use Doxygen style comments for documenting code. 
  - Use /** ... */ for function and class documentation.
  - Use // for inline comments.

- Check that existing comments are still valid and update them if necessary.

- On major changes, check the project documentation and update it if necessary. 
  - This includes the user manual, implementation documentation and known-bugs.md.

- Prefer British spelling over American spelling in comments and documentation as well as variable and function naming.