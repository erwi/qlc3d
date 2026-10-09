# Coding Standards
This file contains coding standards for AI assisted development.

## High Level Guidelines
- Follow widely known best practices for software development.
  - SOLID principle
  - DRY principle
  - KISS principle
- keep code coupling low and well modularized
- code should be easily testable

## Tests
- All new code must have unit tests.
- do not relax test tolerances to make tests pass without permission from human operator. 
- Write tests in BDD style using GIVEN/WHEN/THEN comments interleaved (BDD comment location relative to test code is important for readability) with the code to make it clear what is being tested and how:
```
GINVEN some initial state
<test set-up code>
WHEN some action is performed
<test action code>
THEN some expected result should occur
<test assertion code>
```
- Do not use the SCENARIO/GIVEN/WHEN/THEN from the catch2 framework. Use only the GIVEN/WHEN/THEN comments as described above.

## Error Handling
- The general philosophy is to fail fast and loud. 
  - Don't try to recover from unexpected states, missing arguments, null pointers etc. by falling back to some defaults unless the user explicitly asked for it.
- Use the RUNTIME_ERROR macro for error handling. It automatically fills in the location of the error.

## Comments and Documentation
- Use Doxygen style comments for documenting code.
  - All methods, function, classes, structs, enums etc. in .h files should have descriptive comments
  - functions in .cpp files that don't have declaration in .h file should have descriptive comments
  - Use /** ... */ for function and class documentation.
  - Use // for inline comments.

- Check that existing comments are still valid and update them if necessary.

- On major changes, check the project documentation and update it if necessary. 
  - This includes the user manual, implementation documentation and known-bugs.md.

- Prefer British spelling over American spelling in comments and documentation as well as variable and function naming.

## Rules
- Use namespaces for all functions except for the main() function and test code.
  - no indentation for namespace blocks, but indent everything inside the namespace block.
- Use PascalCase for class names and camelCase for functions, variable names and objects.
- Use braces for all if, else, for, while and do statements even if they are one line.
- use of pointers is discouraged. Use references instead. If you must use pointers, use smart pointers instead of raw pointers to indicate ownership.
  - Raw pointers can be passed to functions to indicate that the function does not take ownership of the object. In this case, use a comment to indicate that the function does not take ownership of the object.
- Const correctness is important, use it.
- For simple conditionals, check before calling. Validate the condition at the call site and only invoke the function when the preconditions are satisfied. Do not hide a simple guard inside the function and have it return early; instead, write the call as an explicit conditional so the function name reads like the action being performed. Example: `if (isReady()) { start(); }` is preferred over `start();` plus a nested `if (!isReady()) return;` inside `start()`.