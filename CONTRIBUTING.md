# PAWpySeed Contributions

## Levels of Access

### Direct Repository Access

If you have direct edit access to the repository,
this implies that you are trusted to make frequent changes
to the code without supervision. Please commit all
changes to dev or an issue-specific or user-specific branch.
When you want to merge with master, merge into dev if you
have not already, run the full unit test suite,
and then merge into master.

### Pull Requests

If you do not have direct edit access to the repository,
please submit all changes as pull requests to the dev
branch after thoroughly testing your code. The maintainer
will accept the pull request, run the unit test suite,
and then merge it with master

## Style Standards

Please follow the PEP 8 Style Guide for Python code
(https://www.python.org/dev/peps/pep-0008/) with one major
exception: As long as the code is in early development
(version 0.\*.\*), use tabs instead of spaces
to make the writing process easier. This is subject to change later
Small deviations from these guidelines will be excepted,
but major deviations will result in rejected code.

GNU has good standards for C
(https://www.gnu.org/prep/standards/html_node/Syntactic-Conventions.html#Syntactic-Conventions).
The C code will likely involve a lot of lines of math. In this case, it may be helpful
to break the occasional rule for readability's sake. Using intuition is fine.
I like to group terms in long lines of math (e.g. (5\*x+2) + (a+b+c) or
(5\*x+2) \* (a+b+c)), which results in fewer spaces than generally recommended
but makes the terms clear. Open brackets should be at the end of a line of C,
and closing brackets should have their own line, at the same indentation
level as the condition, for/while, or function header.

## Releases

PAWpySeed will use the semantic version
number system described at https://semver.org.
All releases will be published on the release branch.
Each time a release is published, the maintainer
will publish the release to PyPI so the code can
be installed by pip.
