# Contributing to spectrum_utils

:+1::tada: First off, thanks for taking the time to contribute! :tada::+1:

The following document provides guidelines for contributing to the documentation and the code of spectrum_utils. These are mostly guidelines, not rules. Use your best judgment, and feel free to propose changes to this document in a pull request.
**No contribution is too small!**
Even fixing a simple typo in the documentation is immensely helpful.

## Contributing to the documentation

We use [sphinx](https://www.sphinx-doc.org/en/master/) to generate our documentation and deploy it to this site.
Most of the pages on the site are created from simple text files written in markdown.
There are two exceptions to this:

1. The API documentation is automatically generated from the documentation contained in the code.

2. The vignettes are created from Jupyter notebooks.

### Editing most documents

The easiest way to edit a document is by clicking the "Edit on GitHub" like in the top right hand corner of each page.
You'll be taken to GitHub where you can click on the pencil to edit the document.

You can then make your changes directly on GitHub.
Once you're finished, fill in a description of what you changed and click the "Propose Changes" button.

Alternatively, these documents live in the `docs/src` directory of the repository and can be edited like code.
See [Contributing to the code](#contributing-to-the-code) below for more details on contributing this way.


## Contributing to the code

We welcome contributions to the source code of spectrum_utils---particularly ones that address discussed [issues](https://github.com/bittremieux/spectrum_utils/issues).

Contributions to spectrum_utils follow a standard GitHub contribution workflow:

1. Create your own fork of the spectrum_utils repository on GitHub.

2. Clone your forked spectrum_utils repository to work on locally.

3. Create a new branch with a descriptive name for your changes:

    ```bash
    git checkout -b fix_x
    ```

4. Make your changes (make sure to read below first).

5. Add, commit, and push your changes to your forked repository.

6. On the GitHub page for you forked repository, click "Pull request" to propose adding your changes to spectrum_utils.

7. We'll review, discuss, and help you make any revisions that are required.
If all goes well, your changes will be added to spectrum_utils in the next release!


### Python code style

The spectrum_utils project follows the [PEP 8 guidelines](https://www.python.org/dev/peps/pep-0008/) for Python code style.
More specifically, we use [black](https://black.readthedocs.io/en/stable/) to format and lint Python code in spectrum_utils.

We highly recommend setting up a pre-commit hook for black.
This will run black on all of the Python source files before the changes can be committed.
Because we run black for code linting as part of our tests, setting up this hook can save you from having to revise code formatting.
Take the following steps to set up the pre-commit hook:

1. Verify that black and pre-commit are installed.
If not, you can install them with pip or conda:

    ```bash
    # Using pip
    pip install black pre-commit
    
    # Using conda
    conda -c conda-forge black pre-commit
    ```

2. Navigate to your local copy of the spectrum_utils repository and activate the hook:

    ```bash
    pre-commit install
    ```

One the hook is installed, black will be run before any commit is made.
If a file is changed by black, then you need to `git add` the file again before finished the commit.
 
