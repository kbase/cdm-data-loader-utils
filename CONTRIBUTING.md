# Contributing to the codebase

:+1: First of all: Thank you for taking the time to contribute!

The following is a set of guidelines for contributing to
cdm-schema. These guidelines are not strict rules.
Use your best judgment, and feel free to propose changes to this document
in a pull request.

## Table Of Contents

- [Contributing to the codebase](#contributing-to-the-codebase)
  - [Table Of Contents](#table-of-contents)
  - [Code of Conduct](#code-of-conduct)
  - [Guidelines for Contributions and Requests](#guidelines-for-contributions-and-requests)
    - [Reporting problems and suggesting changes to with the data model](#reporting-problems-and-suggesting-changes-to-with-the-data-model)
    - [Questions and Discussions](#questions-and-discussions)
    - [Adding new code or bug fixes](#adding-new-code-or-bug-fixes)
  - [GitHub Best Practices](#github-best-practices)
    - [Creating and curating issues](#creating-and-curating-issues)
    - [Pull Requests](#pull-requests)


<a id="code-of-conduct"></a>

## Code of Conduct

The KBase CDM team strives to create a welcoming environment for editors, users and other contributors.
Please carefully read our [Code of Conduct](CODE_OF_CONDUCT.md).

<a id="contributions"></a>

## Guidelines for Contributions and Requests

<a id="reporting-issues"></a>

### Reporting problems and suggesting changes to with the data model

Please use our [Issue Tracker][issues] for any of the following:

- Reporting problems
- Requesting new data loaders

<a id="questions-and-discussions"></a>

### Questions and Discussions

Please use the [Issue Tracker][issues] to ask general questions or contribute to discussions.

<a id="adding-elements"></a>

### Adding new code or bug fixes

Please submit a [Pull Request][pulls] to submit new code or to add a bug fix for an existing issue.

<a id="best-practices"></a>

## GitHub Best Practices

### Creating and curating issues

    - Read ["About Issues"][about-issues]
    - Issues should be focused and actionable
    - Complex issues should be broken down into simpler issues where possible

### Pull Requests

    - Read ["About Pull Requests"][about-pulls]
    - Read ["About Branches"][about-branches]
    - Read [GitHub Pull Requests: 10 Tips to Know](https://blog.mergify.com/github-pull-requests-10-tips-to-know/)
    - Pull Requests (PRs) should be atomic and aim to close a single issue
    - Long running PRs should be avoided where possible
    - PRs should reference issues following standard conventions (e.g. “fixes #123”)
    - Never work on the main branch, always work on an issue/feature branch
    - Always create a PR on a branch to maximize transparency of what you are doing
    - PRs should be reviewed and merged in a timely fashion by the CDM science technical leads
    - PRs that do not pass GitHub actions should never be merged
    - In the case of git conflicts, the contributor should try and resolve the conflict
    - If a PR fails a GitHub action check, the contributor should try and resolve the issue in a timely fashion


[about-branches]: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches
[about-issues]: https://docs.github.com/en/issues/tracking-your-work-with-issues/about-issues
[about-pulls]: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests
[issues]: https://github.com/kbase/cdm-schema/issues/
[pulls]: https://github.com/kbase/cdm-schema/pulls/

We recommend also reading [GitHub Pull Requests: 10 Tips to Know](https://blog.mergify.com/github-pull-requests-10-tips-to-know/)
