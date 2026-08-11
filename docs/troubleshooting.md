# Troubleshooting

## I cannot clone the repository

Check that the repository URL is correct and that you have been added as a collaborator. If the repository is private, sign in to the correct GitHub account. The [GitHub cloning guide](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) and [GitHub Desktop guide](https://docs.github.com/en/desktop/overview/getting-started-with-github-desktop) cover the two supported paths.

## Git asks for a password

GitHub no longer accepts account passwords for Git operations over HTTPS. Use GitHub Desktop or follow GitHub's [authentication guidance](https://docs.github.com/en/authentication).

## My branch is behind

Save your work, switch to the default branch, pull the latest changes, and then return to your branch. Ask a teammate before resolving a conflict if you are unsure which version should remain.

## A command or tool is missing

Record the exact command, your operating system, and the full error message. Ask Copilot to explain the error without sharing credentials or private data. Check the project's Tools and stack section in the README before installing anything.

## My data is too large or sensitive for Git

Do not commit it. Check whether it can be shared at all, record where it is stored and how it was obtained, and add a small README describing the expected file or download step without including the data itself.

