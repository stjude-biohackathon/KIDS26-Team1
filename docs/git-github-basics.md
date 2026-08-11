# Git and GitHub Basics

This is the short path for participants who are new to Git or GitHub. You only need the parts that help your team work today.

## Official Guides

- [Hello World](https://docs.github.com/en/get-started/start-your-journey/hello-world)
- [GitHub flow](https://docs.github.com/en/get-started/using-github/github-flow)
- [Git basics](https://docs.github.com/en/get-started/git-basics)
- [Cloning a repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
- [About branches](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-branches/about-branches)
- [About pull requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/about-pull-requests)
- [About issues](https://docs.github.com/en/issues/tracking-your-work-with-issues/about-issues)
- [About Projects](https://docs.github.com/en/issues/planning-and-tracking-with-projects/learning-about-projects/about-projects)
- [GitHub Desktop documentation](https://docs.github.com/en/desktop)
- [GitHub Desktop download](https://desktop.github.com/download/)

## The Core Idea

Git records versions of files on your computer. GitHub stores a shared copy of the repository and provides collaboration tools. A branch is a safe workspace for one change. A pull request asks the team to review that change before it becomes part of the shared default branch.

Organizers provide the starting repository for each team. Confirm that you can open your assigned repository on GitHub before cloning it. If you cannot access it, ask an organizer or mentor to check your team permissions.

## Terminal Workflow

Copy the repository URL provided by the organizers, then replace the example URL and branch name below. If your team is TeamX and your name is Joe Do, then do the following:

```bash
git clone https://github.com/stjude-biohackathon/KIDS26-TeamX.git
cd KIDS26-TeamX 
git fetch --all
git pull
# Create a new branch for yourself
git switch -c Joe-Do
# Make your changes, say you made changes to a file file1.txt
git status
git add file1.txt
git commit -m "edited file1.txt"
git push --set-upstream origin Joe-Do
```

Then open GitHub, select **Compare & pull request**, describe the change, and ask a teammate to review it. After the pull request is merged, update your local default branch before starting the next task:

```bash
git switch main
git fetch --all
git pull
```

If your repository uses another default branch name, use that name instead of `main`.

## Issues and Pull Requests

Use an issue when the team needs to record a task, idea, question, bug, or blocker. A useful issue answers:

- What are we trying to do or understand?
- Why does it matter for this project?
- What is the smallest next step?
- Who might help?

Use a pull request for a focused change. A useful description includes:

- What changed and why
- How it was checked
- Important assumptions or data changes
- What remains to do

You do not need a large backlog. Three to six clear first tasks are usually enough for Day 1.

## GitHub Desktop Workflow

GitHub Desktop provides a graphical way to clone repositories, create branches, commit changes, push branches, and open pull requests. Use the [official GitHub Desktop guide](https://docs.github.com/en/desktop) if you prefer buttons and menus to terminal commands. The underlying workflow is the same: make a branch, make a focused change, commit, push, and open a pull request.

## Merge Conflicts

A conflict means two branches changed the same part of a file. Do not guess. Read both versions, ask a teammate if the intended result is unclear, keep the correct content, save the file, and test or inspect it before committing the resolution. The [official conflict guide](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/addressing-merge-conflicts/about-merge-conflicts) provides detailed steps.

## Keep the Project Reproducible

Record where inputs came from, which files were changed, what commands or tools were used, and what limitations remain. Do not commit passwords, API keys, private information, or sensitive human or clinical data. Check data licenses before sharing.
