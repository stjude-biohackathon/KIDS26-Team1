# Team Lead Checklist

Use this page to get the team moving. It is intentionally short: a three-day project needs enough structure to coordinate work, not a second project to maintain.

## Start Here

1. **Access your assigned team repository.** Organizers will provide the repository and add team members. Confirm that you can open it on GitHub and clone it to your computer. Follow [Git and GitHub basics](../docs/git-github-basics.md) for the full path.
2. **Complete the project profile.** Agree on the question, inputs, expected output, tools, and team roles before pursuing a large implementation.
3. **Make a small first change.** Create a branch, update this README or document a data source, commit the change, and open a pull request. Use the terminal steps in [Git and GitHub basics](../docs/git-github-basics.md) or [GitHub Desktop](https://desktop.github.com/) if you prefer a graphical interface.
4. **Ask for help early.** Record blockers in an issue, raise them at a check-in, or ask a mentor. See [troubleshooting](../docs/troubleshooting.md) for common recovery steps.

## The Team

- **Members and roles:** Record these in [team.md](team.md).
- **Ways of working:** Use the [team lead checklist](CHECKLIST.md) to agree on branches, reviews, communication, and check-ins.
- **Current plan:** Keep small first tasks and risks in [project-plan.md](project-plan.md).

## Project Structure

Use the folders that fit your project. You do not need to fill every folder. The following is just a suggestion, yours might look different.

```text
data/raw/           Original inputs; do not edit in place
data/processed/     Cleaned or transformed data
models/             Models, predictions, or model notes
project-management/ Team plan, roles, decisions, and check-ins
src/                Reusable code, organized by purpose
docs/               Optional learning and troubleshooting guides
assets/             Images or other supporting project assets
```

## Data, meta data and secrets
**Do not commit passwords, API keys, private information, or identifiable human or clinical data. Check the source and license before sharing external data or media.**

## Resources

- New to Git or GitHub or need to know how to work with git in a shared repo: read [Git and GitHub basics](../docs/git-github-basics.md).
- Using Copilot agents: read [AI assistance](../docs/ai-guidance.md).
- Stuck during setup: open [troubleshooting](../docs/troubleshooting.md).
- Collaborating on changes: see [Contributing to your team](#contributing-to-your-team) below.

## Contributing to your team
### A Simple Workflow

1. Pick a small task or write down a blocker.
2. Create a branch with a clear name, such as `add-project-profile` or `fix-data-path`.
3. Make one focused change and commit it with a short message.
4. Push the branch and open a pull request.
5. Ask another teammate to look at the change before merging.
6. Update the README or project notes when the change affects how someone uses the project.

The [Git and GitHub basics](../docs/git-github-basics.md) guide explains each step, including a GitHub Desktop workflow.

### A Pull Request Is Ready When

- The change has a clear purpose.
- A teammate can understand what changed.
- You have recorded how you checked it, or explained why checking was not possible.
- Relevant assumptions, data sources, and limitations are documented.
- The change does not include credentials or sensitive data.

Small, incomplete pull requests are welcome when they make the current state visible and clearly describe what remains.


## Reproducibility and Attribution

Make work easier to inspect and reuse by keeping inputs, decisions, methods, and limitations visible. Prefer small readable steps over unexplained one-off commands. Cite data, code, models, and external resources that your project depends on. These practices help the next person understand what happened.

This repository is a reusable template. See [LICENSE.md](../LICENSE.md) for the licensing terms and update the project profile and attribution when you create a team project.

## Code of Conduct

In this repository, we use St. Jude's Code of Conduct document, outlining our expectations for all participants. Ask for help early, give feedback about the work rather than the person, and make room for different levels of experience.

For details, please visit: https://issuu.com/sjcrh/docs/st._jude_code_of_conduct.


## Team Leads: Before the Event

- [ ] Complete the [project profile](../README.md#project-profile).
- [ ] Agree on one communication channel and a short check-in rhythm.
- [ ] Create three to six small first tasks in the project board or [project-plan.md](project-plan.md).
- [ ] Use the plan and the expected output to suggest practical roles in [team.md](team.md).


## During the Three Days

### Day 1: Align and start

- Confirm the question, problem, or opportunity.
- Confirm the inputs and expected output.
- Make sure everyone can clone the repository and make a small change.
- Agree on branch, commit, and review habits.

### Day 2: Build and learn

- Keep tasks small enough to finish or review in one sitting.
- Record decisions that change the approach in a decision log (create `decisions.md` if useful).
- Document data sources, assumptions, and unexpected limitations as they appear.
- Check in briefly and redistribute work when someone is blocked.

### Day 3: Explain and hand off

- Decide what the final demo and booth must show as a team.
- Make the main workflow understandable to someone who was not in the room.
- Capture what worked, what did not, and what should happen next.
- Run the available checks and record their results.

# Final Output and Handoff

Use this space for the material that helps someone understand the project after the event.

- **Final demo or report:** [Add a link]
- **Main result:** [Summarize what the team built or learned]
- **How to reproduce or run it:** [Link to the relevant instructions]
- **Data and source notes:** [Link to provenance, citations, and licenses]
- **Known limitations:** [What should a reader be careful about?]
- **Next steps:** [What would be worth doing next?]

Keep generated figures and reports clearly named. Do not commit sensitive data or files that cannot be redistributed.

## Communications

Keep communication easy to find and easy to use during the three-day event.

- **Primary channel:** [Add the team channel or meeting link]
- **Slack team channel:** [Add the team slack channel]
- **Team lead:** [Add name and GitHub handle]
- **Mentor or support contact:** [Add contact]
- **Check-in time:** [Add a recurring time or agree in the team plan]
- **Slack general channel:** [Use this general channel for communication to all teams](https://stjudebiohackathon.slack.com/archives/C04JD4M3TCM)

Use `project-management/check-in.md` for short updates when useful (create the file if needed). Do not store private contact details or sensitive project information in this public repository.

## Optional Templates

- [Team and roles](team.md)
- [Project plan](project-plan.md)

Use only the templates that help. The repository should make progress easier, not require paperwork for its own sake.
