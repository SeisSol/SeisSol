# Contributing to SeisSol
First off, thanks for taking the time to contribute! The following is a set of guidelines for contributing to SeisSol. We sincerely ask you to read and follow our [**Code of Conduct**](CODE_OF_CONDUCT.md).

## Contributing as a user
### Reporting issues

A great way to contribute to SeisSol is to send a detailed report when you encounter a problem. Sufficiently detailed reports are highly appreciate. A good report should allow us reproducing the issue. 

Go [here](https://github.com/SeisSol/SeisSol/issues) and click on "New issue". Select either "Bug report" or "Feature request" and click on the corresponding "Get started" button. A new web page will pop up where you will need to fill in a template and give a descriptive title. Once it is done, click on "Submit new issue".

## Contributing as a developer
### Workflow

### Step 1
You need a local fork of the project. Please, go to our [main GitHub page](https://github.com/SeisSol/SeisSol) and press the “fork” button in GitHub. This will create a copy of the SeisSol repository in your own GitHub account.

### Step 2
Clone the forked SeisSol project from GitHub to your PC or laptop:
```
$ git clone --recurse-submodules https://github.com/<your github account>/SeisSol.git
```

Let's open the new project’s directory:
```
$ cd SeisSol
```

At this point, your local copy of the SeisSol project has a single reference to a remote repository i.e., the one that you've just forked in the previous step.

```
$ git remote -v
origin https://github.com/<your_github_account>/SeisSol.git (fetch)
origin https://github.com/<your_github_account>/SeisSol.git (push)
```

You need to set up a reference to the original remote SeisSol repository (referred to as `upstream`) to be able to grab new changes from the SeisSol master branch. It will allow you to synchronize your contribution with us. 
```
$ git remote add upstream https://github.com/SeisSol/SeisSol.git
$ git remote -v
origin https://github.com/<your_github_account>/SeisSol.git (fetch)
origin https://github.com/<your_github_account>/SeisSol.git (push)
upstream https://github.com/SeisSol/SeisSol.git (fetch)
upstream https://github.com/SeisSol/SeisSol.git (push)
```

### Step 3
We highly recommend cloning the latest master branch of the SeisSol project, creating a new branch out of it with a descriptive name and adding your contribution to SeisSol there.
```
$ git checkout master
$ git pull upstream master
$ git branch <descriptive_branch_name>
$ git checkout <descriptive_branch_name>
```
We also recommend following the following format for your branch names i.e., `<prefix>/<short_name>` where `<prefix>` can be *feature*, *bugfix*, *extension*, etc.

### Step 4
Make a commit once you did some changes in your local copy of SeisSol

```
$ git add <files_to_be_part_of_the_commit>
$ git commit --message <descriptive_message_of_this_particular_commit>
```
Note that code formatting is enforced in src/DynamicRupture and src/tests/DynamicRupture.
If you change files in these folders, you can enforce a suitable formatting by running:

```
.ci/format.sh $(which clang-format) .
```

Push it to your remote repository
```
git push origin <descriptive_branch_name>
```
Now it is time to make a pull request (PR). Open the following web-page
```
https://github.com/<your_github_account>/SeisSol/tree/<descriptive_branch_name>
```
and click on "Compare & pull request". Write a descriptive title and a body of your PR, and click on "Create pull request". 

You can create a `draft pull request` If your contribution is not ready yet but you still want to show it to SeisSol maintainers. In this case, click on a green downward arrow next to "Create pull request", select "Create draft pull request" and proceed.

### Step 5
Once you submit your PR, one or two SeisSol maintainers will review it with you. After that, we may have questions. Please, check back on your PR to keep up with the conversation. Maintainers will start reviewing your contribution if at least the following requirements are fulfilled:

- There is no merge conflict with the latest SeisSol master branch
- All CI tests are passed

### Step 6
Congratulations, your PR is merged! The whole SeisSol community thanks you.
