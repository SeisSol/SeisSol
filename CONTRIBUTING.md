<!--
    SPDX-FileCopyrightText: 2021-2024 SeisSol Group

    SPDX-License-Identifier: BSD-3-Clause
    SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

    SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
-->

# Contributing to SeisSol

Firstly, thank you so much for taking the time to contribute!

The following is a set of guidelines for contributing to SeisSol. We sincerely
ask you to read and follow our [**Code of Conduct**](CODE_OF_CONDUCT.md).

## Contributing as a User

We highly appreciate feedback from our users.

### Reporting issues

Sending us a detailed report when encountering a problem is a great way to
contribute to SeisSol. You are highly encouraged to include as much information
as possible. That will give us a higher chance
to reproduce the problem or maybe enable us to directly track it down.

To raise a new issue, visit [here](https://github.com/SeisSol/SeisSol/issues)
and click on "New issue" there. Next, select either the "Bug report" or "Feature
request" template and click on the corresponding "Get started" button.
You will be given a template to fill in.
Don't forget to give your issue a descriptive title.
Once you are done, click on "Submit new issue".

## Contributing as a Developer

### Workflow

To be able to contribute, you will need to open a Pull Request to the
SeisSol Github repository.

### Step 1

To do so, you will require a local fork of the project in your Github profile.
That can be done by venturing to the
[SeisSol GitHub page](https://github.com/SeisSol/SeisSol)
and pressing the “fork” button there.

For SeisSol-internal developers, you may as well upload a branch to the main
SeisSol repository—that will also automatically enable the tests.
In that case, skip this step.

### Step 2

Clone the forked SeisSol project from GitHub to your PC or laptop via

```bash
git clone --recurse-submodules https://github.com/<your github account>/SeisSol.git
```

(**important:** do not forget the `--recurse-submodules` as the submodules
contain important files for the project. If you did not use it while cloning,
you can do so afterwards by running `get submodule update --init --recursive`)

Next, open the new project directory:

```bash
cd SeisSol
```

At this point, your local copy of the SeisSol project has a single reference to
a remote repository i.e., the one that you've just forked in the previous step.
That can be verified by checking

```bash
$ git remote -v
origin https://github.com/<your_github_account>/SeisSol.git (fetch)
origin https://github.com/<your_github_account>/SeisSol.git (push)
```

You need to set up a reference to the original remote SeisSol repository
(referred to as `upstream` here) to be able to grab new changes from the
SeisSol master branch. It will allow you to synchronize your contribution with us.

```bash
$ git remote add upstream https://github.com/SeisSol/SeisSol.git
$ git remote -v
origin https://github.com/<your_github_account>/SeisSol.git (fetch)
origin https://github.com/<your_github_account>/SeisSol.git (push)
upstream https://github.com/SeisSol/SeisSol.git (fetch)
upstream https://github.com/SeisSol/SeisSol.git (push)
```

For SeisSol-internal developers, you may run (`upstream` and `origin` then conflate)

```bash
git clone --recurse-submodules https://github.com/SeisSol/SeisSol.git
```

### Step 3

We highly recommend cloning the latest master branch of the SeisSol project,
creating a new branch out of it with a descriptive name and adding your
contribution to SeisSol there.

```bash
git checkout master
git pull upstream master
git submodule update --recursive --init
git checkout -b <descriptive_branch_name>
```

We also recommend following the following format for your branch names i.e.,
`<prefix>/<short_name>` where `<prefix>` can be `feature`, `bugfix`, `extension`,
etc.

### Step 4

Once you are done (and happy) with your changes, the next step is to turn them
into a commit. That can be done by running the following commands:

```bash
git add <files_to_be_part_of_the_commit>
git commit --message <descriptive_message_of_this_particular_commit>
```

Note that some folders inside SeisSol adhere to formatting and code standards.

If you change files in these folders, you can enforce a suitable formatting by
running (running `clang-format` in version 18.1.5 in the background;
you can install it e.g. via `pip`):

```bash
.ci/format.sh $(which clang-format) .
```

Similarly, a check for `clang-tidy` is being run. To set it up, go to your build
folder (in our case here, it's called `build` inside your cloned repository),
and run.

```bash
cd build
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON .
```

Then, you can run (in the main folder)

```bash
.ci/tidy.sh ./ build/ -fix
```

In the background, `clang-tidy` in version 18 is called.

Once completed, you can push your changes to your remote repository

```bash
git push origin <descriptive_branch_name>
```

Now it is time to make a pull request (PR). Open the following web-page

```bash
https://github.com/<your_github_account>/SeisSol/tree/<descriptive_branch_name>
```

and click on "Compare & pull request". Write a descriptive title and a body of
your PR, and click on "Create pull request".

You can create a "Draft pull request" If your contribution is not ready yet, but
you still want to show it to SeisSol maintainers. In this case, click on a green
downward arrow next to "Create pull request", select "Create draft pull request"
and proceed.

### Step 5

Once you submit your PR, the SeisSol maintainers will review it. After they are
finished, they may have some comments and/or change requests to your code.
Therefore, please check back on your PR from time to time to keep up with the conversation.

We also have a CI enabled for all branches in the SeisSol repository. If you
submit a contribution from outside, a maintainer will probably clone your branch
into the SeisSol repository to let the tests run.

In particular, the following status checks are mandatory right now:

- There is no merge conflict with the SeisSol `master` branch
- All required CI tests pass. That includes:
  - `clang-format` has been applied (where enabled)
  - `clang-tidy` has been applied (where enabled)
  - All equation systems compile in Debug and Release mode, and the respective
    unit tests pass

There is currently no enforced check for the end-to-end CPU and GPU tests;
however, the SeisSol maintainers may still want them to be respected by your changes.

### Step 6

The SeisSol maintainers approve your PR and add it to the merge queue.
The mandatory status checks (as mentioned in step 5) will need to pass here.

### Step 7

Congratulations, your PR is merged! The whole SeisSol community highly
appreciates your contribution.
