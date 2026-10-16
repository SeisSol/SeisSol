..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _strong-scaling:

Strong scaling series
=====================

``scripts/benchmark/scaling.py`` runs a strong scaling series of SeisSol on a cluster from a single campaign file: a
case, the builds and run-time variants to compare, the node counts, and the machine. It plans the runs, builds
SeisSol, checks every run before submitting it, submits one job per node count, follows the jobs, and attempts again
what did not complete.

It needs Python 3.10 or later, and PyYAML for campaign files and machine profiles in YAML (JSON works without it).
Slurm is supported, and the local machine without a batch system, which runs a job when it is submitted, for trying
a campaign out.

A first campaign
----------------

.. code-block:: yaml

    version: 1
    name: tpv5-strong
    machine:
      profile: supermuc-ng
      slurm: {account: pn12ab}
    case:
      parameter_file: tpv5/parameters.par
      set:
        abortcriteria.endtime: 1.0
        output.wavefieldoutput: 0
      epoch_interval: 0.1
    builds:
      - name: current
        variants: ../SeisSol/scripts/compile/compile-variants.yaml
        preset: quick
    variants:
      - name: default
      - name: two-ranks
        ranks_per_node: 2
    scaling:
      nodes: {min: 4, max: 64, factor: 2}
      repetitions: 3
      min_cells_per_rank: 20000
    job:
      time_limit: "02:00:00"

.. code-block:: bash

    scaling.py plan   tpv5-strong.yaml   # the runs and jobs it would submit
    scaling.py check  tpv5-strong.yaml   # layouts, tools, SeisSol --check
    scaling.py submit tpv5-strong.yaml   # build, check, and submit
    scaling.py status tpv5-strong.yaml   # the state of every run

``scripts/benchmark/examples/campaign.yaml`` is a commented example. Every command can be repeated: ``submit``
submits only the runs that did not complete, and attempts a failed run again up to ``job.max_attempts`` times, as a
new attempt; a run whose check failed is not attempted again. A campaign keeps everything in its directory
(``directory``, by default ``<name>/`` next to the campaign file), and the commands take that directory in place of
the campaign file as well.

The campaign file
-----------------

Paths are relative to the campaign file. Every mapping is checked for the keys it may have, so that a misspelt key is
an error.

``name``
    The name of the campaign, of its directory and of its jobs.
``machine``
    The machine profile: a name of ``scripts/benchmark/machines/`` (e.g. ``supermuc-ng``, ``lumi-g``, ``local``), a
    path, or a mapping with ``profile`` and settings that replace those of the profile, e.g. the account.
``case``
    ``parameter_file``: the parameter file of every run; the runs start in its directory, so that the files it names
    are found as in a run started by hand. ``set``: overrides of every run, as ``--set`` takes them
    (:ref:`parameter-overrides`); the output prefix is set per run, in its directory. ``epoch_interval``: an epoch
    at least every so much simulated time (``EpochInterval``, :ref:`run-summary`). ``env``: environment variables of
    every run. ``cells``: the cells of the mesh, if they should not come from the check.
``builds``
    The executables to compare: each either ``binary``, an executable built before, or ``variants``, a variant-set
    file of ``scripts/compile/compile_variants.py``, with ``preset``, ``toolchain``, ``set`` (its ``--set``
    options) and, for a set that builds several executables, ``select``, the one to take. A build of a variant-set
    file is built anew, in the campaign directory, whenever runs are submitted, so that a series measures the source
    as it is; CMake rebuilds only what changed.
``variants``
    The variants of the runs, crossed with the builds: ``env`` and ``set`` of their own, a layout of their own
    (``ranks_per_node``, ``cpus_per_rank``, ``threads_per_rank``), and ``launcher_args``, appended to the launcher.
    Without variants, there is one, ``default``.
``scaling``
    ``nodes``: a list of node counts, or a range ``{min, max, factor}``; ``min: auto`` takes the smallest node count
    that holds the mesh at ``max_cells_per_node``. ``repetitions`` of every run. The guards ``min_cells_per_rank``
    and ``max_cells_per_node`` leave out the node counts at which a rank would hold too few cells to scale, or a node
    too many to fit; they need the cells of the mesh, which the check finds. ``order``: ``interleaved`` (the
    default) or ``blocked``.
``job``
    ``time_limit`` of a job, ``run_time_limit`` of a run, ``max_attempts`` of a run, and ``check_in_job``
    (whether the job checks every run before running it; on by default).

Machine profiles
----------------

A machine profile says what a node is and how a run starts on it, after a job script that is known to run SeisSol
well:

``scheduler``
    ``slurm`` or ``local``.
``node``
    ``cores``, ``threads_per_core``, ``gpus``, ``memory`` and ``numa_domains`` of a node.
``layout``
    ``ranks_per_node``, ``cpus_per_rank`` (the CPUs of a rank, its workers and its communication thread),
    ``threads_per_rank`` (its OpenMP workers, ``OMP_NUM_THREADS``) and ``gpus_per_rank``.
``slurm``
    ``partition`` (one, or a list of ``{max_nodes, partition}`` by the size of a job), ``account``, ``qos``,
    ``constraint``, ``exclusive``, ``cpus_per_task`` (whether a job asks for the CPUs of a rank) and ``extra``
    sbatch options.
``launcher``
    The command that starts a run, in front of the wrapper and the executable; ``{nodes}``, ``{ranks}``,
    ``{ranks_per_node}``, ``{cpus_per_rank}``, ``{threads_per_rank}`` and ``{gpus_per_rank}`` stand for the values of
    the run, while shell variables such as ``${CPU_BIND}`` are left to the job.
``wrapper``
    A script started per rank in front of the executable, e.g. ``seissol-launch``, which selects the GPU of a rank, or
    one that pins it.
``setup`` and ``env``
    Shell lines that set up a job, e.g. module loads, and the environment of every run.

The runs
--------

Each node count gets a job of its own, which runs all builds, variants and repetitions at that node count, so that
their comparison does not depend on which nodes a job got. The repetitions take turns: the first one of every build
and variant, then the second one, each time in another order. Within a job, each attempt at a run gets a directory,
``runs/<run>/attempt-<k>/``, and a run ID (``SEISSOL_RUN_ID``), which the files of the attempt carry. Before it
runs, the job checks it (``SeisSol --check``, :ref:`parameter-check`); then it runs it, with the environment of its
variant, ``OMP_NUM_THREADS`` after its layout, and its output in ``output/`` of its directory. The attempt leaves
``started``, ``check.json`` and ``check.log``, ``seissol.log``, ``exitcode`` (and ``check-failed`` if its check
failed), and the output of SeisSol, its run summary among it.

``check``, and ``submit`` before it submits, checks what it can before a job waits in a queue: whether the layouts
fit a node, whether the scheduler and the partitions take the jobs, and, with ``SeisSol --check`` of every build and
variant, whether the parameters, the files they name, the mesh and the output directories are fine. That check also
finds the size of the mesh, which the guards of ``scaling`` need. An executable built for the compute nodes may not
run on a login node; then this check is left to the jobs.
