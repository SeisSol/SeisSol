..
  SPDX-FileCopyrightText: 2022-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff


.. _pypi_behind_firewall:

Accessing internet behind a firewall
------------------------------------


Many HPC facilities have rather tight security restrictions, which may include a firewall, preventing any outgoing connections.
To be able to access git or PyPi behind a firewall, we can create a reverse SSH tunnel for internet access.
We describe how to proceed in the following.


1. On your local machine in ~/.ssh/config add the following `RemoteForward` line:

::

    Host supermucNG
        ...
        RemoteForward ddddd localhost:8899

where ddddd is an arbitrary port number with 5 digits.
ddddd is the port on the remote server that is forwarded to your port 8899 of your local machine.
(This number should be different from port number used in other RemoteForward entries.)

2. Install proxy.py on your local machine.

::

    pip install --upgrade --user proxy.py

3. Start proxy.py on your local machine. (And keep it running.)


::

    ~/.local/bin/proxy --port 8899 &

4. Login to the HPC cluster (e.g. ``ssh supermucNG``).
Check that you do not get: `Warning: remote port forwarding failed for listen port ddddd`.
In this case you would need to change ddddd to a different port.
Note that the problem might also be you have already an opened connection to the HPC cluster.

Add to your ``~/.bashrc`` file on the HPC cluster:

::

    export http_proxy=http://localhost:ddddd
    export https_proxy=http://localhost:ddddd

where ddddd is your arbitrary port number.

Then pip or git should be reachable. You can e.g. install pip packages with:

::

    pip install <package name>

In addition, you might need to add the ``--no-build-isolation`` flag to the pip command.

For more information, see also this `link <https://doku.lrz.de/faq-installing-your-own-applications-on-supermug-ng-internet-access-from-supermuc-ng-10746066.html>_`.

.. _git_behind_firewall:

Accessing github behind a firewall (outdated)
---------------------------------------------

Warning: This procedure works, but a much simpler procedure is available at :ref:`pypi_behind_firewall`.

Some HPC clusters (e.g. SuperMUC-NG) restricts access to outside sources and thus does not allow connections to https servers.
Nevertheless, GitHub can be used if remote port forwarding is correctly set.
Here, we described the procedure to set up such port forwarding.


1. On your local machine, add to your ~/.ssh/config the following lines:

::

  Host supermucNG
     Hostname skx.supermuc.lrz.de
     User <Your Login>
     RemoteForward ddddd github.com:22

where ddddd is an arbitrary 5-digital port number, smaller than 65535.

2. Use the following command to login onto the HPC cluster:

.. code-block:: bash

  ssh supermucNG

Add the following lines to your ~/.ssh/config (on the HPC cluster):

::

  Host github.com
     HostName localhost
     User git
     Port ddddd

With ddddd the same port number as before.

3. Create SSH key by typing (use a non-empty passphrase, not too long as you will need to type it often)

.. code-block:: bash

  ssh-keygen -t rsa

4. Go to https://github.com/settings/ssh, add a new SSH key, and paste the public SSH key you just created (the content of ~/.ssh/id_rsa.pub on the HPC cluster).

5. To allow cloning using ssh on SuperMUC-NG with the https address of git repository, add to ``~/.gitconfig``:

::

    [url "git@github.com:"]
        insteadOf = https://github.com/


6. (Optionnally) To avoid having to Enter your passphrase many times during the session, you can execute on the HPC cluster:

.. code-block:: bash

    eval `ssh-agent -s`
    ssh-add ~/.ssh/id_rsa


You should now be able to clone any GitHub repository, e.g. SeisSol including the submodules using:

.. code-block:: bash

  git clone --recursive https://github.com/SeisSol/SeisSol.git


If it works, you will see several lines, for example:

::

  Cloning into 'SeisSol'...
  remote: Enumerating objects: 25806, done.
  remote: Counting objects: 100% (4435/4435), done.
  remote: Compressing objects: 100% (1820/1820), done.
  remote: Total 25806 (delta 2972), reused 3710 (delta 2551), pack-reused 21371
  Receiving objects: 100% (25806/25806), 110.50 MiB | 9.79 MiB/s, done.
  Resolving deltas: 100% (19382/19382), done.




