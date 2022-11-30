
.. _git_behind_firewall:

Accessing github behind a firewall
-----------------------------------

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


.. _pypi_behind_firewall:

Accessing PyPI behind a firewall
--------------------------------

Many post-processing scripts of SeisSol require Python dependencies.
We describe how to use pip on a HPC cluster with restricts access to outside sources in the following.


1. On your local machine in ~/.ssh/config add the following `RemoteForward` line:

::

    Host supermucNG
        ...
        RemoteForward ddddd localhost:8899

where ddddd is an arbitrary port number with 5 digits.
(This number should be different from port number used in other RemoteForward entries.)

2. Install proxy.py on your local machine.

::

    pip install --upgrade --user proxy.py

3. Start proxy.py on your local machine. (And keep it running.)


::

    ~/.local/bin/proxy --port 8899

4. Login to the HPC cluster (e.g. `ssh supermucNG`). Pip can be used with

::

    pip install <package name> --user --proxy localhost:ddddd

where ddddd is your arbitrary port number.

