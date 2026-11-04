=======================
Biofold plugin
=======================

This is a **Scipion** plugin for protein / nucleic-acid structure prediction ("folding").
It can **run** several diffusion / deep-learning models locally and/or **import** results
produced by their public web servers:
`Boltz-2 <https://github.com/jwohlwend/boltz>`_,
`Chai-1 <https://github.com/chaidiscovery/chai-lab>`_,
`IntelliFold <https://github.com/IntelliGen-AI/IntelliFold>`_,
`Protenix <https://github.com/bytedance/Protenix>`_ and
`AlphaFold3 <https://github.com/google-deepmind/alphafold3>`_.


==================================
Supported models and capabilities
==================================

Each model is exposed through its own protocol. **Run** means the prediction is computed
locally (the plugin installs a dedicated conda environment for the model and execution uses
the GPU). **Import** means results generated on the model's public web server are parsed into
Scipion objects.

- **Boltz-2** (v2.2.1) — *run and import*. Protein / DNA / RNA; optional binding-affinity
  estimation; controls for recycling steps, sampling steps, diffusion samples, step scale and
  inference potentials; cyclic entities.
- **Chai-1** (v0.6.1) — *run and import*. Protein / DNA / RNA; optional MSA server.
- **IntelliFold** (v2.0.0) — *run only*. Protein / DNA / RNA; selectable model
  (``v1`` / ``v2`` / ``v2-flash``); MSA toggle; recycling and sampling steps; number of
  output models.
- **Protenix** (v2.0.0) — *run and import*. Selectable model checkpoint
  (``base_default_v1.0.0`` / ``base_20250630_v1.0.0`` / ``base_default_v0.5.0``); MSA toggle;
  sampling steps; number of output models.
- **AlphaFold3** — *import only* (no local run); results imported from the AlphaFold3 server.

**Inputs for local runs:** a Sequence, a SetOfSequences, an AtomStruct (a chain and residue
positions are extracted from it), or a FASTA file. Runs execute on GPU (GPU IDs are
selectable in the form).

**Importing predictions** — the *import structure predictions* protocol parses result
archives into Scipion ``AtomStruct`` objects (each tagged with its origin server) from:

- AlphaFold3 server — https://alphafoldserver.com/
- Protenix server — https://protenix-server.com/
- Chai server — https://lab.chaidiscovery.com/
- Boltz server (Tamarind) — https://app.tamarind.bio/boltz/


==========================
Install this plugin
==========================

You will need to first install
`Scipion3 <https://scipion-em.github.io/docs/release-3.0.0/docs/scipion-modes/how-to-install.html>`_  and
`Scipion-chem <https://github.com/scipion-chem/scipion-chem>`_ to run these protocols.


1. **Install the plugin in Scipion**

Biofold is installed automatically by scipion.

- **Install the stable version (Not available yet)**

    Through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

    or

.. code-block::

    scipion3 installp -p scipion-chem-biofold


- **Developer's version**

    1. **Download repository**:

    .. code-block::

        git clone https://github.com/scipion-chem/scipion-chem-biofold.git

    2. **Switch to the desired branch** (master or devel):

    Scipion-chem-biofold is constantly under development and including new features.
    If you want a relatively older an more stable version, use master branch (default).
    If you want the latest changes and developments, user devel branch.

    .. code-block::

                cd scipion-chem-biofold
                git checkout devel

    3. **Install**:

    .. code-block::

        scipion3 installp -p path_to_scipion-chem-biofold --devel




